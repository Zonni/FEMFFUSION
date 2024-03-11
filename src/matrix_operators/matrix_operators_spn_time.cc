/**
 *
 * @file   matrix_operators_spn_time.cc
 * @brief  Implementation of SystemMatrixTimeSPN and FissionMAtrix classes to handle block matrices.
 *
 */

#include "../../include/matrix_operators/matrix_operators_spn_time.h"
#include "../../include/matrix_operators/matrix_operators_base.h"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <algorithm>
#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  SystemMatrixTimeSPN<dim, n_fe_degree>::SystemMatrixTimeSPN (const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      TransportMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(
        dof_handler),
      constraints(constraints)
  {
    // Silly initializations
    n_groups = 0;
    this->type_approximation = "spn";
    listen_to_material_id = false;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SystemMatrixTimeSPN<dim, n_fe_degree>::reinit (const Materials &materials,
    const unsigned int _n_moments,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const double deltat,
    const std::string timescheme,
    const MatrixFreeType &_matrixfree_type,
    bool _listen_to_material_id)
  {
    this->matrixfree_type = _matrixfree_type;
    n_groups = materials.get_n_groups();
    this->n_moments = _n_moments;
    this->n_blocks = this->n_moments * n_groups;

    std::cout << "n_blocks: " << this->n_blocks << std::endl;
    this->n_dofs_block = dof_handler.n_dofs();

    this->delta_t = deltat;
    this->type_scheme = timescheme;

    this->boundary_conditions = boundary_conditions;
    this->albedo_factors = albedo_factors;

    listen_to_material_id = _listen_to_material_id;

    if (this->matrixfree_type == non_diagonal)
    {
      reinit_non_diagonal(materials, listen_to_material_id,
        boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      reinit_full_matrixfree(materials, listen_to_material_id,
        boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(this->n_blocks,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();

      DynamicSparsityPattern dsp(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->dof_handler.n_locally_owned_dofs_per_processor(),
        this->comm, this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int b1 = 0; b1 < this->n_blocks; b1++)
        for (unsigned int b2 = 0; b2 < this->n_blocks; b2++)
        {
          this->matrix_blocks[b1][b2] =
                                        new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[b1][b2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs, this->sp, this->comm);
        }

      assemble_full_matrices(materials, boundary_conditions, albedo_factors);
    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SystemMatrixTimeSPN<dim, n_fe_degree>::reinit_non_diagonal (
    const Materials &materials,
    bool listen_to_material_id,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler,
      this->locally_relevant_dofs);

    // Build Transport mass coefficients
    coeffs_L.resize(n_groups,
      std::vector<std::vector<std::vector<Vector<double> > > >(n_groups,
        std::vector<std::vector<Vector<double> > >(this->n_moments,
          std::vector<Vector<double> >(this->n_moments,
            Vector<double>(materials.get_n_mats())))));

    std::vector<double> exp_value(materials.get_n_precursors());

    if (this->type_scheme == "implicit-exponential")
    {
      for (unsigned int np = 0; np < materials.get_n_precursors(); np++)
        exp_value[np] = materials.get_delayed_fraction(0, np)
                        * (1
                           - exp(
                             -materials.get_delayed_decay_constant(0, np)
                             * this->delta_t));
    }

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        for (unsigned int mi = 0; mi < this->n_moments; mi++)
          for (unsigned int mj = 0; mj < this->n_moments; mj++)
          {
            coeffs_L[gi][gj][mi][mj] = Vector<double>(materials.get_n_mats());
            for (unsigned int mat = 0; mat < materials.get_n_mats();
                mat++)
            {
              coeffs_L[gi][gj][mi][mj][mat] = 0.0;
              if (gi == gj)
              {
                coeffs_L[gi][gj][mi][mj][mat] = sp_coeff[0][mi][mj]
                                                * materials.get_sigma_r(gi, mat);
                for (unsigned int m3 = 1; m3 < this->n_moments; m3++)
                  coeffs_L[gi][gj][mi][mj][mat] +=
                      sp_coeff[m3][mi][mj]
                      * (materials.get_sigma_r(gi, mat)
                         + materials.get_sigma_s(gi, gi, mat));
                // velocity
                for (unsigned int m3 = 0; m3 < this->n_moments; m3++)
                  coeffs_L[gi][gj][mi][mj][mat] +=
                                                   sp_coeff[m3][mi][mj] / this->delta_t
                                                   / materials.get_velocity(mat,
                                                     gi);

              }
              else
              {
                coeffs_L[gi][gj][mi][mj][mat] = -sp_coeff[0][mi][mj]
                                                * materials.get_sigma_s(gj, gi, mat);
              }

              // fission
              coeffs_L[gi][gj][mi][mj][mat] += -sp_coeff[0][mi][mj]
                  * (1 - materials.get_delayed_fraction_sum(mat))
                  * materials.get_prompt_spectra(mat, gi)
                  * materials.get_nu_sigma_f(gj, mat);

              for (unsigned int np = 0;
                  np < materials.get_n_precursors(); np++)
              {
                coeffs_L[gi][gj][mi][mj][mat] +=
                                                 -sp_coeff[0][mi][mj] * exp_value[np]
                                                 * materials.get_prompt_spectra(mat,
                                                   gi)
                                                 * materials.get_nu_sigma_f(gj, mat);
              }

            }

          }

    //  --------- Matrix-Free Blocks  ---------

    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_JxW_values);
    additional_data.mapping_update_flags_boundary_faces = (update_values
                                                           | update_JxW_values);

    matfree_data.reinit(dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(
        this->n_blocks));
    unsigned int bi, bj;
    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int mi = 0; mi < this->n_moments; mi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mj = 0; mj < this->n_moments; mj++)
          {
            bi = gm_to_b(gi, mi);
            bj = gm_to_b(gj, mj);
            if (bi != bj)
            {
              this->mass_mf_blocks[bi][bj] = new MassOperator<dim,
                  n_fe_degree, double>(matfree_data);
              if (gi == gj)
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  boundary_conditions, marshack_coeff[mi][mj],
                  listen_to_material_id); // With BC
              else
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  listen_to_material_id); // Without BC
            }
//            std::cout << "L(" << bi << ", " << bj << ") = " << std::endl;
          }

    //  --------- Diagonal Blocks  ---------

    this->matrix_blocks.resize(this->n_blocks,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));

    DynamicSparsityPattern dsp(this->locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

    SparsityTools::distribute_sparsity_pattern(dsp,
      this->dof_handler.n_locally_owned_dofs_per_processor(), this->comm,
      this->locally_relevant_dofs);
    this->sp.copy_from(dsp);

    for (unsigned int b = 0; b < this->n_blocks; b++)
    {
      this->matrix_blocks[b][b] = new PETScWrappers::MPI::SparseMatrix;
      this->matrix_blocks[b][b]->reinit(this->locally_owned_dofs,
        this->locally_owned_dofs, this->sp, this->comm);
    }

    assemble_diagonal_matrices(materials, boundary_conditions, albedo_factors);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SystemMatrixTimeSPN<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials,
    bool listen_to_material_id,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler,
      this->locally_relevant_dofs);

    // Build Transport mass coefficients
    coeffs_L.resize(n_groups,
      std::vector<std::vector<std::vector<Vector<double> > > >(n_groups,
        std::vector<std::vector<Vector<double> > >(this->n_moments,
          std::vector<Vector<double> >(this->n_moments,
            Vector<double>(materials.get_n_mats())))));

    std::vector<double> exp_value(materials.get_n_precursors());

    if (this->type_scheme == "implicit-exponential")
    {
      for (unsigned int np = 0; np < materials.get_n_precursors(); np++)
        exp_value[np] = materials.get_delayed_fraction(0, np)
                        * (1
                           - exp(
                             -materials.get_delayed_decay_constant(0, np)
                             * this->delta_t));
    }

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        for (unsigned int mi = 0; mi < this->n_moments; mi++)
          for (unsigned int mj = 0; mj < this->n_moments; mj++)
          {
            coeffs_L[gi][gj][mi][mj] = Vector<double>(materials.get_n_mats());
            for (unsigned int mat = 0; mat < materials.get_n_mats();
                mat++)
            {
              coeffs_L[gi][gj][mi][mj][mat] = 0.0;
              if (gi == gj)
              {
                coeffs_L[gi][gj][mi][mj][mat] = sp_coeff[0][mi][mj]
                                                * materials.get_sigma_r(gi, mat);
                for (unsigned int m3 = 1; m3 < this->n_moments; m3++)
                  coeffs_L[gi][gj][mi][mj][mat] +=
                      sp_coeff[m3][mi][mj]
                      * (materials.get_sigma_r(gi, mat)
                         + materials.get_sigma_s(gi, gi, mat));
                // velocity
                for (unsigned int m3 = 0; m3 < this->n_moments; m3++)
                  coeffs_L[gi][gj][mi][mj][mat] +=
                                                   sp_coeff[m3][mi][mj] / this->delta_t
                                                   / materials.get_velocity(mat,
                                                     gi);

              }
              else
              {
                coeffs_L[gi][gj][mi][mj][mat] = -sp_coeff[0][mi][mj]
                                                * materials.get_sigma_s(gj, gi, mat);
              }

              // fission
              coeffs_L[gi][gj][mi][mj][mat] += -sp_coeff[0][mi][mj]
                  * (1 - materials.get_delayed_fraction_sum(mat))
                  * materials.get_prompt_spectra(mat, gi)
                  * materials.get_nu_sigma_f(gj, mat);

              for (unsigned int np = 0;
                  np < materials.get_n_precursors(); np++)
              {
                coeffs_L[gi][gj][mi][mj][mat] +=
                                                 -sp_coeff[0][mi][mj] * exp_value[np]
                                                 * materials.get_prompt_spectra(mat,
                                                   gi)
                                                 * materials.get_nu_sigma_f(gj, mat);
              }

            }

          }

    //  --------- Matrix-Free Blocks  ---------

    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_JxW_values);
    additional_data.mapping_update_flags_boundary_faces = (update_values
                                                           | update_JxW_values);

    matfree_data.reinit(dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(
        this->n_blocks));
    unsigned int bi, bj;
    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int mi = 0; mi < this->n_moments; mi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mj = 0; mj < this->n_moments; mj++)
          {
            bi = gm_to_b(gi, mi);
            bj = gm_to_b(gj, mj);
            if (bi != bj)
            {
              this->mass_mf_blocks[bi][bj] = new MassOperator<dim,
                  n_fe_degree, double>(matfree_data);
              if (gi == gj)
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  boundary_conditions, marshack_coeff[mi][mj],
                  listen_to_material_id); // With BC
              else
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  listen_to_material_id); // Without BC
            }
//            std::cout << "L(" << bi << ", " << bj << ") = " << std::endl;
          }

    // Diagonal blocks

    coeffs_grad.resize(this->n_blocks);
    coeffs_val.resize(this->n_blocks);
    coeffs_bound.resize(this->n_blocks);

    unsigned int b;
    for (unsigned int g = 0; g < n_groups; g++)
      for (unsigned int m = 0; m < this->n_moments; m++)
      {
        b = gm_to_b(g, m);
        coeffs_grad[b].reinit(materials.get_n_mats());
        coeffs_val[b].reinit(materials.get_n_mats());
        coeffs_bound[b] = marshack_coeff[m][m];
        for (unsigned int mat = 0; mat < materials.get_n_mats(); mat++)
        {
          coeffs_grad[b][mat] = diff_coeff[m]
                                / (materials.get_sigma_tr(g, mat));
          coeffs_val[b][mat] = coeffs_L[g][g][m][m][mat];
        }
      }

    this->poison_mf_blocks.resize(this->n_blocks);

    // Diagonal blocks
    for (unsigned int b = 0; b < this->n_blocks; b++)
    {
      this->poison_mf_blocks[b] =
                                  new PoissonOperator<dim, n_fe_degree, double>(
                                    matfree_data);

      this->poison_mf_blocks[b]->reinit(b, this->constraints, materials,
        materials.get_materials_vector(), coeffs_val[b],
        coeffs_grad[b], boundary_conditions, albedo_factors, true,
        coeffs_bound[b]);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  unsigned int SystemMatrixTimeSPN<dim, n_fe_degree>::gm_to_b (
    const unsigned int group,
    const unsigned int moment) const
  {
    return this->n_moments * group + moment;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SystemMatrixTimeSPN<dim, n_fe_degree>::assemble_diagonal_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma, bc_coeff;
    double factor = 0.0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_bound(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    bool show_materials = false;
    bool show_boundary = false;
    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        const unsigned int mat = materials.get_material_id<dim>(cell);

        cell_grad = 0;
        cell_val = 0;
        cell_bound = 0;

        //
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
            {
              val = fe_values.shape_value(i, q_p)
                    * fe_values.shape_value(j, q_p)
                    * fe_values.JxW(q_p);

              cell_val(i, j) += val;

              grad = fe_values.shape_grad(i, q_p)
                     * fe_values.shape_grad(j, q_p)
                     * fe_values.JxW(q_p);

              cell_grad(i, j) += grad;
            }

        // Take care of albedo Boundary Conditions: Boundary integral
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
          {
            types::boundary_id boundary_id =
                                             cell->face(f)->boundary_id();
            AssertIndexRange(boundary_id, boundary_conditions.size());

            if (boundary_conditions[boundary_id] > 1)
            {
              fe_face_values.reinit(cell, f);
              switch (boundary_conditions[boundary_id])
              {
                case 2: // Vacuum BC
                  factor = 1.0;
                  break;
                default: // Custom Albedo BC
                  factor =
                           albedo_factors[boundary_conditions[boundary_id]
                                          - 3];
                  break;
              }

              for (unsigned int q_p = 0; q_p < n_face_q_points; ++q_p)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < dofs_per_cell;
                      ++j)
                  {
                    val = fe_face_values.shape_value(i, q_p)
                          * fe_face_values.shape_value(j, q_p)
                          * fe_face_values.JxW(q_p);

                    cell_bound(i, j) += factor * val;
                  }
            }
          }

        // Distribute in the Sparse matrix

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
          for (unsigned int mi = 0; mi < this->n_moments; ++mi)
          {
            // Get the material coefficients:
            D = diff_coeff[mi] / (materials.get_sigma_tr(gi, mat));
            sigma = coeffs_L[gi][gi][mi][mi][mat];
            bc_coeff = marshack_coeff[mi][mi];
            cell_L.equ(D, cell_grad, sigma, cell_val, bc_coeff,
              cell_bound);
            constraints.distribute_local_to_global(cell_L,
              local_dof_indices,
              *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(gi,
                mi)]));
          }

      }

    // Compress the matrices
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      this->matrix_blocks[bi][bi]->compress(VectorOperation::add);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SystemMatrixTimeSPN<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma, bc_coeff;
    double factor = 0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_bound(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    bool show_materials = false;
    bool show_boundary = false;
    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      const unsigned int mat = materials.get_material_id<dim>(cell);

      cell_grad = 0;
      cell_val = 0;
      cell_bound = 0;

      //
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
          {
            val = fe_values.shape_value(i, q_p)
                  * fe_values.shape_value(j, q_p)
                  * fe_values.JxW(q_p);

            cell_val(i, j) += val;

            grad = fe_values.shape_grad(i, q_p)
                   * fe_values.shape_grad(j, q_p)
                   * fe_values.JxW(q_p);

            cell_grad(i, j) += grad;
          }

      // Take care of albedo Boundary Conditions: Boundary integral
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
        {
          types::boundary_id boundary_id = cell->face(f)->boundary_id();
          AssertIndexRange(boundary_id, boundary_conditions.size());

          if (boundary_conditions[boundary_id] > 1)
          {
            fe_face_values.reinit(cell, f);
            switch (boundary_conditions[boundary_id])
            {
              case 2: // Vacuum BC
                factor = 1.0;
                break;
              default: // Custom Albedo BC
                factor = albedo_factors[boundary_conditions[boundary_id]
                                        - 3];
                break;
            }

            for (unsigned int q_p = 0; q_p < n_face_q_points; ++q_p)
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  val = fe_face_values.shape_value(i, q_p)
                        * fe_face_values.shape_value(j, q_p)
                        * fe_face_values.JxW(q_p);

                  cell_bound(i, j) += factor * val;
                }
          }
        }

      std::vector<double> exp_value(materials.get_n_precursors());

      if (this->type_scheme == "implicit-exponential")
      {
        for (unsigned int np = 0; np < materials.get_n_precursors(); np++)
          exp_value[np] = materials.get_delayed_fraction(0, np)
                          * (1
                             - exp(
                               -materials.get_delayed_decay_constant(0,
                                 np)
                               * this->delta_t));
      }

      // Distribute in the Sparse matrix
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          for (unsigned int mi = 0; mi < this->n_moments; ++mi)
            for (unsigned int mj = 0; mj < this->n_moments; ++mj)
            {
              // Get the material coefficients:
              if (gi == gj)
              {
                if (mi == mj)
                {
                  D = diff_coeff[mi]
                      / (materials.get_sigma_tr(gi, mat));
                  sigma = sp_coeff[0][mi][mj]
                          * materials.get_sigma_r(gi, mat);
                  for (unsigned int m3 = 1; m3 < this->n_moments; ++m3)
                    sigma += sp_coeff[m3][mi][mj]
                             * (materials.get_sigma_r(gi, mat)
                                + materials.get_sigma_s(gi, gi, mat));

                  for (unsigned int m3 = 0; m3 < this->n_moments; ++m3)
                    sigma += sp_coeff[m3][mi][mj] / this->delta_t
                             / materials.get_velocity(mat, gi);

                  sigma += -sp_coeff[0][mi][mj] *
                           (1 - materials.get_delayed_fraction_sum(mat))
                           * materials.get_prompt_spectra(mat, gi)
                           * materials.get_nu_sigma_f(gi, mat);

                  for (unsigned int np = 0;
                      np < materials.get_n_precursors(); np++)

                    sigma += -sp_coeff[0][mi][mj] * exp_value[np]
                             * materials.get_prompt_spectra(mat,
                               gi)
                             * materials.get_nu_sigma_f(gj, mat);

                  bc_coeff = marshack_coeff[mi][mj];
                  cell_L.equ(D, cell_grad, sigma, cell_val,
                    bc_coeff, cell_bound);
                  constraints.distribute_local_to_global(cell_L,
                    local_dof_indices,
                    *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(
                      gj, mj)]));
                }
                else
                {
                  sigma = sp_coeff[0][mi][mj] * materials.get_sigma_r(gi, mat);
                  for (unsigned int m3 = 1; m3 < this->n_moments; ++m3)
                  {
                    sigma += sp_coeff[m3][mi][mj]
                             * (materials.get_sigma_r(gi, mat)
                                + materials.get_sigma_s(gi, gi, mat));

                  }
                  for (unsigned int m3 = 0; m3 < this->n_moments; ++m3)
                    sigma += sp_coeff[m3][mi][mj] / this->delta_t
                             / materials.get_velocity(mat, gi);

                  sigma += -sp_coeff[0][mi][mj]
                           * (1 - materials.get_delayed_fraction_sum(mat))
                           * materials.get_prompt_spectra(mat, gi)
                           * materials.get_nu_sigma_f(gj, mat);
                  for (unsigned int np = 0;
                      np < materials.get_n_precursors(); np++)
                    // TODO check to change get_prompt_spectra by get_delayed_spectra
                    sigma += -sp_coeff[0][mi][mj] * exp_value[np]
                             * materials.get_prompt_spectra(mat, gi)
                             * materials.get_nu_sigma_f(gj, mat);

                  bc_coeff = marshack_coeff[mi][mj];
                  cell_L.equ(sigma, cell_val, bc_coeff,
                    cell_bound);
                  constraints.distribute_local_to_global(cell_L,
                    local_dof_indices,
                    *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(
                      gj, mj)]));
                }
              }
              else
              {

                sigma = -sp_coeff[0][mi][mj] * materials.get_sigma_s(gj, gi, mat);

                sigma += -sp_coeff[0][mi][mj]
                         * (1 - materials.get_delayed_fraction_sum(mat))
                         * materials.get_prompt_spectra(mat, gi)
                         * materials.get_nu_sigma_f(gj, mat);

                for (unsigned int np = 0;
                    np < materials.get_n_precursors(); np++)
                  // TODO check to change get_prompt_spectra by get_delayed_spectra
                  sigma += -exp_value[np]
                           * materials.get_prompt_spectra(mat, gi)
                           * materials.get_nu_sigma_f(gj, mat);

                cell_L.equ(sigma, cell_val);
                constraints.distribute_local_to_global(cell_L,
                  local_dof_indices,
                  *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(
                    gj, mj)]));

              }
            }
    }

    //
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks; ++bj)
        this->matrix_blocks[bi][bj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t SystemMatrixTimeSPN<dim, n_fe_degree>::memory_consumption () const
  {

    std::size_t memory = 0;
    if (this->matrixfree_type == non_diagonal)
    {
      memory = matfree_data.memory_consumption();
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
      {
        memory += this->matrix_blocks[i][i]->memory_consumption();
      }
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      memory = this->matfree_data.memory_consumption();
    }
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- FissionMatrix    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  FisionMatrixTimeSPN<dim, n_fe_degree>::FisionMatrixTimeSPN (const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(
        dof_handler),
      constraints(constraints)
  {
    n_moments = 0;
  }

/**
 * @brief Star the Operator
 */
template <int dim, int n_fe_degree>
  void FisionMatrixTimeSPN<dim, n_fe_degree>::reinit (const Materials &materials,
    const unsigned int _n_moments,
    const MatrixFreeType &_matrixfree_type,
    bool listen_to_material_id)
  {
    const unsigned int n_groups = materials.get_n_groups();
    n_moments = _n_moments;
    this->n_blocks = _n_moments * n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(
          this->n_blocks));
      coeffs_M.resize(n_groups,
        std::vector<std::vector<std::vector<Vector<double> > > >(
          n_groups,
          std::vector<std::vector<Vector<double> > >(n_moments,
            std::vector<Vector<double> >(n_moments,
              Vector<double>(n_mats)))));
      //
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
          for (unsigned int mi = 0; mi < n_moments; mi++)
            for (unsigned int mj = 0; mj < n_moments; mj++)
              for (unsigned int mat = 0; mat < n_mats; mat++)
              {
                coeffs_M[to_g][from_g][mi][mj][mat] =
                                                      sp_coeff[0][mi][mj]
                                                      * materials.get_xi_nu_sigma_f(
                                                        from_g, to_g, mat);
              }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::partition_partition;
      additional_data.mapping_update_flags = (update_values
                                              | update_JxW_values);
      matfree_data.reinit(dof_handler, constraints,
        QGauss<1>(n_fe_degree + 1), additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mi = 0; mi < n_moments; mi++)
            for (unsigned int mj = 0; mj < n_moments; mj++)
            {
              this->mass_mf_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)] =
                                                                       new MassOperator<
                                                                           dim,
                                                                           n_fe_degree,
                                                                           double>(
                                                                         matfree_data);
              this->mass_mf_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)]->reinit(
                constraints, materials.get_materials_vector(),
                coeffs_M[gi][gj][mi][mj],
                listen_to_material_id);
            }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(this->n_blocks,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));

      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->dof_handler.n_locally_owned_dofs_per_processor(),
        this->comm, this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int bi = 0; bi < this->n_blocks; bi++)
        for (unsigned int bj = 0; bj < this->n_blocks; bj++)
        {
          this->matrix_blocks[bi][bj] =
                                        new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[bi][bj]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs, this->sp, this->comm);
        }

      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void FisionMatrixTimeSPN<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
    double xi_nu_sigma_f;
    unsigned int mat;
    std::vector<types::global_dof_index> local_dof_indices(
      dof_handler.get_fe().dofs_per_cell);
    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_M(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      mat = materials.get_material_id<dim>(cell);

      cell_val = 0;
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
          {
            val = fe_values.shape_value(i, q_p)
                  * fe_values.shape_value(j, q_p)
                  * fe_values.JxW(q_p);

            cell_val(i, j) += val;
          }

      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          for (unsigned int mi = 0; mi < n_moments; ++mi)
            for (unsigned int mj = 0; mj < n_moments; ++mj)
            {
              // Get the material coefficients:
              xi_nu_sigma_f = sp_coeff[0][mi][mj]
                              * materials.get_xi_nu_sigma_f(gj, gi, mat);
              cell_M.equ(xi_nu_sigma_f, cell_val);
              constraints.distribute_local_to_global(cell_M,
                local_dof_indices,
                *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(
                  gj, mj)]));
            }
    }

    // Compress
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks; ++bj)
        this->matrix_blocks[bi][bj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  unsigned int FisionMatrixTimeSPN<dim, n_fe_degree>::gm_to_b (
    const unsigned int group,
    const unsigned int moment) const
  {
    return n_moments * group + moment;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t FisionMatrixTimeSPN<dim, n_fe_degree>::memory_consumption () const
  {
    return 0;
    /*
     std::size_t memory = 0;
     if (this->matrixfree_type == non_diagonal)
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     {
     if (i == j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     else
     memory += this->mass_mf_blocks[i][j]->memory_consumption();
     }
     else
     {
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     }
     return memory;
     */
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- FissionMatrix    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  MassMatrixTimeSPN<dim, n_fe_degree>::MassMatrixTimeSPN (const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(
        dof_handler),
      constraints(constraints)
  {
    n_moments = 0;
  }

/**
 * @brief Star the Operator
 */
template <int dim, int n_fe_degree>
  void MassMatrixTimeSPN<dim, n_fe_degree>::reinit (const Materials &materials,
    const unsigned int _n_moments,
    const MatrixFreeType &_matrixfree_type,
    bool listen_to_material_id)
  {
    const unsigned int n_groups = materials.get_n_groups();
    n_moments = _n_moments;
    this->n_blocks = _n_moments * n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(
          this->n_blocks));
      coeffs_M.resize(n_groups,
        std::vector<std::vector<std::vector<Vector<double> > > >(
          n_groups,
          std::vector<std::vector<Vector<double> > >(n_moments,
            std::vector<Vector<double> >(n_moments,
              Vector<double>(n_mats)))));
      //
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
          for (unsigned int mi = 0; mi < n_moments; mi++)
            for (unsigned int mj = 0; mj < n_moments; mj++)
            {
              if (from_g == to_g)
              {
                coeffs_M[to_g][from_g][mi][mj] = Vector<double>(n_mats);
                for (unsigned int mat = 0; mat < n_mats; mat++)
                {
                  coeffs_M[to_g][from_g][mi][mj][mat] = 0.0;
                  for (unsigned int nm = 0; nm < n_moments; nm++)
                  {
                    coeffs_M[to_g][from_g][mi][mj][mat] +=
                                                           sp_coeff[nm][mi][mj]
                                                           / materials.get_velocity(
                                                             mat, from_g);
                  }
                }

              }
            }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::partition_partition;

      additional_data.mapping_update_flags = (update_values
                                              | update_JxW_values);

      matfree_data.reinit(dof_handler, constraints,
        QGauss<1>(n_fe_degree + 1), additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mi = 0; mi < n_moments; mi++)
            for (unsigned int mj = 0; mj < n_moments; mj++)
            {
              this->mass_mf_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)] =
                                                                       new MassOperator<
                                                                           dim,
                                                                           n_fe_degree,
                                                                           double>(
                                                                         matfree_data);
              this->mass_mf_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)]->reinit(
                constraints, materials.get_materials_vector(),
                coeffs_M[gi][gj][mi][mj],
                listen_to_material_id);
            }

    }
    else if (this->matrixfree_type == full_allocated)
    {

      // Resize matrix_blocks
      this->matrix_blocks.resize(this->n_blocks,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));

      this->locally_owned_dofs = dof_handler.locally_owned_dofs();

      DynamicSparsityPattern dsp(dof_handler.n_dofs());

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->dof_handler.n_locally_owned_dofs_per_processor(),
        this->comm, this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int bi = 0; bi < this->n_blocks; bi++)
        for (unsigned int bj = 0; bj < this->n_blocks; bj++)
        {
          this->matrix_blocks[bi][bj] =
                                        new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[bi][bj]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs, this->sp, this->comm);
        }

      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void MassMatrixTimeSPN<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
    unsigned int mat;
    std::vector<types::global_dof_index> local_dof_indices(
      dof_handler.get_fe().dofs_per_cell);
    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_M(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      mat = materials.get_material_id<dim>(cell);

      cell_val = 0;
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
          {
            val = fe_values.shape_value(i, q_p)
                  * fe_values.shape_value(j, q_p)
                  * fe_values.JxW(q_p);

            cell_val(i, j) += val;
          }

      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          for (unsigned int mi = 0; mi < n_moments; ++mi)
            for (unsigned int mj = 0; mj < n_moments; ++mj)
            {
              // Get the material coefficients:
              double velocity_val = 0.0;
              if (gi == gj)
              {
                for (unsigned int nm = 0; nm < n_moments; nm++)
                  velocity_val += sp_coeff[nm][mi][mj]
                                  / materials.get_velocity(mat, gi);

                cell_M.equ(velocity_val, cell_val);
                constraints.distribute_local_to_global(cell_M,
                  local_dof_indices,
                  *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(
                    gj, mj)]));

              }
            }
    }

    // Compress
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks; ++bj)
        this->matrix_blocks[bi][bj]->compress(VectorOperation::add);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  unsigned int MassMatrixTimeSPN<dim, n_fe_degree>::gm_to_b (
    const unsigned int group,
    const unsigned int moment) const
  {
    return n_moments * group + moment;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t MassMatrixTimeSPN<dim, n_fe_degree>::memory_consumption () const
  {
    return 0;
    /*
     std::size_t memory = 0;
     if (this->matrixfree_type == non_diagonal)
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     {
     if (i == j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     else
     memory += this->mass_mf_blocks[i][j]->memory_consumption();
     }
     else
     {
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     }
     return memory;
     */
  }

// ----------- Explicit Instantations ----------- //

template class SystemMatrixTimeSPN<1, 1> ;
template class SystemMatrixTimeSPN<1, 2> ;
template class SystemMatrixTimeSPN<1, 3> ;
template class SystemMatrixTimeSPN<1, 4> ;
template class SystemMatrixTimeSPN<1, 5> ;

template class SystemMatrixTimeSPN<2, 1> ;
template class SystemMatrixTimeSPN<2, 2> ;
template class SystemMatrixTimeSPN<2, 3> ;
template class SystemMatrixTimeSPN<2, 4> ;
template class SystemMatrixTimeSPN<2, 5> ;

template class SystemMatrixTimeSPN<3, 1> ;
template class SystemMatrixTimeSPN<3, 2> ;
template class SystemMatrixTimeSPN<3, 3> ;
template class SystemMatrixTimeSPN<3, 4> ;
template class SystemMatrixTimeSPN<3, 5> ;

template class FisionMatrixTimeSPN<1, 1> ;
template class FisionMatrixTimeSPN<1, 2> ;
template class FisionMatrixTimeSPN<1, 3> ;
template class FisionMatrixTimeSPN<1, 4> ;
template class FisionMatrixTimeSPN<1, 5> ;

template class FisionMatrixTimeSPN<2, 1> ;
template class FisionMatrixTimeSPN<2, 2> ;
template class FisionMatrixTimeSPN<2, 3> ;
template class FisionMatrixTimeSPN<2, 4> ;
template class FisionMatrixTimeSPN<2, 5> ;

template class FisionMatrixTimeSPN<3, 1> ;
template class FisionMatrixTimeSPN<3, 2> ;
template class FisionMatrixTimeSPN<3, 3> ;
template class FisionMatrixTimeSPN<3, 4> ;
template class FisionMatrixTimeSPN<3, 5> ;

template class MassMatrixTimeSPN<1, 1> ;
template class MassMatrixTimeSPN<1, 2> ;
template class MassMatrixTimeSPN<1, 3> ;
template class MassMatrixTimeSPN<1, 4> ;
template class MassMatrixTimeSPN<1, 5> ;

template class MassMatrixTimeSPN<2, 1> ;
template class MassMatrixTimeSPN<2, 2> ;
template class MassMatrixTimeSPN<2, 3> ;
template class MassMatrixTimeSPN<2, 4> ;
template class MassMatrixTimeSPN<2, 5> ;

template class MassMatrixTimeSPN<3, 1> ;
template class MassMatrixTimeSPN<3, 2> ;
template class MassMatrixTimeSPN<3, 3> ;
template class MassMatrixTimeSPN<3, 4> ;
template class MassMatrixTimeSPN<3, 5> ;
