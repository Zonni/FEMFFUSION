/**
 * @file   matrix_operators_noise_diff.cc
 * @brief
 */

#include "matrix_operators_noise_diff.h"
#include "matrix_operators_base.h"

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

#include <iostream>
#include <fstream>
#include <sstream>

using namespace dealii;

/**
 * @brief Calu
 */
void calculate_precursors_factor (
  const unsigned int mat_id,
  const Materials &materials,
  const double &omega,
  std::vector<complex> &prec_factor)
{
  //std::vector<complex> prec_factor(materials.get_n_groups());

  for (unsigned int g = 0; g < materials.get_n_groups(); g++)
  {
    prec_factor[g] = (1 - materials.get_delayed_fraction_sum(mat_id)) *
                     materials.get_prompt_spectra(mat_id, g);

    for (unsigned int p = 0; p < materials.get_n_precursors(); p++)
    {
      const complex divisor(materials.get_delayed_decay_constant(mat_id, p), omega);
      prec_factor[g] += materials.get_delayed_decay_constant(mat_id, p)
                        * materials.get_delayed_fraction(mat_id, p)
                        / (divisor)
                        * materials.get_delayed_spectra(mat_id, p, g);
    }
  }
}

/**
 *
 */
template <int dim, int n_fe_degree>
  NoiseAMatrix<dim, n_fe_degree>::NoiseAMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double> &constraints) :
      TransportMatrixComplexBase<dim, n_fe_degree>(comm, dh, constraints),
      tria(dh.get_triangulation()),
      dof_handler(dh),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const ComplexPerturbation &pert,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const MatrixFreeType &_matrixfree_type)
  {

    this->matrixfree_type = _matrixfree_type;
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = 2 * n_groups;
    this->n_blocks_real = n_groups;
    this->n_dofs_block = dof_handler.n_dofs();
    this->boundary_conditions = boundary_conditions;
    this->albedo_factors = albedo_factors;

    if (this->matrixfree_type == non_diagonal)
    {
      reinit_non_diagonal(materials, pert, boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      reinit_full_matrixfree(materials, pert, boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks_real.resize(this->n_blocks_real,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks_real));
      this->matrix_blocks_imag.resize(this->n_blocks_real,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks_real));

      DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();

      DynamicSparsityPattern dsp(this->locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->dof_handler.n_locally_owned_dofs_per_processor(),
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int b1 = 0; b1 < this->n_blocks_real; b1++)
        for (unsigned int b2 = 0; b2 < this->n_blocks_real; b2++)
        {
          this->matrix_blocks_real[b1][b2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks_real[b1][b2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);

          this->matrix_blocks_imag[b1][b2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks_imag[b1][b2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }

      assemble_full_matrices(materials, pert, boundary_conditions, albedo_factors);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::calculate_A_factor (
    const unsigned int mat_id,
    const Materials &materials,
    const ComplexPerturbation &pert,
    std::vector<std::vector<complex> > &A_factor)
  {

    std::vector<complex> prec_factor(materials.get_n_groups());
    const double omega = pert.get_frequency() * 2 * M_PI;

    calculate_precursors_factor(mat_id, materials, omega, prec_factor);

    for (unsigned int g1 = 0; g1 < materials.get_n_groups(); g1++)
      for (unsigned int g2 = 0; g2 < materials.get_n_groups(); g2++)
      {
        if (g1 == g2)
        {
          A_factor[g1][g1] = complex(0, omega / materials.get_velocitiy(mat_id, g1))
                             + materials.get_sigma_r(g1, mat_id);
        }
        else
        {
          A_factor[g1][g2] = -materials.get_sigma_s(g2, g1, mat_id);
        }

        A_factor[g1][g2] -= prec_factor[g1]
                            * materials.get_nu_sigma_f(g2, mat_id);

      }

  }

/**
 *
 */

template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::reinit_non_diagonal (
    const Materials &materials,
    const ComplexPerturbation &pert,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {

    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    //  --------- Matrix-Free Blocks  ---------
    this->mass_mf_blocks_real.resize(this->n_blocks_real,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks_real));
    this->mass_mf_blocks_imag.resize(this->n_blocks_real,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks_real));

    coeffs_real.resize(n_groups,
      std::vector<Vector<double> >(n_groups));
    coeffs_imag.resize(n_groups,
      std::vector<Vector<double> >(n_groups));

    std::vector<std::vector<complex>> A_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));

    for (unsigned int g1 = 0; g1 < n_groups; g1++)
      for (unsigned int g2 = 0; g2 < n_groups; g2++)
      {
        coeffs_real[g1][g2].reinit(n_mats);
        coeffs_imag[g1][g2].reinit(n_mats);
      }

    // Fill coeffs
    for (unsigned int mat = 0; mat < n_mats; mat++)
    {
      calculate_A_factor(mat, materials, pert, A_factor);

      // Calculate coeffs
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          coeffs_real[g1][g2][mat] = A_factor[g1][g2].real();
          coeffs_imag[g1][g2][mat] = A_factor[g1][g2].imag();
        }
    }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values
                                            | update_JxW_values);
    matfree_data.reinit(dof_handler, constraints,
      QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
      {
        this->mass_mf_blocks_real[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
          matfree_data);
        this->mass_mf_blocks_real[gi][gj]->reinit(constraints,
          materials.get_materials_vector(),
          coeffs_real[gi][gj],
          materials.listen_to_material_id);

        this->mass_mf_blocks_imag[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
          matfree_data);
        this->mass_mf_blocks_imag[gi][gj]->reinit(constraints,
          materials.get_materials_vector(),
          coeffs_imag[gi][gj],
          materials.listen_to_material_id);
      }

    //  --------- Real Matrices of the diagonal  ---------
    // Resize matrix_blocks
    this->matrix_blocks_real.resize(this->n_blocks_real,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks_real));

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DynamicSparsityPattern dsp(this->locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

    SparsityTools::distribute_sparsity_pattern(dsp,
      dof_handler.n_locally_owned_dofs_per_processor(),
      this->comm,
      this->locally_relevant_dofs);
    this->sp.copy_from(dsp);

    for (unsigned int b = 0; b < this->n_blocks_real; b++)
    {
      this->matrix_blocks_real[b][b] = new PETScWrappers::MPI::SparseMatrix;
      this->matrix_blocks_real[b][b]->reinit(this->locally_owned_dofs,
        this->locally_owned_dofs,
        this->sp,
        this->comm);
    }

    assemble_diagonal_matrices(materials, pert, boundary_conditions, albedo_factors);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials,
    const ComplexPerturbation &pert,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    //  --------- Non -diagonal Blocks ---------
    this->mass_mf_blocks_real.resize(this->n_blocks_real,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks_real));
    this->mass_mf_blocks_imag.resize(this->n_blocks_real,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks_real));

    coeffs_real.resize(n_groups,
      std::vector<Vector<double> >(n_groups));
    coeffs_imag.resize(n_groups,
      std::vector<Vector<double> >(n_groups));

    std::vector<std::vector<complex>> A_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));

    for (unsigned int g1 = 0; g1 < n_groups; g1++)
      for (unsigned int g2 = 0; g2 < n_groups; g2++)
      {
        coeffs_real[g1][g2].reinit(n_mats);
        coeffs_imag[g1][g2].reinit(n_mats);
      }

    // Fill coeffs
    for (unsigned int mat = 0; mat < n_mats; mat++)
    {
      calculate_A_factor(mat, materials, pert, A_factor);

      // Calculate coeffs
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          coeffs_real[g1][g2][mat] = A_factor[g1][g2].real();
          coeffs_imag[g1][g2][mat] = A_factor[g1][g2].imag();
        }
    }

    //  --------- Diagonal Blocks ---------
    this->poison_mf_blocks_real.resize(this->n_blocks_real);

    // Fill Diagonal coeffs -> Posion Operators
    coeffs_grad.resize(n_groups);

    // Fill diagonal
    for (unsigned int g = 0; g < n_groups; g++)
    {
      coeffs_grad[g].reinit(n_mats);

      for (unsigned int mat = 0; mat < n_mats; mat++)
      {
        coeffs_grad[g][mat] = materials.get_diffusion_coefficient(g, mat);
      }
    }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::partition_color;
    additional_data.mapping_update_flags = (update_values
                                            | update_gradients
                                            | update_JxW_values);

    bool boundary = false;
    for (unsigned int bc = 0; bc < boundary_conditions.size(); bc++)
      if (boundary_conditions[bc] > 1)
      {
        boundary = true;
        break;
      }

    if (boundary == true)
      additional_data.mapping_update_flags_boundary_faces = (update_values
                                                             | update_JxW_values
                                                             | update_quadrature_points);

    matfree_data.reinit(dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
      {
        // Imaginary block are always nondiagonal
        this->mass_mf_blocks_imag[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
          matfree_data);
        this->mass_mf_blocks_imag[gi][gj]->reinit(constraints,
          materials.get_materials_vector(),
          coeffs_imag[gi][gj],
          materials.listen_to_material_id);

        if (gi != gj) // Non diagonal
        {
          this->mass_mf_blocks_real[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
            matfree_data);
          this->mass_mf_blocks_real[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs_real[gi][gj],
            materials.listen_to_material_id);

        }
        else // Diagonal
        {
          this->poison_mf_blocks_real[gi] = new PoissonOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->poison_mf_blocks_real[gi]->reinit(gi,
            this->constraints,
            materials,
            materials.get_materials_vector(),
            this->coeffs_real[gi][gi],
            this->coeffs_grad[gi],
            boundary_conditions,
            albedo_factors);
        }
      }

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::assemble_diagonal_matrices (
    const Materials &materials,
    const ComplexPerturbation &pert,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double factor = 0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(),
      quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    std::vector<std::vector<complex> > A_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));
    std::vector<double> grad_value(materials.get_n_groups());

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_A(dofs_per_cell, dofs_per_cell);

    std::vector<FullMatrix<double> > bound(materials.get_n_groups());
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    bool show_materials = true;
    bool show_boundary = false;
    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        const unsigned int mat_id = materials.get_material_id<dim>(cell);

        cell_grad = 0;
        cell_val = 0;
        for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
          bound[g] = 0;

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
            for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
            {
              types::boundary_id boundary_id = cell->face(f)->boundary_id();
              AssertIndexRange(boundary_id, boundary_conditions.size());
              if (show_boundary)
              {
                std::cout << " Cell " << cell->user_index() + 1
                << " face "
                << f << " with mat "
                << materials.get_material_id<dim>(cell) + 1
                << " at boundary " << int(boundary_id)
                << " with bc conditions "
                << boundary_conditions[boundary_id] << std::endl;
              }

              if (boundary_conditions[boundary_id] > 1)
              {
                fe_face_values.reinit(cell, f);
                switch (boundary_conditions[boundary_id])
                {
                  case 2: // Vacuum BC
                    factor = 0.5;
                    break;
                  default: // Custom Albedo BC
                    factor = albedo_factors[(boundary_conditions[boundary_id] - 3)
                                            * materials.get_n_groups()
                                            + g];
                    break;
                }

                for (unsigned int q_p = 0; q_p < n_face_q_points; ++q_p)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      val = fe_face_values.shape_value(i, q_p)
                            * fe_face_values.shape_value(j, q_p)
                            * fe_face_values.JxW(q_p);
                      bound[g](i, j) += factor * val;
                    }
              }
            }
          }

        for (unsigned int g = 0; g < materials.get_n_groups(); g++)
          grad_value[g] = materials.get_diffusion_coefficient(g, mat_id);

        calculate_A_factor(mat_id, materials, pert, A_factor);

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        {

          cell_A.equ(grad_value[gi], cell_grad, A_factor[gi][gi].real(),
            cell_val, 1.0, bound[gi]);

          constraints.distribute_local_to_global(cell_A,
            local_dof_indices,
            *(this->matrix_blocks_real[gi][gi]));
        }
      }

    for (unsigned int bi = 0; bi < materials.get_n_groups(); ++bi)
      this->matrix_blocks_real[bi][bi]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const ComplexPerturbation &pert,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double factor = 0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(),
      quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    std::vector<std::vector<complex> > A_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));
    std::vector<double> grad_value(materials.get_n_groups());

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_A(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups());
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    bool show_materials = false;
    bool show_boundary = false;
    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        const unsigned int mat_id = materials.get_material_id<dim>(cell);

        cell_grad = 0;
        cell_val = 0;
        for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
          bound[g] = 0;

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
            for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
            {
              types::boundary_id boundary_id = cell->face(f)->boundary_id();
              AssertIndexRange(boundary_id, boundary_conditions.size());
              if (show_boundary)
              {
                std::cout << " Cell " << cell->user_index() + 1
                << " face "
                << f << " with mat "
                << materials.get_material_id<dim>(cell) + 1
                << " at boundary " << int(boundary_id)
                << " with bc conditions "
                << boundary_conditions[boundary_id] << std::endl;
              }

              if (boundary_conditions[boundary_id] > 1)
              {
                fe_face_values.reinit(cell, f);
                switch (boundary_conditions[boundary_id])
                {
                  case 2: // Vacuum BC
                    factor = 0.5;
                    break;
                  default: // Custom Albedo BC
                    factor = albedo_factors[(boundary_conditions[boundary_id] - 3)
                                            * materials.get_n_groups()
                                            + g];
                    break;
                }

                for (unsigned int q_p = 0; q_p < n_face_q_points; ++q_p)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      val = fe_face_values.shape_value(i, q_p)
                            * fe_face_values.shape_value(j, q_p)
                            * fe_face_values.JxW(q_p);
                      bound[g](i, j) += factor * val;
                    }
              }
            }
          }

        for (unsigned int g = 0; g < materials.get_n_groups(); g++)
          grad_value[g] = materials.get_diffusion_coefficient(g, mat_id);

        calculate_A_factor(mat_id, materials, pert, A_factor);

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
          for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          {

            if (gi == gj)
            {
              cell_A.equ(grad_value[gi], cell_grad, A_factor[gi][gj].real(),
                cell_val, 1.0, bound[gi]);
              constraints.distribute_local_to_global(cell_A,
                local_dof_indices,
                *(this->matrix_blocks_real[gi][gj]));

            }
            else
            {
              cell_A.equ(A_factor[gi][gj].real(), cell_val);
              constraints.distribute_local_to_global(cell_A,
                local_dof_indices,
                *(this->matrix_blocks_real[gi][gj]));
            }

            cell_A.equ(A_factor[gi][gj].imag(), cell_val);
            constraints.distribute_local_to_global(cell_A,
              local_dof_indices,
              *(this->matrix_blocks_imag[gi][gj]));

          }

      }

    for (unsigned int bi = 0; bi < materials.get_n_groups(); ++bi)
      for (unsigned int bj = 0; bj < materials.get_n_groups(); ++bj)
      {
        this->matrix_blocks_real[bi][bj]->compress(VectorOperation::add);
        this->matrix_blocks_imag[bi][bj]->compress(VectorOperation::add);
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t NoiseAMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;
    if (this->matrixfree_type == non_diagonal)
    {
      memory = matfree_data.memory_consumption();
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks_real; ++i)
      {
        memory += this->matrix_blocks_real[i][i]->memory_consumption();
      }
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      memory = this->matfree_data.memory_consumption();
    }
    else // full_allocated
    {
      memory += this->sp.memory_consumption();

      for (unsigned int i = 0; i < this->n_blocks_real; ++i)
        for (unsigned int j = 0; j < this->n_blocks_real; ++j)
        {
          memory += this->matrix_blocks_real[i][j]->memory_consumption();
          memory += this->matrix_blocks_imag[i][j]->memory_consumption();
        }

    }
    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseAMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    // coeffs.clear();

    this->TransportMatrixComplexBase<dim, n_fe_degree>::clear();

    return;
  }

///**
// *
// */
//template <int dim, int n_fe_degree>
//  void NoiseAMatrix<dim, n_fe_degree>::clear ()
//  {
//  if (this->matrixfree_type == non_diagonal)
//  {
//    for (unsigned int row = 0; row < this->n_blocks; ++row)
//      for (unsigned int col = 0; col < this->n_blocks; ++col)
//        if (row != col)
//        {
//          this->mass_mf_blocks[row][col]->clear();
//          delete this->mass_mf_blocks[row][col];
//        }
//        else
//        {
//          this->matrix_blocks[row][col]->clear();
//          delete this->matrix_blocks[row][col];
//        }
//  }
//  else if (this->matrixfree_type == full_matrixfree)
//  {
//    for (unsigned int row = 0; row < this->n_blocks; ++row)
//      for (unsigned int col = 0; col < this->n_blocks; ++col)
//        if (row != col)
//        {
//          this->mass_mf_blocks[row][col]->clear();
//          delete this->mass_mf_blocks[row][col];
//        }
//        else
//        {
//          this->poison_mf_blocks[row]->clear();
//          delete this->poison_mf_blocks[row];
//        }
//  }
//  else
//    for (unsigned int row = 0; row < this->n_blocks; ++row)
//      for (unsigned int col = 0; col < this->n_blocks; ++col)
//        this->matrix_blocks[row][col]->clear();
//  }

//-------------------------------------------------------------------------------------------//
//                                    NoiseBMatrix
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  NoiseBMatrix<dim, n_fe_degree>::NoiseBMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      MassMatrixComplexBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseBMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const ComplexPerturbation &pert,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = 2 * n_groups;
    this->n_blocks_real = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == full_matrixfree)
    {
      AssertRelease(pert.pert_type == "Cell_Wise",
        "full_matrixfree in B Not implemented for Borders perturbations");
      const unsigned int n_mats = materials.get_n_mats();

      this->mass_mf_blocks_real.resize(this->n_blocks_real,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks_real));
      this->mass_mf_blocks_imag.resize(this->n_blocks_real,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks_real));

      coeffs_real.resize(n_groups, std::vector<Vector<double> >(n_groups));
      coeffs_imag.resize(n_groups, std::vector<Vector<double> >(n_groups));
      std::vector<std::vector<complex>> B_factor(materials.get_n_groups(),
        std::vector<complex>(materials.get_n_groups()));

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          coeffs_real[g1][g2].reinit(n_mats);
          coeffs_imag[g1][g2].reinit(n_mats);
        }

      // Fill coeffs
      for (unsigned int mat = 0; mat < n_mats; mat++)
      {
        calculate_B_factor(mat, mat, materials, pert, B_factor);

        // Calculate coeffs
        for (unsigned int g1 = 0; g1 < n_groups; g1++)
          for (unsigned int g2 = 0; g2 < n_groups; g2++)
          {
            coeffs_real[g1][g2][mat] = B_factor[g1][g2].real();
            coeffs_imag[g1][g2][mat] = B_factor[g1][g2].imag();
          }
      }
      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values
                                              | update_JxW_values);
      matfree_data.reinit(dof_handler, constraints,
        QGauss<1>(n_fe_degree + 1),
        additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks_real[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
            matfree_data);
          this->mass_mf_blocks_real[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs_real[gi][gj],
            materials.listen_to_material_id);

          this->mass_mf_blocks_imag[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
            matfree_data);
          this->mass_mf_blocks_imag[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs_imag[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {

      // Resize matrix_blocks
      this->matrix_blocks_real.resize(this->n_blocks_real,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks_real));
      this->matrix_blocks_imag.resize(this->n_blocks_real,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks_real));

      DoFTools::extract_locally_relevant_dofs(dof_handler,
        this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
      this->local_dofs_per_process =
                                     dof_handler.n_locally_owned_dofs_per_processor();

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->local_dofs_per_process,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int b1 = 0; b1 < this->n_blocks_real; b1++)
        for (unsigned int b2 = 0; b2 < this->n_blocks_real; b2++)
        {
          this->matrix_blocks_real[b1][b2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks_real[b1][b2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);

          this->matrix_blocks_imag[b1][b2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks_imag[b1][b2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }

      if (pert.pert_type == "Cell_Wise")
      {
        assemble_full_rhs_cell_wise(materials, pert);
      }
      else if (pert.pert_type == "Borders")
      {
        assemble_full_rhs_borders(materials, pert);
      }
      else if (pert.pert_type == "BordersHex")
      {
        assemble_full_rhs_bordershex(materials, pert);

      }
      else
        AssertRelease(false,
          "Not valid pert_type, only valid:  Cell_Wise | Borders | BordersHex");

    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseBMatrix<dim, n_fe_degree>::calculate_B_factor (
    const unsigned int mat_id,
    const unsigned int pert_id,
    const Materials &materials,
    const ComplexPerturbation &pert,
    std::vector<std::vector<complex> > &B_factor)
  {
    //B_factor(materials.get_n_groups(),
    //        std::vector<complex>(materials.get_n_groups()));

    std::vector<complex> prec_factor(materials.get_n_groups());
    const double omega = pert.get_frequency() * 2 * M_PI;

    calculate_precursors_factor(mat_id, materials, omega, prec_factor);

    for (unsigned int g1 = 0; g1 < materials.get_n_groups(); g1++)
      for (unsigned int g2 = 0; g2 < materials.get_n_groups(); g2++)
      {
        if (g1 == g2)
        {
          B_factor[g1][g1] = -pert.get_delta_sigma_r(g1, pert_id);
        }
        else
        {
          B_factor[g1][g2] = +pert.get_delta_sigma_s(g2, g1, pert_id);
        }

        B_factor[g1][g2] += prec_factor[g1]
                            * pert.get_delta_sigma_f(g2, pert_id)
                            / materials.keff;
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseBMatrix<dim, n_fe_degree>::assemble_full_rhs_cell_wise (
    const Materials &materials,
    const ComplexPerturbation &pert)
  {

    double val;
    unsigned int mat_id;

    std::vector<std::vector<complex> > B_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));

    std::vector<types::global_dof_index> local_dof_indices(
      dof_handler.get_fe().dofs_per_cell);
    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_B(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        mat_id = materials.get_material_id<dim>(cell);

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

        calculate_B_factor(mat_id, mat_id, materials, pert, B_factor);

        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
          for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          {
            cell_B.equ(B_factor[gi][gj].real(), cell_val); // Real
            constraints.distribute_local_to_global(cell_B,
              local_dof_indices,
              *(this->matrix_blocks_real[gi][gj]));

            cell_B.equ(B_factor[gi][gj].imag(), cell_val); // imag
            constraints.distribute_local_to_global(cell_B,
              local_dof_indices,
              *(this->matrix_blocks_imag[gi][gj]));
          }
      }

    for (unsigned int bi = 0; bi < this->n_blocks_real; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks_real; ++bj)
      {
        this->matrix_blocks_real[bi][bj]->compress(VectorOperation::add);
        this->matrix_blocks_imag[bi][bj]->compress(VectorOperation::add);
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseBMatrix<dim, n_fe_degree>::assemble_full_rhs_borders (
    const Materials &materials,
    const ComplexPerturbation &pert)
  {

    double val;
    unsigned int mat_id;

    std::vector<std::vector<complex>> B_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));

    std::vector<types::global_dof_index> local_dof_indices(
      dof_handler.get_fe().dofs_per_cell);

    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values);
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_face_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_B(dofs_per_cell, dofs_per_cell);

    for (auto cell : dof_handler.cell_iterators_on_level(0))
    {
      unsigned int n_active_children = cell->number_of_children();
      auto active_cells = GridTools::get_active_child_cells<DoFHandler<dim> >(cell);
      if (active_cells.size() == 0) /* get_active_child_cells only selects children */
      {
        active_cells.push_back(cell);
      }

      // -------------------------------------------------------------------------- //
      for (unsigned int ch = 0; ch < n_active_children; ch++)
      {
        //  if (cell_child->active())
        typename DoFHandler<dim>::cell_iterator cell_child_active = active_cells[ch];
        mat_id = materials.get_material_id<dim>(cell_child_active);
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
            ++face)
        {
          // We could have used child_cell_on_face()
          const unsigned int face_pert_id = pert.get_pertubation_face_id(
            mat_id, face, ch, n_active_children);
          if (face_pert_id != static_cast<unsigned int>(-1))
          {
            fe_face_values.reinit(cell_child_active, face);

            cell_B = 0;

            calculate_B_factor(mat_id, face_pert_id, materials, pert, B_factor);

            cell_face_val = 0;
            cell_child_active->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
                {
                  val = fe_face_values.shape_value(i, q_p)
                        * fe_face_values.shape_value(j, q_p)
                        * fe_face_values.JxW(q_p);

                  cell_face_val(i, j) += val;
                }

            for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
              for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
              {
                cell_B.equ(B_factor[gi][gj].real(), cell_face_val);
                constraints.distribute_local_to_global(cell_B,
                  local_dof_indices,
                  *(this->matrix_blocks_real[gi][gj]));

                cell_B.equ(B_factor[gi][gj].imag(), cell_face_val);
                constraints.distribute_local_to_global(cell_B,
                  local_dof_indices,
                  *(this->matrix_blocks_imag[gi][gj]));
              }
          } // End IF there is perturbation
        } // End FOR face
      } // End FOR active children
    } // End FOR cell at level(0)

    for (unsigned int bi = 0; bi < this->n_blocks_real; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks_real; ++bj)
      {
        this->matrix_blocks_real[bi][bj]->compress(VectorOperation::add);
        this->matrix_blocks_imag[bi][bj]->compress(VectorOperation::add);
      }

  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseBMatrix<dim, n_fe_degree>::assemble_full_rhs_bordershex (
    const Materials &materials,
    const ComplexPerturbation &pert)
  {

    double val;
    unsigned int mat_id;

    std::vector<unsigned int> quad_in_hex(materials.get_n_assemblies(), 0);

    std::vector<std::vector<complex>> B_factor(materials.get_n_groups(),
      std::vector<complex>(materials.get_n_groups()));

    std::vector<types::global_dof_index> local_dof_indices(
      dof_handler.get_fe().dofs_per_cell);

    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values);
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_face_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_B(dofs_per_cell, dofs_per_cell);

    // Integrate Interiors Cells
    // -------------------------------------------------------------------------- //
    for (auto cell : dof_handler.cell_iterators_on_level(0))
    {
      unsigned int n_active_children = cell->number_of_children();
      auto active_cells = GridTools::get_active_child_cells<DoFHandler<dim> >(cell);
      if (active_cells.size() == 0) /* get_active_child_cells only selects children */
      {
        active_cells.push_back(cell);
      }

      // -------------------------------------------------------------------------- //
      for (unsigned int ch = 0; ch < n_active_children; ch++)
      {
        typename DoFHandler<dim>::cell_iterator cell_child_active = active_cells[ch];
        const unsigned quad_id = quad_in_hex[cell_child_active->user_index()]
                                 / n_active_children;
        quad_in_hex[cell_child_active->user_index()]++;
        mat_id = materials.get_material_id<dim>(cell_child_active);
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
            ++face)
        {
          const unsigned int face_pert_id = pert.get_pertubation_face_id_hex(
            mat_id,
            quad_id, face,
            ch, n_active_children);

          if (face_pert_id != static_cast<unsigned int>(-1))
          {
            fe_face_values.reinit(cell_child_active, face);

            cell_B = 0;

            calculate_B_factor(mat_id, face_pert_id, materials, pert, B_factor);

            cell_face_val = 0;
            cell_child_active->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
                {
                  val = fe_face_values.shape_value(i, q_p)
                        * fe_face_values.shape_value(j, q_p)
                        * fe_face_values.JxW(q_p);

                  cell_face_val(i, j) += val;
                }

            for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
              for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
              {
                cell_B.equ(B_factor[gi][gj].real(), cell_face_val);
                constraints.distribute_local_to_global(cell_B,
                  local_dof_indices,
                  *(this->matrix_blocks_real[gi][gj]));

                cell_B.equ(B_factor[gi][gj].imag(), cell_face_val);
                constraints.distribute_local_to_global(cell_B,
                  local_dof_indices,
                  *(this->matrix_blocks_imag[gi][gj]));
              }

          } // End IF there is perturbation
        } // End FOR face
      } // End FOR active children
    } // End FOR cell at level(0)

    for (unsigned int bi = 0; bi < this->n_blocks_real; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks_real; ++bj)
      {
        this->matrix_blocks_real[bi][bj]->compress(VectorOperation::add);
        this->matrix_blocks_imag[bi][bj]->compress(VectorOperation::add);
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t NoiseBMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks_real; ++i)
        for (unsigned int j = 0; j < this->n_blocks_real; ++j)
        {
          memory += this->matrix_blocks_real[i][j]->memory_consumption();
          memory += this->matrix_blocks_imag[i][j]->memory_consumption();
        }
    }

    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void NoiseBMatrix<dim, n_fe_degree>::clear ()
  {
    matfree_data.clear();
    coeffs_real.clear();
    coeffs_imag.clear();

    this->MassMatrixComplexBase<dim, n_fe_degree>::clear();

    return;
  }

// ---------------------------------------------- //
// ----------- Explicit Instantations ----------- //
// ---------------------------------------------- //

template class NoiseAMatrix<1, 1> ;
template class NoiseAMatrix<1, 2> ;
template class NoiseAMatrix<1, 3> ;
template class NoiseAMatrix<1, 4> ;
template class NoiseAMatrix<1, 5> ;

template class NoiseAMatrix<2, 1> ;
template class NoiseAMatrix<2, 2> ;
template class NoiseAMatrix<2, 3> ;
template class NoiseAMatrix<2, 4> ;
template class NoiseAMatrix<2, 5> ;

template class NoiseAMatrix<3, 1> ;
template class NoiseAMatrix<3, 2> ;
template class NoiseAMatrix<3, 3> ;
template class NoiseAMatrix<3, 4> ;
template class NoiseAMatrix<3, 5> ;

template class NoiseBMatrix<1, 1> ;
template class NoiseBMatrix<1, 2> ;
template class NoiseBMatrix<1, 3> ;
template class NoiseBMatrix<1, 4> ;
template class NoiseBMatrix<1, 5> ;

template class NoiseBMatrix<2, 1> ;
template class NoiseBMatrix<2, 2> ;
template class NoiseBMatrix<2, 3> ;
template class NoiseBMatrix<2, 4> ;
template class NoiseBMatrix<2, 5> ;

template class NoiseBMatrix<3, 1> ;
template class NoiseBMatrix<3, 2> ;
template class NoiseBMatrix<3, 3> ;
template class NoiseBMatrix<3, 4> ;
template class NoiseBMatrix<3, 5> ;

