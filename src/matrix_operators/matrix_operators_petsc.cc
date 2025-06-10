/**
 * @file   matrix_operators_petsc.cc
 * @brief
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
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

#include "../../include/matrix_operators/matrix_operators_petsc.h"
#include "../../include/matrix_operators/matrix_operators_base.h"

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  TransportMatrix<dim, n_fe_degree>::TransportMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double> &constraints) :
      TransportMatrixBase<dim, n_fe_degree>(comm, dh, constraints),
      tria(dh.get_triangulation()),
      dof_handler(dh),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const MatrixFreeType &_matrixfree_type)
  {

    this->matrixfree_type = _matrixfree_type;
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->n_dofs_block = dof_handler.n_dofs();

    // Resize matrix_blocks
    this->matrix_blocks.resize(n_groups,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();

    if (this->matrixfree_type == non_diagonal)
    {
      reinit_non_diagonal(materials, boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      reinit_full_matrixfree(materials, boundary_conditions,
        albedo_factors);
    }
    else if (this->matrixfree_type == full_allocated)
    {

      DynamicSparsityPattern dsp(this->locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }

      assemble_full_matrices(materials, boundary_conditions, albedo_factors);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrix<dim, n_fe_degree>::reinit_non_diagonal (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
    coeffs.resize(n_groups,
      std::vector<Vector<double> >(n_groups));

    for (unsigned int from_g = 0; from_g < n_groups; from_g++)
      for (unsigned int to_g = 0; to_g < n_groups; to_g++)
      {
        coeffs[to_g][from_g].reinit(n_mats);
        if (to_g != from_g)
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = -materials.get_sigma_s(from_g, to_g, mat);
          }
      }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_JxW_values);
    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        if (gi != gj)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);

          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            this->coeffs[gi][gj],
            materials.listen_to_material_id);
        }

    //  --------- Matrices of the diagonal  ---------
    // Resize matrix_blocks
    this->matrix_blocks.resize(n_groups,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DynamicSparsityPattern dsp(this->locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

    SparsityTools::distribute_sparsity_pattern(dsp,
      this->locally_owned_dofs,
      this->comm,
      this->locally_relevant_dofs);
    this->sp.copy_from(dsp);

    for (unsigned int g1 = 0; g1 < n_groups; g1++)
    {
      this->matrix_blocks[g1][g1] = new PETScWrappers::MPI::SparseMatrix;
      this->matrix_blocks[g1][g1]->reinit(this->locally_owned_dofs,
        this->locally_owned_dofs,
        this->sp, this->comm);
    }

    assemble_diagonal_matrices(materials, boundary_conditions, albedo_factors);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrix<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
    this->poison_mf_blocks.resize(this->n_blocks);

    coeffs.resize(n_groups);
    //

    for (unsigned int to_g = 0; to_g < n_groups; to_g++)
    {
      coeffs[to_g].resize(n_groups);
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
      {
        coeffs[to_g][from_g].reinit(n_mats);
        if (to_g != from_g)
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = -materials.get_sigma_s(from_g,
              to_g,
              mat);
          }
      }
    }

    coeffs_grad.resize(n_groups);
    coeffs_val.resize(n_groups);

    for (unsigned int g = 0; g < n_groups; g++)
    {
      coeffs_grad[g].reinit(n_mats);
      coeffs_val[g].reinit(n_mats);
      for (unsigned int mat = 0; mat < n_mats; mat++)
      {
        coeffs_grad[g][mat] = materials.get_diffusion_coefficient(g, mat);
        coeffs_val[g][mat] = materials.get_sigma_r(g, mat);
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

    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        if (gi != gj)
        {
          this->mass_mf_blocks[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            this->coeffs[gi][gj],
            materials.listen_to_material_id);
        }
        else // Diagonal
        {
          this->poison_mf_blocks[gi] = new PoissonOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->poison_mf_blocks[gi]->reinit(gi,
            this->constraints,
            materials,
            materials.get_materials_vector(),
            this->coeffs_val[gi],
            this->coeffs_grad[gi],
            boundary_conditions,
            albedo_factors);
        }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrix<dim, n_fe_degree>::assemble_diagonal_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma_r;
    double factor = 0.0;
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

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups(),
      FullMatrix<double>(dofs_per_cell, dofs_per_cell));
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<bool> is_done(materials.n_total_assemblies, false);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        const unsigned int mat = materials.get_material_id<dim>(cell);

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

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        {
          // Get the material coefficients:
          D = materials.get_diffusion_coefficient(gi, mat);
          sigma_r = materials.get_sigma_r(gi, mat);
          cell_L.equ(D, cell_grad, sigma_r, cell_val, 1.0, bound[gi]);
          constraints.distribute_local_to_global(cell_L,
            local_dof_indices,
            *(this->matrix_blocks[gi][gi]));
        }
      }

    //
    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      this->matrix_blocks[gi][gi]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma_r, sigma_s;
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

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups());
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

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

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
          for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          {
            // Get the material coefficients:
            if (gi == gj)
            {
              D = materials.get_diffusion_coefficient(gi, mat);
              sigma_r = materials.get_sigma_r(gi, mat);
              cell_L.equ(D, cell_grad, sigma_r, cell_val, 1.0, bound[gi]);
              constraints.distribute_local_to_global(cell_L,
                local_dof_indices,
                *(this->matrix_blocks[gi][gj]));
            }
            else
            {
              sigma_s = -materials.get_sigma_s(gj, gi, mat);
              cell_L.equ(sigma_s, cell_val);
              constraints.distribute_local_to_global(cell_L, local_dof_indices,
                *(this->matrix_blocks[gi][gj]));
            }
          }
      }

    //
    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t TransportMatrix<dim, n_fe_degree>::memory_consumption () const
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

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->TransportMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

///**
// *
// */
//template <int dim, int n_fe_degree>
//  void TransportMatrix<dim, n_fe_degree>::clear ()
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

/**
 *
 */
template <int dim, int n_fe_degree>
  LeackageMatrix<dim, n_fe_degree>::LeackageMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double> &constraints) :
      TransportMatrixBase<dim, n_fe_degree>(comm, dh, constraints),
      tria(dh.get_triangulation()),
      dof_handler(dh),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void LeackageMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const MatrixFreeType &_matrixfree_type)
  {

    this->matrixfree_type = _matrixfree_type;
    if (this->matrixfree_type == non_diagonal)
      this->matrixfree_type = full_allocated;
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->n_dofs_block = dof_handler.n_dofs();

    // Resize matrix_blocks
    this->matrix_blocks.resize(n_groups,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();

    if (this->matrixfree_type == full_matrixfree)
    {
      reinit_full_matrixfree(materials, boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_allocated)
    {

      DynamicSparsityPattern dsp(this->locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }

      assemble_full_matrices(materials, boundary_conditions, albedo_factors);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void LeackageMatrix<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    AssertRelease(false, "this does not work yet");

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
    this->poison_mf_blocks.resize(this->n_blocks);
    coeffs.resize(n_groups,
      std::vector<Vector<double> >(n_groups,
        Vector<double>(n_mats)));
    //
    for (unsigned int from_g = 0; from_g < n_groups; from_g++)
      for (unsigned int to_g = 0; to_g < n_groups; to_g++)
      {
        if (to_g != from_g)
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = 0.0;
          }
      }

    coeffs_grad.resize(n_groups);
    coeffs_val.resize(n_groups);

    for (unsigned int g = 0; g < n_groups; g++)
    {
      coeffs_grad[g].reinit(n_mats);
      coeffs_val[g].reinit(n_mats);
      for (unsigned int mat = 0; mat < n_mats; mat++)
      {
        coeffs_grad[g][mat] = materials.get_diffusion_coefficient(g, mat);
        coeffs_val[g][mat] = materials.get_sigma_r(g, mat);
      }

    }
    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
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

    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        if (gi != gj)
        {
          this->mass_mf_blocks[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            this->coeffs[gi][gj],
            materials.listen_to_material_id);
        }
        else // Diagonal
        {
          this->poison_mf_blocks[gi] = new PoissonOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->poison_mf_blocks[gi]->reinit(gi,
            this->constraints,
            materials, materials.get_materials_vector(),
            this->coeffs_val[gi], this->coeffs_grad[gi],
            boundary_conditions,
            albedo_factors, false);
        }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void LeackageMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma_r;
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

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups());
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

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

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        {
          // Get the material coefficients:
          D = materials.get_diffusion_coefficient(gi, mat);
          sigma_r = materials.get_sigma_r(gi, mat);
          cell_L.equ(D, cell_grad, sigma_r, cell_val, 1.0, bound[gi]);
          constraints.distribute_local_to_global(cell_L,
            local_dof_indices,
            *(this->matrix_blocks[gi][gi]));

        }
      }

    //
    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      this->matrix_blocks[gi][gi]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t LeackageMatrix<dim, n_fe_degree>::memory_consumption () const
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

/**
 *
 */
template <int dim, int n_fe_degree>
  void LeackageMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->TransportMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  AlphaMatrix<dim, n_fe_degree>::AlphaMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double> &constraints) :
      TransportMatrixBase<dim, n_fe_degree>(comm, dh, constraints),
      tria(dh.get_triangulation()),
      dof_handler(dh),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void AlphaMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const MatrixFreeType &_matrixfree_type)
  {

    this->matrixfree_type = _matrixfree_type;
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->n_dofs_block = dof_handler.n_dofs();

    // Resize matrix_blocks
    this->matrix_blocks.resize(n_groups,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();

    if (this->matrixfree_type == non_diagonal)
    {
      reinit_non_diagonal(materials, boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      reinit_full_matrixfree(materials, boundary_conditions,
        albedo_factors);
    }
    else if (this->matrixfree_type == full_allocated)
    {

      DynamicSparsityPattern dsp(this->locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }

      assemble_full_matrices(materials, boundary_conditions, albedo_factors);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void AlphaMatrix<dim, n_fe_degree>::reinit_non_diagonal (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
    coeffs.resize(n_groups,
      std::vector<Vector<double> >(n_groups));

    for (unsigned int from_g = 0; from_g < n_groups; from_g++)
      for (unsigned int to_g = 0; to_g < n_groups; to_g++)
      {
        coeffs[to_g][from_g].reinit(n_mats);
        if (to_g != from_g)
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = materials.get_sigma_s(from_g, to_g, mat);
            coeffs[to_g][from_g][mat] += materials.get_xi_nu_sigma_f(from_g, to_g, mat);
          }
      }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_JxW_values);
    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        if (gi != gj)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);

          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            this->coeffs[gi][gj],
            materials.listen_to_material_id);
        }

    //  --------- Matrices of the diagonal  ---------
    // Resize matrix_blocks
    this->matrix_blocks.resize(n_groups,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

    DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DynamicSparsityPattern dsp(this->locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

    SparsityTools::distribute_sparsity_pattern(dsp,
      this->locally_owned_dofs,
      this->comm,
      this->locally_relevant_dofs);
    this->sp.copy_from(dsp);

    for (unsigned int g1 = 0; g1 < n_groups; g1++)
    {
      this->matrix_blocks[g1][g1] = new PETScWrappers::MPI::SparseMatrix;
      this->matrix_blocks[g1][g1]->reinit(this->locally_owned_dofs,
        this->locally_owned_dofs,
        this->sp, this->comm);
    }

    assemble_diagonal_matrices(materials, boundary_conditions, albedo_factors);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void AlphaMatrix<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    const unsigned int n_mats = materials.get_n_mats();
    const unsigned int n_groups = materials.get_n_groups();

    AssertRelease(false, "this does not work yet");

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
    this->poison_mf_blocks.resize(this->n_blocks);
    coeffs.resize(n_groups,
      std::vector<Vector<double> >(n_groups,
        Vector<double>(n_mats)));
    //
    for (unsigned int from_g = 0; from_g < n_groups; from_g++)
      for (unsigned int to_g = 0; to_g < n_groups; to_g++)
      {
        if (to_g != from_g)
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = materials.get_sigma_s(from_g, to_g, mat);
            coeffs[to_g][from_g][mat] += materials.get_xi_nu_sigma_f(from_g, to_g, mat);
          }
      }

    coeffs_grad.resize(n_groups);
    coeffs_val.resize(n_groups);

    for (unsigned int g = 0; g < n_groups; g++)
    {
      coeffs_grad[g].reinit(n_mats);
      coeffs_val[g].reinit(n_mats);
      for (unsigned int mat = 0; mat < n_mats; mat++)
      {
        coeffs_grad[g][mat] = materials.get_diffusion_coefficient(g, mat);
        coeffs_val[g][mat] = materials.get_sigma_r(g, mat);
      }
    }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
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

    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        if (gi != gj)
        {
          this->mass_mf_blocks[gi][gj] = new MassOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            this->coeffs[gi][gj],
            materials.listen_to_material_id);
        }
        else // Diagonal
        {
          this->poison_mf_blocks[gi] = new PoissonOperator<dim, n_fe_degree, double>(
            matfree_data);

          this->poison_mf_blocks[gi]->reinit(gi,
            this->constraints,
            materials, materials.get_materials_vector(),
            this->coeffs_val[gi], this->coeffs_grad[gi],
            boundary_conditions,
            albedo_factors, true);
        }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void AlphaMatrix<dim, n_fe_degree>::assemble_diagonal_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma_r, nu_sigma_f, coeff_val;
    double factor = 0.0;
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

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups(),
      FullMatrix<double>(dofs_per_cell, dofs_per_cell));
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        const unsigned int mat = materials.get_material_id<dim>(cell);

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

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        {
          // Get the material coefficients:
          D = -materials.get_diffusion_coefficient(gi, mat);
          sigma_r = -materials.get_sigma_r(gi, mat);
          nu_sigma_f = materials.get_xi_nu_sigma_f(gi, gi, mat);
          coeff_val = sigma_r + nu_sigma_f;
          cell_L.equ(D, cell_grad, coeff_val, cell_val, 1.0, bound[gi]);
          constraints.distribute_local_to_global(cell_L,
            local_dof_indices,
            *(this->matrix_blocks[gi][gi]));
        }
      }

    //
    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      this->matrix_blocks[gi][gi]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void AlphaMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma_r, sigma_s;
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

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups());
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

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

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
          for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          {
            // Get the material coefficients:
            if (gi == gj)
            {
              D = -materials.get_diffusion_coefficient(gi, mat);
              sigma_r = -materials.get_sigma_r(gi, mat);
              sigma_r += materials.get_xi_nu_sigma_f(gi, gi, mat);
              cell_L.equ(D, cell_grad, sigma_r, cell_val, 1.0, bound[gi]);
              constraints.distribute_local_to_global(cell_L,
                local_dof_indices,
                *(this->matrix_blocks[gi][gj]));
            }
            else
            {
              sigma_s = materials.get_sigma_s(gj, gi, mat);
              sigma_s += materials.get_xi_nu_sigma_f(gj, gi, mat);
              cell_L.equ(sigma_s, cell_val);
              constraints.distribute_local_to_global(cell_L, local_dof_indices,
                *(this->matrix_blocks[gi][gj]));
            }
          }
      }

    //
    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t AlphaMatrix<dim, n_fe_degree>::memory_consumption () const
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

/**
 *
 */
template <int dim, int n_fe_degree>
  void AlphaMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->TransportMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- FissionMatrix    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  FisionMatrix<dim, n_fe_degree>::FisionMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal or this->matrixfree_type
                                                 == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
      coeffs.resize(n_groups,
        std::vector<Vector<double> >(n_groups));
      // Fill coeffs
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
        {
          coeffs[to_g][from_g].reinit(n_mats);
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = materials.get_xi_nu_sigma_f(from_g, to_g, mat);
          }
        }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
        additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);
          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(n_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

      DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }
      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrix<dim, n_fe_degree>::assemble_full_matrices (
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
      if (cell->is_locally_owned())
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
          {
            // Get the material coefficients:
            xi_nu_sigma_f = materials.get_xi_nu_sigma_f(gj, gi, mat);
            cell_M.equ(xi_nu_sigma_f, cell_val);
            constraints.distribute_local_to_global(cell_M, local_dof_indices,
              *(this->matrix_blocks[gi][gj]));
          }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t FisionMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->FisionMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- FissionMatrix    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  FisionDelayedMatrix<dim, n_fe_degree>::FisionDelayedMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionDelayedMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal or this->matrixfree_type
                                                 == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
      coeffs.resize(n_groups);

      // Fill coeffs
      for (unsigned int to_g = 0; to_g < n_groups; to_g++)
      {
        coeffs[to_g].resize(n_groups);
        for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        {
          coeffs[to_g][from_g].reinit(n_mats);
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = materials.get_prompt_spectra(
                                          mat, to_g)
                                        * materials.get_nu_sigma_f(from_g, mat)
                                        * (1.0 - materials.get_delayed_fraction_sum(mat));
          }
        }
      }
      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
        additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);
          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(n_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

      DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }
      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionDelayedMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
    double coeffs_matrix;
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
      if (cell->is_locally_owned())
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
          {
            // Get the material coefficients:
            coeffs_matrix = materials.get_prompt_spectra(
                              mat, gi)
                            * materials.get_nu_sigma_f(gj, mat)
                            * (1.0 - materials.get_delayed_fraction_sum(mat));
            cell_M.equ(coeffs_matrix, cell_val);
            constraints.distribute_local_to_global(cell_M, local_dof_indices,
              *(this->matrix_blocks[gi][gj]));
          }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t FisionDelayedMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionDelayedMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->FisionMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- GammaMatrix    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  GammaMatrix<dim, n_fe_degree>::GammaMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void GammaMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal or this->matrixfree_type
                                                 == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
      coeffs.resize(n_groups,
        std::vector<Vector<double> >(n_groups));
      // Fill coeffs
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
        {
          coeffs[to_g][from_g].reinit(n_mats);
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = materials.get_xi_nu_sigma_f(from_g, to_g, mat);
            coeffs[to_g][from_g][mat] += materials.get_sigma_s(from_g, to_g, mat);
          }
        }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
        additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);
          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(n_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

      DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }
      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void GammaMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
    double xi_nu_sigma_f, sigma_s;
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
      if (cell->is_locally_owned())
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
          {
            // Get the material coefficients:
            xi_nu_sigma_f = materials.get_xi_nu_sigma_f(gj, gi, mat);
            sigma_s = materials.get_sigma_s(gj, gi, mat);
            double coefficient = xi_nu_sigma_f + sigma_s;
            cell_M.equ(coefficient, cell_val);
            constraints.distribute_local_to_global(cell_M, local_dof_indices,
              *(this->matrix_blocks[gi][gj]));
          }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t GammaMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void GammaMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->FisionMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- MassMatrix    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  MassMatrix<dim, n_fe_degree>::MassMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void MassMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal or this->matrixfree_type
                                                 == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
      coeffs.resize(n_groups,
        std::vector<Vector<double> >(n_groups));
      // Fill coeffs
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
        {
          coeffs[to_g][from_g].reinit(n_mats);
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            if (to_g == from_g)
              coeffs[to_g][from_g][mat] = 1.0;
          }
        }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
        additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);
          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(n_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

      DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }
      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);

    // Reinit de velocities
    this->velocities.resize(materials.get_n_mats());
    for (unsigned int mat = 0; mat < materials.get_n_mats(); mat++)
    {
      this->velocities[mat].resize(materials.get_n_groups());
      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      {
        this->velocities[mat][gi] = materials.get_velocity(mat, gi);
      }
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void MassMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
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
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

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
          {
            // Get the material coefficients:
            constraints.distribute_local_to_global(cell_val, local_dof_indices,
              *(this->matrix_blocks[gi][gj]));
          }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t MassMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void MassMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->FisionMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- MassMatrix    ------------------
//-------------------------------------------------------------------------------------------//

template <int dim, int n_fe_degree>
  VelocityMatrix<dim, n_fe_degree>::VelocityMatrix (
    const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void VelocityMatrix<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal or this->matrixfree_type
                                                 == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
      coeffs.resize(n_groups,
        std::vector<Vector<double> >(n_groups));
      // Fill coeffs
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
        {
          coeffs[to_g][from_g].reinit(n_mats);
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            if (to_g == from_g)
              coeffs[to_g][from_g][mat] = 1.0 / materials.get_velocity(mat, to_g);
          }
        }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
        additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks[gi][gj] =
                                         new MassOperator<dim, n_fe_degree, double>(
                                           matfree_data);
          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(),
            coeffs[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(n_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

      DoFTools::extract_locally_relevant_dofs(dof_handler, this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] = new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs,
            this->sp, this->comm);
        }
      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false, "Invalid matrixfree_type: " + this->matrixfree_type);

    // Reinit de velocities
    // Reinit de velocities
    this->velocities.resize(materials.get_n_mats());
    for (unsigned int mat = 0; mat < materials.get_n_mats(); mat++)
    {
      this->velocities[mat].resize(materials.get_n_groups());
      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      {
        this->velocities[mat][gi] = materials.get_velocity(mat, gi);
      }
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void VelocityMatrix<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
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
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        unsigned int mat = materials.get_material_id<dim>(cell);

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
          {
            if (gi == gj)
            {
              // Get the material coefficients:
              cell_M.equ(1.0 / materials.get_velocity(mat, gi), cell_val);
              constraints.distribute_local_to_global(cell_M, local_dof_indices,
                *(this->matrix_blocks[gi][gj]));
            }
          }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      this->matrix_blocks[gi][gi]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t VelocityMatrix<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void VelocityMatrix<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->FisionMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

//-------------------------------------------------------------------------------------------//
//  ------------------- SpectraBetaFission    ------------------
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints.
 */
template <int dim, int n_fe_degree>
  SpectraBetaFission<dim, n_fe_degree>::SpectraBetaFission (const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(dof_handler),
      constraints(constraints)
  {
    type_precursor = 0;
  }
/**
 *
 */
template <int dim, int n_fe_degree>
  void SpectraBetaFission<dim, n_fe_degree>::reinit (const Materials &materials,
    const MatrixFreeType &_matrixfree_type,
    const unsigned int _type_precursor)
  {
    type_precursor = _type_precursor;
    const unsigned int n_groups = materials.get_n_groups();
    this->n_blocks = n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(
          this->n_blocks));
      coeffs.resize(n_groups, std::vector<Vector<double> >(n_groups));
      // Fill coeffs
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
        {
          coeffs[to_g][from_g].reinit(n_mats);
          for (unsigned int mat = 0; mat < n_mats; mat++)
          {
            coeffs[to_g][from_g][mat] = materials.get_delayed_spectra(
                                          mat, type_precursor, to_g)
                                        * materials.get_delayed_fraction(mat,
                                          type_precursor)
                                        * materials.get_nu_sigma_f(from_g, mat);
          }
        }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values
                                              | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints,
        QGauss<1>(n_fe_degree + 1), additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
        {
          this->mass_mf_blocks[gi][gj] = new MassOperator<dim,
              n_fe_degree, double>(matfree_data);
          this->mass_mf_blocks[gi][gj]->reinit(constraints,
            materials.get_materials_vector(), coeffs[gi][gj],
            materials.listen_to_material_id);
        }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(n_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));

      DoFTools::extract_locally_relevant_dofs(dof_handler,
        this->locally_relevant_dofs);
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(this->locally_relevant_dofs);

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        this->locally_owned_dofs,
        this->comm,
        this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          this->matrix_blocks[g1][g2] =
                                        new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[g1][g2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs, this->sp, this->comm);
        }
      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SpectraBetaFission<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
    double xsec_value;
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
      if (cell->is_locally_owned())
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
          {
            // Get the material coefficients:
            xsec_value = materials.get_delayed_spectra(
                           mat, type_precursor, gi)
                         * materials.get_delayed_fraction(mat, type_precursor)
                         * materials.get_nu_sigma_f(gj, mat);
            cell_M.equ(xsec_value, cell_val);
            constraints.distribute_local_to_global(cell_M,
              local_dof_indices, *(this->matrix_blocks[gi][gj]));
          }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_blocks[gi][gj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SpectraBetaFission<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    this->FisionMatrixBase<dim, n_fe_degree>::clear();

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t SpectraBetaFission<dim, n_fe_degree>::memory_consumption () const
  {
    std::size_t memory = 0;

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
      memory = matfree_data.memory_consumption();
    else
    {
      memory += this->sp.memory_consumption();
      for (unsigned int i = 0; i < this->n_blocks; ++i)
        for (unsigned int j = 0; j < this->n_blocks; ++j)
          memory += this->matrix_blocks[i][j]->memory_consumption();
    }

    return memory;
  }

// ----------- Explicit Instantations ----------- //

template class TransportMatrix<1, 1> ;
template class TransportMatrix<1, 2> ;
template class TransportMatrix<1, 3> ;
template class TransportMatrix<1, 4> ;
template class TransportMatrix<1, 5> ;

template class TransportMatrix<2, 1> ;
template class TransportMatrix<2, 2> ;
template class TransportMatrix<2, 3> ;
template class TransportMatrix<2, 4> ;
template class TransportMatrix<2, 5> ;

template class TransportMatrix<3, 1> ;
template class TransportMatrix<3, 2> ;
template class TransportMatrix<3, 3> ;
template class TransportMatrix<3, 4> ;
template class TransportMatrix<3, 5> ;

template class LeackageMatrix<1, 1> ;
template class LeackageMatrix<1, 2> ;
template class LeackageMatrix<1, 3> ;
template class LeackageMatrix<1, 4> ;
template class LeackageMatrix<1, 5> ;

template class LeackageMatrix<2, 1> ;
template class LeackageMatrix<2, 2> ;
template class LeackageMatrix<2, 3> ;
template class LeackageMatrix<2, 4> ;
template class LeackageMatrix<2, 5> ;

template class LeackageMatrix<3, 1> ;
template class LeackageMatrix<3, 2> ;
template class LeackageMatrix<3, 3> ;
template class LeackageMatrix<3, 4> ;
template class LeackageMatrix<3, 5> ;

template class AlphaMatrix<1, 1> ;
template class AlphaMatrix<1, 2> ;
template class AlphaMatrix<1, 3> ;
template class AlphaMatrix<1, 4> ;
template class AlphaMatrix<1, 5> ;

template class AlphaMatrix<2, 1> ;
template class AlphaMatrix<2, 2> ;
template class AlphaMatrix<2, 3> ;
template class AlphaMatrix<2, 4> ;
template class AlphaMatrix<2, 5> ;

template class AlphaMatrix<3, 1> ;
template class AlphaMatrix<3, 2> ;
template class AlphaMatrix<3, 3> ;
template class AlphaMatrix<3, 4> ;
template class AlphaMatrix<3, 5> ;

template class FisionMatrix<1, 1> ;
template class FisionMatrix<1, 2> ;
template class FisionMatrix<1, 3> ;
template class FisionMatrix<1, 4> ;
template class FisionMatrix<1, 5> ;

template class FisionMatrix<2, 1> ;
template class FisionMatrix<2, 2> ;
template class FisionMatrix<2, 3> ;
template class FisionMatrix<2, 4> ;
template class FisionMatrix<2, 5> ;

template class FisionMatrix<3, 1> ;
template class FisionMatrix<3, 2> ;
template class FisionMatrix<3, 3> ;
template class FisionMatrix<3, 4> ;
template class FisionMatrix<3, 5> ;

template class FisionDelayedMatrix<1, 1> ;
template class FisionDelayedMatrix<1, 2> ;
template class FisionDelayedMatrix<1, 3> ;
template class FisionDelayedMatrix<1, 4> ;
template class FisionDelayedMatrix<1, 5> ;

template class FisionDelayedMatrix<2, 1> ;
template class FisionDelayedMatrix<2, 2> ;
template class FisionDelayedMatrix<2, 3> ;
template class FisionDelayedMatrix<2, 4> ;
template class FisionDelayedMatrix<2, 5> ;

template class FisionDelayedMatrix<3, 1> ;
template class FisionDelayedMatrix<3, 2> ;
template class FisionDelayedMatrix<3, 3> ;
template class FisionDelayedMatrix<3, 4> ;
template class FisionDelayedMatrix<3, 5> ;

template class GammaMatrix<1, 1> ;
template class GammaMatrix<1, 2> ;
template class GammaMatrix<1, 3> ;
template class GammaMatrix<1, 4> ;
template class GammaMatrix<1, 5> ;

template class GammaMatrix<2, 1> ;
template class GammaMatrix<2, 2> ;
template class GammaMatrix<2, 3> ;
template class GammaMatrix<2, 4> ;
template class GammaMatrix<2, 5> ;

template class GammaMatrix<3, 1> ;
template class GammaMatrix<3, 2> ;
template class GammaMatrix<3, 3> ;
template class GammaMatrix<3, 4> ;
template class GammaMatrix<3, 5> ;

template class MassMatrix<1, 1> ;
template class MassMatrix<1, 2> ;
template class MassMatrix<1, 3> ;
template class MassMatrix<1, 4> ;
template class MassMatrix<1, 5> ;

template class MassMatrix<2, 1> ;
template class MassMatrix<2, 2> ;
template class MassMatrix<2, 3> ;
template class MassMatrix<2, 4> ;
template class MassMatrix<2, 5> ;

template class MassMatrix<3, 1> ;
template class MassMatrix<3, 2> ;
template class MassMatrix<3, 3> ;
template class MassMatrix<3, 4> ;
template class MassMatrix<3, 5> ;

template class VelocityMatrix<1, 1> ;
template class VelocityMatrix<1, 2> ;
template class VelocityMatrix<1, 3> ;
template class VelocityMatrix<1, 4> ;
template class VelocityMatrix<1, 5> ;

template class VelocityMatrix<2, 1> ;
template class VelocityMatrix<2, 2> ;
template class VelocityMatrix<2, 3> ;
template class VelocityMatrix<2, 4> ;
template class VelocityMatrix<2, 5> ;

template class VelocityMatrix<3, 1> ;
template class VelocityMatrix<3, 2> ;
template class VelocityMatrix<3, 3> ;
template class VelocityMatrix<3, 4> ;
template class VelocityMatrix<3, 5> ;

template class SpectraBetaFission<1, 1> ;
template class SpectraBetaFission<1, 2> ;
template class SpectraBetaFission<1, 3> ;
template class SpectraBetaFission<1, 4> ;
template class SpectraBetaFission<1, 5> ;

template class SpectraBetaFission<2, 1> ;
template class SpectraBetaFission<2, 2> ;
template class SpectraBetaFission<2, 3> ;
template class SpectraBetaFission<2, 4> ;
template class SpectraBetaFission<2, 5> ;

template class SpectraBetaFission<3, 1> ;
template class SpectraBetaFission<3, 2> ;
template class SpectraBetaFission<3, 3> ;
template class SpectraBetaFission<3, 4> ;
template class SpectraBetaFission<3, 5> ;
