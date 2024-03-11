/**
 * @file   matrix_operators_free.cc
 * @brief  Implementation of TransportMatrix and FissionMAtrix classes to handle block matrices.
 */

#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <algorithm>
#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "../../include/matrix_operators/matrix_operators_free.h"
#include "../../include/materials.h"
#include "../../include/femffusion.h"

using namespace dealii;

// -----------------------------------------------------------------------------------------------//
// Matrix Free Version of the Source Operators

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  MassOperator<dim, n_fe_degree, number>::MassOperator (
    dealii::MatrixFree<dim, number> &_matfree_data) :
      matfree_data(_matfree_data)
  {
    // Dummy initialization
    materials = NULL;

    alpha = NULL;
    constraints = NULL;

    boundary_conditions = NULL;
    bc_coeff = 0.0;

    zero_matrix = false;
    listen_to_material_id = false;
    boundary_req = false;
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  unsigned int MassOperator<dim, n_fe_degree, number>::size () const
  {
    return matfree_data.get_vector_partitioner()->size();
  }

/**
 * m()
 * @brief Return the number of rows in this matrix.
 */
template <int dim, int n_fe_degree, typename number>
  unsigned int MassOperator<dim, n_fe_degree, number>::m () const
  {
    return matfree_data.get_vector_partitioner()->size();
  }

/**
 * @brief Return the number of columns in this matrix.
 */
template <int dim, int n_fe_degree, typename number>
  unsigned int MassOperator<dim, n_fe_degree, number>::n () const
  {
    return matfree_data.get_vector_partitioner()->size();
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::clear ()
  {
    constraints = NULL;

    boundary_req = false;

    this->~MassOperator();
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::reinit (
    const AffineConstraints<double> &_constraints,
    const std::vector<unsigned int> &_materials,
    const Vector<double> &_alpha,
    const bool _listen_to_material_id)
  {
    materials = &_materials;
    alpha = &_alpha;
    listen_to_material_id = _listen_to_material_id;

    // Easy access to problem structures
    constraints = &(_constraints);

    if (alpha->all_zero())
    {
      zero_matrix = true;
      return;
    }

  }

//
// reinit(
//
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::reinit (
    const AffineConstraints<double> &_contraints,
    const std::vector<unsigned int> &_materials,
    const Vector<double> &_alpha,
    const std::vector<unsigned int> &_boundary_conditions,
    const double &_bc_coeff,
    const bool _listen_to_material_id)

  {
    materials = &_materials;
    alpha = &_alpha;
    listen_to_material_id = _listen_to_material_id;

    // Easy access to problem structures
    constraints = &_contraints;

    // Boundary
    boundary_conditions = &_boundary_conditions;
    bc_coeff = _bc_coeff;

    for (unsigned int bc = 0; bc < boundary_conditions->size(); bc++)
    {
      if ((*boundary_conditions)[bc] > 1)
      {
        boundary_req = true;

        break;
      }
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::cell_local_apply (
    const dealii::MatrixFree<dim, number> &matfree_data,
    ParallelVector &dst,
    const ParallelVector &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    VectorizedArray<double> alpha_va;
    // Only first two template parameters always required
    // <dimension, n_fe_degree, quad_degree, n_components, number_type>
    FEEvaluation<dim, n_fe_degree, n_fe_degree + 1, 1, number> fe_eval(
      matfree_data);

    // Loop over the cells range we are working in this kernel
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      alpha_va = 0;

      for (unsigned int v = 0; v < matfree_data.n_active_entries_per_cell_batch(cell);
          ++v)
      {
        typename DoFHandler<dim>::cell_iterator cell_it =
                                                          matfree_data.get_cell_iterator(
                                                            cell, v);
        if (listen_to_material_id)
        {

          alpha_va[v] = (*alpha)[cell_it->material_id()];
        }
        else
        {
          alpha_va[v] = (*alpha)[(*materials)[cell_it->user_index()]];
        }
      }

      // Reinit Values
      fe_eval.reinit(cell);

      // Read in the values of the source vector including the resolution constraints.
      // This stores u_cell.
      fe_eval.read_dof_values(src);

      // Compute only unit cell values  (not gradients or hessians)
      fe_eval.evaluate(true, false, false);

      // Jacobi transformation:
      //  alpha * values(in real space) * JxW
      //AssertIndexRange(cell, materials->size());

      for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
        fe_eval.submit_value(alpha_va * fe_eval.get_value(q), q);

      // Summation over all quadrature points in the cell.
      // Only integrate the values, not the gradients
      fe_eval.integrate(true, false);

      // Distribute the local values in the global vector
      fe_eval.distribute_local_to_global(dst);
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::cell_boundary_apply (
    const dealii::MatrixFree<dim, number> &matfree_data,
    ParallelVector &dst,
    const ParallelVector &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {

    VectorizedArray<double> albedo_factor;
    // Only first two template parameters always required
    // <dimension, n_fe_degree, quad_degree, n_components, number_type>
    FEFaceEvaluation<dim, n_fe_degree, n_fe_degree + 1, 1, number> face_eval(
      matfree_data);

    // Loop over the cells range we are working in this kernel
    for (unsigned int face = face_range.first; face < face_range.second;
        ++face)
    {
      types::boundary_id boundary_id = matfree_data.get_boundary_id(face);
      switch ((*boundary_conditions)[boundary_id])
      {
        case 0: // Vacuum BC
          AssertRelease(false,
            "Zero Flux Bc not possible, fake it with albedo_factor=100.0");
          break;
        case 1: // Zero Current BC
          continue;
        case 2: // Vacuum BC
          albedo_factor = (bc_coeff);
          break;

        default: // Custom Albedo BC
          albedo_factor = (bc_coeff) * 100.0;
          break;
      }

      // Reinit Values
      face_eval.reinit(face);

      // Compute only unit cell values not gradients
      face_eval.gather_evaluate(src, true, false);

      // Jacobi transformation:
      //  alpha * values(in real space) * JxW
      //AssertIndexRange(cell, materials->size());

      for (unsigned int q = 0; q < face_eval.n_q_points; ++q)
        face_eval.submit_value(albedo_factor * face_eval.get_value(q), q);

      // Summation over all quadrature points in the cell.
      // Only integrate the values and gradients
      // and  Distribute the local values in the global vector
      face_eval.integrate_scatter(true, false, dst);

    }
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::cell_face_apply (
    const dealii::MatrixFree<dim, number>&,
    ParallelVector&,
    const ParallelVector&,
    const std::pair<unsigned int, unsigned int>&) const
  {
    return;
  }

// ------------------------------------------------------------------ //
// ---------------------- ParallelVector operations ----------------- //
// ------------------------------------------------------------------ //
/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::vmult (
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    dst = 0;

    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::vmult_add (
    ParallelVector &dst,
    const ParallelVector &src) const
  {

    if (zero_matrix) // If zero matrix do not add anything
      return;

    if (boundary_req)
    {
      matfree_data.loop(&MassOperator::cell_local_apply,
        &MassOperator::cell_face_apply,
        &MassOperator::cell_boundary_apply, this, dst, src,
        /*zero_dst =*/false,
        MatrixFree<dim, number>::DataAccessOnFaces::none,
        MatrixFree<dim, number>::DataAccessOnFaces::none);
    }
    else
      matfree_data.cell_loop(&MassOperator::cell_local_apply, this, dst, src);

  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::Tvmult (
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::Tvmult_add (
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    vmult_add(dst, src);
  }

// ------------------------------------------------------------------ //
// --------------- Now for PETScWrappers::MPI::Vector --------------- //
// ------------------------------------------------------------------ //
/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::vmult (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    dst = 0;

    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::vmult_add (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {

    if (zero_matrix) // If zero matrix do not add anything
      return;

    ParallelVector dst_vec;
    ParallelVector src_vec;
    matfree_data.initialize_dof_vector(dst_vec);
    matfree_data.initialize_dof_vector(src_vec);

    LinearAlgebra::ReadWriteVector<double> rwv(src_vec.locally_owned_elements());
    rwv.import(src, VectorOperation::insert);
    src_vec.import(rwv, VectorOperation::insert);

    rwv.import(dst, VectorOperation::insert);
    dst_vec.import(rwv, VectorOperation::insert);

    if (boundary_req)
    {
      matfree_data.loop(&MassOperator::cell_local_apply,
        &MassOperator::cell_face_apply,
        &MassOperator::cell_boundary_apply, this, dst_vec, src_vec,
        /*zero_dst =*/false,
        MatrixFree<dim, number>::DataAccessOnFaces::values,
        MatrixFree<dim, number>::DataAccessOnFaces::values);
    }
    else
    {
      matfree_data.cell_loop(&MassOperator::cell_local_apply, this, dst_vec,
        src_vec);
    }
    copy_to_Vector(dst, dst_vec);

  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::Tvmult (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::Tvmult_add (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    vmult_add(dst, src);
  }

// ------------------------------------------------------------------ //
// ----------------------- Now for (PETSc) Vec ---------------------- //
// ------------------------------------------------------------------ //
/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::vmult (
    Vec &dst,
    const Vec &src) const
  {
    dst = 0;

    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::vmult_add (
    Vec &dst,
    const Vec &src) const
  {

    if (zero_matrix) // If zero matrix do not add anything
      return;
    ParallelVector src_vec, dst_vec;
    matfree_data.initialize_dof_vector(src_vec);
    matfree_data.initialize_dof_vector(dst_vec);

    copy_to_Vector(src_vec, src);
    copy_to_Vector(dst_vec, dst);

    if (boundary_req)
    {
      matfree_data.loop(&MassOperator::cell_local_apply,
        &MassOperator::cell_face_apply,
        &MassOperator::cell_boundary_apply, this, dst_vec, src_vec,
        /*zero_dst =*/false,
        MatrixFree<dim, number>::DataAccessOnFaces::none,
        MatrixFree<dim, number>::DataAccessOnFaces::none);
    }
    else
    {
      matfree_data.cell_loop(&MassOperator::cell_local_apply, this, dst_vec,
        src_vec);
    }
    copy_to_Vec(dst, dst_vec);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::Tvmult (Vec &dst,
    const Vec &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void MassOperator<dim, n_fe_degree, number>::Tvmult_add (Vec &dst,
    const Vec &src) const
  {
    vmult_add(dst, src);
  }

// ----------- Explicit Instantations ----------- //
template class MassOperator<1, 1, double> ;
template class MassOperator<1, 2, double> ;
template class MassOperator<1, 3, double> ;
template class MassOperator<1, 4, double> ;
template class MassOperator<1, 5, double> ;

template class MassOperator<2, 1, double> ;
template class MassOperator<2, 2, double> ;
template class MassOperator<2, 3, double> ;
template class MassOperator<2, 4, double> ;
template class MassOperator<2, 5, double> ;

template class MassOperator<3, 1, double> ;
template class MassOperator<3, 2, double> ;
template class MassOperator<3, 3, double> ;
template class MassOperator<3, 4, double> ;
template class MassOperator<3, 5, double> ;

// -----------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------//

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  PoissonOperator<dim, n_fe_degree, number>::PoissonOperator (
    dealii::MatrixFree<dim, number> &_matfree_data) :
      MatrixFreeOperators::Base<dim, ParallelVector>(),
      matfree_data(
        _matfree_data)
  {
    // Dummy initialization
    boundary_req = false;
    group = -1;

    constraints = NULL;

    materials = NULL;

    materials_vector = NULL;
    coef_val = NULL;
    coef_grad = NULL;

    boundary_conditions = NULL;
    albedo_factors = NULL;

    alpha_modes = false;
    bc_factor = 0.5;

    inverse_diagonal = 0;
  }

/**
 *
 */

template <int dim, int n_fe_degree, typename number>
  unsigned int PoissonOperator<dim, n_fe_degree, number>::size () const
  {
    return matfree_data.get_vector_partitioner()->size();
  }

/**
 * m()
 * @brief Return the number of rows in this matrix.
 */
template <int dim, int n_fe_degree, typename number>
  unsigned int PoissonOperator<dim, n_fe_degree, number>::m () const
  {
    return matfree_data.get_vector_partitioner()->size();
  }

/**
 * @brief Return the number of columns in this matrix.
 */
template <int dim, int n_fe_degree, typename number>
  unsigned int PoissonOperator<dim, n_fe_degree, number>::n () const
  {
    return matfree_data.get_vector_partitioner()->size();
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::clear ()
  {
    constraints = NULL;

    this->~PoissonOperator();
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::reinit (
    const unsigned int _group,
    const AffineConstraints<double> &_contraints,
    const Materials &_materials,
    const std::vector<unsigned int> &_materials_vec,
    const Vector<double> &_coef_val,
    const Vector<double> &_coef_grad,
    const std::vector<unsigned int> &_boundary_conditions,
    const std::vector<double> &_albedo_factors,
    const bool _alpha_modes,
    const double _bc_factor)

  {
    // Easy access to problem structures
    group = _group;

    constraints = &_contraints;

    materials = &_materials;
    materials_vector = &_materials_vec;
    coef_val = &_coef_val;
    coef_grad = &_coef_grad;
    boundary_conditions = &_boundary_conditions;
    albedo_factors = &_albedo_factors;
    bc_factor = _bc_factor;

    for (unsigned int bc = 0; bc < boundary_conditions->size(); bc++)
      if ((*boundary_conditions)[bc] > 1)
      {
        boundary_req = true;
        break;
      }

    alpha_modes = _alpha_modes;
  }

template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::compute_diagonal ()
  {
    matfree_data.initialize_dof_vector(inverse_diagonal);

    unsigned int dummy = 0;
    this->matfree_data.cell_loop(&PoissonOperator::local_compute_diagonal, this,
      inverse_diagonal, dummy);

    IndexSet index_set(inverse_diagonal.locally_owned_elements());
    for (IndexSet::ElementIterator it = index_set.begin();
        it != index_set.end(); it++)
    {
      if (inverse_diagonal[*it] > 0 or inverse_diagonal[*it] < 0)
        inverse_diagonal[*it] = 1.0 / inverse_diagonal[*it];

    }

    inverse_diagonal.compress(VectorOperation::insert);
  }

template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::local_compute_diagonal (
    const MatrixFree<dim, number> &data,
    ParallelVector &dst,
    const unsigned int&,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, n_fe_degree, n_fe_degree + 1, 1, number> fe_eval(data);
    AlignedVector<VectorizedArray<number> > diagonal(fe_eval.dofs_per_cell);
    VectorizedArray<double> sigma_va, D_va;
    for (unsigned int cell = cell_range.first; cell < cell_range.second;
        ++cell)
    {
      sigma_va = 0;
      D_va = 0;
      for (unsigned int v = 0; v < matfree_data.n_active_entries_per_cell_batch(cell);
          ++v)
      {
        typename DoFHandler<dim>::active_cell_iterator cell_it =
            matfree_data.get_cell_iterator(cell, v);

        D_va[v] = (*coef_grad)[(*materials_vector)[cell_it->user_index()]];
        sigma_va[v] =
                      (*coef_val)[(*materials_vector)[cell_it->user_index()]];

      }

      fe_eval.reinit(cell);
      for (unsigned int i = 0; i < fe_eval.dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < fe_eval.dofs_per_cell; ++j)
          fe_eval.submit_dof_value(VectorizedArray<number>(), j);
        fe_eval.submit_dof_value(make_vectorized_array<number>(1.), i);

        fe_eval.evaluate(true, true);
        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
        {
          fe_eval.submit_value(sigma_va * fe_eval.get_value(q), q);
          fe_eval.submit_gradient(D_va * fe_eval.get_gradient(q), q);
        }
        fe_eval.integrate(true, true);
        diagonal[i] = fe_eval.get_dof_value(i);
      }
      for (unsigned int i = 0; i < fe_eval.dofs_per_cell; ++i)
        fe_eval.submit_dof_value(diagonal[i], i);
      fe_eval.distribute_local_to_global(dst);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::cell_local_apply (
    const dealii::MatrixFree<dim, number> &matfree_data,
    ParallelVector &dst,
    const ParallelVector &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    VectorizedArray<double> cell_va, gradient_va;
    // Only first two template parameters always required
    // <dimension, n_fe_degree, quad_degree, n_components, number_type>
    FEEvaluation<dim, n_fe_degree> fe_eval(matfree_data);

    // Loop over the cells range we are working in this kernel
    for (unsigned int cell = cell_range.first; cell < cell_range.second;
        ++cell)
    {
      cell_va = 0;
      gradient_va = 0;
      for (unsigned int v = 0; v < matfree_data.n_active_entries_per_cell_batch(cell);
          ++v)
      {
        typename DoFHandler<dim>::active_cell_iterator cell_it =
            matfree_data.get_cell_iterator(cell, v);

        gradient_va[v] = (*coef_grad)[(*materials_vector)[cell_it->user_index()]];
        cell_va[v] = (*coef_val)[(*materials_vector)[cell_it->user_index()]];

      }

      // Reinit Values
      fe_eval.reinit(cell);

      // Read in the values of the source vector including the resolution constraints.
      // This stores u_cell.
      fe_eval.read_dof_values(src);

      // Compute only unit cell values and gradients (not hessians)
      fe_eval.evaluate(true, true, false);

      // Jacobi transformation:
      //  alpha * values(in real space) * JxW
      //AssertIndexRange(cell, materials->size());

      for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
        fe_eval.submit_value(cell_va * fe_eval.get_value(q), q);

      for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
        fe_eval.submit_gradient(gradient_va * fe_eval.get_gradient(q), q);

      // Summation over all quadrature points in the cell.
      // Only integrate the values and gradients
      fe_eval.integrate(true, true);

      // Distribute the local values in the global vector
      fe_eval.distribute_local_to_global(dst);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::cell_boundary_apply (
    const dealii::MatrixFree<dim, number> &matfree_data,
    ParallelVector &dst,
    const ParallelVector &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    VectorizedArray<double> albedo_factor;
    // Only first two template parameters always required
    // <dimension, n_fe_degree, quad_degree, n_components, number_type>
    FEFaceEvaluation<dim, n_fe_degree, n_fe_degree + 1, 1, number> face_eval(
      matfree_data);

    // Loop over the cells range we are working in this kernel
    for (unsigned int face = face_range.first; face < face_range.second;
        ++face)
    {
      types::boundary_id boundary_id = matfree_data.get_boundary_id(face);
      switch ((*boundary_conditions)[boundary_id])
      {
        case 0: // Vacuum BC
          AssertRelease(false,
            "Zero Flux Bc not possible, fake it with albedo_factor=100.0");
          break;
        case 1: // Zero Current BC
          continue;
        case 2: // Vacuum BC
          albedo_factor = bc_factor;
          break;

        default: // Custom Albedo BC
          albedo_factor =
                          (*albedo_factors)[((*boundary_conditions)[boundary_id] - 3)
                                            * materials->get_n_groups()
                                            + group];
          break;
      }

      // Reinit Values
      face_eval.reinit(face);

      // Compute only unit cell values not gradients
      face_eval.gather_evaluate(src, true, false);

      // Jacobi transformation:
      //  alpha * values(in real space) * JxW
      //AssertIndexRange(cell, materials->size());

      for (unsigned int q = 0; q < face_eval.n_q_points; ++q)
        face_eval.submit_value(albedo_factor * face_eval.get_value(q), q);

      // Summation over all quadrature points in the cell.
      // Only integrate the values and gradients
      // and  Distribute the local values in the global vector
      face_eval.integrate_scatter(true, false, dst);

    }
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::cell_face_apply (
    const dealii::MatrixFree<dim, number>&,
    ParallelVector&,
    const ParallelVector&,
    const std::pair<unsigned int, unsigned int>&) const
  {
    return;
  }

// ------------------------------------------------------------------ //
// ---------------------- ParallelVector operations ----------------- //
// ------------------------------------------------------------------ //
/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::apply_add (ParallelVector &dst,
    const ParallelVector &src) const
  {
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::vmult (ParallelVector &dst,
    const ParallelVector &src) const
  {
    dst = 0;

    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::vmult_add (ParallelVector &dst,
    const ParallelVector &src) const
  {

    if (boundary_req)
    {
      matfree_data.loop(&PoissonOperator::cell_local_apply,
        &PoissonOperator::cell_face_apply,
        &PoissonOperator::cell_boundary_apply, this, dst, src,
        /*zero_dst =*/false,
        MatrixFree<dim, number>::DataAccessOnFaces::none,
        MatrixFree<dim, number>::DataAccessOnFaces::none);
    }
    else
      matfree_data.cell_loop(&PoissonOperator::cell_local_apply, this, dst,
        src);

  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::Tvmult (
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::Tvmult_add (
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    vmult_add(dst, src);
  }

// ------------------------------------------------------------------ //
// --------------- Now for PETScWrappers::MPI::Vector --------------- //
// ------------------------------------------------------------------ //

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::vmult (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    dst = 0;

    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::vmult_add (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {

    ParallelVector dst_vec;
    ParallelVector src_vec;
    matfree_data.initialize_dof_vector(dst_vec);
    matfree_data.initialize_dof_vector(src_vec);

    // src_vec = src
    LinearAlgebra::ReadWriteVector<double> rwv(src_vec.locally_owned_elements());
    rwv.import(src, VectorOperation::insert);
    src_vec.import(rwv, VectorOperation::insert);
    // dst_vec = dst
    rwv.import(dst, VectorOperation::insert);
    dst_vec.import(rwv, VectorOperation::insert);

    if (boundary_req)
    {
      matfree_data.loop(&PoissonOperator::cell_local_apply,
        &PoissonOperator::cell_face_apply,
        &PoissonOperator::cell_boundary_apply, this, dst_vec, src_vec,
        /*zero_dst =*/false,
        MatrixFree<dim, number>::DataAccessOnFaces::none,
        MatrixFree<dim, number>::DataAccessOnFaces::none);
    }
    else
    {
      matfree_data.cell_loop(&PoissonOperator::cell_local_apply, this,
        dst_vec, src_vec);
    }
    copy_to_Vector(dst, dst_vec);

  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::Tvmult (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::Tvmult_add (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    vmult_add(dst, src);
  }

// ------------------------------------------------------------------ //
// ----------------------- Now for (PETSc) Vec ---------------------- //
// ------------------------------------------------------------------ //
/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::vmult (Vec &dst,
    const Vec &src) const
  {
    VecSet(dst, 0);
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::vmult_add (
    Vec &dst,
    const Vec &src) const
  {

    ParallelVector src_vec, dst_vec;
    matfree_data.initialize_dof_vector(src_vec);
    matfree_data.initialize_dof_vector(dst_vec);

    copy_to_Vector(src_vec, src);
    copy_to_Vector(dst_vec, dst);

    if (boundary_req)
    {
      matfree_data.loop(&PoissonOperator::cell_local_apply,
        &PoissonOperator::cell_face_apply,
        &PoissonOperator::cell_boundary_apply, this, dst_vec, src_vec,
        /*zero_dst =*/false,
        MatrixFree<dim, number>::DataAccessOnFaces::none,
        MatrixFree<dim, number>::DataAccessOnFaces::none);
    }
    else
    {
      matfree_data.cell_loop(&PoissonOperator::cell_local_apply, this,
        dst_vec, src_vec);
    }
    copy_to_Vec(dst, dst_vec);

  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::Tvmult (
    Vec &dst,
    const Vec &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree, typename number>
  void PoissonOperator<dim, n_fe_degree, number>::Tvmult_add (
    Vec &dst,
    const Vec &src) const
  {
    vmult_add(dst, src);
  }

// ----------- Explicit Instantations ----------- //
template class PoissonOperator<1, 1, double> ;
template class PoissonOperator<1, 2, double> ;
template class PoissonOperator<1, 3, double> ;
template class PoissonOperator<1, 4, double> ;
template class PoissonOperator<1, 5, double> ;

template class PoissonOperator<2, 1, double> ;
template class PoissonOperator<2, 2, double> ;
template class PoissonOperator<2, 3, double> ;
template class PoissonOperator<2, 4, double> ;
template class PoissonOperator<2, 5, double> ;

template class PoissonOperator<3, 1, double> ;
template class PoissonOperator<3, 2, double> ;
template class PoissonOperator<3, 3, double> ;
template class PoissonOperator<3, 4, double> ;
template class PoissonOperator<3, 5, double> ;

