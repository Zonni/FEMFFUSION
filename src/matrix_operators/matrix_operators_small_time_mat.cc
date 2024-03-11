/**
 * @file   matrix_operators_small_time_mat.cc
 * @brief
 */

#include "../../include/matrix_operators/matrix_operators_small_time_mat.h"

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
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

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
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */

template <int dim, int n_fe_degree>
  SmallDelayedFissionMatrices<dim, n_fe_degree>::SmallDelayedFissionMatrices (const MPI_Comm &_comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double> &constraints) :
      comm(_comm),
      tria(dh.get_triangulation()),
      dof_handler(dh),
      constraints(
        constraints)
  {
    // Silly initialization

    n_dofs = 0;
    matrixfree_type = full_matrixfree;
    n_prec_groups = 0;
    n_energy_groups = 0;
    listen_to_material_id = false;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const MatrixFreeType &_matrixfree_type,
    bool _listen_to_material_id)
  {

    matrixfree_type = _matrixfree_type;
    n_prec_groups = materials.get_n_precursors();
    n_energy_groups = materials.get_n_groups();
    n_dofs = dof_handler.n_dofs();
    listen_to_material_id = _listen_to_material_id;

    if (matrixfree_type == full_matrixfree)
    {

      reinit_full_matrixfree(materials);

    }
    else if (matrixfree_type == full_allocated)
    {

      matrix_csr.resize(n_energy_groups,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_prec_groups));

      locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        dof_handler.n_locally_owned_dofs_per_processor(), comm,
        locally_relevant_dofs);
      sp.copy_from(dsp);

      for (unsigned ng = 0; ng < n_energy_groups; ng++)
        for (unsigned np = 0; np < n_prec_groups; np++)
        {
          matrix_csr[ng][np] = new PETScWrappers::MPI::SparseMatrix;
          matrix_csr[ng][np]->reinit(locally_owned_dofs,
            locally_owned_dofs, sp, comm);
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
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {

    double val;
    double coeff;
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

      // Get the material coefficients:
      for (unsigned ng = 0; ng < n_energy_groups; ng++)
        for (unsigned np = 0; np < n_prec_groups; np++)
        {
          coeff = materials.get_delayed_fraction(mat, np)
                  * materials.get_nu_sigma_f(ng, mat);
          cell_M.equ(coeff, cell_val);
          constraints.distribute_local_to_global(cell_M,
            local_dof_indices, *(matrix_csr[ng][np]));
        }

    }

    // Compress
    for (unsigned ng = 0; ng < n_energy_groups; ng++)
      for (unsigned np = 0; np < n_prec_groups; np++)
        matrix_csr[ng][np]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials)
  {

    const unsigned int n_mats = materials.get_n_mats();
    mass_mf.resize(n_energy_groups,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(
        n_prec_groups));

    coeffs.resize(n_energy_groups, std::vector<Vector<double> >(n_prec_groups));
    // Fill coeffs
    for (unsigned ng = 0; ng < n_energy_groups; ng++)
      for (unsigned np = 0; np < n_prec_groups; np++)
      {
        coeffs[ng][np].reinit(n_mats);
        for (unsigned int mat = 0; mat < n_mats; mat++)
        {
          coeffs[ng][np][mat] = materials.get_delayed_fraction(mat, np)
                                * materials.get_nu_sigma_f(ng, mat);
        }
      }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_JxW_values);
    matfree_data.reinit(dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    for (unsigned ng = 0; ng < n_energy_groups; ng++)
      for (unsigned np = 0; np < n_prec_groups; np++)
      {
        mass_mf[ng][np] = new MassOperator<dim, n_fe_degree, double>(
          matfree_data);
        mass_mf[ng][np]->reinit(constraints,
          materials.get_materials_vector(), coeffs[ng][np], listen_to_material_id);
      }

  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::clear ()
  {

    matfree_data.clear();
    coeffs.clear();

    if (matrixfree_type == full_matrixfree)
      for (unsigned int g = 0; g < n_energy_groups; ++g)
        for (unsigned int p = 0; p < n_prec_groups; ++p)
        {
          mass_mf[g][p]->clear();
          delete mass_mf[g][p];
        }
    else
      for (unsigned int g = 0; g < n_energy_groups; ++g)
        for (unsigned int p = 0; p < n_prec_groups; ++p)
          matrix_csr[g][p]->clear();
  }

/**
 * @brief Return the total number of columns in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int SmallDelayedFissionMatrices<dim, n_fe_degree>::n () const
  {
    return n_dofs;
  }

/**
 * @brief Return the total number of rows in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int SmallDelayedFissionMatrices<dim, n_fe_degree>::m () const
  {
    return n_dofs;
  }

// ----------------------------------------------------------------- //
// ------- PETScWrappers::MPI::BlockVector Multiplications  -------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = SmallDelayedFissionMatrices * src
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {

    AssertDimension(dst.n_blocks(), n_energy_groups);

    for (unsigned int i = 0; i < n_energy_groups; i++)
      dst.block(i) = 0.0;

    vmult_add(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + SmallDelayedFissionMatrices * src
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult_add (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_energy_groups);
    AssertDimension(src.n_blocks(), n_prec_groups);

    for (unsigned int i = 0; i < n_energy_groups; i++)
      for (unsigned int j = 0; j < n_prec_groups; j++)
      {
        vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Complete Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = SmallDelayedFissionMatrices * src
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_energy_groups);

    for (unsigned int i = 0; i < n_energy_groups; i++)
      dst.block(i) = 0.0;

    vmult_add(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + SmallDelayedFissionMatrices * src
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult_add (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {

    AssertDimension(dst.n_blocks(), n_energy_groups);
    AssertDimension(src.n_blocks(), n_prec_groups);

    for (unsigned int i = 0; i < n_energy_groups; i++)
      for (unsigned int j = 0; j < n_prec_groups; j++)
      {
        vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------ PETScWrappers::MPI::Vector ------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_energy_groups);
    AssertIndexRange(col, n_prec_groups);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf[row][col]->vmult(dst, src);
    else
      matrix_csr[row][col]->vmult(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_energy_groups);
    AssertIndexRange(col, n_prec_groups);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
    {
      mass_mf[row][col]->vmult_add(dst, src);
    }
    else
      matrix_csr[row][col]->vmult_add(dst, src);
  }

/**
 * v1 = v2 + A * v3
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &v1,
    PETScWrappers::MPI::Vector &v2,
    const PETScWrappers::MPI::Vector &v3) const
  {
    AssertIndexRange(row, n_energy_groups);
    AssertIndexRange(col, n_prec_groups);

    v1 = v2;

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf[row][col]->vmult_add(v1, v3);
    else
      matrix_csr[row][col]->vmult_add(v1, v3);
  }

// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------------ ParallelVector ------------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_energy_groups);
    AssertIndexRange(col, n_prec_groups);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf[row][col]->vmult_add(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      matrix_csr[row][col]->vmult_add(dst_vec, src_vec);
      // dst = dst_vec;
      LinearAlgebra::ReadWriteVector<double> rwv(dst_vec.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);

    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SmallDelayedFissionMatrices<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_energy_groups);
    AssertIndexRange(col, n_prec_groups);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf[row][col]->vmult(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      matrix_csr[row][col]->vmult(dst_vec, src_vec);
      // dst = dst_vec;

      LinearAlgebra::ReadWriteVector<double> rwv(dst_vec.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);
    }
  }

// ----------- Explicit Instantations ----------- //

template class SmallDelayedFissionMatrices<1, 1> ;
template class SmallDelayedFissionMatrices<1, 2> ;
template class SmallDelayedFissionMatrices<1, 3> ;
template class SmallDelayedFissionMatrices<1, 4> ;
template class SmallDelayedFissionMatrices<1, 5> ;

template class SmallDelayedFissionMatrices<2, 1> ;
template class SmallDelayedFissionMatrices<2, 2> ;
template class SmallDelayedFissionMatrices<2, 3> ;
template class SmallDelayedFissionMatrices<2, 4> ;
template class SmallDelayedFissionMatrices<2, 5> ;

template class SmallDelayedFissionMatrices<3, 1> ;
template class SmallDelayedFissionMatrices<3, 2> ;
template class SmallDelayedFissionMatrices<3, 3> ;
template class SmallDelayedFissionMatrices<3, 4> ;
template class SmallDelayedFissionMatrices<3, 5> ;





