/**
 * @file   matrix_operators_complex_base.cc
 * @brief  Implementation of TransportMatrixComplexBase and MassMatrixComplexBase
 *  classes to handle complex block matrices.
 */

#include "../../include/matrix_operators/matrix_operators_complex_base.h"
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
 *
 */
template <int dim, int n_fe_degree>
  TransportMatrixComplexBase<dim, n_fe_degree>::TransportMatrixComplexBase (
    const MPI_Comm &_comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double>&) :
      comm(_comm),
      dof_handler(dh)
  {
    // Silly initializations
    n_dofs_block = 0;
    n_blocks = 0;
    n_blocks_real = 0;
    matrixfree_type = non_diagonal;
    n_moments = 1;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::clear ()
  {

    if (matrixfree_type == non_diagonal)
    {
      for (unsigned int row = 0; row < n_blocks_real; ++row)
        for (unsigned int col = 0; col < n_blocks_real; ++col)
        {
          if (row == col)
          {
            matrix_blocks_real[row][col]->clear();
            delete matrix_blocks_real[row][col];
          }
          else
          {
            mass_mf_blocks_real[row][col]->clear();
            delete mass_mf_blocks_real[row][col];

            mass_mf_blocks_imag[row][col]->clear();
            delete mass_mf_blocks_imag[row][col];
          }
        }
    }
    else if (matrixfree_type == full_matrixfree)
    {
      for (unsigned int row = 0; row < n_blocks_real; ++row)
        for (unsigned int col = 0; col < n_blocks_real; ++col)
          if (row != col)
          {
            mass_mf_blocks_real[row][col]->clear();
            delete mass_mf_blocks_real[row][col];

            mass_mf_blocks_imag[row][col]->clear();
            delete mass_mf_blocks_imag[row][col];
          }
          else
          {
            poison_mf_blocks_real[row]->clear();
            delete poison_mf_blocks_real[row];
          }
    }
    else
      // full allocated
      for (unsigned int row = 0; row < n_blocks_real; ++row)
        for (unsigned int col = 0; col < n_blocks_real; ++col)
        {
          matrix_blocks_real[row][col]->clear();
          matrix_blocks_imag[row][col]->clear();
        }

    n_blocks = 0;
    n_blocks_real = 0;
  }

/**
 * @brief Return the total number of columns in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixComplexBase<dim, n_fe_degree>::n () const
  {
    return n_blocks * n_dofs_block;
  }

/**
 * @brief Return the total number of rows in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixComplexBase<dim, n_fe_degree>::m () const
  {
    return n_blocks * n_dofs_block;
  }

/**
 * @brief Return the number of row blocks in a column.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixComplexBase<dim, n_fe_degree>::n_blocks_rows () const
  {
    return n_blocks;
  }

/**
 * @brief Return the number of column blocks in a row.
 * @return The number of column blocks in every row.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixComplexBase<dim, n_fe_degree>::n_blocks_cols () const
  {
    return n_blocks;
  }

/**
 * @brief Return the number of dofs per block.
 * @return The number of dofs per block.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixComplexBase<dim, n_fe_degree>::n_dofs_blocks () const
  {
    return n_dofs_block;
  }

/**
 * @brief Return the number of dofs per block.
 * @return The number of dofs per block.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixComplexBase<dim, n_fe_degree>::get_n_moments () const
  {
    return n_moments;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PETScWrappers::MPI::SparseMatrix&
  TransportMatrixComplexBase<dim, n_fe_degree>::block (unsigned int i,
    unsigned int j)
  {
    AssertRelease(
      matrixfree_type == full_allocated or (i == j and matrixfree_type == non_diagonal),
      "Block not accessible in this matrix-free mode");

    if (is_block_real(i, j))
    {
      return *(matrix_blocks_real[i / 2][j / 2]);

    }
    else if (is_block_imag_negative(i, j))
    {
      AssertRelease(false, "Negative not Imag block not implimented to extract");
      //return -*(matrix_blocks_imag[i/2][j/2]);
    }
    else
    {
      return *(matrix_blocks_imag[i / 2][j / 2]);
    }
    return *(matrix_blocks_imag[0][0]);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::get_inv_diagonal (
    const unsigned int block,
    PETScWrappers::MPI::Vector &inv_diagonal)
  {

    if (matrixfree_type == full_matrixfree)
    {

      if (this->poison_mf_blocks_real[block / 2]->inverse_diagonal.size() == 0)
      {
        this->poison_mf_blocks_real[block / 2]->compute_diagonal();
      }

      if (inv_diagonal.size() != this->n_dofs_block)
        inv_diagonal.reinit(comm, this->n_dofs_block,
          this->locally_owned_dofs.n_elements());

      copy_to_Vector(inv_diagonal,
        this->poison_mf_blocks_real[block / 2]->inverse_diagonal);
    }
    else
    {

      inv_diagonal.reinit(comm, this->n_dofs_block,
        this->locally_owned_dofs.n_elements());

      IndexSet index_set(inv_diagonal.locally_owned_elements());
      int i = 0;
      for (IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++, i++)
      {
        double diag_value = this->matrix_blocks_real[block / 2][block / 2]->diag_element(
          *it);
        inv_diagonal[*it] = 1. / diag_value;
      }

      inv_diagonal.compress(VectorOperation::insert);

    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::get_inv_diagonal (
    PETScWrappers::MPI::BlockVector &inv_mat)
  {
    inv_mat.reinit(n_blocks, comm, this->n_dofs_block,
      this->locally_owned_dofs.n_elements());

    for (unsigned int nb = 0; nb < n_blocks; nb++)
    {
      get_inv_diagonal(nb, inv_mat.block(nb));
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::get_inv_diagonal (
    const unsigned int block,
    DiagonalMatrix<PETScWrappers::MPI::Vector> &inv_mat)
  {

    if (this->poison_mf_blocks_real[block / 2]->inverse_diagonal.size() == 0)
    {
      this->poison_mf_blocks_real[block / 2]->compute_diagonal();
    }

    inv_mat.get_vector().reinit(comm, this->n_dofs_block,
      this->locally_owned_dofs.n_elements());

    copy_to_Vector(inv_mat.get_vector(),
      this->poison_mf_blocks_real[block / 2]->inverse_diagonal);

    inv_mat.get_vector().compress(VectorOperation::insert);

    // Delete diagonal
    this->poison_mf_blocks_real[block / 2]->inverse_diagonal.reinit(0);

  }

// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector Complete Multiplications  --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_k3 (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add_k3(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_k3 (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add_k3(i, j, dst.block(i), src.block(j));
      }
  }

// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector Real and Imaginary parts  --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_real (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks_real; i++)
      dst.block(i) = 0.0;

    this->vmult_add_real(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_real (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks_real);
    AssertDimension(src.n_blocks(), n_blocks_real);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks_real; i++)
      for (unsigned int j = 0; j < n_blocks_real; j++)
      {
        this->vmult_add_real(i, j, dst.block(i), src.block(j));
      }
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_imag (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks_real; i++)
      dst.block(i) = 0.0;

    this->vmult_add_imag(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_imag (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks_real);
    AssertDimension(src.n_blocks(), n_blocks_real);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks_real; i++)
      for (unsigned int j = 0; j < n_blocks_real; j++)
      {
        this->vmult_add_imag(i, j, dst.block(i), src.block(j));
      }
  }

// ----------------------------------------------------------------- //
// ----------- Block Real and Imag Multiplications ----------------- //
// ------------------ PETScWrappers::MPI::Vector ------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_real (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or
        (matrixfree_type == non_diagonal and row != col))
    {
      mass_mf_blocks_real[row][col]->vmult(dst, src);
    }
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks_real[row]->vmult(dst, src);
    else // full_allocated
    {
      matrix_blocks_real[row][col]->vmult(dst, src);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_real (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or
        (matrixfree_type == non_diagonal and row != col))
    {
      mass_mf_blocks_real[row][col]->vmult_add(dst, src);
    }
    else // full_allocated
    {
      matrix_blocks_real[row][col]->vmult_add(dst, src);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_imag (
    const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or matrixfree_type == non_diagonal)
    {
      mass_mf_blocks_imag[row][col]->vmult(dst, src);
    }
    else // full_allocated
    {
      matrix_blocks_imag[row][col]->vmult(dst, src);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_imag (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or matrixfree_type == non_diagonal)
    {
      mass_mf_blocks_imag[row][col]->vmult_add(dst, src);
    }
    else // full_allocated
    {
      matrix_blocks_imag[row][col]->vmult_add(dst, src);
    }
  }

//////////////////////////////////////////////////////////////////////////////

template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_group (
    const unsigned int g1,
    const unsigned int g2,
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertIndexRange(g1, n_blocks);
    AssertIndexRange(g2, n_blocks);
    PETScWrappers::MPI::Vector inter = src.block(1);
    inter *= -1;

    if (matrixfree_type == full_matrixfree or
        (matrixfree_type == non_diagonal and g1 != g2))
    {
      mass_mf_blocks_real[g1][g2]->vmult(dst.block(0), src.block(0));
      mass_mf_blocks_imag[g1][g2]->vmult_add(dst.block(0), inter);
      mass_mf_blocks_imag[g1][g2]->vmult(dst.block(1), src.block(0));
      mass_mf_blocks_real[g1][g2]->vmult_add(dst.block(1), src.block(1));
    }
    else if (matrixfree_type == full_matrixfree and g1 == g2)
    {
      poison_mf_blocks_real[g1]->vmult(dst.block(0), src.block(0));
      mass_mf_blocks_imag[g1][g2]->vmult_add(dst.block(0), inter);
      mass_mf_blocks_imag[g1][g2]->vmult(dst.block(1), src.block(0));
      poison_mf_blocks_real[g1]->vmult_add(dst.block(1), src.block(1));
    }
    else // full_allocated
    {
      matrix_blocks_real[g1][g2]->vmult(dst.block(0), src.block(0));
      matrix_blocks_imag[g1][g2]->vmult_add(dst.block(0), inter);
      matrix_blocks_imag[g1][g2]->vmult(dst.block(1), src.block(0));
      matrix_blocks_real[g1][g2]->vmult_add(dst.block(1), src.block(1));
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_group (
    const unsigned int g1,
    const unsigned int g2,
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertIndexRange(g1, n_blocks);
    AssertIndexRange(g2, n_blocks);
    PETScWrappers::MPI::Vector inter = src.block(1);
    inter *= -1;

    if (matrixfree_type == full_matrixfree or
        (matrixfree_type == non_diagonal and g1 != g2))
    {
      mass_mf_blocks_real[g1][g2]->vmult_add(dst.block(0), src.block(0));
      mass_mf_blocks_imag[g1][g2]->vmult_add(dst.block(0), inter);
      mass_mf_blocks_imag[g1][g2]->vmult_add(dst.block(1), src.block(0));
      mass_mf_blocks_real[g1][g2]->vmult_add(dst.block(1), src.block(1));
    }
    else if (matrixfree_type == full_matrixfree and g1 == g2)
    {
      poison_mf_blocks_real[g1]->vmult_add(dst.block(0), src.block(0));
      mass_mf_blocks_imag[g1][g2]->vmult_add(dst.block(0), inter);
      mass_mf_blocks_imag[g1][g2]->vmult_add(dst.block(1), src.block(0));
      poison_mf_blocks_real[g1]->vmult_add(dst.block(1), src.block(1));
    }
    else // full_allocated
    {
      matrix_blocks_real[g1][g2]->vmult_add(dst.block(0), src.block(0));
      matrix_blocks_imag[g1][g2]->vmult_add(dst.block(0), inter);
      matrix_blocks_imag[g1][g2]->vmult_add(dst.block(1), src.block(0));
      matrix_blocks_real[g1][g2]->vmult_add(dst.block(1), src.block(1));
    }
  }

// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector Transpose  Multiplications  --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_transpose (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(j, i, dst.block(i), src.block(j));
      }
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_transpose (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add_transpose(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::Tvmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add_transpose(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  double TransportMatrixComplexBase<dim, n_fe_degree>::vmult_dot (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    PETScWrappers::MPI::BlockVector inter(src);

    vmult(inter, src);
    double val = dst * inter;

    for (unsigned int g = 0; g < inter.n_blocks(); g++)
      inter.block(g).clear();

    return val;

  }

// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Complete Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Transpose Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_transpose (
    ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(j, i, dst.block(i), src.block(j));
      }
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_transpose (
    ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add_transpose(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  double TransportMatrixComplexBase<dim, n_fe_degree>::vmult_dot (
    ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    ParallelBlockVector inter(src);

    vmult(inter, src);
    double val = dst * inter;

    return val;

  }

// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------ PETScWrappers::MPI::Vector ------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or (matrixfree_type == non_diagonal
                                               and row != col))
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
        dst *= -1.0;
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks_real[row / 2]->vmult(dst, src);
    else // full_allocated
    {
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult(dst, src);
        dst *= -1.0;
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult(dst, src);
      }
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree)
        and row != col)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks_real[row / 2]->vmult_add(dst, src);
    else // full_allocated
    {
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src;
        inter *= -1.0;
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add_k3 (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree)
        and row != col)
    {
      if (is_block_imag_k3(row, col))
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_negative_real_k3(row, col))
      {
        PETScWrappers::MPI::Vector inter = src;
        inter *= -1.0;
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks_real[row / 2]->vmult_add(dst, src);
    else // full_allocated
    {

      if (is_block_imag_k3(row, col))
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_negative_real_k3(row, col))
      {
        PETScWrappers::MPI::Vector inter = src;
        inter *= -1.0;
        matrix_blocks_real[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        matrix_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }

    }
  }
// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------------------ Vec ------------------------------ //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    Vec &dst,
    const Vec &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or (matrixfree_type == non_diagonal
                                               and row != col))
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
        VecScale(dst, -1);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks_real[row / 2]->vmult(dst, src);
    else // full_allocated
    {
      AssertRelease(false, "NOT IMPLEMENTED");
//          if (is_block_real(row, col))
//          {
//            matrix_blocks_real[row / 2][col / 2]->vmult(dst, src);
//          }
//          else if (is_block_imag_negative(row, col))
//          {
//            matrix_blocks_imag[row / 2][col / 2]->vmult(dst, src);
//            dst *= -1.0;
//          }
//          else
//          {
//            matrix_blocks_imag[row / 2][col / 2]->vmult(dst, src);
//          }
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    Vec &dst,
    const Vec &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree or (matrixfree_type == non_diagonal
                                               and row != col))
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        Vec inter;
        VecCopy(src, inter);
        VecScale(inter, -1);
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
        VecDestroy(&inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks_real[row / 2]->vmult_add(dst, src);
    else // full_allocated
    {
      AssertRelease(false, "NOT IMPLEMENTED");
//          if (is_block_real(row, col))
//          {
//            matrix_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
//          }
//          else if (is_block_imag_negative(row, col))
//          {
//            matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
//            dst *= -1.0;
//          }
//          else
//          {
//            matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
//          }
    }

  }

// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------------ Vector<double> ------------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree)
        and row != col)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        ParallelVector inter = src;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks_real[row / 2]->vmult_add(dst, src);
    else // full_allocated
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult_add(dst_vec, src_vec);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src_vec;
        inter *= -1.0;
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst_vec, inter);
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst_vec, src_vec);
      }
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
  void TransportMatrixComplexBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == full_matrixfree
         or matrixfree_type == non_diagonal)
        and row != col)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        ParallelVector inter = src;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
      }
    }
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks_real[row / 2]->vmult(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult(dst_vec, src_vec);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src_vec;
        inter *= -1.0;
        matrix_blocks_imag[row / 2][col / 2]->vmult(dst_vec, inter);
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult(dst_vec, src_vec);
      }

      //dst = dst_vec;

      LinearAlgebra::ReadWriteVector<double> rwv(dst_vec.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);
    }
  }

template <int dim, int n_fe_degree>
  void TransportMatrixComplexBase<dim, n_fe_degree>::print_matrix_in_matlab (
    unsigned int block_i,
    unsigned int block_j,
    const std::string &mat_name,
    std::ostream &out,
    const unsigned int precision)
  {
    out << std::setprecision(precision);

    if (is_block_real(block_i, block_j))
    {
      block_i /= 2;
      block_j /= 2;
      unsigned int nnz = matrix_blocks_real[block_i][block_j]->n_nonzero_elements();
      out << "Nnz = " << nnz << ";" << std::endl;
      std::vector<unsigned int> i_index;
      std::vector<unsigned int> j_index;
      std::vector<double> mat_values;

      i_index.reserve(nnz);
      j_index.reserve(nnz);
      mat_values.reserve(nnz);

      for (unsigned int i = 0; i < matrix_blocks_real[block_i][block_j]->n(); ++i)
        for (unsigned int j = 0; j < matrix_blocks_real[block_i][block_j]->m(); ++j)
          if (std::abs(matrix_blocks_real[block_i][block_j]->el(i, j)) > 1e-15)
          {
            i_index.push_back(i + 1);
            j_index.push_back(j + 1);
            mat_values.push_back(matrix_blocks_real[block_i][block_j]->el(i, j));
          }

      // Print in the vectors in the out
      print_vector_in_matlab(i_index, "i_index", out, precision);
      print_vector_in_matlab(j_index, "j_index", out, precision);
      print_vector_in_matlab(mat_values, "mat_values", out, precision);

      out << "m = " << matrix_blocks_real[block_i][block_j]->m() << ";" << std::endl;
      out << "n = " << matrix_blocks_real[block_i][block_j]->n() << ";" << std::endl;
      out << mat_name << " = sparse(i_index, j_index, mat_values, m, n, Nnz);"
          << std::endl;
      out << "clear  i_index; clear  j_index;" << std::endl;
      out << "clear mat_values; clear m; clear n; clear Nnz;" << std::endl;
      out << std::endl;
    }
    else
    {
      block_i /= 2;
      block_j /= 2;

      unsigned int nnz = matrix_blocks_imag[block_i][block_j]->n_nonzero_elements();
      out << "Nnz = " << nnz << ";" << std::endl;
      std::vector<unsigned int> i_index;
      std::vector<unsigned int> j_index;
      std::vector<double> mat_values;

      i_index.reserve(nnz);
      j_index.reserve(nnz);
      mat_values.reserve(nnz);

      for (unsigned int i = 0; i < matrix_blocks_imag[block_i][block_j]->n(); ++i)
        for (unsigned int j = 0; j < matrix_blocks_imag[block_i][block_j]->m(); ++j)
          if (std::abs(matrix_blocks_imag[block_i][block_j]->el(i, j)) > 1e-15)
          {
            i_index.push_back(i + 1);
            j_index.push_back(j + 1);
            mat_values.push_back(matrix_blocks_imag[block_i][block_j]->el(i, j));
          }

      // Print in the vectors in the out
      print_vector_in_matlab(i_index, "i_index", out, precision);
      print_vector_in_matlab(j_index, "j_index", out, precision);
      print_vector_in_matlab(mat_values, "mat_values", out, precision);

      out << "m = " << matrix_blocks_imag[block_i][block_j]->m() << ";" << std::endl;
      out << "n = " << matrix_blocks_imag[block_i][block_j]->n() << ";" << std::endl;
      out << mat_name << " = sparse(i_index, j_index, mat_values, m, n, Nnz);"
          << std::endl;
      out << "clear  i_index; clear  j_index;" << std::endl;
      out << "clear mat_values; clear m; clear n; clear Nnz;" << std::endl;
      out << std::endl;
    }
  }

/**
 * @brief Is this block Imag in k3 formulation?
 * @return true or false
 */
template <int dim, int n_fe_degree>
  bool TransportMatrixComplexBase<dim, n_fe_degree>::is_block_imag_k3 (
    const unsigned int row,
    const unsigned int col) const
  {

    return (row % 2 == 0 and col % 2 == 0) or (row % 2 == 1 and col % 2 == 1);
  }

/**
 * @brief Is this block real negative in k3 formulation?
 * @return true or false
 */
template <int dim, int n_fe_degree>
  bool TransportMatrixComplexBase<dim, n_fe_degree>::is_block_negative_real_k3 (
    const unsigned int row,
    const unsigned int col) const
  {

    return (row % 2 == 1 and col % 2 == 0);
  }

/**
 * @brief Is this block real?
 * @return true or false
 */
template <int dim, int n_fe_degree>
  bool TransportMatrixComplexBase<dim, n_fe_degree>::is_block_real (
    const unsigned int row,
    const unsigned int col) const
  {

    return (row % 2 == 0 and col % 2 == 0) or (row % 2 == 1 and col % 2 == 1);
  }

/**
 * @brief Is this block real?
 * @return true or false
 */
template <int dim, int n_fe_degree>
  bool TransportMatrixComplexBase<dim, n_fe_degree>::is_block_imag_negative (
    const unsigned int row,
    const unsigned int col) const
  {
    return (row % 2 == 0 and col % 2 == 1);
  }

//-------------------------------------------------------------------------------------------//
//  ---------------------------------------------------------------------------------------- //
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints
 */
template <int dim, int n_fe_degree>
  MassMatrixComplexBase<dim, n_fe_degree>::MassMatrixComplexBase (const MPI_Comm &_comm,
    const DoFHandler<dim>&,
    const AffineConstraints<double>&) :
      comm(_comm)
  {
    // Silly initialization
    n_blocks = 0;
    n_blocks_real = 0;
    n_dofs_block = 0;
    matrixfree_type = non_diagonal;
  }

/**
 * @brief
 */ // TODO
template <int dim, int n_fe_degree>
  PETScWrappers::MPI::SparseMatrix&
  MassMatrixComplexBase<dim, n_fe_degree>::block (unsigned int i,
    unsigned int j)
  {
    AssertRelease(
      matrixfree_type == full_allocated
      or (i == j and matrixfree_type == non_diagonal),
      "Block not accessible in this matrix-free mode");

    if (is_block_real(i, j))
    {
      return *(matrix_blocks_real[i / 2][j / 2]);

    }
    else if (is_block_imag_negative(i, j))
    {
      AssertRelease(false, "Negative not Imag block not implimented to extract");
      //return -*(matrix_blocks_imag[i/2][j/2]);
    }
    else
    {
      return *(matrix_blocks_imag[i / 2][j / 2]);
    }

    return *(matrix_blocks_real[0][0]);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::clear ()
  {
    if (matrixfree_type == full_matrixfree)
    {
      for (unsigned int row = 0; row < n_blocks_real; ++row)
        for (unsigned int col = 0; col < n_blocks_real; ++col)
        {
          mass_mf_blocks_real[row][col]->clear();
          delete mass_mf_blocks_real[row][col];

          mass_mf_blocks_imag[row][col]->clear();
          delete mass_mf_blocks_imag[row][col];
        }
    }
    else
    {
      for (unsigned int row = 0; row < n_blocks_real; ++row)
        for (unsigned int col = 0; col < n_blocks_real; ++col)
        {
          matrix_blocks_real[row][col]->clear();
          matrix_blocks_imag[row][col]->clear();
        }
    }
  }

/**
 * @brief Return the number of row blocks in a column
 */
template <int dim, int n_fe_degree>
  unsigned int MassMatrixComplexBase<dim, n_fe_degree>::n_blocks_rows () const
  {
    return n_blocks;
  }

/**
 * @brief Return the number of column blocks in a row.
 */
template <int dim, int n_fe_degree>
  unsigned int MassMatrixComplexBase<dim, n_fe_degree>::n_blocks_cols () const
  {
    return n_blocks;
  }

/**
 * @brief Return the total number of columns in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int MassMatrixComplexBase<dim, n_fe_degree>::n () const
  {
    return n_blocks * n_dofs_block;
  }

/**
 * @brief Return the total number of rows in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int MassMatrixComplexBase<dim, n_fe_degree>::m () const
  {
    return n_blocks * n_dofs_block;
  }

// ----------------------------------------------------------------- //
// ------- PETScWrappers::MPI::BlockVector Multiplications  -------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector Transpose Multiplications --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_transpose (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add_transpose(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add_transpose (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(j, i, dst.block(i), src.block(j));
      }
  }

/**
 * @brief Complete matrix scalar product
 *  val = dst * A * src
 */
template <int dim, int n_fe_degree>
  double MassMatrixComplexBase<dim, n_fe_degree>::vmult_dot (
    const PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    PETScWrappers::MPI::BlockVector inter(src);

    vmult(inter, src);
    double val = dst * inter;

    for (unsigned int g = 0; g < inter.n_blocks(); g++)
      inter.block(g).clear();

    return val;

  }

// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Complete Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Complete Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_transpose (
    ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add_transpose(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + MassMatrixComplexBase * src
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add_transpose (
    ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(dst.size(), this->m());
    AssertDimension(src.size(), this->m());

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        this->vmult_add(j, i, dst.block(i), src.block(j));
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
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
        dst *= -1.0;
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult(dst, src);
      }
    }
    else
    {
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult(dst, src);
        dst *= -1.0;
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult(dst, src);
      }
    }
  }

/**
 * @brief Is this block real?
 * @return true or false
 */
template <int dim, int n_fe_degree>
  bool MassMatrixComplexBase<dim, n_fe_degree>::is_block_real (
    const unsigned int row,
    const unsigned int col) const
  {

    return (row % 2 == 0 and col % 2 == 0) or (row % 2 == 1 and col % 2 == 1);
  }

/**
 * @brief Is this block real?
 * @return true or false
 */
template <int dim, int n_fe_degree>
  bool MassMatrixComplexBase<dim, n_fe_degree>::is_block_imag_negative (
    const unsigned int row,
    const unsigned int col) const
  {
    return (row % 2 == 0 and col % 2 == 1);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else
    {
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src;
        inter *= -1.0;
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
  }

/**
 * v1 = v2 + A * v3
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &v1,
    PETScWrappers::MPI::Vector &v2,
    const PETScWrappers::MPI::Vector &v3) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    v1 = v2;

    if (matrixfree_type == full_matrixfree)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(v1, v3);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = v3;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(v1, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(v1, v3);
      }
    }
    else
    {
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult_add(v1, v3);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = v3;
        inter *= -1.0;
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(v1, inter);
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(v1, v3);
      }
    }
  }

// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------------ ParallelVector ------------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult_add (
    const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        ParallelVector inter = src;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      if (is_block_real(row, col))
      {
        matrix_blocks_real[row / 2][col / 2]->vmult_add(dst_vec, src_vec);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src_vec;
        inter *= -1.0;
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst_vec, inter);
      }
      else
      {
        matrix_blocks_imag[row / 2][col / 2]->vmult_add(dst_vec, src_vec);
      }
      //      dst = dst_vec;
      LinearAlgebra::ReadWriteVector<double> rwv(dst.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);

    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void MassMatrixComplexBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == full_matrixfree)
    {
      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst, src);
      }
      else if (is_block_imag_negative(row, col))
      {
        ParallelVector inter = src;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst, src);
      }
    }
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);

      if (is_block_real(row, col))
      {
        mass_mf_blocks_real[row / 2][col / 2]->vmult_add(dst_vec, src_vec);
      }
      else if (is_block_imag_negative(row, col))
      {
        PETScWrappers::MPI::Vector inter = src_vec;
        inter *= -1.0;
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst_vec, inter);
      }
      else
      {
        mass_mf_blocks_imag[row / 2][col / 2]->vmult_add(dst_vec, src_vec);
      }
      //dst = dst_vec;
      LinearAlgebra::ReadWriteVector<double> rwv(dst_vec.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);
    }

  }

// ----------- Explicit Instantations ----------- //

template class TransportMatrixComplexBase<1, 1> ;
template class TransportMatrixComplexBase<1, 2> ;
template class TransportMatrixComplexBase<1, 3> ;
template class TransportMatrixComplexBase<1, 4> ;
template class TransportMatrixComplexBase<1, 5> ;

template class TransportMatrixComplexBase<2, 1> ;
template class TransportMatrixComplexBase<2, 2> ;
template class TransportMatrixComplexBase<2, 3> ;
template class TransportMatrixComplexBase<2, 4> ;
template class TransportMatrixComplexBase<2, 5> ;

template class TransportMatrixComplexBase<3, 1> ;
template class TransportMatrixComplexBase<3, 2> ;
template class TransportMatrixComplexBase<3, 3> ;
template class TransportMatrixComplexBase<3, 4> ;
template class TransportMatrixComplexBase<3, 5> ;

template class MassMatrixComplexBase<1, 1> ;
template class MassMatrixComplexBase<1, 2> ;
template class MassMatrixComplexBase<1, 3> ;
template class MassMatrixComplexBase<1, 4> ;
template class MassMatrixComplexBase<1, 5> ;

template class MassMatrixComplexBase<2, 1> ;
template class MassMatrixComplexBase<2, 2> ;
template class MassMatrixComplexBase<2, 3> ;
template class MassMatrixComplexBase<2, 4> ;
template class MassMatrixComplexBase<2, 5> ;

template class MassMatrixComplexBase<3, 1> ;
template class MassMatrixComplexBase<3, 2> ;
template class MassMatrixComplexBase<3, 3> ;
template class MassMatrixComplexBase<3, 4> ;
template class MassMatrixComplexBase<3, 5> ;

