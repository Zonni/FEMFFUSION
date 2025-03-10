/**
 * @file   matrix_operators.cc
 * @brief  Implementation of TransportMatrixBase and FissionMAtrix classes to handle block matrices.
 */

#include "../../include/matrix_operators/matrix_operators_petsc.h"
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

MatrixFreeType string_to_enum (std::string matrix_free_string)
{
  lower_case(matrix_free_string);
  if (matrix_free_string == "non_diagonal")
    return non_diagonal;
  else if (matrix_free_string == "full_matrixfree")
    return full_matrixfree;
  else if (matrix_free_string == "full_allocated")
    return full_allocated;
  else
    AssertRelease(false, "Not valid MatrixFreeType " + matrix_free_string);

  return non_diagonal; // Error
}

std::string enum_to_string (MatrixFreeType matrix_free)
{
  if (matrix_free == non_diagonal)
    return "non_diagonal";
  else if (matrix_free == full_matrixfree)
    return "full_matrixfree";
  else if (matrix_free == full_allocated)
    return "full_allocated";

  return "";
}

/**
 * @brief Get a MatrixFreeType from the options commands with the given keyword.
 *    If the keyword is not given the default behavior is not to change the
 *    given result string.
 */
PetscErrorCode get_enum_from_options (const std::string &keyword,
  MatrixFreeType &result)
{
  char result_char[PETSC_MAX_PATH_LEN];
  PetscBool flg;
  PetscErrorCode ierr;

#if PETSC_VERSION_MINOR > 6
  ierr = PetscOptionsGetString(NULL, NULL, keyword.c_str(), result_char,
    sizeof(result_char), &flg);
#else
  ierr = PetscOptionsGetString(NULL, keyword.c_str(), result_char, sizeof(result_char), &flg);
#endif
  CHKERRQ(ierr);

  if (flg)
  {
    std::string result_string = result_char;
    result = string_to_enum(result_string);
  }

  return ierr;
}

/**
 *
 */
template <int dim, int n_fe_degree>
  TransportMatrixBase<dim, n_fe_degree>::TransportMatrixBase (
    const MPI_Comm &_comm,
    const DoFHandler<dim> &dh,
    const AffineConstraints<double>&) :
      comm(_comm),
      dof_handler(dh)
  {
    // Silly initializations
    n_dofs_block = 0;
    n_blocks = 0;
    matrixfree_type = non_diagonal;
    delta_t = 0;
    n_moments = 1;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::clear ()
  {

    if (matrixfree_type == non_diagonal)
    {
      for (unsigned int row = 0; row < n_blocks; ++row)
        for (unsigned int col = 0; col < n_blocks; ++col)
        {
          if (row == col)
          {
            matrix_blocks[row][col]->clear();
            delete matrix_blocks[row][col];

          }
          else
          {
            mass_mf_blocks[row][col]->clear();
            delete mass_mf_blocks[row][col];
          }
        }
    }
    else if (matrixfree_type == full_matrixfree)
    {
      for (unsigned int row = 0; row < n_blocks; ++row)
        for (unsigned int col = 0; col < n_blocks; ++col)
          if (row != col)
          {
            mass_mf_blocks[row][col]->clear();
            delete mass_mf_blocks[row][col];
          }
          else
          {
            poison_mf_blocks[row]->clear();
            delete poison_mf_blocks[row];
          }
    }
    else
      for (unsigned int row = 0; row < n_blocks; ++row)
        for (unsigned int col = 0; col < n_blocks; ++col)
          matrix_blocks[row][col]->clear();

    n_blocks = 0;
  }

/**
 * @brief Return the total number of columns in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixBase<dim, n_fe_degree>::n () const
  {
    return n_blocks * n_dofs_block;
  }

/**
 * @brief Return the total number of rows in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixBase<dim, n_fe_degree>::m () const
  {
    return n_blocks * n_dofs_block;
  }

/**
 * @brief Return the number of row blocks in a column.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixBase<dim, n_fe_degree>::n_blocks_rows () const
  {
    return n_blocks;
  }

/**
 * @brief Return the number of column blocks in a row.
 * @return The number of column blocks in every row.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixBase<dim, n_fe_degree>::n_blocks_cols () const
  {
    return n_blocks;
  }

/**
 * @brief Return the number of dofs per block.
 * @return The number of dofs per block.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixBase<dim, n_fe_degree>::n_dofs_blocks () const
  {
    return n_dofs_block;
  }

/**
 * @brief Return the number of dofs per block.
 * @return The number of dofs per block.
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixBase<dim, n_fe_degree>::get_n_moments () const
  {
    return n_moments;
  }
/**
 *
 */
template <int dim, int n_fe_degree>
  PETScWrappers::MPI::SparseMatrix&
  TransportMatrixBase<dim, n_fe_degree>::block (unsigned int i,
    unsigned int j)
  {

    Assert(
      matrixfree_type == full_allocated or ( i == j and matrixfree_type == non_diagonal),
      ExcMessage("Block not accessible in this matrix-free mode"));
    return *(matrix_blocks[i][j]);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::get_inv_diagonal (
    const unsigned int group,
    PETScWrappers::MPI::Vector &inv_diagonal)
  {

    if (matrixfree_type == full_matrixfree)
    {

      if (this->poison_mf_blocks[group]->inverse_diagonal.size() == 0)
      {
        this->poison_mf_blocks[group]->compute_diagonal();
      }

      if (inv_diagonal.size() != this->n_dofs_block)
        inv_diagonal.reinit(comm, this->n_dofs_block,
          this->locally_owned_dofs.n_elements());

      copy_to_Vector(inv_diagonal,
        this->poison_mf_blocks[group]->inverse_diagonal);
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
        double diag_value = this->matrix_blocks[group][group]->diag_element(
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
  void TransportMatrixBase<dim, n_fe_degree>::get_inv_diagonal (
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
  void TransportMatrixBase<dim, n_fe_degree>::get_inv_diagonal (
    const unsigned int group,
    DiagonalMatrix<PETScWrappers::MPI::Vector> &inv_mat)
  {

    if (this->poison_mf_blocks[group]->inverse_diagonal.size() == 0)
    {
      this->poison_mf_blocks[group]->compute_diagonal();
    }

    inv_mat.get_vector().reinit(comm, this->n_dofs_block,
      this->locally_owned_dofs.n_elements());

    copy_to_Vector(inv_mat.get_vector(),
      this->poison_mf_blocks[group]->inverse_diagonal);

    inv_mat.get_vector().compress(VectorOperation::insert);

    // Delete diagonal
    this->poison_mf_blocks[group]->inverse_diagonal.reinit(0);

  }

// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector Complete Multiplications  --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (
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
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector ADJOINT Multiplications  --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add_transpose (
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
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_transpose (
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
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::Tvmult (
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
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  double TransportMatrixBase<dim, n_fe_degree>::vmult_dot (
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
 * dst = dst + TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (ParallelBlockVector &dst,
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
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult (ParallelBlockVector &dst,
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
 * dst = dst + TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add_transpose (
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
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_transpose (
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
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  double TransportMatrixBase<dim, n_fe_degree>::vmult_dot (
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
  void TransportMatrixBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult(dst, src);
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks[row]->vmult(dst, src);
    else
      matrix_blocks[row][col]->vmult(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult_add(dst, src);
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks[row]->vmult_add(dst, src);
    else
      matrix_blocks[row][col]->vmult_add(dst, src);
  }

/**
 * v1 = v2 + A * v3
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &v1,
    PETScWrappers::MPI::Vector &v2,
    const PETScWrappers::MPI::Vector &v3) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    v1 = v2;

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult_add(v1, v3);
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks[row]->vmult_add(v1, v3);
    else
      matrix_blocks[row][col]->vmult_add(v1, v3);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_row (
    double &dst,
    const PETScWrappers::MPI::BlockVector &src,
    unsigned int row) const
  {

    dst = 0.0;

    this->vmult_add_row(dst, src, row);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add_row (
    double &dst,
    const PETScWrappers::MPI::BlockVector &src,
    unsigned int row) const
  {
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(src.size(), this->m());

    const unsigned int row_block = static_cast<unsigned int>(row / this->n_dofs_block);
    const unsigned int row_index_on_block = row %  this->n_dofs_block;

    for (unsigned int j = 0; j < n_blocks; j++)
    {
      this->vmult_add_row(row_block, j, dst, src.block(j), row_index_on_block);
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
  void TransportMatrixBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    Vec &dst,
    const Vec &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult(dst, src);
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks[row]->vmult(dst, src);
    else
    {
      AssertRelease(false, "NOT IMPLEMENTED");
      //matrix_blocks[row][col]->vmult(dst, src);
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    Vec &dst,
    const Vec &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult_add(dst, src);
    else if (matrixfree_type == full_matrixfree)
      poison_mf_blocks[row]->vmult_add(dst, src);
    else
    {
      AssertRelease(false, "NOT IMPLEMENTED");
      //matrix_blocks[row][col]->vmult_add(dst, src);
    }
  }

/**
 * v1 = v2 + A * v3
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    Vec &v1,
    Vec &v2,
    const Vec &v3) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    v1 = v2;

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult_add(v1, v3);
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks[row]->vmult_add(v1, v3);
    else
    {
      AssertRelease(false, "NOT IMPLEMENTED");
      //matrix_blocks[row][col]->vmult_add(v1, v3);
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
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == non_diagonal
         or matrixfree_type == full_matrixfree)
        and row != col)
      mass_mf_blocks[row][col]->vmult_add(dst, src);
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks[row]->vmult_add(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      matrix_blocks[row][col]->vmult_add(dst_vec, src_vec);
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
  void TransportMatrixBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if ((matrixfree_type == full_matrixfree
         or matrixfree_type == non_diagonal)
        and row != col)
      mass_mf_blocks[row][col]->vmult(dst, src);
    else if (matrixfree_type == full_matrixfree and row == col)
      poison_mf_blocks[row]->vmult(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      matrix_blocks[row][col]->vmult_add(dst_vec, src_vec);
      //dst = dst_vec;

      LinearAlgebra::ReadWriteVector<double> rwv(dst_vec.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixBase<dim, n_fe_degree>::vmult_add_row (
    const unsigned int row,
    const unsigned int col,
    double &dst,
    const PETScWrappers::MPI::Vector &src,
    unsigned int row_index_on_block) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    //if (matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree)
    //  mass_mf_blocks[row][col]->vmult_add_row(dst, src, dst_row);
    //else FIXME
    //  matrix_blocks[row][col]->vmult_add(dst, src, dst_row);
    AssertRelease(matrixfree_type == full_allocated,
      "LUPOD ONLY IMPLEMENTED FOR full_allocated");

    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    MatGetRow((*matrix_blocks[row][col]), row_index_on_block, &ncols, &cols, &vals);
    for (PetscInt j = 0; j < ncols; j++)
      dst += src[cols[j]] * vals[j];

    MatRestoreRow((*matrix_blocks[row][col]), row_index_on_block, &ncols, &cols, &vals);
  }

//-------------------------------------------------------------------------------------------//
//  ---------------------------------------------------------------------------------------- //
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints
 */
template <int dim, int n_fe_degree>
  FisionMatrixBase<dim, n_fe_degree>::FisionMatrixBase (const MPI_Comm &_comm,
    const DoFHandler<dim>&,
    const AffineConstraints<double>&) :
      comm(_comm)
  {
    // Silly initialization
    n_blocks = 0;
    n_dofs_block = 0;
    matrixfree_type = non_diagonal;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PETScWrappers::MPI::SparseMatrix&
  FisionMatrixBase<dim, n_fe_degree>::block (unsigned int i,
    unsigned int j)
  {
    AssertRelease(
      matrixfree_type == full_allocated
      or (i == j and matrixfree_type == non_diagonal),
      "Block not accessible in this matrix-free mode");
    return *(matrix_blocks[i][j]);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::clear ()
  {
    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      for (unsigned int row = 0; row < n_blocks; ++row)
        for (unsigned int col = 0; col < n_blocks; ++col)
        {
          mass_mf_blocks[row][col]->clear();
          delete mass_mf_blocks[row][col];
        }
    else
      for (unsigned int row = 0; row < n_blocks; ++row)
        for (unsigned int col = 0; col < n_blocks; ++col)
          matrix_blocks[row][col]->clear();
  }

/**
 * @brief Return the number of row blocks in a column
 */
template <int dim, int n_fe_degree>
  unsigned int FisionMatrixBase<dim, n_fe_degree>::n_blocks_rows () const
  {
    return n_blocks;
  }

/**
 * @brief Return the number of column blocks in a row.
 */
template <int dim, int n_fe_degree>
  unsigned int FisionMatrixBase<dim, n_fe_degree>::n_blocks_cols () const
  {
    return n_blocks;
  }

/**
 * @brief Return the total number of columns in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int FisionMatrixBase<dim, n_fe_degree>::n () const
  {
    return n_blocks * n_dofs_block;
  }

/**
 * @brief Return the total number of rows in this matrix.
 */
template <int dim, int n_fe_degree>
  unsigned int FisionMatrixBase<dim, n_fe_degree>::m () const
  {
    return n_blocks * n_dofs_block;
  }

// ----------------------------------------------------------------- //
// ------- PETScWrappers::MPI::BlockVector Multiplications  -------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_row (
    double &dst,
    const PETScWrappers::MPI::BlockVector &src,
    unsigned int row) const
  {
    AssertIndexRange(row, src.size());

    dst = 0.0;
    this->vmult_add_row(dst, src, row);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add_row (
    double &dst,
    const PETScWrappers::MPI::BlockVector &src,
    unsigned int row) const
  {
    AssertDimension(src.n_blocks(), n_blocks);
    AssertDimension(src.size(), this->m());
    AssertIndexRange(row, src.size());

    const unsigned int row_block = static_cast<unsigned int>(row / this->n_dofs_block);
    const unsigned int row_index_on_block = row % this->n_dofs_block;
       for (unsigned int j = 0; j < n_blocks; j++)
    {
      this->vmult_add_row(row_block, j, dst, src.block(j), row_index_on_block);
    }
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add (
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
 * dst = TimeMassMatrix * src
 // */
//template <int dim, int n_fe_degree>
//  void FisionMatrixBase<dim, n_fe_degree>::vmult_velocities (
//    PETScWrappers::MPI::BlockVector &dst,
//    const PETScWrappers::MPI::BlockVector &src) const
//  {
//    AssertDimension(dst.n_blocks(), n_blocks);
//
//    for (unsigned int i = 0; i < n_blocks; i++)
//      dst.block(i) = 0.0;
//
//    for (unsigned int i = 0; i < n_blocks; i++)
//    {
//      this->vmult_add(i, i, dst.block(i), src.block(i));
//      dst.block(i) /= velocities[i];
//    }
//
//  }
// ----------------------------------------------------------------- //
// --- PETScWrappers::MPI::BlockVector Transpose Multiplications --- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_transpose (
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
 * dst = dst + FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add_transpose (
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
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  double FisionMatrixBase<dim, n_fe_degree>::vmult_dot (
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

/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
//template <int dim, int n_fe_degree>
//  double FisionMatrixBase<dim, n_fe_degree>::vmult_dot_velocities (
//    PETScWrappers::MPI::BlockVector &dst,
//    const PETScWrappers::MPI::BlockVector &src) const
//  {
//    AssertDimension(dst.n_blocks(), n_blocks);
//
//    PETScWrappers::MPI::BlockVector inter(src);
//    vmult(inter, src);
//    for (unsigned int nb = 0; nb < n_blocks; nb++)
//      inter.block(nb) /= velocities[nb];
//
//    double val = dst * inter;
//
//    for (unsigned int g = 0; g < inter.n_blocks(); g++)
//      inter.block(g).clear();
//
//    return val;
//
//  }
// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Complete Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    for (unsigned int i = 0; i < n_blocks; i++)
      dst.block(i) = 0.0;

    this->vmult_add(dst, src);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = dst + FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add (ParallelBlockVector &dst,
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
 * dst = FisionMatrixBase * src
 */
//template <int dim, int n_fe_degree>
//  void FisionMatrixBase<dim, n_fe_degree>::vmult_velocities (
//    ParallelBlockVector &dst,
//    const ParallelBlockVector &src) const
//  {
//    AssertDimension(dst.n_blocks(), n_blocks);
//
//    for (unsigned int i = 0; i < n_blocks; i++)
//      dst.block(i) = 0.0;
//
//    for (unsigned int i = 0; i < n_blocks; i++)
//    {
//      this->vmult_add(i, i, dst.block(i), src.block(i));
//      dst.block(i) /= velocities[i];
//    }
//  }
/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  double FisionMatrixBase<dim, n_fe_degree>::vmult_dot (ParallelBlockVector &dst,
    const ParallelBlockVector &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);
    ParallelBlockVector inter(src);

    vmult(inter, src);
    double val = dst * inter;

    return val;

  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
//template <int dim, int n_fe_degree>
//  double FisionMatrixBase<dim, n_fe_degree>::vmult_dot_velocities (
//    ParallelBlockVector &dst,
//    const ParallelBlockVector &src) const
//  {
//    AssertDimension(dst.n_blocks(), n_blocks);
//
//    ParallelBlockVector inter(src);
//    vmult(inter, src);
//    for (unsigned int nb = 0; nb < n_blocks; nb++)
//      inter.block(nb) /= velocities[nb];
//
//    double val = dst * inter;
//    return val;
//
//  }
// ----------------------------------------------------------------- //
// --------- ParallelBlockVector Complete Multiplications  --------- //
// ----------------------------------------------------------------- //
/**
 * @brief Complete matrix-vector multiplication.
 * dst = FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_transpose (
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
 * dst = dst + FisionMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add_transpose (
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
  void FisionMatrixBase<dim, n_fe_degree>::vmult (
    const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf_blocks[row][col]->vmult(dst, src);
    else
      matrix_blocks[row][col]->vmult(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add (
    const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf_blocks[row][col]->vmult_add(dst, src);
    else
      matrix_blocks[row][col]->vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add_row (
    const unsigned int row,
    const unsigned int col,
    double &dst,
    const PETScWrappers::MPI::Vector &src,
    unsigned int dst_row_on_block) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    //if (matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree)
    //  mass_mf_blocks[row][col]->vmult_add_row(dst, src, dst_row);
    //else FIXME
    //matrix_blocks[row][col]->vmult_add_row(dst, src, dst_row);
    AssertRelease(matrixfree_type == full_allocated,
      "LUPOD ONLY IMPLEMENTED FOR full_allocated");

    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    MatGetRow((*matrix_blocks[row][col]), dst_row_on_block, &ncols, &cols, &vals);
    for (PetscInt j = 0; j < ncols; j++)
      dst += src[cols[j]] * vals[j];

    MatRestoreRow((*matrix_blocks[row][col]), dst_row_on_block, &ncols, &cols, &vals);
  }

/**
 * v1 = v2 + A * v3
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add (
    const unsigned int row,
    const unsigned int col,
    PETScWrappers::MPI::Vector &v1,
    PETScWrappers::MPI::Vector &v2,
    const PETScWrappers::MPI::Vector &v3) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    v1 = v2;

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf_blocks[row][col]->vmult_add(v1, v3);
    else
      matrix_blocks[row][col]->vmult_add(v1, v3);
  }

// ----------------------------------------------------------------- //
// -------------------- Block Multiplications ---------------------- //
// ------------------------ ParallelVector ------------------------- //
// ----------------------------------------------------------------- //
/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixBase<dim, n_fe_degree>::vmult_add (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf_blocks[row][col]->vmult_add(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      matrix_blocks[row][col]->vmult_add(dst_vec, src_vec);
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
  void FisionMatrixBase<dim, n_fe_degree>::vmult (const unsigned int row,
    const unsigned int col,
    ParallelVector &dst,
    const ParallelVector &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    if (matrixfree_type == non_diagonal
        or matrixfree_type == full_matrixfree)
      mass_mf_blocks[row][col]->vmult(dst, src);
    else
    {
      PETScWrappers::MPI::Vector dst_vec(comm, dst.size(), dst.local_size());
      PETScWrappers::MPI::Vector src_vec(comm, src.size(), src.local_size());
      copy_to_Vector(src_vec, src);
      matrix_blocks[row][col]->vmult_add(dst_vec, src_vec);
      //dst = dst_vec;
      LinearAlgebra::ReadWriteVector<double> rwv(dst_vec.locally_owned_elements());
      rwv.import(dst_vec, VectorOperation::insert);
      dst.import(rwv, VectorOperation::insert);
    }

  }

// ----------- Explicit Instantations ----------- //

template class TransportMatrixBase<1, 1> ;
template class TransportMatrixBase<1, 2> ;
template class TransportMatrixBase<1, 3> ;
template class TransportMatrixBase<1, 4> ;
template class TransportMatrixBase<1, 5> ;

template class TransportMatrixBase<2, 1> ;
template class TransportMatrixBase<2, 2> ;
template class TransportMatrixBase<2, 3> ;
template class TransportMatrixBase<2, 4> ;
template class TransportMatrixBase<2, 5> ;

template class TransportMatrixBase<3, 1> ;
template class TransportMatrixBase<3, 2> ;
template class TransportMatrixBase<3, 3> ;
template class TransportMatrixBase<3, 4> ;
template class TransportMatrixBase<3, 5> ;

template class FisionMatrixBase<1, 1> ;
template class FisionMatrixBase<1, 2> ;
template class FisionMatrixBase<1, 3> ;
template class FisionMatrixBase<1, 4> ;
template class FisionMatrixBase<1, 5> ;

template class FisionMatrixBase<2, 1> ;
template class FisionMatrixBase<2, 2> ;
template class FisionMatrixBase<2, 3> ;
template class FisionMatrixBase<2, 4> ;
template class FisionMatrixBase<2, 5> ;

template class FisionMatrixBase<3, 1> ;
template class FisionMatrixBase<3, 2> ;
template class FisionMatrixBase<3, 3> ;
template class FisionMatrixBase<3, 4> ;
template class FisionMatrixBase<3, 5> ;

