/**
 * @file   matrix_operators_complex_base.h
 * @brief  Class to handle complex block matrices.
 */

#ifndef MATRIX_OPERATORS_COMPLEX_BASE_H_
#define MATRIX_OPERATORS_COMPLEX_BASE_H_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_data.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/diagonal_matrix.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "../io/materials.h"
#include "matrix_operators_free.h"
#include "matrix_operators_base.h"

using namespace dealii;

typedef LinearAlgebra::distributed::Vector<double, MemorySpace::Host> ParallelVector;
typedef LinearAlgebra::distributed::BlockVector<double> ParallelBlockVector;

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

template <int dim, int n_fe_degree>
  class TransportMatrixComplexBase : public Subscriptor
  {

    public:

    /**
     * @brief Constructor
     */
    TransportMatrixComplexBase (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dh,
      const AffineConstraints<double> &constraints);

    /**
     *
     */
    void clear ();

    /**
     * @brief Return the total number of columns in this matrix.
     */
    unsigned int n () const;

    /**
     * @brief  Return the total number of rows in this matrix.
     */
    unsigned int m () const;

    /**
     * @brief  Return the number of row blocks in a column.
     */
    unsigned int n_blocks_rows () const;

    /**
     * @brief  Return the number of column blocks in a row.
     */
    unsigned int n_blocks_cols () const;

    /**
     * @brief Return the number of dofs per block.
     */
    unsigned int n_dofs_blocks () const;

    /**
     * @brief Return the number of dofs per block.
     */
    unsigned int get_n_moments () const;

    /**
     *
     */
    PETScWrappers::MPI::SparseMatrix& block (
      unsigned int i,
      unsigned int j);

    /**
     *
     */
    void get_inv_diagonal (const unsigned int group,
      PETScWrappers::MPI::Vector &inv_diagonal);

    /*
     *
     */
    void get_inv_diagonal (
      PETScWrappers::MPI::BlockVector &inv_mat);
    /**
     *
     */
    void get_inv_diagonal (const unsigned int group,
      DiagonalMatrix<PETScWrappers::MPI::Vector> &inv_diagonal);

    /**
     *
     */
    bool is_block_real (
      const unsigned int row,
      const unsigned int col) const;

    /**
     *
     */
    bool is_block_imag_negative (
      const unsigned int row,
      const unsigned int col) const;

    // ----------------------------------------------------------------- //
    // --- PETScWrappers::MPI::BlockVector Complete Multiplications  --- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication. k3 real formulation.
     * dst = TransportMatrix * src
     */
    void vmult_k3 (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication. k3 real formulation.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add_k3 (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --- PETScWrappers::MPI::BlockVector Real and Imaginary parts  --- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrixComplexBase * src
     */
    void vmult_real (
      PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrixComplexBase * src
     */
    void vmult_add_real (
      PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrixComplexBase * src
     */
    void vmult_imag (
      PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrixComplexBase * src
     */
    void vmult_add_imag (
      PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // ----------- Block Real and Imag Multiplications ----------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult_real (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add_real (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_imag (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add_imag (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;


    // ----------------------------------------------------------------- //
    // ----------- Block Real and Imag Multiplications ----------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult_group (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     *
     */
    void vmult_add_group (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;


    // ----------------------------------------------------------------- //
    // --- PETScWrappers::MPI::BlockVector Adjoint Multiplications  --- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    void Tvmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * dst
     */
    double vmult_dot (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Complete Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Transpose Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * dst
     */
    double vmult_dot (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * v1 = v2 + A * v3
     */
//    void vmult_add (
//      const unsigned int row,
//      const unsigned int col,
//      PETScWrappers::MPI::Vector &v1,
//      PETScWrappers::MPI::Vector &v2,
//      const PETScWrappers::MPI::Vector &v3) const;
    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ----------------------------- Vec ------------------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      Vec &dst,
      const Vec &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      Vec &dst,
      const Vec &src) const;

//    /**
//     * v1 = v2 + A * v3
//     */
//    void vmult_add (
//      const unsigned int row,
//      const unsigned int col,
//      Vec &v1,
//      Vec &v2,
//      const Vec &v3) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------------ Vector<double> ------------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    void print_matrix_in_matlab (
      unsigned int block_i,
      unsigned int block_j,
      const std::string &mat_name,
      std::ostream &out,
      const unsigned int precision);

    // ----------------------------------------------------------------- //
    // --------------------  k3 Formultaion       ---------------------- //
    // ------------------------ Vector<double> ------------------------- //
    // ----------------------------------------------------------------- //

    void vmult_add_k3 (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    bool is_block_imag_k3 (
      const unsigned int row,
      const unsigned int col) const;

    bool is_block_negative_real_k3 (
      const unsigned int row,
      const unsigned int col) const;

    MPI_Comm comm;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    MatrixFreeType matrixfree_type;

    unsigned int n_dofs_block;
    unsigned int n_blocks;
    unsigned int n_blocks_real;
    unsigned int n_moments;

    const DoFHandler<dim> &dof_handler;

    std::vector<unsigned int> boundary_conditions;
    std::vector<double> albedo_factors;
    // double delta_t;
    // std::string type_scheme;

    protected:

    // Matrices Structures
    SparsityPattern sp;
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*> > matrix_blocks_real;
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*> > matrix_blocks_imag;
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*> > mass_mf_blocks_real;
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*> > mass_mf_blocks_imag;

    std::vector<PoissonOperator<dim, n_fe_degree, double>*> poison_mf_blocks_real;
  };

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template <int dim, int n_fe_degree>
  class MassMatrixComplexBase : public Subscriptor
  {

    public:

    /**
     * @brief Constructor associate the object to a DoFHandler and a AffineConstraints.
     * It stores references to this objects
     */
    MassMatrixComplexBase (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     * @brief Clear the matrix.
     */
    void clear ();

    /**
     * @brief Return the total number of columns in this matrix.
     */
    unsigned int n () const;

    /**
     * @Brief Return the total number of rows in this matrix.
     */
    unsigned int m () const;

    /**
     * @Brief Return the number of row blocks in a column.
     */
    unsigned int n_blocks_rows () const;

    /**
     * @Brief Return the number of column blocks in a row.
     */
    unsigned int n_blocks_cols () const;

    /**
     * @brief Get a  constant reference  SparseMatrix<double> to the block matrix.
     */
    PETScWrappers::MPI::SparseMatrix& block (unsigned int i,
      unsigned int j);

    /**
     *
     */
    bool is_block_real (
      const unsigned int row,
      const unsigned int col) const;

    /**
     *
     */
    bool is_block_imag_negative (
      const unsigned int row,
      const unsigned int col) const;

    // ----------------------------------------------------------------- //
    // ------- PETScWrappers::MPI::BlockVector Multiplications  -------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * A * dst
     */
    double vmult_dot (const PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --PETScWrappers::MPI::BlockVector Transpose Multiplications  ---- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Complete Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Transpose Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //

    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * v1 = v2 + A * v3
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &v1,
      PETScWrappers::MPI::Vector &v2,
      const PETScWrappers::MPI::Vector &v3) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------------ Vector<double> ------------------------- //
    // ----------------------------------------------------------------- //

    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    MPI_Comm comm;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    // Members
    MatrixFreeType matrixfree_type;
    unsigned int n_blocks;
    unsigned int n_blocks_real;
    unsigned int n_dofs_block;

    // Allocate Matrices Structures
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*> > matrix_blocks_real;
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*> > matrix_blocks_imag;
    SparsityPattern sp;
    // Matrix Free Structures
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*> > mass_mf_blocks_real;
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*> > mass_mf_blocks_imag;
  };

#endif /* MATRIX_OPERATORS_COMPLEX_BASE_H_ */
