/**
 * @file   matrix_operators_small_time_mat.h
 * @brief
 */
#ifndef MATRIX_OPERATORS_SMALL_TIME_MAT_H_
#define MATRIX_OPERATORS_SMALL_TIME_MAT_H_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_base.h>

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

#include "../materials.h"
#include "matrix_operators_free.h"
#include "matrix_operators_base.h"

using namespace dealii;

typedef LinearAlgebra::distributed::Vector<double, MemorySpace::Host> ParallelVector;
typedef LinearAlgebra::distributed::BlockVector<double> ParallelBlockVector;

/**
 *
 */
template <int dim, int n_fe_degree>
  class SmallDelayedFissionMatrices : public Subscriptor
  {

    public:

    /**
     * @brief Constructor associate the object to a DoFHandler and a AffineConstraints<double>.
     * It stores references to this objects
     */
    SmallDelayedFissionMatrices (const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     *
     */
    void reinit (const Materials &materials,
      const MatrixFreeType &matrix_free_type = non_diagonal,
      bool listen_to_material_id = false);

    /**
     *
     */
    std::size_t memory_consumption () const;

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
     * dst = dst + SmallMatrixBase * src
     */
    void vmult_add (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Complete Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = SmallMatrixBase * src
     */
    void vmult (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + SmallMatrixBase * src
     */
    void vmult_add (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */

    void vmult_add (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * v1 = v2 + A * v3
     */
    void vmult_add (const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &v1,
      PETScWrappers::MPI::Vector &v2,
      const PETScWrappers::MPI::Vector &v3) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------------ ParallelVector ------------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */

    void vmult_add (const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    void vmult (const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    dealii::MatrixFree<dim, double> matfree_data;

    MPI_Comm comm;

    // Easy access to problem structures
    const Triangulation<dim> &tria;
    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    std::vector<types::global_dof_index> local_dofs_per_process;

    // Members
    MatrixFreeType matrixfree_type;
    unsigned int n_prec_groups;
    unsigned int n_energy_groups;
    unsigned int n_dofs;
    bool listen_to_material_id;

    // Allocate Matrices Structures
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*>> matrix_csr;
    SparsityPattern sp;
    // Matrix Free Structures
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*>> mass_mf;

    private:

    /**
     *
     */
    void reinit_full_matrixfree (const Materials &materials);

    /**
     *
     */
    void assemble_full_matrices (const Materials &materials);

    std::vector<std::vector<Vector<double> > > coeffs;

  };

#endif /* MATRIX_OPERATORS_SMALL_TIME_MAT_H_ */
