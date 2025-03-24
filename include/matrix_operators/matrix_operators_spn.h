/**
 *
 * @file   matrix_operators.h
 * @brief  Class to handle block matrices.
 *
 */

#ifndef MATFREE_OPERATORS_SPN_H_
#define MATFREE_OPERATORS_SPN_H_

#include <deal.II/base/quadrature_lib.h>

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

static const unsigned int max_mom = 4;

static const double phi_to_u_coeff[max_mom][max_mom] =
                                                         {
                                                             { 1., 2., 0., 0. },
                                                             { 0., 3., 4., 0. },
                                                             { 0., 0., 5., 6. },
                                                             { 0., 0., 0., 7. }
                                                         };

static const double u_to_phi_coeff[max_mom][max_mom] =
      {
          { 1, -2. / 3, 8. / 15, -16.
                                 / 35 },
          { 0, 1. / 3, -4. / 15, +8.
                                 / 35 },
          { 0, 0, 1. / 5, -6. / 35 },
          { 0, 0, 0, +1. / 7 }
      };

static const double sp_coeff[max_mom][max_mom][max_mom] =
      {
          {
              { 1., -2. / 3, 8. / 15, -16. / 35 },
              { -2. / 3, 4. / 9, -16. / 45, 32. / 105 },
              { 8. / 15, -16. / 45, 64. / 225, -128. / 525 },
              { -16. / 35, 32. / 105, -128. / 525, 256. / 1225 },
          },
          {
              { 0.0, 0.0, 0.0, 0.0 },
              { 0.0, 5. / 9, -4. / 9, 8. / 21 },
              { 0.0, -4. / 9, 16. / 45, -32. / 105 },
              { 0.0, 8. / 21, -32. / 105, 64. / 245 },
          },
          {
              { 0.0, 0.0, 0.0, 0.0 },
              { 0.0, 0.0, 0.0, 0.0 },
              { 0.0, 0.0, 9. / 25, -54. / 175 },
              { 0.0, 0.0, -54. / 175, 324. / 1225 },
          },
          {
              { 0.0, 0.0, 0.0, 0.0 },
              { 0.0, 0.0, 0.0, 0.0 },
              { 0.0, 0.0, 0.0, 0.0 },
              { 0.0, 0.0, 0.0, 13. / 49 },
          }

      };

static const double marshack_coeff[max_mom][max_mom] =
                                                         {
                                                             { +1. / 2, -1. / 8, 1. / 16,  -5. / 128 },
                                                             { -1. / 8, +7. / 24, -41./ 384,     1. / 16 },
                                                             { +1. / 16, -41. / 384, 407./ 1920, -233. / 2560 },
                                                             { -5. / 128, 1. / 16, -233./ 2560,   3023. / 17920 }
                                                         };

static const double diff_coeff[max_mom] =
                                            { 1. / 3, 1. / 7, 1. / 11., 1. / 15 };


// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

typedef LinearAlgebra::distributed::Vector<double> ParallelVector;

template <int dim, int n_fe_degree>
  class TransportMatrixSPN : public TransportMatrixBase<dim, n_fe_degree>
  {

    public:

    /**
     * @brief Constructor
     */
    TransportMatrixSPN (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     *
     */
    void reinit (
      const Materials &materials,
      const unsigned int spn,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors,
      const MatrixFreeType &matrix_free = non_diagonal,
      bool listen_to_material_id = false);

    /**
     *
     */
    void reinit_non_diagonal (
      const Materials &materials,
      bool listen_to_material_id,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     *
     */
    void reinit_full_matrixfree (
      const Materials &materials,
      bool listen_to_material_id,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     *
     */
    std::size_t memory_consumption () const;

    std::vector<Vector<double> > coeffs_grad;
    std::vector<Vector<double> > coeffs_val;
    std::vector<double> coeffs_bound;

    protected:

    /**
     *
     */
    void assemble_full_matrices (
      const Materials &materials,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     *
     */
    void assemble_diagonal_matrices (
      const Materials &materials,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     *
     */
    void make_boundary_sp (
      const DoFHandler<dim> &dof_handler,
      DynamicSparsityPattern &sparsity);

    /**
     *
     */
    unsigned int gm_to_b (const unsigned int group,
      const unsigned int moment) const;

    unsigned int n_moments;
    unsigned int n_groups;

    // Easy access to problem structures
    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;

    // Matrix Free Structures
    std::vector<std::vector<std::vector<std::vector<Vector<double> > > > > coeffs_L;
    dealii::MatrixFree<dim, double> matfree_data;
  };

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template <int dim, int n_fe_degree>
  class FisionMatrixSPN : public FisionMatrixBase<dim, n_fe_degree>
  {
    public:

    /**
     * @brief Constructor associate the object to a DoFHandler and a AffineConstraints<double>.
     * It stores references to this objects
     */
    FisionMatrixSPN (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     * @brief Allocate and assemble the matrix associated to this materials and fission
     * cross sections.
     */
    void reinit (const Materials &materials,
      const unsigned int spn,
      const MatrixFreeType &matrix_free = non_diagonal,
      bool listen_to_material_id = false);

    /**
     *
     *
     */
    std::size_t memory_consumption () const;

    std::vector<Vector<double> > coeffs_grad;
    std::vector<Vector<double> > coeffs_val;
    std::vector<Vector<double> > coeffs_bound;

    private:

    /**
     *
     */
    void assemble_full_matrices (const Materials &materials);

    unsigned int gm_to_b (const unsigned int group,
      const unsigned int moment) const;

    unsigned int n_moments;

    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;

    std::vector<std::vector<std::vector<std::vector<Vector<double> > > > > coeffs_M;
    dealii::MatrixFree<dim, double> matfree_data;
  };

#endif /* MATFREE_OPERATORS_SPN_H_ */
