/**
 * @file   matrix_operators_noise_diff.h
 * @brief  Class to handle complex block matrices.
 */

#ifndef MATFREE_OPERATORS_NOISE_DIFF_H_
#define MATFREE_OPERATORS_NOISE_DIFF_H_

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_base.h>

#include <deal.II/distributed/shared_tria.h>

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
#include "../complex_perturbation.h"
#include "matrix_operators_base.h"
#include "matrix_operators_complex_base.h"
#include "matrix_operators_free.h"
#include "matrix_operators_noise_diff.h"

using namespace dealii;
typedef LinearAlgebra::distributed::Vector<double> ParallelVector;

/**
 * @brief Calculate the precursors factor depending on a materials
 */
void calculate_precursors_factor (
  const unsigned int mat_id,
  const Materials &materials,
  const double &omega,
  std::vector<complex> &prec_factor);

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

template <int dim, int n_fe_degree>
  class NoiseAMatrix : public TransportMatrixComplexBase<dim, n_fe_degree>
  {
    public:

    /**
     * @brief Constructor
     */
    NoiseAMatrix (
      const MPI_Comm &_comm,
      const DoFHandler<dim> &dh,
      const AffineConstraints<double> &constraints);

    /**
     *
     */
    void reinit (
      const Materials &materials,
      const ComplexPerturbation &pert,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors,
      const MatrixFreeType &matrix_free_type = non_diagonal);

    /**
     *
     */
    std::size_t memory_consumption () const;

    void clear ();

    dealii::MatrixFree<dim, double> matfree_data;

    std::vector<std::vector<Vector<double> > > coeffs_real;
    std::vector<std::vector<Vector<double> > > coeffs_imag;

    std::vector<Vector<double> > coeffs_grad;
    //std::vector<std::vector<Vector<double> > > coeffs_val;
    //std::vector<Vector<double> >  coeffs_bound;

    // Easy access to problem structures
    const Triangulation<dim> &tria;
    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;

    private:
    /**
     *
     */
    void reinit_non_diagonal (
      const Materials &materials,
      const ComplexPerturbation &pert,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     * @Calculate A factor
     */
    void calculate_A_factor (
      const unsigned int mat_id,
      const Materials &materials,
      const ComplexPerturbation &pert,
      std::vector<std::vector<complex> > &B_factor);

    /**
     *
     */
    void reinit_full_matrixfree (
      const Materials &materials,
      const ComplexPerturbation &pert,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     *
     */
    void assemble_full_matrices (
      const Materials &materials,
      const ComplexPerturbation &pert,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

    /**
     *
     */
    void assemble_diagonal_matrices (
      const Materials &materials,
      const ComplexPerturbation &pert,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors);

  };

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template <int dim, int n_fe_degree>
  class NoiseBMatrix : public MassMatrixComplexBase<dim, n_fe_degree>
  {
    public:

    /**
     * @brief Constructor associate the object to a DoFHandler and a AffineConstraints<double>.
     * It stores references to this objects
     */
    NoiseBMatrix (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     * @Calculate B factor
     */
    void calculate_B_factor (
      const unsigned int mat_id,
      const unsigned int pert_id,
      const Materials &materials,
      const ComplexPerturbation &pert,
      std::vector<std::vector<complex> > &B_factor);

    /**
     * @brief Allocate and assemble the matrix associated to this materials and fission
     * cross sections.
     */
    void reinit (const Materials &materials,
      const ComplexPerturbation &pert,
      const MatrixFreeType &matrix_free_type = non_diagonal);

    /**
     * @brief
     */
    std::size_t memory_consumption () const;

    void clear ();

    dealii::MatrixFree<dim, double> matfree_data;

    private:

    /**
     * @brief
     */
    void assemble_full_rhs_cell_wise (
      const Materials &materials,
      const ComplexPerturbation &pert);

    /**
     *
     */
    void assemble_full_rhs_borders (
      const Materials &materials,
      const ComplexPerturbation &pert);

    /**
     *
     */
    void assemble_full_rhs_bordershex (
      const Materials &materials,
      const ComplexPerturbation &pert);

    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;

    std::vector<std::vector<Vector<double> > > coeffs_real;
    std::vector<std::vector<Vector<double> > > coeffs_imag;
  };
#endif /* MATFREE_OPERATORS_NOISE_DIFF_H_ */
