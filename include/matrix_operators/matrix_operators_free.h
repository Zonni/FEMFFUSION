/**
 * @file   matrix_operators.h
 * @brief  Class to handle block matrices.
 */

#ifndef MATFREE_OPERATORS_FREE_H_
#define MATFREE_OPERATORS_FREE_H_

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/operators.h>

#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "../io/materials.h"
#include "../femffusion.h"

using namespace dealii;

typedef LinearAlgebra::distributed::Vector<double> ParallelVector;

// TODO Make a common base class that initializes all the common parameters that use the
// block-matrices
template <int dim, int n_fe_degree, typename number>
  class MassOperator : public Subscriptor
  {
    public:

    /**
     * @brief Constructor
     */
    MassOperator (dealii::MatrixFree<dim, number> &matfree_data);

    /**
     *
     */
    void clear ();

    /**
     * @brief Reinit and Assemble the Mass matrix.
     * This version differs from the previous one because it does not implement
     * boundary conditions.
     */
    void reinit (
      const AffineConstraints<double> &constraints,
      const std::vector<unsigned int> &materials,
      const Vector<double> &alpha,
      const bool listen_to_material_id);

    /**
     * @brief
     */
    void reinit (
      const AffineConstraints<double> &_contraints,
      const std::vector<unsigned int> &_materials,
      const Vector<double> &_alpha,
      const std::vector<unsigned int> &_boundary_conditions,
      const double &_bc_coeff,
      const bool _listen_to_material_id);

    /**
     * @brief
     */
    unsigned int size () const;

    /**
     * @brief Return the number of columns in this matrix.
     */
    unsigned int n () const;

    /**
     * @brief Return the number of rows in this matrix.
     */
    unsigned int m () const;

    /**
     * @brief Vector Matrix multiplication. dst = Mat * src
     * The dst vector is deleted before the multiplication.
     */
    void vmult (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief Vector Matrix multiplication and addition of the previous results. dst += Mat * src
     */
    void vmult_add (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief
     */
    void Tvmult (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief
     */
    void Tvmult_add (ParallelVector &dst,
      const ParallelVector &src) const;

    // The same but for PETSc Vectors
    /**
     * @brief
     */
    void vmult (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * @brief
     */
    void Tvmult (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * @brief
     */
    void vmult_add (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * @brief
     */
    void vmult_add_row (
      double &dst,
      const PETScWrappers::MPI::Vector &src,
      unsigned int row_dst);

    /**
     * @brief
     */
    void Tvmult_add (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    // The same but for (PETSc) Vec
    /**
     * @brief
     */
    void vmult (Vec &dst,
      const Vec &src) const;

    /**
     * @brief
     */
    void Tvmult (Vec &dst,
      const Vec &src) const;

    /**
     * @brief
     */
    void vmult_add (Vec &dst,
      const Vec &src) const;

    /**
     * @brief
     */
    void Tvmult_add (Vec &dst,
      const Vec &src) const;

    private:
    /**
     * @brief
     */
    void cell_local_apply (const dealii::MatrixFree<dim, number> &data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    /**
     * @brief
     */
    void cell_local_apply_row (
      const dealii::MatrixFree<dim, number> &matfree_data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    void cell_boundary_apply (const dealii::MatrixFree<dim, number> &matfree_data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &face_range) const;

    void cell_face_apply (const dealii::MatrixFree<dim, number> &matfree_data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &face_range) const;

    dealii::MatrixFree<dim, number> &matfree_data;
    const std::vector<unsigned int> *materials;
    const Vector<double> *alpha;

    // Easy access to problem structures
    const AffineConstraints<double> *constraints;

    bool listen_to_material_id;
    bool zero_matrix; // Is a zero matrix?

    // Boundary
    bool boundary_req;
    const std::vector<unsigned int> *boundary_conditions;
    double bc_coeff;

    // ROM - LUPOD
    unsigned int apply_on_row_index;

  };

// TODO Make a common base class that initializes all the common parameters that use the
// block-matrices
template <int dim, int n_fe_degree, typename number>
  class PoissonOperator : public MatrixFreeOperators::Base<dim, ParallelVector>
  {
    public:

    /**
     * @brief Constructor
     */
    PoissonOperator (dealii::MatrixFree<dim, number> &matfree_data);

    /**
     *
     */
    void clear ();

    /**
     * @brief
     */
    void reinit (
      const unsigned int _group,
      const AffineConstraints<double> &_contraints,
      const Materials &_materials,
      const std::vector<unsigned int> &_materials_vec,
      const Vector<double> &_coef_val,
      const Vector<double> &_coef_grad,
      const std::vector<unsigned int> &_boundary_conditions,
      const std::vector<double> &_albedo_factors,
      const bool _alpha_modes = false,
      const double bc_factor = 0.5);

    unsigned int size () const;

    /**
     * @brief Return the number of columns in this matrix.
     */
    unsigned int n () const;

    /**
     * @brief Return the number of rows in this matrix.
     */
    unsigned int m () const;

    /**
     * @brief Vector Matrix multiplication. dst = Mat * src
     * The dst vector is deleted before the multiplication.
     */
    void vmult (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief Vector Matrix multiplication and addition of the previous results. dst += Mat * src
     */
    void vmult_add (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief Vector Matrix multiplication only 1 row
     */
    void vmult_add_row (
      double &dst,
      const PETScWrappers::MPI::Vector &src,
      unsigned int row_dst);

    /**
     * @brief
     */
    void Tvmult (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief
     */
    void Tvmult_add (ParallelVector &dst,
      const ParallelVector &src) const;

    // The same but for PETSc Vectors
    /**
     * @brief
     */
    void vmult (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * @brief
     */
    void Tvmult (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * @brief
     */
    void vmult_add (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * @brief
     */
    void Tvmult_add (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    // The same but for (PETSc) Vec
    /**
     * @brief
     */
    void vmult (Vec &dst,
      const Vec &src) const;

    /**
     * @brief
     */
    void Tvmult (Vec &dst,
      const Vec &src) const;

    /**
     * @brief
     */
    void vmult_add (Vec &dst,
      const Vec &src) const;

    /**
     * @brief
     */
    void Tvmult_add (Vec &dst,
      const Vec &src) const;

    /**
     *
     */
    void compute_diagonal ();

    ParallelVector inverse_diagonal;

    private:

    virtual void apply_add (ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     * @brief
     */
    void cell_local_apply (const dealii::MatrixFree<dim, number> &data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    /**
     *
     */
    void cell_local_apply_row (
      const dealii::MatrixFree<dim, number> &matfree_data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    void cell_boundary_apply (const dealii::MatrixFree<dim, number> &matfree_data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &face_range) const;

    void cell_face_apply (const dealii::MatrixFree<dim, number> &matfree_data,
      ParallelVector &dst,
      const ParallelVector &src,
      const std::pair<unsigned int, unsigned int> &face_range) const;

    /**
     *
     */
    void local_compute_diagonal (const MatrixFree<dim, number> &data,
      ParallelVector &dst,
      const unsigned int&,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    dealii::MatrixFree<dim, number> &matfree_data;
    const std::vector<unsigned int> *materials_vector;
    const Vector<double> *coef_val, *coef_grad;

    bool boundary_req;
    unsigned int group;
    bool alpha_modes;
    double bc_factor;

    // Materials
    const Materials *materials;
    // Boundary
    const std::vector<unsigned int> *boundary_conditions;
    const std::vector<double> *albedo_factors;
    const AffineConstraints<double> *constraints;

    unsigned int apply_on_row_index;

  };

#endif /* MATFREE_OPERATORS_FREE_H_ */
