/**
 * @file   matrix_operators_petsc_time.cc
 * @brief  Implementation of SystemMatrixTime to handle block matrices.
 */


#ifndef MATFREE_OPERATORS_TIME_H_
#define MATFREE_OPERATORS_TIME_H_

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

#include "../io/materials.h"
#include "matrix_operators_base.h"
#include "matrix_operators_free.h"

using namespace dealii;
typedef LinearAlgebra::distributed::Vector<double> ParallelVector;

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

template<int dim, int n_fe_degree>
class SystemMatrixTime: public TransportMatrixBase<dim, n_fe_degree> {
public:

	/**
	 * @brief Constructor
	 */
	SystemMatrixTime(const MPI_Comm &_comm, const DoFHandler<dim> &dh,
			const AffineConstraints<double> &constraints);

	/**
	 *
	 */
	void reinit(const Materials &materials,
			const std::vector<unsigned int> &boundary_conditions,
			const std::vector<double> &albedo_factors, double deltat,
			const std::string &_scheme_type,
			const MatrixFreeType &matrix_free_type = non_diagonal,
			bool listen_to_material_id = false);

	/**
	 *
	 */
	std::size_t memory_consumption() const;

	dealii::MatrixFree<dim, double> matfree_data;

	// Easy access to problem structures
	const Triangulation<dim> &tria;
	const DoFHandler<dim> &dof_handler;
	const AffineConstraints<double> &constraints;

//	std::vector<unsigned int > boundary_conditions;
//	std::vector<double> albedo_factors;
//
//	double delta_t;
//
//	std::string type_scheme;

private:
	/**
	 *
	 */
	void reinit_non_diagonal(const Materials &materials,
			bool listen_to_material_id,
			const std::vector<unsigned int> &boundary_conditions,
			const std::vector<double> &albedo_factors);

	/**
	 *
	 */
	void reinit_full_matrixfree(const Materials &materials,
			bool listen_to_material_id,
			const std::vector<unsigned int> &boundary_conditions,
			const std::vector<double> &albedo_factors);

	/**
	 *
	 */
	void assemble_full_matrices(const Materials &materials,
			const std::vector<unsigned int> &boundary_conditions,
			const std::vector<double> &albedo_factors);

	/**
	 *
	 */
	void assemble_diagonal_matrices(const Materials &materials,
			const std::vector<unsigned int> &boundary_conditions,
			const std::vector<double> &albedo_factors);

	std::vector<std::vector<Vector<double> > > coeffs;
    std::vector<Vector<double>  > coeffs_grad;
    std::vector<Vector<double>  > coeffs_val;



};

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template<int dim, int n_fe_degree>
class FisionMatrixTime: public FisionMatrixBase<dim, n_fe_degree> {
public:

	/**
	 * @brief Constructor associate the object to a DoFHandler and a AffineConstraints<double>.
	 * It stores references to this objects
	 */
	FisionMatrixTime(const MPI_Comm &comm, const DoFHandler<dim> &dof_handler,
			const AffineConstraints<double> &constraints);

	/**
	 * @brief Allocate and assemble the matrix associated to this materials and fission
	 * cross sections.
	 */
	void reinit(const Materials &materials,
			const MatrixFreeType &matrix_free_type = non_diagonal,
			bool listen_to_material_id = false);

	/**
	 * @brief
	 */
	std::size_t memory_consumption() const;

	dealii::MatrixFree<dim, double> matfree_data;

private:

	/**
	 * @brief
	 */
	void assemble_full_matrices(const Materials &materials);

	const DoFHandler<dim> &dof_handler;
	const AffineConstraints<double> &constraints;

	std::vector<std::vector<Vector<double> > > coeffs;

};

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template<int dim, int n_fe_degree>
class MassMatrixTime: public FisionMatrixBase<dim, n_fe_degree> {
public:

	/**
	 * @brief Constructor associate the object to a DoFHandler and a AffineConstraints<double>.
	 * It stores references to this objects
	 */
	MassMatrixTime(const MPI_Comm &comm, const DoFHandler<dim> &dof_handler,
			const AffineConstraints<double> &constraints);

	/**
	 * @brief Allocate and assemble the matrix associated to this materials and fission
	 * cross sections.
	 */
	void reinit(const Materials &materials,
			const MatrixFreeType &matrix_free_type = non_diagonal,
			bool listen_to_material_id = false);

	/**
	 * @brief
	 */
	std::size_t memory_consumption() const;

	dealii::MatrixFree<dim, double> matfree_data;

private:

	/**
	 * @brief
	 */
	void assemble_full_matrices(const Materials &materials);

	const DoFHandler<dim> &dof_handler;
	const AffineConstraints<double> &constraints;

	std::vector<std::vector<Vector<double> > > coeffs;

};

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//



#endif /* EPS_MATFREE_SOLVER_H_ */
