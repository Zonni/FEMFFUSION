/**
 * @file   eps_solver_2g.h
 * @brief  SolverEPS2G class. Solver optimized for two energy groups without up-scattering and
 * the energy fission spectrum all for the fast flux. Thus, the eigenvalues problem is yields,
 *  \begin{equation}
 *    \begin{pmatrix}
 *      \bf L_{11} & 0 \\  -\bf L_{21} & \bf L_{22} \\
 *    \end{pmatrix} \binom{\boldsymbol{\tilde{ \phi}_1}}{\boldsymbol{\tilde{ \phi}_2}}
 *    = \frac{1}{\lambda}
 *        \begin{pmatrix} \bf M_{11} & \bf M_{12} \\ 0 & 0\\ \end{pmatrix}
 *        \binom{\boldsymbol{\tilde{ \phi}_1}}{\boldsymbol{\tilde{ \phi}_2}} \,.
 * \end{equation}
 * That can be solved as
 *  \bf L_{11}^{\rm -1} \paren{M_{11} + M_{12} L_{22} ^{\rm -1} L_{21}} \boldsymbol{\tilde{ \phi}_1} = {\rm \lambda } \boldsymbol{\tilde{ \phi}_1} \,,
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/slepc_solver.h>


#include <slepceps.h>
#include <petscksp.h>
#include "../matrix_operators/matrix_operators_base.h"
#include "../matrix_operators/matrix_operators_petsc.h"


#ifndef EPSSOLVER_2G_H_
#define EPSSOLVER_2G_H_

using namespace dealii;

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template<int dim, int n_fe_degree>
void shell_vmult_2g(Mat PescMat, Vec src, Vec dst);

template<int dim, int n_fe_degree>
void shell_vmult_2g_adj(Mat PescMat, Vec src, Vec dst);


/**
 *   @brief This Class solves the eigenvalue problem associated with the two groups
 *   diffusion equation called the Lambda modes problem of a nuclear reactor.
 *   It uses PETSC and SLEPc in order to solve it.
 *   By default it uses Krylov-Schur eigenvalue solver from SLEPC with a shell Matrix.
 *   Every Matrix Vector multiplication is defined in the function MatxVec.
 *   This function needs to calculate multiplication by the inverse of global_A_11 and
 *   global_A_00. This is done with the PETCS routine of Conjugate Gradient
 *   and the Incomplete Cholesky preconditioner.
 */
template<int dim, int n_fe_degree>
class SolverEPS2G
{
public:

  /**
   * @brief Constructor of the EPS solver.
   */
  SolverEPS2G(
    TransportMatrixBase<dim, n_fe_degree> &L,
    FisionMatrixBase<dim, n_fe_degree> &M,
    unsigned int _n_eigenvalues);

  /**
   * @brief Destroy the eps and all related objects.
   */
  ~SolverEPS2G();

  /**
   * @brief Solve the eigenvalue problem.
   * It also computes the adjoint eigenvectors if they are requested at the constructor.
   */
  void solve(std::vector<double> &eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector> &phi);

  PetscErrorCode ksp_solve(KSP ksp, Vec b, Vec x);

  /**
   * @brief Get the mean number of iterations for the linear system related with group one,
   * <tt>its_1</tt>,  and  with group two, <tt>its_2</tt>.
   */
  void get_n_inner_iterations(unsigned int& its_1, unsigned int& its_2) const;

  /**
   * @brief Get the total number iterations for one-group linear systems.
   */
  unsigned int get_n_inner_iterations() const;

  /**
   * @brief Compute the flux related with the thermal flux, <tt>phi_1</tt>,
   * from the fast flux  <tt>phi_0</tt>.
   * It need to solve a linear system related with L_11 block.
   */
  PetscErrorCode compute_phi_1(
    const PETScWrappers::MPI::Vector &phi_0,
    PETScWrappers::MPI::Vector &phi_1);


  PetscErrorCode compute_phi_1_adj(double eigenvalue,
      const PETScWrappers::MPI::Vector &phi_0,
      PETScWrappers::MPI::Vector &phi_1);
  /**
   * @brief Computes the adjoint problem from the standard eigenvector.
   */
  void compute_adjoint_slepc (
    std::vector<double>& eigenvalues,
    std::vector<
        PETScWrappers::MPI::BlockVector>& phi_adj);

  void compute_adjoint_fix_point(
      std::vector<double>& eigenvalues,
	    std::vector<
	        PETScWrappers::MPI::BlockVector>& phi,
      std::vector<
          PETScWrappers::MPI::BlockVector>& phi_adj);

  PetscErrorCode apply_fix_point (
    double eigen,
    const PETScWrappers::MPI::Vector& phi0,
    PETScWrappers::MPI::Vector& phi0_adj,
    PETScWrappers::MPI::Vector& phi1_adj,
	double &norm);

  PetscErrorCode gram_schmidt_mod(
    std::vector<PETScWrappers::MPI::Vector>& phi1);

  /**
   * @brief Compute the actual absolute error norm of the eps problem.
   */
  PetscErrorCode validate(
    const PETScWrappers::MPI::Vector& phi0,
    const PETScWrappers::MPI::Vector& phi1,
    double eigen,
    bool adj,
    double& norm_sum);

  PetscErrorCode gram_schmidt_bior (
    std::vector<
        PETScWrappers::MPI::BlockVector>& phi,
    std::vector<
        PETScWrappers::MPI::BlockVector>& phiadj);

  // Operators
  TransportMatrixBase<dim, n_fe_degree>& L;
  FisionMatrixBase<dim, n_fe_degree>& M;

  // PETSc objects
  EPS eps;
  Mat shell_mat, shell_mat_adj;
  KSP ksp_00, ksp_11;
  MPI_Comm& comm;

  // Adjoint problem
  bool adjoint;
  std::string adjoint_solver;
  std::vector<PETScWrappers::MPI::BlockVector> phi_adj;

  int n_iterations;
  int n_multiplications;
  double tol_eps;
  double tol_ksp;

  unsigned int max_iterations_eps;
  unsigned int max_iterations_ksp;
  unsigned int n_iterations_solver_00;
  unsigned int n_iterations_solver_11;
  unsigned int n_uses_solver_00;
  unsigned int n_uses_solver_11;

  const unsigned int n_dofs;
  const unsigned int n_dofs_local;
  const unsigned int n_eigenvalues;
  unsigned int n_groups;
};

/**
 * @brief Setup a KSP before using it.
 */
PetscErrorCode ksp_setup(
  KSP &ksp,
  PETScWrappers::MPI::SparseMatrix &Matrix,
  double tol_ksp,
  unsigned int max_iterations_ksp);

/**
 * @brief Change the tolerance of a ksp.
 */
PetscErrorCode ksp_change_tol(
  KSP &ksp,
  double tol_ksp,
  unsigned int max_iterations_ksp);

#endif /* EPSSOLVER_2G_H_ */
