/**
 * @file   pc_multilevel.h
 * @brief Implementation of the neutron noise
 */

#ifndef PC_MULTILEVEL_H_
#define PC_MULTILEVEL_H_

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/types.h>
#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/la_vector.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include "matrix_operators/matrix_operators_base.h"
#include "matrix_operators/matrix_operators_complex_base.h"
#include "matrix_operators/matrix_operators_noise_diff.h"
#include "matrix_operators/matrix_operators_petsc.h"
#include "matrix_operators/matrix_operators_spn.h"
#include "matrix_operators/matrix_operators_petsc_time.h"
#include "matrix_operators/matrix_operators_spn_time.h"

#include "complex_perturbation.h"

#include "utils.h"

using namespace dealii;
/**
 * Preconditioning with a Chebyshev polynomial for symmetric positive definite
 * matrices. This preconditioner is based on an iteration of an inner
 * preconditioner of type @p DiagonalMatrix<PETScWrappers::MPI::Vector> with coefficients that are
 * adapted to optimally cover an eigenvalue range between the largest
 * eigenvalue down to a given lower eigenvalue specified by the optional
 * parameter @p smoothing_range. The typical use case for the preconditioner
 * is a Jacobi preconditioner specified through DiagonalMatrix, which is also
 * the default value for the preconditioner. Note that if the degree variable
 * is set to zero, the Chebyshev iteration corresponds to a Jacobi
 * preconditioner (or the underlying preconditioner type) with relaxation
 * parameter according to the specified smoothing range.
 *
 * Besides the default choice of a pointwise Jacobi preconditioner, this class
 * also allows for more advanced types of preconditioners, for example
 * iterating block-Jacobi preconditioners in DG methods.
 *
 * Apart from the inner preconditioner object, this iteration does not need
 * access to matrix entries, which makes it an ideal ingredient for
 * matrix-free computations. In that context, this class can be used as a
 * multigrid smoother that is trivially %parallel (assuming that matrix-vector
 * products are %parallel and the inner preconditioner is %parallel). Its use
 * is demonstrated in the step-37 tutorial program.
 *
 * <h4>Algorithm execution</h4>
 *
 * The Chebyshev method relies on an estimate of the eigenvalues of the matrix
 * which are computed during the first invocation of vmult(). The algorithm
 * invokes a conjugate gradient solver so symmetry and positive definiteness
 * of the (preconditioned) matrix system are strong requirements. The
 * computation of eigenvalues needs to be deferred until the first vmult()
 * invocation because temporary vectors of the same layout as the source and
 * destination vectors are necessary for these computations and this
 * information gets only available through vmult().
 *
 * The estimation of eigenvalues can also be bypassed by setting
 * SmootherChebyshev::AdditionalData::eig_cg_n_iterations to zero and
 * providing sensible values for the largest eigenvalues in the field
 * SmootherChebyshev::AdditionalData::max_eigenvalue. If the range
 * <tt>[max_eigenvalue/smoothing_range, max_eigenvalue]</tt> contains all
 * eigenvalues of the preconditioned matrix system and the degree (i.e.,
 * number of iterations) is high enough, this class can also be used as a
 * direct solver. For an error estimation of the Chebyshev iteration that can
 * be used to determine the number of iteration, see Varga (2009).
 *
 * In order to use Chebyshev as a solver, set the degree to
 * numbers::invalid_unsigned_int to force the automatic computation of the
 * number of iterations needed to reach a given target tolerance. In this
 * case, the target tolerance is read from the variable
 * SmootherChebyshev::AdditionalData::smoothing_range (it needs to be a
 * number less than one to force any iterations obviously).
 *
 * For details on the algorithm, see section 5.1 of
 * @code{.bib}
 * @Book{Varga2009,
 *   Title       = {Matrix iterative analysis},
 *   Author      = {Varga, R. S.},
 *   Publisher   = {Springer},
 *   Address     = {Berlin},
 *   Edition     = {2nd},
 *   Year        = {2009},
 * }
 * @endcode
 *
 * <h4>Requirements on the templated classes</h4>
 *
 * The class MatrixType must be derived from Subscriptor because a
 * SmartPointer to MatrixType is held in the class. In particular, this means
 * that the matrix object needs to persist during the lifetime of
 * SmootherChebyshev. The preconditioner is held in a shared_ptr that is
 * copied into the AdditionalData member variable of the class, so the
 * variable used for initialization can safely be discarded after calling
 * initialize(). Both the matrix and the preconditioner need to provide @p
 * vmult functions for the matrix-vector product and @p m functions for
 * accessing the number of rows in the (square) matrix. Furthermore, the
 * matrix must provide <tt>el(i,i)</tt> methods for accessing the matrix
 * diagonal in case the preconditioner type is a diagonal matrix. Even though
 * it is highly recommended to pass the inverse diagonal entries inside a
 * separate preconditioner object for implementing the Jacobi method (which is
 * the only possible way to operate this class when computing in %parallel
 * with MPI because there is no knowledge about the locally stored range of
 * entries that would be needed from the matrix alone), there is a backward
 * compatibility function that can extract the diagonal in case of a serial
 * computation.
 *
 * @author Martin Kronbichler, 2009, 2016; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
template <class MatrixType>
  class SmootherChebyshev
  {
    public:
    /**
     * Declare type for container size.
     */
    typedef unsigned int size_type;

    // avoid warning about use of deprecated variables

    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData (const unsigned int degree = 0,
        const double smoothing_range = 0.,
        const bool nonzero_starting =
                                      false,
        const unsigned int eig_cg_n_iterations = 8,
        const double eig_cg_residual = 1e-2,
        const double max_eigenvalue = 1);

      /**
       * This determines the degree of the Chebyshev polynomial. The degree of
       * the polynomial gives the number of matrix-vector products to be
       * performed for one application of the vmult() operation. Degree zero
       * corresponds to a damped Jacobi method.
       *
       * If the degree is set to numbers::invalid_unsigned_int, the algorithm
       * will automatically determine the number of necessary iterations based
       * on the usual Chebyshev error formula as mentioned in the discussion of
       * the main class.
       */
      unsigned int degree;

      /**
       * This sets the range between the largest eigenvalue in the matrix and
       * the smallest eigenvalue to be treated. If the parameter is set to a
       * number less than 1, an estimate for the largest and for the smallest
       * eigenvalue will be calculated internally. For a smoothing range larger
       * than one, the Chebyshev polynomial will act in the interval
       * $[\lambda_\mathrm{max}/ \tt{smoothing\_range}, \lambda_\mathrm{max}]$,
       * where $\lambda_\mathrm{max}$ is an estimate of the maximum eigenvalue
       * of the matrix. A choice of <tt>smoothing_range</tt> between 5 and 20 is
       * useful in case the preconditioner is used as a smoother in multigrid.
       */
      double smoothing_range;

      /**
       * When this flag is set to <tt>true</tt>, it enables the method
       * <tt>vmult(dst, src)</tt> to use non-zero data in the vector
       * <tt>dst</tt>, appending to it the Chebyshev corrections. This can be
       * useful in some situations (e.g. when used for high-frequency error
       * smoothing in a multigrid algorithm), but not the way the solver classes
       * expect a preconditioner to work (where one ignores the content in
       * <tt>dst</tt> for the preconditioner application).
       *
       * @deprecated For non-zero starting, use the step() and Tstep()
       * interfaces, whereas vmult() provides the preconditioner interface.
       */
      bool nonzero_starting;

      /**
       * Maximum number of CG iterations performed for finding the maximum
       * eigenvalue. If set to zero, no computations are performed and the
       * eigenvalues according to the given input are used instead.
       */
      unsigned int eig_cg_n_iterations;

      /**
       * Tolerance for CG iterations performed for finding the maximum
       * eigenvalue.
       */
      double eig_cg_residual;

      /**
       * Maximum eigenvalue to work with. Only in effect if @p
       * eig_cg_n_iterations is set to zero, otherwise this parameter is
       * ignored.
       */
      double max_eigenvalue;

      /**
       * Stores the inverse of the diagonal of the underlying matrix.
       *
       * @deprecated Set the variable @p preconditioner defined below instead.
       */

      PETScWrappers::MPI::Vector matrix_diagonal_inverse;
      /**
       * Stores the preconditioner object that the Chebyshev is wrapped around.
       */
      DiagonalMatrix<PETScWrappers::MPI::Vector> *preconditioner;
    };

    SmootherChebyshev ();

    /**
     * Initialize function. Takes the matrix which is used to form the
     * preconditioner, and additional flags if there are any. This function
     * works only if the input matrix has an operator <tt>el(i,i)</tt> for
     * accessing all the elements in the diagonal. Alternatively, the diagonal
     * can be supplied with the help of the AdditionalData field.
     *
     * This function calculates an estimate of the eigenvalue range of the
     * matrix weighted by its diagonal using a modified CG iteration in case the
     * given number of iterations is positive.
     */
    void initialize (MatrixType *matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Compute the action of the preconditioner on <tt>src</tt>, storing the
     * result in <tt>dst</tt>.
     */
    void vmult (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src);

    /**
     * Compute the action of the transposed preconditioner on <tt>src</tt>,
     * storing the result in <tt>dst</tt>.
     */
    void Tvmult (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src);

    /**
     * Perform one step of the preconditioned Richardson iteration.
     */
    void step (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src);

    /**
     * Perform one transposed step of the preconditioned Richardson iteration.
     */
    void Tstep (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src);

    /**
     * Resets the preconditioner.
     */
    void clear ();

    /**
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is of dimension $m \times n$.
     */
    size_type m () const;

    /**
     * Return the dimension of the domain space. Note that the matrix is of
     * dimension $m \times n$.
     */
    size_type n () const;

    /**
     *
     */
    void set_initial_guess (PETScWrappers::MPI::Vector &vector) const;

    private:

    /**
     * A pointer to the underlying matrix.
     */
    MatrixType *matrix_ptr;

    /**
     * Internal vector used for the <tt>vmult</tt> operation.
     */
    mutable PETScWrappers::MPI::Vector update1;

    /**
     * Internal vector used for the <tt>vmult</tt> operation.
     */
    mutable PETScWrappers::MPI::Vector update2;

    /**
     * Internal vector used for the <tt>vmult</tt> operation.
     */
    mutable PETScWrappers::MPI::Vector update3;

    /**
     * Stores the additional data passed to the initialize function, obtained
     * through a copy operation.
     */
    AdditionalData data;

    /**
     * Average of the largest and smallest eigenvalue under consideration.
     */
    double theta;

    /**
     * Half the interval length between the largest and smallest eigenvalue
     * under consideration.
     */
    double delta;

    /**
     * Stores whether the preconditioner has been set up and eigenvalues have
     * been computed.
     */
    bool eigenvalues_are_initialized;

    /**
     * A mutex to avoid that multiple vmult() invocations by different threads
     * overwrite the temporary vectors.
     */
    mutable Threads::Mutex mutex;

    /**
     * Runs the inner loop of the Chebyshev preconditioner that is the same for
     * vmult() and step() methods.
     */
    void do_chebyshev_loop (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * Runs the inner loop of the Chebyshev preconditioner that is the same for
     * vmult() and step() methods. Uses a separate function to not force users
     * to provide both vmult() and Tvmult() in case only one variant is
     * requested in subsequent calls.
     */
    void do_transpose_chebyshev_loop (PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * Initializes the factors theta and delta based on an eigenvalue
     * computation. If the user set provided values for the largest eigenvalue
     * in AdditionalData, no computation is performed and the information given
     * by the user is used.
     */
    void estimate_eigenvalues (const PETScWrappers::MPI::Vector &src);
  };

template <class MatrixType>
  class FullSmootherChebyshev
  {
    public:
    /**
     * Declare type for container size.
     */
    typedef unsigned int size_type;

    // avoid warning about use of deprecated variables

    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData (const unsigned int degree = 0,
        const double smoothing_range = 0.,
        const bool nonzero_starting =
                                      false,
        const unsigned int eig_cg_n_iterations = 8,
        const double eig_cg_residual = 1e-2,
        const double max_eigenvalue = 1);

      /**
       * This determines the degree of the Chebyshev polynomial. The degree of
       * the polynomial gives the number of matrix-vector products to be
       * performed for one application of the vmult() operation. Degree zero
       * corresponds to a damped Jacobi method.
       *
       * If the degree is set to numbers::invalid_unsigned_int, the algorithm
       * will automatically determine the number of necessary iterations based
       * on the usual Chebyshev error formula as mentioned in the discussion of
       * the main class.
       */
      unsigned int degree;

      /**
       * This sets the range between the largest eigenvalue in the matrix and
       * the smallest eigenvalue to be treated. If the parameter is set to a
       * number less than 1, an estimate for the largest and for the smallest
       * eigenvalue will be calculated internally. For a smoothing range larger
       * than one, the Chebyshev polynomial will act in the interval
       * $[\lambda_\mathrm{max}/ \tt{smoothing\_range}, \lambda_\mathrm{max}]$,
       * where $\lambda_\mathrm{max}$ is an estimate of the maximum eigenvalue
       * of the matrix. A choice of <tt>smoothing_range</tt> between 5 and 20 is
       * useful in case the preconditioner is used as a smoother in multigrid.
       */
      double smoothing_range;

      /**
       * When this flag is set to <tt>true</tt>, it enables the method
       * <tt>vmult(dst, src)</tt> to use non-zero data in the vector
       * <tt>dst</tt>, appending to it the Chebyshev corrections. This can be
       * useful in some situations (e.g. when used for high-frequency error
       * smoothing in a multigrid algorithm), but not the way the solver classes
       * expect a preconditioner to work (where one ignores the content in
       * <tt>dst</tt> for the preconditioner application).
       *
       * @deprecated For non-zero starting, use the step() and Tstep()
       * interfaces, whereas vmult() provides the preconditioner interface.
       */
      bool nonzero_starting;

      /**
       * Maximum number of CG iterations performed for finding the maximum
       * eigenvalue. If set to zero, no computations are performed and the
       * eigenvalues according to the given input are used instead.
       */
      unsigned int eig_cg_n_iterations;

      /**
       * Tolerance for CG iterations performed for finding the maximum
       * eigenvalue.
       */
      double eig_cg_residual;

      /**
       * Maximum eigenvalue to work with. Only in effect if @p
       * eig_cg_n_iterations is set to zero, otherwise this parameter is
       * ignored.
       */
      double max_eigenvalue;

      /**
       * Stores the inverse of the diagonal of the underlying matrix.
       *
       * @deprecated Set the variable @p preconditioner defined below instead.
       */

      PETScWrappers::MPI::BlockVector matrix_diagonal_inverse;
      /**
       * Stores the preconditioner object that the Chebyshev is wrapped around.
       */
      DiagonalMatrix<PETScWrappers::MPI::BlockVector> *preconditioner;
    };

    FullSmootherChebyshev ();

    /**
     * Initialize function. Takes the matrix which is used to form the
     * preconditioner, and additional flags if there are any. This function
     * works only if the input matrix has an operator <tt>el(i,i)</tt> for
     * accessing all the elements in the diagonal. Alternatively, the diagonal
     * can be supplied with the help of the AdditionalData field.
     *
     * This function calculates an estimate of the eigenvalue range of the
     * matrix weighted by its diagonal using a modified CG iteration in case the
     * given number of iterations is positive.
     */
    void initialize (MatrixType *matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Compute the action of the preconditioner on <tt>src</tt>, storing the
     * result in <tt>dst</tt>.
     */
    void vmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src);

    /**
     * Compute the action of the transposed preconditioner on <tt>src</tt>,
     * storing the result in <tt>dst</tt>.
     */
    void Tvmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src);

    /**
     * Perform one step of the preconditioned Richardson iteration.
     */
    void step (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src);

    /**
     * Perform one transposed step of the preconditioned Richardson iteration.
     */
    void Tstep (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src);

    /**
     * Resets the preconditioner.
     */
    void clear ();

    /**
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is of dimension $m \times n$.
     */
    size_type m () const;

    /**
     * Return the dimension of the domain space. Note that the matrix is of
     * dimension $m \times n$.
     */
    size_type n () const;

    /**
     *
     */
    void set_initial_guess (PETScWrappers::MPI::BlockVector &vector) const;

    /**
     * A pointer to the underlying matrix.
     */
    MatrixType *matrix_ptr;

    /**
     * Internal vector used for the <tt>vmult</tt> operation.
     */
    mutable PETScWrappers::MPI::BlockVector update1;

    /**
     * Internal vector used for the <tt>vmult</tt> operation.
     */
    mutable PETScWrappers::MPI::BlockVector update2;

    /**
     * Internal vector used for the <tt>vmult</tt> operation.
     */
    mutable PETScWrappers::MPI::BlockVector update3;

    /**
     * Stores the additional data passed to the initialize function, obtained
     * through a copy operation.
     */
    AdditionalData data;

    /**
     * Average of the largest and smallest eigenvalue under consideration.
     */
    double theta;

    /**
     * Half the interval length between the largest and smallest eigenvalue
     * under consideration.
     */
    double delta;

    /**
     * Stores whether the preconditioner has been set up and eigenvalues have
     * been computed.
     */
    bool eigenvalues_are_initialized;

    /**
     * A mutex to avoid that multiple vmult() invocations by different threads
     * overwrite the temporary vectors.
     */
    mutable Threads::Mutex mutex;

    /**
     * Runs the inner loop of the Chebyshev preconditioner that is the same for
     * vmult() and step() methods.
     */
    void do_chebyshev_loop (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * Runs the inner loop of the Chebyshev preconditioner that is the same for
     * vmult() and step() methods. Uses a separate function to not force users
     * to provide both vmult() and Tvmult() in case only one variant is
     * requested in subsequent calls.
     */
    void do_transpose_chebyshev_loop (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * Initializes the factors theta and delta based on an eigenvalue
     * computation. If the user set provided values for the largest eigenvalue
     * in AdditionalData, no computation is performed and the information given
     * by the user is used.
     */
    void estimate_eigenvalues (const PETScWrappers::MPI::BlockVector &src);

    private:
  };

template <class MatrixType>
  void initialize_preconditioner (const MatrixType &matrix,
    DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
    PETScWrappers::MPI::Vector &diagonal_inverse);

template <class MatrixType>
  void initialize_preconditioner (const MatrixType &matrix,
    DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
    PETScWrappers::MPI::BlockVector &diagonal_inverse);

//---------------------------------------------------------------------------

// for deal.II vectors, perform updates for Chebyshev preconditioner all
// at once to reduce memory transfer. Here, we select between general
// vectors and deal.II vectors where we expand the loop over the (local)
// size of the vector

// generic part for non-deal.II vectors
inline
void vector_updates (const PETScWrappers::MPI::Vector &src,
  const DiagonalMatrix<PETScWrappers::MPI::Vector> &preconditioner,
  const bool start_zero,
  const double factor1,
  const double factor2,
  PETScWrappers::MPI::Vector &update1,
  PETScWrappers::MPI::Vector &update2,
  PETScWrappers::MPI::Vector &update3,
  PETScWrappers::MPI::Vector &dst)
{
  if (start_zero)
  {
    update1.equ(factor2, src);
    preconditioner.vmult(dst, update1);
    update1.equ(-1., dst);
  }
  else
  {
    update2 -= src;
    preconditioner.vmult(update3, update2);
    update2 = update3;
    if (factor1 == 0.)
      update1.equ(factor2, update2);
    else
      update1.sadd(factor1, factor2, update2);
    dst -= update1;
  }
}

inline
void vector_updates (const PETScWrappers::MPI::BlockVector &src,
  const DiagonalMatrix<PETScWrappers::MPI::BlockVector> &preconditioner,
  const bool start_zero,
  const double factor1,
  const double factor2,
  PETScWrappers::MPI::BlockVector &update1,
  PETScWrappers::MPI::BlockVector &update2,
  PETScWrappers::MPI::BlockVector &update3,
  PETScWrappers::MPI::BlockVector &dst)
{
  if (start_zero)
  {
    update1.equ(factor2, src);
//    preconditioner.vmult(dst, update1);
    update1.equ(-1., dst);
  }
  else
  {
    update2 -= src;
    preconditioner.vmult(update3, update2);
    update2 = update3;
    if (factor1 == 0.)
      update1.equ(factor2, update2);
    else
      update1.sadd(factor1, factor2, update2);
    dst -= update1;
  }
}

//// worker routine for deal.II vectors. Because of vectorization, we need
//// to put the loop into an extra structure because the virtual function of
//// VectorUpdatesRange prevents the compiler from applying vectorization.
//template <typename Number>
//  struct VectorUpdater
//  {
//    VectorUpdater (const Number *src,
//      const Number *matrix_diagonal_inverse,
//      const bool start_zero,
//      const Number factor1,
//      const Number factor2,
//      Number *update1,
//      Number *update2,
//      Number *dst) :
//        src(src),
//        matrix_diagonal_inverse(matrix_diagonal_inverse),
//        do_startup(
//          factor1 == Number()),
//        start_zero(start_zero),
//        factor1(
//          factor1),
//        factor2(factor2),
//        update1(update1),
//        update2(
//          update2),
//        dst(dst)
//    {
//    }
//
//    void apply_to_subrange (const std::size_t begin,
//      const std::size_t end) const
//      {
//      // To circumvent a bug in gcc
//      // (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63945), we create copies
//      // of the variables factor1 and factor2 and do not check based on
//      // factor1.
//      const Number factor1 = this->factor1;
//      const Number factor2 = this->factor2;
//      if (do_startup)
//      {
//        if (start_zero)
//          DEAL_II_OPENMP_SIMD_PRAGMA
//          for (std::size_t i = begin; i < end; ++i)
//          {
//            dst[i] = factor2 * src[i] * matrix_diagonal_inverse[i];
//            update1[i] = -dst[i];
//          }
//        else
//          DEAL_II_OPENMP_SIMD_PRAGMA
//          for (std::size_t i = begin; i < end; ++i)
//          {
//            update1[i] = ((update2[i] - src[i]) * factor2
//                          * matrix_diagonal_inverse[i]);
//            dst[i] -= update1[i];
//          }
//      }
//      else
//        DEAL_II_OPENMP_SIMD_PRAGMA
//        for (std::size_t i = begin; i < end; ++i)
//        {
//          const Number update = factor1 * update1[i]
//                                + factor2
//                                  * ((update2[i] - src[i])
//                                     * matrix_diagonal_inverse[i]);
//          update1[i] = update;
//          dst[i] -= update;
//        }
//    }
//
//    const Number *src;
//    const Number *matrix_diagonal_inverse;
//    const bool do_startup;
//    const bool start_zero;
//    const Number factor1;
//    const Number factor2;
//    mutable Number *update1;
//    mutable Number *update2;
//    mutable Number *dst;
//  };
//
//template <typename Number>
//  struct VectorUpdatesRange : public parallel::ParallelForInteger
//  {
//    VectorUpdatesRange (const VectorUpdater<Number> &updater,
//      const std::size_t size) :
//        updater(updater)
//    {
//      if (size < internal::VectorImplementation::minimum_parallel_grain_size)
//        apply_to_subrange(0, size);
//      else
//        apply_parallel(0, size,
//          internal::VectorImplementation::minimum_parallel_grain_size);
//    }
//
//    ~VectorUpdatesRange () = default;
//
//    virtual void apply_to_subrange (const std::size_t begin,
//      const std::size_t end) const
//      {
//      updater.apply_to_subrange(begin, end);
//    }
//
//    const VectorUpdater<Number> &updater;
//  };

/**
 *
 */
struct EigenvalueTracker
{
  public:
  void slot (const std::vector<double> &eigenvalues)
  {
    values = eigenvalues;
  }

  std::vector<double> values;
};

/**
 *
 */
struct EigenvalueTrackerComplex
{
  public:
  void slot (const std::vector<std::complex<double>> &eigenvalues)
  {
    values_complex = eigenvalues;
  }

  std::vector<std::complex<double>> values_complex;
};

// ------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------- //
/*
 * **********************
 *     PC_GMG Class     *
 * **********************
 */
template <int dim, int n_fe_degree>
  void shell_D_coarse_gmg (Mat shell_mat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  class PC_MLFE
  {
    public:

    PC_MLFE (const Triangulation<dim> &tria,
      unsigned int p_coarse,
      const DoFHandler<dim> &dof_handler,
      const FE_Q<dim> &fe,
      TransportMatrixBase<dim, n_fe_degree> &_L,
      FisionMatrixBase<dim, n_fe_degree> &_M);

    void reinit ();

    void clear ();

    void apply_gmg (unsigned int ngr,
      unsigned int ngc,
      PETScWrappers::MPI::Vector &in,
      PETScWrappers::MPI::Vector &out);

    void apply_gmg (double rho,
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &out);

    void ksp_setup (KSP &ksp);

    void ksp_block_setup (KSP &ksp,
      unsigned int nb);

    const Triangulation<dim> &tria;
    unsigned int n_fe_coarse;

    const DoFHandler<dim> &dof_handler;
    const FE_Q<dim> &fe;

    TransportMatrixBase<dim, n_fe_degree> &L;
    FisionMatrixBase<dim, n_fe_degree> &M;

    FE_Q<dim> fe_coarse;
    DoFHandler<dim> dof_handler_coarse;

    FullMatrix<double> transfer_mat_restrict;
    FullMatrix<double> transfer_mat_prolongate;

    unsigned int ndofs_coarse;
    unsigned int n_eigenvalues;
    unsigned int n_blocks;

    KSP ksp;
    double lambda;
    Mat shellD;

    std::vector<KSP> ksp_blocks;
  };

template <int dim, int n_fe_degree>
  void shell_T_coarse (Mat shell_mat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_T_coarse (PC shell_mat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_T_coarse_gs_cgilu (PC shell_mat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  class PC_MLFE_Time
  {
    public:

    PC_MLFE_Time (const MPI_Comm &_comm,
      const TransportMatrixBase<dim, n_fe_degree> &_T,
      const Materials &_materials);

    void reinit ();

    void clear ();

    void apply_gmg (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void restrict_solution (
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &in_coarse);

    void prolongate_solution (
      PETScWrappers::MPI::BlockVector &out_coarse,
      PETScWrappers::MPI::BlockVector &out);

    void ksp_setup ();

    void pc_setup ();

    void pc_gs_cgilu_setup ();

    const MPI_Comm comm;

    const TransportMatrixBase<dim, n_fe_degree> &T;
    const Materials &materials;

    const DoFHandler<dim> &dof_handler;
    DoFHandler<dim> dof_handler_coarse;

    // Parameters for T_coarse
    std::string transport_fine;
    unsigned int n_fe_coarse;
    unsigned int n_blocks_coarse;
    std::string transport_coarse;
    SystemMatrixTime<dim, 1> T_coarse_1;
    SystemMatrixTime<dim, 2> T_coarse_2;
    SystemMatrixTime<dim, 3> T_coarse_3;

    SystemMatrixTimeSPN<dim, 1> T_coarse_spn_1;
    SystemMatrixTimeSPN<dim, 2> T_coarse_spn_2;

    AffineConstraints<double> constraints_coarse;

    std::vector<unsigned int> boundary_conditions;
    std::vector<double> albedo_factors;
    double delta_t;
    std::string time_scheme;
    MatrixFreeType matrixfree_type_coarse;
    std::string pc_coarse;

    FullMatrix<double> transfer_mat_restrict;
    FullMatrix<double> transfer_mat_prolongate;

    unsigned int n_dofs_coarse;
    unsigned int n_dofs_coarse_local;

    unsigned int n_blocks;

    KSP ksp;
    double lambda;
    Mat shellT_coarse;

    bool first_computation;

    DiagonalMatrix<PETScWrappers::MPI::BlockVector> prec_diag;
    std::vector<PC> pc_blocks;
    std::vector<KSP> ksp_blocks;

    unsigned int total_its;
    unsigned int n_applications;

  };

template <int dim, int n_fe_degree>
  class PC_MLFE_Complex
  {
    public:

    PC_MLFE_Complex (const MPI_Comm &_comm,
      const TransportMatrixComplexBase<dim, n_fe_degree> &_T,
      const Materials &_materials,
      const ComplexPerturbation &_pert);

    void reinit ();

    void clear ();

    void apply_gmg (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void restrict_solution (
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &in_coarse);

    void prolongate_solution (
      PETScWrappers::MPI::BlockVector &out_coarse,
      PETScWrappers::MPI::BlockVector &out);

    void ksp_setup ();

    void pc_setup ();

    void pc_gs_cgilu_setup ();

    const MPI_Comm comm;

    const TransportMatrixComplexBase<dim, n_fe_degree> &A;
    const Materials &materials;
    const ComplexPerturbation &pert;

    const DoFHandler<dim> &dof_handler;
    DoFHandler<dim> dof_handler_coarse;

    // Parameters for T_coarse
    //std::string transport_fine;
    unsigned int n_fe_coarse;
    unsigned int n_blocks_coarse;
    std::string transport_coarse;

    NoiseAMatrix<dim, 1> A_coarse_1;
    NoiseAMatrix<dim, 2> A_coarse_2;
    NoiseAMatrix<dim, 3> A_coarse_3;

//    SystemMatrixTimeSPN<dim, 1> A_coarse_spn_1;
//    SystemMatrixTimeSPN<dim, 2> A_coarse_spn_2;

    AffineConstraints<double> constraints_coarse;

    std::vector<unsigned int> boundary_conditions;
    std::vector<double> albedo_factors;

    MatrixFreeType matrixfree_type_coarse;
    std::string pc_coarse;

    FullMatrix<double> transfer_mat_restrict;
    FullMatrix<double> transfer_mat_prolongate;

    unsigned int n_dofs_coarse;
    unsigned int n_dofs_coarse_local;

    unsigned int n_blocks;

    KSP ksp;
    double lambda;
    Mat shellA_coarse;

    bool first_computation;

    DiagonalMatrix<PETScWrappers::MPI::BlockVector> prec_diag;
    std::vector<PC> pc_blocks;
    std::vector<KSP> ksp_blocks;

    unsigned int total_its;
    unsigned int n_applications;

  };

template <int dim, int n_fe_degree>
  void shell_A_coarse (Mat shell_mat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_A_coarse (PC shell_mat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_A_coarse_gs_cgilu (PC shell_mat,
    Vec src,
    Vec dst);

#endif /* PRECONDITION_CHEBYSHEV_H_ */
