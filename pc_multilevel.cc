/**
 * @file   pc_multilevel.cc
 * @brief Implementation of the neutron noise
 */

#include "pc_multilevel.h"

using namespace dealii;

//template struct VectorUpdater<double> ;

/**
 *
 */
template <class MatrixType>
  SmootherChebyshev<MatrixType>::AdditionalData::AdditionalData (
    const unsigned int degree,
    const double smoothing_range,
    const bool nonzero_starting,
    const unsigned int eig_cg_n_iterations,
    const double eig_cg_residual,
    const double max_eigenvalue) :
      degree(degree),
      smoothing_range(smoothing_range),
      nonzero_starting(
        nonzero_starting),
      eig_cg_n_iterations(eig_cg_n_iterations),
      eig_cg_residual(
        eig_cg_residual),
      max_eigenvalue(max_eigenvalue),
      preconditioner(
        nullptr)
  {
  }

/**
 *
 */
template <class MatrixType>
  SmootherChebyshev<MatrixType>::SmootherChebyshev () :
      matrix_ptr(nullptr),
      theta(1.),
      delta(1.),
      eigenvalues_are_initialized(
        false)
  {
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::initialize (MatrixType *matrix,
    const AdditionalData &additional_data)
  {
    matrix_ptr = matrix;
    data = additional_data;
    initialize_preconditioner<MatrixType>(*matrix, data.preconditioner,
      data.matrix_diagonal_inverse);

    eigenvalues_are_initialized = false;
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::clear ()
  {
    eigenvalues_are_initialized = false;
    theta = delta = 1.0;
    matrix_ptr = nullptr;
    {
      PETScWrappers::MPI::Vector empty_vector;
      data.matrix_diagonal_inverse.reinit(empty_vector);
      update1.reinit(empty_vector);
      update2.reinit(empty_vector);
      update3.reinit(empty_vector);
    }
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::estimate_eigenvalues (
    const PETScWrappers::MPI::Vector &src)
  {
    Assert(eigenvalues_are_initialized == false, ExcInternalError());
    // Assert(data.preconditioner.get() != nullptr, ExcNotInitialized());

    update1.reinit(src);
    update2.reinit(src, true);

    // calculate largest eigenvalue using a hand-tuned CG iteration on the
    // matrix weighted by its diagonal. we start with a vector that consists of
    // ones only, weighted by the length.
    double max_eigenvalue, min_eigenvalue;
    if (data.eig_cg_n_iterations > 0)
    {
      Assert(data.eig_cg_n_iterations > 2,
        ExcMessage ("Need to set at least two iterations to find eigenvalues."));

      // set a very strict tolerance to force at least two iterations
      ReductionControl control(data.eig_cg_n_iterations,
        std::sqrt(std::numeric_limits<double>::epsilon()), 1e-10, false,
        false);

      EigenvalueTracker eigenvalue_tracker;
      SolverCG<PETScWrappers::MPI::Vector> solver(control);
      solver.connect_eigenvalues_slot(
        std::bind(&EigenvalueTracker::slot, &eigenvalue_tracker,
          std::placeholders::_1));

      // set an initial guess which is close to the constant vector but where
      // one entry is different to trigger high frequencies
      set_initial_guess(update2);

      update2.compress(VectorOperation::insert);
      update1.compress(VectorOperation::insert);
      try
      {
        solver.solve(*matrix_ptr, update1, update2, *data.preconditioner);
      }
      catch (SolverControl::NoConvergence&)
      {
      }

      // read the eigenvalues from the attached eigenvalue tracker
      if (eigenvalue_tracker.values.empty())
        min_eigenvalue = max_eigenvalue = 1;
      else
      {
        min_eigenvalue = eigenvalue_tracker.values.front();

        // include a safety factor since the CG method will in general not
        // be converged
        max_eigenvalue = 1.2 * eigenvalue_tracker.values.back();
      }
    }
    else
    {
      max_eigenvalue = data.max_eigenvalue;
      min_eigenvalue = data.max_eigenvalue / data.smoothing_range;
    }

    const double alpha = (
                         data.smoothing_range > 1. ?
                             max_eigenvalue / data.smoothing_range :
                             std::min(0.9 * max_eigenvalue, min_eigenvalue));

    // in case the user set the degree to invalid unsigned int, we have to
    // determine the number of necessary iterations from the Chebyshev error
    // estimate, given the target tolerance specified by smoothing_range. This
    // estimate is based on the error formula given in section 5.1 of
    // R. S. Varga, Matrix iterative analysis, 2nd ed., Springer, 2009
    if (data.degree == numbers::invalid_unsigned_int)
    {
      const double actual_range = max_eigenvalue / alpha;
      const double sigma = (1. - std::sqrt(1. / actual_range))
                           / (1. + std::sqrt(1. / actual_range));
      const double eps = data.smoothing_range;
      const_cast<SmootherChebyshev<MatrixType>*>(this)->data.degree = 1
          + std::log(1. / eps + std::sqrt(1. / eps / eps - 1))
            / std::log(1. / sigma);
    }

    const_cast<SmootherChebyshev<MatrixType>*>(this)->delta = (max_eigenvalue
                                                               - alpha)
                                                              * 0.5;
    const_cast<SmootherChebyshev<MatrixType>*>(this)->theta = (max_eigenvalue
                                                               + alpha)
                                                              * 0.5;

    // We do not need the third auxiliary vector in case we have a
    // DiagonalMatrix as preconditioner and use deal.II's own vectors
    update3.reinit(src, true);

    const_cast<SmootherChebyshev<MatrixType>*>(this)->eigenvalues_are_initialized =
                                                                                    true;
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::do_chebyshev_loop (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    // if delta is zero, we do not need to iterate because the updates will be
    // zero
    if (std::abs(delta) < 1e-40)
      return;

    double rhok = delta / theta, sigma = theta / delta;
    for (unsigned int k = 0; k < data.degree; ++k)
    {
      matrix_ptr->vmult(update2, dst);
      const double rhokp = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      vector_updates(src, *data.preconditioner, false, factor1, factor2,
        update1, update2, update3, dst);
    }
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::do_transpose_chebyshev_loop (
    PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src) const
  {
    double rhok = delta / theta, sigma = theta / delta;
    for (unsigned int k = 0; k < data.degree; ++k)
    {
      matrix_ptr->Tvmult(update2, dst);
      const double rhokp = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      vector_updates(src, *data.preconditioner, false, factor1, factor2,
        update1, update2, update3, dst);
    }
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::set_initial_guess (
    PETScWrappers::MPI::Vector &vector) const
  {
    vector = 1. / std::sqrt(static_cast<double>(vector.size()));
    if (vector.locally_owned_elements().is_element(0))
      vector(0) = 0.;
  }

template <class MatrixType>
  void SmootherChebyshev<MatrixType>::vmult (PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    vector_updates(src, *data.preconditioner, true, 0., 1. / theta, update1,
      update2, update3, dst);

    do_chebyshev_loop(dst, src);
  }

template <class MatrixType>
  void SmootherChebyshev<MatrixType>::Tvmult (PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    vector_updates(src, *data.preconditioner, true, 0., 1. / theta, update1,
      update2, update3, dst);

    do_transpose_chebyshev_loop(dst, src);
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::step (PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    matrix_ptr->vmult(update2, dst);
    vector_updates(src, *data.preconditioner, false, 0., 1. / theta, update1,
      update2, update3, dst);

    do_chebyshev_loop(dst, src);
  }

/**
 *
 */
template <class MatrixType>
  void SmootherChebyshev<MatrixType>::Tstep (PETScWrappers::MPI::Vector &dst,
    const PETScWrappers::MPI::Vector &src)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    matrix_ptr->Tvmult(update2, dst);
    vector_updates(src, *data.preconditioner, false, 0., 1. / theta, update1,
      update2, update3, dst);

    do_transpose_chebyshev_loop(dst, src);
  }

/**
 *
 */
template <class MatrixType>
  typename SmootherChebyshev<MatrixType>::size_type SmootherChebyshev<MatrixType>::m () const
  {
    Assert(matrix_ptr!=nullptr, ExcNotInitialized());
    return matrix_ptr->m();
  }

/**
 *
 */
template <class MatrixType>
  typename SmootherChebyshev<MatrixType>::size_type SmootherChebyshev<MatrixType>::n () const
  {
    Assert(matrix_ptr!=nullptr, ExcNotInitialized());
    return matrix_ptr->n();
  }

template class SmootherChebyshev<PETScWrappers::MPI::SparseMatrix> ;
template class SmootherChebyshev<PoissonOperator<1, 1, double> > ;
template class SmootherChebyshev<PoissonOperator<1, 2, double> > ;
template class SmootherChebyshev<PoissonOperator<1, 3, double> > ;
template class SmootherChebyshev<PoissonOperator<1, 4, double> > ;
template class SmootherChebyshev<PoissonOperator<1, 5, double> > ;

template class SmootherChebyshev<PoissonOperator<2, 1, double> > ;
template class SmootherChebyshev<PoissonOperator<2, 2, double> > ;
template class SmootherChebyshev<PoissonOperator<2, 3, double> > ;
template class SmootherChebyshev<PoissonOperator<2, 4, double> > ;
template class SmootherChebyshev<PoissonOperator<2, 5, double> > ;

template class SmootherChebyshev<PoissonOperator<3, 1, double> > ;
template class SmootherChebyshev<PoissonOperator<3, 2, double> > ;
template class SmootherChebyshev<PoissonOperator<3, 3, double> > ;
template class SmootherChebyshev<PoissonOperator<3, 4, double> > ;
template class SmootherChebyshev<PoissonOperator<3, 5, double> > ;

/**
 *
 */
template <class MatrixType>
  FullSmootherChebyshev<MatrixType>::AdditionalData::AdditionalData (
    const unsigned int degree,
    const double smoothing_range,
    const bool nonzero_starting,
    const unsigned int eig_cg_n_iterations,
    const double eig_cg_residual,
    const double max_eigenvalue) :
      degree(degree),
      smoothing_range(smoothing_range),
      nonzero_starting(
        nonzero_starting),
      eig_cg_n_iterations(eig_cg_n_iterations),
      eig_cg_residual(
        eig_cg_residual),
      max_eigenvalue(max_eigenvalue),
      preconditioner(
        nullptr)
  {
  }

/**
 *
 */
template <class MatrixType>
  FullSmootherChebyshev<MatrixType>::FullSmootherChebyshev () :
      matrix_ptr(nullptr),
      theta(1.),
      delta(1.),
      eigenvalues_are_initialized(
        false)
  {
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::initialize (MatrixType *matrix,
    const AdditionalData &additional_data)
  {
    matrix_ptr = matrix;
    data = additional_data;
    initialize_preconditioner<MatrixType>(*matrix, data.preconditioner,
      data.matrix_diagonal_inverse);

    eigenvalues_are_initialized = false;
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::clear ()
  {
    eigenvalues_are_initialized = false;
    theta = delta = 1.0;
    matrix_ptr = nullptr;
    {
      PETScWrappers::MPI::BlockVector empty_vector;
      data.matrix_diagonal_inverse.reinit(empty_vector);
      update1.reinit(empty_vector);
      update2.reinit(empty_vector);
      update3.reinit(empty_vector);
    }
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::estimate_eigenvalues (
    const PETScWrappers::MPI::BlockVector &src)
  {
    Assert(eigenvalues_are_initialized == false, ExcInternalError());
    // Assert(data.preconditioner.get() != nullptr, ExcNotInitialized());

    update1.reinit(src);

    // calculate largest eigenvalue using a hand-tuned CG iteration on the
    // matrix weighted by its diagonal. we start with a vector that consists of
    // ones only, weighted by the length.
    double max_eigenvalue, min_eigenvalue;
    if (data.eig_cg_n_iterations > 0)
    {
      Assert(data.eig_cg_n_iterations > 2,
        ExcMessage ("Need to set at least two iterations to find eigenvalues."));

      // set a very strict tolerance to force at least two iterations
      ReductionControl control(data.eig_cg_n_iterations,
        std::sqrt(std::numeric_limits<double>::epsilon()), 1e-10, false,
        false);

      EigenvalueTrackerComplex eigenvalue_tracker;
      SolverGMRES<PETScWrappers::MPI::BlockVector> solver(control);
      solver.connect_eigenvalues_slot(
        std::bind(&EigenvalueTrackerComplex::slot, &eigenvalue_tracker,
          std::placeholders::_1));

      // set an initial guess which is close to the constant vector but where
      // one entry is different to trigger high frequencies

      if (data.nonzero_starting == false)
        set_initial_guess(update2);

      update2.compress(VectorOperation::insert);
      update1.compress(VectorOperation::insert);
      try
      {
        solver.solve(*matrix_ptr, update1, update2, *data.preconditioner);
      }
      catch (SolverControl::NoConvergence&)
      {
      }

      // read the eigenvalues from the attached eigenvalue tracker
      if (eigenvalue_tracker.values_complex.empty())
        min_eigenvalue = max_eigenvalue = 1;
      else
      {
        min_eigenvalue = std::real(
          eigenvalue_tracker.values_complex.front());

        // include a safety factor since the CG method will in general not
        // be converged
        max_eigenvalue = 1.2
                         * std::real(eigenvalue_tracker.values_complex.back());
      }
    }
    else
    {
      max_eigenvalue = data.max_eigenvalue;
      min_eigenvalue = data.max_eigenvalue / data.smoothing_range;
    }

    const double alpha = (
                         data.smoothing_range > 1. ?
                             max_eigenvalue / data.smoothing_range :
                             std::min(0.9 * max_eigenvalue, min_eigenvalue));

    // in case the user set the degree to invalid unsigned int, we have to
    // determine the number of necessary iterations from the Chebyshev error
    // estimate, given the target tolerance specified by smoothing_range. This
    // estimate is based on the error formula given in section 5.1 of
    // R. S. Varga, Matrix iterative analysis, 2nd ed., Springer, 2009
    if (data.degree == numbers::invalid_unsigned_int)
    {
      const double actual_range = max_eigenvalue / alpha;
      const double sigma = (1. - std::sqrt(1. / actual_range))
                           / (1. + std::sqrt(1. / actual_range));
      const double eps = data.smoothing_range;
      const_cast<FullSmootherChebyshev<MatrixType>*>(this)->data.degree = 1
          + std::log(1. / eps + std::sqrt(1. / eps / eps - 1))
            / std::log(1. / sigma);

    }

    const_cast<FullSmootherChebyshev<MatrixType>*>(this)->delta =
        (max_eigenvalue - alpha) * 0.5;
    const_cast<FullSmootherChebyshev<MatrixType>*>(this)->theta =
        (max_eigenvalue + alpha) * 0.5;

    // We do not need the third auxiliary vector in case we have a
    // DiagonalMatrix as preconditioner and use deal.II's own vectors
    update3.reinit(src, true);

    const_cast<FullSmootherChebyshev<MatrixType>*>(this)->eigenvalues_are_initialized =
        true;

  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::do_chebyshev_loop (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    // if delta is zero, we do not need to iterate because the updates will be
    // zero
    if (std::abs(delta) < 1e-40)
      return;

    double rhok = delta / theta, sigma = theta / delta;
    for (unsigned int k = 0; k < data.degree; ++k)
    {
      matrix_ptr->vmult(update2, dst);
      const double rhokp = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      vector_updates(src, *data.preconditioner, false, factor1, factor2,
        update1, update2, update3, dst);
    }
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::do_transpose_chebyshev_loop (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    double rhok = delta / theta, sigma = theta / delta;
    for (unsigned int k = 0; k < data.degree; ++k)
    {
      matrix_ptr->Tvmult(update2, dst);
      const double rhokp = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      vector_updates(src, *data.preconditioner, false, factor1, factor2,
        update1, update2, update3, dst);
    }
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::set_initial_guess (
    PETScWrappers::MPI::BlockVector &vector) const
  {

    vector = 1. / std::sqrt(static_cast<double>(vector.size()));
    if (vector.locally_owned_elements().is_element(0))
      vector(0) = 0.;
  }

template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::vmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src)
  {

    std::lock_guard<std::mutex> lock(mutex);
    update2.reinit(dst, true);

    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    vector_updates(src, *data.preconditioner, true, 0., 1. / theta, update1,
      update2, update3, dst);

    do_chebyshev_loop(dst, src);
  }

template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::Tvmult (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src)
  {

    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    vector_updates(src, *data.preconditioner, true, 0., 1. / theta, update1,
      update2, update3, dst);

    do_transpose_chebyshev_loop(dst, src);
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::step (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    matrix_ptr->vmult(update2, dst);
    vector_updates(src, *data.preconditioner, false, 0., 1. / theta, update1,
      update2, update3, dst);

    do_chebyshev_loop(dst, src);
  }

/**
 *
 */
template <class MatrixType>
  void FullSmootherChebyshev<MatrixType>::Tstep (
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (eigenvalues_are_initialized == false)
      estimate_eigenvalues(src);

    matrix_ptr->Tvmult(update2, dst);
    vector_updates(src, *data.preconditioner, false, 0., 1. / theta, update1,
      update2, update3, dst);

    do_transpose_chebyshev_loop(dst, src);
  }

/**
 *
 */
template <class MatrixType>
  typename FullSmootherChebyshev<MatrixType>::size_type FullSmootherChebyshev<
      MatrixType>::m () const
  {
    Assert(matrix_ptr!=nullptr, ExcNotInitialized());
    return matrix_ptr->m();
  }

/**
 *
 */
template <class MatrixType>
  typename FullSmootherChebyshev<MatrixType>::size_type FullSmootherChebyshev<
      MatrixType>::n () const
  {
    Assert(matrix_ptr!=nullptr, ExcNotInitialized());
    return matrix_ptr->n();
  }

template class FullSmootherChebyshev<TransportMatrixBase<1, 1> > ;
template class FullSmootherChebyshev<TransportMatrixBase<1, 2> > ;
template class FullSmootherChebyshev<TransportMatrixBase<1, 3> > ;
template class FullSmootherChebyshev<TransportMatrixBase<1, 4> > ;
template class FullSmootherChebyshev<TransportMatrixBase<1, 5> > ;

template class FullSmootherChebyshev<TransportMatrixBase<2, 1> > ;
template class FullSmootherChebyshev<TransportMatrixBase<2, 2> > ;
template class FullSmootherChebyshev<TransportMatrixBase<2, 3> > ;
template class FullSmootherChebyshev<TransportMatrixBase<2, 4> > ;
template class FullSmootherChebyshev<TransportMatrixBase<2, 5> > ;

template class FullSmootherChebyshev<TransportMatrixBase<3, 1> > ;
template class FullSmootherChebyshev<TransportMatrixBase<3, 2> > ;
template class FullSmootherChebyshev<TransportMatrixBase<3, 3> > ;
template class FullSmootherChebyshev<TransportMatrixBase<3, 4> > ;
template class FullSmootherChebyshev<TransportMatrixBase<3, 5> > ;

template class FullSmootherChebyshev<TransportMatrixComplexBase<1, 1> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<1, 2> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<1, 3> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<1, 4> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<1, 5> > ;

template class FullSmootherChebyshev<TransportMatrixComplexBase<2, 1> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<2, 2> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<2, 3> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<2, 4> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<2, 5> > ;

template class FullSmootherChebyshev<TransportMatrixComplexBase<3, 1> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<3, 2> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<3, 3> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<3, 4> > ;
template class FullSmootherChebyshev<TransportMatrixComplexBase<3, 5> > ;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 *
 */
template <class MatrixType>
  void initialize_preconditioner (const MatrixType &matrix,
    DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
    PETScWrappers::MPI::Vector &diagonal_inverse)
  {

    if (preconditioner == nullptr)
    {
      // Check if we can initialize from vector that then gets set to zero
      // as the matrix will own the memory
      preconditioner = new DiagonalMatrix<PETScWrappers::MPI::Vector>();
      preconditioner->reinit(diagonal_inverse);
      {
        PETScWrappers::MPI::Vector empty_vector;
        diagonal_inverse.reinit(empty_vector);
      }

      // This part only works in serial
      if (preconditioner->m() != matrix.m())
      {
        preconditioner->get_vector().reinit(MPI_COMM_WORLD, matrix.m(),
          matrix.m());

        for (unsigned int i = 0; i < matrix.m(); ++i)
          preconditioner->get_vector()(i) = 1. / matrix.el(i, i);

        preconditioner->get_vector().compress(VectorOperation::insert);
      }

    }
    Assert(preconditioner->m() == matrix.m(),
      ExcMessage("Preconditioner appears to be initialized but not sized correctly"));

  }

template void initialize_preconditioner (
  const PETScWrappers::MPI::SparseMatrix &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<1, 1, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<1, 2, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<1, 3, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<1, 4, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<1, 5, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<2, 1, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<2, 2, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<2, 3, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<2, 4, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<2, 5, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<3, 1, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<3, 2, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<3, 3, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<3, 4, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);
template void initialize_preconditioner (
  const PoissonOperator<3, 5, double> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::Vector> *&preconditioner,
  PETScWrappers::MPI::Vector &diagonal_inverse);

/**
 *
 */
template <class MatrixType>
  void initialize_preconditioner (const MatrixType&,
    DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
    PETScWrappers::MPI::BlockVector &diagonal_inverse)
  {
    if (preconditioner == nullptr)
    {
      // Check if we can initialize from vector that then gets set to zero
      // as the matrix will own the memory
      preconditioner = new DiagonalMatrix<PETScWrappers::MPI::BlockVector>();
      preconditioner->reinit(diagonal_inverse);
      {
        PETScWrappers::MPI::BlockVector empty_vector;
        diagonal_inverse.reinit(empty_vector);
      }

    }

  }

template void initialize_preconditioner (const SystemMatrixTime<1, 1> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<1, 2> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<1, 3> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<1, 4> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<1, 5> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<2, 1> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<2, 2> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<2, 3> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<2, 4> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<2, 5> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<3, 1> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<3, 2> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<3, 3> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<3, 4> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
template void initialize_preconditioner (const SystemMatrixTime<3, 5> &matrix,
  DiagonalMatrix<PETScWrappers::MPI::BlockVector> *&preconditioner,
  PETScWrappers::MPI::BlockVector &diagonal_inverse);
// ------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------- //

/**
 *
 */
template <int dim, int n_fe_degree>
  PC_MLFE<dim, n_fe_degree>::PC_MLFE (const Triangulation<dim> &_tria,
    unsigned int p_coarse,
    const DoFHandler<dim> &_dof_handler,
    const FE_Q<dim> &_fe,
    TransportMatrixBase<dim, n_fe_degree> &_L,
    FisionMatrixBase<dim, n_fe_degree> &_M) :
      tria(_tria),
      n_fe_coarse(p_coarse),
      dof_handler(_dof_handler),
      fe(_fe),
      L(
        _L),
      M(_M),
      fe_coarse(QGaussLobatto<1>(p_coarse + 1)),
      dof_handler_coarse(
        tria)
  {
    n_blocks = 0;
    ksp = NULL;
    ndofs_coarse = 0;
    lambda = 1.00;
    n_eigenvalues = 0;
    shellD = NULL;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE<dim, n_fe_degree>::reinit ()
  {

    n_blocks = L.n_blocks_rows();

    ksp_blocks.resize(n_blocks);

    ndofs_coarse = L.block(0, 0).m();

    dof_handler_coarse.initialize(tria, fe_coarse);
    dof_handler_coarse.distribute_dofs(fe_coarse);

    transfer_mat_restrict.reinit(fe_coarse.dofs_per_cell, fe.dofs_per_cell);
    transfer_mat_prolongate.reinit(fe.dofs_per_cell, fe_coarse.dofs_per_cell);

    set_transfer_matrix(fe, fe_coarse, transfer_mat_restrict);
    set_transfer_matrix(fe_coarse, fe, transfer_mat_prolongate);

    MatCreateShell(PETSC_COMM_WORLD, n_blocks * ndofs_coarse,
      n_blocks * ndofs_coarse, n_blocks * ndofs_coarse,
      n_blocks * ndofs_coarse, this, &shellD);
    MatShellSetOperation(shellD, MATOP_MULT,
      (void (*) ()) shell_D_coarse_gmg<dim, n_fe_degree>);
    ;
    ksp_setup(ksp);

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      ksp_block_setup(ksp_blocks[nb], nb);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE<dim, n_fe_degree>::clear ()
  {

    KSPDestroy(&ksp);

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      KSPDestroy(&ksp_blocks[nb]);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE<dim, n_fe_degree>::apply_gmg (double rho,
    PETScWrappers::MPI::BlockVector &in,
    PETScWrappers::MPI::BlockVector &out)
  {
    lambda = rho;

    PETScWrappers::MPI::BlockVector in_coarse;
    PETScWrappers::MPI::BlockVector out_coarse;
    in_coarse.reinit(n_blocks, MPI_COMM_WORLD, ndofs_coarse, ndofs_coarse);
    out_coarse.reinit(n_blocks, MPI_COMM_WORLD, ndofs_coarse, ndofs_coarse);

    Vec incoarse_vec, outcoarse_vec;
    VecCreateSeq(MPI_COMM_WORLD, ndofs_coarse * n_blocks, &incoarse_vec);
    VecDuplicate(incoarse_vec, &outcoarse_vec);

    // Restrict the solution
    solution_interpolate_dof(dof_handler, dof_handler_coarse,
      transfer_mat_restrict, in, in_coarse);

    // Solve the linear system in the linear elements
    copy_to_Vec(incoarse_vec, in_coarse);
    KSPSolve(ksp, incoarse_vec, outcoarse_vec);
    copy_to_BlockVector(out_coarse, outcoarse_vec);

    // Extend the solution
    solution_interpolate_dof(dof_handler_coarse, dof_handler,
      transfer_mat_prolongate, out_coarse, out);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE<dim, n_fe_degree>::apply_gmg (unsigned int ngr,
    unsigned int ngc,
    PETScWrappers::MPI::Vector &in,
    PETScWrappers::MPI::Vector &out)
  {
    Assert(ngr==ngc, ExcMessage("ngr must be equal to ngc"));

    (void) ngc;

    PETScWrappers::MPI::Vector in_coarse;
    PETScWrappers::MPI::Vector out_coarse;
    in_coarse.reinit(MPI_COMM_WORLD, ndofs_coarse, ndofs_coarse);
    out_coarse.reinit(MPI_COMM_WORLD, ndofs_coarse, ndofs_coarse);
    out.reinit(MPI_COMM_WORLD, dof_handler.n_dofs(), dof_handler.n_dofs());

    Vec incoarse_vec, outcoarse_vec;
    VecCreateSeq(MPI_COMM_WORLD, ndofs_coarse, &incoarse_vec);
    VecDuplicate(incoarse_vec, &outcoarse_vec);

    // Restrict the solution
    solution_interpolate_dof(dof_handler, dof_handler_coarse,
      transfer_mat_restrict, in, in_coarse);

    // Solve the linear system in the linear elements
    VecAssemblyBegin(in_coarse);
    VecAssemblyEnd(in_coarse);
    VecCopy(in_coarse, incoarse_vec);
    VecAssemblyBegin(incoarse_vec);
    VecAssemblyEnd(incoarse_vec);
    KSPSolve(ksp_blocks[ngr], incoarse_vec, outcoarse_vec);
    VecCopy(outcoarse_vec, out_coarse);

    // Prolongate the solution
    solution_interpolate_dof(dof_handler_coarse, dof_handler,
      transfer_mat_prolongate, out_coarse, out);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE<dim, n_fe_degree>::ksp_setup (KSP &ksp)
  {
    PC pc;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, shellD, shellD);
    KSPSetType(ksp, KSPGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
    KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 100);
    //KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE<dim, n_fe_degree>::ksp_block_setup (KSP &ksp,
    unsigned int nb)
  {
    PC pc;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetTolerances(ksp, 1e-5, PETSC_DEFAULT,
    PETSC_DEFAULT, 15);
    KSPSetOperators(ksp, L.block(nb, nb), L.block(nb, nb));
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCICC);
    PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
    PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
    //KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
  }

// ----------- Explicit Instantations ----------- //
template class PC_MLFE<1, 1> ;
template class PC_MLFE<1, 2> ;
template class PC_MLFE<1, 3> ;
template class PC_MLFE<1, 4> ;
template class PC_MLFE<1, 5> ;

template class PC_MLFE<2, 1> ;
template class PC_MLFE<2, 2> ;
template class PC_MLFE<2, 3> ;
template class PC_MLFE<2, 4> ;
template class PC_MLFE<2, 5> ;

template class PC_MLFE<3, 1> ;
template class PC_MLFE<3, 2> ;
template class PC_MLFE<3, 3> ;
template class PC_MLFE<3, 4> ;
template class PC_MLFE<3, 5> ;

template <int dim, int n_fe_degree>
  void shell_D_coarse_gmg (Mat shell_mat,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    PC_MLFE<dim, n_fe_degree> *EPSobject = (PC_MLFE<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(EPSobject->n_blocks, MPI_COMM_WORLD,
      EPSobject->ndofs_coarse, EPSobject->ndofs_coarse);
    dstblock.reinit(EPSobject->n_blocks, MPI_COMM_WORLD,
      EPSobject->ndofs_coarse, EPSobject->ndofs_coarse);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)
    EPSobject->M.vmult(dstblock, srcblock);
    dstblock *= -EPSobject->lambda;
    EPSobject->L.vmult_add(dstblock, srcblock);

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < EPSobject->n_blocks; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PC_MLFE_Time<dim, n_fe_degree>::PC_MLFE_Time (const MPI_Comm &_comm,
    const TransportMatrixBase<dim, n_fe_degree> &_T,
    const Materials &_materials) :
      comm(_comm),
      T(_T),
      materials(_materials),
      dof_handler(T.dof_handler),
      dof_handler_coarse(dof_handler.get_triangulation()),
      T_coarse_1(_comm, dof_handler_coarse, constraints_coarse),
      T_coarse_2(_comm, dof_handler_coarse, constraints_coarse),
      T_coarse_3(_comm, dof_handler_coarse, constraints_coarse),
      T_coarse_spn_1(_comm, dof_handler_coarse, constraints_coarse),
      T_coarse_spn_2(_comm, dof_handler_coarse, constraints_coarse)
  {
    // Silly initializations
    n_blocks = 0;
    ksp = NULL;
    n_dofs_coarse = 0;
    n_blocks_coarse = 0;
    n_dofs_coarse_local = 0;
    lambda = 1.00;
    shellT_coarse = NULL;

    total_its = 0;
    n_applications = 0;

    time_scheme = "implicit";
    matrixfree_type_coarse = full_matrixfree;
    delta_t = 0.1;

    first_computation = true;
    n_fe_coarse = 1;
    get_uint_from_options("-n_fe_coarse", n_fe_coarse);
    transport_coarse = T.type_approximation;
    get_string_from_options("-transport_coarse", transport_coarse);
    transport_fine = T.type_approximation;
    pc_coarse = "diagonal";
    get_string_from_options("-pc_coarse", pc_coarse);
    get_enum_from_options("-matrixfree_type_coarse", matrixfree_type_coarse);

    if (pc_coarse == "gs-cgilu")
      matrixfree_type_coarse = non_diagonal;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::reinit ()
  {
    delta_t = T.delta_t;
    time_scheme = T.type_scheme;

    if (!first_computation)
    {
      if (transport_coarse == "diffusion")
      {
        if (n_fe_coarse == 1)
          T_coarse_1.clear();
        else if (n_fe_coarse == 2)
          T_coarse_2.clear();
        else if (n_fe_coarse == 3)
          T_coarse_3.clear();
      }

      else if (transport_coarse == "spn")
      {
        if (n_fe_coarse == 1)
          T_coarse_spn_1.clear();
        else if (n_fe_coarse == 2)
          T_coarse_spn_2.clear();
      }
    }

    if (first_computation)
    {
      boundary_conditions = T.boundary_conditions;
      albedo_factors = T.albedo_factors;
      time_scheme = T.type_scheme;

      n_blocks = T.n_blocks_rows();
      if (transport_coarse == transport_fine)
        n_blocks_coarse = T.n_blocks_rows();
      else
        n_blocks_coarse = T.n_blocks_rows() / 2;

      FE_Q<dim> fe_coarse(QGaussLobatto<1>(n_fe_coarse + 1));

      dof_handler_coarse.distribute_dofs(fe_coarse);

      n_dofs_coarse = dof_handler_coarse.n_dofs();
      n_dofs_coarse_local =
                            dof_handler_coarse.locally_owned_dofs().n_elements();

      constraints_coarse.clear();
      for (unsigned int c = 0; c < boundary_conditions.size(); c++)
        if (boundary_conditions[c] == 0)
        {
          DoFTools::make_zero_boundary_constraints(dof_handler_coarse, c,
            constraints_coarse);
        }
      constraints_coarse.close();

      // Set the transfer matrices
      transfer_mat_restrict.reinit(fe_coarse.dofs_per_cell,
        dof_handler.get_fe().dofs_per_cell);
      transfer_mat_prolongate.reinit(dof_handler.get_fe().dofs_per_cell,
        fe_coarse.dofs_per_cell);

      set_transfer_matrix(dof_handler.get_fe(), fe_coarse,
        transfer_mat_restrict);
      set_transfer_matrix(fe_coarse, dof_handler.get_fe(),
        transfer_mat_prolongate);

      MatCreateShell(comm, n_blocks_coarse * n_dofs_coarse_local,
        n_blocks_coarse * n_dofs_coarse_local,
        n_blocks_coarse * n_dofs_coarse,
        n_blocks_coarse * n_dofs_coarse, this, &shellT_coarse);

      MatShellSetOperation(shellT_coarse, MATOP_MULT,
        (void (*) ()) shell_T_coarse<dim, n_fe_degree>);
      ;
      first_computation = false;
    }

    // Set the shell coarse T
    if (transport_coarse == "diffusion")
    {
      if (n_fe_coarse == 1)
        T_coarse_1.reinit(materials, boundary_conditions, albedo_factors,
          delta_t, time_scheme, matrixfree_type_coarse, false);
      else if (n_fe_coarse == 2)
        T_coarse_2.reinit(materials, boundary_conditions, albedo_factors,
          delta_t, time_scheme, matrixfree_type_coarse, false);
      else if (n_fe_coarse == 3)
        T_coarse_3.reinit(materials, boundary_conditions, albedo_factors,
          delta_t, time_scheme, matrixfree_type_coarse, false);
    }

    else if (transport_coarse == "spn")
    {
      if (n_fe_coarse == 1)
        T_coarse_spn_1.reinit(materials, T.get_n_moments(),
          boundary_conditions, albedo_factors, delta_t, time_scheme,
          matrixfree_type_coarse, false);
      else if (n_fe_coarse == 2)
        T_coarse_spn_2.reinit(materials, T.get_n_moments(),
          boundary_conditions, albedo_factors, delta_t, time_scheme,
          matrixfree_type_coarse, false);
    }

    // Setup the solver
    KSPDestroy(&ksp);

    if (pc_coarse == "diagonal")
      pc_setup();
    else if (pc_coarse == "gs-cgilu" and ksp_blocks.size() < 1)
      pc_gs_cgilu_setup();

    ksp_setup();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::apply_gmg (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {
    AssertRelease(Utilities::MPI::n_mpi_processes(comm) == 1,
      "The function solution_interpolate_dof must be fixed for parallel computing");

    PETScWrappers::MPI::BlockVector in_coarse;
    PETScWrappers::MPI::BlockVector out_coarse;
    in_coarse.reinit(n_blocks_coarse, comm, n_dofs_coarse, n_dofs_coarse_local);
    out_coarse.reinit(n_blocks_coarse, comm, n_dofs_coarse,
      n_dofs_coarse_local);

    Vec incoarse_vec, outcoarse_vec;
    VecCreateMPI(comm, n_dofs_coarse_local * n_blocks_coarse,
      n_dofs_coarse * n_blocks_coarse, &incoarse_vec);
    VecCreateMPI(comm, n_dofs_coarse_local * n_blocks_coarse,
      n_dofs_coarse * n_blocks_coarse, &outcoarse_vec);

    // Restrict the solution
    restrict_solution(in, in_coarse);
    //	solution_interpolate_dof(dof_handler, dof_handler_coarse,
    //			transfer_mat_restrict, in, in_coarse);

    // Solve the linear system in the linear elements
    copy_to_Vec(incoarse_vec, in_coarse);
    KSPSolve(ksp, incoarse_vec, outcoarse_vec);
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    total_its += its;
    n_applications++;
    std::cout << "its" << its << std::endl;
    copy_to_BlockVector(out_coarse, outcoarse_vec);

    // Extend the solution
    prolongate_solution(out_coarse, out);
    //	solution_interpolate_dof(dof_handler_coarse, dof_handler,
    //			transfer_mat_prolongate, out_coarse, out);

    VecDestroy(&incoarse_vec);
    VecDestroy(&outcoarse_vec);
    for (unsigned nb = 0; nb < n_blocks_coarse; nb++)
    {
      in_coarse.block(nb).clear();
      out_coarse.block(nb).clear();
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::restrict_solution (
    PETScWrappers::MPI::BlockVector &in,
    PETScWrappers::MPI::BlockVector &in_coarse)
  {

    if (transport_fine == transport_coarse)
    {
      solution_interpolate_dof(dof_handler, dof_handler_coarse,
        transfer_mat_restrict, in, in_coarse);
    }
    else if (n_fe_degree == n_fe_coarse)
    {
      for (unsigned int b = 0; b < n_blocks_coarse; b++)
        in_coarse.block(b) = in.block(2 * b);
    }
    else
    {
      PETScWrappers::MPI::BlockVector innervec;
      innervec.reinit(n_blocks_coarse, comm, dof_handler.n_dofs(), dof_handler.n_dofs());

      in.compress(VectorOperation::insert);

      // Restrict the approximation
      for (unsigned int b = 0; b < n_blocks_coarse; b++)
      {
        innervec.block(b) = in.block(2 * b);
        innervec.block(b).compress(VectorOperation::insert);
      }

      // Restrict the fed
      solution_interpolate_dof(dof_handler, dof_handler_coarse,
        transfer_mat_restrict, innervec, in_coarse);

      for (unsigned int b = 0; b < n_blocks_coarse; b++)
        innervec.block(b).clear();
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::prolongate_solution (
    PETScWrappers::MPI::BlockVector &out_coarse,
    PETScWrappers::MPI::BlockVector &out)
  {

    if (transport_fine == transport_coarse)
    {
      solution_interpolate_dof(dof_handler_coarse, dof_handler,
        transfer_mat_prolongate, out_coarse, out);
    }
    else if (n_fe_degree == n_fe_coarse)
    {

      out *= 0.0;
      for (unsigned int b = 0; b < n_blocks_coarse; b++)
        out.block(2 * b) = out_coarse.block(b);
    }
    else
    {
      PETScWrappers::MPI::BlockVector innervec;
      innervec.reinit(n_blocks, comm, dof_handler_coarse.n_dofs(),
        dof_handler_coarse.n_dofs());

      // Restrict the approximation
      for (unsigned int b = 0; b < n_blocks_coarse; b++)
        innervec.block(2 * b) = out_coarse.block(b);

      // Restrict the fed
      solution_interpolate_dof(dof_handler_coarse, dof_handler,
        transfer_mat_prolongate, innervec, out);

      for (unsigned int b = 0; b < n_blocks_coarse; b++)
        innervec.block(b).clear();

    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::ksp_setup ()
  {
    PC pc;
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, shellT_coarse, shellT_coarse);
    KSPSetType(ksp, KSPGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);
    if (pc_coarse == "diagonal")
      PCShellSetApply(pc, shell_pc_T_coarse<dim, n_fe_degree>);
    else if (pc_coarse == "gs-cgilu")
      PCShellSetApply(pc, shell_pc_T_coarse_gs_cgilu<dim, n_fe_degree>);
    PCShellSetContext(pc, this);
    KSPSetTolerances(ksp, 1e-5, 1e-5, PETSC_DEFAULT, 100);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::pc_setup ()
  {
    PETScWrappers::MPI::BlockVector diagonal_inverse;
    diagonal_inverse.reinit(n_blocks_coarse, comm, n_dofs_coarse,
      n_dofs_coarse_local);

    if (transport_coarse == "diffusion")
    {
      if (n_fe_coarse == 1)
        T_coarse_1.get_inv_diagonal(diagonal_inverse);
      else if (n_fe_coarse == 2)
        T_coarse_2.get_inv_diagonal(diagonal_inverse);
      else if (n_fe_coarse == 3)
        T_coarse_3.get_inv_diagonal(diagonal_inverse);
    }

    else if (transport_coarse == "spn")
    {
      if (n_fe_coarse == 1)
        T_coarse_spn_1.get_inv_diagonal(diagonal_inverse);
      else if (n_fe_coarse == 2)
        T_coarse_spn_2.get_inv_diagonal(diagonal_inverse);
    }

    prec_diag.reinit(diagonal_inverse);

    for (unsigned int nb = 0; nb < n_blocks_coarse; nb++)
      diagonal_inverse.block(nb).clear();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Time<dim, n_fe_degree>::pc_gs_cgilu_setup ()
  {

    ksp_blocks.resize(n_blocks_coarse);
    pc_blocks.resize(n_blocks_coarse);
    KSP *subksp;
    PC subpc;
    PetscInt n_local, p;

    // Set up the GS Preconditioner
    for (unsigned int i = 0; i < n_blocks_coarse; ++i)
    {
      KSPCreate(comm, &(ksp_blocks[i]));
      KSPSetType(ksp_blocks[i], KSPCG);
      KSPSetTolerances(ksp_blocks[i], 1e-3, 1e-3,
      PETSC_DEFAULT, 10);

      if (transport_coarse == "diffusion")
      {
        if (n_fe_coarse == 1)
          KSPSetOperators(ksp_blocks[i], T_coarse_1.block(i, i), T_coarse_1.block(i, i));
        else if (n_fe_coarse == 2)
          KSPSetOperators(ksp_blocks[i], T_coarse_2.block(i, i), T_coarse_2.block(i, i));
        else if (n_fe_coarse == 3)
          KSPSetOperators(ksp_blocks[i], T_coarse_3.block(i, i), T_coarse_3.block(i, i));
      }

      else if (transport_coarse == "spn")
      {
        if (n_fe_coarse == 1)
          KSPSetOperators(ksp_blocks[i], T_coarse_spn_1.block(i, i),
            T_coarse_spn_1.block(i, i));
        else if (n_fe_coarse == 2)
          KSPSetOperators(ksp_blocks[i], T_coarse_spn_2.block(i, i),
            T_coarse_spn_2.block(i, i));
      }

      KSPGetPC(ksp_blocks[i], &pc_blocks[i]);
      PCSetType(pc_blocks[i], PCBJACOBI);

      PCFactorSetShiftType(pc_blocks[i], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pc_blocks[i], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(ksp_blocks[i], PETSC_TRUE);
      KSPSetFromOptions(ksp_blocks[i]);
      KSPSetNormType(ksp_blocks[i], KSP_NORM_UNPRECONDITIONED);
      KSPSetUp(ksp_blocks[i]);

      PCBJacobiGetSubKSP(pc_blocks[i], &n_local, NULL, &subksp);
      for (p = 0; p < n_local; p++)
      {
        KSPSetType(subksp[p], KSPPREONLY);
        KSPGetPC(subksp[p], &subpc);
        PCSetType(subpc, PCILU);
        PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
        PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
      }
    }

    return;
  }

// ----------- Explicit Instantations ----------- //
template class PC_MLFE_Time<1, 1> ;
template class PC_MLFE_Time<1, 2> ;
template class PC_MLFE_Time<1, 3> ;
template class PC_MLFE_Time<1, 4> ;
template class PC_MLFE_Time<1, 5> ;

template class PC_MLFE_Time<2, 1> ;
template class PC_MLFE_Time<2, 2> ;
template class PC_MLFE_Time<2, 3> ;
template class PC_MLFE_Time<2, 4> ;
template class PC_MLFE_Time<2, 5> ;

template class PC_MLFE_Time<3, 1> ;
template class PC_MLFE_Time<3, 2> ;
template class PC_MLFE_Time<3, 3> ;
template class PC_MLFE_Time<3, 4> ;
template class PC_MLFE_Time<3, 5> ;

template <int dim, int n_fe_degree>
  void shell_T_coarse (Mat shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    PC_MLFE_Time<dim, n_fe_degree> *PCobject =
                                               (PC_MLFE_Time<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    dstblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)

    if (PCobject->transport_coarse == "diffusion")
    {
      if (PCobject->n_fe_coarse == 1)
        PCobject->T_coarse_1.vmult(dstblock, srcblock);
      else if (PCobject->n_fe_coarse == 2)
        PCobject->T_coarse_2.vmult(dstblock, srcblock);
      else if (PCobject->n_fe_coarse == 3)
        PCobject->T_coarse_3.vmult(dstblock, srcblock);
    }

    else if (PCobject->transport_coarse == "spn")
    {
      if (PCobject->n_fe_coarse == 1)
        PCobject->T_coarse_spn_1.vmult(dstblock, srcblock);
      else if (PCobject->n_fe_coarse == 2)
        PCobject->T_coarse_spn_2.vmult(dstblock, srcblock);
    }

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < PCobject->n_blocks_coarse; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_T_coarse (PC shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    PCShellGetContext(shell_mat, &ctx);
    PC_MLFE_Time<dim, n_fe_degree> *PCobject =
                                               (PC_MLFE_Time<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    dstblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)
    PCobject->prec_diag.vmult(dstblock, srcblock);

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < PCobject->n_blocks_coarse; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return 0;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_T_coarse_gs_cgilu (PC shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    PCShellGetContext(shell_mat, &ctx);
    PC_MLFE_Time<dim, n_fe_degree> *PCobject =
                                               (PC_MLFE_Time<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    dstblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)
    ///////////////////////////
    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(PCobject->comm, PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    vecacc.reinit(PCobject->comm, PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    // Compute x1
    KSPSolve(PCobject->ksp_blocks[0], srcblock.block(0), dstblock.block(0));
    // Compute x2..
    for (unsigned int ng = 1; ng < PCobject->n_blocks_coarse; ng++)
    {
      vecacc = srcblock.block(ng);
      for (unsigned int subg = 0; subg < ng; subg++)
      {
        if (PCobject->transport_coarse == "diffusion")
        {
          if (PCobject->n_fe_coarse == 1)
            PCobject->T_coarse_1.vmult(ng, subg, inter1, dstblock.block(subg));
          else if (PCobject->n_fe_coarse == 2)
            PCobject->T_coarse_2.vmult(ng, subg, inter1, dstblock.block(subg));
          else if (PCobject->n_fe_coarse == 3)
            PCobject->T_coarse_3.vmult(ng, subg, inter1, dstblock.block(subg));
        }

        else if (PCobject->transport_coarse == "spn")
        {
          if (PCobject->n_fe_coarse == 1)
            PCobject->T_coarse_spn_1.vmult(ng, subg, inter1, dstblock.block(subg));
          else if (PCobject->n_fe_coarse == 2)
            PCobject->T_coarse_spn_2.vmult(ng, subg, inter1, dstblock.block(subg));
        }
        VecAXPY(vecacc, -1.0, inter1);
      }
      KSPSolve(PCobject->ksp_blocks[ng], vecacc, dstblock.block(ng));
    }

    inter1.clear();
    vecacc.clear();
    //////////////////////////////

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < PCobject->n_blocks_coarse; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return 0;
  }

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// FIXME
/**
 *
 */
template <int dim, int n_fe_degree>
  PC_MLFE_Complex<dim, n_fe_degree>::PC_MLFE_Complex (const MPI_Comm &_comm,
    const TransportMatrixComplexBase<dim, n_fe_degree> &_A,
    const Materials &_materials,
    const ComplexPerturbation &_pert) :
      comm(_comm),
      A(_A),
      materials(_materials),
      pert(_pert),
      dof_handler(A.dof_handler),
      dof_handler_coarse(dof_handler.get_triangulation()),
      A_coarse_1(_comm, dof_handler_coarse, constraints_coarse),
      A_coarse_2(_comm, dof_handler_coarse, constraints_coarse),
      A_coarse_3(_comm, dof_handler_coarse, constraints_coarse)
  //      A_coarse_spn_1(_comm, dof_handler_coarse, constraints_coarse),
//      A_coarse_spn_2(_comm, dof_handler_coarse, constraints_coarse)
  {
    // Silly initializations
    n_blocks = 0;
    ksp = NULL;
    n_dofs_coarse = 0;
    n_blocks_coarse = 0;
    n_dofs_coarse_local = 0;
    lambda = 1.00;
    shellA_coarse = NULL;

    total_its = 0;
    n_applications = 0;

    first_computation = true;
    n_fe_coarse = 1;
    get_uint_from_options("-n_fe_coarse", n_fe_coarse);
    transport_coarse = "diffusion";
    //get_string_from_options("-transport_coarse", transport_coarse);
    //transport_fine = A.type_approximation;
    pc_coarse = "diagonal";
    get_string_from_options("-pc_coarse", pc_coarse);
    get_enum_from_options("-matrixfree_type_coarse", matrixfree_type_coarse);

    matrixfree_type_coarse = full_matrixfree;
    if (pc_coarse == "gs-cgilu")
      matrixfree_type_coarse = non_diagonal;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::reinit ()
  {

    if (!first_computation)
    {
      // if (transport_coarse == "diffusion")
      // {
      if (n_fe_coarse == 1)
        A_coarse_1.clear();
      else if (n_fe_coarse == 2)
        A_coarse_2.clear();
      else if (n_fe_coarse == 3)
        A_coarse_3.clear();
      //}

      //      else if (transport_coarse == "spn")
      //      {
      //        if (n_fe_coarse == 1)
      //          A_coarse_spn_1.clear();
      //        else if (n_fe_coarse == 2)
      //          A_coarse_spn_2.clear();
      //      }
    }

    if (first_computation)
    {
      boundary_conditions = A.boundary_conditions;
      albedo_factors = A.albedo_factors;

      n_blocks = A.n_blocks_rows();
      //if (transport_coarse == transport_fine)
      n_blocks_coarse = A.n_blocks_rows();

      FE_Q<dim> fe_coarse(QGaussLobatto<1>(n_fe_coarse + 1));

      dof_handler_coarse.distribute_dofs(fe_coarse);

      n_dofs_coarse = dof_handler_coarse.n_dofs();
      n_dofs_coarse_local = dof_handler_coarse.locally_owned_dofs().n_elements();

      constraints_coarse.clear();
      for (unsigned int c = 0; c < boundary_conditions.size(); c++)
        if (boundary_conditions[c] == 0)
        {
          DoFTools::make_zero_boundary_constraints(dof_handler_coarse, c,
            constraints_coarse);
        }
      constraints_coarse.close();

      // Set the transfer matrices
      transfer_mat_restrict.reinit(fe_coarse.dofs_per_cell,
        dof_handler.get_fe().dofs_per_cell);
      transfer_mat_prolongate.reinit(dof_handler.get_fe().dofs_per_cell,
        fe_coarse.dofs_per_cell);

      set_transfer_matrix(dof_handler.get_fe(), fe_coarse,
        transfer_mat_restrict);
      set_transfer_matrix(fe_coarse, dof_handler.get_fe(),
        transfer_mat_prolongate);

      MatCreateShell(comm, n_blocks_coarse * n_dofs_coarse_local,
        n_blocks_coarse * n_dofs_coarse_local,
        n_blocks_coarse * n_dofs_coarse,
        n_blocks_coarse * n_dofs_coarse, this, &shellA_coarse);

      MatShellSetOperation(shellA_coarse, MATOP_MULT,
        (void (*) ()) shell_A_coarse<dim, n_fe_degree>);
      ;
      first_computation = false;

      std::cout << "size small problem: " << n_blocks_coarse * n_dofs_coarse_local
                << std::endl;
    }

    // Set the shell coarse A
    //if (transport_coarse == "diffusion")
    //{
    if (n_fe_coarse == 1)
      A_coarse_1.reinit(materials, pert, boundary_conditions, albedo_factors,
        matrixfree_type_coarse);
    else if (n_fe_coarse == 2)
      A_coarse_2.reinit(materials, pert, boundary_conditions, albedo_factors,
        matrixfree_type_coarse);
    else if (n_fe_coarse == 3)
      A_coarse_3.reinit(materials, pert, boundary_conditions, albedo_factors,
        matrixfree_type_coarse);

    //}

    //    else if (transport_coarse == "spn")
    //    {
    //      if (n_fe_coarse == 1)
    //        A_coarse_spn_1.reinit(materials, T.get_n_moments(),
    //          boundary_conditions, albedo_factors, delta_t, time_scheme,
    //          matrixfree_type_coarse, false);
    //      else if (n_fe_coarse == 2)
    //        A_coarse_spn_2.reinit(materials, T.get_n_moments(),
    //          boundary_conditions, albedo_factors, delta_t, time_scheme,
    //          matrixfree_type_coarse, false);
    //    }

    // Setup the solver
    KSPDestroy(&ksp);

    if (pc_coarse == "diagonal")
      pc_setup();
    else if (pc_coarse == "gs-cgilu" and ksp_blocks.size() < 1)
      pc_gs_cgilu_setup();

    ksp_setup();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::apply_gmg (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {
    AssertRelease(Utilities::MPI::n_mpi_processes(comm) == 1,
      "The function solution_interpolate_dof must be fixed for parallel computing");

    PETScWrappers::MPI::BlockVector in_coarse;
    PETScWrappers::MPI::BlockVector out_coarse;
    in_coarse.reinit(n_blocks_coarse, comm, n_dofs_coarse, n_dofs_coarse_local);
    out_coarse.reinit(n_blocks_coarse, comm, n_dofs_coarse, n_dofs_coarse_local);

    Vec incoarse_vec, outcoarse_vec;
    VecCreateMPI(comm, n_dofs_coarse_local * n_blocks_coarse,
      n_dofs_coarse * n_blocks_coarse, &incoarse_vec);
    VecCreateMPI(comm, n_dofs_coarse_local * n_blocks_coarse,
      n_dofs_coarse * n_blocks_coarse, &outcoarse_vec);

    // Restrict the solution
    in_coarse.compress(VectorOperation::insert);
    in.compress(VectorOperation::insert);
    VecAssemblyBegin(incoarse_vec);
    VecAssemblyEnd(incoarse_vec);

    restrict_solution(in, in_coarse);

    //  solution_interpolate_dof(dof_handler, dof_handler_coarse,
    //      transfer_mat_restrict, in, in_coarse);

    // Solve the linear system in the linear elements
    copy_to_Vec(incoarse_vec, in_coarse);
    KSPSolve(ksp, incoarse_vec, outcoarse_vec);
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    std::cout << "    Coarse problem n its: " << its << std::endl;

    total_its += its;
    n_applications++;

    copy_to_BlockVector(out_coarse, outcoarse_vec);

    // Extend the solution

    prolongate_solution(out_coarse, out);
    //  solution_interpolate_dof(dof_handler_coarse, dof_handler,
    //      transfer_mat_prolongate, out_coarse, out);

    VecDestroy(&incoarse_vec);
    VecDestroy(&outcoarse_vec);
    for (unsigned nb = 0; nb < n_blocks_coarse; nb++)
    {
      in_coarse.block(nb).clear();
      out_coarse.block(nb).clear();
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::restrict_solution (
    PETScWrappers::MPI::BlockVector &in,
    PETScWrappers::MPI::BlockVector &in_coarse)
  {

    /// if (transport_fine == transport_coarse)
    // {

    solution_interpolate_dof(dof_handler, dof_handler_coarse,
      transfer_mat_restrict, in, in_coarse);
    //}
//    else if (n_fe_degree == n_fe_coarse)
//    {
//      for (unsigned int b = 0; b < n_blocks_coarse; b++)
//        in_coarse.block(b) = in.block(2 * b);
//    }
//    else
//    {
//      PETScWrappers::MPI::BlockVector innervec;
//      innervec.reinit(n_blocks_coarse, comm, dof_handler.n_dofs(), dof_handler.n_dofs());
//
//      in.compress(VectorOperation::insert);
//
//      // Restrict the approximation
//      for (unsigned int b = 0; b < n_blocks_coarse; b++)
//      {
//        innervec.block(b) = in.block(2 * b);
//        innervec.block(b).compress(VectorOperation::insert);
//      }
//
//      // Restrict the fed
//      solution_interpolate_dof(dof_handler, dof_handler_coarse,
//        transfer_mat_restrict, innervec, in_coarse);
//
//      for (unsigned int b = 0; b < n_blocks_coarse; b++)
//        innervec.block(b).clear();
//    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::prolongate_solution (
    PETScWrappers::MPI::BlockVector &out_coarse,
    PETScWrappers::MPI::BlockVector &out)
  {

    //if (transport_fine == transport_coarse)
    //{
    solution_interpolate_dof(dof_handler_coarse, dof_handler,
      transfer_mat_prolongate, out_coarse, out);
    //}
//    else if (n_fe_degree == n_fe_coarse)
//    {
//
//      out *= 0.0;
//      for (unsigned int b = 0; b < n_blocks_coarse; b++)
//        out.block(2 * b) = out_coarse.block(b);
//    }
//    else
//    {
//      PETScWrappers::MPI::BlockVector innervec;
//      innervec.reinit(n_blocks, comm, dof_handler_coarse.n_dofs(),
//        dof_handler_coarse.n_dofs());
//
//      // Restrict the approximation
//      for (unsigned int b = 0; b < n_blocks_coarse; b++)
//        innervec.block(2 * b) = out_coarse.block(b);
//
//      // Restrict the fed
//      solution_interpolate_dof(dof_handler_coarse, dof_handler,
//        transfer_mat_prolongate, innervec, out);
//
//      for (unsigned int b = 0; b < n_blocks_coarse; b++)
//        innervec.block(b).clear();
//
//    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::ksp_setup ()
  {
    PC pc;
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, shellA_coarse, shellA_coarse);
    KSPSetType(ksp, KSPFGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);
    if (pc_coarse == "diagonal")
      PCShellSetApply(pc, shell_pc_A_coarse<dim, n_fe_degree>);
    else if (pc_coarse == "gs-cgilu")
      PCShellSetApply(pc, shell_pc_A_coarse_gs_cgilu<dim, n_fe_degree>);
    PCShellSetContext(pc, this);
    KSPSetTolerances(ksp, 1e-3, 1e-3, PETSC_DEFAULT, 100);
    //KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::pc_setup ()
  {
    PETScWrappers::MPI::BlockVector diagonal_inverse;
    diagonal_inverse.reinit(n_blocks_coarse, comm, n_dofs_coarse,
      n_dofs_coarse_local);

    //if (transport_coarse == "diffusion")
    //{
    if (n_fe_coarse == 1)
      A_coarse_1.get_inv_diagonal(diagonal_inverse);
    else if (n_fe_coarse == 2)
      A_coarse_2.get_inv_diagonal(diagonal_inverse);
    else if (n_fe_coarse == 3)
      A_coarse_3.get_inv_diagonal(diagonal_inverse);
    //}
    //    else if (transport_coarse == "spn")
    //    {
    //      if (n_fe_coarse == 1)
    //        A_coarse_spn_1.get_inv_diagonal(diagonal_inverse);
    //      else if (n_fe_coarse == 2)
    //        A_coarse_spn_2.get_inv_diagonal(diagonal_inverse);
    //    }

    prec_diag.reinit(diagonal_inverse);

    for (unsigned int nb = 0; nb < n_blocks_coarse; nb++)
      diagonal_inverse.block(nb).clear();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void PC_MLFE_Complex<dim, n_fe_degree>::pc_gs_cgilu_setup ()
  {

    // FIXME Se necesitan la mitat de ksp blocks
    ksp_blocks.resize(n_blocks_coarse);
    pc_blocks.resize(n_blocks_coarse);
    KSP *subksp;
    PC subpc;
    PetscInt n_local, p;

    // Set up the GS Preconditioner
    for (unsigned int i = 0; i < n_blocks_coarse; ++i)
    {
      KSPCreate(comm, &(ksp_blocks[i]));
      KSPSetType(ksp_blocks[i], KSPCG);
      KSPSetTolerances(ksp_blocks[i], 1e-3, 1e-3,
      PETSC_DEFAULT, 10);

//      if (transport_coarse == "diffusion")
//      {
      if (n_fe_coarse == 1)
        KSPSetOperators(ksp_blocks[i], A_coarse_1.block(i, i), A_coarse_1.block(i, i));
      else if (n_fe_coarse == 2)
        KSPSetOperators(ksp_blocks[i], A_coarse_2.block(i, i), A_coarse_2.block(i, i));
      else if (n_fe_coarse == 3)
        KSPSetOperators(ksp_blocks[i], A_coarse_3.block(i, i), A_coarse_3.block(i, i));
//      }
//
//      else if (transport_coarse == "spn")
//      {
//        if (n_fe_coarse == 1)
//          KSPSetOperators(ksp_blocks[i], A_coarse_spn_1.block(i, i),
//            A_coarse_spn_1.block(i, i));
//        else if (n_fe_coarse == 2)
//          KSPSetOperators(ksp_blocks[i], A_coarse_spn_2.block(i, i),
//            A_coarse_spn_2.block(i, i));
//      }

      KSPGetPC(ksp_blocks[i], &pc_blocks[i]);
      PCSetType(pc_blocks[i], PCBJACOBI);

      PCFactorSetShiftType(pc_blocks[i], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pc_blocks[i], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(ksp_blocks[i], PETSC_TRUE);
      (ksp_blocks[i]);
      KSPSetNormType(ksp_blocks[i], KSP_NORM_UNPRECONDITIONED);
      KSPSetUp(ksp_blocks[i]);

      PCBJacobiGetSubKSP(pc_blocks[i], &n_local, NULL, &subksp);
      for (p = 0; p < n_local; p++)
      {
        KSPSetType(subksp[p], KSPPREONLY);
        KSPGetPC(subksp[p], &subpc);
        PCSetType(subpc, PCILU);
        PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
        PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
      }
    }

    return;
  }

// ----------- Explicit Instantations ----------- //
template class PC_MLFE_Complex<1, 1> ;
template class PC_MLFE_Complex<1, 2> ;
template class PC_MLFE_Complex<1, 3> ;
template class PC_MLFE_Complex<1, 4> ;
template class PC_MLFE_Complex<1, 5> ;

template class PC_MLFE_Complex<2, 1> ;
template class PC_MLFE_Complex<2, 2> ;
template class PC_MLFE_Complex<2, 3> ;
template class PC_MLFE_Complex<2, 4> ;
template class PC_MLFE_Complex<2, 5> ;

template class PC_MLFE_Complex<3, 1> ;
template class PC_MLFE_Complex<3, 2> ;
template class PC_MLFE_Complex<3, 3> ;
template class PC_MLFE_Complex<3, 4> ;
template class PC_MLFE_Complex<3, 5> ;

template <int dim, int n_fe_degree>
  void shell_A_coarse (Mat shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    PC_MLFE_Complex<dim, n_fe_degree> *PCobject =
                                                  (PC_MLFE_Complex<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    dstblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)

    if (PCobject->transport_coarse == "diffusion")
    {
      if (PCobject->n_fe_coarse == 1)
        PCobject->A_coarse_1.vmult(dstblock, srcblock);
      else if (PCobject->n_fe_coarse == 2)
        PCobject->A_coarse_2.vmult(dstblock, srcblock);
      else if (PCobject->n_fe_coarse == 3)
        PCobject->A_coarse_3.vmult(dstblock, srcblock);
    }

//    else if (PCobject->transport_coarse == "spn")
//    {
//      if (PCobject->n_fe_coarse == 1)
//        PCobject->A_coarse_spn_1.vmult(dstblock, srcblock);
//      else if (PCobject->n_fe_coarse == 2)
//        PCobject->A_coarse_spn_2.vmult(dstblock, srcblock);
//    }

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < PCobject->n_blocks_coarse; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_A_coarse (PC shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    PCShellGetContext(shell_mat, &ctx);
    PC_MLFE_Complex<dim, n_fe_degree> *PCobject =
                                                  (PC_MLFE_Complex<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    dstblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)
    PCobject->prec_diag.vmult(dstblock, srcblock);

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < PCobject->n_blocks_coarse; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return 0;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode shell_pc_A_coarse_gs_cgilu (PC shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS object
    // so we can access the data of the problem.

    void *ctx;
    PCShellGetContext(shell_mat, &ctx);
    PC_MLFE_Complex<dim, n_fe_degree> *PCobject =
                                                  (PC_MLFE_Complex<dim, n_fe_degree>*) ctx;

    // Create Vectors
    PETScWrappers::MPI::BlockVector srcblock, dstblock;
    srcblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    dstblock.reinit(PCobject->n_blocks_coarse, PCobject->comm,
      PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    copy_to_BlockVector(srcblock, src);

    // Compute matrix (M-lambdaL)
    ///////////////////////////
    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(PCobject->comm, PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);
    vecacc.reinit(PCobject->comm, PCobject->n_dofs_coarse, PCobject->n_dofs_coarse_local);

    // Compute x1
    KSPSolve(PCobject->ksp_blocks[0], srcblock.block(0), dstblock.block(0));
    // Compute x2..
    for (unsigned int ng = 1; ng < PCobject->n_blocks_coarse; ng++)
    {
      vecacc = srcblock.block(ng);
      for (unsigned int subg = 0; subg < ng; subg++)
      {
//        if (PCobject->transport_coarse == "diffusion")
//        {
        if (PCobject->n_fe_coarse == 1)
          PCobject->A_coarse_1.vmult(ng, subg, inter1, dstblock.block(subg));
        else if (PCobject->n_fe_coarse == 2)
          PCobject->A_coarse_2.vmult(ng, subg, inter1, dstblock.block(subg));
        else if (PCobject->n_fe_coarse == 3)
          PCobject->A_coarse_3.vmult(ng, subg, inter1, dstblock.block(subg));
//        }
//
//        else if (PCobject->transport_coarse == "spn")
//        {
//          if (PCobject->n_fe_coarse == 1)
//            PCobject->A_coarse_spn_1.vmult(ng, subg, inter1, dstblock.block(subg));
//          else if (PCobject->n_fe_coarse == 2)
//            PCobject->A_coarse_spn_2.vmult(ng, subg, inter1, dstblock.block(subg));
//        }
        VecAXPY(vecacc, -1.0, inter1);
      }
      KSPSolve(PCobject->ksp_blocks[ng], vecacc, dstblock.block(ng));
    }

    inter1.clear();
    vecacc.clear();
    //////////////////////////////

    copy_to_Vec(dst, dstblock);

    for (unsigned int ng = 0; ng < PCobject->n_blocks_coarse; ng++)
    {
      srcblock.block(ng).clear();
      dstblock.block(ng).clear();
    }

    return 0;
  }

