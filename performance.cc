/**
 * @file   performance.cc
 * @brief Implementation of the neutron noise
 */

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/petsc_vector.h>

#include "femffusion.h"
#include "static_diffusion.h"
//#include "static_spn.h"
#include "test.h"
#include "matrix_operators/matrix_operators_petsc.h"
#include "performance.h"

using namespace dealii;

/**
 * @brief Run all test asserting the right values.
 */

int run_performance_matmult ()
{
  ConditionalOStream cout(std::cout,
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
  std::string input_file;

  const unsigned int n_vect_doubles = VectorizedArray<double>::n_array_elements;
  const unsigned int n_vect_bits = 8 * sizeof(double) * n_vect_doubles;
  cout << "  Vectorization over " << n_vect_doubles
       << " doubles = "
       << n_vect_bits << " bits ("
       << Utilities::System::get_current_vectorization_level()
       << "), VECTORIZATION_LEVEL=" << DEAL_II_COMPILER_VECTORIZATION_LEVEL
       << std::endl
       << std::endl;

  // --------------------------------------------------------------- //
  //  CHECK MATRIX_FREE MULTIPLICATION
  double norm_allocated, norm_matfree;


  cout << "CHECK MATRIX_FREE MULTIPLICATION... " << std::flush;
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_fullmatrixfree.prm";
  norm_matfree = valid_mf_2D<3>(input_file, 0, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  norm_allocated = valid_mf_2D<3>(input_file, 0, 0);


  AssertRelease(std::abs(norm_matfree - norm_allocated) < 1e-10,
    "Revise the fullmatrix free computation");

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_fullmatrixfree.prm";
  norm_matfree = valid_mf_2D<3>(input_file, 1, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  norm_allocated = valid_mf_2D<3>(input_file, 1, 0);

  AssertRelease(std::abs(norm_matfree - norm_allocated) < 1e-10,
    "Revise the fullmatrix free computation 2");


  cout << "DONE! " << std::endl;

  // --------------------------------------------------------------- //
  //  TEST MULTIPLICATION OF MATRICES
  double t_allocated, t_matfree;
  cout << "PERFORMANCE OF MATRIX-FREE MULTIPLICATIONS " << std::endl;

#ifdef DEBUG
  std::cout<<"Compiled in DEBUG mode!!"<<std::endl;
  exit(0);
#endif

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_fullmatrixfree.prm";
  t_matfree = test_mf_performance_2D<1>(input_file, 1e3, 0, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<1>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << std::fixed << "   L(0, 0) 2D BIBLIS with FE deg "
       << 1
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_fullmatrixfree.prm";
  t_matfree = test_mf_performance_2D<2>(input_file, 1e3, 0, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<2>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << "   L(0, 0) 2D BIBLIS with FE deg "
       << 2
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_fullmatrixfree.prm";
  t_matfree = test_mf_performance_2D<3>(input_file, 1e3, 0, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<3>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << "   L(0, 0) 2D BIBLIS with FE deg "
       << 3
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_fullmatrixfree.prm";
  t_matfree = test_mf_performance_2D<5>(input_file, 1e3, 0, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<5>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << "   L(0, 0) 2D BIBLIS with FE deg "
       << 5
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  cout << std::endl;
  // --------------------------------------------------------------- //
  // --------------------------------------------------------------- //
  // --------------------------------------------------------------- //
  // Test Matrix L_10
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_nondiagonal.prm";
  t_matfree = test_mf_performance_2D<1>(input_file, 1e3, 1, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<1>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 2D BIBLIS with FE deg "
       << 1
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_nondiagonal.prm";
  t_matfree = test_mf_performance_2D<2>(input_file, 1e3, 1, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<2>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 2D BIBLIS with FE deg "
       << 2
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_nondiagonal.prm";
  t_matfree = test_mf_performance_2D<3>(input_file, 1e3, 1, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<3>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 2D BIBLIS with FE deg "
       << 3
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE3_nondiagonal.prm";
  t_matfree = test_mf_performance_2D<5>(input_file, 1e3, 1, 0);

  input_file = "test/2D_biblis_matrixfree/biblis_FE3_allocated.prm";
  t_allocated = test_mf_performance_2D<5>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 2D BIBLIS with FE deg "
       << 5
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  cout << std::endl;

  // --------------------------------------------------------------- //
  // --------------------------------------------------------------- //
  // --------------------------------------------------------------- //
  // TEST 3D
  input_file = "test/3D_IAEA_test/IAEA_FE1_fullmatrixfree.prm";
  t_matfree = test_mf_performance_3D<1>(input_file, 1e3, 0, 0);
  input_file = "test/3D_IAEA_test/IAEA_FE1_allocated.prm";
  t_allocated = test_mf_performance_3D<1>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << "   L(0, 0) 3D IAEA with FE deg "
       << 1
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  // --------------------------------------------------------------- //
  input_file = "test/3D_IAEA_test/IAEA_FE1_fullmatrixfree.prm";
  t_matfree = test_mf_performance_3D<2>(input_file, 1e3, 0, 0);
  input_file = "test/3D_IAEA_test/IAEA_FE1_allocated.prm";
  t_allocated = test_mf_performance_3D<2>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << "   L(0, 0) 3D IAEA with FE deg "
       << 2
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/3D_IAEA_test/IAEA_FE1_fullmatrixfree.prm";
  t_matfree = test_mf_performance_3D<3>(input_file, 1e3, 0, 0);
  input_file = "test/3D_IAEA_test/IAEA_FE1_allocated.prm";
  t_allocated = test_mf_performance_3D<3>(input_file, 1e3, 0, 0);

  cout << std::setprecision(2) << "   L(0, 0) 3D IAEA with FE deg "
       << 3
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  cout << std::endl;
  //
  // --------------------------------------------------------------- //
  // TEST 3D L_10
  input_file = "test/3D_IAEA_test/IAEA_FE1_fullmatrixfree.prm";
  t_matfree = test_mf_performance_3D<1>(input_file, 1e3, 1, 0);
  input_file = "test/3D_IAEA_test/IAEA_FE1_allocated.prm";
  t_allocated = test_mf_performance_3D<1>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 3D IAEA with FE deg "
       << 1
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  // --------------------------------------------------------------- //
  input_file = "test/3D_IAEA_test/IAEA_FE1_fullmatrixfree.prm";
  t_matfree = test_mf_performance_3D<2>(input_file, 1e3, 1, 0);
  input_file = "test/3D_IAEA_test/IAEA_FE1_allocated.prm";
  t_allocated = test_mf_performance_3D<2>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 3D IAEA with FE deg "
       << 2
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;

  // --------------------------------------------------------------- //
  input_file = "test/3D_IAEA_test/IAEA_FE1_fullmatrixfree.prm";
  t_matfree = test_mf_performance_3D<3>(input_file, 1e3, 1, 0);
  input_file = "test/3D_IAEA_test/IAEA_FE1_allocated.prm";
  t_allocated = test_mf_performance_3D<3>(input_file, 1e3, 1, 0);

  cout << std::setprecision(2) << "   L(1, 0) 3D IAEA with FE deg "
       << 3
       << ": MF-> "
       << t_matfree
       << "s vs SpM-> "
       << t_allocated
       << "s. "
       << t_allocated / t_matfree
       << " times. " << std::endl;
  // --------------------------------------------------------------- //
  cout << std::endl;

  return 0;
}

int run_performance_solvers ()
{

  ConditionalOStream cout(std::cout,
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  cout << "PERFORMANCE OF SOLVERS NOT IMPLEMENTED" << std::endl;

  return 0;
}

/**
 *
 */
template <int fe_degree>
  double test_mf_performance_2D (
    const std::string& input_file,
    const unsigned int n_mults,
    const unsigned int block_row,
    const unsigned int block_col)
  {
    ParameterHandler prm;
    prm_declare_entries(prm);
    prm.parse_input(input_file);
    AssertRelease(2 == prm.get_integer("Dimension"), "Only dim==2 is allowed");
    AssertRelease(3 == prm.get_integer("FE_Degree"), "Only p==3 is allowed");
    std::string transport = prm.get("Transport_Appr");
    lower_case(transport);
    AssertRelease(transport == "diffusion", "Only diffusion is allowed");
    StaticDiffusion<2, fe_degree> problem(prm, input_file, false, true, true);

    problem.make_dofs();
    problem.assemble_system_lambda();

    PETScWrappers::MPI::Vector dst(problem.locally_owned_dofs, problem.comm);
    PETScWrappers::MPI::Vector src(problem.locally_owned_dofs, problem.comm);

    Timer timer;
    timer.start();

    for (unsigned int mult = 0; mult < n_mults; mult++)
    {
      problem.T.vmult(block_row, block_col, dst, src);
    }

    double time = timer.cpu_time();

    problem.F.clear();
    problem.T.clear();

    return time;
  }

/**
 *
 */
template <int fe_degree>
  double valid_mf_2D (
    const std::string& input_file,
    const unsigned int block_row,
    const unsigned int block_col)
  {
    ParameterHandler prm;
    prm_declare_entries(prm);
    prm.parse_input(input_file);
    AssertRelease(2 == prm.get_integer("Dimension"), "Only dim==2 is allowed");
    AssertRelease(3 == prm.get_integer("FE_Degree"), "Only p==3 is allowed");
    std::string transport = prm.get("Transport_Appr");
    lower_case(transport);
    AssertRelease(transport == "diffusion", "Only diffusion is allowed");
    StaticDiffusion<2, fe_degree> problem(prm, input_file, false, true, true);

    problem.make_dofs();
    problem.assemble_system_lambda();

    PETScWrappers::MPI::Vector dst(problem.locally_owned_dofs, problem.comm);
    PETScWrappers::MPI::Vector src(problem.locally_owned_dofs, problem.comm);

    src.add(1.0);

    problem.T.vmult(block_row, block_col, dst, src);

    double norm_L = dst.l2_norm();

    problem.F.clear();
    problem.T.clear();

    return norm_L;
  }

/**
 *
 */
template <int fe_degree>
  double test_mf_performance_3D (
    const std::string& input_file,
    const unsigned int n_mults,
    const unsigned int block_row,
    const unsigned int block_col)
  {
    ParameterHandler prm;
    prm_declare_entries(prm);
    prm.parse_input(input_file);
    AssertRelease(3 == prm.get_integer("Dimension"), "Only dim==3 is allowed");
    AssertRelease(1 == prm.get_integer("FE_Degree"), "Only p==1 is allowed");
    std::string transport = prm.get("Transport_Appr");
    lower_case(transport);
    AssertRelease(transport == "diffusion", "Only diffusion is allowed");

    StaticDiffusion<3, fe_degree> problem(prm, input_file, false, true, true);
    problem.make_dofs();
    problem.assemble_system_lambda();

    PETScWrappers::MPI::Vector dst(problem.locally_owned_dofs, problem.comm);
    PETScWrappers::MPI::Vector src(problem.locally_owned_dofs, problem.comm);

    Timer timer;
    timer.start();
    for (unsigned int mult = 0; mult < n_mults; mult++)
    {
      problem.T.vmult(block_row, block_col, dst, src);
    }
    double time = timer.cpu_time();
    problem.F.clear();
    problem.T.clear();

    return time;
  }

