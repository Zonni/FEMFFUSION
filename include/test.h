/**
 * @file   test.h
 * @brief  Automated Testing of the program in some  defined inputs.
 */

#ifndef TEST_H
#define TEST_H

using namespace dealii;

/**
 *  @brief Run all tests.
 */
int
run_tests ();

/**
 *  @brief Run a lambda modes test, only keff is checked.
 */
void
test_keff_problem (const std::string &input_file,
  const double keff,
  const double tol = 2e-5);

/*
 * @brief
 */
void
test_power_evolution (const std::string &input_file,
  const std::vector<double> &power,
  const double tol = 1e-4);

/*
 *  @brief Test noise problem through the mean and the l2 norm of the delta_phi.
 */
void test_noise_problem (
  const std::string &input_file,
  const double ref_l2norm_delta_phi,
  const double tol = 1e-2);

/**
 * @brief
 */
void chek_output_vector (
  const std::string &output_file,
  const std::string &headline,
  const std::vector<double> &vector_correct,
  const unsigned int n_lines = 1,
  const double tol = 1e-4);

/**
 * @brief test_sum_vector()
 */
void test_sum_vector ();

/**
 * @brief test_compute_eigenvalues()
 */
void test_compute_eigenvalues ();

/**
 * @brief  test_lower_case()
 */
void test_lower_case ();

/**
 *
 */
void test_copy_to ();
/**
 * @brief test_round
 */
void test_round ();

/**
 * run_tests
 */
void run_tests_utils ();

/**
 * @brief run_test_static_rom
 */
void run_test_static_rom (const std::string &input_file,
  unsigned int n_tests);

#endif

