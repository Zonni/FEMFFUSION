/**
 * @file   pc_multilevel.h
 * @brief
 */

#ifndef PERFORMANCE_H_
#define PERFORMANCE_H_

using namespace dealii;

/**
 *  Run all tests.
 */
int run_performance_matmult ();

int run_performance_solvers ();

/**
 * @brief
 */
template<int fe_degree>
double test_mf_performance_2D(
  const std::string& input_file,
  const unsigned int n_mults,
  const unsigned int block_row,
  const unsigned int block_col);

/*
 *
 */

template <int fe_degree>
  double valid_mf_2D (
    const std::string& input_file,
    const unsigned int block_row,
    const unsigned int block_col);

/**
 * @brief
 */
template<int fe_degree>
double test_mf_performance_3D(
  const std::string& input_file,
  const unsigned int n_mults,
  const unsigned int block_row,
  const unsigned int block_col);


#endif /* PERFORMANCE_H_ */
