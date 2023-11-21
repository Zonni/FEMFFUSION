/**
 * @file   printing.h
 * @brief  Functions for pinting the results in output files.
 */

#ifndef PRINTING_H_
#define PRINTING_H_

#include "materials.h"

/**
 * @brief Print a cell wise distribution in a file.
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<double>& matrix,
  const std::vector<unsigned int>& n_cells_per_dim,
  std::ofstream &out,
  const Materials & materials,
  const std::string &introduction,
  bool fill_with_zeros = false,
  unsigned int precision = 6);

/**
 * @brief The same as before but with a filename string.
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<double>& matrix,
  const std::vector<unsigned int>& n_cells_per_dim,
  const std::string &filename,
  const Materials & materials,
  const std::string &introduction,
  bool fill_with_zeros = false,
  unsigned int precision = 6);

/**
 * @brief rint a cell wise complex-valued distribution in a file.
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<std::complex<double> >& matrix,
  const std::vector<unsigned int>& n_cells_per_dim,
  std::ofstream &out,
  const Materials & materials,
  const std::string &introduction,
  bool fill_with_zeros = false,
  unsigned int precision = 6);

/**
 * @brief The same as before but with a filename string.
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<std::complex<double> >& matrix,
  const std::vector<unsigned int>& n_cells_per_dim,
  const std::string &filename,
  const Materials & materials,
  const std::string &introduction,
  bool fill_with_zeros = false,
  unsigned int precision = 6);

/**
 * @brief Print a vector in a file with an introduction.
 */
void print_vector_in_file (const std::vector<double> &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a vector in a file with an introduction.
 */
void print_vector_in_file (const std::vector<unsigned int> &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a vector in a file with an introduction.
 */
void print_vector_in_file (const PETScWrappers::MPI::Vector &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a vector in a file with an introduction.
 */
void print_vector_in_file (const PETScWrappers::MPI::BlockVector &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a std::vector<double> in a file with an introduction.
 */
void print_vector_in_file (const std::vector<double> &vect,
  std::string filename,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a std::vector<double> in a file with an introduction.
 */
void print_vector_in_file (const std::vector<unsigned int> &vect,
  std::string filename,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a std::vector<double> in a file with an introduction.
 */
void print_vector_in_file (const PETScWrappers::MPI::Vector &vect,
  std::string filename,
  std::string introduction,
  bool horizontal = false,
  int precision = 6);

/**
 * @brief Print a std::vector<Point<dim> > in a File with an introduction.
 */
template <int dim>
  void print_vector_in_file (const std::vector<Point<dim> > &vect,
    std::string filename,
    std::string introduction,
    bool horizontal = false,
    int precision = 6);

/**
 *
 */
void print_matrix_in_file (const std::vector<std::vector<double> > &mat,
  std::string filename,
  std::string introduction = "\n \n",
  int precision = 6);

/**als) :
 comm(_comm), T(_T), mate
 *
 */
void print_matrix_in_file (const std::vector<std::vector<double> > &mat,
  std::ofstream& out,
  std::string introduction = "\n \n",
  int precision = 6);

/**
 * @brief Print a std::vector<double> in a File with an introduction.
 */
template <typename num>
  void print_in_file (num numb,
    std::string filename,
    std::string introduction = "\n",
    unsigned int precision = 6);

/**
 * @brief The same as before but in a ofstream
 */
template <typename num>
  void print_in_file (num numb,
    std::ofstream &out_stream,
    std::string introduction = "\n",
    unsigned int precision = 6);

/**
 * @brief Print a list of points in a file
 */
template <int dim>
  void print_points_in_file (const std::vector<Point<dim> > &vect,
    std::string filename,
    std::string introduction = "\n \n",
    bool horizontal = false,
    int precision = 6);

#endif /* PRINTING_H_ */
