/**
 * @file   utils.h
 * @brief
 */
#ifndef UTILS_H
#define UTILS_H

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/table.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/fe/fe_q.h>

#include <petscsys.h>
#include <slepceps.h>
#include <petscksp.h>
#include <petscerror.h>

#include <cstddef>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

using namespace dealii;
typedef std::complex<double> complex;

DeclException1(ExcCantConvertString, std::string,
  << "Can't convert the string " << arg1 << " to the desired type");

//
// bool petscBool_to_bool(PetscBool flag)
//   Convert a PetscBool into a C bool.
bool petscBool_to_bool (PetscBool flag);

//
// is_similar(double num1, double num2, double tol)
//
bool is_similar (const double num1,
  const double num2,
  const double tol = 1e-5);

//
/**
 * Assert the condition also in  Release Mode. Useful to check input parameters.
 */
void AssertRelease (bool condition,
  const std::string &message);

/**
 *  Assert if two vectors are similarly equal.
 */
template <typename NUMBER>
  void assert_vectors_similar (std::vector<NUMBER> test_vector,
    std::vector<NUMBER> ref_vector,
    double tol = 1e-5);

/**
 *    Get a boolean from the options commands with the given name, the default answer is
 *     not to change the given flag.
 */
PetscErrorCode get_bool_from_options (const std::string &name,
  bool &flag);

/**
 * Get an unsigned int from the options commands with the given name, the default answer
 *  is not to change the given unsigned int.
 */
PetscErrorCode get_uint_from_options (const std::string &name,
  unsigned int &result);

/**
 *  Get a double from the options commands with the given name, the default answer is not to change the
 *  given number.
 */
PetscErrorCode get_double_from_options (const std::string &name,
  double &result);

/**
 *  Get a std::string from the options commands with the given name, the default answer is not
 *  to change the given string.
 */
PetscErrorCode get_string_from_options (const std::string &name,
  std::string &result);

/**
 * Check if the file exist and it can be open.
 */
bool fexists (std::string filename);

/*
 * Trim string from start (in place)
 */
static inline void ltrim (std::string &s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [] (int ch)
    {
      return !std::isspace(ch);
    }));
}

/**
 * Trim string from end (in place)
 */
static inline void rtrim (std::string &s)
{
  s.erase(std::find_if(s.rbegin(), s.rend(), [] (int ch)
    {
      return !std::isspace(ch);
    }).base(), s.end());
}

/*
 * Trim string from both ends (in place)
 */
static inline void trim (std::string &s)
{
  ltrim(s);
  rtrim(s);
}

/**
 * Trim string from start (copying)
 */
static inline std::string ltrim_copy (std::string s)
{
  ltrim(s);
  return s;
}

/**
 * Trim string from end (copying)
 */
static inline std::string rtrim_copy (std::string s)
{
  rtrim(s);
  return s;
}

/**
 * @brief Trim from both ends (copying)
 */
static inline std::string trim_copy (std::string s)
{
  trim(s);
  return s;
}

/**
 * Return true if the string starts with a commentary symbol.
 */
bool is_commentary (std::string &str);

/**
 * @bref Round number to a number of digits
 */
double round (double x,
  const unsigned int digits);

/**
 *
 */
template <typename T>
  std::string num_to_str (T Number);

/**
 *
 */
template <typename T>
  std::string vec_to_str (std::vector<T> vector);

/**
 *
 */
template <typename T>
  T str_to_num (std::string &Text);

/**
 * Convert a string to a vector in a comma separated format.
 */
template <typename Number>
  void
  str_to_vector (std::string &in,
    std::vector<Number> &out);

/**
 * Convert a string to a vector in a comma separated format.
 */
template <typename Number>
  void
  str_to_vector (std::string &in,
    std::vector<std::vector<Number> > &out);

/**
 * Convert a string to a vector in a comma separated format.
 */
template <typename Number>
  void
  vector_to_str (const std::vector<Number> &in,
    std::string &out,
    const std::string &indent = std::string(""),
    const unsigned int n_indent = 0);

/**
 * Convert a string to a vector in a comma separated format.
 * We only use the indent for the pretty output
 */
template <typename Number>
  void
  vector_to_str (const std::vector<std::vector<Number> > &in,
    std::string &out,
    const std::string &indent = std::string(""),
    const unsigned int n_indent = 0);

/** */
void str_to_vector_complex (std::string &in,
  std::vector<complex> &out);

/** */
void str_to_vector_complex (std::string &in,
  std::vector<std::vector<complex> > &out);

template <int dim>
  RefinementCase<dim> getRefCase (char *str);

/**
 *
 */
std::string get_new_valid_line (std::ifstream &input,
  std::string &line);

/**
 * Get a new valid line_stream from the input stream.
 */
std::ifstream get_new_valid_line_stream (std::ifstream&,
  std::string&);

/**
 * Return the value_points of a a solution along a line,
 *  it needs the triangulation for to check and future uses of this function.
 */
template <int dim>
  void make_line (DoFHandler<dim> &dof_handler,
    Vector<double> &solution,
    Point<dim> from,
    Point<dim> to,
    unsigned int steps,
    std::vector<double> &out,
    std::vector<double> &out_x);

/**
 *  Sum a vector of numbers and returns the result.
 */
template <class number>
  number sum_vector (const std::vector<number> &vector);

/**
 * Return the average of a std::vector<double>.
 */
double average (const std::vector<double> &vectortoaverage);

/**
 * Return the average of a std::vector<double> making
 * the absolute value of each value.
 */
double average_abs (const std::vector<double> &vectortoaverage);

/**
 * Returns the l2 norm of a std::vector<doubl i.e the square root of the elements squared.
 */
double l2_norm (const std::vector<double> vector);

/**
 * Returns the l2 norm of a std::vector<doubl i.e the square root of the elements squared.
 */
double l2_norm (const std::vector<std::vector<double> > vector);

/**
 * Normalize a vector with a given factor.
 */
template <class VECTOR>
  void normalize_vector (VECTOR &vect,
    const double factor);

/**
 * Normalize a vector with a given vector of factors.
 */
template <class number>
  void normalize_vector (
    std::vector<number> &vect,
    const std::vector<double> &factor);

/**
 * Normalize a vector with a given vector of factors.
 */
void normalize_vector (
  std::vector<std::vector<double> > &vect,
  std::vector<std::vector<double> > &factor);

/**
 * Normalize a vector of vectors by dividing every element by a given factor.
 */
void normalize_vector (
  std::vector<std::vector<double> > &vect,
  const double factor);

/**
 * Normalize a vector of vectors of vectors by dividing every element by a given factor.
 */
void normalize_vector (
  std::vector<std::vector<std::vector<double> > > &vect,
  const double factor);

/**
 * Normalize a vector of vector with a given vector of factors.
 */
void normalize_vector (
  std::vector<std::vector<double> > &vect,
  const std::vector<double> &factor);

/**
 * Normalize a vector of vector with a given vector of factors.
 */
void normalize_vector (std::vector<std::vector<std::vector<double> > > &vect,
  const std::vector<double> &factor);

/**
 * Normalize a vector with a given factor
 */
void normalize_vector (Table<2, double> &vect,
  const double factor);

/**
 * Normalize a vector with a given vector of factors
 */
void normalize_vector (Table<2, double> &vect,
  const Table<2, double> &factor);

/**
 * @brief Normalize a vector with a given vector of factors
 */
void normalize_vector (Table<1, double> &vect,
  const double factor);

/**
 * @brief Normalize a vector with a given vector of factors
 */
void normalize_vector (Table<1, double> &vect,
  const Table<1, double> &factor);

/**
 * @brief Assignment sol = a*u
 */
template <typename Number>
  void equ (const Number a,
    const std::vector<Number> &u,
    std::vector<Number> &sol);

/**
 * @brief Assignment sol = a*u + b*v.
 */
template <typename Number>
  void equ (const Number a,
    const std::vector<Number> &u,
    const Number b,
    const std::vector<Number> &v,
    std::vector<Number> &sol);

/**
 * @brief Assignment sol = a*u + b*v + b*w.
 */
template <typename Number>
  void equ (const Number a,
    const std::vector<Number> &u,
    const Number b,
    const std::vector<Number> &v,
    const Number c,
    const std::vector<Number> &w,
    std::vector<Number> &sol);

/**
 *
 */
template <class number>
  void parse_multiline_vector (std::ifstream &input_file,
    unsigned int nlines,
    std::vector<number> &out,
    bool deprecateFirst = true);

/**
 *
 */
template <class number>
  void parse_multiline_vector (std::ifstream &input_file,
    unsigned int nlines,
    std::vector<std::vector<number> > &out,
    bool deprecateFirst = true);

/**
 * Convert string to lower case
 */
std::string lower_case (std::string &data);

/**
 * Convert const char* to lower case
 */
std::string lower_case (const char *data);

/**
 *
 */
void getVec (const std::string filename,
  const std::string header,
  std::vector<double> &out,
  const unsigned int length);

/**
 *
 */
double getDouble (const std::string &filename,
  const std::string &header,
  double def = 0.0);

/**
 *
 */
template <int dim>
  std::vector<unsigned int> vectorize (const Point<dim> p1,
    const Point<dim> p2);

/**
 *
 */
template <int dim>
  void addRectangle (Triangulation<dim> &tri_in,
    Point<dim> p1,
    Point<dim> p2,
    double Lx,
    double Ly,
    double Lz);

/**
 * @brief Print a triangulation to a .eps file.
 */
template <int dim>
  void print_grid (Triangulation<dim> &tri,
    std::string filename);

/**
 * @brief  Prints a std::vector to the console.
 * Space separated values. One line.
 */
template <class num>
  void print_vector (std::vector<num> vect,
    bool end_line = true);

/**
 * @brief Prints a std::vector to the console.
 *  Space separated values. One line.
 */
template <class VECTOR>
  void print_vector (VECTOR vect,
    bool end_line = true);

/**
 * @brief Prints a std::vector to the console.
 * Space separated values. One line.
 */
void print_vector (PETScWrappers::MPI::Vector vect,
  bool end_line = true);

/**
 * @brieff Calculate the percentile of an unsorted vector v. In other words,
 * value that is placed in the percent place of a sorted vector.
 */
float percentile (Vector<float> v,
  const float percent);

/**
 * @brief Some utilities to change the indexes when it is used Gauss4 quadrature.
 *  It computes the row and the colum given...
 */
unsigned int findRow (const std::vector<unsigned int> nCol,
  const int index,
  int &col,
  int &sum);

/**
 * @brief Make the geometry points from a specific triangulation, and its vector of cell
 * centers.
 */
template <int dim>
  void makeGeometyPoints (
    const std::vector<std::pair<Point<dim>, unsigned int> > &centers,
    const Triangulation<dim> &tria,
    std::vector<unsigned int> &geometry_points,
    std::vector<unsigned int> &nRods_out);

/**
 *
 */
template <int dim, class num>
  void print_table (const Table<dim, num> tab);

/**
 *  Parse a vector in a file after a headline.
 */
void parse_vector_in_file (const std::string &file,
  const std::string &headline,
  std::vector<double> &vector_out,
  const unsigned int n_lines = 1,
  const unsigned int expected_vector_size = -1);

/**
 *  Parse a vector in a file after a headline.
 */
void parse_vector_in_file (const std::string &file,
  const std::string &headline,
  std::vector<unsigned int> &vector_out,
  const unsigned int n_lines = 1,
  const unsigned int expected_vector_size = -1);

/**
 * @brief This function construct the default rectangular geometry_points.
 * E.G 2D, 9 by 4 reactor:
 *  1 9 1 9 1 9 1 9 1 9
 */
std::vector<unsigned int> default_geometry_points (
  std::vector<unsigned int> n_cells_per_dim);

/*
 * Parse a matrix from a space and '\n' separated string.
 * n_rows and n_cols can be set to ensure a matrix size
 */
template <typename num>
  void parse_matrix (const std::string &input,
    std::vector<std::vector<num> > &out_matrix,
    const unsigned int n_rows = 0,
    const unsigned int n_cols = 0);

/**
 * @brief Parse a vector from an string. Check at the end if the given vector have
 * the expected length. If the string is empty return the default input.
 * @input The input string.
 * @input The output vector.
 * @input The expected length. If 0, do not check the length.
 * @input The default vector in case the string is empty.
 */
void parse_vector (std::string input,
  std::vector<unsigned int> &out,
  unsigned int length = 0,
  std::vector<unsigned int> def = std::vector<unsigned int>());

void parse_vector (std::string input,
  std::vector<int> &out,
  unsigned int length = 0,
  std::vector<int> def = std::vector<int>());

/**
 * @brief Same as before but for double vector. There is not default values in this case.
 * It exist an special character in order to repeat number:
 * 4* 10.0 = 10.0 10.0 10.0 10.0
 * 2.0 2*1.0 2.0 = 2.0 1.0 1.0 2.0
 *
 * @input The input string.
 * @input The output vector.
 * @input The expected length. If 0, do not check the length.
 */
void parse_vector (
  std::string input,
  std::vector<double> &out,
  unsigned int length = 0);

/**
 *
 */
void parse_vector (std::string input,
  PETScWrappers::MPI::BlockVector &out,
  unsigned int n_blocks,
  unsigned int n_dofs_per_block);

/**
 *
 */
void parse_vector (std::string input,
  PETScWrappers::MPI::Vector &out,
  unsigned int n_dofs);

/**
 *
 */

void parse_vector (std::string input,
  std::vector<std::string> &out,
  unsigned int length = 0);

/**
 *
 */
void parse_materials (std::string XECSFile,
  std::vector<unsigned int> &n_cells_per_dim,
  unsigned int n_assemblies,
  std::vector<unsigned int> &materials);

/**
 *
 */
template <int dim>
  void extrude_triangulation (const Triangulation<2, 2> &input,
    const unsigned int n_slices,
    const double height,
    Triangulation<dim> &result,
    const bool changeMat);

/**
 * This function does a binary search and returns the index i such that
 * x is between wl[i] and wl[i+1], except that i is restricted to the
 * range from 0 to n-2 inclusive.
 */
unsigned int binary_search (double x,
  std::vector<double> wl,
  unsigned int n);

/**
 * interp_linear(double x, std::vector<double> wl, std::vector<double> f,
 * unsigned int n, unsigned int *starti)
 *
 *  This routine uses linear interpolation.
 * Extracted from
 * https://trac.stsci.edu/ssb/stsci_python/browser/hstcal/trunk/pkg/stis/calstis/lib/interp1d.c
 *   The function value is the interpolated value.
 *
 *   The independent variable array (wl) is assumed to be monotonically
 *   increasing.  If x is outside the range of wl values, the nearest
 *   element (i.e. the first or last) in the dependent variable array (f)
 *   will be returned.
 *
 *   Before the first call to this routine, starti can be initialized to
 *   some value such as one.  This is the starting index for searching for
 *   the minimum difference between x and elements of wl; starti will be
 *   updated by this function.
 *
 * double x               : the value at which the function is to be evaluated
 * std::vector<double> wl : the vector of independent variable values
 * double f[]             : the array of dependent variable values
 * unsigned int n         : size of wl and f arrays
 * unsigned int *starti   : begin here for finding nearest wl to x
 */
double interp_linear (double x,
  std::vector<double> wl,
  std::vector<double> f,
  unsigned int n,
  unsigned int &starti);

/*
 * @brief Compute maximum and minimum vertex of cell in the selected coordinate.
 */
template <int dim>
  void getMaxMinVertex (
    TriaIterator<DoFCellAccessor<dim, dim, false> > cell,
    unsigned int coord,
    double &maxp,
    double &minp);

/**
 *
 */
template <class DH, class SparsityPattern>
  void make_partial_sparsity_pattern (const DH &dof_handler,
    SparsityPattern &sparsity,
    const std::vector<types::global_dof_index> &map_global_to_partial,
    unsigned int i_start,
    unsigned int i_end);

/**
 *
 */
template <class DH, class SparsityPattern>
  void make_partial_sparsity_pattern (const DH &dof_handler,
    SparsityPattern &sparsity,
    const std::vector<types::global_dof_index> &map_global_to_partial_i,
    unsigned int i_start,
    unsigned int i_end,
    const std::vector<types::global_dof_index> &map_global_to_partial_j,
    unsigned int j_start,
    unsigned int j_end);

/**
 *
 */
template <int dim>
  void map_global_to_partial (const DoFHandler<dim> &dof_handler,
    std::vector<types::global_dof_index> &map_global_to_partial,
    const unsigned int i_start,
    const unsigned int i_end);

/**
 *
 */
template <int dim>
  void map_partial_to_global (const DoFHandler<dim> &dof_handler,
    std::vector<types::global_dof_index> &map_partial_to_global,
    const unsigned int n_partial_dofs,
    const unsigned int i_start,
    const unsigned int i_end);

/**
 * Translate a subpartial_to_global map  to a subpartial_to_partial map using the corresponding
 * map_global_to_partial dictionary.
 */
void map_subpartial_to_partial (
  const std::vector<types::global_dof_index> &map_subpartial_to_global,
  const std::vector<types::global_dof_index> &map_global_to_partial,
  std::vector<types::global_dof_index> &map_subpartial_to_partial);

/**
 *
 */
void map_subpartials_to_partial (
  const std::vector<std::vector<types::global_dof_index> > &map_subpartials_to_global,
  const std::vector<types::global_dof_index> &map_global_to_partial,
  std::vector<std::vector<types::global_dof_index> > &map_subpartials_to_partial,
  std::vector<unsigned int> &partition_of_unity_global,
  std::vector<std::vector<double> > &partition_of_unity);

/**
 *
 */
template <int dim>
  void map_edges_to_global (const DoFHandler<dim> &dof_handler,
    std::vector<std::vector<types::global_dof_index> > &partials_global_map,
    std::vector<unsigned int> &partition_of_unity_vector,
    unsigned int &n_subdomains,
    unsigned int &n_dofs_per_subdomain);

/**
 *
 */
template <class VECTOR>
  void print_vector_in_matlab (const VECTOR &vec,
    const std::string &vec_name = "vec",
    std::ostream &out = std::cout,
    const unsigned int precision = 8);

/**
 *
 */
template <class BLOCKVECTOR>
  void print_block_vector_in_matlab (const BLOCKVECTOR &vec,
    const std::string &vec_name = "vec",
    std::ostream &out = std::cout,
    const unsigned int precision = 8);

/**
 *
 */
template <class MATRIX>
  void print_matrix_in_matlab (const MATRIX &mat,
    const std::string &mat_name = "mat",
    std::ostream &out = std::cout,
    const unsigned int precision = 8);

/**
 *
 */
template <class VECTOR>
  void print_vector_in_python (const VECTOR &vec,
    const std::string &vec_name = "vec",
    std::ostream &out = std::cout,
    const unsigned int precision = 8);

/**
 *
 */
template <class MATRIX>
  void print_matrix_in_python (const MATRIX &mat,
    const std::string &mat_name = "mat",
    std::ostream &out = std::cout,
    const unsigned int precision = 8);

/**
 * Copy a Petsc Vec to deal.ii Vector<double>.
 *    Be careful this function involves to copy a vector entry by entry.
 */
void copy_to_Vector (
  Vector<double> &dst,
  Vec src);

/**
 * Copy a Petsc Vec to deal.ii Vector<double>. Be careful this function involves copy a vector.
 */
void copy_to_BlockVector (
  BlockVector<double> &dst,
  Vec src);

/**
 * @brief Copy a LinearAlgebra::distributed::Vector to PETScWrappers::MPI::Vector.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vector (PETScWrappers::MPI::Vector &dst,
  const LinearAlgebra::distributed::Vector<double, MemorySpace::Host> &src);

/**
 * @brief Copy a Petsc Vec to LinearAlgebra::distributed::Vector.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vector (LinearAlgebra::distributed::Vector<double, MemorySpace::Host> &dst,
  Vec src);

/**
 * @brief Copy a Petsc Vec to std::vector<PETScWrappers::MPI::BlockVector>
 * Be careful this function involves copy a vector.
 */
void copy_to_stdBlockVector (std::vector<PETScWrappers::MPI::BlockVector> &dst,
  Vec src);

/**
 * @brief Copy a Petsc Vec to ETScWrappers::MPI::BlockVector.
 * Be careful this function involves copy a vector.
 */
void copy_to_BlockVector (PETScWrappers::MPI::BlockVector &dst,
  Vec src);

/**
 * @brief Copy a  BlockVector<double> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const BlockVector<double> &src);

/**
 * @brief Copy a  std::vector<PETScWrappers::MPI::BlockVector> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const std::vector<PETScWrappers::MPI::BlockVector> &src);

/**
 * @brief Copy a  PETScWrappers::MPI::BlockVector to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const PETScWrappers::MPI::BlockVector &src);

/**
 * @brief Copy a   LinearAlgebra::distributed::Vector to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const LinearAlgebra::distributed::Vector<double, MemorySpace::Host> &src);

/**
 * @brief Copy a  Vector<double> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const Vector<double> &src);

/**
 *
 */
void compute_max_eig (
  const LAPACKFullMatrix<double> &matrix,
  double &eig,
  std::vector<double> &eigenvector);

/**
 *
 */
void compute_eigs (
  const LAPACKFullMatrix<double> &matrix,
  std::vector<double> &wr,
  LAPACKFullMatrix<double> &eigenvectors,
  bool sort = false);

/**
 * @brief Join and order two vector into one.
 */
template <int dim>
  void join_vectors (
    const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const PETScWrappers::MPI::Vector &phi_1,
    const PETScWrappers::MPI::Vector &phi_2,
    PETScWrappers::MPI::Vector &phi_joint);

/**
 * @brief Join and order two vector into one.
 */
template <int dim>
  void join_vectors (const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const std::vector<PETScWrappers::MPI::Vector> &phi_sep,
    PETScWrappers::MPI::Vector &phi_joint);

/**
 * @brief Separate a vector into two ones.
 */
template <int dim>
  void separate_vectors (const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const PETScWrappers::MPI::Vector &phi_sol,
    PETScWrappers::MPI::Vector &phi_1,
    PETScWrappers::MPI::Vector &phi_2);

/**
 *
 */
template <int dim>
  void separate_vectors (const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const PETScWrappers::MPI::Vector &phi_joint,
    std::vector<PETScWrappers::MPI::Vector> &phi_sep);

/**
 * @brief Count the real number of non_zero_elements
 */
template <class MATRIX>
  unsigned int n_real_nonzero_elements (const MATRIX &mat,
    double tol = 1e-12);

/**
 * @brief Return the flag that corresponds to cutting a cell along the axis given as argument
 */
template <int dim>
  RefinementCase<dim>
  cut_axis (const unsigned int);

template <>
  RefinementCase<1>
  cut_axis<1> (const unsigned int);

template <>
  RefinementCase<2>
  cut_axis<2> (const unsigned int);

template <>
  RefinementCase<3>
  cut_axis<3> (const unsigned int);

/**
 * @brief Interpolates the discrete function in, which is a vector on the grid before the refinement, to the function out which then is a vector on the refined grid
 */
template <int dim>
  void interpolate_operator (
    const DoFHandler<dim> &dof_handler,
    const std::vector<std::vector<unsigned int> > &cell_map,
    const PETScWrappers::MPI::Vector &vec_in,
    PETScWrappers::MPI::Vector &vec_out);

template <int dim>
  void restriction_operator (
    const DoFHandler<dim> &dof_handler,
    const std::vector<std::vector<unsigned int> > &cell_map,
    const PETScWrappers::MPI::Vector &vec_in,
    PETScWrappers::MPI::Vector &vec_out);

/**
 * @brief This function transfers the solution phi_old of fe coarse
 * with n_fe_degree=p_coarse into the dof_handler
 */
template <int dim>
  void solution_transfer_dof (
    const Triangulation<dim> &tria,
    unsigned int p_coarse,
    const DoFHandler<dim> &dof_handler,
    const FE_Q<dim> &fe,
    std::vector<PETScWrappers::MPI::BlockVector> &phi_old,
    std::vector<PETScWrappers::MPI::BlockVector> &phi_new);

/**
 * @brief This function transfers the solution phi_old of fe coarse
 * with n_fe_degree=p_coarse into the dof_handler
 */
template <int dim>
  void solution_transfer_dof (
    const Triangulation<dim> &tria,
    unsigned int p_coarse,
    const DoFHandler<dim> &dof_handler,
    const FE_Q<dim> &fe,
    std::vector<BlockVector<double> > &phi_old,
    std::vector<BlockVector<double> > &phi_new);

/**
 * @brief This function create the transfer matrix
 * to use the solution_interpolate_dof
 */
template <int dim>
  void
  set_transfer_matrix (
    const FE_Q<dim> &fe_in,
    const FE_Q<dim> &fe_out,
    FullMatrix<double> &transfer_mat);

/**
 * @brief This function create the transfer matrix
 * to use the solution_interpolate_dof
 */
template <int dim>
  void set_transfer_matrix (
    const FiniteElement<dim> &fe_in,
    const FiniteElement<dim> &fe_out,
    FullMatrix<double> &transfer_mat);

/**
 * @brief This function transfers the solution in
 * into the solution out. It is needed the dofhandler_in,
 * the dof_handler_out and the transfer mat
 */
template <int dim>
  void solution_interpolate_dof (
    const DoFHandler<dim> &dof_handler_in,
    const DoFHandler<dim> &dof_handler_out,
    const FullMatrix<double> &transfer_mat,
    PETScWrappers::MPI::Vector &sol_in,
    PETScWrappers::MPI::Vector &sol_out);

/**
 * @brief This function transfers the solution in
 * into the solution out. It is needed the dofhandler_in,
 * the dof_handler_out and the transfer mat
 */
template <int dim>
  void solution_interpolate_dof (
    const DoFHandler<dim> &dof_handler_in,
    const DoFHandler<dim> &dof_handler_out,
    const FullMatrix<double> &transfer_mat,
    PETScWrappers::MPI::BlockVector &sol_in,
    PETScWrappers::MPI::BlockVector &sol_out);

/**
 *
 */
const std::vector<unsigned int> child_pos (
  const unsigned int child,
  const unsigned int n_levels,
  unsigned int dim);

/**
 *
 */
bool child_at_face (
  const unsigned int face,
  const unsigned int child,
  const unsigned int n_children,
  const unsigned int dim);

#endif

