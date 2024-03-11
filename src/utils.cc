/**
 * @file   utils.cc
 * @brief
 */
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/table.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/index_set.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/fe/fe_tools.h>

#include <petscsys.h>
#include <slepceps.h>
#include <petscksp.h>
#include <petscvec.h>
#include <petscoptions.h>

#include <cstddef>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

#include <boost/algorithm/string.hpp>

#include "../include/utils.h"

using namespace dealii;

//
// bool petscBool_to_bool(PetscBool flag)
//   Convert a PetscBool into a C bool.
bool petscBool_to_bool (PetscBool flag)
{
  if (flag == PETSC_TRUE)
    return true;
  else
    return false;
}

//
// is_similar(double num1, double num2, double tol)
//
bool is_similar (const double num1,
  const double num2,
  const double tol)
{
  return (std::abs(num1 - num2) < tol);
}

//
// AssertRealse(bool condition, std::string message);
// Assert the condition also in  Release Mode. Useful to check input parameters.
//
void AssertRelease (bool condition,
  const std::string &message)
{
  if (!condition)
  {
    std::cerr << std::endl;
    std::cerr << "----------------------------------------------------"
              << std::endl;
    std::cerr << "ERROR! " << std::endl;
    std::cerr << "Exception on release processing: " << std::endl;
    std::cerr << "   " << message << std::endl << "Aborting!" << std::endl;
    std::cerr << "----------------------------------------------------"
              << std::endl;
    exit(1);
  }
}

/**
 *  Assert if two vectors are similarly equal up to a tolerance.
 */
template <typename number>
  void assert_vectors_similar (std::vector<number> test_vector,
    std::vector<number> ref_vector,
    double tol)
  {

    AssertRelease(test_vector.size() == ref_vector.size(),
      "Vectors do not have the same size");
    for (unsigned int i = 0; i < ref_vector.size(); ++i)
      AssertRelease(is_similar(test_vector[i], ref_vector[i], tol),
        "Error in " + num_to_str(i) + " entry. "
        + num_to_str(test_vector[i])
        + " vs "
        + num_to_str(ref_vector[i]));
  }

template void assert_vectors_similar (
  std::vector<double> test_vector,
  std::vector<double> ref_vector,
  double tol);
template void assert_vectors_similar (
  std::vector<unsigned int> test_vector,
  std::vector<unsigned int> ref_vector,
  double tol);

/**
 *    Get a boolean from the options commands with the given name, the default answer is
 *     not to change the given flag.
 */
PetscErrorCode get_bool_from_options (const std::string &name,
  bool &flag)
{
  PetscBool flg;
  PetscBool set;
  PetscErrorCode ierr;
#if PETSC_VERSION_MINOR > 6
  ierr = PetscOptionsGetBool(NULL, NULL, name.c_str(), &flg, &set);
#else
  ierr = PetscOptionsGetBool(NULL, name.c_str(), &flg, &set);
#endif
  CHKERRQ(ierr);
  if (set == PETSC_TRUE)
    flag = petscBool_to_bool(flg);

  return ierr;
}

/**
 * @brief Get an unsigned int from the options commands with the given name, the default answer is not to
 * change the given string.
 */
PetscErrorCode get_uint_from_options (const std::string &name,
  unsigned int &result)
{
  PetscErrorCode ierr;
  int i = (int) result;
#if PETSC_VERSION_MINOR > 6
  ierr = PetscOptionsGetInt(NULL, NULL, name.c_str(), &i, NULL);
#else
  ierr = PetscOptionsGetInt(NULL, name.c_str(), &i, NULL);
#endif
  CHKERRQ(ierr);
  result = (unsigned int) i;
  return ierr;
}

/**
 * @brief Get an unsigned int from the options commands with the given name, the default answer is not to
 * change the given string.
 */
PetscErrorCode get_double_from_options (const std::string &name,
  double &result)
{
  PetscErrorCode ierr;

#if PETSC_VERSION_MINOR > 6
  ierr = PetscOptionsGetReal(NULL, NULL, name.c_str(), &result, NULL);
#else
  ierr = PetscOptionsGetReal(NULL, name.c_str(), &result, NULL);
#endif
  CHKERRQ(ierr);
  return ierr;
}

/**
 * @brief Get a std::string from the options commands with the given keyword.
 *    If the keyword is not given the default behavior is not to change the
 *    given result string.
 */
PetscErrorCode get_string_from_options (const std::string &keyword,
  std::string &result)
{
  char result_char[PETSC_MAX_PATH_LEN];
  PetscBool flg;
  PetscErrorCode ierr;

#if PETSC_VERSION_MINOR > 6
  ierr = PetscOptionsGetString(NULL, NULL, keyword.c_str(), result_char,
    sizeof(result_char), &flg);
#else
  ierr = PetscOptionsGetString(NULL, keyword.c_str(), result_char, sizeof(result_char), &flg);
#endif
  CHKERRQ(ierr);
  if (flg)
    result = result_char;

  return ierr;
}

//
// bool fexists(std::string filename)
//   Check if the file exist and it can be open.
bool fexists (std::string filename)
{
  std::ifstream ifile(filename.c_str());
  return ifile.is_open();
}

//
// bool is_commentary(std::string &str)
//   Return true if the string starts with a commentary symbol.
bool is_commentary (std::string &str)
{
  if (str == "") // Blank line
    return true;
  if (str.compare("#") == 0)
    return true; // Compare == 0 if equal
  if (str.compare("!") == 0)
    return true; // Compare == 0 if equal
  if (str.compare(0, 2, "//") == 0)
    return true; // Compare == 0 if equal
  else
    return false;
}

/**
 *
 */
double round (double x,
  const unsigned int digits)
{
  return round(x * pow(10, (double) digits)) / pow(10, (double) digits);
}
//
// std::string num_to_str(T Number)
//   Convert a number to a string.
template <typename T>
  std::string num_to_str (T Number)
  {
    std::ostringstream ss;
    ss << Number;
    return ss.str();
  }

template std::string num_to_str (double);
template std::string num_to_str (unsigned int);
template std::string num_to_str (int);
template std::string num_to_str (float);
template std::string num_to_str (unsigned char);

//
// std::string vec_to_str(T Number)
//
template <typename T>
  std::string vec_to_str (std::vector<T> vector)
  {
    std::ostringstream ss;
    for (unsigned int i = 0; i < vector.size(); i++)
      ss << vector[i] << " ";
    return ss.str();
  }

template std::string vec_to_str (std::vector<double>);
template std::string vec_to_str (std::vector<unsigned int>);
template std::string vec_to_str (std::vector<int>);
template std::string vec_to_str (std::vector<float>);

//
// str_to_num(std::string &Text)
//   Converts a string to a number
//     E.g    str_to_num<int>('23') -> 23
//     E.g    str_to_num<double>('23.2') -> 23.2
//     E.g    str_to_num<double>('23') -> 23.0
template <typename T>
  T str_to_num (std::string &text)
  {
    trim(text);
    std::istringstream iss(text);
    T result;
    iss >> result;
    AssertRelease(!iss.fail(), "ExcCantConvertString " + text);
    return result;
  }

template double str_to_num (std::string&);
template unsigned int str_to_num (std::string&);
template int str_to_num (std::string&);

template <typename Number>
  void str_to_vector (std::string &in,
    std::vector<Number> &out)
  {
    std::vector<std::string> strs;
    boost::trim(in);
    boost::split(strs, in, boost::is_any_of(", "), boost::token_compress_on);
    for (unsigned int i = 0; i < strs.size(); ++i)
    {
      out.push_back(str_to_num<Number>(strs[i]));
    }
  }

template void str_to_vector (std::string &in,
  std::vector<double> &out);
template void str_to_vector (std::string &in,
  std::vector<int> &out);
template void str_to_vector (std::string &in,
  std::vector<unsigned int> &out);
template void str_to_vector (std::string &in,
  std::vector<float> &out);

/** @todo document me */
template <typename Number>
  void str_to_vector (std::string &in,
    std::vector<std::vector<Number> > &out)
  {
    boost::trim_if(in, boost::is_any_of("\n "));
    boost::trim_right_if(in, boost::is_any_of(";"));
    std::vector<std::string> strs;
    boost::split(strs, in, boost::is_any_of(";"), boost::token_compress_on);
    out.resize(strs.size());
    for (unsigned int i = 0; i < strs.size(); ++i)
    {
      str_to_vector(strs[i], out[i]);
    }
  }

template void str_to_vector (std::string &in,
  std::vector<std::vector<double> > &out);
template void str_to_vector (std::string &in,
  std::vector<std::vector<int> > &out);
template void str_to_vector (std::string &in,
  std::vector<std::vector<unsigned int> > &out);
template void str_to_vector (std::string &in,
  std::vector<std::vector<float> > &out);

/**
 *
 */
void str_to_vector_complex (std::string &in,
  std::vector<complex> &out)
{
  std::vector<std::string> strs;
  boost::trim(in);
  boost::split(strs, in, boost::is_any_of(", "), boost::token_compress_on);

  for (unsigned int i = 0; i < strs.size() / 2; ++i)
  {
    out.push_back(
      complex(str_to_num<double>(strs[2 * i]), str_to_num<double>(strs[2 * i + 1])));
  }
}

/**
 *
 */
void str_to_vector_complex (std::string &in,
  std::vector<std::vector<complex> > &out)
{
  boost::trim_if(in, boost::is_any_of("\n "));
  boost::trim_right_if(in, boost::is_any_of(";"));
  std::vector<std::string> strs;
  boost::split(strs, in, boost::is_any_of(";"), boost::token_compress_on);
  out.resize(strs.size());
  for (unsigned int i = 0; i < strs.size(); ++i)
  {
    str_to_vector_complex(strs[i], out[i]);
  }
}

/** @todo document me */
template <typename Number>
  void vector_to_str (const std::vector<Number> &in,
    std::string &out,
    const std::string &indent,
    const unsigned int n_indent)
  {
    out.clear();

    // We should check that the vector is not empty
    assert(in.size() != 0);

    // We put the numbers inside the string
    std::stringstream ss;
    ss << std::fixed;
    ss.precision(5);
    // ss.setf( std::ios::fixed, std::ios::floatfield);

    for (unsigned int i = 0; i < in.size(); ++i)
    {
      ss << std::setw(7) << in[i] << " ";
    }

    //std::copy(in.begin(), in.end(), std::ostream_iterator<Number>(ss, " "));
    for (unsigned int j = 0; j < n_indent; ++j)
      out += indent;
    out += ss.str();
    out = out.substr(0, out.length() - 1);

  }

template void
vector_to_str (const std::vector<double> &in,
  std::string &out,
  const std::string &indent,
  const unsigned int n_indent);

template void
vector_to_str (const std::vector<unsigned int> &in,
  std::string &out,
  const std::string &indent,
  const unsigned int n_indent);

/**
 *  We only use the indent for the pretty output
 */
template <typename Number>
  void vector_to_str (const std::vector<std::vector<Number> > &in,
    std::string &out,
    const std::string &indent,
    const unsigned int n_indent)
  {
    out.clear();
    out += "\n";
    std::string bin;
    for (unsigned int i = 0; i < in.size(); ++i)
    {
      vector_to_str(in[i], bin, indent, n_indent);
      out += bin + ";\n";
    }
    for (unsigned int j = 0; j < n_indent - 1; ++j)
      out += indent;
  }

template void vector_to_str (
  const std::vector<std::vector<double> > &in,
  std::string &out,
  const std::string &indent,
  const unsigned int n_indent);
template void vector_to_str (
  const std::vector<std::vector<unsigned int> > &in,
  std::string &out,
  const std::string &indent,
  const unsigned int n_indent);

/**
 *
 */
template <int dim>
  RefinementCase<dim> getRefCase (char *str)
  {
    if (*str == 'x')
      return RefinementCase<dim>::cut_x;
    else if (*str == 'y')
      return RefinementCase<dim>::cut_y;
    else if (*str == 'z')
      return RefinementCase<dim>::cut_z;
    else if (*str == 'i')
      return RefinementCase<dim>::isometric_refinement;
    else
      AssertRelease(false, "Refinement type not recognized");
    return RefinementCase<dim>::isometric_refinement;

  }

/**
 * Get a new valid line from the input stream.
 */
std::string get_new_valid_line (std::ifstream &input,
  std::string &line)
{
  std::string str;
  getline(input, line);
  trim(line);

  std::istringstream iss(line);
  iss >> str;

  // If this line is a commentary try with the next line
  if (is_commentary(str))
  {
    return get_new_valid_line(input, line);
  }

  return line;
}

/**
 *  Calculate the value_points of a solution along a line,
 *   it needs the triangulation for to check and future uses of this function.
 */
template <int dim>
  void make_line (DoFHandler<dim> &dof_handler,
    Vector<double> &solution,
    Point<dim> from,
    Point<dim> to,
    unsigned int steps,
    std::vector<double> &out,
    std::vector<double> &out_x)
  {
// TODO This function can be optimized using fe_values.
    Point<dim> p;
    out.resize(steps);
    out_x.resize(steps);
    for (unsigned int i = 0; i < steps; i++)
    {
      p = from + ((double) i / (double) (steps - 1)) * (to - from);
      out_x[i] = p(0);
      out[i] = VectorTools::point_value(dof_handler, solution, p);
    }
  }

template void make_line (DoFHandler<1> &dof_handler,
  Vector<double> &solution,
  Point<1> from,
  Point<1> to,
  unsigned int steps,
  std::vector<double> &out,
  std::vector<double> &out_x);

/**
 * @brief Sum a vector of numbers and throws the result.
 */
template <class number>
  number sum_vector (const std::vector<number> &vector)
  {
    number sum = 0;
    for (unsigned int i = 0; i < vector.size(); i++)
      sum += vector[i];
    return sum;
  }

template float sum_vector (const std::vector<float> &vector);
template double sum_vector (const std::vector<double> &vector);
template int sum_vector (const std::vector<int> &vector);
template unsigned int sum_vector (const std::vector<unsigned int> &vector);

/**
 * @brief Return the average of a std::vector<double>.
 */
double average (const std::vector<double> &vector)
{
  double solution = 0.0;
  for (unsigned int i = 0; i < vector.size(); i++)
    solution += vector[i];

  return (solution / vector.size());
}

/**
 *  Return the average of a std::vector<double> making
 *  the absolute value of each value.
 */
double average_abs (const std::vector<double> vectortoaverage)
{
  double solution = 0.0;
  for (unsigned int i = 0; i < vectortoaverage.size(); i++)
    solution += numbers::NumberTraits<double>::abs(vectortoaverage[i]);

  return (solution / vectortoaverage.size());
}

/**
 * Returns the l2 norm of a std::vector<doubl i.e the square root of the elements squared.
 */
double l2_norm (const std::vector<double> vector)
{
  double solution = 0.0;
  for (unsigned int i = 0; i < vector.size(); i++)
    solution += vector[i] * vector[i];

  return (sqrt(solution));
}

/**
 * @brief Returns the l2 norm of a std::vector<doubl i.e the square root of the elements squared.
 *
 */
double l2_norm (const std::vector<std::vector<double> > vector)
{
  double solution = 0.0;
  for (unsigned int i = 0; i < vector.size(); i++)
    for (unsigned int j = 0; j < vector[i].size(); j++)
      solution += vector[i][j] * vector[i][j];

  return (sqrt(solution));
}

/**
 * Normalize a vector by dividing every element by a given factor.
 */
template <class VECTOR>
  void normalize_vector (VECTOR &vect,
    const double factor)
  {
    const double factor_inverse = 1.0 / factor;
    for (unsigned int i = 0; i < vect.size(); i++)
    {
      vect[i] *= factor_inverse;
    }
  }

template void normalize_vector (std::vector<double> &vect,
  const double factor);
template void normalize_vector (std::vector<complex> &vect,
  const double factor);
template void normalize_vector (Vector<double> &vect,
  const double factor);
template void normalize_vector (PETScWrappers::MPI::Vector &vect,
  const double factor);

/**
 * Normalize a vector with a given vector of factors
 */
template <class number>
  void normalize_vector (std::vector<number> &vect,
    const std::vector<double> &factor)
  {
    assert(vect.size() == factor.size());
    for (unsigned int i = 0; i < vect.size(); i++)
      if (factor[i] != 0)
        vect[i] /= factor[i];
  }

template void normalize_vector (std::vector<double> &vect,
  const std::vector<double> &factor);
template void normalize_vector (std::vector<complex> &vect,
  const std::vector<double> &factor);

/**
 * Normalize a vector of vectors by dividing every element by a given factor.
 */
void normalize_vector (std::vector<std::vector<double> > &vect,
  const double factor)
{
  const double factor_inverse = 1.0 / factor;
  for (unsigned int i = 0; i < vect.size(); i++)
    for (unsigned int j = 0; j < vect[i].size(); j++)
      vect[i][j] *= factor_inverse;

}

//
// normalize_vector(std::vector<std::vector<std::vector<double> > > vect, double factor)
//   Normalize a vector of vectors of vectors by dividing every element by a given factor.
void normalize_vector (std::vector<std::vector<std::vector<double> > > &vect,
  const double factor)
{
  double factor_inverse = 1.0 / factor;
  for (unsigned int i = 0; i < vect.size(); i++)
    for (unsigned int j = 0; j < vect[i].size(); j++)
      for (unsigned int k = 0; k < vect[i][j].size(); k++)
        vect[i][j][k] *= factor_inverse;

}

//
// void normalize_vector(std::vector<std::vector<double> > &vect, std::vector<std::vector<double> > &factor)
//   Normalize a vector with a given vector of factors.
void normalize_vector (std::vector<std::vector<double> > &vect,
  const std::vector<std::vector<double> > &factor)
{
  assert(vect.size() == factor.size());

  for (unsigned int i = 0; i < vect.size(); i++)
  {
    assert(vect[i].size() == factor[i].size());
    for (unsigned int j = 0; j < vect[i].size(); j++)
      if (factor[i][j] != 0)
        vect[i][j] /= factor[i][j];
  }
}

//
// void normalize_vector(std::vector<std::vector<double> > &vect, std::vector<double> &factor)
//   Normalize a vector of vector with a given vector of factors.
void normalize_vector (std::vector<std::vector<double> > &vect,
  const std::vector<double> &factor)
{
  assert(vect.size() == factor.size());

  for (unsigned int i = 0; i < vect.size(); i++)
  {
    if (factor[i] != 0.0)
      for (unsigned int j = 0; j < vect[i].size(); j++)
        vect[i][j] /= factor[i];
  }
}

//
// void normalize_vector(std::vector<std::vector<double> > &vect, std::vector<double> &factor)
//   Normalize a vector of vector with a given vector of factors.
void normalize_vector (std::vector<std::vector<std::vector<double> > > &vect,
  const std::vector<double> &factor)
{
  assert(vect.size() == factor.size());
  double factor_inv;
  for (unsigned int i = 0; i < vect.size(); i++)
  {
    if (factor[i] != 0.0)
    {
      factor_inv = 1.0 / factor[i];
      for (unsigned int j = 0; j < vect[i].size(); j++)
        for (unsigned int k = 0; k < vect[i][k].size(); k++)
          vect[i][j][k] *= factor_inv;
    }
  }
}

/**
 * Normalize a vector with a given vector of factors
 */
void normalize_vector (PETScWrappers::MPI::Vector &vect,
  const double factor)
{
  const double factorInv = 1.0 / factor;
  for (unsigned int i = 0; i < vect.size(); i++)
    vect[i] *= factorInv;
}

/**
 * Normalize a vector with a given vector of factors
 */
void normalize_vector (Table<2, double> &vect,
  const double factor)
{
  const double factorInv = 1 / factor;
  for (unsigned int i = 0; i < vect.size(0); i++)
    for (unsigned int j = 0; j < vect.size(1); j++)
    {
      vect[i][j] *= factorInv;
    }
}

/**
 *
 */
void normalize_vector (Table<2, double> &vect,
  const Table<2, double> &factor)
{

  for (unsigned int i = 0; i < vect.size(0); i++)
    for (unsigned int j = 0; j < vect.size(1); j++)
    {
      if (factor[i][j] != 0.0 and vect[i][j] != 0.0)
        vect[i][j] /= factor[i][j];
    }
}

/**
 *  Normalize a vector with the given factor
 */
void normalize_vector (Table<1, double> &vect,
  const double factor)
{
  const double factorInv = 1 / factor;
  for (unsigned int i = 0; i < vect.size(0); i++)
  {
    vect[i] *= factorInv;
  }
}

/**
 *  Normalize a vector with a given vector of factors
 */
void normalize_vector (Table<1, double> &table,
  const Table<1, double> &factor)
{

  for (unsigned int i = 0; i < table.size(0); i++)
  {
    if (factor[i] != 0.0 and table[i] != 0.0)
      table[i] /= factor[i];
  }
}

/**
 * Assignment sol = a*u
 */
template <typename Number>
  void equ (const Number a,
    const std::vector<Number> &u,
    std::vector<Number> &sol)
  {
    if (sol.size() != u.size())
      sol.resize(u.size());

    for (unsigned int i = 0; i < sol.size(); i++)
      sol[i] = a * u[i];
  }

template void equ (const double a,
  const std::vector<double> &u,
  std::vector<double> &sol);

//
//void equ(const Number a, const std::vector<Number> &u, const Number b,
//  const std::vector<Number> &v, std::vector<Number> &sol)
//  Assignment sol = a*u + b*v
template <typename Number>
  void equ (const Number a,
    const std::vector<Number> &u,
    const Number b,
    const std::vector<Number> &v,
    std::vector<Number> &sol)
  {
    AssertDimension(u.size(), v.size());

    if (sol.size() != u.size())
      sol.resize(u.size());

    for (unsigned int i = 0; i < sol.size(); i++)
      sol[i] = a * u[i] + b * v[i];
  }

template void equ (const double a,
  const std::vector<double> &u,
  const double b,
  const std::vector<double> &v,
  std::vector<double> &sol);

/**
 * Assignment sol = a*u + b*v + b*w
 */
template <typename Number>
  void equ (const Number a,
    const std::vector<Number> &u,
    const Number b,
    const std::vector<Number> &v,
    const Number c,
    const std::vector<Number> &w,
    std::vector<Number> &sol)
  {
    AssertDimension(u.size(), v.size());
    AssertDimension(u.size(), w.size());

    if (sol.size() != u.size())
      sol.resize(u.size());

    for (unsigned int i = 0; i < sol.size(); i++)
      sol[i] = a * u[i] + b * v[i] + c * w[i];
  }

template void equ (const double a,
  const std::vector<double> &u,
  const double b,
  const std::vector<double> &v,
  const double c,
  const std::vector<double> &w,
  std::vector<double> &sol);

template <int dim>
  std::pair<unsigned int, RefinementCase<dim> > getRefPair (std::string str)
  {
    RefinementCase<dim> refCase;
    if (str[1] == 'x')
      refCase = RefinementCase<dim>::cut_x;
    else if (str[1] == 'y')
      refCase = RefinementCase<dim>::cut_y;
    else if (str[1] == 'z' and dim == 3)
      refCase = RefinementCase<dim>::cut_axis(2);
    else if (str[1] == 'z' and dim == 2)
      Assert(false, ExcMessage("Cannot refine in z in 2 dimensions!"))
    else if (str[1] == 'i')
      refCase = RefinementCase<dim>::isotropic_refinement;
    else
      Assert(false, ExcMessage("Refinement type not recognized"));
    return std::make_pair(std::atoi(&str[0]), refCase);
  }

template std::pair<unsigned int, RefinementCase<3> > getRefPair (
  std::string str);

template <>
  std::pair<unsigned int, RefinementCase<1> > getRefPair (std::string str)
  {
    RefinementCase<1> refCase;
    if (str[1] == 'x')
      refCase = RefinementCase<1>::cut_x;
    else
      Assert(false, ExcMessage("Refinement type not recognized"));
    return std::make_pair(std::atoi(&str[0]), refCase);
  }

template <>
  std::pair<unsigned int, RefinementCase<2> > getRefPair (std::string str)
  {
    RefinementCase<2> refCase;
    if (str[1] == 'x')
      refCase = RefinementCase<2>::cut_x;
    else if (str[1] == 'y')
      refCase = RefinementCase<2>::cut_y;
    else if (str[1] == 'z')
      Assert(false, ExcMessage("Cannot refine in z in 2 dimensions!"))
    else if (str[1] == 'i')
      refCase = RefinementCase<2>::isotropic_refinement;
    else
      Assert(false, ExcMessage("Refinement type not recognized"));
    return std::make_pair(std::atoi(&str[0]), refCase);
  }

template <int dim>
  void getRefinement (const std::string &filename,
    const std::string &header,
    std::vector<std::pair<unsigned int, RefinementCase<dim> > > &out,
    std::vector<unsigned int> &nCol,
    const int height)
  {
    std::ifstream myfile(filename.c_str(), std::ios::in);
    std::string line;
    unsigned int n_assemblies_per_plane = sum_vector(nCol);
    unsigned int i;

    std::string coeff, coeff2;
    if (myfile.is_open())
    {
      while (myfile.good())
      {
        myfile >> line;
        if (line == header)
        {
          for (int h = 0; h < height; h++)
          {
            myfile >> line; // Throw the word PLANE
            myfile >> line; // Throw the plane number
            for (i = 0; i < n_assemblies_per_plane; i++)
            {
              myfile >> coeff;

              if (*coeff.rbegin() == '*') // Last char is a *
              {
                coeff = coeff.substr(0, coeff.size() - 1);
                myfile >> coeff2;
                unsigned int len = Utilities::string_to_int(coeff);

                for (unsigned int j = 0; j < len; j++)
                  out.push_back(getRefPair<dim>(coeff2));

                i += len - 1;

              }
              else
                out.push_back(getRefPair<dim>(coeff));

            }

          }
          myfile.close();
          return;
        };
      }
    }
  }

template void getRefinement (const std::string &filename,
  const std::string &header,
  std::vector<std::pair<unsigned int, RefinementCase<1> > > &out,
  std::vector<unsigned int> &nCol,
  const int height);
template void getRefinement (const std::string &filename,
  const std::string &header,
  std::vector<std::pair<unsigned int, RefinementCase<2> > > &out,
  std::vector<unsigned int> &nCol,
  const int height);
template void getRefinement (const std::string &filename,
  const std::string &header,
  std::vector<std::pair<unsigned int, RefinementCase<3> > > &out,
  std::vector<unsigned int> &nCol,
  const int height);

/* Function to parse a multi-line vector
 *  Also completes in nCol how many numbers per line there are.
 */
void parse_multiline_vector (std::ifstream &input_file,
  unsigned int nlines,
  std::map<types::material_id, double> &out,
  bool deprecateFirst)
{
  std::cout << "we are using pareseMultilineVector without template..."
            << std::endl;
// TODO Allow blank lines and comments here!
  std::string line;
  std::string bin;
  double num;
  int k = 0;

  for (unsigned int j = 0; j < nlines; j++)
  {
    getline(input_file, line);
    trim(line);
    std::istringstream iss(line);

    if (deprecateFirst == true)
    {
      iss >> bin;
      if (bin == "PLANE" or bin == "plane" or bin == "Plane")
      {
        j--;
        continue;
      }
      else if ((bin == "#") or (bin == "!") or (bin == "//")
               or bin == "")
      {
        j--;
        continue;
      }
      else
      {
        std::istringstream is2(bin);
        unsigned int line_num = 9999999;
        is2 >> line_num;
        AssertIndexRange(line_num - 1, nlines);
        Assert(line_num > 0,
          ExcMessage("Not valid header in XSEC File: " + bin));
      }
    }

    while (iss.good())
    {
      iss >> num;
      if (num != 999999)
        out[k] = num;
      num = 999999;
      k++;
    }

  }
  return;
}

template <class number>
  void parse_multiline_vector (std::ifstream &input_file,
    unsigned int nlines,
    std::vector<number> &out,
    bool deprecateFirst)
  {
    std::string line;
    std::string bin;
    std::string bin2;
    number num;
    unsigned int line_num;

    for (unsigned int j = 0; j < nlines; j++)
    {
      get_new_valid_line(input_file, line);
      trim(line);
      std::istringstream iss(line);
      iss >> bin;
      trim(bin);
      if (lower_case(bin) == "plane")
      {
        j--;
        continue;
      }

      if (deprecateFirst == true)
      {
        std::istringstream is2(bin);
        is2 >> bin2;
        Assert(!is2.fail(),
          ExcMessage("Error in line numbers, no line numbers? "));
        line_num = Utilities::string_to_int(bin2);
        AssertIndexRange(line_num - 1, nlines);
        AssertRelease(line_num > 0,
          "Not valid header in XSEC File: " + bin);
        iss >> bin;
        Assert(!iss.fail(), ExcMessage("Error in XSEC File"));
      }
      while (iss.good())
      {
        num = str_to_num<number>(bin);
        out.push_back(num);

        iss >> bin;
        trim(bin);
      }
      num = str_to_num<number>(bin);
      out.push_back(num);
    }
    return;
  }

template void parse_multiline_vector<unsigned int> (std::ifstream &input_file,
  unsigned int nlines,
  std::vector<std::vector<unsigned int> > &out,
  bool deprecateFirst);

template void parse_multiline_vector<double> (std::ifstream &input_file,
  unsigned int nlines,
  std::vector<std::vector<double> > &out,
  bool deprecateFirst);

/*
 * Only defined specializations of these:
 */
// TODO change to parse_multiline_vector
template <class number>
  void parse_multiline_vector (std::ifstream &input_file,
    unsigned int nlines,
    std::vector<std::vector<number> > &out,
    bool deprecateFirst)
  {
    std::string line;
    std::string bin;
    std::string bin2;
    number num;
    std::vector<number> vec_num;
    unsigned int line_num;

    for (unsigned int j = 0; j < nlines; j++)
    {
      getline(input_file, line);
      trim(line);
      std::istringstream iss(line);
      iss >> bin;
      trim(bin);

      vec_num.clear();

      if ((bin == "#") or (bin == "!") or (bin == "//") or bin == "")
      {
        j--;
        continue;
      }

      if (deprecateFirst == true)
      {
        std::istringstream is2(bin);
        is2 >> bin2;
        Assert(!is2.fail(),
          ExcMessage("Error in line numbers, no line numbers? "));
        line_num = Utilities::string_to_int(bin2);
        AssertIndexRange(line_num - 1, nlines);
        AssertRelease(line_num > 0,
          "Not valid header in XSEC File: " + bin);
        iss >> bin;
        Assert(!iss.fail(), ExcMessage("Error in XSEC File"));
      }
      while (iss.good())
      {
        num = str_to_num<number>(bin);
        vec_num.push_back(num);
        iss >> bin;
        trim(bin);
      }
      num = str_to_num<number>(bin);
      vec_num.push_back(num);
      out.push_back(vec_num);
    }
    return;
  }

template void parse_multiline_vector<unsigned int> (std::ifstream &input_file,
  unsigned int nlines,
  std::vector<unsigned int> &out,
  bool deprecateFirst);

template void parse_multiline_vector<double> (std::ifstream &input_file,
  unsigned int nlines,
  std::vector<double> &out,
  bool deprecateFirst);

/**
 * Convert string to lower case
 */
std::string lower_case (std::string &data)
{
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  return data;
}

/**
 * Convert const char* to lower case
 */
std::string lower_case (const char *data)
{
  std::string str(data);
  return lower_case(str);
}

/**
 *
 */
void getVec (const std::string filename,
  const std::string header,
  std::vector<double> &out,
  const unsigned int length)
{
  std::ifstream myfile(filename.c_str(), std::ios::in);
  std::string line;
  std::string number;
  double n;
  if (myfile.is_open())
  {
    while (myfile.good())
    {
      myfile >> line;
      if (line == header)
      {
        for (unsigned int i = 0; i < length; i++)
        {
          myfile >> number;
          if (*number.rbegin() == '*')
          {
            number = number.substr(0, number.size() - 1);
            myfile >> n;
            unsigned int len = Utilities::string_to_int(number);
            assert(i + len <= length);
            for (unsigned int j = 0; j < len; j++)
              out.push_back(n);

            i += len - 1;

          }
          else
            out.push_back(Utilities::string_to_double(number));
        }
        myfile.close();
        return;
      }
    }
    return;
  }
  return;
}

double getDouble (const std::string &filename,
  const std::string &header,
  double def)
{
  std::vector<double> out;
  getVec(filename, header, out, 1);
  if (out.size() == 1)
    return out[0];
  else
    return def;

}

template <int dim>
  std::vector<unsigned int> vectorize (const Point<dim> p1,
    const Point<dim> p2)
  {
// make a vector  = Points<dim> p1 - p2
// useful for addRectangle.
    std::vector<unsigned int> n;
    for (unsigned int i = 0; i < dim; i++)
    {
      n.push_back(uint(std::abs(p2(i) - p1(i))));
    }

    return (n);
  }

//
// addRectangle()
//   Add a rectangle to a triangulation. It works for 2D and 3D.
template <int dim>
  void addRectangle (Triangulation<dim> &tri_in,
    Point<dim> p1,
    Point<dim> p2,
    double Lx,
    double Ly,
    double Lz)
  {
    Triangulation<dim> tri_new;
    std::vector<unsigned int> n_grid = vectorize(p1, p2);
    p1[0] *= Lx;
    p2[0] *= Lx;
    p1[1] *= Ly;
    p2[1] *= Ly;
    if (dim == 3)
    {
      p1[2] *= Lz;
      p2[2] *= Lz;
    }

    if (tri_in.n_active_cells() != 0)
    {
      GridGenerator::subdivided_hyper_rectangle(tri_new, n_grid, p1, p2);
      GridGenerator::merge_triangulations(tri_in, tri_new, tri_in);
    }
    else
      GridGenerator::subdivided_hyper_rectangle(tri_in, n_grid, p1, p2);
  }

//
// print_grid(Triangulation<dim> &tri, std::string filename)
//   Print a Grid to a .eps file.
template <int dim>
  void print_grid (Triangulation<dim> &tri,
    std::string filename)
  {
    Assert(dim !=1, ExcImpossibleInDim(dim));
    std::ofstream out(filename.c_str());
    GridOut grid_out;

    grid_out.write_eps(tri, out);
  }

template void print_grid<1> (Triangulation<1> &tri,
  std::string filename);
template void print_grid<2> (Triangulation<2> &tri,
  std::string filename);
template void print_grid<3> (Triangulation<3> &tri,
  std::string filename);

// Prints a std::vector to the console.
// Space separted values. One line.
template <class num>
  void print_vector (std::vector<num> vect,
    bool end_line)
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      std::cout << vect[i] << " ";
    if (end_line and vect.size() > 0)
      std::cout << std::endl;
  }
template void print_vector (std::vector<double> vect,
  bool end_line);
template void print_vector (std::vector<unsigned int> vect,
  bool end_line);

// Prints a std::vector to the console.
// Space separted values. One line.
template <class VECTOR>
  void print_vector (VECTOR vect,
    bool end_line)
  {

    for (unsigned int i = 0; i < vect.size(); i++)
      std::cout << vect[i] << " ";
    if (end_line and vect.size() > 0)
      std::cout << std::endl;
  }

template void print_vector (Vector<double> vect,
  bool end_line);
template void print_vector (Vector<float> vect,
  bool end_line);
template void print_vector (BlockVector<double> vect,
  bool end_line);

/**
 * @brief Prints a std::vector to the console.
 * Space separated values. One line.
 */
void print_vector (PETScWrappers::MPI::Vector vect,
  bool end_line)
{
  for (unsigned int i = 0; i < vect.size(); i++)
    std::cout << vect[i] << " ";
  if (end_line and vect.size() > 0)
    std::cout << std::endl;
}

float percentile (Vector<float> v,
  const float percent)
// Calculate the percentile of an unsorted vector v. In other words,
// value that is placed in the percent place of a sorted vector.
{

  int place = int(v.size() * percent);
  std::nth_element(v.begin(), v.begin() + place, v.end());

  return v[place];
}

unsigned int findRow (const std::vector<unsigned int> nCol,
  const int index,
  int &col,
  int &sum)
// Some utilities to change the indexes when it is used Gauss4 quadrature.
// It computes the row and the colum given...
{
  for (unsigned int row = 0; row < nCol.size(); row++)
  {
    sum += nCol[row];
    if (index < sum)
    {
      if (row == 0)
        col = index;
      else
        col = index - (sum - nCol[row - 1]);
      sum -= nCol[row];
      return row;
    }
  }
  Assert(false, ExcInternalError());
  return 1;
}

//
// direction. The cut direction indicates the in which dimenension it is defined
//  the side
//
template <int dim>
  unsigned int getSide (const unsigned int cut_dir,
    const Point<dim> &point)
  {
    if (point[cut_dir] > 0.5)
      return 1;
    return 0;
  }

/*
 * Returns the side 0 or 1 of the point in the reference cell depending of the cut
 * direction. The cut direction indicates the in which dimenension it is defined
 *  the side
 */
template <int dim>
  unsigned int getSide (const unsigned int cut_dir,
    const Point<dim> &point,
    const unsigned int N)
  {
    double inv = 1.0 / N;
    for (unsigned int i = N - 1; i != 0; i--)
      if (point[cut_dir] > (i * inv))
        return i;
    return 0;

  }

template <int dim>
  void makeGeometyPoints (
    const std::vector<std::pair<Point<dim>, unsigned int> > &centers,
    const Triangulation<dim> &tria,
    std::vector<unsigned int> &geometry_points,
    std::vector<unsigned int> &nRods_out)
  /*
   * makeGeometyPoints(const std::vector<double> centers,
   * const Triangulation<dim> &tria, std::vector<unsigned int> &geometry_points)
   *
   * Make the geometry points from a specific triangulation, and its vector of cell
   * centers.
   * Todo: change the set like implementation to a vector with set entrance
   * implementation. Or maybe a copy-it to a vector.
   *
   */
  {
// Produce an ordered set of the x_center
    std::set<double> x_center, y_center, z_center;

    for (unsigned int k = 0; k < centers.size(); ++k)
    {
      x_center.insert(round(centers[k].first[0]));
      y_center.insert(round(centers[k].first[1]));
      if (dim == 3)
        z_center.insert(round(centers[k].first[2]));
    }

    std::set<double>::iterator it;
    double y, y_ant, z_ant;
    unsigned int in, fin;

// Get the first center
    it = x_center.find(rint(centers[0].first[0]));
    in = std::distance(x_center.begin(), it) + 1;
    y_ant = tria.begin_active()->center()[1];
    if (dim == 3)
      z_ant = tria.begin_active()->center()[2];
    geometry_points.push_back(in);

    for (unsigned int i = 0; i < centers.size(); i++)
    {
      y = centers[i].first[1];

      // If there is a change in the plane break!
      if ((dim == 3))
        if (z_ant != centers[i].first[2])
          break;
      //Look for changes in y
      if (fabs(y - y_ant) > 0.1)
      {

        it = x_center.find(round(centers[i - 1].first[0]));
        fin = std::distance(x_center.begin(), it) + 1;
        geometry_points.push_back(fin);

        it = x_center.find(round(centers[i].first[0]));
        in = std::distance(x_center.begin(), it) + 1;
        geometry_points.push_back(in);

      }
      y_ant = y;
    }
    // Don't forget the last one
    it = x_center.find(round(centers[centers.size() - 1].first[0]));
    fin = std::distance(x_center.begin(), it) + 1;
    geometry_points.push_back(fin);

    // Creation of nRods_out:
    nRods_out.resize(3);
    nRods_out[0] = x_center.size();
    nRods_out[1] = y_center.size();
    if (dim == 2)
      nRods_out[2] = 1;
    else
      nRods_out[2] = z_center.size();
  }

/**
 *
 */
template <class num>
  void print_table (const Table<2, num> tab)
  {
    for (unsigned int j = 0; j < tab.size(1); ++j)
      for (unsigned int i = 0; i < tab.size(0); ++i)
      {
        std::cout << (num) (tab[i][j]) << " ";
        std::cout << "\n";
      }

  }

template void print_table<double> (const Table<2, double> tab);

/**
 *  Parse a vector in a file after a headline.
 *  The expected size of the vector is not checked
 */
void parse_vector_in_file (const std::string &file,
  const std::string &headline,
  std::vector<double> &vector_out,
  const unsigned int n_lines,
  const unsigned int expected_vector_size)
{
  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;
  bool flag = false;

  if (expected_vector_size != static_cast<unsigned int>(-1))
    vector_out.reserve(expected_vector_size);

  // for every line
  for (std::string line; getline(input, line);)
  {
    trim(line);

    // First definition Material and XSecs:
    if (line == headline)
    {
      parse_multiline_vector<double>(input, n_lines, vector_out, false);
      flag = true;
      break;
    }
  }
  input.close();

  AssertRelease(flag, "Headline: " + headline + " Not Found");
  if (expected_vector_size != static_cast<unsigned int>(-1))
    AssertRelease(expected_vector_size == vector_out.size(),
      "Error in vector expected size");
}

/**
 *  Parse a vector in a file after a headline.
 *  The expected size of the vector is not checked
 */
void parse_vector_in_file (const std::string &file,
  const std::string &headline,
  std::vector<unsigned int> &vector_out,
  const unsigned int n_lines,
  const unsigned int expected_vector_size)
{
  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;
  bool flag = false;

  if (expected_vector_size != static_cast<unsigned int>(-1))
    vector_out.reserve(expected_vector_size);

  // for every line
  for (std::string line; getline(input, line);)
  {
    trim(line);

    // First definition Material and XSecs:
    if (line == headline)
    {
      parse_multiline_vector<unsigned int>(input, n_lines, vector_out, false);
      flag = true;
      break;
    }
  }
  input.close();

  AssertRelease(flag, "Headline: " + headline + " Not Found");
  if (expected_vector_size != static_cast<unsigned int>(-1))
    AssertRelease(expected_vector_size == vector_out.size(),
      "Error in vector expected size");
}

void parse_file (unsigned int mat,
  std::string file,
  std::string intro,
  std::vector<double> &D1v,
  std::vector<double> &D2v,
  std::vector<double> &SigmaA1v,
  std::vector<double> &SigmaA2v,
  std::vector<double> &Sigma12v,
  std::vector<double> &SigmaF1v,
  std::vector<double> &SigmaF2v)
{
  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;

  // for every line
  for (std::string line; getline(input, line);)
  {

    trim(line);

    // First definition Material and XSecs:
    if (line == intro)
    {
      std::vector<double> Xsec;
      Xsec.reserve(7);
      parse_multiline_vector(input, 1, Xsec, false);
      Assert(Xsec.size()==7, ExcMessage("problem"));

      D1v[mat] = Xsec[0];
      D2v[mat] = Xsec[1];
      SigmaA1v[mat] = Xsec[2];
      SigmaA2v[mat] = Xsec[3];
      Sigma12v[mat] = Xsec[4];
      SigmaF1v[mat] = Xsec[5];
      SigmaF2v[mat] = Xsec[6];

      break;
    }
  }
  input.close();
}

void parse_file (unsigned int mat,
  std::string file,
  std::string intro,
  std::vector<double> &dfs_left_0,
  std::vector<double> &dfs_right_0,
  std::vector<double> &dfs_left_1,
  std::vector<double> &dfs_right_1)
{
// TODO Hacer esta funcion un poco mas general no tan especializada.

  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;

// for every line
  for (std::string line; getline(input, line);)
  {
    trim(line);
    // First definition Material and XSecs:
    if (line == intro)
    {
      std::vector<double> dfs;
      dfs.reserve(4);
      parse_multiline_vector(input, 1, dfs, false);
      Assert(dfs.size()==4, ExcMessage("problem"))
      dfs_left_0[mat] = dfs[0];
      dfs_right_0[mat] = dfs[1];
      dfs_left_1[mat] = dfs[2];
      dfs_right_1[mat] = dfs[3];

      break;
    }
  }
  input.close();
}

void parse_file (unsigned int,
  std::string file,
  std::string intro,
  std::vector<double> &df_faces,
  unsigned int n)
{
  Assert(fexists(file),
    ExcMessage("Required file " + file + " does not exist."));
  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;

  // for every line
  for (std::string line; getline(input, line);)
  {
    trim(line);
    // First definition Material and XSecs:
    if (line == intro)
    {
      std::vector<double> dfs;
      dfs.reserve(n);
      parse_multiline_vector(input, 1, dfs, false);
      Assert(dfs.size() == n,
        ExcMessage("In file " + file +" , "+ intro + " there are not " + Utilities::int_to_string(n) + " elements, there are " + Utilities::int_to_string(dfs.size()) + "."));

      for (unsigned int i = 0; i < df_faces.size(); i++)
        df_faces[i] = dfs[i];

      return;
    }
  }
  input.close();
  Assert(false, ExcMessage(intro + " does not found in file " + file));
}

void parse_file (std::string file,
  std::string intro,
  std::vector<std::vector<double> > &dfs_poly,
  unsigned int n)
{
  Assert(fexists(file),
    ExcMessage("Required file " + file + " does not exist."));
  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;

// for every line
  for (std::string line; getline(input, line);)
  {
    trim(line);
    // First definition Material and XSecs:
    if (line == intro)
    {
      parse_multiline_vector(input, n, dfs_poly, false);
      Assert(dfs_poly.size() == n,
        ExcMessage("In file " + file +" , "+ intro + " there are not " + Utilities::int_to_string(n) + " elements, there are " + Utilities::int_to_string(dfs_poly.size()) + "."));
      return;
    }
  }
  input.close();
  Assert(false, ExcMessage(intro + " does not found in file " + file));
}

void parse_file (unsigned int,
  std::string file,
  std::string intro,
  std::vector<double> &shape_fun,
  unsigned int n,
  unsigned int nlines)
{
  Assert(fexists(file),
    ExcMessage("Required file " + file + " does not exist."));
  std::ifstream input(file.c_str(), std::ios::in);
  std::string sub, Name;

// for every line
  for (std::string line; getline(input, line);)
  {
    trim(line);
    // First definition Material and XSecs:
    if (line == intro)
    {
      std::vector<double> vect;
      vect.reserve(n);
      parse_multiline_vector(input, nlines, vect, false);
      Assert(vect.size() == n,
        ExcMessage("In file " + file +", "+ intro + " there are not " + Utilities::int_to_string(n) + " elements, there are " + Utilities::int_to_string(vect.size()) + "."));

      for (unsigned int i = 0; i < shape_fun.size(); i++)
        shape_fun[i] = vect[i];

      return;
    }
  }
  input.close();
  Assert(false, ExcMessage(intro + " does not found in file " + file));
}

/* std::vector<unsigned int> default_geometry_points(td::vector<unsigned int> n_cells_per_dim)
 * This function construct the default rectangular geometry_points.
 * E.G 2D, 9 by 4 reactor:
 * 1 9 1 9 1 9 1 9 1 9
 */
std::vector<unsigned int> default_geometry_points (
  std::vector<unsigned int> n_cells_per_dim)
{
  std::vector<unsigned int> out(n_cells_per_dim[1] * 2);
  for (unsigned int i = 0; i < n_cells_per_dim[1]; i++)
  {
    out[2 * i] = 1;
    out[2 * i + 1] = n_cells_per_dim[0];
  }
  return out;
}

/*
 * Parse a matrix from a space and '\n' separated string.
 * n_rows and n_cols can be set to ensure a matrix size
 */
template <typename num>
  void parse_matrix (const std::string &input,
    std::vector<std::vector<num> > &out_matrix,
    const unsigned int n_rows,
    const unsigned int n_cols)
  {
    std::string in = input;
    boost::trim_if(in, boost::is_any_of(";\n "));
    std::vector<std::string> strs;
    boost::split(strs, in, boost::is_any_of(";\n"), boost::token_compress_on);

    AssertRelease(n_rows == 0 or strs.size() == n_rows,
      "Invalid number of rows");
    out_matrix.resize(strs.size());
    for (unsigned int row = 0; row < strs.size(); ++row)
    {
      str_to_vector(strs[row], out_matrix[row]);
      AssertRelease(n_cols == 0 or out_matrix[row].size() == n_cols,
        "Invalid number of cols");
    }
  }

template void parse_matrix (const std::string &input,
  std::vector<std::vector<int> > &out_matrix,
  const unsigned int n_rows,
  const unsigned int n_cols);

template void parse_matrix (const std::string &input,
  std::vector<std::vector<unsigned int> > &out_matrix,
  const unsigned int n_rows,
  const unsigned int n_cols);

template void parse_matrix (const std::string &input,
  std::vector<std::vector<double> > &out_matrix,
  const unsigned int n_rows,
  const unsigned int n_cols);
/**
 * @brief Parse a vector from an string. Check at the end if the given vector have
 * the expected length. If the string is empty return the default input.
 */
void parse_vector (std::string input,
  std::vector<unsigned int> &out,
  unsigned int length,
  std::vector<unsigned int> def)
{
  Assert(out.size()==0, ExcMessage("Vector not empty"));
  trim(input);
  std::istringstream iss(input);
  unsigned int number;
  std::string str;
  out.reserve(length);

// if input is empty return default parameter
  if (input.empty())
  {
    out = def;
    return;
  }

  while (iss.good())
  {
    iss >> str;

    // Allow "4* 10" notation for repetitions
    if (*str.rbegin() == '*')
    {
      str = str.substr(0, str.size() - 1);
      iss >> number;
      unsigned int len = Utilities::string_to_int(str);
      for (unsigned int j = 0; j < len; j++)
        out.push_back(number);
    }
    else
      out.push_back(Utilities::string_to_int(str));
  }

  if (length != 0)
    AssertRelease(out.size() == length,
      "There aren't the number of arguments it should be \n");

  return;
}

/**
 * Same as before but for double vector. There is not default values in this case.
 * It exist an special character in order to repeat number:
 * 4* 10.0 = 10.0 10.0 10.0 10.0
 * 2.0 2*1.0 2.0 = 2.0 1.0 1.0 2.0
 */
void parse_vector (std::string input,
  std::vector<double> &out,
  unsigned int length)
{
  Assert(out.size()==0, ExcMessage("Vector not empty"));
  trim(input);
  std::istringstream iss(input);

  double number;
  std::string str;
  out.reserve(length);
  while (iss.good())
  {
    iss >> str;

    // Allow "4* 10.0" notation for repetition
    if (*str.rbegin() == '*')
    {
      str = str.substr(0, str.size() - 1);
      iss >> number;
      unsigned int len = Utilities::string_to_int(str);
      for (unsigned int j = 0; j < len; j++)
        out.push_back(number);
    }
    else
      out.push_back(Utilities::string_to_double(str));
  }

  if (length != 0)
    AssertRelease(out.size() == length,
      "There aren't the number of arguments it should be \n");

  return;
}

/**
 * Same as before but for double vector. There is not default values in this case.
 * It exist an special character in order to repeat number:
 * 4* 10.0 = 10.0 10.0 10.0 10.0
 * 2.0 2*1.0 2.0 = 2.0 1.0 1.0 2.0
 */
void parse_vector (std::string input,
  PETScWrappers::MPI::BlockVector &out,
  unsigned int n_blocks,
  unsigned int n_dofs_per_block)
{
  Assert(out.size()==0, ExcMessage("Vector not empty"));
  trim(input);
  std::istringstream iss(input);

  double num;
  out.reinit(n_blocks, PETSC_COMM_WORLD,  n_dofs_per_block, n_dofs_per_block);
  for (unsigned int i=0; i< n_blocks *n_dofs_per_block; i++)
  {
    iss >> num;
     out[i] = num;
  }


  return;
}

/**
 * Same as before but for double vector. There is not default values in this case.
 * It exist an special character in order to repeat number:
 * 4* 10.0 = 10.0 10.0 10.0 10.0
 * 2.0 2*1.0 2.0 = 2.0 1.0 1.0 2.0
 */
void parse_vector (std::string input,
  PETScWrappers::MPI::Vector &out,
  unsigned int n_dofs)
{
  Assert(out.size()==0, ExcMessage("Vector not empty"));
  trim(input);
  std::istringstream iss(input);

  double num;
  out.reinit(PETSC_COMM_WORLD,  n_dofs, n_dofs);
  for (unsigned int i=0; i< n_dofs; i++)
  {
    iss >> num;
     out[i] = num;
  }


  return;
}

template <int dim>
  void extrude_triangulation (const Triangulation<2, 2>&,
    const unsigned int,
    const double,
    Triangulation<dim>&,
    const bool)
  {
    AssertRelease(false,
      "Impossible dim" + num_to_str(dim) + "for extrude_triangulation");
  }

//
//
//
template <>
  void extrude_triangulation (const Triangulation<2, 2> &input,
    const unsigned int n_slices,
    const double height,
    Triangulation<3> &result,
    const bool changeMat)
  {
    Assert(input.n_levels() == 1,
      ExcMessage ("The input triangulation must be coarse meshes."));
    Assert(result.n_cells()==0,
      ExcMessage("Result triangulation need to be empty upon calling extrude_triangulation."));
    Assert(height>0,
      ExcMessage("The height in extrude_triangulation needs to be positive."));
    Assert(n_slices>=2,
      ExcMessage("The number of slices in extrude_triangulation needs to be at least 2."));

    std::vector<Point<3> > points(n_slices * input.n_vertices());
    std::vector<CellData<3> > cells;
    cells.reserve((n_slices - 1) * input.n_active_cells());

// Get the maximum material number
    unsigned int n_materials = 0;
    for (Triangulation<2, 2>::cell_iterator cell = input.begin();
        cell != input.end(); ++cell)
    {
      if (cell->material_id() > n_materials)
        n_materials = cell->material_id();
    }

    for (unsigned int slice = 0; slice < n_slices; ++slice)
    {
      for (unsigned int i = 0; i < input.n_vertices(); ++i)

      {
        const Point<2> &v = input.get_vertices()[i];
        points[i + slice * input.n_vertices()](0) = v(0);
        points[i + slice * input.n_vertices()](1) = v(1);
        points[i + slice * input.n_vertices()](2) = height * slice
                                                    / (n_slices - 1);
      }
    }

    for (Triangulation<2, 2>::cell_iterator cell = input.begin();
        cell != input.end(); ++cell)
    {
      for (unsigned int slice = 0; slice < n_slices - 1; ++slice)
      {
        CellData<3> this_cell;
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell;
            ++v)
        {
          this_cell.vertices[v] = cell->vertex_index(v)
                                  + slice * input.n_vertices();
          this_cell.vertices[v + GeometryInfo<2>::vertices_per_cell] =
              cell->vertex_index(v)
              + (slice + 1) * input.n_vertices();
        }

        if (changeMat == true)
          this_cell.material_id = cell->material_id()
                                  + slice * n_materials;
        else
          this_cell.material_id = cell->material_id();

        cells.push_back(this_cell);
      }
    }

    SubCellData s;
// Put every boundary indicator to 1
    types::boundary_id bid = 1;
    s.boundary_quads.reserve(
      input.n_active_lines() * (n_slices - 1)
      + input.n_active_cells() * 2);
    for (Triangulation<2, 2>::cell_iterator cell = input.begin();
        cell != input.end(); ++cell)
    {
      CellData<2> quad;
      for (unsigned int f = 0; f < 4; ++f)
        if (cell->at_boundary(f))
        {
          quad.boundary_id = bid;
          for (unsigned int slice = 0; slice < n_slices - 1; ++slice)
          {
            quad.vertices[0] = cell->face(f)->vertex_index(0)
                               + slice * input.n_vertices();
            quad.vertices[1] = cell->face(f)->vertex_index(1)
                               + slice * input.n_vertices();
            quad.vertices[2] = cell->face(f)->vertex_index(0)
                               + (slice + 1) * input.n_vertices();
            quad.vertices[3] = cell->face(f)->vertex_index(1)
                               + (slice + 1) * input.n_vertices();
            s.boundary_quads.push_back(quad);
          }
        }
    }

    for (Triangulation<2, 2>::cell_iterator cell = input.begin();
        cell != input.end(); ++cell)
    {
      CellData<2> quad;
      quad.boundary_id = bid;
      quad.vertices[0] = cell->vertex_index(0);
      quad.vertices[1] = cell->vertex_index(1);
      quad.vertices[2] = cell->vertex_index(2);
      quad.vertices[3] = cell->vertex_index(3);
      s.boundary_quads.push_back(quad);

      quad.boundary_id = bid;
      for (int i = 0; i < 4; ++i)
        quad.vertices[i] += (n_slices - 1) * input.n_vertices();
      s.boundary_quads.push_back(quad);
    }

    result.create_triangulation(points, cells, s);
  }

template void extrude_triangulation (const Triangulation<2, 2> &input,
  const unsigned int n_slices,
  const double height,
  Triangulation<1> &result,
  const bool changeMat);
template void extrude_triangulation (const Triangulation<2, 2> &input,
  const unsigned int n_slices,
  const double height,
  Triangulation<2> &result,
  const bool changeMat);

/**
 *
 */
void parse_materials (std::string xsec_file,
  std::vector<unsigned int> &n_cells_per_dim,
  unsigned int,
  std::vector<unsigned int> &materials)
{
  Assert(fexists(xsec_file), ExcMessage("XECSFile doesn't exist"));
  Assert(n_cells_per_dim.size() == 3, ExcMessage("Error in n_cells_per_dim"));

  std::ifstream input(xsec_file.c_str(), std::ios::in);
  std::string sub;
  std::string keyword;

  // for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    keyword.clear();
    iss >> keyword;
    trim(keyword);

    if (is_commentary(keyword))
      continue;

    // First definition Material and XSecs:
    else if (keyword == "Materials")
      parse_multiline_vector(input, n_cells_per_dim[1] * n_cells_per_dim[2],
        materials, true);
  }
}

/**
 * This function does a binary search and returns the index i such that
 *  x is between wl[i] and wl[i+1], except that i is restricted to the
 *   range from 0 to n-2 inclusive
 */
unsigned int binary_search (double x,
  std::vector<double> wl,
  unsigned int n)
{
  unsigned int low, high; /* range of elements to consider */
  unsigned int k; /* middle element between low and high */

  low = 0;
  high = n - 2;

  while (high - low > 1)
  {
    k = (low + high) / 2;
    if (x < wl[k])
    {
      high = k;
    }
    else if (x > wl[k + 1])
    {
      low = k;
    }
    else
    {
      high = k;
      low = k;
    }
  }

  if (low < high)
  {
    if (wl[high] <= x)
      return (high);
    else
      return (low);
  }
  else
  {
    return (low);
  }
}

/**
 *  This routine implements a linear interpolation.
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
 *
 */
double interp_linear (double x,
  std::vector<double> wl,
  std::vector<double> f,
  unsigned int n,
  unsigned int &starti)
{

  double p0, p1; /* differences between x and wl elements */
  double y; /* interpolated value */
  unsigned int i; /* index for which wl[i] is closest to x */
  const double tol_round = 1e-12;

  if (n == 1 || std::abs(x - wl[0]) < tol_round)
  {
    y = f[0];
  }
  else if (std::abs(x - wl[n - 1]) < tol_round)
  {
    y = f[n - 1];
  }
  else if (x > wl[n - 1])
  {
    AssertRelease(false,
      "interp_linear is extrapolating in the upper bound\n x: "
      + num_to_str(x)
      + " bound: " + num_to_str(wl[n - 1]));
    y = 0.0;
  }
  else if (x < wl[0])
  {
    AssertRelease(false,
      "interp_linear is extrapolating in the lower bound\n x: "
      + num_to_str(x)
      + " bound: " + num_to_str(wl[0]));
    y = 0.0;
  }
  else
  {
    //Find the index i such that x is between wl[i] and wl[i+1];
    // i can have any value from zero to n-2, inclusive.

    i = starti;

    if (i > n - 2)
    { /* i is out of range? */

      i = binary_search(x, wl, n);

    }
    else if (i > 0 && wl[i - 1] <= x && x < wl[i])
    {

      i--; /* x is in the previous interval */

    }
    else if (i < n - 2 && wl[i + 1] <= x && x <= wl[i + 2])
    {

      i++; /* x is in the next interval */

    }
    else if (x < wl[i] || x > wl[i + 1])
    {

      i = binary_search(x, wl, n);

    }

    // Update the starting index, in case we call this function
    // again.
    starti = i;

    // linear interpolation
    p0 = (wl[i + 1] - x) / (wl[i + 1] - wl[i]);
    p1 = 1. - p0;
    y = f[i] * p0 + f[i + 1] * p1;
  }

  return (y);
}

/**
 * Compute max and min vertex of cell in the coord coordinate.
 */
template <int dim>
  void getMaxMinVertex (
    TriaIterator<DoFCellAccessor<dim, dim, false> > cell,
    unsigned int coord,
    double &maxp,
    double &minp)
  {
    Assert(coord < dim, ExcImpossibleInDim(dim));
    maxp = cell->vertex(0)(coord);
    minp = cell->vertex(0)(coord);

    for (unsigned int v = 1; v < GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      maxp = std::max(maxp, cell->vertex(v)[coord]);
      minp = std::min(minp, cell->vertex(v)[coord]);
    }
    Assert(maxp >= 0.0, ExcMessage(" "))
    Assert(minp >= 0.0, ExcMessage(" "))
  }

template void getMaxMinVertex (
  TriaIterator<DoFCellAccessor<1, 1, false>  > cell,
  unsigned int coord,
  double &maxp,
  double &minp);
template void getMaxMinVertex (
  TriaIterator<DoFCellAccessor<2, 2, false>  > cell,
  unsigned int coord,
  double &maxp,
  double &minp);
template void getMaxMinVertex (
  TriaIterator<DoFCellAccessor<3, 3, false>  > cell,
  unsigned int coord,
  double &maxp,
  double &minp);

// make_partial_sparsity_pattern
//
//
template <class DH, class SparsityPattern>
  void make_partial_sparsity_pattern (const DH &dof_handler,
    SparsityPattern &sparsity,
    const std::vector<types::global_dof_index> &map_global_to_partial_i,
    unsigned int i_start,
    unsigned int i_end,
    const std::vector<types::global_dof_index> &map_global_to_partial_j,
    unsigned int j_start,
    unsigned int j_end)
  {

    AssertDimension(map_global_to_partial_i.size(), dof_handler.n_dofs());
    AssertDimension(map_global_to_partial_j.size(), dof_handler.n_dofs());

    const unsigned int n_dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> cell_dofs(n_dofs_per_cell);

    // Loop over all cells. And allocate the sparsity memory.
    // Note that if a place in the sparsity is allocate twice nothing bad happens.
    typename DH::active_cell_iterator cell = dof_handler.begin_active(),
        end_cell = dof_handler.end();
    for (; cell != end_cell; ++cell)
    {
      cell->get_dof_indices(cell_dofs);

      // Make sparsity pattern for this cell
      for (unsigned int i = i_start; i < i_end; ++i)
        for (unsigned int j = j_start; j < j_end; ++j)
          sparsity.add(map_global_to_partial_i[cell_dofs[i]],
            map_global_to_partial_j[cell_dofs[j]]);
    }
  }

// make_partial_sparsity_pattern
//
//
template <class DH, class SparsityPattern>
  void make_partial_sparsity_pattern (const DH &dof_handler,
    SparsityPattern &sparsity,
    const std::vector<types::global_dof_index> &map_global_to_partial_i,
    unsigned int i_start,
    unsigned int i_end)
  {

    AssertDimension(sparsity.n_rows(), sparsity.n_cols());

    // Use the other function
    make_partial_sparsity_pattern(dof_handler, sparsity,
      map_global_to_partial_i, i_start, i_end, map_global_to_partial_i,
      i_start, i_end);
  }

template void make_partial_sparsity_pattern (const DoFHandler<1> &dof_handler,
  DynamicSparsityPattern &sparsity,
  const std::vector<types::global_dof_index> &map_global_to_partial_i,
  unsigned int i_start,
  unsigned int i_end,
  const std::vector<types::global_dof_index> &map_global_to_partial_j,
  unsigned int j_start,
  unsigned int j_end);

template void make_partial_sparsity_pattern (const DoFHandler<2> &dof_handler,
  DynamicSparsityPattern &sparsity,
  const std::vector<types::global_dof_index> &map_global_to_partial_i,
  unsigned int i_start,
  unsigned int i_end,
  const std::vector<types::global_dof_index> &map_global_to_partial_j,
  unsigned int j_start,
  unsigned int j_end);

template void make_partial_sparsity_pattern (const DoFHandler<3> &dof_handler,
  DynamicSparsityPattern &sparsity,
  const std::vector<types::global_dof_index> &map_global_to_partial_i,
  unsigned int i_start,
  unsigned int i_end,
  const std::vector<types::global_dof_index> &map_global_to_partial_j,
  unsigned int j_start,
  unsigned int j_end);

template void make_partial_sparsity_pattern (const DoFHandler<1> &dof_handler,
  DynamicSparsityPattern &sparsity,
  const std::vector<types::global_dof_index> &global_to_partial_map,
  unsigned int i_start,
  unsigned int i_end);

template void make_partial_sparsity_pattern (const DoFHandler<2> &dof_handler,
  DynamicSparsityPattern &sparsity,
  const std::vector<types::global_dof_index> &global_to_partial_map,
  unsigned int i_start,
  unsigned int i_end);

template void make_partial_sparsity_pattern (const DoFHandler<3> &dof_handler,
  DynamicSparsityPattern &sparsity,
  const std::vector<types::global_dof_index> &global_to_partial_map,
  unsigned int i_start,
  unsigned int i_end);

template <int dim>
  void map_global_to_partial (const DoFHandler<dim> &dof_handler,
    std::vector<types::global_dof_index> &global_to_partial_map,
    const unsigned int i_start,
    const unsigned int i_end)
  {
    // global_to_partial_map should be empty
    Assert(global_to_partial_map.empty(), ExcMessage("The map should be empty"));

    // Construct the edges_global_maps
    unsigned int local_num = 0;

    global_to_partial_map.resize(dof_handler.n_dofs());
    std::vector<bool> is_done(dof_handler.n_dofs(), false);
    std::vector<types::global_dof_index> cell_dofs(
      dof_handler.get_fe().dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        end_cell = dof_handler.end();
    for (; cell != end_cell; ++cell)
    {
      cell->get_dof_indices(cell_dofs);

      for (unsigned int i = i_start; i < i_end; ++i)
      {
        if (is_done[cell_dofs[i]] == false)
        {
          global_to_partial_map[cell_dofs[i]] = local_num;

          ++local_num;
          is_done[cell_dofs[i]] = true;
        }
      }
    }
  }

template
void map_global_to_partial (const DoFHandler<1> &dof_handler,
  std::vector<types::global_dof_index> &global_to_partial_map,
  const unsigned int i_start,
  const unsigned int i_end);

template
void map_global_to_partial (const DoFHandler<2> &dof_handler,
  std::vector<types::global_dof_index> &global_to_partial_map,
  const unsigned int i_start,
  const unsigned int i_end);

template
void map_global_to_partial (const DoFHandler<3> &dof_handler,
  std::vector<types::global_dof_index> &global_to_partial_map,
  const unsigned int i_start,
  const unsigned int i_end);

template <int dim>
  void map_partial_to_global (const DoFHandler<dim> &dof_handler,
    std::vector<types::global_dof_index> &partial_global_map,
    const unsigned int n_partial_dofs,
    const unsigned int i_start,
    const unsigned int i_end)
  {
    // global_to_partial_map should be empty
    Assert(partial_global_map.empty(), ExcMessage("The map should be empty"));

    // Construct the edges_global_maps
    unsigned int local_num = 0;

    partial_global_map.resize(n_partial_dofs);
    std::vector<bool> is_done(dof_handler.n_dofs(), false);
    std::vector<types::global_dof_index> cell_dofs(
      dof_handler.get_fe().dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        end_cell = dof_handler.end();
    for (; cell != end_cell; ++cell)
    {
      cell->get_dof_indices(cell_dofs);

      for (unsigned int i = i_start; i < i_end; ++i)
      {
        if (is_done[cell_dofs[i]] == false)
        {
          partial_global_map[local_num] = cell_dofs[i];

          ++local_num;
          is_done[cell_dofs[i]] = true;
        }
      }
    }
    AssertDimension(local_num, n_partial_dofs);
  }

template
void map_partial_to_global (const DoFHandler<1> &dof_handler,
  std::vector<types::global_dof_index> &global_to_partial_map,
  const unsigned int n_partial_dofs,
  const unsigned int i_start,
  const unsigned int i_end);

template
void map_partial_to_global (const DoFHandler<2> &dof_handler,
  std::vector<types::global_dof_index> &global_to_partial_map,
  const unsigned int n_partial_dofs,
  const unsigned int i_start,
  const unsigned int i_end);

template
void map_partial_to_global (const DoFHandler<3> &dof_handler,
  std::vector<types::global_dof_index> &global_to_partial_map,
  const unsigned int n_partial_dofs,
  const unsigned int i_start,
  const unsigned int i_end);

//
// map_subpartial_to_partial(std::vector partial_to_global, std::vector partial_to_global,
//    std::vector result partial_to_subpartial)
//    Translate a subpartial_to_global map  to a subpartial_to_partial map using the corresponding
//    map_global_to_partial dictionary.
void map_subpartial_to_partial (
  const std::vector<types::global_dof_index> &map_subpartial_to_global,
  const std::vector<types::global_dof_index> &map_global_to_partial,
  std::vector<types::global_dof_index> &map_subpartial_to_partial)
{
  // The result map should be empty
  Assert(map_subpartial_to_partial.empty(),
    ExcMessage("The result map should be empty"));

  const unsigned int n_subpartials = map_subpartial_to_global.size();
  map_subpartial_to_partial.resize(n_subpartials);
  for (unsigned int i = 0; i < n_subpartials; ++i)
    map_subpartial_to_partial[i] =
                                   map_global_to_partial[map_subpartial_to_global[i]];
}

/**
 * Translate a subpartial_to_global map  to a subpartial_to_partial map using the corresponding
 * map_global_to_partial dictionary..
 */
void map_subpartials_to_partial (
  const std::vector<std::vector<types::global_dof_index> > &map_subpartials_to_global,
  const std::vector<types::global_dof_index> &map_global_to_partial,
  std::vector<std::vector<types::global_dof_index> > &map_subpartials_to_partial,
  std::vector<unsigned int> &partition_of_unity_global,
  std::vector<std::vector<double> > &partition_of_unity)
{
  // The result map should be empty
  Assert(map_subpartials_to_partial.empty(),
    ExcMessage("The result map should be empty"));

  const unsigned int n_subdomains = map_subpartials_to_global.size();
  const unsigned int dofs_per_subdomain = map_subpartials_to_global[0].size();
  map_subpartials_to_partial.resize(n_subdomains,
    std::vector<types::global_dof_index>(dofs_per_subdomain));

  for (unsigned int s = 0; s < n_subdomains; ++s)
    for (unsigned int i = 0; i < dofs_per_subdomain; ++i)
      map_subpartials_to_partial[s][i] =
          map_global_to_partial[map_subpartials_to_global[s][i]];

  // Resizing partitions of unity
  std::vector<unsigned int> partition_of_unity_schur(
    partition_of_unity_global.size(), -1);

  partition_of_unity.resize(n_subdomains);
  for (unsigned int s = 0; s < n_subdomains; ++s)
    partition_of_unity[s].resize(dofs_per_subdomain);

  for (unsigned int i = 0; i < partition_of_unity_global.size(); ++i)
    partition_of_unity_schur[map_global_to_partial[i]] =
                                                         partition_of_unity_global[i];

  // Invert and redistribute
  for (unsigned int s = 0; s < n_subdomains; ++s)
    for (unsigned int i = 0; i < dofs_per_subdomain; ++i)
      partition_of_unity[s][i] =
          1.0
          / partition_of_unity_schur[map_subpartials_to_partial[s][i]];
}

//
// map_edges_to_global ( ... )
//
template <int dim>
  void map_edges_to_global (const DoFHandler<dim> &dof_handler,
    std::vector<std::vector<types::global_dof_index> > &partials_global_map,
    std::vector<unsigned int> &partition_of_unity_vector,
    unsigned int &n_subdomains,
    unsigned int &n_dofs_per_subdomain)
  {
    AssertRelease(false, "Not Valid dim");
  }

//
// map_edges_to_global ( ... )
//
template <>
  void map_edges_to_global (const DoFHandler<2> &dof_handler,
    std::vector<std::vector<types::global_dof_index> > &partials_global_map,
    std::vector<unsigned int> &partition_of_unity_vector,
    unsigned int &n_subdomains,
    unsigned int &n_dofs_per_subdomain)
  {
    // global_to_partial_map should be empty
    Assert(partials_global_map.empty(), ExcMessage("The map should be empty"));

    n_subdomains = dof_handler.get_triangulation().n_active_lines();
    n_dofs_per_subdomain = GeometryInfo<2>::vertices_per_face
                           + dof_handler.get_fe().dofs_per_line;

    // Construct the partials_global_map
    partials_global_map.resize(n_subdomains);
    for (unsigned int s = 0; s < n_subdomains; ++s)
      partials_global_map[s].resize(n_dofs_per_subdomain);

    partition_of_unity_vector.resize(dof_handler.n_dofs(), 0.0);

    std::vector<bool> is_done(n_subdomains, false);
    unsigned int line_index;

    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
        end_cell = dof_handler.end();
    for (; cell != end_cell; ++cell)
    {
      for (unsigned int l = 0; l < GeometryInfo<2>::lines_per_cell; ++l)
      {
        line_index = cell->line(l)->index();

        if (is_done[line_index] == false)
        {
          is_done[line_index] = true;

          cell->line(l)->get_dof_indices(partials_global_map[line_index]);

          for (unsigned int i = 0; i < n_dofs_per_subdomain; ++i)
            partition_of_unity_vector[partials_global_map[line_index][i]] +=
                                                                             1;
        }
      }
    }
  }

//
// void map_edges_to_global ( ... )
//
template <>
  void map_edges_to_global (const DoFHandler<3> &dof_handler,
    std::vector<std::vector<types::global_dof_index> > &partials_global_map,
    std::vector<unsigned int> &partition_of_unity_vector,
    unsigned int &n_subdomains,
    unsigned int &n_dofs_per_subdomain)
  {

    // global_to_partial_map should be empty
    Assert(partials_global_map.empty(), ExcMessage("The map should be empty"));

    n_subdomains = dof_handler.get_triangulation().n_active_quads();
    n_dofs_per_subdomain = GeometryInfo<3>::vertices_per_face
                           + +GeometryInfo<3>::lines_per_face
                             * dof_handler.get_fe().dofs_per_line
                           + dof_handler.get_fe().dofs_per_quad;

    // Construct the partials_global_map
    partials_global_map.resize(n_subdomains);
    for (unsigned int s = 0; s < n_subdomains; ++s)
      partials_global_map[s].resize(n_dofs_per_subdomain, -1);

    partition_of_unity_vector.resize(dof_handler.n_dofs(), 0.0);
    std::vector<bool> is_done(n_subdomains, false);
    unsigned int quad_index;

    DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(),
        end_cell = dof_handler.end();
    for (; cell != end_cell; ++cell)
    {
      for (unsigned int q = 0; q < GeometryInfo<3>::quads_per_cell; ++q)
      {
        quad_index = cell->quad(q)->index();

        if (is_done[quad_index] == false)
        {
          is_done[quad_index] = true;

          cell->quad(q)->get_dof_indices(partials_global_map[quad_index]);

          for (unsigned int i = 0; i < n_dofs_per_subdomain; ++i)
            partition_of_unity_vector[partials_global_map[quad_index][i]] +=
                                                                             1;
        }
      }
    }
  }

template <class VECTOR>
  void print_vector_in_matlab (const VECTOR &vec,
    const std::string &vec_name,
    std::ostream &out,
    const unsigned int precision)
  {
    out << std::setprecision(precision);

    out << vec_name << " = [";
    for (unsigned int i = 0; i < vec.size(); ++i)
      out << vec[i] << " ";

    out << "]';" << std::endl;
  }

template void print_vector_in_matlab (const std::vector<int> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_matlab (const std::vector<unsigned int> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_matlab (const std::vector<double> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_matlab (const Vector<double> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_matlab (const PETScWrappers::MPI::Vector &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template <class BLOCKVECTOR>
  void print_block_vector_in_matlab (const BLOCKVECTOR &vec,
    const std::string &vec_name,
    std::ostream &out,
    const unsigned int precision)
  {
    out << std::setprecision(precision);

    for (unsigned int b = 0; b < vec.n_blocks(); b++)
    {
      std::string name = vec_name + std::to_string(b + 1);
      print_vector_in_matlab(vec.block(b), name, out, 12);
    }

    out << vec_name + "=[";

    for (unsigned int b = 0; b < vec.n_blocks(); b++)
    {
      out << vec_name + std::to_string(b + 1) + ";";
    }

    out << "];" << std::endl;

    for (unsigned int b = 0; b < vec.n_blocks(); b++)
    {
      out << "clear " + vec_name + std::to_string(b + 1) + ";";
    }

  }

template void print_block_vector_in_matlab (const BlockVector<double> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_block_vector_in_matlab (const PETScWrappers::MPI::BlockVector &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template <class VECTOR>
  void print_vector_in_python (const VECTOR &vec,
    const std::string &vec_name,
    std::ostream &out,
    const unsigned int precision)
  {
    out << std::setprecision(precision);

    out << vec_name << " = np.array([";
    for (unsigned int i = 0; i < vec.size() - 1; ++i)
      out << vec[i] << ", ";

    out << vec[vec.size() - 1] << " ";
    out << "]);" << std::endl;
  }

template void print_vector_in_python (const std::vector<int> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_python (const std::vector<unsigned int> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_python (const std::vector<double> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_python (const Vector<double> &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template void print_vector_in_python (const PETScWrappers::MPI::Vector &vec,
  const std::string &vec_name,
  std::ostream &out,
  const unsigned int precision);

template <class MATRIX>
  void print_matrix_in_matlab (const MATRIX &mat,
    const std::string &mat_name,
    std::ostream &out,
    const unsigned int precision)
  {

    out << std::setprecision(precision);

    unsigned int nnz = mat.n_nonzero_elements();
    out << "Nnz = " << nnz << ";" << std::endl;
    std::vector<unsigned int> i_index;
    std::vector<unsigned int> j_index;
    std::vector<double> mat_values;

    i_index.reserve(nnz);
    j_index.reserve(nnz);
    mat_values.reserve(nnz);

    for (unsigned int i = 0; i < mat.n(); ++i)
      for (unsigned int j = 0; j < mat.m(); ++j)
        if (std::abs(mat.el(i, j)) > 1e-15)
        {
          i_index.push_back(i + 1);
          j_index.push_back(j + 1);
          mat_values.push_back(mat.el(i, j));
        }

    // Print in the vectors in the out
    print_vector_in_matlab(i_index, "i_index", out, precision);
    print_vector_in_matlab(j_index, "j_index", out, precision);
    print_vector_in_matlab(mat_values, "mat_values", out, precision);

    out << "m = " << mat.m() << ";" << std::endl;
    out << "n = " << mat.n() << ";" << std::endl;
    out << mat_name << " = sparse(i_index, j_index, mat_values, m, n, Nnz);"
        << std::endl;
    out << "clear  i_index; clear  j_index;" << std::endl;
    out << "clear mat_values; clear m; clear n; clear Nnz;" << std::endl;
    out << std::endl;
  }

template void print_matrix_in_matlab (const SparseMatrix<double> &mat,
  const std::string &mat_name,
  std::ostream &out,
  const unsigned int precision);

template void print_matrix_in_matlab (const PETScWrappers::SparseMatrix &mat,
  const std::string &mat_name,
  std::ostream &out,
  const unsigned int precision);

template void print_matrix_in_matlab (
  const PETScWrappers::MPI::SparseMatrix &mat,
  const std::string &mat_name,
  std::ostream &out,
  const unsigned int precision);

template <class MATRIX>
  void print_matrix_in_python (const MATRIX &mat,
    const std::string &mat_name,
    std::ostream &out,
    const unsigned int precision)
  {

    out << std::setprecision(precision);

    unsigned int nnz = mat.n_nonzero_elements();
    std::vector<unsigned int> i_index;
    std::vector<unsigned int> j_index;
    std::vector<double> mat_values;

    i_index.reserve(nnz);
    j_index.reserve(nnz);
    mat_values.reserve(nnz);

    for (unsigned int i = 0; i < mat.n(); ++i)
      for (unsigned int j = 0; j < mat.m(); ++j)
        if (std::abs(mat.el(i, j)) > 1e-15)
        {
          i_index.push_back(i);
          j_index.push_back(j);
          mat_values.push_back(mat.el(i, j));
        }

    // Print in the vectors in the out
    print_vector_in_python(i_index, "i_index", out, precision);
    print_vector_in_python(j_index, "j_index", out, precision);
    print_vector_in_python(mat_values, "mat_values", out, precision);

    out << "m = " << mat.m() << ";" << std::endl;
    out << "n = " << mat.n() << ";" << std::endl;
    out << mat_name
        << " = csr_matrix((mat_values, (i_index, j_index)), shape=(m,n), dtype='double');"
        << std::endl;
    out << "del  i_index, j_index " << std::endl;
    out << "del mat_values, m, n " << std::endl;
    out << std::endl;
  }

template void print_matrix_in_python (const SparseMatrix<double> &mat,
  const std::string &mat_name,
  std::ostream &out,
  const unsigned int precision);

template void print_matrix_in_python (const PETScWrappers::SparseMatrix &mat,
  const std::string &mat_name,
  std::ostream &out,
  const unsigned int precision);

template void print_matrix_in_python (
  const PETScWrappers::MPI::SparseMatrix &mat,
  const std::string &mat_name,
  std::ostream &out,
  const unsigned int precision);

/*
 *
 */
void copy_to_Vector (Vector<double> &dst,
  Vec src)
{

  const double *aray_vec;
  VecGetArrayRead(src, &aray_vec);
  for (unsigned int i = 0; i < dst.size(); i++)
  {
    dst[i] = aray_vec[i];
  }
  VecRestoreArrayRead(src, &aray_vec);

  dst.compress(VectorOperation::insert);
}

/*
 *
 */
void copy_to_Vector (LinearAlgebra::distributed::Vector<double> &dst,
  Vec src)
{

  const double *aray_vec;
  VecGetArrayRead(src, &aray_vec);
  IndexSet index_set(dst.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {
    dst[*it] = aray_vec[i];
  }
  VecRestoreArrayRead(src, &aray_vec);

  dst.compress(VectorOperation::insert);
  index_set.clear();
}

/*
 *
 */
void copy_to_BlockVector (BlockVector<double> &dst,
  Vec src)
{
  // Assert Dimensions
  const double *aray_vec;
  VecGetArrayRead(src, &aray_vec);
  IndexSet index_set(dst.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {
    dst[*it] = aray_vec[i];
  }
  VecRestoreArrayRead(src, &aray_vec);

  dst.compress(VectorOperation::insert);
  index_set.clear();
}

/**
 *
 */
void copy_to_Vector (PETScWrappers::MPI::Vector &dst,
  const LinearAlgebra::distributed::Vector<double> &src)
{
  unsigned int i = 0;
  IndexSet index_set(src.locally_owned_elements());
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {
    dst[*it] = src.local_element(i);
  }

  dst.compress(VectorOperation::insert);
  index_set.clear();
}
/**
 * @brief Copy a Petsc Vec to std::vector<PETScWrappers::MPI::BlockVector>
 * Be careful this function involves copy a vector.
 */
void copy_to_stdBlockVector (std::vector<PETScWrappers::MPI::BlockVector> &dst,
  Vec src)
{
  const double *aray_vec;
  VecGetArrayRead(src, &aray_vec);

  int i = 0;
  IndexSet index_set;
  for (unsigned int b = 0; b < dst.size(); b++)
    for (unsigned int c = 0; c < dst[b].n_blocks(); c++)
    {
      //const unsigned int block_size = dst[b].block(c).size();
      index_set = dst[b].block(c).locally_owned_elements();

      for (IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++, i++)
      {
        //std::cout << "num " << *it + c * block_size << std::endl;
        dst[b].block(c)[*it] = aray_vec[i];

      }
    }
  VecRestoreArrayRead(src, &aray_vec);

  for (unsigned int b = 0; b < dst.size(); b++)
    dst[b].compress(VectorOperation::insert);

  index_set.clear();

  return;
}

/**
 * @brief Copy a Petsc Vec to ETScWrappers::MPI::BlockVector.
 * Be careful this function involves copy a vector.
 */
void copy_to_BlockVector (PETScWrappers::MPI::BlockVector &dst,
  Vec src)

{
  const double *aray_vec;
  VecGetArrayRead(src, &aray_vec);
  IndexSet index_set(dst.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {
    dst[*it] = aray_vec[i];
  }
  VecRestoreArrayRead(src, &aray_vec);

  dst.compress(VectorOperation::insert);
  index_set.clear();

  return;
}

/**
 * @brief Copy a  BlockVector<double> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const BlockVector<double> &src)
{

  // Assert Dimensions
  int s;
  VecGetSize(dst, &s);
  AssertDimension(src.size(), (unsigned int) s);
  double *aray_vec;
  VecGetArray(dst, &aray_vec);
  IndexSet index_set(src.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {

    aray_vec[i] = src[*it];
  }
  VecRestoreArray(dst, &aray_vec);

  VecAssemblyBegin(dst);
  VecAssemblyEnd(dst);
  index_set.clear();

  return;
}

/**
 * @brief Copy a  BlockVector<double> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const Vector<double> &src)
{
  // Assert Dimensions
  int s;
  VecGetSize(dst, &s);
  AssertDimension(src.size(), (unsigned int) s);
  double *aray_vec;
  VecGetArray(dst, &aray_vec);
  for (int i = 0; i < s; i++)
  {
    aray_vec[i] = src[i];
  }
  VecRestoreArray(dst, &aray_vec);

  VecAssemblyBegin(dst);
  VecAssemblyEnd(dst);

  return;
}

/**
 * @brief Copy a  BlockVector<double> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const LinearAlgebra::distributed::Vector<double> &src)
{

  // Assert Dimensions
  int s;
  VecGetSize(dst, &s);
  AssertDimension(src.size(), (unsigned int) s);
  double *aray_vec;
  VecGetArray(dst, &aray_vec);
  IndexSet index_set(src.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {

    aray_vec[i] = src[*it];
  }
  VecRestoreArray(dst, &aray_vec);

  VecAssemblyBegin(dst);
  VecAssemblyEnd(dst);

  return;
}

/**
 * @brief Copy a  BlockVector<double> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (PETScWrappers::MPI::Vector dst,
  const LinearAlgebra::distributed::Vector<double> &src)
{

  // Assert Dimensions
  int s;
  VecGetSize(dst, &s);
  AssertDimension(src.size(), (unsigned int) s);
  IndexSet index_set(dst.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {

    dst[*it] = src[i];
  }

  dst.compress(VectorOperation::insert);

  return;
}

/**
 * @brief Copy a  PETScWrappers::MPI::BlockVector to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const PETScWrappers::MPI::BlockVector &src)
{

  // Assert Dimensions
  int s;
  VecGetSize(dst, &s);
  AssertDimension(src.size(), (unsigned int) s);

  double *aray_vec;
  VecGetArray(dst, &aray_vec);
  IndexSet index_set(src.locally_owned_elements());
  int i = 0;
  for (IndexSet::ElementIterator it = index_set.begin();
      it != index_set.end(); it++, i++)
  {
    aray_vec[i] = src[*it];
  }
  VecRestoreArray(dst, &aray_vec);

  VecAssemblyBegin(dst);
  VecAssemblyEnd(dst);

  return;
}

/**
 * @brief Copy a std::vector<PETScWrappers::MPI::BlockVector> to Petsc Vec.
 * Be careful this function involves copy a vector.
 */
void copy_to_Vec (Vec dst,
  const std::vector<PETScWrappers::MPI::BlockVector> &src)
{
  // Assert Dimensions
  int s;
  VecGetSize(dst, &s);
  AssertDimension(src.size(), (unsigned int) s);

  double *aray_vec;
  VecGetArray(dst, &aray_vec);
  IndexSet index_set;
  const unsigned int block_size_b = src[0].size();
  const unsigned int block_size_c = src[0].block(0).size();

  for (unsigned int b = 0; b < src.size(); b++)
    for (unsigned int c = 0; c < src[b].n_blocks(); c++)
    {
      index_set = src[b].block(c).locally_owned_elements();
      int i = 0;
      for ( IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++, i++)
      {
        aray_vec[i + block_size_b * b + c * block_size_c] = src[b].block(c)[*it];
      }
    }
  VecRestoreArray(dst, &aray_vec);
  VecAssemblyBegin(dst);
  VecAssemblyEnd(dst);

  return;
}

/**
 *
 */
void compute_max_eig (const LAPACKFullMatrix<double> &matrix,
  double &eig,
  std::vector<double> &eigenvector)
{
  AssertDimension(eigenvector.size(), matrix.m());
  AssertDimension(matrix.n(), matrix.m());

  std::vector<double> wr;
  std::vector<double> wi;
  std::vector<double> vr;
  std::vector<double> vl;
  std::vector<double> work;

  const int n = matrix.m();
  const char N = 'N';
  const char V = 'V';

  const bool right = true;
  const bool left = false;

  wr.resize(n);
  wi.resize(n);

  if (right)
    vr.resize(n * n);
  if (left)
    vl.resize(n * n);

  // Copy matrix values
  // In LAPACK and in LAPACKFullMatrix the values are stored in column order (as in FORTRAN)
  AlignedVector<double> values_;
  values_.resize(n * n, 0.0);
  for (int i = 0; i < n * n; i++)
    values_[i] = matrix(i % n, i / n);
  double *values = &values_[0];

  int info = 0;
  int lwork = 1;
  const char *const jobvr = (right) ? (&V) : (&N);
  const char *const jobvl = (left) ? (&V) : (&N);

  // The LAPACK routine xGEEV :
  lwork = -1;
  work.resize(1);

  geev(jobvl, jobvr, &n, values, &n, &wr[0], &wi[0], &vl[0], &n, &vr[0], &n,
    &work[0], &lwork, &info);
  // geev returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert(info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<int>(work[0] + 1);

  // resize workspace array
  work.resize((unsigned int) lwork);

  // Finally compute the eigenvalues.
  geev(jobvl, jobvr, &n, values, &n, &wr[0], &wi[0], &vl[0], &n, &vr[0], &n,
    &work[0], &lwork, &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert(info >=0, ExcInternalError());

  if (info != 0)
    std::cerr << "LAPACK error in geev" << std::endl;

  eig = 0.0;
  for (int i = 0; i < n; i++)
    if (wr[i] > eig)
    {
      eig = wr[i];
      for (int j = 0; j < n; j++)
        eigenvector[j] = vr[i * n + j];
    }
}

/**
 *
 */
void compute_eigs (const LAPACKFullMatrix<double> &matrix,
  std::vector<double> &wr,
  LAPACKFullMatrix<double> &eigenvectors,
  bool sort)
{

  AssertDimension(wr.size(), matrix.m());
  AssertDimension(matrix.n(), matrix.m());
  AssertDimension(eigenvectors.m(), matrix.m());

  std::vector<double> wi;
  std::vector<double> vr;
  std::vector<double> vl;
  std::vector<double> work;

  int n = matrix.m();
  const char N = 'N';
  const char V = 'V';

  bool right = true;
  bool left = false;

  wr.resize(n);
  wi.resize(n);

  if (right)
    vr.resize(n * n);
  if (left)
    vl.resize(n * n);

  // Copy matrix values
  // In LAPACK and in LAPACKFullMatrix the values are stored in column order (as in FORTRAN)
  AlignedVector<double> values_;
  values_.resize(n * n, 0.0);
  for (int i = 0; i < n * n; i++)
    values_[i] = matrix(i % n, i / n);
  double *values = &values_[0];

  int info = 0;
  int lwork = 1;
  const char *const jobvr = (right) ? (&V) : (&N);
  const char *const jobvl = (left) ? (&V) : (&N);

  // The LAPACK routine xGEEV :
  lwork = -1;
  work.resize(1);

  geev(jobvl, jobvr, &n, values, &n, &wr[0], &wi[0], &vl[0], &n, &vr[0], &n,
    &work[0], &lwork, &info);
  // geev returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert(info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<int>(work[0] + 1);

  // resize workspace array
  work.resize((unsigned int) lwork);

  // Finally compute the eigenvalues.
  geev(jobvl, jobvr, &n, values, &n, &wr[0], &wi[0], &vl[0], &n, &vr[0], &n,
    &work[0], &lwork, &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert(info >=0, ExcInternalError());

  if (info != 0)
    std::cerr << "LAPACK error in geev" << std::endl;

  //  for ( int j = 0; j < n; ++j)
  //      std::cout<<"real: "<<wr[j]<<"img: "<<wi[j]<<std::endl;

  // Sort the eigenvalues by largest magnitude

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      eigenvectors(i, j) = vr[j * n + i];

  if (sort == true)
  {
    double aux;
    std::vector<double> vecaux;
    vecaux.resize(n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n - i - 1; j++)
        if (std::abs(wr[j]) < std::abs(wr[j + 1]))
        {
          aux = wr[j];
          wr[j] = wr[j + 1];
          wr[j + 1] = aux;
          for (int l = 0; l < n; l++)
          {
            vecaux[l] = eigenvectors(l, j);
            eigenvectors(l, j) = eigenvectors(l, j + 1);
            eigenvectors(l, j + 1) = vecaux[l];
          }
        }
  }

  // Checking
  LAPACKFullMatrix<double> check_mat;
  check_mat.reinit(n);

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
    {
      for (int k = 0; k < n; ++k)
      {
        check_mat(i, j) += matrix(i, k) * eigenvectors(k, j);
      }
      check_mat(i, j) -= eigenvectors(i, j) * wr[j];
      Assert(check_mat(i,j)<1e-14,
        ExcMessage("Error solving with LAPACK."));
    }

}

/**
 *
 */
template <int dim>
  void join_vectors (const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const PETScWrappers::MPI::Vector &phi_1,
    const PETScWrappers::MPI::Vector &phi_2,
    PETScWrappers::MPI::Vector &phi_joint)
  {

    unsigned int i1, i2;
    const unsigned int dofs_per_cell_joint =
                                             dof_handler_joint.get_fe().dofs_per_cell;
    const unsigned int dofs_per_cell_sep =
                                           dof_handler_separated.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices_joint(
      dofs_per_cell_joint);
    std::vector<types::global_dof_index> local_dof_indices_sep(
      dofs_per_cell_sep);

    // Pre-built list to which component a given dof on a cell should go.
    std::vector<unsigned int> local_component_list(dofs_per_cell_joint);
    for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      local_component_list[i] =
                                dof_handler_joint.get_fe().system_to_component_index(
                                  i).first;

    // Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell_sep =
        dof_handler_separated.begin_active();
    // endc_out = dof_handler_separated.end();
    typename DoFHandler<dim>::active_cell_iterator cell_joint =
        dof_handler_joint.begin_active(), endc = dof_handler_joint.end();
    for (; cell_joint != endc; ++cell_joint, ++cell_sep)
    {
      // On each cell: get dof indices and insert them into the global
      // list using their component;
      cell_joint->get_dof_indices(local_dof_indices_joint);
      cell_sep->get_dof_indices(local_dof_indices_sep);
      i1 = 0;
      i2 = 0;
      for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      {
        if (local_component_list[i] == 0)
        {
          // std::cout << "i " << i << " i1 " << i1 << "  global "
          //   << local_dof_indices_joint[i] << std::endl;
          phi_joint[local_dof_indices_joint[i]] = phi_1[local_dof_indices_sep[i1]];
          i1++;

        }
        else if (local_component_list[i] == 1)
        {
          // std::cout << "i " << i << " i2 " << i2 << "  global "
          //   << local_dof_indices_joint[i] << std::endl;
          phi_joint[local_dof_indices_joint[i]] = phi_2[local_dof_indices_sep[i2]];
          i2++;
        }
      }
    }
  }

template void join_vectors<1> (
  const DoFHandler<1> &dof_handler_separated,
  const DoFHandler<1> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_1,
  const PETScWrappers::MPI::Vector &phi_2,
  PETScWrappers::MPI::Vector &phi_joint);
template void join_vectors<2> (
  const DoFHandler<2> &dof_handler_separated,
  const DoFHandler<2> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_1,
  const PETScWrappers::MPI::Vector &phi_2,
  PETScWrappers::MPI::Vector &phi_joint);
template void join_vectors<3> (
  const DoFHandler<3> &dof_handler_separated,
  const DoFHandler<3> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_1,
  const PETScWrappers::MPI::Vector &phi_2,
  PETScWrappers::MPI::Vector &phi_joint);

/**
 *
 */
template <int dim>
  void join_vectors (const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const std::vector<PETScWrappers::MPI::Vector> &phi_sep,
    PETScWrappers::MPI::Vector &phi_joint)
  {

    unsigned int n_components = dof_handler_joint.get_fe().n_components();
    std::vector<unsigned int> it(n_components);

    const unsigned int dofs_per_cell_joint = dof_handler_joint.get_fe().dofs_per_cell;
    const unsigned int dofs_per_cell_sep =
                                           dof_handler_separated.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices_joint(dofs_per_cell_joint);
    std::vector<types::global_dof_index> local_dof_indices_sep(dofs_per_cell_sep);

    // Pre-built list to which component a given dof on a cell should go.
    std::vector<unsigned int> local_component_list(dofs_per_cell_joint);
    for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      local_component_list[i] =
                                dof_handler_joint.get_fe().system_to_component_index(
                                  i).first;

    // Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell_sep =
        dof_handler_separated.begin_active();
    typename DoFHandler<dim>::active_cell_iterator cell_joint =
        dof_handler_joint.begin_active(), endc = dof_handler_joint.end();
    for (; cell_joint != endc; ++cell_joint, ++cell_sep)
    {
      // On each cell: get dof indices and insert them into the global
      // list using their component;
      cell_joint->get_dof_indices(local_dof_indices_joint);
      cell_sep->get_dof_indices(local_dof_indices_sep);

      for (unsigned int nc = 0; nc < n_components; nc++)
        it[nc] = 0;

      for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      {
        for (unsigned int nc = 0; nc < n_components; nc++)
        {
          if (local_component_list[i] == nc)
          {
            // std::cout << "i " << i << " i1 " << i1 << "  global "
            //   << local_dof_indices_joint[i] << std::endl;
            phi_joint[local_dof_indices_joint[i]] =
                phi_sep[nc][local_dof_indices_sep[it[nc]]];
            it[nc]++;

          }
        }
      }
    }

    phi_joint.compress(VectorOperation::insert);
  }

template void join_vectors<1> (const DoFHandler<1> &dof_handler_separated,
  const DoFHandler<1> &dof_handler_joint,
  const std::vector<PETScWrappers::MPI::Vector> &phi_sep,
  PETScWrappers::MPI::Vector &phi_joint);
template void join_vectors<2> (const DoFHandler<2> &dof_handler_separated,
  const DoFHandler<2> &dof_handler_joint,
  const std::vector<PETScWrappers::MPI::Vector> &phi_sep,
  PETScWrappers::MPI::Vector &phi_joint);
template void join_vectors<3> (const DoFHandler<3> &dof_handler_separated,
  const DoFHandler<3> &dof_handler_joint,
  const std::vector<PETScWrappers::MPI::Vector> &phi_sep,
  PETScWrappers::MPI::Vector &phi_joint);

/**
 *
 */
template <int dim>
  void separate_vectors (const DoFHandler<dim> &dof_handler_separated,
    const DoFHandler<dim> &dof_handler_joint,
    const PETScWrappers::MPI::Vector &phi_sol,
    PETScWrappers::MPI::Vector &phi_1,
    PETScWrappers::MPI::Vector &phi_2)
  {

    unsigned int i1, i2;
    const unsigned int dofs_per_cell_joint =
                                             dof_handler_joint.get_fe().dofs_per_cell;
    const unsigned int dofs_per_cell_sep =
                                           dof_handler_separated.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices_joint(
      dofs_per_cell_joint);
    std::vector<types::global_dof_index> local_dof_indices_sep(
      dofs_per_cell_sep);

    // Pre-built list to which component a given dof on a cell should go.
    std::vector<unsigned int> local_component_list(dofs_per_cell_joint);
    for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      local_component_list[i] =
                                dof_handler_joint.get_fe().system_to_component_index(
                                  i).first;

    // Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell_sep =
        dof_handler_separated.begin_active();
    //typename DoFHandler<dim>::active_cell_iterator endc_out = dof_handler_separated.end();
    typename DoFHandler<dim>::active_cell_iterator cell_joint =
        dof_handler_joint.begin_active(), endc = dof_handler_joint.end();
    for (; cell_joint != endc; ++cell_joint, ++cell_sep)
    {
      // On each cell: get dof indices and insert them into the global
      // list using their component;
      cell_joint->get_dof_indices(local_dof_indices_joint);
      cell_sep->get_dof_indices(local_dof_indices_sep);
      i1 = 0;
      i2 = 0;
      for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      {
        if (local_component_list[i] == 0)
        {
          // std::cout << "i " << i << " i1 " << i1 << "  global "
          //   << local_dof_indices_joint[i] << std::endl;
          phi_1[local_dof_indices_sep[i1]] =
                                             phi_sol[local_dof_indices_joint[i]];
          i1++;
        }
        else if (local_component_list[i] == 1)
        {
          // std::cout << "i " << i << " i2 " << i2 << "  global "
          //   << local_dof_indices_joint[i] << std::endl;
          phi_2[local_dof_indices_sep[i2]] =
                                             phi_sol[local_dof_indices_joint[i]];
          i2++;
        }
      }
    }
  }

template void separate_vectors<1> (const DoFHandler<1> &dof_handler_separated,
  const DoFHandler<1> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_sol,
  PETScWrappers::MPI::Vector &phi_1,
  PETScWrappers::MPI::Vector &phi_2);
template void separate_vectors<2> (const DoFHandler<2> &dof_handler_separated,
  const DoFHandler<2> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_sol,
  PETScWrappers::MPI::Vector &phi_1,
  PETScWrappers::MPI::Vector &phi_2);
template void separate_vectors<3> (const DoFHandler<3> &dof_handler_separated,
  const DoFHandler<3> &dof_handler_joint,
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
    std::vector<PETScWrappers::MPI::Vector> &phi_sep)
  {

    unsigned int n_components = dof_handler_joint.get_fe().n_components();
    std::vector<unsigned int> it(n_components);

    const unsigned int dofs_per_cell_joint = dof_handler_joint.get_fe().dofs_per_cell;
    const unsigned int dofs_per_cell_sep =
                                           dof_handler_separated.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices_joint(dofs_per_cell_joint);
    std::vector<types::global_dof_index> local_dof_indices_sep(dofs_per_cell_sep);

    // Pre-built list to which component a given dof on a cell should go.
    std::vector<unsigned int> local_component_list(dofs_per_cell_joint);
    for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      local_component_list[i] =
                                dof_handler_joint.get_fe().system_to_component_index(
                                  i).first;

    // Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell_sep =
        dof_handler_separated.begin_active();
    typename DoFHandler<dim>::active_cell_iterator cell_joint =
        dof_handler_joint.begin_active(), endc = dof_handler_joint.end();
    for (; cell_joint != endc; ++cell_joint, ++cell_sep)
    {
      // On each cell: get dof indices and insert them into the global
      // list using their component;
      cell_joint->get_dof_indices(local_dof_indices_joint);
      cell_sep->get_dof_indices(local_dof_indices_sep);

      for (unsigned int nc = 0; nc < n_components; nc++)
        it[nc] = 0;

      for (unsigned int i = 0; i < dofs_per_cell_joint; ++i)
      {
        for (unsigned int nc = 0; nc < n_components; nc++)
        {
          if (local_component_list[i] == nc)
          {
            phi_sep[nc][local_dof_indices_sep[it[nc]]] =
                phi_joint[local_dof_indices_joint[i]];
            it[nc]++;
          }
        }
      }

    }

    for (unsigned int nc = 0; nc < n_components; nc++)
      phi_sep[nc].compress(VectorOperation::insert);
  }

template void separate_vectors<1> (const DoFHandler<1> &dof_handler_separated,
  const DoFHandler<1> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_joint,
  std::vector<PETScWrappers::MPI::Vector> &phi_sep);
template void separate_vectors<2> (const DoFHandler<2> &dof_handler_separated,
  const DoFHandler<2> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_joint,
  std::vector<PETScWrappers::MPI::Vector> &phi_sep);
template void separate_vectors<3> (const DoFHandler<3> &dof_handler_separated,
  const DoFHandler<3> &dof_handler_joint,
  const PETScWrappers::MPI::Vector &phi_joint,
  std::vector<PETScWrappers::MPI::Vector> &phi_sep);

/**
 * Count the real number of non_zero_elements
 */
template <class MATRIX>
  unsigned int n_real_nonzero_elements (const MATRIX &mat,
    double tol)
  {

    unsigned int nnz = 0;
    for (auto it = mat.begin(); it != mat.end(); it++)
    {
      if (std::abs(it->value()) > tol)
        nnz++;
    }
    return nnz;

  }

template
unsigned int n_real_nonzero_elements (const PETScWrappers::SparseMatrix &mat,
  double tol);

template <int dim>
  RefinementCase<dim> cut_axis (const unsigned int)
  {
    Assert(false, ExcInternalError());
    return static_cast<unsigned char>(-1);
  }

template <>
  RefinementCase<1> cut_axis<1> (const unsigned int i)
  {

    Assert(i < 2, ExcIndexRange(i, 0, 2));

    static const RefinementCase<1> options[2] =
          {
            RefinementPossibilities<1>::no_refinement,
            RefinementPossibilities<1>::cut_x };
    return options[i];
  }

template <>
  RefinementCase<2> cut_axis<2> (const unsigned int i)
  {

    Assert(i < 4, ExcIndexRange(i, 0, 4));

    static const RefinementCase<2> options[4] =
          {
            RefinementPossibilities<2>::no_refinement,
            RefinementPossibilities<2>::cut_x,
            RefinementPossibilities<2>::cut_y,
            RefinementPossibilities<2>::cut_xy };
    return options[i];
  }

template <>
  RefinementCase<3> cut_axis<3> (const unsigned int i)
  {
    Assert(i < 8, ExcIndexRange(i, 0, 8));

    static const RefinementCase<3> options[8] =
          {
            RefinementPossibilities<3>::no_refinement,
            RefinementPossibilities<3>::cut_x,
            RefinementPossibilities<3>::cut_y,
            RefinementPossibilities<3>::cut_xy,
            RefinementPossibilities<3>::cut_z,
            RefinementPossibilities<3>::cut_xz,
            RefinementPossibilities<3>::cut_yz,
            RefinementPossibilities<3>::cut_xyz };
    return options[i];
  }

//
// interpolate_operator() only for 1 level of coarsening
//
template <int dim>
  void interpolate_operator (const DoFHandler<dim> &dof_handler,
    const std::vector<std::vector<unsigned int> > &cell_map,
    const PETScWrappers::MPI::Vector &vec_in,
    PETScWrappers::MPI::Vector &vec_out)
  {
    unsigned int dofs_per_cell = cell_map[0].size();
    Vector<double> local_values(dofs_per_cell);
    for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(0);
        cell != dof_handler.end(0); ++cell)
    {
      for (unsigned int i = 0; i < dofs_per_cell; i++)
        local_values(i) = vec_in(cell_map[cell->index()][i]);
      cell->set_dof_values_by_interpolation(local_values, vec_out);
    }
    VecAssemblyBegin(vec_out);
    VecAssemblyEnd(vec_out);
  }

template void interpolate_operator<1> (const DoFHandler<1> &dof_handler,
  const std::vector<std::vector<unsigned int> > &cell_map,
  const PETScWrappers::MPI::Vector &vec_in,
  PETScWrappers::MPI::Vector &vec_out);
template void interpolate_operator<2> (const DoFHandler<2> &dof_handler,
  const std::vector<std::vector<unsigned int> > &cell_map,
  const PETScWrappers::MPI::Vector &vec_in,
  PETScWrappers::MPI::Vector &vec_out);
template void interpolate_operator<3> (const DoFHandler<3> &dof_handler,
  const std::vector<std::vector<unsigned int> > &cell_map,
  const PETScWrappers::MPI::Vector &vec_in,
  PETScWrappers::MPI::Vector &vec_out);

//
// interpolate_operator() only for 1 level of coarsening
//
template <int dim>
  void restriction_operator (const DoFHandler<dim> &dof_handler,
    const std::vector<std::vector<unsigned int> > &cell_map,
    const PETScWrappers::MPI::Vector &vec_in,
    PETScWrappers::MPI::Vector &vec_out)
  {

    unsigned int dofs_per_cell = cell_map[0].size();
    Vector<double> local_values(dofs_per_cell);
    for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(0);
        cell != dof_handler.end(0); ++cell)
    {

      cell->get_interpolated_dof_values(vec_in, local_values);

      for (unsigned int i = 0; i < dofs_per_cell; i++)
      {
        vec_out(cell_map[cell->index()][i]) = local_values(i);
      }
    }

    VecAssemblyBegin(vec_out);
    VecAssemblyEnd(vec_out);
  }

template void restriction_operator<1> (const DoFHandler<1> &dof_handler,
  const std::vector<std::vector<unsigned int> > &cell_map,
  const PETScWrappers::MPI::Vector &vec_in,
  PETScWrappers::MPI::Vector &vec_out);
template void restriction_operator<2> (const DoFHandler<2> &dof_handler,
  const std::vector<std::vector<unsigned int> > &cell_map,
  const PETScWrappers::MPI::Vector &vec_in,
  PETScWrappers::MPI::Vector &vec_out);
template void restriction_operator<3> (const DoFHandler<3> &dof_handler,
  const std::vector<std::vector<unsigned int> > &cell_map,
  const PETScWrappers::MPI::Vector &vec_in,
  PETScWrappers::MPI::Vector &vec_out);

/**
 * @brief This function transfers the solution phi_old of fe coarse
 * with n_fe_degree=p_coarse into the dof_handler
 */
template <int dim>
  void solution_transfer_dof (const Triangulation<dim> &tria,
    unsigned int p_coarse,
    const DoFHandler<dim> &dof_handler_fine,
    const FE_Q<dim> &fe_fine,
    std::vector<PETScWrappers::MPI::BlockVector> &phi_coarse,
    std::vector<PETScWrappers::MPI::BlockVector> &phi_fine)
  {

    unsigned int n_dofs = dof_handler_fine.n_dofs();
    unsigned int n_eigenvalues = phi_coarse.size();
    unsigned int n_groups = phi_coarse[0].n_blocks();
    // Reinit
    phi_fine.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      phi_fine[eig].reinit(n_groups, MPI_COMM_WORLD, n_dofs, n_dofs);
      phi_fine[eig].compress(VectorOperation::insert);
    }

    // Coarse_fe
    FE_Q<dim> fe_coarse(QGaussLobatto<1>(p_coarse + 1));
    DoFHandler<dim> dof_handler_coarse(tria);

    //dof_handler_out.initialize(tria, fe_new);
    dof_handler_coarse.distribute_dofs(fe_coarse);
    FullMatrix<double> transfer_mat(fe_fine.dofs_per_cell,
      fe_coarse.dofs_per_cell);

    FETools::get_interpolation_matrix(fe_coarse, fe_fine, transfer_mat);

#ifdef DEBUG
    Vector<double> data1(dof_handler_coarse.n_dofs());
    Vector<double> data2(dof_handler_fine.n_dofs());
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      for (unsigned int g = 0; g < n_groups; g++)
      {
        copy_to_Vector(data1, phi_coarse[eig].block(g));

        VectorTools::interpolate(
            dof_handler_coarse,
            dof_handler_fine,
            transfer_mat,
            data1,
            data2);

        copy_to_Vec(phi_fine[eig].block(g), data2);
      }
    }
#else
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      for (unsigned int g = 0; g < n_groups; g++)
      {
        VectorTools::interpolate(dof_handler_coarse, dof_handler_fine,
          transfer_mat, phi_coarse[eig].block(g),
          phi_fine[eig].block(g));
      }
    }
#endif

  }

template void solution_transfer_dof<1> (const Triangulation<1> &tria,
  unsigned int p_coarse,
  const DoFHandler<1> &dof_handler,
  const FE_Q<1> &fe,
  std::vector<PETScWrappers::MPI::BlockVector> &phi_old,
  std::vector<PETScWrappers::MPI::BlockVector> &phi_new);

template void solution_transfer_dof<2> (const Triangulation<2> &tria,
  unsigned int p_coarse,
  const DoFHandler<2> &dof_handler,
  const FE_Q<2> &fe,
  std::vector<PETScWrappers::MPI::BlockVector> &phi_old,
  std::vector<PETScWrappers::MPI::BlockVector> &phi_new);

template void solution_transfer_dof<3> (const Triangulation<3> &tria,
  unsigned int p_coarse,
  const DoFHandler<3> &dof_handler,
  const FE_Q<3> &fe,
  std::vector<PETScWrappers::MPI::BlockVector> &phi_old,
  std::vector<PETScWrappers::MPI::BlockVector> &phi_new);

/**
 * @brief This function transfers the solution phi_old of fe coarse
 * with n_fe_degree=p_coarse into the dof_handler
 */
template <int dim>
  void solution_transfer_dof (const Triangulation<dim> &tria,
    unsigned int p_coarse,
    const DoFHandler<dim> &dof_handler_fine,
    const FE_Q<dim> &fe_fine,
    std::vector<BlockVector<double> > &phi_coarse,
    std::vector<BlockVector<double> > &phi_fine)
  {

    unsigned int n_dofs = dof_handler_fine.n_dofs();
    unsigned int n_eigenvalues = phi_coarse.size();
    unsigned int n_groups = phi_coarse[0].n_blocks();
    // Reinit
    phi_fine.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      phi_fine[eig].reinit(n_groups, n_dofs);

    // Coarse_fe
    FE_Q<dim> fe_coarse(QGaussLobatto<1>(p_coarse + 1));
    DoFHandler<dim> dof_handler_coarse(tria);

    //dof_handler_out.initialize(tria, fe_new);
    dof_handler_coarse.distribute_dofs(fe_coarse);
    FullMatrix<double> transfer_mat(fe_fine.dofs_per_cell,
      fe_coarse.dofs_per_cell);

    FETools::get_interpolation_matrix(fe_coarse, fe_fine, transfer_mat);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      for (unsigned int g = 0; g < n_groups; g++)
      {
        VectorTools::interpolate(dof_handler_coarse, dof_handler_fine,
          transfer_mat, phi_coarse[eig].block(g),
          phi_fine[eig].block(g));
      }
    }

  }

template void solution_transfer_dof<1> (const Triangulation<1> &tria,
  unsigned int p_coarse,
  const DoFHandler<1> &dof_handler,
  const FE_Q<1> &fe,
  std::vector<BlockVector<double> > &phi_old,
  std::vector<BlockVector<double> > &phi_new);

template void solution_transfer_dof<2> (const Triangulation<2> &tria,
  unsigned int p_coarse,
  const DoFHandler<2> &dof_handler,
  const FE_Q<2> &fe,
  std::vector<BlockVector<double> > &phi_old,
  std::vector<BlockVector<double> > &phi_new);

template void solution_transfer_dof<3> (const Triangulation<3> &tria,
  unsigned int p_coarse,
  const DoFHandler<3> &dof_handler,
  const FE_Q<3> &fe,
  std::vector<BlockVector<double> > &phi_old,
  std::vector<BlockVector<double> > &phi_new);

/**
 * @brief This function create the transfer matrix
 * to use the solution_interpolate_dof
 */
template <int dim>
  void set_transfer_matrix (const FE_Q<dim> &fe_in,
    const FE_Q<dim> &fe_out,
    FullMatrix<double> &transfer_mat)
  {

    Assert(transfer_mat.m()==fe_out.dofs_per_cell, ExcInternalError());
    Assert(transfer_mat.n()==fe_in.dofs_per_cell, ExcInternalError());

    if (fe_in.dofs_per_cell > fe_out.dofs_per_cell)
    {
      for (unsigned int nc = 0; nc < fe_out.dofs_per_cell; nc++)
        transfer_mat.set(nc, nc, 1.0);
    }

    else if (fe_in.dofs_per_cell < fe_out.dofs_per_cell)
    {
      fe_out.get_interpolation_matrix(fe_in, transfer_mat);
    }

  }

template void set_transfer_matrix<1> (const FE_Q<1> &fe_in,
  const FE_Q<1> &fe_out,
  FullMatrix<double> &transfer_mat);

template void set_transfer_matrix<2> (const FE_Q<2> &fe_in,
  const FE_Q<2> &fe_out,
  FullMatrix<double> &transfer_mat);

template void set_transfer_matrix<3> (const FE_Q<3> &fe_in,
  const FE_Q<3> &fe_out,
  FullMatrix<double> &transfer_mat);

/**
 * @brief This function create the transfer matrix
 * to use the solution_interpolate_dof
 */
template <int dim>
  void set_transfer_matrix (const FiniteElement<dim> &fe_in,
    const FiniteElement<dim> &fe_out,
    FullMatrix<double> &transfer_mat)
  {

    Assert(transfer_mat.m()==fe_out.dofs_per_cell, ExcInternalError());
    Assert(transfer_mat.n()==fe_in.dofs_per_cell, ExcInternalError());

//    if (fe_in.dofs_per_cell > fe_out.dofs_per_cell)
//    {
//      fe_out.get_interpolation_matrix(fe_in, transfer_mat);
//    }
//
//    else if (fe_in.dofs_per_cell < fe_out.dofs_per_cell)
//    {
//      fe_out.get_interpolation_matrix(fe_in, transfer_mat);
//    }
    fe_out.get_interpolation_matrix(fe_in, transfer_mat);

//    transfer_mat.print(std::cout,10,5);

  }

template void set_transfer_matrix<1> (const FiniteElement<1> &fe_in,
  const FiniteElement<1> &fe_out,
  FullMatrix<double> &transfer_mat);

template void set_transfer_matrix<2> (const FiniteElement<2> &fe_in,
  const FiniteElement<2> &fe_out,
  FullMatrix<double> &transfer_mat);

template void set_transfer_matrix<3> (const FiniteElement<3> &fe_in,
  const FiniteElement<3> &fe_out,
  FullMatrix<double> &transfer_mat);

/**
 * @brief This function transfers the solution in
 * into the solution out. It is needed the dofhandler_in,
 * the dof_handler_out and the transfer mat
 */
template <int dim>
  void solution_interpolate_dof (const DoFHandler<dim> &dof_handler_in,
    const DoFHandler<dim> &dof_handler_out,
    const FullMatrix<double> &transfer_mat,
    PETScWrappers::MPI::Vector &sol_in,
    PETScWrappers::MPI::Vector &sol_out)
  {

#ifdef DEBUG
    Vector<double> data1(dof_handler_in.n_dofs());
    Vector<double> data2(dof_handler_out.n_dofs());
    copy_to_Vector(data1, sol_in);

    VectorTools::interpolate(
        dof_handler_in,
        dof_handler_out,
        transfer_mat,
        data1,
        data2);

    copy_to_Vec(sol_out, data2);

#else

    VectorTools::interpolate(dof_handler_in, dof_handler_out, transfer_mat,
      sol_in, sol_out);

#endif

  }

template void solution_interpolate_dof<1> (const DoFHandler<1> &dof_handler_in,
  const DoFHandler<1> &dof_handler_out,
  const FullMatrix<double> &transfer_mat,
  PETScWrappers::MPI::Vector &sol_in,
  PETScWrappers::MPI::Vector &sol_out);

template void solution_interpolate_dof<2> (const DoFHandler<2> &dof_handler_in,
  const DoFHandler<2> &dof_handler_out,
  const FullMatrix<double> &transfer_mat,
  PETScWrappers::MPI::Vector &sol_in,
  PETScWrappers::MPI::Vector &sol_out);

template void solution_interpolate_dof<3> (const DoFHandler<3> &dof_handler_in,
  const DoFHandler<3> &dof_handler_out,
  const FullMatrix<double> &transfer_mat,
  PETScWrappers::MPI::Vector &sol_in,
  PETScWrappers::MPI::Vector &sol_out);

/**
 * @brief This function transfers the solution in
 * into the solution out. It is needed the dofhandler_in,
 * the dof_handler_out and the transfer mat
 */
template <int dim>
  void solution_interpolate_dof (const DoFHandler<dim> &dof_handler_in,
    const DoFHandler<dim> &dof_handler_out,
    const FullMatrix<double> &transfer_mat,
    PETScWrappers::MPI::BlockVector &sol_in,
    PETScWrappers::MPI::BlockVector &sol_out)
  {

    unsigned int nblocks = sol_in.n_blocks();
//
//#ifdef DEBUG

//	LinearAlgebra::distributed::Vector<double> data_in(dof_handler_in.n_dofs());
//	LinearAlgebra::distributed::Vector<double> data_out(dof_handler_out.n_dofs());
//
//    for (unsigned int g = 0; g < nblocks; g++)
//    {
//
//      copy_to_Vector(data_in, sol_in.block(g));
//      VectorTools::interpolate(
//          dof_handler_in,
//          dof_handler_out,
//          transfer_mat,
//		  data_in,
//		  data_out);
//
//      copy_to_Vec(sol_out.block(g), data_out);
//    }

//#else
//    LinearAlgebra::distributed::Vector<double> data_in(dof_handler_in.n_dofs());
//    LinearAlgebra::distributed::Vector<double> data_out(dof_handler_out.n_dofs());
//
//    for (unsigned int g = 0; g < nblocks; g++){
//    	copy_to_Vector(data_in, sol_in.block(g));
////    	data_in.compress(VectorOperation::insert);
//      VectorTools::interpolate(dof_handler_in, dof_handler_out, transfer_mat,
//    		  data_in, data_out);
//      copy_to_Vec(sol_out.block(g), data_out);
////      sol_out.block(g).compress(VectorOperation::insert);
//    }

    //FIXME FOR PARALLEL
    AssertRelease(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) < 2,
      "Fixme for parallel computation");
    sol_in.compress(VectorOperation::insert);
    for (unsigned int g = 0; g < nblocks; g++)
      VectorTools::interpolate(dof_handler_in, dof_handler_out, transfer_mat,
        sol_in.block(g), sol_out.block(g));

    sol_out.compress(VectorOperation::insert);

//#endif

  }

template void solution_interpolate_dof<1> (const DoFHandler<1> &dof_handler_in,
  const DoFHandler<1> &dof_handler_out,
  const FullMatrix<double> &transfer_mat,
  PETScWrappers::MPI::BlockVector &sol_in,
  PETScWrappers::MPI::BlockVector &sol_out);

template void solution_interpolate_dof<2> (const DoFHandler<2> &dof_handler_in,
  const DoFHandler<2> &dof_handler_out,
  const FullMatrix<double> &transfer_mat,
  PETScWrappers::MPI::BlockVector &sol_in,
  PETScWrappers::MPI::BlockVector &sol_out);

template void solution_interpolate_dof<3> (const DoFHandler<3> &dof_handler_in,
  const DoFHandler<3> &dof_handler_out,
  const FullMatrix<double> &transfer_mat,
  PETScWrappers::MPI::BlockVector &sol_in,
  PETScWrappers::MPI::BlockVector &sol_out);

/**
 * @brief Return the position of a child in a refined cell.
 * @return The position vector in (x, y, z) order.
 */
const std::vector<unsigned int> child_pos (const unsigned int child,
  const unsigned int n_levels,
  const unsigned int dim)
{
  std::vector<unsigned int> pos(dim);
  if (dim == 1)
    pos[0] = child;
  else if (dim == 2)
  {
    std::vector<unsigned int> level_pos(n_levels);
    unsigned int rem = child;
    for (unsigned int level = 0; level < n_levels; level++)
    {
      level_pos[level] = rem / pow(exp2(dim), n_levels - level - 1);
      rem = rem % static_cast<unsigned int>(pow(exp2(dim), n_levels - level - 1));
    }

    pos[0] = 0;
    pos[1] = 0;
    for (unsigned int level = 0; level < n_levels; level++)
    {
      pos[0] += (level_pos[level] % 2) * exp2(n_levels - level - 1);
      pos[1] += (level_pos[level] / 2) * exp2(n_levels - level - 1);
    }

  }
  else if (dim == 3)
  {
    std::vector<unsigned int> level_pos(n_levels);
    unsigned int rem = child;
    for (unsigned int level = 0; level < n_levels; level++)
    {
      level_pos[level] = rem / pow(exp2(dim), n_levels - level - 1);
      rem = rem % static_cast<unsigned int>(pow(exp2(dim), n_levels - level - 1));
    }

    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    for (unsigned int level = 0; level < n_levels; level++)
    {
      pos[0] += (level_pos[level] % 2) * exp2(n_levels - level - 1);
      pos[1] += ((level_pos[level] % 4) / 2) * exp2(n_levels - level - 1);
      pos[2] += (level_pos[level] / 4) * exp2(n_levels - level - 1);
    }
  }

  return pos;
}

/**
 * @brief
 * @return If the child is adjacent to that cell
 */
bool child_at_face (const unsigned int face,
  const unsigned int child,
  const unsigned int n_children,
  const unsigned int dim)
{
  const unsigned int n_levels = log10(n_children) / log10(exp2(dim));
  const unsigned int max_index_per_dim = exp2(n_levels) - 1; // Indeed max_child_index_per_dim
  std::vector<unsigned int> pos;

  pos = child_pos(child, n_levels, dim);

  if (dim == 1)
  {
    if (face == 0 and pos[0] == 0)
      return true;
    else if (face == 1 and pos[0] == max_index_per_dim)
      return true;
    else
      return false;
  }
  else if (dim == 2)
  {
    if (face == 0 and pos[0] == 0)
      return true;
    else if (face == 1 and pos[0] == max_index_per_dim)
      return true;
    else if (face == 2 and pos[1] == 0)
      return true;
    else if (face == 3 and pos[1] == max_index_per_dim)
      return true;
    else
      return false;
  }
  else if (dim == 3)
  {
    if (face == 0 and pos[0] == 0)
      return true;
    else if (face == 1 and pos[0] == max_index_per_dim)
      return true;
    else if (face == 2 and pos[1] == 0)
      return true;
    else if (face == 3 and pos[1] == max_index_per_dim)
      return true;
    else if (face == 4 and pos[2] == 0)
      return true;
    else if (face == 5 and pos[2] == max_index_per_dim)
      return true;
    else
      return false;
  }

  return false;
}
