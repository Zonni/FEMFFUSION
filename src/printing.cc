/**
 * @file   printing.cc
 * @brief  Printing functions.
 */

#include "../include/printing.h"

/**
 *
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<double> &matrix,
  const std::vector<unsigned int> &n_cells_per_dim,
  std::ofstream &out,
  const Materials &materials,
  const std::string &introduction,
  bool fill_with_zeros,
  unsigned int precision)
{
  const double tol = pow(10, -(int) precision);
  // Set scientific notation
  out.setf(std::ios::scientific);
  out.precision(precision);

  out << introduction;

  unsigned int idx_mat = 0;
  for (unsigned int plane = 0; plane < n_cells_per_dim[2]; plane++)
  {
    if (dim == 3)
      out << "Plane " << plane + 1 << "\n";

    for (unsigned int j = 0; j < n_cells_per_dim[1]; j++)
    {
      for (unsigned int i = 0; i < n_cells_per_dim[0]; i++)
      {
        if (materials.exist(i, j, plane))
        {
          // Avoid negative zeros
          if ((matrix[idx_mat] < tol) & (matrix[idx_mat] > -tol))
          {
            out << " " << 0.0 << " ";
          }
          else
          {
            if (matrix[idx_mat] > 0.0)
              out << " ";
            out << matrix[idx_mat] << " ";
          }
          idx_mat++;
        }
        // Print spaces or zeros when there is no material?
        else
        {
          if (fill_with_zeros)
            out << " " << 0.0 << " ";

          else
          {
            out << "       ";
            for (unsigned int p = 0; p < precision; p++)
              out << " ";
            out << " ";
          }
        }
      }
      out << "\n";
    }
  }
  out << "\n\n";
}

/**
 *
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<double> &matrix,
  const std::vector<unsigned int> &n_cells_per_dim,
  const std::string &filename,
  const Materials &materials,
  const std::string &introduction,
  bool fill_with_zeros,
  unsigned int precision)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  print_cell_distribution_in_file(
    dim,
    matrix,
    n_cells_per_dim,
    out,
    materials,
    introduction,
    fill_with_zeros,
    precision);
  out.close();
}

/**
 *
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<std::complex<double> > &matrix,
  const std::vector<unsigned int> &n_cells_per_dim,
  std::ofstream &out,
  const Materials &materials,
  const std::string &introduction,
  bool fill_with_zeros,
  unsigned int precision)
{
  // Set scientific notation
  out.setf(std::ios::scientific);
  out.precision(precision);
  out << std::showpos; // show + sign before positive numbers
  out << introduction;

  unsigned int idx_mat = 0;
  for (unsigned int plane = 0; plane < n_cells_per_dim[2]; plane++)
  {
    if (dim >= 3)
      out << "Plane " << plane + 1 << "\n";

    for (unsigned int j = 0; j < n_cells_per_dim[1]; j++)
    {
      for (unsigned int i = 0; i < n_cells_per_dim[0]; i++)
      {
        if (materials.exist(i, j, plane))
        {
          out << matrix[idx_mat].real() << matrix[idx_mat].imag() << "j"
              << " ";

          idx_mat++;
        }
        // Print spaces or zeros when there is no material?
        else
        {
          if (fill_with_zeros)
            out << " " << 0.0 << " ";

          else
          {
            out << "       ";
            for (unsigned int p = 0; p < precision; p++)
              out << " ";
            out << " ";
          }
        }
      }
      out << "\n";
    }
  }
  out << "\n\n";
  out << std::noshowpos;
}

/**
 *
 */
void print_cell_distribution_in_file (
  const unsigned int dim,
  const std::vector<std::complex<double> > &matrix,
  const std::vector<unsigned int> &n_cells_per_dim,
  const std::string &filename,
  const Materials &materials,
  const std::string &introduction,
  bool fill_with_zeros,
  unsigned int precision)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  print_cell_distribution_in_file(
    dim,
    matrix,
    n_cells_per_dim,
    out,
    materials,
    introduction,
    fill_with_zeros,
    precision);
  out.close();
}

/**
 * Print a vector in a File with an introduction.
 */
void print_vector_in_file (const std::vector<double> &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  print_vector_in_file(vect, out, introduction, horizontal, precision);
  out.close();
}

/**
 * Print a vector in a File with an introduction.
 */
void print_vector_in_file (const std::vector<unsigned int> &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  print_vector_in_file(vect, out, introduction, horizontal, precision);
  out.close();
}

/**
 * Print a vector in a File with an introduction.
 */
void print_vector_in_file (const PETScWrappers::MPI::Vector &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  print_vector_in_file(vect, out, introduction, horizontal, precision);
  out.close();
}

/**
 * Print a std::vector<double> in a File with an introduction.
 */
void print_vector_in_file (const std::vector<double> &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal,
  int precision)
{
  out.setf(std::ios::scientific);
  out.precision(precision);
  out << introduction;
  if (horizontal == true)
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " ";
  }
  else
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " \n";
  }
  out << "\n \n";
}

/**
 * Print a std::vector<double> in a File with an introduction.
 */
void print_vector_in_file (const std::vector<unsigned int> &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal,
  int precision)
{
  out.setf(std::ios::scientific);
  out.precision(precision);
  out << introduction;
  if (horizontal == true)
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " ";
  }
  else
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " \n";
  }
  out << "\n \n";
}

/**
 * Print a std::vector<double> in a File with an introduction.
 */
void print_vector_in_file (const PETScWrappers::MPI::Vector &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal,
  int precision)
{
  out.setf(std::ios::scientific);
  out.precision(precision);
  out << introduction;
  if (horizontal == true)
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " ";
  }
  else
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " \n";
  }
  out << "\n";
}

/**
 * @brief Print a VPETScWrappers::MPI::BlockVector in a file with an introduction.
 */
void print_vector_in_file (const PETScWrappers::MPI::BlockVector &vect,
  std::ofstream &out,
  std::string introduction,
  bool horizontal,
  int precision)
{
  out.setf(std::ios::scientific);
  out.precision(precision);
  out << introduction;
  if (horizontal == true)
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " ";
  }
  else
  {
    for (unsigned int i = 0; i < vect.size(); i++)
      out << vect[i] << " \n";
  }
  out << "\n \n";
}


/**
 * Print a std::vector<double> in a File with an introduction.
 */
template <int dim>
  void print_vector_in_file (const std::vector<Point<dim> > &vect,
    std::string filename,
    std::string introduction,
    bool horizontal,
    int precision)
  {
    std::ofstream out(filename.c_str(), std::ios::app);
    out.setf(std::ios::scientific);
    out.precision(precision);
    out << introduction;
    if (horizontal == true)
    {
      for (unsigned int d = 0; d < dim; ++d)
      {
        for (unsigned int i = 0; i < vect.size(); i++)
          out << vect[i][d] << " ";

        out << "\n";
      }
    }
    else
    {
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int i = 0; i < vect.size(); i++)
          out << vect[i][d] << " \n";
    }
    out << "\n \n";
    out.close();
  }

template void print_vector_in_file (const std::vector<Point<1> > &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision);
template void print_vector_in_file (const std::vector<Point<2> > &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision);
template void print_vector_in_file (const std::vector<Point<3> > &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision);

/**
 *
 */
void print_matrix_in_file (const std::vector<std::vector<double> > &mat,
  std::string filename,
  std::string introduction,
  int precision)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  print_matrix_in_file(mat, out, introduction, precision);
  out.close();
}

/**
 *
 */
void print_matrix_in_file (const std::vector<std::vector<double> > &mat,
  std::ofstream &out,
  std::string introduction,
  int precision)
{
  out.setf(std::ios::scientific);
  out.precision(precision);
  out << introduction;
  const double tol = pow(10, -(int) precision);

  for (unsigned int i = 0; i < mat.size(); i++)
  {
    for (unsigned int j = 0; j < mat[i].size(); j++)
    {
      // Avoid negative zeros
      if ((mat[i][j] < tol) and (mat[i][j] > -tol))
      {
        out << " " << 0.0 << " ";
      }
      else
      {
        if (mat[i][j] > 0.0)
          out << " ";
        out << mat[i][j] << " ";
      }
    }
    out << "\n";
  }
  out << "\n \n";
}

/**
 * @brief Print a number in a file with an introduction.
 */
template <class num>
  void print_in_file (num numb,
    std::string filename,
    std::string introduction,
    unsigned int precision)
  {
    std::ofstream out(filename.c_str(), std::ios::app);
    out.precision(precision);
    out << introduction << numb << " \n";
    out.close();
  }

template void print_in_file (double numb,
  std::string filename,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (unsigned int numb,
  std::string filename,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (int numb,
  std::string filename,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (std::string numb,
  std::string filename,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (const char *numb,
  std::string filename,
  std::string introduction = "\n",
  unsigned int precision = 8);

/**
 *
 */
template <class num>
  void print_in_file (num numb,
    std::ofstream &out,
    std::string introduction,
    unsigned int precision)
  // Print a std::vector<double> in a ofstream with an introduction.
  {
    out.precision(precision);
    out << introduction << numb << " \n";
  }

template void print_in_file (double numb,
  std::ofstream &out,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (unsigned int numb,
  std::ofstream &out,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (int numb,
  std::ofstream &out,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (std::string numb,
  std::ofstream &out,
  std::string introduction = "\n",
  unsigned int precision = 8);
template void print_in_file (const char *numb,
  std::ofstream &out,
  std::string introduction = "\n",
  unsigned int precision = 8);

/**
 *
 */
template <int dim>
  void print_points_in_file (const std::vector<Point<dim> > &vect,
    std::string filename,
    std::string introduction,
    bool horizontal,
    int precision)
  {

    std::ofstream out(filename.c_str(), std::ios::app);
    out.precision(precision);
    out << introduction;
    std::vector<std::string> introDim;
    introDim.push_back("x :");
    introDim.push_back("y :");
    introDim.push_back("z :");

    for (int d = 0; d < dim; d++)
    {
      out << introDim[d] << " \n";
      if (horizontal == true)
      {
        for (unsigned int i = 0; i < vect.size(); i++)
          out << vect[i][d] << " ";
      }
      else
      {
        for (unsigned int i = 0; i < vect.size(); i++)
          out << vect[i][d] << " \n";
      }

      out << "\n \n";
    }
    out.close();
  }

template void print_points_in_file<1> (const std::vector<Point<1> > &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision);
template void print_points_in_file<2> (const std::vector<Point<2> > &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision);
template void print_points_in_file<3> (const std::vector<Point<3> > &vect,
  std::string filename,
  std::string introduction,
  bool horizontal,
  int precision);

