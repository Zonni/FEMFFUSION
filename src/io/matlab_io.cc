/**
 * @file   matlab_io.cc
 * @brief  Implementation to read and write variables in a .mat from matlab using matio
 */
#ifdef MATIO
#include <iostream>

#include "../include/femffusion.h"
#include "matlab_io.h"

using namespace dealii;


/**
 *
 */
void load_matlab_matrix (
  mat_t* matfp,
  const std::string var_name,
  std::vector<double>& data)
{
  matvar_t *matvar;
  matvar = Mat_VarRead(matfp, var_name.c_str());
  AssertRelease(NULL != matvar, "Variable " + var_name +
                                " not found, or error reading MAT file");
  //AssertRelease(!matvar->isComplex, "Complex XS " + var_name);
  //AssertRelease(matvar->rank == 3, var_name + " must have rank 3");

  double* data_array = static_cast<double*>(matvar->data);
  unsigned data_array_size = sizeof(data_array) / sizeof(double);

  data.insert(data.end(), &data_array[0], &data_array[data_array_size]);
  Mat_VarFree(matvar);
}

/**
 *
 */
void load_matlab_num (
  mat_t* matfp,
  const std::string var_name,
  double& num)
{
  matvar_t *matvar;
  matvar = Mat_VarRead(matfp, var_name.c_str());

  AssertRelease(NULL != matvar,
    "Variable " + var_name + " not found, or error reading MAT file");
  AssertRelease(!matvar->isComplex, "Complex XS " + var_name);
  AssertRelease(matvar->rank == 2, var_name + " must have rank 2 to be a number");

  double* data_array = static_cast<double*>(matvar->data);
  //unsigned data_array_size = sizeof(data_array) / sizeof(int);
  //AssertRelease(data_array_size == 1, var_name + " must be a number");
  num = data_array[0];

  Mat_VarFree(matvar);
}

/**
 *
 */
void load_matlab_xs (
  mat_t* matfp,
  const std::string xs_name,
  const unsigned int size,
  std::vector<double>& xs)
{
  Assert(xs.size() == 0, ExcMessage("The xs must be empty"));
  matvar_t *matvar;
  matvar = Mat_VarRead(matfp, xs_name.c_str());
  AssertRelease(NULL != matvar, "Variable " + xs_name +
                                " not found, or error reading MAT file");

  AssertRelease(!matvar->isComplex, "Complex XS " + xs_name);
  //AssertRelease(matvar->rank == 3, xs_name + " must have rank 3");

  double* dataArray = static_cast<double*>(matvar->data);

  unsigned int size_test = 1;
  for (int r = 0; r < matvar->rank; r++)
    size_test *= matvar->dims[r];
  AssertRelease(size_test == size, "Error in size of " + xs_name);

  xs.insert(xs.begin(), &dataArray[0], &dataArray[size]);
  Mat_VarFree(matvar);
}

/**
 *
 */
void load_matlab_dimension (
  mat_t* matfp,
  const std::string xs_name,
  std::vector<unsigned int>& assem_per_dim)
{
  matvar_t *matvar;
  matvar = Mat_VarReadInfo(matfp, xs_name.c_str());

  AssertRelease(NULL != matvar, "Error reading variable " + xs_name);
  // TODO
  //AssertRelease(!matvar->isComplex, "Complex number " + xs_name);
  //AssertRelease(matvar->rank == 3, xs_name + " must have rank 3");
  assem_per_dim.resize(3);
  for (int d = 0; d < 3; d++)
  {
    if (d < matvar->rank)
      assem_per_dim[d] = matvar->dims[d];
    else
      assem_per_dim[d] = 1;
  }

  Mat_VarFree(matvar);
}

/**
 *
 */
void load_matlab_dxs (
  mat_t* matfp,
  const std::string xs_name,
  const unsigned int size,
  std::vector<std::complex<double> >& dxs)
{
  Assert(dxs.size() == 0, ExcMessage("The xs must be empty"));
  matvar_t *matvar;
  matvar = Mat_VarRead(matfp, xs_name.c_str());

  if (NULL == matvar)
  {
    // Variable not found
    dxs.resize(size, std::complex<double>(0.0));
    return;
  }
  else
  {
    //    std::cout << " matvar->class_type " << matvar->class_type << std::endl;
    //    std::cout << "IS THE MATVAR COMPLEX??  " << matvar->isComplex << std::endl;
    AssertRelease(matvar->isComplex == false,
      "Do not work (reading) the complex perturbation values.");

    double* dataArray = static_cast<double*>(matvar->data);
    //    mat_complex_split_t z;
    //    double *ptr = new double[size];
    //    double *pti = new double[size];
    //Mat_VarReadData(matfp, matvar, &z, start, stride, edge);

    unsigned int size_test = 1;
    for (int r = 0; r < matvar->rank; r++)
      size_test *= matvar->dims[r];
    AssertRelease(size_test == size, "Error in size of " + xs_name);

    dxs.insert(dxs.begin(), &dataArray[0], &dataArray[size]);

    //    std::cout << xs_name << std::endl;
    //    for (unsigned int i = 0; i < size; i++)
    //      std::cout << i + 1 << " ->  " << dxs[i] << std::endl;
    Mat_VarFree(matvar);
  }
}

/**
 *
 */
void write_matlab_matrix (
  mat_t* &matfp,
  const std::string &var_name,
  const std::vector<unsigned int>& dimensions,
  const Vector<double>& data)
{
  matvar_t *matvar;
  int rank = dimensions.size();
  size_t* dims = new size_t[rank];
  unsigned int n_data = 1;

  for (unsigned int d = 0; d < dimensions.size(); d++)
  {
    n_data *= dimensions[d];
    dims[d] = static_cast<size_t>(dimensions[d]);
  }

  AssertRelease((data.size() == n_data),
    "Invalid size of real data in the output to Matlab");

  // Copy to array
  double *ptr = new double[n_data];
  double *pti = new double[n_data];
  for (unsigned int i = 0; i < n_data; i++)
  {
    ptr[i] = data(i);
    pti[i] = 0.0;
  }

  // Create Complex Structure
  struct mat_complex_split_t z;
  z.Re = ptr;
  z.Im = pti;

  matvar = Mat_VarCreate(var_name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, rank, dims, &z,
    MAT_F_COMPLEX);
  AssertRelease(matvar != NULL,
    "Error creating the variable" + var_name + " in .mat file.");

  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
  AssertRelease(matvar != NULL,
    "Error creating the variable " + var_name + " in .mat file.");
}

/**
 *
 */
void write_matlab_matrix (
  mat_t* &matfp,
  const std::string &var_name,
  const std::vector<unsigned int>& dimensions,
  const Vector<double>& data_real,
  const Vector<double>& data_imag)
{
  matvar_t *matvar;
  int rank = dimensions.size();
  size_t* dims = new size_t[rank];
  unsigned int n_data = 1;

  for (unsigned int d = 0; d < dimensions.size(); d++)
  {
    n_data *= dimensions[d];
    dims[d] = static_cast<size_t>(dimensions[d]);
  }

  AssertRelease((data_real.size() == n_data),
    "Invalid size of real data in the output to Matlab");
  AssertRelease((data_imag.size() == n_data),
    "Invalid size of imag data in the output to Matlab");

  // Copy to array
  double *ptr = new double[n_data];
  double *pti = new double[n_data];
  for (unsigned int i = 0; i < n_data; i++)
  {
    ptr[i] = data_real(i);
    pti[i] = data_imag(i);
  }

  // Create Complex Structure
  struct mat_complex_split_t z;
  z.Re = ptr;
  z.Im = pti;

  matvar = Mat_VarCreate(var_name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, rank, dims, &z,
    MAT_F_COMPLEX);
  AssertRelease(matvar != NULL,
    "Error creating the variable" + var_name + " in .mat file.");

  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
  AssertRelease(matvar != NULL,
    "Error creating the variable " + var_name + " in .mat file.");
}

/**
 *
 */
void write_matlab_number (
  mat_t* &matfp,
  const std::string &var_name,
  const std::complex<double>& num)
{
  matvar_t *matvar;
  const int rank = 2;
  size_t* dims = new size_t[rank];
  const unsigned int n_data = 1;
  dims[0] = 1;
  dims[1] = 1;

  // Prepare data
  double *ptr = new double[n_data];
  double *pti = new double[n_data];

  ptr[0] = num.real();
  pti[0] = num.imag();

  // Create Complex Structure
  struct mat_complex_split_t z;
  z.Re = ptr;
  z.Im = pti;

  matvar = Mat_VarCreate(var_name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, rank, dims, &z,
    MAT_F_COMPLEX);
  AssertRelease(matvar != NULL,
    "Error creating the variable" + var_name + " in .mat file.");

  // MAke a function of this
  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
  AssertRelease(matvar != NULL,
    "Error creating the variable " + var_name + " in .mat file.");
}

/**
 *
 */
void write_matlab_number (
  mat_t* &matfp,
  const std::string &var_name,
  const double& num)
{
  matvar_t *matvar;
  const int rank = 2;
  size_t* dims = new size_t[rank];
  dims[0] = 1;
  dims[1] = 1;
  const unsigned int n_data = 1;

  // Prepare data
  double *ptr = new double[n_data];
  double *pti = new double[n_data];

  ptr[0] = num;
  pti[0] = 0.0;

  // Create Complex Structure
  struct mat_complex_split_t z;
  z.Re = ptr;
  z.Im = pti;

  matvar = Mat_VarCreate(var_name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, rank, dims, &z,
    MAT_F_COMPLEX);
  AssertRelease(matvar != NULL,
    "Error creating the variable" + var_name + " in .mat file.");

  // Make a function of this
  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
  AssertRelease(matvar != NULL,
    "Error creating the variable " + var_name + " in .mat file.");
}
#endif
