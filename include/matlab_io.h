/**
 * @file   matlab_io.h
 * @brief  Implementation to read and write variables in a .mat from matlab using matio.
 */

#ifndef MATLAB_IO_H_
#define MATLAB_IO_H_

#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector.h>
#include "femffusion.h"

#ifdef MATIO
#include  <matio.h>

/**
 *
 */
void load_matlab_matrix (
                         mat_t *matlab_fp,
                         const std::string var_name,
                         std::vector<double> &data);

/**
 *
 */
void load_matlab_xs (
                     mat_t *matlab_fp,
                     const std::string xs_name,
                     const unsigned int size,
                     std::vector<double> &xs);

/**
 *
 */
void load_matlab_dxs (
                      mat_t *matfp,
                      const std::string xs_name,
                      const unsigned int size,
                      std::vector<std::complex<double> > &xs);

/**
 *
 */
void load_matlab_num (
                      mat_t *matfp,
                      const std::string var_name,
                      double &num);

/**
 *
 */
void load_matlab_dimension (
                            mat_t *matfp,
                            const std::string xs_name,
                            std::vector<unsigned int> &assem_per_dim);

/**
 *
 */
void write_matlab_matrix (
                          mat_t *&matfp,
                          const std::string &var_name,
                          const std::vector<unsigned int> &dimensions,
                          const Vector<double> &data);

/**
 *
 */
void write_matlab_matrix (
                          mat_t *&matfp,
                          const std::string &var_name,
                          const std::vector<unsigned int> &dimensions,
                          const Vector<double> &data_real,
                          const Vector<double> &data_imag);

/**
 *
 */
void write_matlab_number (
                          mat_t *&matfp,
                          const std::string &var_name,
                          const std::complex<double> &num);

/**
 *
 */
void write_matlab_number (
                          mat_t *&matfp,
                          const std::string &var_name,
                          const double &num);

#endif /* MATLAB_IO_H_ */
#endif /* MATIO */
