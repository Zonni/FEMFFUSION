/**
 * @brief  Implementation of class Materials
 */

#ifndef COMPLEX_PERTURBATION_H_
#define COMPLEX_PERTURBATION_H_

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/table.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <iostream>
#include <fstream>
#include <sstream>

// This is for the .xml parser from boost
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <complex>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>

#include "matlab_io.h"
#include "materials.h"

using namespace dealii;
typedef std::complex<double> complex;

/**
 * @class ComplexPerturbation
 *
 * @brief Class to store the Complex pertubation (Fourier Transformed)
 * on the .
 *
 *
 */
class ComplexPerturbation
{
  public:

  /**
   * @brief Constructor
   */
  ComplexPerturbation (const Materials &materials,
    const unsigned int dim,
    const bool verbose = false);

  /**
   *
   */
  void reinit (const std::string &dxs_file,
    const std::string &xsec_type,
    const std::string &pert_type);

  /**
   * @brief Get the number of energy groups.
   * @return n_groups
   */
  unsigned int get_n_groups () const;

  /**
   * @brief Get the number of materials defined.
   * @return n_mats
   */
  unsigned int get_n_mats () const;

  /**
   *
   */
  void set_frequency (const double freq);

  /**
   *
   */
  double get_frequency () const;

  /**
   * @brief
   */
  unsigned int get_pertubation_face_id (const unsigned int mat_id,
    const unsigned int face,
    const unsigned int child,
    const unsigned int n_children) const;

  /**
   * @brief
   */
  unsigned int get_pertubation_face_id_hex (const unsigned int mat_id,
    const unsigned int face,
    const unsigned int quad_in_hex,
    const unsigned int child,
    const unsigned int n_children) const;

  /**
   * @brief
   */
  complex get_delta_sigma_f (const unsigned int group,
    const unsigned int mat) const;

  /**
   * @brief
   */
  complex get_delta_sigma_f (
    const unsigned int group_i,
    const unsigned int group_j,
    const unsigned int pert_mat,
    const unsigned int mat) const;

  /**
   * @brief
   */
  complex get_delta_sigma_t (const unsigned int group,
    const unsigned int mat) const;

  /**
   * @brief
   */
  complex get_delta_sigma_r (const unsigned int group,
    const unsigned int mat) const;

  /**
   *
   */
  complex get_delta_sigma_s (const unsigned int groupi,
    const unsigned int groupj,
    const unsigned int mat) const;

  /**
   * Perturbation Type
   */
  std::string pert_type;

  private:

  /**
   * @brief
   */
  void parse_dxs_file (const std::string &dxs_file);

  /**
   * @brief
   */
  void parse_dxs_XSEC (const std::string &dxs_file);

  /**
   * @brief
   */
  void parse_forest_dxs (const std::string &xml_file);

  /**
   * @brief
   */
  void parse_borders_file (const std::string &dxs_file);

  /**
   * @brief
   */
  void parse_borders_hex_file (const std::string &dxs_file);

  const Materials &materials;
  ConditionalOStream verbose_cout;
  const unsigned int dim;
  double freq;

  // Perturbation of Cross Sections
  // DXS[group][material]
  std::vector<std::vector<complex> > delta_sigma_f,
      delta_sigma_a, delta_sigma_r, delta_sigma_t; // TODO Remove sigma_a (make it local)
  std::vector<std::vector<std::vector<complex>>> delta_sigma_s;

  // Face Perturbation
  std::vector<std::vector<unsigned int> > faces_id;
  const unsigned int faces_map[3][2] =
                                         {
                                             { 2, 0 },
                                             { 1, 3 },
                                             { 5, 4 } };
};

#endif /* COMPLEX_PERTURBATION_H_ */
