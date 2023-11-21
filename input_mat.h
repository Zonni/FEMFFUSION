/**
 * @file   input/input_mat.h
 * @brief  InputMat class template declarations
 */
#ifndef FOREST_INPUT_MAT_H
#define FOREST_INPUT_MAT_H

// This is for the xml parser from boost

#include <map>
#include <string>
#include <vector>
#include <utility>

namespace XMLInput
{

/**
 * @brief Type of cross sections.
 *
 * @p XS_type is to know the cross sections needed for
 * this particular problem.
 *
 */
enum class XS_type
  : unsigned int
  {
    diffussion = 0, /**<  Diffusion cross sections. */
    transport = 1 /**<  Transport cross sections. */
};

/**
 * @brief Structure containing the data of a single material.
 *
 * Different data defining the properties of a single material,
 * and booleans variable to check which data has been provided
 * to this material or not.
 */
struct XS_single
{
  /** @todo document me */
  std::string name;
  /** @todo document me */
  std::vector<double> sigma_t;
  /** @todo document me */
  std::vector<double> nu_sigma_f;
  /** @todo document me */
  std::vector<double> chi;
  /** @todo document me */
  std::vector<double> lambda;
  /** @todo document me */
  std::vector<double> beta;
  /** @todo document me */
  double beta_eff;
  /** @todo document me */
  std::vector<double> velocities;
  /** @todo document me */
  std::vector<std::vector<double> > sigma_s;
  /** @todo document me */
  std::vector<std::vector<double> > chi_d;
  /** @todo document me */
  std::vector<double> chi_p;

  /** @todo document me */
  bool exist_sigma_t = false;
  /** @todo document me */
  bool exist_nu_sigma_f = false;
  /** @todo document me */
  bool exist_chi = false;
  /** @todo document me */
  bool exist_lambda = false;
  /** @todo document me */
  bool exist_beta = false;
  /** @todo document me */
  bool exist_velocities = false;
  /** @todo document me */
  bool exist_chi_d = false;
  /** @todo document me */
  bool exist_sigma_s = false;

  /** @todo document me */
  std::vector<double> sigma_a;
  /** @todo document me */
  std::vector<double> sigma_r;
  /** @todo document me */
  std::vector<double> sigma_tr;
  /** @todo document me */
  std::vector<double> nu;
  /** @todo document me */
  std::vector<double> sigma_f;

  /** @todo document me */
  bool exist_sigma_a = false;
  /** @todo document me */
  bool exist_nu = false;
  /** @todo document me */
  bool exist_sigma_f = false;

  /** @todo document me */
  unsigned int id;
};

/**
 * @class InputMat
 *
 * @brief Class for the cross sections.
 *
 * Here we have the cross sections for all the materials,
 * as well as some functions to read and write this data
 * to xml format, and some functions to check the data
 * when enough cross sections are available, i.e.,
 * we check that
 * \f$ \nu\Sigma_{f} = \nu*\Sigma_{f} \f$
 *  and we check
 * \f$ \Sigma_{t,g} = \sum_{h}(\Sigma_{s,g,h}) +
 * \Sigma_{a,g} \f$
 *
 * @todo The cross sections for the diffusion equation
 * should be added here.
 *
 */
class InputMat
{
public:

  /** @todo document me */
  typedef std::map<unsigned int, XS_single> XS_map;
  /** @todo document me */
  typedef std::pair<unsigned int, XS_single> XS_pair;

  /** @todo document me */
  XS_map xs;

  /**
   @brief Load the material data from the file @p filename
   @param filename
   */
  void
  load(const std::string &filename);

  /**
   @brief Print the material data to the file @p filename
   @param filename
   */
  void
  save(const std::string &filename);

  /**
   @brief Check Consistency of the materials data.
   @details check the data
   when enough cross sections are available, i.e.,
   we check that
   \f$ \nu\Sigma_{f} = \nu*\Sigma_{f} \f$
   and we check
   \f$ \Sigma_{t,g} = \sum_{h}(\Sigma_{s,g,h}) + \Sigma_{a,g} \f$
   */
  void
  check();

  /**
   * @brief Get the number of energy groups.
   * @return n_groups
   */
  unsigned int
  get_n_groups() const
  {
    return n_groups;
  }

  /**
   * @brief Get the number of different materials defined.
   * @return n_mat
   */
  unsigned int
  get_n_mat() const
  {
    return n_mat;
  }

  /**
   * @brief Get the number of different materials defined.
   * @return n_mat
   */
  unsigned int
  get_n_precursors() const
  {
    return n_precursors;
  }

  /**
   * @brief Get the number of different materials defined.
   * @return n_mat
   */
  std::vector<unsigned int>&
  get_materials_vector()
  {
    return materials_vector;
  }


private:

  unsigned int n_groups;
  unsigned int n_mat;
  unsigned int n_precursors;
  std::vector<unsigned int> materials_vector;

  /** @todo document me */
  void
  check_nusigf(XS_single & xs_) const;
  /** @todo document me */
  void
  calc_nusigf(XS_single & xs_);
  /** @todo document me */
  void
  check_sigmat(XS_single & xs_) const;
  /** @todo document me */
  void
  calc_sigmat(XS_single & xs_);
  /** @todo document me */
  void
  calc_chip(XS_single & xs_);
  /** @todo document me */
  void
  calc_betaeff (XS_single &xs_);
  /** @todo document me */
  void
  norm_chi(XS_single & xs_);

};

} // end of namespace XMLInput

#endif
