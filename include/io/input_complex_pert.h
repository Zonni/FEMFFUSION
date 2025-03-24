/**
 * @file
 * @brief
 */

#ifndef FOREST_COMPLEX_INPUT_PERT_H
#define FOREST_COMPLEX_INPUT_PERT_H

// This is for the xml parser from boost

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <complex>

typedef std::complex<double> complex;

namespace XMLInput
{

/**
 * @brief Type of cross sections.
 *
 * @p XS_type is to know the cross sections needed for
 * this particular problem.
 *
 */
enum class DXS_type
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
struct DXS_single
{
  /** @todo document me */
  std::string name;
  /** @todo document me */
  std::vector<complex> sigma_t;
  /** @todo document me */
  std::vector<complex> nu_sigma_f;
  /** @todo document me */
  std::vector<complex> chi;
  /** @todo document me */
  std::vector<std::vector<complex> > sigma_s;

  /** @todo document me */
  bool exist_sigma_t = false;
  /** @todo document me */
  bool exist_nu_sigma_f = false;
  /** @todo document me */
  bool exist_chi = false;
  /** @todo document me */
  bool exist_sigma_s = false;

  /** @todo document me */
  std::vector<complex> sigma_a;
  /** @todo document me */
  std::vector<complex> nu;
  /** @todo document me */
  std::vector<complex> sigma_f;

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
class InputPert
{
  public:

  typedef std::map<unsigned int, DXS_single> DXS_map;
  typedef std::pair<unsigned int, DXS_single> DXS_pair;

  DXS_map xs;

  /**
   @brief Load the material data from the file @p filename
   @param filename
   */
  void load (const std::string &filename);

  /**
   @brief Check Consistency of the materials data.
   @details check the data
   when enough cross sections are available, i.e.,
   we check that
   \f$ \nu\Sigma_{f} = \nu*\Sigma_{f} \f$
   and we check
   \f$ \Sigma_{t,g} = \sum_{h}(\Sigma_{s,g,h}) + \Sigma_{a,g} \f$
   */
  void check ();

  /**
   * @brief Get the number of energy groups.
   * @return n_groups
   */
  unsigned int get_n_groups () const
  {
    return n_groups;
  }

  /**
   * @brief Get the number of different materials defined.
   * @return n_mat
   */
  unsigned int get_n_mat () const
  {
    return n_mat;
  }

  /**
   * @brief Get the number of different materials defined.
   * @return n_mat
   */
  std::vector<unsigned int>& get_materials_vector ()
  {
    return materials_vector;
  }

  private:

  unsigned int n_groups;
  unsigned int n_mat;
  std::vector<unsigned int> materials_vector;

  /** Check if NuSigma_f = Nu*Sigma_f */
  void check_nusigf (DXS_single &xs_) const;

  /**  document me */
  void calc_nusigf (DXS_single &xs_);

  /** @todo document me */
  void check_sigmat (DXS_single &xs_) const;

  /** @todo document me */
  void calc_sigmat (DXS_single &xs_);
};

} // end of namespace XMLInput

#endif
