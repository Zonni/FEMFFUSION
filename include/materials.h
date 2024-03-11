/**
 * @file   materials.h
 * @brief  InputMat class template Materials
 */

#ifndef MATERIALS_H
#define MATERIALS_H

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

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>

#include "../include/femffusion.h"

using namespace dealii;

/**
 * @class Materials
 *
 * @brief Class for the cross sections.
 *
 *
 */
class Materials
{
  public:

  Materials (ConditionalOStream &verbose_cout);

  /**
   *
   */
  void reinit (const std::string &xsec_file,
    const std::string &xsec_type,
    const unsigned int n_groups,
    const std::vector<unsigned int> &n_assemblies_per_dim,
    unsigned int &n_assemblies,
    bool listen_to_material_id = false,
    std::vector<unsigned int> geo_ps = std::vector<unsigned int>(),
    const std::string precursors_file = "none");

  /**
   *
   */
  void reinit_xsec (std::string xs_name,
    unsigned int _n_mats,
    unsigned int _n_groups,
    std::vector<unsigned int> _materials_vector,
    std::vector<std::vector<double> > _nu_sigma_f);

  /**
   *
   */
  void reinit (unsigned int _n_mats,
    unsigned int _n_groups,
    std::vector<unsigned int> _materials_vector,
    std::vector<std::vector<std::vector<double>> > _new_xsec,
    std::vector<std::vector<std::vector<double>> > _new_xsec_s);

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
   * @brief Get the number of assemblies
   */
  unsigned int get_n_assemblies () const;

  /**
   * @brief Get the number of precursors groups.
   * @return n_precursors
   */
  unsigned int get_n_precursors () const;

  /**
   *
   */
  void set_n_precursors (const unsigned int n_prec);

  /**
   *
   */
  void set_velocity (
    const unsigned int group,
    double vel);

  /**
   *
   */
  void set_lambda_prec (
    const unsigned int group,
    double lambda);

  /**
   *
   */
  void set_beta_prec (
    const unsigned int group,
    double beta);

  /**
   * Set the default 2 group spectra: 1.0, 0.0
   */
  void set_default_spectra ();
  /**
   *
   */
  double get_lambda_prec (const unsigned int group) const;

  /**
   *
   */
  double get_beta_prec (const unsigned int group) const;

  /**
   *
   */
  double get_beta_total () const;

  /**
   *
   */
  std::vector<double> get_delayed_fraction_sum ();

  /**
   *
   */
  double get_delayed_fraction_sum (unsigned int mat) const;

  /**
   *
   */
  double get_velocity (const unsigned int mat,
    const unsigned int group) const;

  /**
   *
   */
  double get_delayed_decay_constant (
    const unsigned int mat,
    const unsigned int precursor) const;

  /**
   *
   */
  double get_delayed_fraction (
    const unsigned int mat,
    const unsigned int precursor) const;

  /**
   *
   */
  double get_delayed_spectra (
    const unsigned int mat,
    const unsigned int precursor,
    const unsigned int group) const;

  /**
   *
   */
  double get_prompt_spectra (
    const unsigned int mat,
    const unsigned int group) const;

  /**
   *
   */
  void remove_precursors ();

  /**
   * @brief Get a constant reference to materials_vector.
   */
  const std::vector<unsigned int>& get_materials_vector () const;

  /**
   * @brief Get a constant reference to get_geometry_matrix
   */
  const std::vector<std::vector<unsigned int> >& get_geometry_matrix () const;

  /**
   * @brief Get the material id of a cell given by its reference.
   * @return material_id
   */
  template <int dim>
    unsigned int get_material_id (
      typename DoFHandler<dim>::cell_iterator &cell) const;

  /**
   * @brief Get the original (before any vibration) material id of a cell given by its reference.
   * @return material_id
   */
  template <int dim>
    unsigned int get_original_material_id (
      typename DoFHandler<dim>::active_cell_iterator &cell) const;

  /**
   *
   */
  void set_materials_id (
    const unsigned int cell_user_index,
    const unsigned int mat_id);

  /**
   *
   */
  template <int dim>
    void get_materials_table (
      Table<dim, types::material_id> &materials_table,
      const std::vector<unsigned int> &n_assemblies_per_dim);

  /**
   * @brief Check if the cell in position (pos_x, pos_y, pos_z) is fuel
   * it also advances an index if the cell is a hole or the reflector.
   */
  bool is_fuel (
    const unsigned int pos_x,
    const unsigned int pos_y,
    const unsigned int pos_z,
    unsigned int &index) const;

  /**
   * @brief Check if the cell in position (pos_x, pos_y, pos_z) is reflector.
   */
  bool is_reflector (
    const unsigned int pos_x,
    const unsigned int pos_y,
    const unsigned int pos_z) const;

  /**
   * @brief Return the plane number where the cell_index pertains
   */
  unsigned int plane (int cell_index) const;

  /**
   * @brief Set geometry_matrix from the  geometry_points structure.
   */
  void set_geometry_matrix (
    const std::vector<unsigned int> &assem_per_dim,
    const std::vector<unsigned int> &geo_ps);

  /**
   * @brief Set geometry_matrix from a string (from the input file)
   */
  void set_geometry_matrix (
    const std::vector<unsigned int> &assem_per_dim,
    const std::string &str);

  /**
   * @brief Get the transport cross section given energy group and the material id.
   * @return sigma_tr
   */
  double get_sigma_tr (
    const unsigned int group,
    const unsigned int mat) const;

  /**
   * @brief Get the total cross section given energy group and the material id.
   * @return sigma_t
   */
  double get_sigma_t (
    const unsigned int group,
    const unsigned int mat) const;

  /**
   *
   */
  double get_chi (
    const unsigned int group,
    const unsigned int mat) const;

  /**
   *
   */
  double get_xi_nu_sigma_f (
    const unsigned int from_group,
    const unsigned int to_group,
    const unsigned int mat) const;

  /**
   *
   */
  double get_nu_sigma_f (
    const unsigned int group,
    const unsigned int mat) const;

  /**
   *
   */
  double get_sigma_f (
    const unsigned int from_group,
    const unsigned int mat) const;

  /**
   *
   */
  std::vector<double> get_nu_sigma_f (unsigned int group) const;

  /**
   *
   */
  std::vector<double> get_xi_nu_sigma_f (
    const unsigned int from_group,
    const unsigned int to_group) const;

  /**
   *
   */
  double get_sigma_s (
    const unsigned int from_group,
    const unsigned int to_group,
    const unsigned int mat) const;

  /**
   *
   */
  double get_sigma_r (
    const unsigned int group,
    unsigned int mat) const;

  /**
   *
   */
  std::vector<double> get_sigma_tr (unsigned int group) const;

  /**
   *
   */
  std::vector<double> get_sigma_r (unsigned int group) const;

  /**
   *
   */
  std::vector<double> get_sigma_f (const unsigned int group) const;

  /**
   *
   */
  std::vector<double> get_sigma_s (
    const unsigned int group_i,
    const unsigned int group_j) const;

  /**
   *
   */
  double get_diffusion_coefficient (
    const unsigned int group,
    const unsigned int mat) const;

  // ***********************************************************************************//
  // SOME SETTERS

  /**
   * @brief Set the set_sigma_f cross section of a given energy group and a material id.
   */
  void set_sigma_f (
    const double sigma_f,
    const unsigned int group,
    const unsigned int mat);

  /**
   *
   */
  void set_nu_sigma_f (
    const double nu_sigma_f,
    const unsigned int group,
    const unsigned int mat);

  /**
   *
   */
  void set_chi (
    const double chi,
    const unsigned int group,
    const unsigned int mat);

  /**
   *
   */
  void set_sigma_s (
    const double set_sigma_s,
    const unsigned int from_group,
    const unsigned int to_group,
    const unsigned int mat);

  /**
   *
   */
  void set_sigma_r (
    const double sigma_r,
    unsigned int group,
    unsigned int mat);

  /**
   *
   */
  void set_sigma_tr (
    const double sigma_r,
    unsigned int group,
    unsigned int mat);

  /**
   *
   */
  void set_sigma_t (
    const double sigma_r,
    unsigned int group,
    unsigned int mat);

  /*
   *
   */
  void change_mat_value (
    unsigned int old_material,
    unsigned int new_material);

  /**
   *
   */
  void create_new_mixed_mat (
    const unsigned int averaged_mat,
    const double frac,
    const unsigned int mat_bar,
    const unsigned int mat_no_bar,
    const unsigned int user_index);

  /*
   *
   */

  void create_new_mixed_mat_flux (
    const unsigned int new_mat,
    const double frac,
    const unsigned int mat_bar,
    const unsigned int mat_no_bar,
    const double flux_bar,
    const double flux_nobar,
    const unsigned int user_index);

  /**
   *
   */
  void make_critical (const double &keff);

  /*
   *
   */
  void modify_xsec (
    std::string xsec_type,
    unsigned int mat,
    std::vector<double> &delta_xsec);

  /**
   *
   */
  void modify_xsec_all (
    std::string xsec_type,
    unsigned int mat,
    std::vector<std::vector<double>> &new_xsec);

  /**
   *
   */
  void modify_xsec_7g (
    std::string xsec_type,
    double sim_time,
    std::vector<double> amplitudes,
    unsigned int n_mat);

  /**
   *
   */
  void modify_xsec_c5g7_td11 (double sim_time);

  /**
   *
   */
  void create_new_added_mat (
    const unsigned int new_mat,
    const double frac,
    double inc_xsec_r,
    const unsigned int mat_no_bar,
    const unsigned int user_index);

  /*
   *
   */
  void add_perturbation_xsec (
    std::string xsec_type,
    unsigned int mat,
    std::vector<double> &delta_xsec);

  /*
   *
   */
  void add_perturbation_xsec (
    std::string xsec_type,
    double &coeff);

  /*
   *
   *
   */

  void add_perturbation_xsec (
    std::string xsec_type,
    std::vector<double> &delta_xsec);

  /*
   *
   */
  void remove_perturbation_xsec (
    std::string xsec_type,
    unsigned int mat,
    std::vector<double> &delta_xsec);

  /*
   *
   */
  void remove_perturbation_xsec (
    std::string xsec_type,
    double &coeff);

  /*
   *
   *
   */

  void remove_perturbation_xsec (
    std::string xsec_type,
    std::vector<double> &delta_xsec);

  /*
   *
   */
  void save_initial_xsec ();

  /*
   *
   */
  void save_n_mats_init ();

  /**
   *
   */
  bool exist (
    const unsigned int idx_x,
    const unsigned int idx_y,
    const unsigned int idx_z) const;

  unsigned int n_total_assemblies;
  std::vector<unsigned int> assem_per_dim;
  unsigned int n_assemblies;

  double keff;
  unsigned int n_mats_init;
  bool transient;
  bool listen_to_material_id;

  std::vector<std::vector<double> > init_sigma_tr, init_sigma_t, init_nu_sigma_f,
      init_sigma_f, init_chi, init_sigma_r;
  // Group-to-group assemblies  XS[from_group][to_group][material]
  std::vector<std::vector<std::vector<double> > > init_sigma_s;

  private:

  /*
   * @brief Builds the deal.ii Table<dim, types::material_id> materials_table.
   */
  template <int dim>
    Table<dim, types::material_id> build_materials_table (
      const std::vector<std::vector<unsigned int> > &geometry_matrix,
      const std::vector<unsigned int> &n_cells_per_dim,
      const std::vector<unsigned int> &materials_vector);

  /**
   *
   */
  void parse_xsec (const std::string &xsec_file,
    const std::vector<unsigned int> &n_assemblies_per_dim,
    const unsigned int n_assemblies);

  /**
   *
   */
  void parse_precursors_file (const std::string &prec_file);

  /**
   *
   */
  void parse_xsec_2g (const std::string &xsec_file,
    const std::vector<unsigned int> &n_assemblies_per_dim,
    const unsigned int n_assemblies);

  /**
   *
   */
  void parse_forest_xs (const std::string &xml_file);

  /**
   *
   */
  void parse_Valkin_file (std::string xs_file,
    const std::vector<unsigned int> &n_assemblies_per_dim,
    const unsigned int n_assemblies,
    const std::vector<unsigned int> geo_ps);

  bool diffusions_coefficient_computed;

  unsigned int n_groups;
  unsigned int n_mats;
  unsigned int n_precursors;
  ConditionalOStream verbose_cout;

  // Cross Sections
  std::vector<std::vector<double> > nu_sigma_f, sigma_f, chi, sigma_r, sigma_tr, sigma_t; // XS[group][material]
  // Group-to-group assemblies  XS[from_group][to_group][material]
  std::vector<std::vector<std::vector<double> > > sigma_s;

  unsigned int n_prec_mat;
  std::vector<unsigned int> precursors_materials;
  std::vector<std::vector<double>> delayed_fractions, delayed_decay_constants;
  std::vector<std::vector<std::vector<double> > > delayed_spectra;
  std::vector<std::vector<double>> prompt_spectra;
  std::vector<double> delayed_fraction_sum;
  std::vector<std::vector<double>> velocities_vector;

  std::vector<unsigned int> materials_vector;
  std::vector<unsigned int> materials_vector_init;
  std::vector<unsigned int> materials_vector_no_bar;
  std::vector<unsigned int> materials_vector_with_holes;

  // 2D map of fuel cells, reflector and dummy cells
  std::vector<std::vector<unsigned int> > geometry_matrix;

  // Number of assemblies per plane
  std::vector<unsigned int> assemblies_per_plane;

};

#endif
