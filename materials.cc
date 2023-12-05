/**
 * @file   materials.cc
 * @brief  Implementation of class Materials
 */

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include <deal.II/lac/full_matrix.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

#include "materials.h"
#include "input_mat.h"

using namespace dealii;

/**
 *
 */
Materials::Materials (ConditionalOStream &verbose_cout) :
    verbose_cout(verbose_cout)
{
  n_groups = 0;
  n_mats = 0;
  diffusions_coefficient_computed = false;
  listen_to_material_id = false;
  n_precursors = 0;
//  beta_total = 0.0;
  n_mats_init = 0;
  n_total_assemblies = 0;
  n_prec_mat = 0;
  keff = 1.0;
  n_assemblies = 0;
  transient = false;
}

/**
 *
 */
void Materials::reinit (const std::string &xsec_file,
  const std::string &xsec_type,
  const unsigned int n_groups,
  const std::vector<unsigned int> &n_assemblies_per_dim,
  unsigned int &n_assem,
  bool listen_to_material_id,
  std::vector<unsigned int> geo_ps,
  const std::string precursors_file)
{
  this->n_groups = n_groups;
  this->listen_to_material_id = listen_to_material_id;

  assem_per_dim = n_assemblies_per_dim;
  n_assemblies = n_assem;
  AssertRelease(fexists(xsec_file), "XSECS_Filename does not exist");

  if (xsec_type == "XS2G")
  {
    parse_xsec_2g(xsec_file, n_assemblies_per_dim, n_assem);
    AssertRelease(n_groups == 2, "Only 2 energy groups valid for this XS");
  }
  else if (xsec_type == "XSEC")
    parse_xsec(xsec_file, n_assemblies_per_dim, n_assem);
  else if (xsec_type == "XML")
  {
    parse_forest_xs(xsec_file);
    if (precursors_file != "none")
      parse_precursors_file(precursors_file);
  }
  else if (xsec_type == "Valkin")
    parse_Valkin_file(xsec_file, n_assemblies_per_dim, n_assem, geo_ps);
  else
    AssertRelease(false, "Wrong xsec_type " + xsec_type);

  compute_diffusions_coefficients();

  const unsigned int per_plane = n_assem / n_assemblies_per_dim[2];
  assemblies_per_plane.resize(n_assemblies_per_dim[2], per_plane);

  n_total_assemblies = n_assemblies;

  // Delete dummy materials
  std::vector<unsigned int> temp;
  temp.reserve(materials_vector.size());
  for (unsigned int i = 0; i < materials_vector.size(); i++)
  {
    if (materials_vector[i] != static_cast<unsigned int>(-1))
      temp.push_back(materials_vector[i]);
    else
    {
      verbose_cout << "Hole in position: " << i << " in plane "
                   << i / per_plane
                   << std::endl;
      assemblies_per_plane[i / per_plane]--;
      n_assemblies--;
    }
  }

  // Remove materials with holes
  materials_vector_with_holes = materials_vector;
  materials_vector_no_bar = materials_vector;
  materials_vector = temp;
}

/**
 *
 */
void Materials::reinit_xsec (std::string xs_name,
  unsigned int _n_mats,
  unsigned int _n_groups,
  std::vector<unsigned int> _materials_vector,
  std::vector<std::vector<double> > _new_xsec)
{

  n_mats = _n_mats;
  n_groups = _n_groups;
  materials_vector = _materials_vector;

  if (xs_name == "sigma_f")
  {
    nu_sigma_f = _new_xsec;
    chi.resize(n_groups, std::vector<double>(n_mats));
    chi[0].assign(n_mats, 1.0);
    chi[1].assign(n_mats, 0.0);
  }
  else if (xs_name == "sigma_a")
  {
    sigma_r = _new_xsec;
  }
  else
    AssertRelease(false, "Invalid type of xs_name");
}

/**
 *
 */
void Materials::reinit (unsigned int _n_mats,
  unsigned int _n_groups,
  std::vector<unsigned int> _materials_vector,
  std::vector<std::vector<std::vector<double>> > _new_xsec,
  std::vector<std::vector<std::vector<double>> > _new_xsec_s)
{

  n_mats = _n_mats;
  n_groups = _n_groups;
  materials_vector = _materials_vector;

  // 0: diffusion, 1: sigma_r, 2: nusigmaf, 3: sigma_f,
  D = _new_xsec[0];
  sigma_r = _new_xsec[1];

  nu_sigma_f = _new_xsec[2];
  sigma_f = _new_xsec[3];
  sigma_s = _new_xsec_s;

  chi.resize(n_groups, std::vector<double>(n_mats));
  AssertRelease(n_groups == 2, "Change this function for more groups");
  chi[0].assign(n_mats, 1.0);
  chi[1].assign(n_mats, 0.0);

}

/**
 * @brief Get the number of energy groups.
 * @return n_groups
 */
unsigned int Materials::get_n_groups () const
{
  return n_groups;
}

/**
 * @brief Get the number of materials defined.
 * @return n_mats
 */
unsigned int Materials::get_n_mats () const
{
  return n_mats;
}

/**
 * @brief Get a constant reference to geometry_matrix
 */
unsigned int Materials::get_n_assemblies () const
{
  return materials_vector.size();
}

/**
 * @brief Get the number of energy groups.
 * @return n_precursors
 */
unsigned int Materials::get_n_precursors () const
{
  return n_precursors;
}

///**
// *
// */
//double Materials::get_velocitiy (const unsigned int group) const
//{
//  AssertIndexRange(group, velocities.size());
//  return velocities[group];
//}
//
///**
// *
// */
//double Materials::get_lambda_prec (const unsigned int group) const
//{
//  AssertIndexRange(group, lambda_prec.size());
//  return lambda_prec[group];
//}
//
///**
// *
// */
//double Materials::get_beta_prec (const unsigned int group) const
//{
//  AssertIndexRange(group, beta_prec.size());
//  return beta_prec[group];
//}
//
///**
// *
// */
//double Materials::get_beta_total () const
//{
//
//  return beta_total;
//}

/**
 *  Generalization of precursors
 */
std::vector<double> Materials::get_delayed_fraction_sum ()
{

  return delayed_fraction_sum;
}

/**
 *
 */
double Materials::get_delayed_fraction_sum (unsigned int mat) const
{

  return delayed_fraction_sum[mat];
}

/**
 *
 */
double Materials::get_velocitiy (const unsigned int mat,
  const unsigned int group) const
{
  AssertIndexRange(mat, velocities_vector.size());
  return velocities_vector[mat][group];
}

/**
 *
 */
void Materials::set_n_precursors (const unsigned int n_prec)
{
  n_precursors = n_prec;
//  lambda_prec.resize(n_prec);
//  beta_prec.resize(n_prec);

  delayed_decay_constants.resize(n_mats);
  for (unsigned int mat = 0; mat < n_mats; mat++)
    delayed_decay_constants[mat].resize(n_prec);

  delayed_fraction_sum.resize(n_mats);

  delayed_fractions.resize(n_mats);
  for (unsigned int mat = 0; mat < n_mats; mat++)
    delayed_fractions[mat].resize(n_prec);

}

/**
 *
 */
void Materials::set_velocity (const unsigned int group,
  double vel)
{
//  velocities.resize(n_groups);
  AssertIndexRange(group, n_groups);
//  velocities[group] = vel;

  velocities_vector.resize(n_mats);
  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    velocities_vector[mat].resize(n_groups);
    velocities_vector[mat][group] = vel;
  }

}

/**
 *
 */
void Materials::set_lambda_prec (const unsigned int group,
  double lambda)
{
  AssertIndexRange(group, lambda_prec.size());
//  lambda_prec[group] = lambda;

  for (unsigned int mat = 0; mat < n_mats; mat++)
    delayed_decay_constants[mat][group] = lambda;
}

void Materials::set_default_spectra ()
{
  // resize
  delayed_spectra.resize(n_mats);
  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    delayed_spectra[mat].resize(n_precursors);
    for (unsigned int p = 0; p < n_precursors; p++)
      delayed_spectra[mat][p].resize(n_groups);
  }
  prompt_spectra.resize(n_mats);
  for (unsigned int mat = 0; mat < n_mats; mat++)
    prompt_spectra[mat].resize(n_groups);

  // Set it
  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    for (unsigned int p = 0; p < n_precursors; p++)
      delayed_spectra[mat][p][0] = 1.0;
    // Next are 0.0
  }

  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    prompt_spectra[mat][0] = 1.0;
    // Next are 0.0
  }
}

/**
 *
 */
void Materials::set_beta_prec (const unsigned int group,
  double beta)
{
  AssertIndexRange(group, beta_prec.size());
//  beta_prec[group] = beta;
//  beta_total = sum_vector(beta_prec);

  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    delayed_fractions[mat][group] = beta;
    delayed_fraction_sum[mat] = sum_vector(delayed_fractions[mat]);
  }
}
/**
 *
 */
double Materials::get_delayed_decay_constant (const unsigned int mat,
  const unsigned int precursor) const
{
  AssertIndexRange(mat, delayed_decay_constants.size());
  return delayed_decay_constants[mat][precursor];
}

/**
 *
 */
double Materials::get_delayed_fraction (const unsigned int mat,
  const unsigned int precursor) const
{
  AssertIndexRange(mat, delayed_fractions.size());
  return delayed_fractions[mat][precursor];
}

/**
 *
 */
double Materials::get_delayed_spectra (const unsigned int mat,
  const unsigned int precursor,
  const unsigned int group) const
{
  AssertIndexRange(mat, delayed_spectra.size());
  return delayed_spectra[mat][precursor][group];
}

/**
 *
 */
double Materials::get_prompt_spectra (const unsigned int mat,
  const unsigned int group) const
{
  AssertIndexRange(mat, prompt_spectra.size());
  return prompt_spectra[mat][group];
}

/**
 *
 */
void Materials::remove_precursors ()
{

  n_precursors = 0;
//  beta_total = 0;

  for (unsigned int mat = 0; mat < n_mats; mat++)
    delayed_fraction_sum[mat] = 0.0;

}

/**
 *
 */
const std::vector<unsigned int>& Materials::get_materials_vector () const
{
  return materials_vector;
}

/**
 *
 */
template <int dim>
  unsigned int Materials::get_material_id (
    typename DoFHandler<dim>::cell_iterator &cell) const
  {
    if (listen_to_material_id)
    {
      return cell->material_id();
    }
    else
    {
      AssertIndexRange(cell->user_index(), materials_vector.size());
      return materials_vector[cell->user_index()];
    }
  }

template unsigned int Materials::get_material_id<1> (
  typename DoFHandler<1>::cell_iterator&) const;
template unsigned int Materials::get_material_id<2> (
  typename DoFHandler<2>::cell_iterator&) const;
template unsigned int Materials::get_material_id<3> (
  typename DoFHandler<3>::cell_iterator&) const;

/**
 *
 */
template <int dim>
  unsigned int Materials::get_original_material_id (
    typename DoFHandler<dim>::active_cell_iterator &cell) const
  {
    AssertIndexRange(cell->user_index(), materials_vector_no_bar.size());
    return materials_vector_no_bar[cell->user_index()];
  }

template
unsigned int Materials::get_original_material_id<1> (
  typename DoFHandler<1>::active_cell_iterator&) const;
template
unsigned int Materials::get_original_material_id<2> (
  typename DoFHandler<2>::active_cell_iterator&) const;
template
unsigned int Materials::get_original_material_id<3> (
  typename DoFHandler<3>::active_cell_iterator&) const;

/**
 *
 */
void Materials::set_materials_id (const unsigned int cell_user_index,
  const unsigned int mat_id)
{
  AssertIndexRange(cell_user_index, materials_vector.size());
  AssertIndexRange(mat_id, n_mats);
  materials_vector[cell_user_index] = mat_id;
}

/**
 *
 */
double Materials::get_sigma_t (unsigned int group,
  unsigned int mat) const
{
  AssertIndexRange(group, sigma_t.size());
  AssertIndexRange(mat, sigma_t[group].size());
  return sigma_t[group][mat];
}

/**
 *
 */
double Materials::get_sigma_tr (unsigned int group,
  unsigned int mat) const
{
  AssertIndexRange(group, sigma_t.size());
  AssertIndexRange(mat, sigma_t[group].size());
  return sigma_t[group][mat]; // TODO
}

/**
 *
 */
double Materials::get_chi (unsigned int group,
  unsigned int mat) const
{
  AssertIndexRange(group, nu_sigma_f.size());
  AssertIndexRange(mat, nu_sigma_f[group].size());
  return nu_sigma_f[group][materials_vector[mat]];
}

/**
 *
 */
double Materials::get_sigma_r (const unsigned int group,
  const unsigned int mat) const
{
  AssertIndexRange(group, sigma_r.size());
  AssertIndexRange(mat, sigma_r[group].size());
  return sigma_r[group][mat];
}

/**
 *
 */
std::vector<double> Materials::get_sigma_r (const unsigned int group) const
{
  AssertIndexRange(group, sigma_r.size());
  return sigma_r[group];
}

/**
 *
 */
std::vector<double> Materials::get_sigma_f (const unsigned int group) const
{
  AssertIndexRange(group, sigma_f.size());
  return sigma_f[group];
}

/**
 *
 */
std::vector<double> Materials::get_sigma_t (const unsigned int group) const
{
  AssertIndexRange(group, sigma_t.size());
  return sigma_t[group];
}

/**
 *
 */
std::vector<double> Materials::get_sigma_s (const unsigned int group_i,
  const unsigned int group_j) const
{
  AssertIndexRange(group_i, sigma_s.size());
  return sigma_s[group_i][group_j];
}

/**
 *
 */
double Materials::get_xi_nu_sigma_f (const unsigned int from_group,
  const unsigned int to_group,
  const unsigned int mat) const
{
  AssertIndexRange(to_group, chi.size());
  AssertIndexRange(from_group, nu_sigma_f.size());
  AssertIndexRange(mat, chi[to_group].size());
  AssertIndexRange(mat, nu_sigma_f[from_group].size());

  return chi[to_group][mat] * nu_sigma_f[from_group][mat];
}

/**
 *
 */
double Materials::get_nu_sigma_f (const unsigned int group,
  const unsigned int mat) const
{

  AssertIndexRange(group, nu_sigma_f.size());
  AssertIndexRange(mat, nu_sigma_f[group].size());

  return nu_sigma_f[group][mat];
}

/**
 *
 */
double Materials::get_sigma_f (const unsigned int group,
  const unsigned int mat) const
{
  AssertIndexRange(group, sigma_f.size());
  AssertIndexRange(mat, sigma_f[group].size());
  return sigma_f[group][mat];
}

/**
 *
 */
std::vector<double> Materials::get_xi_nu_sigma_f (const unsigned int from_group,
  const unsigned int to_group) const
{
  AssertIndexRange(to_group, chi.size());
  AssertIndexRange(from_group, nu_sigma_f.size());

  std::vector<double> chinusigmaf(n_mats);
  for (unsigned int nm = 0; nm < n_mats; nm++)
    chinusigmaf[nm] = chi[to_group][nm] * nu_sigma_f[from_group][nm];

  return chinusigmaf;
}

/**
 *
 */
std::vector<double> Materials::get_nu_sigma_f (const unsigned int group) const
{

  std::vector<double> chinusigmaf(n_mats);
  for (unsigned int nm = 0; nm < n_mats; nm++)
    chinusigmaf[nm] = nu_sigma_f[group][nm];

  return chinusigmaf;
}

/**
 *
 */
double Materials::get_sigma_s (const unsigned int from_group,
  const unsigned int to_group,
  const unsigned int mat) const
{
  AssertIndexRange(from_group, sigma_s.size());
  AssertIndexRange(to_group, sigma_s[from_group].size());
  AssertIndexRange(mat, sigma_s[from_group][to_group].size());
  return sigma_s[from_group][to_group][mat];
}

/**
 *
 */
double Materials::get_diffusion_coefficient (const unsigned int group,
  const unsigned int mat) const
{
  AssertIndexRange(group, D.size());
  AssertIndexRange(mat, D[group].size());
  return D[group][mat];
}

void Materials::set_sigma_f (const double sigma_f_coeff,
  const unsigned int group,
  const unsigned int mat)
{
  AssertIndexRange(group, sigma_f.size());
  AssertIndexRange(mat, sigma_f[group].size());
  sigma_f[group][mat] = sigma_f_coeff;
}

/**
 *
 */
void Materials::set_nu_sigma_f (const double nu_sigma_f_coeff,
  const unsigned int group,
  const unsigned int mat)
{
  AssertIndexRange(group, nu_sigma_f.size());
  AssertIndexRange(mat, nu_sigma_f[group].size());
  nu_sigma_f[group][mat] = nu_sigma_f_coeff;
}

/**
 *
 */
void Materials::set_chi (const double chi_coeff,
  const unsigned int group,
  const unsigned int mat)
{
  AssertIndexRange(group, chi.size());
  AssertIndexRange(mat, chi[group].size());
  chi[group][mat] = chi_coeff;
}

/**
 *
 */
void Materials::set_sigma_s (const double set_sigma_s_coeff,
  const unsigned int from_group,
  const unsigned int to_group,
  const unsigned int mat)
{
  AssertIndexRange(from_group, sigma_s.size());
  AssertIndexRange(to_group, sigma_s[from_group].size());
  AssertIndexRange(mat, sigma_s[from_group][to_group].size());
  sigma_s[from_group][to_group][mat] = set_sigma_s_coeff;
}

/**
 *
 */
void Materials::set_sigma_r (const double sigma_r_coeff,
  unsigned int group,
  unsigned int mat)
{
  AssertIndexRange(group, sigma_r.size());
  AssertIndexRange(mat, sigma_r[group].size());
  sigma_r[group][mat] = sigma_r_coeff;
}

/**
 *
 */
void Materials::set_diffusion_coefficient (const double diff_coeff,
  const unsigned int group,
  const unsigned int mat)
{
  AssertIndexRange(group, D.size());
  AssertIndexRange(mat, D[group].size());
  D[group][mat] = diff_coeff;
}

void Materials::change_mat_value (unsigned int old_material,
  unsigned int new_material)
{

  for (unsigned int c = 0; c < materials_vector.size(); c++)
    if (materials_vector[c] == old_material)
      materials_vector[c] = new_material;
}

/**
 *
 */
void Materials::create_new_mixed_mat (const unsigned int new_mat,
  const double frac,
  const unsigned int mat_bar,
  const unsigned int mat_no_bar,
  const unsigned int user_index)
{
  Assert(n_groups == 2, ExcMessage("Implemented to n_groups"));
  Assert(frac >= 0.0 and frac <= 1.0,
    ExcMessage("Invalid frac " + num_to_str(frac)));

  if (n_mats == new_mat)
  {
    n_mats++;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      D[g].resize(n_mats);
      sigma_t[g].resize(n_mats);
      sigma_r[g].resize(n_mats);
      nu_sigma_f[g].resize(n_mats);
      chi[g].resize(n_mats);
      sigma_f[g].resize(n_mats);
      for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
        sigma_s[g][to_g].resize(n_mats);
    }

    delayed_fractions.resize(n_mats, std::vector<double>(n_precursors));
    delayed_decay_constants.resize(n_mats,
      std::vector<double>(n_precursors));
    prompt_spectra.resize(n_mats, std::vector<double>(n_precursors));
    velocities_vector.resize(n_mats, std::vector<double>(n_groups));
    delayed_fraction_sum.resize(n_mats);
    delayed_spectra.resize(n_mats,
      std::vector<std::vector<double> >(n_precursors,
        std::vector<double>(n_groups)));
  }

  for (unsigned int g = 0; g < n_groups; ++g)
  {
    double tr_bar = 0.0, tr_nobar = 0.0;
    tr_bar = 1 / (3 * D[g][mat_bar]);
    tr_nobar = 1 / (3 * D[g][mat_no_bar]);

    sigma_t[g][new_mat] = frac * tr_bar + (1 - frac) * tr_nobar;

    D[g][new_mat] = 1 / (3 * sigma_t[g][new_mat]);

//    D[g][new_mat] = frac * D[g][mat_bar]
//                    + (1 - frac) * D[g][mat_no_bar];

    sigma_r[g][new_mat] = frac * sigma_r[g][mat_bar]
                          + (1 - frac) * sigma_r[g][mat_no_bar];

    nu_sigma_f[g][new_mat] = frac * nu_sigma_f[g][mat_bar]
                             + (1 - frac) * nu_sigma_f[g][mat_no_bar];

    sigma_f[g][new_mat] = frac * sigma_f[g][mat_bar]
                          + (1 - frac) * sigma_f[g][mat_no_bar];

    chi[g][new_mat] = frac * chi[g][mat_bar]
                      + (1 - frac) * chi[g][mat_no_bar];

    for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
      sigma_s[g][to_g][new_mat] = frac * sigma_s[g][to_g][mat_bar]
                                  + (1 - frac) * sigma_s[g][to_g][mat_no_bar];

    for (unsigned int p = 0; p < n_precursors; p++)
      delayed_spectra[new_mat][p][g] = frac
                                       * delayed_spectra[mat_bar][p][g]
                                       + (1 - frac) * delayed_spectra[mat_no_bar][p][g];

    velocities_vector[new_mat][g] = frac * velocities_vector[mat_bar][g]
                                    + (1 - frac) * velocities_vector[mat_no_bar][g];

  }

  for (unsigned int p = 0; p < n_precursors; ++p)
  {
    delayed_fractions[new_mat][p] = frac * delayed_fractions[mat_bar][p]
                                    + (1 - frac) * delayed_fractions[mat_no_bar][p];
    delayed_decay_constants[new_mat][p] =
        frac
        * delayed_decay_constants[mat_bar][p]
        + (1 - frac) * delayed_decay_constants[mat_no_bar][p];
    prompt_spectra[new_mat][p] = frac * prompt_spectra[mat_bar][p]
                                 + (1 - frac) * prompt_spectra[mat_no_bar][p];
  }

  delayed_fraction_sum[new_mat] = frac * delayed_fraction_sum[mat_bar]
                                  + (1 - frac) * delayed_fraction_sum[mat_no_bar];

//  AssertIndexRange(user_index, materials_vector.size());
  materials_vector[user_index] = new_mat;

}

/**
 *
 */
void Materials::create_new_mixed_mat_flux (const unsigned int new_mat,
  const double frac,
  const unsigned int mat_bar,
  const unsigned int mat_no_bar,
  const double flux_bar,
  const double flux_nobar,
  const unsigned int user_index)
{
  Assert(n_groups == 2, ExcMessage("Implemented to n_groups"));
  Assert(frac >= 0.0 and frac <= 1.0,
    ExcMessage("Invalid frac " + num_to_str(frac)));

  if (n_mats == new_mat)
  {
    n_mats++;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      D[g].resize(n_mats);
      sigma_t[g].resize(n_mats);
      sigma_r[g].resize(n_mats);
      nu_sigma_f[g].resize(n_mats);
      chi[g].resize(n_mats);
      sigma_f[g].resize(n_mats);
      for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
        sigma_s[g][to_g].resize(n_mats);
    }

    delayed_fractions.resize(n_mats, std::vector<double>(n_precursors));
    delayed_decay_constants.resize(n_mats,
      std::vector<double>(n_precursors));
    prompt_spectra.resize(n_mats, std::vector<double>(n_precursors));
    velocities_vector.resize(n_mats, std::vector<double>(n_groups));
    delayed_fraction_sum.resize(n_mats);
    delayed_spectra.resize(n_mats,
      std::vector<std::vector<double> >(n_precursors,
        std::vector<double>(n_groups)));

  }

  double den_flux = frac * flux_bar + (1 - frac) * flux_nobar;

  for (unsigned int g = 0; g < n_groups; ++g)
  {

    //		double tr_bar = 0.0, tr_nobar = 0.0;
    //		tr_bar = 1 / (3 * D[g][mat_bar]);
    //		tr_nobar = 1 / (3 * D[g][mat_no_bar]);

    sigma_t[g][new_mat] = frac * sigma_t[g][mat_bar] * flux_bar
                          + (1 - frac) * sigma_t[g][mat_no_bar] * flux_nobar;
    sigma_t[g][new_mat] /= den_flux;

    D[g][new_mat] = 1 / (3 * sigma_t[g][new_mat]);

    sigma_r[g][new_mat] = frac * sigma_r[g][mat_bar] * flux_bar
                          + (1 - frac) * sigma_r[g][mat_no_bar] * flux_nobar;
    sigma_r[g][new_mat] /= den_flux;

    nu_sigma_f[g][new_mat] = frac * nu_sigma_f[g][mat_bar] * flux_bar
                             + (1 - frac) * nu_sigma_f[g][mat_no_bar] * flux_nobar;
    nu_sigma_f[g][new_mat] /= den_flux;

    sigma_f[g][new_mat] = frac * sigma_f[g][mat_bar] * flux_bar
                          + (1 - frac) * sigma_f[g][mat_no_bar] * flux_nobar;
    sigma_f[g][new_mat] /= den_flux;

    chi[g][new_mat] = frac * chi[g][mat_bar] * flux_bar
                      + (1 - frac) * chi[g][mat_no_bar] * flux_nobar;
    chi[g][new_mat] /= den_flux;

    for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
    {
      sigma_s[g][to_g][new_mat] = frac * sigma_s[g][to_g][mat_bar]
                                  * flux_bar
                                  + (1 - frac) * sigma_s[g][to_g][mat_no_bar]
                                    * flux_nobar;
      sigma_s[g][to_g][new_mat] /= den_flux;

    }

    for (unsigned int p = 0; p < n_precursors; p++)
    {
      delayed_spectra[new_mat][p][g] = frac
                                       * delayed_spectra[mat_bar][p][g]
                                       * flux_bar
                                       + (1 - frac) * delayed_spectra[mat_no_bar][p][g]
                                         * flux_nobar;
      delayed_spectra[new_mat][p][g] /= den_flux;
    }

    velocities_vector[new_mat][g] = frac * velocities_vector[mat_bar][g]
                                    * flux_bar
                                    + (1 - frac) * velocities_vector[mat_no_bar][g]
                                      * flux_nobar;
    velocities_vector[new_mat][g] /= den_flux;

  }

  for (unsigned int p = 0; p < n_precursors; ++p)
  {
    delayed_fractions[new_mat][p] = frac * delayed_fractions[mat_bar][p]
                                    * flux_bar
                                    + (1 - frac) * delayed_fractions[mat_no_bar][p]
                                      * flux_nobar;
    delayed_fractions[new_mat][p] /= den_flux;

    delayed_decay_constants[new_mat][p] =
        frac
        * delayed_decay_constants[mat_bar][p]
        * flux_bar
        + (1 - frac) * delayed_decay_constants[mat_no_bar][p]
          * flux_nobar;
    delayed_decay_constants[new_mat][p] /= den_flux;

    prompt_spectra[new_mat][p] = frac * prompt_spectra[mat_bar][p]
                                 * flux_bar
                                 + (1 - frac) * prompt_spectra[mat_no_bar][p]
                                   * flux_nobar;
    prompt_spectra[new_mat][p] /= den_flux;
  }

  delayed_fraction_sum[new_mat] = frac * delayed_fraction_sum[mat_bar]
                                  * flux_bar
                                  + (1 - frac) * delayed_fraction_sum[mat_no_bar]
                                    * flux_nobar;
  delayed_fraction_sum[new_mat] /= den_flux;

  //  AssertIndexRange(user_index, materials_vector.size());
  materials_vector[user_index] = new_mat;

}

/**
 *
 */
void Materials::create_new_mixed_mat_toni (const unsigned int new_mat,
  const double frac,
  const unsigned int mat_bar,
  const unsigned int mat_no_bar,
  const unsigned int user_index)
{

  Assert(n_groups == 2, ExcMessage("Implemented to n_groups"));
  Assert(frac >= 0.0 and frac <= 1.0,
    ExcMessage("Invalid fraction " + num_to_str(frac)));
  AssertIndexRange(mat_bar, n_mats);
  AssertIndexRange(mat_no_bar, n_mats);
  AssertIndexRange(new_mat, n_mats+1);
  // If necessary add a reserve space for a new material

  if (n_mats == new_mat)
  {
    n_mats++;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      D[g].resize(n_mats);
      sigma_t[g].resize(n_mats);
      sigma_r[g].resize(n_mats);
      nu_sigma_f[g].resize(n_mats);
      chi[g].resize(n_mats);
      sigma_f[g].resize(n_mats);
      for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
        sigma_s[g][to_g].resize(n_mats);
    }

    delayed_fractions.resize(n_mats, std::vector<double>(n_precursors));
    delayed_decay_constants.resize(n_mats,
      std::vector<double>(n_precursors));
    prompt_spectra.resize(n_mats, std::vector<double>(n_precursors));
    velocities_vector.resize(n_mats, std::vector<double>(n_groups));
    delayed_fraction_sum.resize(n_mats);
    delayed_spectra.resize(n_mats,
      std::vector<std::vector<double> >(n_precursors,
        std::vector<double>(n_groups)));

  }

  for (unsigned int g = 0; g < n_groups; ++g)
  {

    sigma_t[g][new_mat] = frac * sigma_t[g][mat_bar]
                          + (1 - frac) * sigma_t[g][mat_no_bar];
    D[g][new_mat] = 1 / (3 * sigma_t[g][new_mat]);
    //		D[g][new_mat] = frac * D[g][mat_bar] + (1 - frac) * D[g][mat_no_bar];
    sigma_r[g][new_mat] = frac * sigma_r[g][mat_bar]
                          + (1 - frac) * sigma_r[g][mat_no_bar];
    nu_sigma_f[g][new_mat] = frac * nu_sigma_f[g][mat_bar]
                             + (1 - frac) * nu_sigma_f[g][mat_no_bar];
    chi[g][new_mat] = frac * chi[g][mat_bar]
                      + (1 - frac) * chi[g][mat_no_bar];
    sigma_f[g][new_mat] = frac * sigma_f[g][mat_bar]
                          + (1 - frac) * sigma_f[g][mat_no_bar];
    for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
      sigma_s[g][to_g][new_mat] = frac * sigma_s[g][to_g][mat_bar]
                                  + (1 - frac) * sigma_s[g][to_g][mat_no_bar];

    for (unsigned int p = 0; p < n_precursors; p++)
      delayed_spectra[new_mat][p][g] = frac
                                       * delayed_spectra[mat_bar][p][g]
                                       + (1 - frac) * delayed_spectra[mat_no_bar][p][g];

    velocities_vector[new_mat][g] = frac * velocities_vector[mat_bar][g]
                                    + (1 - frac) * velocities_vector[mat_no_bar][g];
  }

  for (unsigned int p = 0; p < n_precursors; ++p)
  {
    delayed_fractions[new_mat][p] = frac * delayed_fractions[mat_bar][p]
                                    + (1 - frac) * delayed_fractions[mat_no_bar][p];
    delayed_decay_constants[new_mat][p] =
        frac
        * delayed_decay_constants[mat_bar][p]
        + (1 - frac) * delayed_decay_constants[mat_no_bar][p];
    prompt_spectra[new_mat][p] = frac * prompt_spectra[mat_bar][p]
                                 + (1 - frac) * prompt_spectra[mat_no_bar][p];
  }

  delayed_fraction_sum[new_mat] = frac * delayed_fraction_sum[mat_bar]
                                  + (1 - frac) * delayed_fraction_sum[mat_no_bar];

  // Create new position in materials_vector
  AssertIndexRange(user_index, materials_vector.size()+1);
  if (user_index == materials_vector.size())
    materials_vector.push_back(new_mat);
  else
    materials_vector[user_index] = new_mat;

}

/**
 * @brief Make the reactor critical.
 */
void Materials::make_critical (const double &keffective)
{

  // Compute diffusion Coefficients
  for (unsigned int g = 0; g < n_groups; ++g)
    for (unsigned int mat = 0; mat < n_mats; ++mat)
    {
      nu_sigma_f[g][mat] /= keffective;
      sigma_f[g][mat] /= keffective;
    }

  keff = keffective;

}

/**
 *
 */
void Materials::modify_xsec (std::string xsec_type,
  unsigned int mat,
  std::vector<double> &new_xsec)
{

  for (unsigned int ng = 0; ng < n_groups; ng++)
  {
    if (xsec_type == "sigma_f")
    {
      double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
      nu_sigma_f[ng][mat] = new_xsec[ng];
      sigma_f[ng][mat] = new_xsec[ng] / nu;
    }
    else if (xsec_type == "sigma_a")
    {
      sigma_r[ng][mat] = new_xsec[ng];
    }
    else
    {
      AssertRelease(false, "Invalid type of xsec");
    }
  }

  return;
}

/**
 *
 */
void Materials::modify_xsec_all (std::string xsec_type,
  unsigned int mat,
  std::vector<std::vector<double>> &new_xsec)
{

  AssertRelease(xsec_type == "all", "This function is for change all xsecs");

  for (unsigned int ng = 0; ng < n_groups; ng++)
  {
    // sigma t
    sigma_t[ng][mat] = new_xsec[0][ng];
    D[ng][mat] = 1 / (3 * new_xsec[0][ng]);
    // sigma_a
    sigma_r[ng][mat] = new_xsec[1][ng];
    // nu*sigma_f
    double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
    sigma_f[ng][mat] = new_xsec[2][ng];
    nu_sigma_f[ng][mat] = nu * sigma_f[ng][mat];

  }
  // sigma_s
  sigma_s[0][1][mat] = new_xsec[3][0];

  return;
}

/**
 *
 */
void Materials::modify_xsec_7g (std::string xsec_type,
  double sim_time,
  std::vector<double> amplitudes,
  unsigned int n_mat)
{

  double frequency = 1.0;

  if (xsec_type == "capture")
  {
    // Copy xsec

    for (unsigned int g = 0; g < n_groups; ++g)
    {

      sigma_t[g][n_mat] = init_sigma_t[g][n_mat]
                          + amplitudes[g]
                            * sin(2 * M_PI * frequency * sim_time);
      D[g][n_mat] = 1 / (3 * sigma_t[g][n_mat]);

      sigma_r[g][n_mat] = init_sigma_r[g][n_mat]
                          + amplitudes[g]
                            * sin(2 * M_PI * frequency * sim_time);
    }
  }

  return;
}

/**
 *
 */
void Materials::modify_xsec_c5g7_td11 (
  double sim_time)
{

  const unsigned int changing_mat = 7;
  const unsigned int rodded_mat = 8;

  if (sim_time < 1.0)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      sigma_t[g][changing_mat] = init_sigma_t[g][changing_mat]
          + 0.01 * (init_sigma_t[g][rodded_mat] - init_sigma_t[g][changing_mat])
            * sim_time;

      D[g][changing_mat] = 1 / (3 * sigma_t[g][changing_mat]);

      sigma_r[g][changing_mat] = init_sigma_r[g][changing_mat]
          + 0.01 * (init_sigma_r[g][rodded_mat] - init_sigma_r[g][changing_mat])
            * sim_time;

      // TODO Init Velocities?
      velocities_vector[g][changing_mat] = velocities_vector[g][changing_mat]
          + 0.01 * (velocities_vector[g][rodded_mat] - velocities_vector[g][changing_mat])
            * sim_time;

      for (unsigned int g2 = 0; g2 < n_groups; g2++)
      {

        sigma_s[g][g2][changing_mat] = init_sigma_s[g][g2][changing_mat]
            + 0.01 * (init_sigma_s[g][g2][rodded_mat] - init_sigma_s[g][g2][changing_mat])
              * sim_time;
      }
    }
  }
  else if (sim_time < 2.0)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      sigma_t[g][changing_mat] = init_sigma_t[g][changing_mat]
          + 0.01 * (init_sigma_t[g][rodded_mat] - init_sigma_t[g][changing_mat])
            * (2.0 - sim_time);

      D[g][changing_mat] = 1 / (3 * sigma_t[g][changing_mat]);

      sigma_r[g][changing_mat] = init_sigma_r[g][changing_mat]
          + 0.01 * (init_sigma_r[g][rodded_mat] - init_sigma_r[g][changing_mat])
            * (2.0 - sim_time);

      // FIXME initi velocities
      velocities_vector[g][changing_mat] = velocities_vector[g][changing_mat]
          + 0.01 * (velocities_vector[g][rodded_mat] - velocities_vector[g][changing_mat])
            * (2.0 - sim_time);

      for (unsigned int g2 = 0; g2 < n_groups; ++g2)
        sigma_s[g][g2][changing_mat] = init_sigma_s[g][g2][changing_mat]
            + 0.01 * (init_sigma_s[g][g2][rodded_mat] - init_sigma_s[g][g2][changing_mat])
              * (2.0 - sim_time);
    }
  }
  else
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      sigma_t[g][changing_mat] = init_sigma_t[g][changing_mat];
      D[g][changing_mat] = 1 / (3 * sigma_t[g][changing_mat]);

      sigma_r[g][changing_mat] = init_sigma_r[g][changing_mat];

      //velocities_vector[g][changing_mat] = velocities_vector[g][changing_mat];

      for (unsigned int g2 = 0; g2 < n_groups; ++g2)
        sigma_s[g][g2][changing_mat] = init_sigma_s[g][g2][changing_mat];
    }
  }
  return;
}

/**
 *
 */
void Materials::create_new_added_mat (const unsigned int new_mat,
  const double frac,
  double inc_xsec_r,
  const unsigned int mat_no_bar,
  const unsigned int user_index)
{
  Assert(n_groups == 2, ExcMessage("Implemented to n_groups"));
  Assert(frac >= 0.0 and frac <= 1.0,
    ExcMessage("Invalid fraction " + num_to_str(frac)));
//  AssertIndexRange(mat_bar, n_mats);
//  AssertIndexRange(mat_no_bar, n_mats);
//  AssertIndexRange(new_mat, n_mats+1);

// If necessary add a reserve space for a new material
  if (n_mats == new_mat)
  {
    n_mats++;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      D[g].resize(n_mats);
      sigma_r[g].resize(n_mats);
      nu_sigma_f[g].resize(n_mats);
      chi[g].resize(n_mats);
      sigma_f[g].resize(n_mats);
      for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
        sigma_s[g][to_g].resize(n_mats);
    }
  }

  for (unsigned int g = 0; g < n_groups; ++g)
  {
    D[g][new_mat] = D[g][mat_no_bar];
    nu_sigma_f[g][new_mat] = nu_sigma_f[g][mat_no_bar];
    chi[g][new_mat] = chi[g][mat_no_bar];
    sigma_f[g][new_mat] = sigma_f[g][mat_no_bar];
    for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
      sigma_s[g][to_g][new_mat] = sigma_s[g][to_g][mat_no_bar];
  }

  sigma_r[0][new_mat] = sigma_r[0][mat_no_bar];
  sigma_r[1][new_mat] = sigma_r[1][mat_no_bar] + frac * inc_xsec_r;

  // Create new position in materials_vector
  AssertIndexRange(user_index, materials_vector.size());
  materials_vector[user_index] = new_mat;
}

/**
 *
 */
void Materials::add_perturbation_xsec (std::string xsec_type,
  unsigned int mat,
  std::vector<double> &delta_xsec)
{

  //
  for (unsigned int ng = 0; ng < n_groups; ng++)
  {
    if (xsec_type == "sigma_f")
    {
      double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
      nu_sigma_f[ng][mat] += delta_xsec[ng];
      sigma_f[ng][mat] += delta_xsec[ng] / nu;
    }
    else if (xsec_type == "sigma_a")
    {
      sigma_r[ng][mat] += delta_xsec[ng];
    }
    else
    {
      AssertRelease(false, "Invalid type of xsec");
    }

  }

  return;
}

/**
 *
 */
void Materials::add_perturbation_xsec (std::string xsec_type,
  std::vector<double> &delta_xsec)
{
  // FIXME nu is assumed equal to 1

  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    for (unsigned int ng = 0; ng < n_groups; ng++)
    {
      if (xsec_type == "sigma_f")
      {
        double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
        nu_sigma_f[ng][mat] += delta_xsec[ng];
        sigma_f[ng][mat] += delta_xsec[ng] / nu;
      }
      else if (xsec_type == "sigma_a")
      {
        sigma_r[ng][mat] += delta_xsec[ng];
      }
      else
      {
        AssertRelease(false, "Invalid type of xsec");
      }

    }
  }

  return;
}

/**
 *
 */
void Materials::add_perturbation_xsec (std::string xsec_type,
  double &coeff)
{
  double nu;
  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    for (unsigned int ng = 0; ng < n_groups; ng++)
    {
      if (xsec_type == "sigma_f")
      {
        nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
        nu_sigma_f[ng][mat] += coeff;
        sigma_f[ng][mat] += coeff / nu;
      }
      else if (xsec_type == "sigma_a")
      {
        sigma_r[ng][mat] += coeff;
      }
      else
      {
        AssertRelease(false, "Invalid type of xsec");
      }

    }
  }

  return;
}

/**
 *
 */
void Materials::remove_perturbation_xsec (std::string xsec_type,
  unsigned int mat,
  std::vector<double> &delta_xsec)
{

  // FIXME nu is assumed equal to 1
  for (unsigned int ng = 0; ng < n_groups; ng++)
  {
    if (xsec_type == "sigma_f")
    {
      double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
      nu_sigma_f[ng][mat] -= delta_xsec[ng];
      sigma_f[ng][mat] -= delta_xsec[ng] / nu;
    }
    else if (xsec_type == "sigma_a")
    {
      sigma_r[ng][mat] -= delta_xsec[ng];
    }
    else
    {
      AssertRelease(false, "Invalid type of xsec");
    }

  }

  return;
}

/**
 *
 */
void Materials::remove_perturbation_xsec (std::string xsec_type,
  std::vector<double> &delta_xsec)
{

  // FIXME nu is assumed equal to 1

  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    for (unsigned int ng = 0; ng < n_groups; ng++)
    {
      if (xsec_type == "sigma_f")
      {
        double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
        nu_sigma_f[ng][mat] -= delta_xsec[ng];
        sigma_f[ng][mat] -= delta_xsec[ng] / nu;
      }
      else if (xsec_type == "sigma_a")
      {
        sigma_r[ng][mat] -= delta_xsec[ng];
      }
      else
      {
        AssertRelease(false, "Invalid type of xsec");
      }

    }
  }

  return;
}

/**
 *
 */
void Materials::remove_perturbation_xsec (std::string xsec_type,
  double &coeff)
{

  // FIXME nu is assumed equal to 1
  for (unsigned int mat = 0; mat < n_mats; mat++)
  {
    for (unsigned int ng = 0; ng < n_groups; ng++)
    {
      if (xsec_type == "sigma_f")
      {
        double nu = nu_sigma_f[ng][mat] / sigma_f[ng][mat];
        nu_sigma_f[ng][mat] -= coeff;
        sigma_f[ng][mat] -= coeff / nu;
      }
      else if (xsec_type == "sigma_a")
      {
        sigma_r[ng][mat] -= coeff;
      }
      else
      {
        AssertRelease(false, "Invalid type of xsec");
      }

    }
  }

  return;
}

/**
 *
 */
void Materials::save_initial_xsec ()
{

  // Resize xsec
  init_D.resize(n_groups, std::vector<double>(n_mats));
  init_sigma_t.resize(n_groups, std::vector<double>(n_mats));
  init_sigma_r.resize(n_groups, std::vector<double>(n_mats));
  init_nu_sigma_f.resize(n_groups, std::vector<double>(n_mats));
  init_chi.resize(n_groups, std::vector<double>(n_mats));
  init_sigma_f.resize(n_groups, std::vector<double>(n_mats));
  init_sigma_s.resize(n_groups,
    std::vector<std::vector<double> >(n_groups,
      std::vector<double>(n_mats)));

  // Copy xsec
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    init_D[g] = D[g];
    init_sigma_t[g] = sigma_t[g];
    init_sigma_r[g] = sigma_r[g];
    init_nu_sigma_f[g] = nu_sigma_f[g];
    init_chi[g] = chi[g];
    init_sigma_f[g] = sigma_f[g];
    for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
      init_sigma_s[g][to_g] = sigma_s[g][to_g];
  }

  materials_vector_init = materials_vector;
}

/**
 *
 */
void Materials::save_n_mats_init ()
{

  // Resize xsec
  n_mats_init = n_mats;

}

/**
 *
 */
void Materials::compute_xsec_perturbation (unsigned int ncell,
  std::vector<std::vector<std::vector<double>>> &xsec_pert,
  std::vector<std::vector<std::vector<double>>> &xsec_pert_s,
  std::vector<unsigned int> &materials_vector_pert)
{

  unsigned int mat_id = materials_vector[ncell];
  materials_vector_pert.push_back(ncell);

  xsec_pert.resize(4, std::vector<std::vector<double> >(n_groups));

  xsec_pert_s.resize(n_groups, std::vector<std::vector<double> >(n_groups));

  for (unsigned int ng = 0; ng < n_groups; ng++)
  {

    // 0: diffusion, 1: sigma_r, 2: nusigmaf, 3: sigma_f,
    xsec_pert[0][ng].push_back(
      D[ng][mat_id] - init_D[ng][materials_vector_init[ncell]]);

    xsec_pert[1][ng].push_back(
      sigma_r[ng][mat_id]
      - init_sigma_r[ng][materials_vector_init[ncell]]);

    xsec_pert[2][ng].push_back(
      nu_sigma_f[ng][mat_id]
      - init_nu_sigma_f[ng][materials_vector_init[ncell]]);
    xsec_pert[3][ng].push_back(
      sigma_f[ng][mat_id]
      - init_sigma_f[ng][materials_vector_init[ncell]]);

    for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
    {
      xsec_pert_s[ng][to_g].push_back(
        sigma_s[ng][to_g][mat_id]
        - init_sigma_s[ng][to_g][materials_vector_init[ncell]]);
    }

  }

}

/**
 *
 */
void Materials::parse_xsec_2g (const std::string &xs_file,
  const std::vector<unsigned int> &n_assemblies_per_dim,
  const unsigned int n_assemblies)
{
  Assert(fexists(xs_file), ExcMessage("xs_file doesn't exist"));
  std::ifstream input(xs_file.c_str(), std::ios::in);
  std::string sub, keyword;
  unsigned int mat;
  double num;


  // for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    keyword.clear();
    iss >> keyword;

    if (is_commentary(keyword))
      continue;
    // First definition Material and XSecs:
    else if (keyword == "Materials")
    {
      verbose_cout << "  Materials " << std::flush;

      materials_vector.reserve(n_assemblies); // Reserve space for allocation
      parse_multiline_vector(input,
        n_assemblies_per_dim[1] * n_assemblies_per_dim[2],
        materials_vector, true);

      // Subtract 1 because we use vectors of C
      for (unsigned int i = 0; i < n_assemblies; i++)
        materials_vector[i]--;

      verbose_cout << " Done! " << std::endl;
    }
    else if (keyword == "XSecs")
    {
      verbose_cout << "  XSecs..." << std::endl;

      std::string str;
      iss >> n_mats;
      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of materials defined!"));
      Assert(n_mats>0,
        ExcMessage("It must be defined the number of materials defined!"));

      verbose_cout << "  n_mats " << n_mats << std::endl;
      // Resize
      sigma_t.resize(n_groups, std::vector<double>(n_mats));
      sigma_r.resize(n_groups, std::vector<double>(n_mats));
      nu_sigma_f.resize(n_groups, std::vector<double>(n_mats));
      chi.resize(n_groups, std::vector<double>(n_mats));
      sigma_f.resize(n_groups, std::vector<double>(n_mats));
      sigma_s.resize(n_groups,
        std::vector<std::vector<double> >(n_groups,
          std::vector<double>(n_mats)));

      for (unsigned int j = 0; j < n_mats; j++)
      {
        str = get_new_valid_line(input, line);
        std::istringstream iss(line);
        iss >> str;
        mat = str_to_num<unsigned int>(str) - 1;

        AssertRelease(mat == j,
          "Error in the XS in line " + num_to_str(j + 1));
        verbose_cout << "    mat " << mat + 1 << std::endl;

        // sigma_tr1
        iss >> num;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        sigma_t[0][mat] = num;

        // sigma_a1
        iss >> num;
        double sigma_a1 = num;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));

        // nu_sigma_f1
        iss >> num;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        nu_sigma_f[0][mat] = num;

        // sigma_f1
        iss >> num;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        sigma_f[0][mat] = num;

        // sigma_12
        iss >> num;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        sigma_s[0][1][mat] = num;
        sigma_r[0][mat] = sigma_a1 + num;

        // new line
        get_new_valid_line(input, line);
        std::istringstream iss2(line);

        // sigma_tr2
        iss2 >> num;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        sigma_t[1][mat] = num;

        // sigma_a2
        iss2 >> num;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        sigma_r[1][mat] = num;

        //  nu_sigma_f2
        iss2 >> num;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        nu_sigma_f[1][mat] = num;

        // sigma_f2
        iss2 >> num;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        sigma_f[1][mat] = num;

        // Possible upscattering sigma_21
        if (!iss2.eof())
        {
          iss2 >> num;
          sigma_s[1][0][mat] = num;
          sigma_r[1][mat] += num;
        }

        // chi
        chi[0][mat] = 1.0;
        chi[1][mat] = 0.0;
      }
      verbose_cout << " Done!" << std::endl;
    }
    // Neutron precursors and velocities
    else if (keyword == "Precursors")
    {

      verbose_cout << "WARNING!" << std::endl;
      verbose_cout << "Precursors are not needed" << std::endl;

      std::string str;
      unsigned int prec_name;
      iss >> n_precursors;
      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of materials defined!"));
      Assert(n_precursors>0,
        ExcMessage("It must be defined the number of precursors defined!"));
//      beta_prec.resize(n_precursors);
//      lambda_prec.resize(n_precursors);

      delayed_fractions.resize(n_mats, std::vector<double>(n_precursors));
      delayed_decay_constants.resize(n_mats,
        std::vector<double>(n_precursors));

      for (unsigned int j = 0; j < n_precursors; j++)
      {
        getline(input, line);
        std::istringstream iss(line);

        iss >> str;
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) precursors specified"));
        trim(str);
        // Allow commentaries
        if ((str == "#") or (str == "!") or (str == "//")
            or str == "")
        {
          j--;
          continue;
        }

        prec_name = Utilities::string_to_int(str) - 1;
        AssertRelease(prec_name == j,
          "Precursors should be defined in order!");
        Assert(n_mats>0,
          ExcMessage("It must be defined the number of materials defined!"));

        iss >> num;
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) precursors specified"));
//        beta_prec[prec_name] = num;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          delayed_fractions[mat][prec_name] = num;

        iss >> num;
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) precursors specified"));
//        lambda_prec[prec_name] = num;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          delayed_decay_constants[mat][prec_name] = num;
      }

//      beta_total = sum_vector(beta_prec);

      delayed_fraction_sum.resize(n_mats);
      for (unsigned int mat = 0; mat < n_mats; mat++)
        for (unsigned int p = 0; p < n_precursors; p++)
          delayed_fraction_sum[mat] += delayed_fractions[mat][p];

      verbose_cout << "Precursors... ok" << std::endl;

    }
    else if (keyword == "Velocity")
    {
      verbose_cout << "  Parsing velocities..." << std::flush;
//      velocities.resize(n_groups);
      velocities_vector.resize(n_mats, std::vector<double>(n_groups));
      getline(input, line);
      std::istringstream iss(line);
      iss >> num;
//      velocities[0] = num;
      for (unsigned int mat = 0; mat < n_mats; mat++)
              velocities_vector[mat][0] = num;

      Assert(! iss.fail(), ExcMessage("Velocity v1 not (well) defined"));
      iss >> num;
//      velocities[1] = num;
      for (unsigned int mat = 0; mat < n_mats; mat++)
                   velocities_vector[mat][1] = num;

      Assert(! iss.fail(), ExcMessage("Velocity v2 not (well) defined"));
      verbose_cout << " Done!" << std::endl;
    }
    else
      // Error
      Assert(false, ExcMessage("Invalid Header in XSEC_file: " + keyword));
  }

  if (!velocities_vector.empty())
  {
    delayed_spectra.resize(n_mats,
      std::vector<std::vector<double> >(n_precursors,
        std::vector<double>(n_groups)));

    for (unsigned int mat = 0; mat < n_mats; mat++)
      for (unsigned int p = 0; p < n_precursors; p++)
      {
        delayed_spectra[mat][p][0] = 1.0;
        delayed_spectra[mat][p][1] = 0.0;
      }

    prompt_spectra.resize(n_mats, std::vector<double>(n_groups));
    for (unsigned int mat = 0; mat < n_mats; mat++)
      for (unsigned int g = 0; g < n_groups; g++)
      {

        double sum_del_beta = 0.0;
        for (unsigned int p = 0; p < n_precursors; p++)
        {
          sum_del_beta += delayed_spectra[mat][p][g]
                          * delayed_fractions[mat][p];
        }

        prompt_spectra[mat][g] = (chi[g][mat] - sum_del_beta)
                                 / (1 - delayed_fraction_sum[mat]);
      }
  }

  std::string precursors_file(xs_file.begin(), xs_file.end() - 4);
  precursors_file += "prec";

  if (velocities_vector.empty() and transient)
    parse_precursors_file(precursors_file);

}

/**
 * parse_xs_and_dfs(std::string xsec_file)
 *  This Functions parses the XSEC file completing
 *  D1v, D2v, sigma_a1v,sigma_a2v, nu_sigma_f1v,  nu_sigma_f2v, nu_sigma_f1v.
 *  and materials_table
 */
void Materials::parse_xsec (const std::string &xsec_file,
  const std::vector<unsigned int> &n_assemblies_per_dim,
  const unsigned int n_assemblies)
{

  AssertRelease(fexists(xsec_file), "xsec_file doesn't exist");
  std::ifstream input(xsec_file.c_str(), std::ios::in);
  std::string sub, word, str;
  std::string assembly_file;
  double num;

// for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    word.clear();
    iss >> word;
    trim(word);

    if (is_commentary(word))
      continue;

    // First definition Material and XSecs:
    else if (word == "Materials")
    {
      materials_vector.reserve(n_assemblies); // Reserve space for allocation
      parse_multiline_vector(input,
        n_assemblies_per_dim[1] * n_assemblies_per_dim[2],
        materials_vector, true);

      // Subtract 1 because we use C-style vectors
      for (unsigned int i = 0; i < n_assemblies; i++)
        materials_vector[i] -= 1;
      verbose_cout << "materials_vector: " << std::flush;
      if (verbose_cout.is_active())
        print_vector(materials_vector);
    }
    else if (word == "XSecs" or word == "XSs")
    {
      verbose_cout << "parsing XSecs..." << std::endl;

      iss >> n_mats;
      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of materials defined!"));
      Assert(n_mats>0,
        ExcMessage("It must be defined the number of materials defined!"));

      verbose_cout << " n_mats " << n_mats << std::endl;
      // Resize cross sections
      sigma_s.resize(n_groups,
        std::vector<std::vector<double> >(n_groups,
          std::vector<double>(n_mats)));
      sigma_t.resize(n_groups, std::vector<double>(n_mats));
      chi.resize(n_groups, std::vector<double>(n_mats));
      nu_sigma_f.resize(n_groups, std::vector<double>(n_mats));
      sigma_f.resize(n_groups, std::vector<double>(n_mats));
      sigma_r.resize(n_groups, std::vector<double>(n_mats));

      for (unsigned int mat = 0; mat < n_mats; mat++)
      {
        verbose_cout << " Material " << mat + 1 << std::endl;
        for (unsigned int g = 0; g < n_groups; g++)
        {
          // Get index of the Material defined
          get_new_valid_line(input, line);

          std::istringstream iss(line);
          Assert(! iss.fail(),
            ExcMessage("There are not enough (well) XSEC specified"));

          // Get the redundant material number
          if (g == 0)
          {
            iss >> str;
            trim(str);
            Assert(mat == (str_to_num<unsigned int>(str) - 1),
              ExcMessage("Materials should be defined in order!"));
          }
          //verbose_cout << "   Group " << g + 1 << std::endl;

          // sigma_t
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) XSEC specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(g + 1)
            + ".");
          sigma_t[g][mat] = num;
          verbose_cout << "      sigma_t_" << g + 1 << " = " << num
                       << std::endl;

          // chi
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) XSEC specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(g + 1)
            + ".");
          chi[g][mat] = num;
          verbose_cout << "      chi_" << g + 1 << " = " << num
                       << std::endl;

          // nu_sigma_f
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) XSEC specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(g + 1)
            + ".");
          nu_sigma_f[g][mat] = num;

          verbose_cout << "      nu_sigma_fv_" << g + 1 << " = "
                       << num
                       << std::endl;

          // sigma_f
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) XSEC specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(g + 1)
            + ".");
          sigma_f[g][mat] = num;
          verbose_cout << "      sigma_fv_" << g + 1 << " = " << num
                       << std::endl;

          for (unsigned int to_group = 0; to_group < n_groups;
              to_group++)
          {
            iss >> num;
            AssertRelease(!iss.fail(),
              "There are not enough (well) XSEC specified in mat "
              + num_to_str(mat + 1)
              + " group "
              + num_to_str(g + 1)
              + ".");
            // Be careful because we define sigma_s negative!
            sigma_s[g][to_group][mat] = num;
            verbose_cout << "      sigma_s" << g + 1 << "->"
                         << to_group + 1
                         << " = " << num << std::endl;

          }

          // sigma_r
          sigma_r[g][mat] = sigma_t[g][mat] - sigma_s[g][g][mat];
        }

      }

    }
    else if (word == "Precursors")
    {

      verbose_cout << "Precursors... " << std::endl;
      std::string str;
      unsigned int prec_name;
      iss >> n_precursors;
      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of materials defined!"));
//      beta_prec.resize(n_precursors);
//      lambda_prec.resize(n_precursors);

      delayed_fractions.resize(n_mats, std::vector<double>(n_precursors));
      delayed_decay_constants.resize(n_mats,
        std::vector<double>(n_precursors));

      for (unsigned int j = 0; j < n_precursors; j++)
      {
        getline(input, line);
        std::istringstream iss(line);

        iss >> str;
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) precursors specified"));
        trim(str);
        // Allow commentaries
        if ((str == "#") or (str == "!") or (str == "//")
            or str == "")
        {
          j--;
          continue;
        }

        prec_name = Utilities::string_to_int(str) - 1;
        AssertRelease(prec_name == j,
          "Precursors should be defined in order!");
        Assert(n_mats>0,
          ExcMessage("It must be defined the number of materials defined!"));

        iss >> num;

        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) precursors specified"));
//        beta_prec[prec_name] = num;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          delayed_fractions[mat][prec_name] = num;

        iss >> num;
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) precursors specified"));
//        lambda_prec[prec_name] = num;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          delayed_decay_constants[mat][prec_name] = num;
      }

//      beta_total = sum_vector(beta_prec);

      delayed_fraction_sum.resize(n_mats);
      for (unsigned int mat = 0; mat < n_mats; mat++)
        for (unsigned int p = 0; p < n_precursors; p++)
          delayed_fraction_sum[mat] += delayed_fractions[mat][p];

      verbose_cout << "Precursors... ok" << std::endl;

    }
    else if (word == "Velocity")
    {
//      velocities.resize(n_groups);
      velocities_vector.resize(n_mats, std::vector<double>(n_groups));
      getline(input, line);
      std::istringstream iss(line);
      for (unsigned int g = 0; g < n_groups; g++)
      {
        iss >> num;
//        velocities[g] = num;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          velocities_vector[mat][g] = num;
        Assert(! iss.fail(),
          ExcMessage("Velocity v not (well) defined"));
      }

      verbose_cout << "Velocity... ok" << std::endl;

    }
    else if (word == "Decay_Spectrum")
    {
      delayed_spectra.resize(n_mats,
        std::vector<std::vector<double> >(n_precursors,
          std::vector<double>(n_groups)));
      getline(input, line);
      std::istringstream iss(line);
      for (unsigned int g = 0; g < n_groups; g++)
      {
        iss >> num;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          for (unsigned int p = 0; p < n_precursors; p++)
            delayed_spectra[mat][p][g] = num;
        Assert(! iss.fail(),
          ExcMessage("Decay_Spectrum not (well) defined"));
      }

      verbose_cout << " Decay_Spectrum... Done!" << std::endl;

    }
    else
      // Error
      Assert(false, ExcMessage("Invalid Header in XSEC_file: " + word))

  }

// Compute sigma_r
  sigma_r.resize(n_groups, std::vector<double>(n_mats));
  for (unsigned int mat = 0; mat < n_mats; mat++)
    for (unsigned int g = 0; g < n_groups; g++)
    {
      sigma_r[g][mat] = sigma_t[g][mat] - sigma_s[g][g][mat];
    }
  verbose_cout << "Done!" << std::endl;

  if (!velocities_vector.empty())
  {

    if (delayed_spectra.size() < 1)
    {
      // Create the delayed spectra vector
      delayed_spectra.resize(n_mats,
        std::vector<std::vector<double> >(n_precursors,
          std::vector<double>(n_groups)));

      for (unsigned int mat = 0; mat < n_mats; mat++)
        for (unsigned int p = 0; p < n_precursors; p++)
          for (unsigned int g = 0; g < n_groups; g++)
          {
            delayed_spectra[mat][p][g] = chi[g][mat];

          }
    }

    prompt_spectra.resize(n_mats, std::vector<double>(n_groups));
    for (unsigned int mat = 0; mat < n_mats; mat++)
      for (unsigned int g = 0; g < n_groups; g++)
      {

        double sum_del_beta = 0.0;
        for (unsigned int p = 0; p < n_precursors; p++)
        {
          sum_del_beta += delayed_spectra[mat][p][g]
                          * delayed_fractions[mat][p];
        }

        prompt_spectra[mat][g] = (chi[g][mat] - sum_del_beta)
                                 / (1 - delayed_fraction_sum[mat]);
      }

  }

  std::string precursors_file(xsec_file.begin(), xsec_file.end() - 4);
  precursors_file += "prec";

  if (velocities_vector.empty() and transient)
    parse_precursors_file(precursors_file);

}

/**
 * parse_precursors_file(std::string prec)
 *  This Functions parses the Precursors file
 */
void Materials::parse_precursors_file (const std::string &prec_file)
{

  AssertRelease(fexists(prec_file), "prec_file doesn't exist");
  std::ifstream input(prec_file.c_str(), std::ios::in);
  std::string sub, word, str;
  std::string assembly_file;
  double num;
  unsigned int numint;

// for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    word.clear();
    iss >> word;
    trim(word);

    if (is_commentary(word))
      continue;

    // First definition Material and XSecs:
    else if (word == "Velocities")
    {
      velocities_vector.resize(n_mats, std::vector<double>(n_groups));

      for (unsigned int g = 0; g < n_groups; g++)
      {
        // Get index of the Material defined
        get_new_valid_line(input, line);

        std::istringstream iss(line);
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) XSEC specified"));

        // Get the redundant energy group number
        iss >> str;
        trim(str);

        verbose_cout << "   Group " << g + 1 << std::endl;

        for (unsigned int mat = 0; mat < n_mats; mat++)
        {

          verbose_cout << " Material " << mat + 1 << std::endl;

          // velocities
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) XSEC specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(g + 1)
            + ".");
          velocities_vector[mat][g] = num;

        }

        verbose_cout << "velocities_vector: " << std::flush;
        for (unsigned int mat = 0; mat < n_mats; mat++)
          if (verbose_cout.is_active())
            verbose_cout << velocities_vector[mat][g] << " "
                         << std::flush;
      }

      verbose_cout << "  " << std::endl;

    }
    else if (word == "N_PrecMat")
    {
      verbose_cout << "parsing N_PrecMat..." << std::endl;

      iss >> n_prec_mat;
      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of materials with precursors!"));
      Assert(n_prec_mat>0,
        ExcMessage("The number of materials with precursors must be greater than 0"));

      verbose_cout << " n_prec_mat " << n_prec_mat << std::endl;
      precursors_materials.resize(n_prec_mat);

    }
    else if (word == "PrecMat")
    {
      verbose_cout << "parsing PrecMat..." << std::endl;

      for (unsigned int nm = 0; nm < n_prec_mat; nm++)
      {
        iss >> numint;

        Assert(!iss.fail(),
          ExcMessage("The number of material with precursors is not correct"));
        Assert(numint>0 and numint<n_mats+1,
          ExcMessage("The number of material with precursors is not correct"));

        precursors_materials[nm] = numint - 1;

      }

      if (verbose_cout.is_active())
        print_vector(precursors_materials);

    }
    else if (word == "N_Delayed_Groups")
    {
      verbose_cout << "   N_Delayed_Groups..." << std::flush;

      iss >> n_precursors;
      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of delayed groups!"));
      Assert(n_precursors>0,
        ExcMessage("The number of delayed groups must be greater than 0"));

      verbose_cout << n_precursors << std::endl;

    }
    else if (word == "Delayed_Fractions")
    {
      verbose_cout << "   Delayed_Fractions" << std::endl;

      delayed_fractions.resize(n_mats, std::vector<double>(n_precursors, 0.0));
      for (unsigned int mat = 0; mat < n_mats; mat++)
        delayed_fractions[mat].resize(n_precursors);

      for (unsigned int p = 0; p < n_precursors; p++)
      {
        // Get index of the Material defined
        get_new_valid_line(input, line);

        std::istringstream iss(line);
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) Delayed_Fractions specified"));

        // Get the redundant energy group number
        iss >> str;
        trim(str);

        for (unsigned int mat = 0; mat < n_prec_mat; mat++)
        {
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) Delayed_Fractions specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(p + 1)
            + ".");
          delayed_fractions[precursors_materials[mat]][p] = num;

        }

        verbose_cout << "delayed_fractions: " << std::flush;
        for (unsigned int mat = 0; mat < n_prec_mat; mat++)
          if (verbose_cout.is_active())
            verbose_cout
            << delayed_fractions[precursors_materials[mat]][p]
            << " "
            << std::flush;
      }
      verbose_cout << "  " << std::endl;

    }
    else if (word == "Delayed_Decay_Constants")
    {
      verbose_cout << "   Delayed_Decay_Constants... " << std::endl;
      // Resize
      delayed_decay_constants.resize(n_mats, std::vector<double>(n_precursors, 0.0));
      for (unsigned int mat = 0; mat < n_mats; mat++)
        delayed_decay_constants[mat].resize(n_precursors);

      for (unsigned int p = 0; p < n_precursors; p++)
      {
        // Get index of the Material defined
        get_new_valid_line(input, line);

        std::istringstream iss(line);
        Assert(! iss.fail(),
          ExcMessage("There are not enough (well) delayed_decay_constants specified"));

        // Get the redundant energy group number
        iss >> str;
        trim(str);

        for (unsigned int mat = 0; mat < n_prec_mat; mat++)
        {
          iss >> num;
          AssertRelease(!iss.fail(),
            "There are not enough (well) delayed_decay_constants specified in mat "
            + num_to_str(mat + 1)
            + " group "
            + num_to_str(p + 1)
            + ".");
          delayed_decay_constants[precursors_materials[mat]][p] = num;
        }

        verbose_cout << "delayed_decay_constants: " << std::flush;
        for (unsigned int mat = 0; mat < n_prec_mat; mat++)
          if (verbose_cout.is_active())
            verbose_cout
            << delayed_decay_constants[precursors_materials[mat]][p]
            << " "
            << std::flush;
      }

      verbose_cout << "  " << std::endl;

    }
    else if (word == "Delayed_Spectra")
    {

      delayed_spectra.resize(n_mats);

      for (unsigned int mat = 0; mat < n_mats; mat++)
        delayed_spectra[mat].resize(n_precursors,
          std::vector<double>(n_groups));

      for (unsigned int mat = 0; mat < n_prec_mat; mat++)
      {

        // Get index of the Material defined
        get_new_valid_line(input, line);

        std::istringstream iss(line);

        word.clear();
        iss >> word;
        trim(word);

        iss >> numint;

        AssertRelease(numint - 1 == precursors_materials[mat],
          "The materials must be equal to precursors materials");

//        delayed_spectra[precursors_materials[mat]].resize(n_precursors,
//          std::vector<double>(n_groups));

        for (unsigned int g = 0; g < n_groups; g++)
        {
          // Get index of the Material defined
          get_new_valid_line(input, line);

          std::istringstream iss(line);
          Assert(! iss.fail(),
            ExcMessage("There are not enough (well) Delayed_Spectra specified"));

          // Get the redundant energy group number
          iss >> str;
          trim(str);

          verbose_cout << "   Energy group " << g + 1 << std::endl;

          for (unsigned int p = 0; p < n_precursors; p++)
          {

            verbose_cout << " Material " << mat + 1 << std::endl;

            // velocities
            iss >> num;
            AssertRelease(!iss.fail(),
              "There are not enough (well) Delayed_Spectra specified in mat "
              + num_to_str(g + 1)
              + " group "
              + num_to_str(p + 1)
              + ".");
            delayed_spectra[precursors_materials[mat]][p][g] = num;
            verbose_cout << "      sigma_t_" << p + 1 << " = "
                         << num
                         << std::endl;

          }

          verbose_cout << "Delayed_Spectra: " << std::flush;
          for (unsigned int p = 0; p < n_precursors; p++)
            if (verbose_cout.is_active())
              verbose_cout
              << delayed_spectra[precursors_materials[mat]][p][g]
              << " "
              << std::flush;

          verbose_cout << "  " << std::endl;
        }

      }

    }
    else
      // Error
      Assert(false, ExcMessage("Invalid Header in XSEC_file: " + word))
  }

  delayed_fraction_sum.resize(n_mats);
  for (unsigned int mat = 0; mat < n_mats; mat++)
    for (unsigned int p = 0; p < n_precursors; p++)
      delayed_fraction_sum[mat] += delayed_fractions[mat][p];

  prompt_spectra.resize(n_mats, std::vector<double>(n_groups));
  for (unsigned int mat = 0; mat < n_mats; mat++)
    for (unsigned int g = 0; g < n_groups; g++)
    {

      double sum_del_beta = 0.0;
      for (unsigned int p = 0; p < n_precursors; p++)
      {
        sum_del_beta += delayed_spectra[mat][p][g]
                        * delayed_fractions[mat][p];
      }

      prompt_spectra[mat][g] = (chi[g][mat] - sum_del_beta)
                               / (1 - delayed_fraction_sum[mat]);
    }

  verbose_cout << "Done!" << std::endl;
}

/**
 *  @brief It parses the XS.xml file with a XML format.
 */
void Materials::parse_forest_xs (const std::string &xml_file)
{
  verbose_cout << "parse_forest_xs...  " << xml_file << std::endl;
  AssertRelease(fexists(xml_file), "forest_file doesn't exist");

  XMLInput::InputMat input;
  input.load(xml_file);
  AssertRelease(input.get_n_groups() == n_groups,
    "n_groups in xml file does not match " + num_to_str(n_groups)
    + " vs "
    + num_to_str(input.get_n_groups()));

  // Resize containers
  n_mats = input.get_n_mat();
  verbose_cout << "n_mats: " << n_mats << std::endl;

  // Resize Containers
  sigma_t.resize(n_groups, std::vector<double>(n_mats));
  sigma_s.resize(n_groups,
    std::vector<std::vector<double> >(n_groups,
      std::vector<double>(n_mats)));
  chi.resize(n_groups, std::vector<double>(n_mats));
  sigma_r.resize(n_groups, std::vector<double>(n_mats));
  nu_sigma_f.resize(n_groups, std::vector<double>(n_mats));
  sigma_f.resize(n_groups, std::vector<double>(n_mats));

  // Resize time material data
  n_precursors = input.get_n_precursors();
  //velocities.resize(n_groups);
  velocities_vector.resize(n_mats, std::vector<double>(n_groups));
  delayed_fractions.resize(n_mats, std::vector<double>(n_precursors));
  delayed_fraction_sum.resize(n_mats);
  delayed_decay_constants.resize(n_mats,
    std::vector<double>(n_precursors));
  prompt_spectra.resize(n_mats, std::vector<double>(n_groups));
  delayed_spectra.resize(n_mats,
    std::vector<std::vector<double> >(n_precursors,
      std::vector<double>(n_groups)));

  // Fill materials Vector
  materials_vector = input.get_materials_vector();
  if (materials_vector.empty())
    for (unsigned int mat = 0; mat < n_mats; ++mat)
      materials_vector.push_back(mat);

  for (unsigned int mat = 0; mat < n_mats; ++mat)
  {
    verbose_cout << "  Material " << mat << std::endl;
    AssertRelease(input.xs[mat].id == mat, "Error in mat ids");
    for (unsigned int from_g = 0; from_g < n_groups; ++from_g)
    {
      verbose_cout << "    Group " << from_g + 1 << std::endl;
      AssertRelease(input.xs[mat].exist_sigma_t,
        "Sigma_t does not exist");
      AssertRelease(input.xs[mat].exist_sigma_s,
        "Sigma_s does not exist");
      AssertRelease(input.xs[mat].exist_chi, "Chi does not exist");
      AssertRelease(input.xs[mat].exist_nu_sigma_f,
        "Nu Sigma_f does not exist");

      // sigma_t
      sigma_t[from_g][mat] = input.xs[mat].sigma_tr[from_g]; // TODO
      verbose_cout << "        sigma_t_" << from_g + 1 << " = "
                   << sigma_t[from_g][mat]
                   << std::endl;

      // chi
      chi[from_g][mat] = input.xs[mat].chi[from_g];
      verbose_cout << "        chi_" << from_g + 1 << " = "
                   << chi[from_g][mat]
                   << std::endl;

      // nu_sigma_f
      nu_sigma_f[from_g][mat] = input.xs[mat].nu_sigma_f[from_g];
      verbose_cout << "        nu_sigma_fv_" << from_g + 1 << " = "
                   << nu_sigma_f[from_g][mat]
                   << std::endl;

      // sigma_r
      sigma_r[from_g][mat] = input.xs[mat].sigma_r[from_g];
      verbose_cout << "        nu_sigma_fv_" << from_g + 1 << " = "
                   << nu_sigma_f[from_g][mat]
                   << std::endl;

      // sigma_f
      // If sigma_f exists get it if not use nusigf instead
      if (!input.xs[mat].exist_sigma_f)
      {
        input.xs[mat].sigma_f = input.xs[mat].nu_sigma_f;
        input.xs[mat].exist_sigma_f = true;
      }
      sigma_f[from_g][mat] = input.xs[mat].sigma_f[from_g];
      verbose_cout << "        sigma_f_" << from_g + 1 << " = "
                   << sigma_f[from_g][mat]
                   << std::endl;

      // sigma_s
      for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
      {
        // Be careful because in input.xs sigma_s is in a different way
        // Also be careful because we define sigma_s negative!
        sigma_s[from_g][to_g][mat] =
                                     input.xs[mat].sigma_s[to_g][from_g];

        verbose_cout << "        sigma_s_intergroup_" << from_g + 1
                     << "->"
                     << to_g + 1 << " = "
                     << sigma_s[from_g][to_g][mat]
                     << std::endl;
      }

      // velocities
      if (input.xs[mat].exist_velocities)
      {
        verbose_cout << "       velocities " << from_g + 1 << " = " << std::flush;
        velocities_vector[mat][from_g] = input.xs[mat].velocities[from_g];
        verbose_cout << velocities_vector[mat][from_g] << std::endl;
      }
      // prompt spectra - chi_p
      if (input.xs[mat].exist_chi_d)
      {
        verbose_cout << "       prompt_spectra " << from_g + 1 << " = " << std::flush;
        prompt_spectra[mat][from_g] = input.xs[mat].chi_p[from_g];
      }
    }

    double sum = 0;
    for (unsigned int to_g = 0; to_g < n_groups; ++to_g){
          sum += chi[to_g][mat];
    }

    // Transient data
    for (unsigned int p = 0; p < n_precursors; ++p)
    {
      // lambda
      if (input.xs[mat].exist_lambda)
      {
        verbose_cout << "       delayed_decay_constants " << p + 1 << " = " << std::flush;
        delayed_decay_constants[mat][p] = input.xs[mat].lambda[p];
        verbose_cout << delayed_decay_constants[mat][p] << std::endl;
      }
      // beta
      if (input.xs[mat].exist_beta)
      {
        verbose_cout << "       delayed_fractions " << p + 1 << " = " << std::flush;
        delayed_fractions[mat][p] = input.xs[mat].beta[p];
        verbose_cout << delayed_fractions[mat][p] << std::endl;
      }

      // prompt spectra
      if (input.xs[mat].exist_chi_d)
      {
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          verbose_cout << "       delayed_spectra (g=" << g + 1 << ") p=" << p + 1
                       << ")  = "
                       << std::flush;
          delayed_spectra[mat][p][g] = input.xs[mat].chi_d[p][g];
          verbose_cout << delayed_spectra[mat][p][g] << std::endl;
        }
      }
    }

    delayed_fraction_sum[mat] = input.xs[mat].beta_eff;

    // Check spectra relation (1-beta)*chip+sum_k betak*chidk =chi
    for (unsigned int g = 0; g < n_groups; ++g){
    double chi_computed=0;
    for (unsigned int p = 0; p < n_precursors; ++p)
    	chi_computed+=delayed_spectra[mat][p][g]*delayed_fractions[mat][p];
    chi_computed+=prompt_spectra[mat][g]*(1.-delayed_fraction_sum[mat]);
    AssertRelease(chi_computed-chi[g][mat]<1e-10, "The relation between the spectral is not satisfied.");
    }
  }


  verbose_cout << "materials_vector: " << std::flush;
  if (verbose_cout.is_active())
    print_vector(materials_vector, false);
  verbose_cout << " Done!" << std::endl;
}

void Materials::parse_Valkin_file (std::string xs_file,
  const std::vector<unsigned int> &n_assemblies_per_dim,
  const unsigned int n_assemblies,
  const std::vector<unsigned int> geo_ps)
{
  n_mats = n_assemblies;
  Assert(fexists(xs_file), ExcMessage("xs_file doesn't exist"));
  std::ifstream input(xs_file.c_str(), std::ios::in);
  std::string keyword, str;
  double num;
  std::vector<unsigned int> n_cell_per_line;
  std::vector<double> sigma_a1;
  unsigned int line_num;

  sigma_a1.reserve(n_mats);
  D.resize(n_groups);
  sigma_t.resize(n_groups, std::vector<double>(n_mats));
  sigma_r.resize(n_groups);
  nu_sigma_f.resize(n_groups);
  chi.resize(n_groups, std::vector<double>(n_mats));
  sigma_f.resize(n_groups, std::vector<double>(n_mats));
  sigma_s.resize(n_groups);
  sigma_s[0].resize(n_groups);
  sigma_s[1].resize(n_groups);

  n_cell_per_line.resize(n_assemblies_per_dim[1]);
  unsigned int n_lines = n_assemblies_per_dim[1] * n_assemblies_per_dim[2];
  for (unsigned int i = 0; i < n_assemblies_per_dim[1]; i++)
  {
    n_cell_per_line[i] = geo_ps[2 * i + 1] - geo_ps[2 * i] + 1;
  }

// For every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);

    keyword.clear();
    iss >> keyword;

    if (is_commentary(keyword))
      continue;

    else if (keyword == "TIME=" || keyword == "TIME=0"
             || keyword == "TIME=0.0"
             || keyword == "time=0"
             || keyword == "Time=0")
    {
      continue;
    }
    // First definition Material and XSecs:
    else if (keyword == "D1")
    {
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;
        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;
          D[0].push_back(num);

        }
      }

    }

    else if (keyword == "D2")
    {
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;
        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;
          D[1].push_back(num);
        }
      }

    }

    else if (keyword == "SIGA1")
    {
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;
        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;
          sigma_a1.push_back(num);
        }
      }

    }

    else if (keyword == "SIGA2")
    {
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;
        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;
          sigma_r[1].push_back(num);
        }
      }

    }
    else if (keyword == "SIGR1")
    {
      // Remember! In this kind of input SIGR1= sigma_s12
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;
        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;
          sigma_s[0][1].push_back(num);
        }
      }
    }

    else if (keyword == "SGNUF1")
    {
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;
        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;
          nu_sigma_f[0].push_back(num);
        }
      }
    }

    else if (keyword == "SGNUF2")
    {
      for (unsigned int j = 0; j < n_lines; j++)
      {
        get_new_valid_line(input, line);
        std::istringstream iss2(line);
        iss2 >> str;

        if (lower_case(str) == "plane" || lower_case(str) == "plano"
            || lower_case(str) == "PLANE")
        {
          j--;
          continue;
        }
        line_num = str_to_num<unsigned int>(str);
        for (unsigned int k = 0; k < n_cell_per_line[line_num - 1];
            k++)
        {
          iss2 >> num;

          nu_sigma_f[1].push_back(num);
        }
      }
    }

    else
      // Error
      Assert(false, ExcMessage("Invalid Header in XSEC_file: " + keyword));
  }

// Define sigma_r1
  sigma_r[0].resize(n_mats);
  sigma_s[1][0].resize(n_mats);
  sigma_s[0][0].resize(n_mats);
  sigma_s[1][1].resize(n_mats);
  for (unsigned int i = 0; i < n_mats; i++)
  {
    sigma_r[0][i] = sigma_s[0][1][i] + sigma_a1[i];
    // chi
    chi[0][i] = 1.0;
    sigma_t[0][i] = 1.0 / (3.0 * D[0][i]);
    sigma_t[1][i] = 1.0 / (3.0 * D[1][i]);
  }

  sigma_f[0] = nu_sigma_f[0];
  sigma_f[1] = nu_sigma_f[1];

// Set materials vector sequentially
  materials_vector.resize(n_mats);
  for (unsigned int i = 0; i < n_assemblies; i++)
    materials_vector[i] = i;

}

/**
 * compute_diffusions_coefficient
 */
void Materials::compute_diffusions_coefficients ()
{
  diffusions_coefficient_computed = true;
// Compute diffusion Coefficients
  D.resize(n_groups, std::vector<double>(n_mats));
  for (unsigned int g = 0; g < n_groups; ++g)
    for (unsigned int mat = 0; mat < n_mats; ++mat)
      D[g][mat] = 1.0 / (3.0 * sigma_t[g][mat]);
}

/**
 * @brief Get the materials table as materials_table valid to
 */
template <int dim>
  void Materials::get_materials_table (
    Table<dim, types::material_id> &materials_table,
    const std::vector<unsigned int> &n_assemblies_per_dim)
  {
    materials_table = build_materials_table<dim>(
      geometry_matrix,
      n_assemblies_per_dim,
      materials_vector_with_holes);
  }

template void Materials::get_materials_table<1> (
  Table<1, types::material_id> &materials_table,
  const std::vector<unsigned int> &n_assemblies_per_dim);

template void Materials::get_materials_table<2> (
  Table<2, types::material_id> &materials_table,
  const std::vector<unsigned int> &n_assemblies_per_dim);

template void Materials::get_materials_table<3> (
  Table<3, types::material_id> &materials_table,
  const std::vector<unsigned int> &n_assemblies_per_dim);

/*
 * This should NOT be used, only defined to specialization below:
 */
template <int dim>
  Table<dim, types::material_id> Materials::build_materials_table (
    const std::vector<std::vector<unsigned int> > &geometry_matrix,
    const std::vector<unsigned int> &n_cells_per_dim,
    const std::vector<unsigned int> &materials)

  {
    AssertRelease(false, "setMaterialsTable ExcImpossibleInDim 1");
    Table<dim, types::material_id> materials_table;
    return materials_table;
  }

/*
 * Template Specialization for dim=2
 *  Fill materials_table a matrix that indicate which cell exist (1) or not (-1)
 */
template <>
  Table<1, types::material_id> Materials::build_materials_table<1> (
    const std::vector<std::vector<unsigned int> > &geometry_matrix,
    const std::vector<unsigned int> &n_cells_per_dim,
    const std::vector<unsigned int> &materials)
  {
    Table<1, types::material_id> materials_table(n_cells_per_dim[0]);
    unsigned int mat = 0;
    for (unsigned int i = 0; i < n_cells_per_dim[0]; ++i)
    {
      if (geometry_matrix[0][i] != 0)
      {
        materials_table[i] = static_cast<types::material_id>(materials[mat]);
        mat++;
      }
      else
        materials_table[i] = static_cast<types::material_id>(-1); // Hole
    }

    return materials_table;
  }

/*
 * Template Specialization for dim=2
 *  Fill materials_table a matrix that indicate which cell exist (1) or not (-1)
 */
template <>
  Table<2, types::material_id> Materials::build_materials_table<2> (
    const std::vector<std::vector<unsigned int> > &geometry_matrix,
    const std::vector<unsigned int> &n_cells_per_dim,
    const std::vector<unsigned int> &materials)
  {
    Table<2, types::material_id> materials_table(n_cells_per_dim[0],
      n_cells_per_dim[1]);

    unsigned int mat = 0;
    for (unsigned int j = 0; j < n_cells_per_dim[1]; j++)
      for (unsigned int i = 0; i < n_cells_per_dim[0]; i++)
      {
        if (geometry_matrix[j][i] != 0)
        {
          materials_table[i][j] = static_cast<types::material_id>(materials[mat]);
          mat++;
        }
        else
          materials_table[i][j] = static_cast<types::material_id>(-1); // Hole
      }
    return materials_table;
  }

/*
 *  Template Specialization for dim=3
 *  Fill materials_table a matrix that indicate which cell exist (1) or not (-1)
 */
template <>
  Table<3, types::material_id> Materials::build_materials_table<3> (
    const std::vector<std::vector<unsigned int> > &geometry_matrix,
    const std::vector<unsigned int> &n_cells_per_dim,
    const std::vector<unsigned int> &materials)
  {
    Table<3, types::material_id> materials_table(n_cells_per_dim[0],
      n_cells_per_dim[1], n_cells_per_dim[2]);
    unsigned int mat = 0;
    for (unsigned int k = 0; k < n_cells_per_dim[2]; k++)
      for (unsigned int j = 0; j < n_cells_per_dim[1]; j++)
        for (unsigned int i = 0; i < n_cells_per_dim[0]; i++)
          if (geometry_matrix[j][i] != 0)
          {
            materials_table[i][j][k] = static_cast<types::material_id>(materials[mat]);
            mat++;
          }
          else
            materials_table[i][j][k] = static_cast<types::material_id>(-1); // Hole

    return materials_table;
  }

/**
 * @brief Check if the cell in position (pos_x, pos_y, pos_z) is fuel.
 * It also advances an index if the cell is a hole or the reflector.
 */
bool Materials::is_fuel (const unsigned int pos_x,
  const unsigned int pos_y,
  const unsigned int,
  unsigned int &index) const
{
  if (listen_to_material_id)
    return true;
  if (geometry_matrix[pos_y][pos_x] == 0)
  {
    return false;
  }
// else if (geometry_matrix[pos_y][pos_x] == 2)
// {
//   index++;
//   return false;
// }
  else if (materials_vector_with_holes[index]
           == static_cast<unsigned int>(-1))
  {
    index++;
    return false;
  }

  return true;
}

/**
 * @brief Check if the cell in position (pos_x, pos_y, pos_z) is reflector.
 */
bool Materials::is_reflector (const unsigned int pos_x,
  const unsigned int pos_y,
  const unsigned int) const
{
  return (geometry_matrix[pos_y][pos_x] == 2);
}

/**
 * @brief Return the plane number where the cell_index pertains
 */
unsigned int Materials::plane (int cell_index) const
{
  unsigned int plane = 0;
  while (cell_index >= 0)
  {
    cell_index -= assemblies_per_plane[plane];
    plane++;
  }
  return plane - 1;
}

/**
 * @brief Set geometry_matrix from the geometry_points structure.
 */
void Materials::set_geometry_matrix (
  const std::vector<unsigned int> &assem_per_dim,
  const std::vector<unsigned int> &geo_ps)
{
  geometry_matrix.resize(assem_per_dim[1],
    std::vector<unsigned int>(assem_per_dim[0], 0));
  for (unsigned int i = 0; i < assem_per_dim[1]; ++i)
    for (unsigned int j = geo_ps[2 * i] - 1; j < geo_ps[2 * i + 1]; ++j)
      geometry_matrix[i][j] = 1;
}

/**
 * @brief Set geometry_matrix from a string (from the input file)
 */
void Materials::set_geometry_matrix (
  const std::vector<unsigned int> &assem_per_dim,
  const std::string &str)
{
  parse_matrix(str, geometry_matrix, assem_per_dim[1], assem_per_dim[0]);
}

/**
 * @brief Get a constant reference to geometry_matrix
 */
const std::vector<std::vector<unsigned int> >& Materials::get_geometry_matrix () const
{
  return geometry_matrix;
}

/**
 * @brief Check if the material in cell position (idx_x, idx_y, idx_z) exists.
 * @return true or false
 */
bool Materials::exist (
  const unsigned int idx_x,
  const unsigned int idx_y,
  const unsigned int) const
{
  return (geometry_matrix[idx_y][idx_x] != 0);
}
