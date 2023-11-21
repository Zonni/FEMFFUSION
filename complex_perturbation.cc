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

#include "complex_perturbation.h"
#include "input_complex_pert.h"

#include <stdlib.h>
#include <stdio.h>

#include "femffusion.h"
#include "matlab_io.h"

using namespace dealii;

/**
 *
 */
ComplexPerturbation::ComplexPerturbation (
  // @suppress("Class members should be properly initialized")
  const Materials &_materials,
  const unsigned int _dim,
  const bool verbose)
:
    materials(_materials),
    verbose_cout(std::cout, verbose),
    dim(_dim)
{
  freq = 0;
}

/**
 *
 */
void ComplexPerturbation::reinit (const std::string &dxs_file,
  const std::string &xsec_type,
  const std::string &_pert_type)
{
  pert_type = _pert_type;
  if (xsec_type == "XS2G")
  {
    if (pert_type == "Cell_Wise")
      parse_dxs_file(dxs_file);
    else if (pert_type == "Borders")
      parse_borders_file(dxs_file);
    else if (pert_type == "BordersHex")
      parse_borders_hex_file(dxs_file);
    else
      AssertRelease(false,
        "Not valid pert_type, only valid:  Cell_Wise | Borders | BordersHex");
  }
  else if (xsec_type == "XSEC") // FIXME CAMBIAR NOMBRE a XS_Multigroup
  {
    if (pert_type == "Cell_Wise")
      parse_dxs_XSEC(dxs_file);
    else
      AssertRelease(false, "Wrong pert_type " + pert_type);
  }
  else if (xsec_type == "XML")
  {
    if (pert_type == "Cell_Wise")
      parse_forest_dxs(dxs_file);
    else
      AssertRelease(false, "Wrong pert_type " + pert_type);
  }
  else
    AssertRelease(false, "Wrong xsec_type " + xsec_type);
}

/**
 * @brief Get the number of materials defined.
 * @return n_mats
 */
unsigned int ComplexPerturbation::get_n_mats () const
{
  return materials.get_n_mats();
}

/**
 *
 */
void ComplexPerturbation::set_frequency (const double freq_)
{
  freq = freq_;
}

/**
 *
 */
double ComplexPerturbation::get_frequency () const
{
  return freq;
}

/**
 * @brief
 * @return face_id
 */
unsigned int ComplexPerturbation::get_pertubation_face_id (
  const unsigned int mat_id,
  const unsigned int face,
  const unsigned int child,
  const unsigned int n_children) const
{
  if (n_children == 1)
  {
    return faces_id[mat_id][face];
  }
  else if (child_at_face(face, child, n_children, dim))
  {
    return faces_id[mat_id][face];
  }
  return -1;
}

/**
 * @brief
 * @return face_id
 */
unsigned int ComplexPerturbation::get_pertubation_face_id_hex (
  const unsigned int mat_id,
  const unsigned int quad_in_hex,
  const unsigned int face,
  const unsigned int child,
  const unsigned int n_children) const
{
  if (face == 0 or face == 2)
    return -1;

  if (n_children == 1)
  {
    return faces_id[mat_id][faces_map[quad_in_hex][face / 2]];
  }
  else if (child_at_face(face, child, n_children, dim))
  {
    return faces_id[mat_id][faces_map[quad_in_hex][face / 2]];
  }
  return -1;
}

/**
 * @brief
 * @return delta_sigma_a
 */
complex ComplexPerturbation::get_delta_sigma_t (const unsigned int group,
  const unsigned int mat) const
{
  AssertIndexRange(group, delta_sigma_t.size());
  AssertIndexRange(mat, delta_sigma_t[group].size());
  return delta_sigma_t[group][mat];
}

/**
 * @brief
 * @return delta_sigma_a
 */
complex ComplexPerturbation::get_delta_sigma_r (const unsigned int group,
  const unsigned int mat) const
{
  AssertIndexRange(group, delta_sigma_r.size());
  AssertIndexRange(mat, delta_sigma_r[group].size());
  return delta_sigma_r[group][mat];
}

/**
 * @brief Get \delta\Sigma_f.
 * Usually you will want to use materials.chi(to_group, mat) * get_delta_sigma_f (from_group, int mat)
 * @return delta_sigma_f
 */
complex ComplexPerturbation::get_delta_sigma_f (const unsigned int group,
  const unsigned int mat) const
{
  AssertIndexRange(group, delta_sigma_f.size());
  AssertIndexRange(mat, delta_sigma_f[group].size());
  return delta_sigma_f[group][mat]; //
}

/**
 * @brief Get \chi\delta\Sigma_f.
 * materials.chi(group_j, mat) * get_delta_sigma_f (group_i, mat)
 */
complex ComplexPerturbation::get_delta_sigma_f (
  const unsigned int group_i,
  const unsigned int group_j,
  const unsigned int pert_mat,
  const unsigned int mat
  ) const
{
  AssertIndexRange(group_i, delta_sigma_f.size());
  AssertIndexRange(pert_mat, delta_sigma_f[group_i].size());
  return materials.get_chi(group_j, mat) * delta_sigma_f[group_i][pert_mat]; //
}

/**
 *
 */
complex ComplexPerturbation::get_delta_sigma_s (const unsigned int group_i,
  const unsigned int group_j,
  const unsigned int mat) const
{
  AssertIndexRange(mat, delta_sigma_s[group_i][group_j].size());
  return delta_sigma_s[group_i][group_j][mat];
}

/**
 *
 */
#ifdef MATIO
void ComplexPerturbation::parse_coresim_file (
    const std::string& dxs_file)
{
  Assert(fexists(dxs_file), ExcMessage("dxs_file doesn't exist"));

  const unsigned int n_groups = materials.get_n_groups();
  const unsigned int n_mats = materials.get_n_mats();

  AssertRelease(n_groups == 1 or n_groups == 2,
      "Invalid number of groups: n_groups = " + std::to_string(n_groups));
  AssertRelease(n_mats > 0,
      "Invalid number of materials: n_mats = " + std::to_string(n_mats));

  // Resize structures
  delta_sigma_r.resize(n_groups);
  delta_sigma_f.resize(n_groups);
  delta_sigma_a.resize(n_groups);
  delta_sigma_t.resize(n_groups);
  delta_sigma_s.resize(n_groups,std::vector<std::vector<complex>>(n_groups));

  mat_t* matfp;
  matfp = Mat_Open(dxs_file.c_str(), MAT_ACC_RDONLY);
  AssertRelease(NULL != matfp, "Error opening MAT file " + dxs_file);

  verbose_cout << "     n_mats " << n_mats << std::endl;
  verbose_cout << "     n_groups " << n_groups << std::endl;

  if (n_groups == 1)
  {
    load_matlab_dxs(matfp, "dABS1", n_mats, delta_sigma_a[0]);
    load_matlab_dxs(matfp, "dNUFIS1", n_mats, delta_sigma_f[0]);
    delta_sigma_r[0] = delta_sigma_a[0];
    delta_sigma_t[0] = std::vector<std::complex<double> >(n_mats, 0.0);

  }
  else if (n_groups == 2)
  {
    load_matlab_dxs(matfp, "dABS1", n_mats, delta_sigma_a[0]);
    load_matlab_dxs(matfp, "dABS2", n_mats, delta_sigma_a[1]);
    load_matlab_dxs(matfp, "dNUFIS1", n_mats, delta_sigma_f[0]);
    load_matlab_dxs(matfp, "dNUFIS2", n_mats, delta_sigma_f[1]);
    load_matlab_dxs(matfp, "dREM", n_mats, delta_sigma_s[0][1]);

    delta_sigma_r[0].resize(n_mats);
    delta_sigma_s[0][0].resize(n_mats);
    delta_sigma_s[1][1].resize(n_mats);
    delta_sigma_s[1][0].resize(n_mats);
    for (unsigned int nm=0; nm<n_mats; nm++)
    {
      delta_sigma_r[0][nm]=delta_sigma_a[0][nm]+delta_sigma_s[0][1][nm];
    }
    delta_sigma_r[1]=delta_sigma_a[1];

    delta_sigma_t[0] = std::vector<std::complex<double> >(n_mats, 0.0);
    delta_sigma_t[1] = std::vector<std::complex<double> >(n_mats, 0.0);

  }
  else

  Mat_Close(matfp);

}
#endif

/**
 *
 */
void ComplexPerturbation::parse_dxs_file (const std::string &dxs_file)
{
  AssertRelease(fexists(dxs_file), "dxs_file doesn't exist " + dxs_file);

  const unsigned int n_groups = materials.get_n_groups();
  const unsigned int n_mats = materials.get_n_mats();

  // Parse dsx_file
  std::ifstream input(dxs_file.c_str(), std::ios::in);
  std::string sub, keyword;
  unsigned int mat;
  std::string str;
  double num_real, num_imag;

  // for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    keyword.clear();
    iss >> keyword;

    if (is_commentary(keyword))
      continue;
    else if (keyword == "dXS")
    {
      verbose_cout << std::endl;
      verbose_cout << "  dXS..." << std::endl;

      // Resize
      delta_sigma_a.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_r.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_t.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_f.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_s.resize(n_groups,
        std::vector<std::vector<complex>>(n_groups));

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          delta_sigma_s[g1][g2].resize(n_mats);

      for (unsigned int j = 0; j < n_mats; j++)
      {
        str = get_new_valid_line(input, line);
        std::istringstream iss(line);
        iss >> str;
        mat = str_to_num<unsigned int>(str) - 1;

        AssertRelease(mat == j,
          "Error in the dXS in line " + num_to_str(j + 1));
        verbose_cout << "    mat " << mat + 1 << std::endl;

        // delta_sigma_t1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_t[0][mat] = complex(num_real, num_imag);

        // delta_sigma_a1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_a[0][mat] = complex(num_real, num_imag);

        // delta_sigma_f1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_f[0][mat] = complex(num_real, num_imag);

        // delta_sigma_12
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_s[0][1][mat] = complex(num_real, num_imag);

        // ---------------------------------------------------------------------------------------------------
        // Second line
        get_new_valid_line(input, line);
        std::istringstream iss2(line);

        // delta_sigma_t2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_t[1][mat] = complex(num_real, num_imag);

        // delta_sigma_a2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_a[1][mat] = complex(num_real, num_imag);

        // delta_sigma_f2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(j + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(j + 1));
        delta_sigma_f[1][mat] = complex(num_real, num_imag);

        //------------------------------------
        delta_sigma_r[0][mat] = delta_sigma_s[0][1][mat]
                                + delta_sigma_a[0][mat];
        delta_sigma_r[1][mat] = delta_sigma_a[1][mat];

        verbose_cout << "    " << delta_sigma_r[0][mat] << "  "
                     << delta_sigma_f[0][mat]
                     << "  "
                     << delta_sigma_s[0][1][mat]
                     << std::endl;
        verbose_cout << "    " << delta_sigma_r[1][mat] << "  "
                     << delta_sigma_f[1][mat]
                     << std::endl;
      }
      verbose_cout << "  Done!" << std::endl;
    }
  }
}

/**
 *
 */
void ComplexPerturbation::parse_dxs_XSEC (const std::string &dxs_file)
{
  AssertRelease(fexists(dxs_file), "dxs_file doesn't exist " + dxs_file);

  const unsigned int n_groups = materials.get_n_groups();
  const unsigned int n_mats = materials.get_n_mats();

  // Parse dsx_file
  std::ifstream input(dxs_file.c_str(), std::ios::in);
  std::string sub, keyword;
  unsigned int mat = 0;
  std::string str;
  double num_real, num_imag;

  // for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    keyword.clear();
    iss >> keyword;

    if (is_commentary(keyword))
      continue;
    else if (keyword == "dXS")
    {
      verbose_cout << std::endl;
      verbose_cout << "  dXS..." << std::endl;

      // Resize
      delta_sigma_a.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_r.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_f.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_t.resize(n_groups, std::vector<complex>(n_mats));
      delta_sigma_s.resize(n_groups,
        std::vector<std::vector<complex>>(n_groups));

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          delta_sigma_s[g1][g2].resize(n_mats);

      for (unsigned int j = 0; j < n_mats; j++)
      {
        for (unsigned int g = 0; g < n_groups; g++)
        {
          str = get_new_valid_line(input, line);

          std::istringstream iss(line);

          if (g == 0)
          {
            iss >> str;
            mat = str_to_num<unsigned int>(str) - 1;
            AssertRelease(mat == j,
              "Error in the dXS in line " + num_to_str(j + 1));
            verbose_cout << "    mat " << mat + 1 << std::endl;
          }

          // delta_sigma_t
          iss >> num_real;
          AssertRelease(!iss.fail(),
            "There are not enough (well) dXS specified in line "
            + num_to_str(j + 1));
          iss >> num_imag;
          AssertRelease(!iss.fail(),
            "There are not enough (well) dXS specified in line "
            + num_to_str(j + 1));
          delta_sigma_t[g][mat] = complex(num_real, num_imag);

          // delta_sigma_f
          iss >> num_real;
          AssertRelease(!iss.fail(),
            "There are not enough (well) XSEC specified in line "
            + num_to_str(j + 1));
          iss >> num_imag;
          AssertRelease(!iss.fail(),
            "There are not enough (well) dXS specified in line "
            + num_to_str(j + 1));
          delta_sigma_f[g][mat] = complex(num_real, num_imag);

          for (unsigned int g2 = 0; g2 < n_groups; g2++)
          {
            // delta_sigma_s
            iss >> num_real;
            AssertRelease(!iss.fail(),
              "There are not enough (well) XSEC specified in line "
              + num_to_str(j + 1));
            iss >> num_imag;
            AssertRelease(!iss.fail(),
              "There are not enough (well) dXS specified in line "
              + num_to_str(j + 1));
            delta_sigma_s[g][g2][mat] = complex(num_real, num_imag);

          }

          /////////////////////////////////////////////////
          // Calculated DXS from above
          delta_sigma_a[g][mat] = delta_sigma_t[g][mat];
          for (unsigned int g2 = 0; g2 < n_groups; g2++)
            delta_sigma_a[g][mat] -= delta_sigma_s[g][g2][mat];

          delta_sigma_r[g][mat] = delta_sigma_t[g][mat] - delta_sigma_s[g][g][mat];

        }
      }
      verbose_cout << "  Done!" << std::endl;
    }
  }
}

/**
 *  @brief It parses the XS.xml file with a XML format.
 */
void ComplexPerturbation::parse_forest_dxs (const std::string &xml_file)
{
  verbose_cout << "parse_forest_dxs...  " << xml_file << std::endl;
  AssertRelease(fexists(xml_file), "forest_file doesn't exist");

  XMLInput::InputPert input;
  input.load(xml_file);
  const unsigned int n_groups = materials.get_n_groups();
  const unsigned int n_mats = materials.get_n_mats();

  AssertRelease(input.get_n_groups() == n_groups,
    "n_groups in xml file does not match "
    + num_to_str(n_groups)
    + " vs " + num_to_str(input.get_n_groups()));

  AssertRelease(input.get_n_mat() == n_mats,
    "n_mats in xml files does not match "
    + num_to_str(n_groups)
    + " vs " + num_to_str(input.get_n_groups()));

  // Resize Containers
  delta_sigma_t.resize(n_groups, std::vector<complex>(n_mats));
  delta_sigma_s.resize(n_groups,
    std::vector<std::vector<complex> >(n_groups, std::vector<complex>(n_mats)));
//  delta_chi.resize(n_groups, std::vector<double>(n_mats));
  delta_sigma_r.resize(n_groups, std::vector<complex>(n_mats));
//  delta_nu_sigma_f.resize(n_groups, std::vector<double>(n_mats));
  delta_sigma_f.resize(n_groups, std::vector<complex>(n_mats));

  for (unsigned int mat = 0; mat < n_mats; ++mat)
  {
    verbose_cout << "  Material " << mat << std::endl;
    AssertRelease(input.xs[mat].id == mat, "Error in mat ids");
    for (unsigned int from_g = 0; from_g < n_groups; ++from_g)
    {
      verbose_cout << "    Group " << from_g + 1 << std::endl;
      AssertRelease(input.xs[mat].exist_sigma_t, "Sigma_t does not exist");
      AssertRelease(input.xs[mat].exist_sigma_s, "Sigma_s does not exist");
      AssertRelease(input.xs[mat].exist_chi, "Chi does not exist");
      AssertRelease(input.xs[mat].exist_nu_sigma_f,
        "Nu Sigma_f does not exist");

      // sigma_t
      delta_sigma_t[from_g][mat] = input.xs[mat].sigma_t[from_g];
      verbose_cout << "        sigma_t_" << from_g + 1 << " = "
                   << delta_sigma_t[from_g][mat]
                   << std::endl;

//      // chi
//      delta_chi[from_g][mat] = input.xs[mat].chi[from_g];
//      verbose_cout << "        chi_" << from_g + 1 << " = "
//                   << chi[from_g][mat]
//                   << std::endl;

//      // nu_sigma_f
//      delta_nu_sigma_f[from_g][mat] = input.xs[mat].nu_sigma_f[from_g];
//      verbose_cout << "        nu_sigma_fv_" << from_g + 1 << " = "
//                   << nu_sigma_f[from_g][mat]
//                   << std::endl;

      // sigma_f
      // If sigma_f exists get it if not use nusigf instead
      if (!input.xs[mat].exist_sigma_f)
      {
        input.xs[mat].sigma_f = input.xs[mat].nu_sigma_f;
        input.xs[mat].exist_sigma_f = true;
      }
      delta_sigma_f[from_g][mat] = input.xs[mat].sigma_f[from_g];
      verbose_cout << "        sigma_f_" << from_g + 1 << " = "
                   << delta_sigma_f[from_g][mat]
                   << std::endl;

      // sigma_s
      for (unsigned int to_g = 0; to_g < n_groups; ++to_g)
      {
        // Be careful because in input.xs sigma_s is in a different way
        // Also be careful because we define sigma_s negative!
        delta_sigma_s[from_g][to_g][mat] =
                                           input.xs[mat].sigma_s[to_g][from_g];

        verbose_cout << "        sigma_s_intergroup_" << from_g + 1 << "->"
                     << to_g + 1
                     << " = " << delta_sigma_s[from_g][to_g][mat]
                     << std::endl;
      }
    }
  }

  // Calculate sigma_r
  for (unsigned int mat = 0; mat < n_mats; ++mat)
    for (unsigned int g = 0; g < n_groups; ++g)
      delta_sigma_r[g][mat] = delta_sigma_t[g][mat] - delta_sigma_s[g][g][mat];

  verbose_cout << "... Done!" << std::endl;
}

/**
 *
 */
void ComplexPerturbation::parse_borders_file (const std::string &dxs_file)
{
  AssertRelease(fexists(dxs_file), "dxs_file doesn't exist " + dxs_file);

  const unsigned int n_groups = materials.get_n_groups();
  const unsigned int n_mats = materials.get_n_mats();
  unsigned int n_faces_pert;
  AssertRelease(n_groups == 2, "Only implemented for 2 groups");

  // Parse dsx_file
  std::ifstream input(dxs_file.c_str(), std::ios::in);
  std::string sub, keyword;
  unsigned int mat, face_id, face;
  std::string str;
  double num_real, num_imag;
  unsigned int faces_per_cell = 0;
  if (dim == 1)
    faces_per_cell = GeometryInfo<1>::faces_per_cell;
  else if (dim == 2)
    faces_per_cell = GeometryInfo<2>::faces_per_cell;
  else
    faces_per_cell = GeometryInfo<3>::faces_per_cell;
  faces_id.resize(n_mats, std::vector<unsigned int>(faces_per_cell));
  // for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    keyword.clear();
    iss >> keyword;

    if (is_commentary(keyword))
      continue;
    else if (keyword == "Perturbed_Faces")
    {
      for (unsigned int mt = 0; mt < n_mats; mt++)
      {
        str = get_new_valid_line(input, line);
        std::istringstream iss(line);
        iss >> str;
        mat = str_to_num<unsigned int>(str) - 1;

        AssertRelease(mat == mt,
          "Error in the Perturbed_Faces in line "
          + num_to_str(mat + 1));
        verbose_cout << "    Faces " << mat + 1 << std::endl;

        for (unsigned int f = 0; f < faces_per_cell; f++)
        {
          iss >> face_id;
          faces_id[mat][f] = face_id - 1;

          AssertRelease(!iss.fail(),
            "There are not enough (well) Perturbed_Faces specified in line "
            + num_to_str(mat + 1));
        }
      }
    }
    else if (keyword == "dXS_Faces")
    {
      verbose_cout << std::endl;
      verbose_cout << "  dXS_Faces..." << std::endl;
      iss >> n_faces_pert;

      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of perturbed Faces defined!"));
      Assert(n_mats>0,
        ExcMessage("It must be defined the number of perturbed Faces defined!"));

      verbose_cout << "  n_faces_pert " << n_faces_pert << std::endl;
      // Resize
      delta_sigma_a.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_r.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_f.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_t.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_s.resize(n_groups,
        std::vector<std::vector<complex>>(n_faces_pert));

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          delta_sigma_s[g1][g2].resize(n_faces_pert);

      for (unsigned int fc = 0; fc < n_faces_pert; fc++)
      {
        str = get_new_valid_line(input, line);
        std::istringstream iss(line);
        iss >> str;
        face = str_to_num<unsigned int>(str) - 1;

        AssertRelease(face == fc,
          "Error in the dXS in line " + num_to_str(face + 1));
        verbose_cout << "    mat " << face + 1 << std::endl;

        // delta_sigma_t1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_t[0][face] = complex(num_real, num_imag);

        // delta_sigma_a1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_a[0][face] = complex(num_real, num_imag);

        // delta_sigma_f1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_f[0][face] = complex(num_real, num_imag);

        // delta_sigma_12
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_s[0][1][face] = complex(num_real, num_imag);
        // ---------------------------------------------------------------------------------------------------
        // Second line
        get_new_valid_line(input, line);
        std::istringstream iss2(line);

        // delta_sigma_t2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_t[1][face] = complex(num_real, num_imag);

        // delta_sigma_a2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_a[1][face] = complex(num_real, num_imag);

        // delta_sigma_f2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(fc + 1));
        delta_sigma_f[1][face] = complex(num_real, num_imag);

        //------------------------------------
        delta_sigma_r[0][face] = delta_sigma_s[0][1][face] + delta_sigma_a[0][face];
        delta_sigma_r[1][face] = delta_sigma_a[1][face];

        verbose_cout << "    " << delta_sigma_t[0][face] << "  "
                     << delta_sigma_r[0][face]
                     << "  "
                     << delta_sigma_f[0][face]
                     << "  " << delta_sigma_s[0][1][face]
                     << std::endl;
        verbose_cout << "    " << delta_sigma_t[1][face] << "  "
                     << delta_sigma_r[1][face]
                     << "  "
                     << delta_sigma_f[1][face]
                     << std::endl;

      }
      verbose_cout << "  Done!" << std::endl;
    }
  }
}

/**
 * @brief
 */
void ComplexPerturbation::parse_borders_hex_file (const std::string &dxs_file)
{
  AssertRelease(fexists(dxs_file), "dxs_file doesn't exist " + dxs_file);

  const unsigned int n_groups = materials.get_n_groups();
  const unsigned int n_mats = materials.get_n_mats();
  unsigned int n_faces_pert;

  // Parse dsx_file
  std::ifstream input(dxs_file.c_str(), std::ios::in);
  std::string sub, keyword;
  unsigned int mat, face_id, face;
  std::string str;
  double num_real, num_imag;
  unsigned int faces_per_cell;

  AssertRelease(n_groups == 2, "Only Valid for n_groups");
  AssertRelease(dim != 1, "BordersHex not valid in dim==1");
  if (dim == 2)
    faces_per_cell = 6;
  else
    faces_per_cell = 8;

  faces_id.resize(n_mats, std::vector<unsigned int>(faces_per_cell));
  // for every line
  for (std::string line; getline(input, line);)
  {
    std::istringstream iss(line);
    keyword.clear();
    iss >> keyword;

    if (is_commentary(keyword))
      continue;
    else if (keyword == "Perturbed_Faces")
    {
      for (unsigned int mt = 0; mt < n_mats; mt++)
      {
        str = get_new_valid_line(input, line);
        std::istringstream iss(line);
        iss >> str;
        mat = str_to_num<unsigned int>(str) - 1;

        AssertRelease(mat == mt,
          "Error in the Perturbed_Faces in line "
          + num_to_str(mat + 1));
        verbose_cout << "    Faces " << mat + 1 << std::endl;

        for (unsigned int f = 0; f < faces_per_cell; f++)
        {
          iss >> face_id;
          faces_id[mat][f] = face_id - 1;

          AssertRelease(!iss.fail(),
            "There are not enough (well) Perturbed_Faces specified in line "
            + num_to_str(mat + 1));
        }
      }
    }
    else if (keyword == "dXS_Faces")
    {
      verbose_cout << std::endl;
      verbose_cout << "  dXS_Faces..." << std::endl;
      iss >> n_faces_pert;

      Assert(!iss.fail(),
        ExcMessage("It must be defined the number of perturbed Faces defined!"));
      Assert(n_mats>0,
        ExcMessage("It must be defined the number of perturbed Faces defined!"));

      verbose_cout << "  n_faces_pert " << n_faces_pert << std::endl;
      // Resize
      delta_sigma_a.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_r.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_f.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_t.resize(n_groups, std::vector<complex>(n_faces_pert));
      delta_sigma_s.resize(n_groups,
        std::vector<std::vector<complex>>(n_faces_pert));

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          delta_sigma_s[g1][g2].resize(n_faces_pert);

      for (unsigned int fc = 0; fc < n_faces_pert; fc++)
      {
        str = get_new_valid_line(input, line);
        std::istringstream iss(line);
        iss >> str;
        face = str_to_num<unsigned int>(str) - 1;

        AssertRelease(face == fc,
          "Error in the dXS in line " + num_to_str(face + 1));
        verbose_cout << "    mat " << face + 1 << std::endl;

        // delta_sigma_t1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_t[0][face] = complex(num_real, num_imag);

        // delta_sigma_a1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_a[0][face] = complex(num_real, num_imag);

        // delta_sigma_f1
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_f[0][face] = complex(num_real, num_imag);

        // delta_sigma_12
        iss >> num_real;
        AssertRelease(!iss.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss >> num_imag;
        AssertRelease(!iss.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_s[0][1][face] = complex(num_real, num_imag);

        // ---------------------------------------------------------------------------------------------------
        // Second line
        get_new_valid_line(input, line);
        std::istringstream iss2(line);

        // delta_sigma_t2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_t[1][face] = complex(num_real, num_imag);

        // delta_sigma_a2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(face + 1));
        delta_sigma_a[1][face] = complex(num_real, num_imag);

        // delta_sigma_f2
        iss2 >> num_real;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) XSEC specified in line "
          + num_to_str(face + 1));
        iss2 >> num_imag;
        AssertRelease(!iss2.fail(),
          "There are not enough (well) dXS specified in line "
          + num_to_str(fc + 1));
        delta_sigma_f[1][face] = complex(num_real, num_imag);

        //------------------------------------
        delta_sigma_r[0][face] = delta_sigma_s[0][1][face] + delta_sigma_a[0][face];
        delta_sigma_r[1][face] = delta_sigma_a[1][face];

        verbose_cout << "    " << delta_sigma_t[0][face] << "  "
                     << delta_sigma_r[0][face]
                     << "  "
                     << delta_sigma_f[0][face]
                     << "  " << delta_sigma_s[0][1][face]
                     << std::endl;
        verbose_cout << "    " << delta_sigma_t[1][face] << "  "
                     << delta_sigma_r[1][face]
                     << "  "
                     << delta_sigma_f[1][face]
                     << std::endl;
      }
      verbose_cout << "  Done!" << std::endl;
    }
  }
}
