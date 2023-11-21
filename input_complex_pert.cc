/**
 * @file   input_mat.cc
 * @brief  Implementation of class InputPert
 */

#include "input_complex_pert.h"

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>                      // for basic_filebuf<>::int_type, etc
#include <sstream>                      // for basic_stringbuf<>::int_type, etc

#include <cassert>
#include <cmath>
#include <complex>

#include "utils.h"

namespace XMLInput
{

/**
 *
 */
void InputPert::load (const std::string &filename)
{
  // Create empty property tree object
  using boost::property_tree::ptree;
  ptree pt;
  std::string bin; // bin is a container to store temporal data

  // Load XML file and put its contents in property tree.
  // No namespace qualification is needed, because of Koenig
  // lookup on the second argument. If reading fails, exception
  // is thrown.
  read_xml(filename, pt, boost::property_tree::xml_parser::trim_whitespace);

  n_groups = pt.get<unsigned int>("materials.<xmlattr>.ngroups");

  // Get Composition
  boost::optional<std::string> val = pt.get_optional<std::string>(
    "Composition");
  if (val.is_initialized())
  {
    bin = val.get();
    str_to_vector<unsigned int>(bin, materials_vector);
    bin.clear();
  }
  // cross sections ----------------------------------------------------
  BOOST_FOREACH(ptree::value_type & v, pt.get_child("materials") )
  { if (v.first == "mix")
  {
    DXS_single mat;
    mat.id = v.second.get<unsigned int>("<xmlattr>.id");
    mat.name = v.second.get<std::string>("name");

    bin = v.second.get<std::string>("SigmaT", std::string(""));
    if (bin != std::string(""))
    {
      std::vector<double> v;
      str_to_vector_complex(bin, mat.sigma_t);
      bin.clear();
      mat.exist_sigma_t = true;
    }

    bin = v.second.get<std::string>("SigmaS", std::string(""));
    if (bin != std::string(""))
    {
      str_to_vector_complex(bin, mat.sigma_s);
      bin.clear();
      mat.exist_sigma_s = true;
    }

    bin = v.second.get<std::string>("SigmaA", std::string(""));
    if (bin != std::string(""))
    {
      str_to_vector_complex(bin, mat.sigma_a);
      bin.clear();
      mat.exist_sigma_a = true;
    }

    bin = v.second.get<std::string>("NuSigF", std::string(""));
    if (bin != std::string(""))
    {
      str_to_vector_complex(bin, mat.nu_sigma_f);
      bin.clear();
      mat.exist_nu_sigma_f = true;
    }

    bin = v.second.get<std::string>("Nu", std::string(""));
    if (bin != std::string(""))
    {
      str_to_vector_complex(bin, mat.nu);
      bin.clear();
      mat.exist_nu = true;
    }

    bin = v.second.get<std::string>("SigF", std::string(""));
    if (bin != std::string(""))
    {
      str_to_vector_complex(bin, mat.sigma_f);
      bin.clear();
      mat.exist_sigma_f = true;
    }

    bin = v.second.get<std::string>("Chi", std::string(""));
    if (bin != std::string(""))
    {
      str_to_vector_complex(bin, mat.chi);
      bin.clear();
      mat.exist_chi = true;
    }

    xs.insert(DXS_pair(mat.id, mat));
  }
}
  n_mat = xs.size();
}

/**
 *
 */
void InputPert::check_nusigf (DXS_single &xs_) const
  {
  const double nusigf_threshold_rel = 1.e-3;
  const double nusigf_threshold_abs = 1.e-5;
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    complex calc_nusigf = xs_.nu[g] * xs_.sigma_f[g];
    if (std::abs(calc_nusigf) < nusigf_threshold_rel)
    {
      if (std::abs(calc_nusigf - xs_.nu_sigma_f[g]) > nusigf_threshold_abs)
      {
        std::cout
        << "Material id = "
        << xs_.id
        << " has a different nu_sigma_f cross-section "
        << "than the product of its nu and fission cross-sections for group "
        << g
        << " : nu_sigma_f = "
        << std::setprecision(6)
        << xs_.nu_sigma_f[g]
        << " calc_nusigf = " << std::setprecision(6) << calc_nusigf
        << "\n";
        assert(false);
      }
    }
    else
    {
      if (std::abs((calc_nusigf - xs_.nu_sigma_f[g]) / calc_nusigf)
          > nusigf_threshold_rel)
      {
        std::cout
        << "Material id = "
        << xs_.id
        << " has a different nu_sigma_f cross-section "
        << "than the product of its nu and fission cross-sections for group "
        << g
        << " : nu_sigma_f = "
        << std::setprecision(6)
        << xs_.nu_sigma_f[g]
        << " calc_nusigf = " << std::setprecision(6) << calc_nusigf
        << "\n";
        assert(false);
      }
    }
  }
}

/**
 *
 */
void InputPert::calc_nusigf (DXS_single &xs_)
{
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    complex calc_nusigf = xs_.nu[g] * xs_.sigma_f[g];
    xs_.nu_sigma_f[g] = calc_nusigf;
  }
}

/**
 *
 */
void InputPert::check_sigmat (DXS_single &xs_) const
  {
  const double sigmat_threshold = 1.e-4;
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    complex calc_sigmat = xs_.sigma_a[g];
    for (unsigned int h = 0; h < n_groups; ++h)
    {
      calc_sigmat += xs_.sigma_s[h][g];
    }
    if (std::abs(calc_sigmat - xs_.sigma_t[g]) > sigmat_threshold)
    {
      std::cout
      << "Material id = "
      << xs_.id
      << " has a different total cross-section "
      << "than the sum of its scattering and absorption cross-sections for group "
      << g
      << " : sigma_t = "
      << std::setprecision(6) << xs_.sigma_t[g]
      << " calc_sigma_t = "
      << std::setprecision(6) << calc_sigmat
      << "\n";
      assert(false);
    }
  }
}

void InputPert::calc_sigmat (DXS_single &xs_)
{
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    complex calc_sigmat = xs_.sigma_a[g];
    for (unsigned int h = 0; h < n_groups; ++h)
    {
      calc_sigmat += xs_.sigma_s[h][g];
    }
    xs_.sigma_t[g] = calc_sigmat;
  }
}

void InputPert::check ()
{
  for (typename DXS_map::iterator mat = xs.begin(); mat != xs.end(); ++mat)
  {

    if (mat->second.exist_nu_sigma_f && mat->second.exist_nu
        && mat->second.exist_sigma_f)
    {
      check_nusigf(xs[mat->second.id]);
      calc_nusigf(xs[mat->second.id]);
    }

    if (mat->second.exist_sigma_t && mat->second.exist_sigma_s
        && mat->second.exist_sigma_a)
    {
      check_sigmat(xs[mat->second.id]);
      calc_sigmat(xs[mat->second.id]);
    }

#ifdef DEBUG
    bool problem_defined = false;
    if (mat->second.exist_sigma_t && mat->second.exist_sigma_s
        && mat->second.exist_chi && mat->second.exist_nu_sigma_f)
    {
      problem_defined = true;
    }
    assert(problem_defined);
#endif
  }
}

} // end of namespace XMLInput
