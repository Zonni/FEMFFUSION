/**
 * @file   input_mat.cc
 * @brief  Implementation of class InputMat
 */

#include "input_mat.h"
#include "femffusion.h"

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

namespace XMLInput
{

void InputMat::load (const std::string &filename)
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
  n_precursors = pt.get<unsigned int>("materials.<xmlattr>.nprecursors", 0);

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
  {
    if (v.first == "mix")
    {
      XS_single mat;
      mat.id = v.second.get<unsigned int>("<xmlattr>.id");
      mat.name = v.second.get<std::string>("name");

      bin = v.second.get<std::string>("SigmaT", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.sigma_t);
        bin.clear();
        mat.exist_sigma_t = true;
      }

      bin = v.second.get<std::string>("SigmaS", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector(bin, mat.sigma_s);
        bin.clear();
        mat.exist_sigma_s = true;
      }

      bin = v.second.get<std::string>("SigmaA", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.sigma_a);
        bin.clear();
        mat.exist_sigma_a = true;
      }

      /**************************************************************/
      bin = v.second.get<std::string>("SigmaTR", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.sigma_tr);
        bin.clear();

      }
      bin = v.second.get<std::string>("SigmaR", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.sigma_r);
        bin.clear();
      }
      /**************************************************************/

      bin = v.second.get<std::string>("NuSigF", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.nu_sigma_f);
        bin.clear();
        mat.exist_nu_sigma_f = true;
      }

      bin = v.second.get<std::string>("Nu", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.nu);
        bin.clear();
        mat.exist_nu = true;
      }

      bin = v.second.get<std::string>("SigF", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.sigma_f);
        bin.clear();
        mat.exist_sigma_f = true;
      }

      bin = v.second.get<std::string>("Chi", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.chi);
        bin.clear();
        mat.exist_chi = true;
      }

      bin = v.second.get<std::string>("Beta", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.beta);
        bin.clear();
        mat.exist_beta = true;
      }

      bin = v.second.get<std::string>("Lambda", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.lambda);
        bin.clear();
        mat.exist_lambda = true;
      }

      bin = v.second.get<std::string>("ChiD", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector(bin, mat.chi_d);
        bin.clear();
        mat.exist_chi_d = true;
      }

      bin = v.second.get<std::string>("Velocities", std::string(""));
      if (bin != std::string(""))
      {
        str_to_vector<double>(bin, mat.velocities);
        bin.clear();
        mat.exist_velocities = true;
      }

      AssertRelease(mat.sigma_tr.size() > 0 or mat.sigma_t.size() > 0,
        "Sigma Tr or Sigma T must be defined");
      if (mat.sigma_tr.size() == 0)
      {

        mat.sigma_tr = mat.sigma_t;
      }

      if (mat.sigma_r.size() == 0)
      {
        mat.sigma_r.resize(n_groups);
        // Calculate sigma_r
        for (unsigned int g = 0; g < n_groups; ++g)
          (mat.sigma_r)[g] = (mat.sigma_t)[g] - (mat.sigma_s)[g][g];
      }
      xs.insert(XS_pair(mat.id, mat));
    }
  }
  n_mat = xs.size();

  check();

}

void InputMat::save (const std::string &filename)
{
  // Create empty property tree object
  using boost::property_tree::ptree;
  ptree pt;

  std::string bin;  // bin is a container to save temporal data
  std::string indent("\t");

  // materials ---------------------------------------------------------
  pt.add("materials", "");
  pt.put("materials.<xmlattr>.ngroups", n_groups);
  pt.put("materials.<xmlattr>.nprecursors", n_precursors);

  for (typename XS_map::iterator mat = xs.begin(); mat != xs.end(); ++mat)
  {
    // here we generate the node
    ptree &node = pt.add("materials.mix", "");

    node.put("<xmlattr>.id", mat->first);
    node.put("name", mat->second.name);

    if (mat->second.exist_sigma_t)
    {
      vector_to_str(mat->second.sigma_t, bin);
      node.put("SigmaT", bin);
      bin.clear();
    }

    if (mat->second.exist_sigma_a)
    {
      vector_to_str(mat->second.sigma_a, bin);
      node.put("SigmaA", bin);
      bin.clear();
    }

    if (mat->second.exist_nu_sigma_f)
    {
      vector_to_str(mat->second.nu_sigma_f, bin);
      node.put("NuSigF", bin);
      bin.clear();
    }

    if (mat->second.exist_nu)
    {
      vector_to_str(mat->second.nu, bin);
      node.put("Nu", bin);
      bin.clear();
    }

    if (mat->second.exist_sigma_f)
    {
      vector_to_str(mat->second.sigma_f, bin);
      node.put("SigF", bin);
      bin.clear();
    }

    if (mat->second.exist_chi)
    {
      vector_to_str(mat->second.chi, bin);
      node.put("Chi", bin);
      bin.clear();
    }

    if (mat->second.exist_sigma_s)
    {
      vector_to_str(mat->second.sigma_s, bin, indent, 3);
      node.put("SigmaS", bin);
      bin.clear();
    }

  }

  // Write property tree to XML file (\t is a tabulator)
#if BOOST_VERSION < 105500
  boost::property_tree::xml_writer_settings<char> settings('\t', 1);
#else
  boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);
#endif
  write_xml(filename, pt, std::locale(), settings);
}

void InputMat::check_nusigf (XS_single &xs_) const
{
  const double nusigf_threshold_rel = 1.e-3;
  const double nusigf_threshold_abs = 1.e-5;
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    double calc_nusigf = xs_.nu[g] * xs_.sigma_f[g];
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

void InputMat::calc_nusigf (XS_single &xs_)
{
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    double calc_nusigf = xs_.nu[g] * xs_.sigma_f[g];
    xs_.nu_sigma_f[g] = calc_nusigf;
  }
}

void InputMat::check_sigmat (XS_single &xs_) const
{
  const double sigmat_threshold = 1.e-2;
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    double calc_sigmat = xs_.sigma_a[g];
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

void InputMat::calc_sigmat (XS_single &xs_)
{
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    double calc_sigmat = xs_.sigma_a[g];
    for (unsigned int h = 0; h < n_groups; ++h)
    {
      calc_sigmat += xs_.sigma_s[h][g];
    }
    xs_.sigma_t[g] = calc_sigmat;
  }
}

void InputMat::calc_chip (XS_single &xs_)
{

  double betaeff;
  xs_.chi_p.resize(n_groups);

  for (unsigned int g = 0; g < n_groups; g++)
  {
    betaeff = 0;
    double sum_del_beta = 0.0;
    for (unsigned int p = 0; p < n_precursors; p++)
    {
      sum_del_beta += xs_.chi_d[p][g] * xs_.beta[p];
      betaeff += xs_.beta[p];
    }

    xs_.chi_p[g] = (xs_.chi[g] - sum_del_beta) / (1 - betaeff);


  }

}

void InputMat::calc_betaeff (XS_single &xs_)
{

  double betaeff = 0.0;

  for (unsigned int p = 0; p < n_precursors; p++)
  {
    betaeff += xs_.beta[p];
  }

  xs_.beta_eff = betaeff;

}

void InputMat::norm_chi (XS_single &xs_)
{
  double sum_chi = 0.;
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    sum_chi += xs_.chi[g];
  }

  if (sum_chi != 0.0)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      xs_.chi[g] = xs_.chi[g] / sum_chi;
    }
  }
}

void InputMat::check ()
{
  for (typename XS_map::iterator mat = xs.begin(); mat != xs.end(); ++mat)
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

    if (mat->second.exist_chi)
    {
      norm_chi(xs[mat->second.id]);
      calc_chip(xs[mat->second.id]);
      calc_betaeff(xs[mat->second.id]);
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
