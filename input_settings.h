
// Include guard
#ifndef INPUT_SETTINGS_H
#define INPUT_SETTINGS_H

// This is for the xml parser from boost
#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "utils_base.h"

// different eigenvalue solvers:
// PI: power iteration
// AM: Arnoldi method
enum class Eig_solver {PI, AM};

class InputSettings
{
public:

  struct Files
  {
    std::string geom;
    std::string mat;
  };

  struct Problem
  {
    std::string type;
    std::string approximation;
    unsigned int sn;
  };

  struct Algebra
  {
    Eig_solver type;
    double tol;
    unsigned int max_it;
  };

  struct FE_settings
  {
    unsigned int degree;
  };
  
  unsigned int dim;

  //Files files;
  Problem problem;
  Algebra algebra;
  FE_settings fe_settings;
  
  void load(const std::string &filename);
  void save(const std::string &filename);
  
};

void InputSettings::load(const std::string &filename)
{
  // Create empty property tree object
  using boost::property_tree::ptree;
  ptree pt;

  // Load XML file and put its contents in property tree. 
  // No namespace qualification is needed, because of Koenig 
  // lookup on the second argument. If reading fails, exception
  // is thrown.
  read_xml(filename, pt, boost::property_tree::xml_parser::trim_whitespace);

  std::string bin; // bin is a container to save temporal data

  dim = pt.get<unsigned int>("settings.dim");
  
  // files -------------------------------------------------------------
  /*
  files.geom = pt.get<std::string>("settings.input_files.geom", filename+"geom.xml");
  files.mat = pt.get<std::string>("settings.input_files.mat", filename+"mat.xml");
  */

  // problem -----------------------------------------------------------
  problem.type = pt.get<std::string>("settings.problem.<xmlattr>.type");
  problem.approximation = pt.get<std::string>("settings.problem.approximation.<xmlattr>.type");
  problem.sn = pt.get<unsigned int>("settings.problem.approximation.sn");

  // algebra -----------------------------------------------------------
  bin = pt.get<std::string>("settings.algebra.eig_solver.<xmlattr>.type");
  if (bin.compare(std::string("PI"))==0)
    algebra.type = Eig_solver::PI;
  else
    algebra.type = Eig_solver::AM;
        
  algebra.tol = pt.get<double>("settings.algebra.eig_solver.tol");
  algebra.max_it = pt.get<unsigned int>("settings.algebra.eig_solver.max_it");


  // fe_settings -------------------------------------------------------
  fe_settings.degree = pt.get<unsigned int>("settings.fe_settings.degree");

}

void InputSettings::save(const std::string &filename)
{
  // Create empty property tree object
  using boost::property_tree::ptree;
  ptree pt;

  std::string bin; // bin is a container to save temporal data
  std::string indent("\t");
  
  pt.add("settings", "");
  
  pt.add("settings.<xmlcomment>"," \n \
    Here we define the different materials associated to the \n \
    material id specified when defining the geometry");
  
  pt.add("settings.<xmlcomment>", "dimension");
  pt.put("settings.dim", dim);
  
  // files -------------------------------------------------------------
  /*
  pt.add("settings.<xmlcomment>", "input files");
  pt.add("settings.input_files", "");
  pt.put("settings.input_files.geom", files.geom);
  pt.put("settings.input_files.mat", files.mat);
  */

  // problem -----------------------------------------------------------
  pt.add("settings.<xmlcomment>", "problem");
  pt.add("settings.problem", "");
  pt.put("settings.problem.<xmlattr>.type",problem.type);
  pt.put("settings.problem.approximation.<xmlattr>.type", problem.approximation);
  pt.put("settings.problem.approximation.sn", problem.sn);

  // algebra -----------------------------------------------------------
  pt.add("settings.<xmlcomment>", "algebra");
  pt.add("settings.algebra", "");
  //pt.put("settings.algebra.eig_solver.<xmlattr>.type", algebra.type);
  pt.put("settings.algebra.eig_solver.tol", algebra.tol);
  pt.put("settings.algebra.eig_solver.max_it", algebra.max_it);

  // fe_settings -------------------------------------------------------
  pt.add("settings.<xmlcomment>", "fe_settings");
  pt.add("settings.fe_settings", "");
  pt.put("settings.fe_settings.degree", fe_settings.degree);
  
  // Write property tree to XML file (\t is a tabulator)
#if BOOST_VERSION < 105500
    boost::property_tree::xml_writer_settings<char> settings ('\t', 1);
#else
    boost::property_tree::xml_writer_settings<std::string> settings ('\t', 1);
#endif
  write_xml(filename, pt, std::locale(), settings);
}

#endif
