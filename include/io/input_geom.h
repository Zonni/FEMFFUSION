/**
 * @file   input_geom.h
 * @brief
 */

#ifndef INPUT_GEOM_H
#define INPUT_GEOM_H

// This is for the xml parser from boost
#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "../../include/femffusion.h"

enum class Pin_type
{
  pin,
  box
};

class InputGeom
{
  public:

  struct Core
  {
    unsigned int dimension;
    std::string name;
    double n_planes;
    std::vector<unsigned int> components;
    std::vector<std::vector<unsigned int> > boundary;
    std::vector<double> length;
  };

  struct Plane
  {
    std::string name;
    unsigned int id;
    std::vector<unsigned int> n_nodes;
    std::vector<std::vector<unsigned int> > components;
  };

  struct Pin
  {
    std::string name;
    unsigned int id;
    Pin_type pin_type; // pin or  box
    double pitch;
    double fuel_radius;
    std::vector<unsigned int> materials;

  };

  Core core;
  std::vector<Plane> planes;
  std::vector<Pin> pins;

  void load (const std::string &filename);
};

#endif

