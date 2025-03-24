/**
 * @file   input_geom.cc
 * @brief
 */

// This is for the xml parser from boost
#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "../../include/utils_base.h"
#include "../../include/io/input_geom.h"

void InputGeom::load (const std::string &filename)
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

  // core --------------------------------------------------------------
  // core.composed = pt.get<bool>("geometry.core.<xmlattr>.composed", false);
  core.name = pt.get<std::string>("geometry.core.name");
  core.dimension = pt.get<unsigned int>("geometry.core.dimension");

  bin = pt.get<std::string>("geometry.core.components");
  string_to_vector(bin, core.components);
  bin.clear();

  bin = pt.get<std::string>("geometry.core.boundary");
  string_to_vector(bin, core.boundary);
  bin.clear();

  if (core.dimension == 3)
  {
    core.n_planes = pt.get<double>("geometry.core.n_planes");
    bin = pt.get<std::string>("geometry.core.length");
    string_to_vector(bin, core.length);
    bin.clear();
  }
  else
  {
    core.n_planes = 1;
    core.length = std::vector<double>(1, 0.0);
  }

  // -------------------------------------------------------------------
  // Assert given values
  AssertRelease(core.length.size() == core.n_planes,
    "NOT a correct number of core.length");
  AssertRelease(core.components.size() == core.n_planes,
    "NOT a correct number of core.components");
  AssertRelease(core.boundary.size() == core.dimension,
    "NOT a correct number of Boundary Conditions");
  for (unsigned int d = 0; d < core.dimension; ++d)
    AssertRelease(core.boundary[d].size() == 2,
      "NOT a correct number of Boundary Conditions");

  // Planes Defined ------------------------------------------------------
  BOOST_FOREACH(ptree::value_type & v, pt.get_child("geometry.planes") )
  {
    if (v.first == "plane")
    {
      Plane plane;
      plane.id = v.second.get<unsigned int>("<xmlattr>.id");
      plane.name = v.second.get<std::string>("name");

      bin = v.second.get<std::string>("n_nodes");
      string_to_vector(bin, plane.n_nodes);

      bin.clear();

      bin = v.second.get<std::string>("components");
      string_to_vector(bin, plane.components);
      bin.clear();

      // -------------------------------------------------------------------
      // Assert given values
      AssertRelease(plane.n_nodes.size() == core.dimension,
        "NOT a correct number of plane.n_nodes");
      AssertRelease(plane.components.size() == plane.n_nodes[1],
        "NOT a correct number of plane.components");
      for (unsigned int i = 0; i < plane.n_nodes[1]; ++i)
        AssertRelease(plane.components.size() == plane.n_nodes[0],
          "NOT a correct number of plane.components");
      // -------------------------------------------------------------------

      planes.push_back(plane);
    }
  }

  // Pins  --------------------------------------------------------------
  BOOST_FOREACH(ptree::value_type & v,
      pt.get_child("geometry.pins") )
  {
    if (v.first == "pin")
    {
      Pin pin;
      pin.name = v.second.get<std::string>("name");
      pin.id = v.second.get<unsigned int>("<xmlattr>.id");

      std::string bin = v.second.get<std::string>("type");
      if (bin.compare(std::string("pin")) == 0)
      {
        pin.pin_type = Pin_type::pin;

        pin.pitch = v.second.get<double>("pitch");
        pin.fuel_radius = v.second.get<double>("fuel_radius");
        std::string str = v.second.get<std::string>("materials");
        str_to_vector<unsigned int>(str, pin.materials);
        bin.clear();

        // -------------------------------------------------------------------
        // Assert given values
        AssertRelease(pin.materials.size() == 2,
          " Only 2 material must be specified pin type geometry");
        AssertRelease(pin.pitch >= 2 * pin.fuel_radius,
          "pin.pitch must be bigger than 2* pin.fuel_radius in pin id "
          + pin.id);
      }
      else if (bin.compare(std::string("box")) == 0)
      {
        pin.pin_type = Pin_type::box;
        pin.pitch = v.second.get<double>("pitch");
        std::string str = v.second.get<std::string>("materials");
        str_to_vector<unsigned int>(str, pin.materials);
        AssertRelease(pin.materials.size() == 1,
          " Only 2 material must be specified box type geometry");
        bin.clear();
      }
      else
        AssertRelease(false, "Not a valid pin type, only valid pin or box");

      // -------------------------------------------------------------------

      pins.push_back(pin);
    }
  }

}
