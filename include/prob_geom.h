/**
 * @file   input_geom.cc
 * @brief
 */

#ifndef PROB_GEOM_H
#define PROB_GEOM_H

// all include files you need here
#include <deal.II/base/types.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_reordering.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/base/parameter_handler.h>

#include <cmath>
#include <math.h>

#include "input_geom.h"

using namespace dealii;

/**
 * Overload of the function merge_triangulation to deal with
 * a vector of triangulations (not just two).
 */
template <int dim, int spacedim>
  void merge_triangulations (
    const std::vector<Triangulation<dim, spacedim> > &tria,
    Triangulation<dim, spacedim> &result);

/**
 * Here the pin function. We assume the pin cell is a square of
 * length @p pitch with a circular subdomain of radius
 * @p pin_radius, centered at @p center.
 */
template <int dim>
  void pin_cell (Triangulation<dim> & /* tria */,
    const Point<2> & /* center */,
    const double /* pitch */,
    const double /* effective_pin_radius */,
    const unsigned int /*inner_material*/,
    const unsigned int /*outer_material*/);

/**
 *
 */
template <int dim>
  void print_pin (const Triangulation<dim> & tria);

/**
 *
 */
template <int dim>
  void
  assign_pin_manifolds (Triangulation<dim> & tria,
    const std::vector<Point<2>> pin_center,
    const std::vector<double> & pin_radius,
    const std::vector<Pin_type> & pin_types,
    const unsigned int first_manifold = 10);

/**
 *
 */
template <int dim>
  void
  refine_assembly (Triangulation<dim> & tria,
    const std::string refinement_model,
    const unsigned int n_ref = 0);

/**
 * @brief Calculate the effective radius to pass to the GridGenerator::pin_cell()
 * to maintain the fuel area.
 */
double calculate_effective_radius (
  const double radius,
  const unsigned int n_refinements);

/**
 *
 */
template <int dim>
  void create_assembly (Triangulation<dim> & tria_out,
    std::vector<Triangulation<2> > & tria_in,
    const std::vector<double> & z_lengths,
    const unsigned int n_assemblies_plane);

/**
 *
 */
template <int dim>
  void set_manifold_and_refine (Triangulation<dim> & tria,
    const std::vector<Point<2> > & pin_center,
    const std::vector<Pin_type> & pin_types,
    const unsigned int n_ref,
    const unsigned int first_manifold,
    const std::string refinement_model);

/**
 *
 */
template <int dim>
  void colorize_rectangular_triangulation (Triangulation<dim> & tria,
    const Point<dim> & p1,
    const Point<dim> & p2,
    double epsilon);

/**
 *
 */
//template <int dim>
//  void extrude (
//    const std::vector<Triangulation<2> > & /* tria_plane */,
//    Triangulation<dim> & /* tria_results */,
//    const std::vector<std::vector<double> > & /* z_lengths */);
/**
 *
 */
template <int dim>
  void pin_area (const Triangulation<dim> & tria);

/**
 *
 */
template <int dim>
  void make_composed_geometry (InputGeom & geom,
    Triangulation<dim> &tria,
    const unsigned int n_ref_radial,
    const unsigned int n_ref_axial,
    const std::string refinement_model = "local");

#endif
