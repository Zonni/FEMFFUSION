/**
 * @file   prob_geom.cc
 * @brief
 */

// all include files you need here
#include <deal.II/base/types.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.templates.h>
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

#include <cmath>
#include <math.h>

#include "../../include/io/input_geom.h"
#include "../../include/io/prob_geom.h"

using namespace dealii;
/**
 * Overload of the function merge_triangulation to deal with
 * a vector of triangulations (not just two).
 */
template <int dim, int spacedim>
  void
  merge_triangulations (
    const std::vector<Triangulation<dim, spacedim> > &tria,
    Triangulation<dim, spacedim> &result)
  {

    const unsigned int n_triangulations = tria.size();

    // Check that the number of triangulation is at least two
    Assert(n_triangulations >= 2,
      ExcMessage ("The number triangulations must be at least two."));

    // Check that the triangulations are coarse meshes
    for (unsigned int i = 0; i < n_triangulations; ++i)
    {
      Assert(tria[i].n_levels() == 1,
        ExcMessage ("The input triangulations must be coarse meshes."));
    }

    // get the union of the set of vertices
    std::vector<Point<spacedim> > vertices = tria[0].get_vertices();
    for (unsigned int i = 1; i < n_triangulations; ++i)
    {
      vertices.insert(vertices.end(),
        tria[i].get_vertices().begin(),
        tria[i].get_vertices().end());
    }

    // now form the union of the set of cells. note that we have to
    // translate the vertex indices
    std::vector<CellData<dim> > cells;
    unsigned int n_cells = 0;
    for (unsigned int i = 0; i < n_triangulations; ++i)
    {
      n_cells += tria[i].n_cells();
    }
    cells.reserve(n_cells);
    std::vector<unsigned int> user_indices;
    user_indices.reserve(n_cells);
    unsigned int n_vertices = 0;
    for (unsigned int i = 0; i < n_triangulations; ++i)
    {
      for (typename Triangulation<dim, spacedim>::cell_iterator
      cell = tria[i].begin(); cell != tria[i].end(); ++cell)
      {
        CellData<dim> this_cell;
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          this_cell.vertices[v] = cell->vertex_index(v) + n_vertices;
        this_cell.material_id = cell->material_id();
        user_indices.push_back(i);
        cells.push_back(this_cell);

      }
      n_vertices += tria[i].n_vertices();
    }

    // we do the loop because we have more than two triangulations, and
    // delete_duplicated_vertices only works for 2 triangulations.
    for (unsigned int i = 0; i < n_triangulations - 1; ++i)
    {
      // throw out duplicated vertices from the two meshes, reorder vertices as
      // necessary and create the triangulationeffective_pin_radius
      SubCellData subcell_data;
      std::vector<unsigned int> considered_vertices;
      GridTools::delete_duplicated_vertices(vertices, cells,
        subcell_data,
        considered_vertices);
    }
    // One more time
    SubCellData subcell_data;
    std::vector<unsigned int> considered_vertices;
    GridTools::delete_duplicated_vertices(vertices, cells,
      subcell_data,
      considered_vertices);

    // reorder the cells to ensure that they satisfy the convention for
    // edge and face directions
    GridReordering<dim, spacedim>::reorder_cells(cells, true);
    result.clear();
    result.create_triangulation(vertices, cells, subcell_data);

    // Set the pertinent Rod (user_index) for each cell.
    typename Triangulation<dim>::active_cell_iterator cell = result.begin_active();
    for (unsigned int i = 0; cell != result.end(); ++cell, ++i)
    {
      cell->set_user_index(user_indices[i]);
    }

  }

/**
 * Here the pin function. We assume the pin cell is a square of
 * length @p pitch with a circular subdomain of radius
 * @p pin_radius, centered at @p center.
 */
template <int dim>
  void
  pin_cell (Triangulation<dim> & /* tria */,
    const Point<2> & /* center */,
    const double /* pitch */,
    const double /* effective_pin_radius */,
    const unsigned int /*inner_material*/,
    const unsigned int /*outer_material*/)
  {
    Assert(false, ExcMessage ("Not implemented."));
  }

/**
 * TODO Cache this function
 */
template <>
  void
  pin_cell (Triangulation<2> & tria,
    const Point<2> & center,
    const double pitch,
    const double effective_pin_radius,
    const unsigned int inner_material,
    const unsigned int outer_material)
  {
    const unsigned int n_verts = 12;
    const unsigned int n_cells = 9;

    const double p = pitch / 2.0;
    // effective_pin_radius is the diagonal and r is the cathetus for x y axis
    const double r = effective_pin_radius / M_SQRT2;
    const double r2 = r / M_SQRT2;

    const double xx = center(0);
    const double yy = center(1);

    /**
     *  Defining special vertices manually
     *
     *  2-------------3
     *  | \         / |
     *  |  6-------7  |
     *  |  | \   / |  |
     *  |  | 10-11 |  |
     *  |  |  8-9  |  |
     *  |  | /   \ |  |
     *  |  4-------5  |
     *  | /         \ |
     *  0-------------1
     */
    const Point<2> vertices[n_verts] =
                                         {
                                           // box
                                           Point<2>(xx - p, yy - p),
                                           Point<2>(xx + p, yy - p),
                                           Point<2>(xx - p, yy + p),
                                           Point<2>(xx + p, yy + p),
                                           // circle
                                           Point<2>(xx - r, yy - r),
                                           Point<2>(xx + r, yy - r),
                                           Point<2>(xx - r, yy + r),
                                           Point<2>(xx + r, yy + r),
                                           // circle2
                                           Point<2>(xx - r2, yy - r2),
                                           Point<2>(xx + r2, yy - r2),
                                           Point<2>(xx - r2, yy + r2),
                                           Point<2>(xx + r2, yy + r2) };

    /**
     *  Define connections following the deal.II standards
     *       3
     *    2-->--3
     *    |     |
     *   0^     ^1
     *    |     |
     *    0-->--1
     *       2
     *
     * This cell is {0, 1, 2, 3}
     */
    const int cell_vertices[n_cells][GeometryInfo<2>::vertices_per_cell] =
          {
              { 0, 1, 4, 5 },
              { 0, 4, 2, 6 },
              { 2, 6, 3, 7 },
              { 1, 3, 5, 7 },
              { 4, 5, 8, 9 },
              { 4, 8, 6, 10 },
              { 6, 10, 7, 11 },
              { 5, 7, 9, 11 },
              { 8, 9, 10, 11 } };

    std::vector<CellData<2> > cells(n_cells, CellData<2>());

    // creating arrays with vertex info
    for (unsigned int i = 0; i < n_cells; ++i)
    {
      for (unsigned int j = 0; j < GeometryInfo<2>::vertices_per_cell; ++j)
      {
        cells[i].vertices[j] = cell_vertices[i][j];
      }
      if (i > 3) // inner cells
        cells[i].material_id = inner_material;
      else
        // outer cells
        cells[i].material_id = outer_material;
    }

    // Create the triangulation with no boundary information
    tria.create_triangulation(
      std::vector<Point<2> >(&vertices[0], &vertices[n_verts]),
      cells,
      SubCellData());

  }

/**
 * Here the pin function. We assume the pin cell is a square of
 * length @p pitch with a circular subdomain of radius
 * @p pin_radius, centered at @p center.
 */
template <int dim>
  void pin_box (Triangulation<dim> & /* tria */,
    const Point<2> & /* center */,
    const double /* pitch */,
    const unsigned int /* first_mat_id = 0 */)
  {
    Assert(false, ExcMessage ("Not implemented."));
  }

/**
 *
 */
template <>
  void pin_box (Triangulation<2> & tria,
    const Point<2> & center,
    const double pitch,
    const unsigned int first_mat_id)
  {
    const unsigned int n_verts = 4;
    const unsigned int n_cells = 1;

    const double p = pitch / 2.0;

    const double xx = center(0);
    const double yy = center(1);

    /**
     *  Defining special vertices manually
     *
     *  2---3
     *  |   |
     *  0---1
     */
    const Point<2> vertices[n_verts] =
                                         {
                                           // box
                                           Point<2>(xx - p, yy - p),
                                           Point<2>(xx + p, yy - p),
                                           Point<2>(xx - p, yy + p),
                                           Point<2>(xx + p, yy + p) };

    /**
     *  Define connections following the deal.II standards
     *       3
     *    2-->--3
     *    |     |
     *   0^     ^1
     *    |     |
     *    0-->--1
     *       2
     *
     * This cell is {0, 1, 2, 3}
     */
    const int cell_vertices[n_cells][GeometryInfo<2>::vertices_per_cell] =
          {
              { 0, 1, 2, 3 } };

    std::vector<CellData<2> > cells(n_cells, CellData<2>());

    // creating arrays with vertex info
    for (unsigned int i = 0; i < n_cells; ++i)
    {
      for (unsigned int j = 0;
          j < GeometryInfo<2>::vertices_per_cell;
          ++j)
      {
        cells[i].vertices[j] = cell_vertices[i][j];
      }
      cells[i].material_id = first_mat_id;
    }

    // Create the triangulation with no boundary information
    tria.create_triangulation(
      std::vector<Point<2> >(&vertices[0], &vertices[n_verts]),
      cells,
      SubCellData());
  }

/**
 *
 */
template <int dim>
  void print_pin (const Triangulation<dim> & tria)
  {
    std::ofstream eps_out("grid_pin.eps");
    GridOut gridout;
    gridout.write_eps(tria, eps_out);

    std::ofstream vtk_out("grid_pin.vtk");

    DoFHandler<dim> dof_handler(tria);
    static const FE_DGQ<dim> fe(1);
    dof_handler.distribute_dofs(fe);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    Vector<double> mat_id_cells(tria.n_active_cells());
    unsigned int j = 0;
    typename Triangulation<dim>::active_cell_iterator cell;
    for (cell = tria.begin_active(); cell != tria.end(); ++cell, ++j)
    {
      unsigned int cell_mat = cell->material_id();
      mat_id_cells[j] = static_cast<double>(cell_mat);
    }
    data_out.add_data_vector(mat_id_cells, "material");
    data_out.build_patches();
    data_out.write_vtk(vtk_out);
  }

/**
 *
 */
template <int dim>
  void
  assign_pin_manifolds (Triangulation<dim> & tria,
    const std::vector<Point<2>> pin_center,
    const std::vector<double> & pin_radius,
    const std::vector<Pin_type> & pin_types,
    const unsigned int first_manifold)
  {
    Assert(dim >= 2, ExcMessage("dim must be at least 2."));
    // how many pins?
    const unsigned int n_total_pins = pin_center.size();
    // After extruding the triangulation we have to specify the interior
    // boundaries (manifolds). We have to do it here because
    // "extrude_triangulation" does not take them into account yet.
    for (typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
        cell != tria.end(); ++cell)
    {
      const Point<dim> cell_center = cell->center();
      for (unsigned int i = 0; i < n_total_pins; ++i)
      {

        if (pin_types[i] == Pin_type::pin)
        {
          // we only use x and y coordinates to see if we are inside the pin
          /*
           double distance = std::sqrt(
           std::pow(pin_center[i][0]-cell_center(0),2)
           + std::pow(pin_center[i][1]-cell_center(1),2));
           double cell_diameter = cell->diameter();
           if ((distance + cell_diameter/2.0 > pin_radius[i]) and
           (distance < pin_radius[i]))
           */
          if ((std::abs(pin_center[i][0] - cell_center(0)) < 0.9 * pin_radius[i])
              and
              (std::abs(pin_center[i][1] - cell_center(1)) < 0.9 * pin_radius[i]))
          {
            const unsigned int manifold_ind = first_manifold + i;
            cell->set_all_manifold_ids(manifold_ind);
            /*
             cell->set_manifold_id (manifold_ind);
             if (cell->has_children())
             for (unsigned int c=0; c<cell->n_children(); ++c)
             cell->child(c)->set_all_manifold_ids (manifold_ind);
             // for hexes also set manifold_id of bounding quads and lines
             // Six bonding quads (only four faces, so we use dim = 2 for this).
             for (unsigned int i=0; i< GeometryInfo<2>::faces_per_cell; ++i)
             cell->face(i)->set_manifold_id(manifold_ind);
             // Twelve bounding lines (or four, depending on the dimension
             for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_cell; ++i)
             cell->line(i)->set_manifold_id(manifold_ind);
             */
            if ((std::abs(pin_center[i][0] - cell_center(0))
                 < 0.01 * pin_radius[i])
                and
                (std::abs(pin_center[i][1] - cell_center(1)) < 0.01 * pin_radius[i]))
            {
              cell->set_all_manifold_ids(-1);
            }
          }
        }
      }
    }
  }

/**
 *
 */
template <int dim>
  void refine_assembly (
    Triangulation<dim> & tria,
    const std::string refinement_model,
    const unsigned int n_ref_radial)
  {
    if (refinement_model == "uniform")
      tria.refine_global(n_ref_radial);

    else if (refinement_model == "local")
    {
      const unsigned int max_loc_refinement = 3;
      Assert(max_loc_refinement <= 10, ExcMessage(
        std::string("Here we set max_loc_refinement = 3 \n") +
        std::string("until we deal with deal with hanging nodes\n") ));
      // Here we iterate for the refinements
      for (unsigned int step = 0; step < n_ref_radial; step++)
      {
        // Here we look for the faces over the internal boundary, and tag
        // this cells to be refined. (if two faces belong to the internal
        // manifold, will be isotropically refined, otherwise will be refined
        // only in one direction)
        if (step < max_loc_refinement)
        {
          typename Triangulation<dim>::active_cell_iterator cell;
          for (cell = tria.begin_active(); cell != tria.end(); ++cell)
          {
            std::vector<unsigned int> direction;
            // we only run over the radial faces (thus, dim = 2 always)
            for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
            {
              if (cell->face(f)->at_boundary())
                continue;
              if (cell->material_id() != cell->neighbor(f)->material_id())
              {
                // Here we use 1-f/2 because:
                // f/2 will be 1 if we want to refine faces f={2,3} (x-axis),
                // then 1-f/2 will be 0, and
                // f/2 will be 0 if we want to refine faces f={0,1} (y-axis),
                // then 1-f/2 will be 1.
                direction.push_back(1 - f / 2);
              }
              /*if (cell->face(f)->manifold_id() != dealii::types::manifold_id(-1))
               {
               // Here we use 1-f/2 because:
               // f/2 will be 1 if we want to refine faces f={2,3} (x-axis),
               // then 1-f/2 will be 0, and
               // f/2 will be 0 if we want to refine faces f={0,1} (y-axis),
               // then 1-f/2 will be 1.
               direction.push_back(1-f/2);
               }*/
            }
            if (direction.size() >= 2)
            {
              cell->set_refine_flag(RefinementCase<dim>::cut_xy);
            }
            else if (direction.size() == 1)
            {
              cell->set_refine_flag(RefinementCase<dim>::cut_axis(direction[0]));
            }
            // Added by Toni to refine the central cell
            else if (direction.size() == 0)
            {
              cell->set_refine_flag(RefinementCase<dim>::cut_xy);
            }

            direction.clear();
          }
          // No more than max_loc_refinement anisotropic localized refinements.
          // After that, we do global isotropic refinement.
          tria.execute_coarsening_and_refinement();
        }
        else
          tria.refine_global(1);
      }
    }
    else
      AssertRelease(false, "Error in: refine_assembly(refinement_model)");
  }

/**
 * Calculate the effective radius to pass to the GridGenerator::pin_cell()
 * to maintain the fuel area.
 */
double calculate_effective_radius (const double radius,
  const unsigned int n_refinements)
{
  const double n_edges = std::pow(2, 2 + n_refinements);
  const double alpha = 2 * M_PI / n_edges; // interior angle of the polygon
  const double circunscribed_radious = radius
      * std::sqrt(2 * M_PI / (n_edges * std::sin(alpha)));

  return circunscribed_radious;
}

/**
 *
 */
template <int dim>
  void
  create_assembly (Triangulation<dim> & /* tria_out */,
    std::vector<Triangulation<2> > & /* tria_in */,
    const std::vector<double> & /* z_lengths */,
    const unsigned int /* n_assemblies_plane */)
  {
    Assert(false, ExcNotImplemented());
  }

/**
 *
 */
template <>
  void create_assembly (Triangulation<2> & tria_out,
    std::vector<Triangulation<2> > & tria_in,
    const std::vector<double> & /* z_lengths */,
    const unsigned int /* n_assemblies_plane */)
  {
    if (tria_in.size() > 1)
      merge_triangulations(tria_in, tria_out);
    else
    {
      AssertDimension(tria_in.size(), 1);
      tria_out.copy_triangulation(tria_in[0]);
    }
  }

/**
 *
 */
template <>
  void
  create_assembly (Triangulation<3> & tria_out,
    std::vector<Triangulation<2> > & tria_in,
    const std::vector<double> & z_lengths,
    const unsigned int n_assemblies_plane)
  {
    Assert(z_lengths.size() > 0, ExcMessage("z_lengths must be minimum 1"));

    Triangulation<2, 2> tria_plane;
    if (tria_in.size() > 1)
      merge_triangulations(tria_in, tria_plane);
    else
    {
      AssertDimension(tria_in.size(), 1);
      tria_plane.copy_triangulation(tria_in[0]);
    }

    // Calculate slice_coordinates
    unsigned int n_slices = z_lengths.size() + 1;
    std::vector<double> slice_coordinates;
    slice_coordinates.reserve(n_slices);
    slice_coordinates.push_back(0.0);
    for (unsigned int pln = 0; pln < z_lengths.size(); pln++)
      slice_coordinates.push_back(slice_coordinates.back() + z_lengths[pln]);

    GridGenerator::extrude_triangulation(
      tria_plane,
      slice_coordinates,
      tria_out);

    // Set user indices
    // extrude_triangulation reorders cell in such a way that cell2d 0 is extrudded
    // into cell3d 0, 1, 2...
    unsigned int n_planes = z_lengths.size();
    //const double tol = 1e-4;
    typename Triangulation<3>::active_cell_iterator cell3d;
    typename Triangulation<2>::active_cell_iterator cell2d;
    for (cell2d = tria_plane.begin_active(), cell3d = tria_out.begin_active();
        cell2d != tria_plane.end(); ++cell2d)
    {
      for (unsigned int pln = 0; pln < n_planes; pln++, ++cell3d)
      {
        cell3d->set_user_index(cell2d->user_index() + n_assemblies_plane * pln);
//        std::cout << "Cell " << cell3d << " height " << cell3d->vertex(0)[2] << "  user "
//                  << cell3d->user_index()
//                  << std::endl;
      }
    }
  }

/**
 *
 */
template <int dim>
  void set_manifold_and_refine (Triangulation<dim> & /* tria */,
    const std::vector<Point<2> > & /* pin_center */,
    const std::vector<Pin_type> & /* pin_types */,
    const unsigned int /* n_ref_radial = 0 */,
    const unsigned int /* first_manifold = 10 */,
    const std::string /* refinent_model */)
  {
    Assert(false, ExcNotImplemented());
  }

/**
 *
 */
template <>
  void set_manifold_and_refine (Triangulation<2> & tria,
    const std::vector<Point<2> > & pin_center,
    const std::vector<Pin_type> & pin_types,
    const unsigned int n_ref_radial,
    const types::manifold_id first_manifold,
    const std::string refinement_model)
  {
    const unsigned int n_total_pins = pin_center.size();

    // we associate the internal manifolds with boundary entities
    std::vector<SphericalManifold<2, 2> > manifold;
    for (unsigned int i = 0; i < n_total_pins; ++i)
    {
      const Point<2> aux_center = Point<2>(pin_center[i]);
      SphericalManifold<2, 2> manifold_aux(aux_center);
      manifold.push_back(manifold_aux);
    }

    for (unsigned int i = 0; i < n_total_pins; ++i)
      if (pin_types[i] == Pin_type::pin)
        tria.set_manifold(first_manifold + i,
          manifold[i]);

    refine_assembly(tria, refinement_model, n_ref_radial);

    // Detaching the manifold from the triangulation (so we can return
    // the triangulation with a function).
    for (unsigned int i = 0; i < n_total_pins; ++i)
      if (pin_types[i] == Pin_type::pin)
        tria.reset_manifold(first_manifold + i);
  }

/**
 *
 */
template <>
  void set_manifold_and_refine (Triangulation<3> & tria,
    const std::vector<Point<2>> & pin_center,
    const std::vector<Pin_type> & pin_types,
    const unsigned int n_ref_radial,
    const unsigned int first_manifold,
    const std::string refinement_model)
  {
    const unsigned int n_total_pins = pin_center.size();

    // we associate the internal manifolds with boundary entities
    std::vector<CylindricalManifold<3>> manifold;
    for (unsigned int i = 0; i < n_total_pins; ++i)
    {
      const Point<3> aux_center = Point<3>(pin_center[i](0), pin_center[i](1), 0);
      const Point<3> aux_dir = Point<3>(0, 0, 1);
      CylindricalManifold<3> manifold_aux(aux_dir,
        aux_center);
      manifold.push_back(manifold_aux);
    }

    for (unsigned int i = 0; i < n_total_pins; ++i)
      if (pin_types[i] == Pin_type::pin)
        tria.set_manifold(first_manifold + i, manifold[i]);

    refine_assembly(tria, refinement_model, n_ref_radial);

    // Detaching the manifold from the triangulation (so we can return
    // the triangulation with a function).
    for (unsigned int i = 0; i < n_total_pins; ++i)
      if (pin_types[i] == Pin_type::pin)
        tria.reset_manifold(first_manifold + i);
  }

/**
 *
 */
template <int dim>
  void colorize_rectangular_triangulation (Triangulation<dim> & tria,
    const Point<dim> & p1,
    const Point<dim> & p2,
    double epsilon)
  {
    // run through all faces and check
    // if one of their center coordinates matches
    // one of the corner points. Comparisons
    // are made using an epsilon which
    // should be smaller than the smallest cell
    // diameter.
    typename Triangulation<dim>::face_iterator face = tria.begin_face(),
        endface = tria.end_face();
    for (; face != endface; ++face)
    {
      if (face->at_boundary())
      {
        const Point<dim> center(face->center());
        if (std::abs(center(0) - p1[0]) < epsilon)
          face->set_boundary_id(0);
        else if (std::abs(center(0) - p2[0]) < epsilon)
          face->set_boundary_id(1);
        else if (dim > 1 && std::abs(center(1) - p1[1]) < epsilon)
          face->set_boundary_id(2);
        else if (dim > 1 && std::abs(center(1) - p2[1]) < epsilon)
          face->set_boundary_id(3);
        else if (dim > 2 && std::abs(center(2) - p1[2]) < epsilon)
          face->set_boundary_id(4);
        else if (dim > 2 && std::abs(center(2) - p2[2]) < epsilon)
          face->set_boundary_id(5);
        else
        {
          // triangulation says it is on the boundary,
          // but we could not find on which boundary.
          AssertRelease(false, "Error with colorize_rectangular_triangulation()");
        }
      }
    }
  }

/**
 *
 */
template <int dim>
  void pin_area (const Triangulation<dim> & tria)
  {
    double total_area = 0;
    double pin_area = 0;
    typename Triangulation<dim>::active_cell_iterator cell;
    for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      total_area += cell->measure();
      if (cell->material_id() == 1)
        pin_area += cell->measure();
    }
    std::cout << "area pin   = " << pin_area << std::endl;
    std::cout << "area total = " << total_area << std::endl;
    std::cout << "ratio = " << 100. * pin_area / total_area << " %" << std::endl;
  }

/**
 *
 */
template <int dim>
  void make_composed_geometry (
    InputGeom & geom,
    Triangulation<dim> & tria_result,
    const unsigned int n_ref_radial,
    const unsigned int n_ref_axial,
    const std::string refinement_model)
  {
    assert(dim == 2 or dim == 3);
    unsigned int counter;
    const unsigned int nothing = static_cast<unsigned int>(-1);

    assert(dim == 3 or geom.core.n_planes == 1);

    std::vector<double> z_lengths = geom.core.length;
    const types::material_id inner_mat = 0;
    const types::material_id outter_mat = 1;
    const types::material_id box_mat = 2;

    if (dim == 3)
      for (unsigned pln = 0; pln < geom.planes.size(); pln++)
        AssertRelease(geom.planes[pln].n_nodes[2] == 1,
          "Only implemented with planes[pln].n_nodes[2] == 1");

    // pin_map
    const std::vector<std::vector<unsigned int> > pins_map = geom.planes[0].components;
    const std::vector<unsigned int> n_pins =
          { static_cast<unsigned int>(pins_map.size()),
            static_cast<unsigned int>(pins_map[0].size()) };

    // How many not empty pins do we have
    unsigned int n_total_pins = 0;
    for (unsigned int i = 0; i < n_pins[0]; ++i)
      for (unsigned int j = 0; j < n_pins[1]; ++j)
        if (pins_map[i][j] != nothing)
          n_total_pins++;

    std::vector<Triangulation<2> > tria_pin(n_total_pins);
    Triangulation<2> tria_plane;

    // Now we work with dim = 2 to generate the geometry before extrusion
    std::vector<Point<2> > pin_center(n_total_pins);
    std::vector<double> pin_radius(n_total_pins);
    std::vector<Pin_type> pin_types(n_total_pins);

    // We generate isolated pins in a vector of triangulations, in order
    // to merge them later
    // It is assumed that all pins have the same pitch
    counter = 0;
    for (unsigned int j = 0; j < n_pins[1]; ++j)
    {
      for (unsigned int i = 0; i < n_pins[0]; ++i)
      {
        if (pins_map[i][j] != nothing)
        {
          InputGeom::Pin & pin = geom.pins[pins_map[i][j]];

          switch (pin.pin_type)
          {
            case Pin_type::pin:
              {
              pin_center[counter][0] = double(i) * pin.pitch + pin.pitch / 2;
              pin_center[counter][1] = double(j) * pin.pitch + pin.pitch / 2;
              pin_radius[counter] = calculate_effective_radius(pin.fuel_radius,
                n_ref_radial);
              pin_types[counter] = Pin_type::pin;

              // Make pin triangulation
              if (dim == 2)
                pin_cell(tria_pin[counter],
                  pin_center[counter],
                  pin.pitch,
                  pin_radius[counter],
                  pin.materials[0],
                  pin.materials[1]);
              else
                pin_cell(tria_pin[counter],
                  pin_center[counter],
                  pin.pitch,
                  pin_radius[counter],
                  inner_mat,
                  outter_mat);

              break;
            }
            case Pin_type::box:
              {
              pin_center[counter][0] = double(i) * pin.pitch + pin.pitch / 2;
              pin_center[counter][1] = double(j) * pin.pitch + pin.pitch / 2;
              pin_radius[counter] = 0.0;
              pin_types[counter] = Pin_type::box;

              // Make pin triangulation
              if (dim == 2)
                pin_box(tria_pin[counter],
                  pin_center[counter],
                  pin.pitch,
                  pin.materials[0]);
              else
                pin_box(tria_pin[counter],
                  pin_center[counter],
                  pin.pitch,
                  box_mat);
              break;
            }
            default:
              {
              Assert(false, ExcNotImplemented())
              break;
            }
          }
          counter++;
        }
      }
    }

    for (unsigned int c = 0; c < n_total_pins; c++)
    {
      typename Triangulation<2>::active_cell_iterator cell = tria_pin[c].begin_active();
      for (; cell != tria_pin[c].end(); ++cell)
        cell->set_user_index(c);
    }

    // merge the pins in a unique mesh
    create_assembly<dim>(tria_result, tria_pin, z_lengths, n_total_pins);

    if (dim == 3) // Set materials Id for 3d geometries
    {
      unsigned int pln, i, j, user_index;
      typename Triangulation<dim>::active_cell_iterator cell;
      for (cell = tria_result.begin_active(); cell != tria_result.end(); ++cell)
      {
        user_index = cell->user_index();
        pln = user_index / n_total_pins;
        j = (user_index - pln * n_total_pins) / n_pins[1];
        i = (user_index - pln * n_total_pins - n_pins[1] * j);

        if (geom.planes[geom.core.components[pln]].components[i][j] == nothing)
          AssertRelease(false, "Error in prob_geom.cc");

        InputGeom::Pin & pin =
            geom.pins[geom.planes[geom.core.components[pln]].components[i][j]];

        if (cell->material_id() == inner_mat)
          cell->set_material_id(pin.materials[0]);
        else if (cell->material_id() == outter_mat)
          cell->set_material_id(pin.materials[1]);
        else if (cell->material_id() == box_mat)
          cell->set_material_id(pin.materials[0]);
        else
          AssertRelease(false, "Error setting materials in prob_geom.cc");
      }

      std::vector<unsigned int> user_indices(tria_result.n_active_cells());
      unsigned int idx = 0;
      // Get user indices
      for (cell = tria_result.begin_active(); cell != tria_result.end(); ++cell, idx++)
      {
        user_indices[idx] = cell->user_index();
      }

      // Refine axially
      for (unsigned int ref = 0; ref < n_ref_axial; ref++)
      {
        for (cell = tria_result.begin_active(); cell != tria_result.end(); ++cell)
        {
          cell->set_refine_flag(RefinementCase<dim>::cut_axis(2));
        }
        tria_result.execute_coarsening_and_refinement();
      }

      // Copy back user_indices
      idx = 0;
      typename Triangulation<dim>::cell_iterator cell_it;
      for (cell_it = tria_result.begin(0); cell_it != tria_result.end(); ++cell_it, ++idx)
      {
        cell_it->recursively_set_user_index(user_indices[idx]);
      }

    }
    // TODO Allow that all pins are not equals
    InputGeom::Pin & pin = geom.pins[pins_map[0][0]];
    const double epsilon = 0.001 * pin.pitch;
    Point<dim> p1, p2;
    for (auto d = 0; d < dim; d++)
    {
      p1[d] = 0.0;
      if (d < 2)
        p2[d] = pin.pitch * n_pins[d];
      else
        p2[d] = sum_vector(z_lengths);

    }

    colorize_rectangular_triangulation(tria_result, p1, p2, epsilon);

    // set the different identifiers to the mesh, in order to attach
    // different manifolds to different identifiers
    const unsigned int first_manifold = 10;
    assign_pin_manifolds(tria_result,
      pin_center,
      pin_radius,
      pin_types,
      first_manifold);

    // attach manifolds to the mesh and refine
    set_manifold_and_refine(tria_result,
      pin_center,
      pin_types,
      n_ref_radial,
      first_manifold,
      refinement_model);

    // Set user_indexes through refinement
    typename Triangulation<dim>::cell_iterator cell = tria_result.begin(0);
    for (; cell != tria_result.end(0); ++cell)
    {
      cell->recursively_set_user_index(cell->user_index());
    }
    // pin area
    //pin_area<dim>(tria);
  }

template void make_composed_geometry (
  InputGeom & geom,
  Triangulation<1> &tria,
  const unsigned int n_ref_radial,
  const unsigned int n_ref_axial,
  const std::string refinement_model);
template void make_composed_geometry (
  InputGeom & geom,
  Triangulation<2> &tria,
  const unsigned int n_ref_radial,
  const unsigned int n_ref_axial,
  const std::string refinement_model);
template void make_composed_geometry (
  InputGeom & geom,
  Triangulation<3> &tria,
  const unsigned int n_ref_radial,
  const unsigned int n_ref_axial,
  const std::string refinement_model);
