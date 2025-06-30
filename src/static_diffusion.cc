/**
 * @file   static_diffusion.cc
 * @brief  Implementation of the class StaticDiffusion and the main functions of
 *  the FEMFFUSION program.
 */

#include <deal.II/lac/solver_selector.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_tools.h>
#include <filesystem>

#include "../include/femffusion.h"
#include "../include/static_diffusion.h"
#include "../include/io/prob_geom.h"
#include "../include/io/input_geom.h"
#include "../include/io/printing.h"
#include "../include/matrix_operators/matrix_operators_petsc.h"
#include "../include/eps_solvers/eps_solver.h"

#include <petscsys.h>

#include <map>
#include <set>
#include <algorithm>

using namespace dealii;

/**
 * @brief Constructor of the main class StaticDiffusion.
 * Reads the input file and it reads or builds the grid.
 */
template <int dim, int n_fe_degree>
  StaticDiffusion<dim, n_fe_degree>::StaticDiffusion (
    ParameterHandler &prm,
    std::string input_file,
    const bool verbose,
    const bool silent,
    const bool _to_init,
    std::string modes) :
      comm(PETSC_COMM_WORLD),
      n_mpi_processes(Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
      n_local_cells(numbers::invalid_unsigned_int),
      verbose_cout(std::cout, verbose and this_mpi_process == 0),
      materials_cout(std::cout, false), // No verbose for materials
      cout(std::cout, !silent and this_mpi_process == 0),
      to_init(_to_init),
      n_groups(prm.get_integer("Energy_Groups")),
      fe(QGaussLobatto<1>(n_fe_degree + 1)),
      tria(comm),
      dof_handler(tria),
      T(comm, dof_handler, constraints),
      F(comm, dof_handler, constraints),
      L(comm, dof_handler, constraints),
      G(comm, dof_handler, constraints),
      A(comm, dof_handler, constraints),
      V(comm, dof_handler, constraints),
      input_file(input_file),
      materials(materials_cout),
      perturbation(prm, materials, dof_handler, verbose_cout)
  {
    verbose_cout << "Start of the program " << std::endl;
    AssertRelease(n_fe_degree > 0, "FE can not be 0");

    // General Options
    n_eigenvalues = prm.get_integer("N_Eigenvalues");
    n_refinements = prm.get_integer("N_Refinements");
    refine_y = false;

    materials.transient = prm.get_bool("Transient");

    // Output
    out_file = prm.get("Output_Filename");
    out_flag = prm.get_bool("Output_Flag");
    print_grid_flag = prm.get_bool("Print_Grid_Flag");
    n_out_ref = prm.get_integer("Out_Refinements");

    // Solver
    show_eps_convergence = not silent;
    residual_norm = true;
    solver_type = prm.get("Solver_Type");
    tol_eps = prm.get_double("EPS_Tolerance");
    tol_ksp = prm.get_double("KSP_Tolerance");
    adjoint = prm.get_bool("Adjoint");
    static_ksp_tol = prm.get_bool("Static_KSP_Tolerance");
    spectral_index = prm.get_bool("Spectral_index");
    matrixfree_type = string_to_enum(prm.get("Matrix_Free_Type"));
    p_init = prm.get_bool("P_Init");
    init_type = "sp1";
    std::string refinement_model = prm.get("Refinement_Model");
    t_modes = prm.get("Spatial_Modes");

    // Albedo Factors
    std::string str = prm.get("Albedo_Factors");
    trim(str);
    if (str.length() > 0)
      parse_vector(str, albedo_factors, 0);

    // Get changes in the parameters through the command line
    get_parameters_from_command_line();

    if (to_init == false)
      adjoint = true;

    init_adjoint_from_direct = true;

    if (modes != "none")
      t_modes = "lambda";

    // ---------------------------------------------------------------------------------
    // GEOMETRY SHAPE
    geo_type = prm.get("Geometry_Type");
    verbose_cout << "geo_type  " << geo_type << std::endl;
    bool listen_to_material_id = (geo_type == "Composed");
    if (geo_type != "Composed")
      get_mesh_shape(prm);

    // ---------------------------------------------------------------------------------
    // XSEC
    xs_file = prm.get("XSECS_Filename");
    std::string xs_type = prm.get("XSEC_Type");
    precursors_file = prm.get("PREC_Filename");

    verbose_cout << "xs_file  " << xs_file << std::endl;
    verbose_cout << "xs_type  " << xs_type << std::endl;

    std::vector<unsigned int> geo_ps;
    if (xs_type == "Valkin")
    {
      parse_vector(prm.get("Geometry_Points"), geo_ps, 2 * assem_per_dim[1],
        default_geometry_points(assem_per_dim));
    }

    if (geo_type != "Composed")
    {
      // Parse  Materials Cross sections
      materials.reinit(xs_file, xs_type, n_groups, assem_per_dim,
        n_assemblies, listen_to_material_id, geo_ps, precursors_file);
    }

    // ---------------------------------------------------------------------------- //
    // GEOMETRY
    // Mesh File
    if (geo_type == "Unstructured")
    {
      std::string mesh_file = prm.get("Mesh_Filename");
      verbose_cout << "mesh_file " << mesh_file << std::endl;

      verbose_cout << "parsing Boundary_Conditions... " << std::flush;
      parse_vector(prm.get("Boundary_Conditions"), boundary_conditions);
      if (verbose_cout.is_active())
        print_vector(boundary_conditions, false);
      verbose_cout << "Done!" << std::endl;

      // Get the z information to extrude
      assembly_pitch.resize(dim);
      if (dim >= 3)
        parse_vector(prm.get("Cell_Pitch_z"), assembly_pitch[2],
          assem_per_dim[2]);
      // TODO Make the possibility to read a 3D gmsh mesh.

      get_unstructured_grid(mesh_file);

    }
    else if (geo_type == "Hexagonal")
    {
      std::string mesh_file = prm.get("Mesh_Filename");
      verbose_cout << "mesh_file " << mesh_file << std::endl;

      verbose_cout << "parsing assembly_pitch... " << std::endl;
      // Pin Pitch definition
      assembly_pitch.resize(dim);
      parse_vector(prm.get("Cell_Pitch_x"), assembly_pitch[0],
        assem_per_dim[0]);
      if (dim >= 2)
        parse_vector(prm.get("Cell_Pitch_y"), assembly_pitch[1],
          assem_per_dim[1]);
      if (dim >= 3)
        parse_vector(prm.get("Cell_Pitch_z"), assembly_pitch[2],
          assem_per_dim[2]);

      get_unstructured_grid(mesh_file);

      // Check every the mesh
      AssertRelease(n_assemblies == tria.n_cells(0) / 3,
        "Not valid Hexagonal Grid");

      // Hexagonal meshes have only one boundary indicator
      verbose_cout << "parsing Boundary_Conditions... " << std::flush;
      parse_vector(prm.get("Boundary_Conditions"), boundary_conditions);
      if (verbose_cout.is_active())
        print_vector(boundary_conditions, false);

      verbose_cout << " Done!" << std::endl;

    }
    else if (geo_type == "Composed")
    {
      n_refinements_radial = prm.get_integer("N_Refs_Radial");
      n_refinements_axial = prm.get_integer("N_Refs_Axial");

      std::string tria_file = prm.get("Triangulation_Filename");
      std::string vec_file = tria_file + ".vec";

      AssertRelease(xs_type == "XML", "Composed geo_type needs XML xs_type ");
      std::string geom_file = prm.get("Geometry_Filename");
      verbose_cout << "geom_file " << geom_file << std::endl;

      verbose_cout << "load geometry file... " << std::flush;
      InputGeom input_geometry;
      input_geometry.load(geom_file);
      verbose_cout << "Done!" << std::endl;

      // Fill assem_per_dim
      assem_per_dim.resize(3, 1);
      for (unsigned int d = 0; d < dim; ++d)
      {
        if (d < 2)
          assem_per_dim[d] = input_geometry.planes[0].n_nodes[d];
        else
          assem_per_dim[d] = input_geometry.core.n_planes;
      }
      n_assemblies = assem_per_dim[0] * assem_per_dim[1] * assem_per_dim[2];

      materials.reinit(xs_file, xs_type, n_groups, assem_per_dim,
        n_assemblies, listen_to_material_id, geo_ps, precursors_file);

      // Load Triangulation
      if (fexists(tria_file))
      {
        verbose_cout << "load_triangualtion... " << std::flush;
        GridIn<dim> gridin;
        std::ifstream f(tria_file.c_str());
        gridin.attach_triangulation(tria);
        gridin.read_ucd(f);

        std::vector<unsigned int> user_indices;
        parse_vector_in_file(vec_file, "User_indices", user_indices,
          tria.n_active_cells(), tria.n_active_cells());

        unsigned int i = 0;
        typename Triangulation<dim>::active_cell_iterator cell =
                                                                 tria.begin_active();
        for (; cell != tria.end(); ++cell, ++i)
        {
          cell->set_user_index(user_indices[i]);
        }

        verbose_cout << "Done!" << std::endl;
      }
      else
      {
        AssertRelease(n_refinements == 0,
          "Not valid n_refinements for geo_type == Composed");

        verbose_cout << "make_composed_geometry... " << std::flush;
        make_composed_geometry(input_geometry, tria, n_refinements_radial,
          n_refinements_axial, refinement_model);
        verbose_cout << "Done!" << std::endl;
      }
      // Set Geometry Matrix
      verbose_cout << "  reading geometry_matrix... " << std::flush;
      std::vector<unsigned int> geo_ps = default_geometry_points(
        assem_per_dim);
      materials.set_geometry_matrix(assem_per_dim, geo_ps);
      verbose_cout << " Done!" << std::endl;

      // Copy Boundary Conditions
      boundary_conditions = input_geometry.core.boundary[0];
      for (unsigned int d = 1; d < dim; ++d)
        boundary_conditions.insert(boundary_conditions.end(),
          input_geometry.core.boundary[0].begin(),
          input_geometry.core.boundary[0].end());
      AssertRelease(prm.get("Boundary_Conditions") == "",
        "In Composed Geometry_Type is the boundary condition must be set in the geo.xml file");

      // Save Triangulation
      if (not tria_file.empty() and not fexists(tria_file))
      {
        verbose_cout << "save_triangualtion... " << std::flush;
        AssertRelease(tria.has_hanging_nodes() == false,
          "The triangulation must not have hanging nodes \n"
            "and this triangulation has them");
        std::ofstream tria_file_stream(tria_file.c_str(), std::ios::out);
        GridOut grid_out;
        GridOutFlags::Ucd gridout_flags(true, true, true);
        grid_out.set_flags(gridout_flags);
        grid_out.write_ucd<dim>(tria, tria_file_stream);

        std::vector<unsigned int> user_indices;
        user_indices.reserve(tria.n_active_cells());
        unsigned int i = 0;
        typename Triangulation<dim>::active_cell_iterator cell =
                                                                 tria.begin_active();
        for (; cell != tria.end(); ++cell, ++i)
        {
          user_indices.push_back(cell->user_index());
        }

        std::ofstream vec_file_stream(vec_file.c_str(), std::ios::out);
        print_vector_in_file(user_indices, vec_file_stream,
          "User_indices\n", false);
        verbose_cout << " Done!" << std::endl;
      }
    }
    else if (geo_type == "Rectangular")
    {
      verbose_cout << "parsing assembly_pitch..." << std::endl;
      // Pin Pitch definition
      assembly_pitch.resize(dim);
      parse_vector(prm.get("Cell_Pitch_x"), assembly_pitch[0],
        assem_per_dim[0]);
      if (dim >= 2)
        parse_vector(prm.get("Cell_Pitch_y"), assembly_pitch[1],
          assem_per_dim[1]);
      if (dim >= 3)
        parse_vector(prm.get("Cell_Pitch_z"), assembly_pitch[2],
          assem_per_dim[2]);

      if (dim >= 1 and verbose_cout.is_active())
      {
        verbose_cout << "Cell_Pitch_x " << std::flush;
        print_vector(assembly_pitch[0]);
      }
      if (dim >= 2 and verbose_cout.is_active())
      {
        verbose_cout << "Cell_Pitch_y " << std::flush;
        print_vector(assembly_pitch[1]);
      }
      if (dim >= 3 and verbose_cout.is_active())
      {
        verbose_cout << "Cell_Pitch_z " << std::flush;
        print_vector(assembly_pitch[2]);
      }

      // Build the rectangular grid
      make_rectangular_grid();
      verbose_cout << "parsing Boundary_Conditions..." << std::flush;
      parse_vector(prm.get("Boundary_Conditions"), boundary_conditions,
        2 * dim);
      if (verbose_cout.is_active())
        print_vector(boundary_conditions, false);
      verbose_cout << " Done!" << std::endl;

    }
    else
      AssertRelease(false, "Unrecognized geo_type: " + geo_type);

    // ---------------------------------------------------------------------------------

    adjoint = false;

    // ---------------------------------------------------------------------------------
    type_perturbation = prm.get("Type_Perturbation");

    if (prm.get("Bar_Filename") != "no.bar" or type_perturbation == "Rods")
    {
      std::string bar_file = prm.get("Bar_Filename");
      AssertRelease(type_perturbation == "Rods",
        "Type_Perturbation must be Rods");
      AssertRelease(bar_file != "no.bar", "It is necessary a bar file");
      verbose_cout << "parsing bar_file... " << bar_file << std::flush;
      perturbation.parse_bar_file(bar_file);
      verbose_cout << "Done!" << std::endl;
    }

    // ---------------------------------------------------------------------------------
    bool save_static = prm.get_bool("Save_Static");
    std::string static_file = prm.get("STA_Filename");

    if (to_init == false)
    {
      if (static_file == "none.sta" or save_static == true)
        run(prm);

      if (static_file != "none.sta" and save_static == false)
        load_static_calculation(static_file);

      if (save_static == true)
        save_static_calculation(static_file);
    }

  }

/**
 * @brief Get the information about the Mesh Size and and Shape.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::get_mesh_shape (ParameterHandler &prm)
  {

    verbose_cout << "parsing the mesh size..." << std::flush;
    parse_vector(prm.get("Mesh_Size"), assem_per_dim, dim);
    verbose_cout << "Done!" << std::endl;

    // Complete assem_per_dim with 1 if dim is not 3.
    assem_per_dim.reserve(3);
    if (dim < 3)
      assem_per_dim.push_back(1);
    if (dim < 2)
      assem_per_dim.push_back(1);

    // Geometry of the reactor
    if (prm.get("Geometry_Matrix") != "")
    {
      verbose_cout << "parsing Geometry_Matrix... " << std::flush;
      materials.set_geometry_matrix(assem_per_dim,
        prm.get("Geometry_Matrix"));
      AssertRelease(prm.get("Geometry_Points") == "",
        "Geometry_Matrix and Geometry_Points cannot be defined at the same time");
      verbose_cout << "Done!" << std::endl;
    }
    else // if geometry point exists and default complete rectangular problem
    {
      verbose_cout << "parsing Geometry_Points... " << std::flush;
      std::vector<unsigned int> geo_ps;
      parse_vector(prm.get("Geometry_Points"), geo_ps, 2 * assem_per_dim[1],
        default_geometry_points(assem_per_dim));

      materials.set_geometry_matrix(assem_per_dim, geo_ps);
      verbose_cout << " Done!" << std::endl;
    }

    // Calculation n_assemblies
    const std::vector<std::vector<unsigned int> > &geometry_matrix =
        materials.get_geometry_matrix();
    n_assemblies = 0;
    for (unsigned int r = 0; r < assem_per_dim[1]; ++r)
      for (unsigned int c = 0; c < assem_per_dim[0]; ++c)
        if (geometry_matrix[r][c] != 0)
          n_assemblies++;
    n_assemblies *= assem_per_dim[2];
    verbose_cout << "n_assemblies: " << n_assemblies << std::endl;
  }

/**
 * @brief Construct a rectangular grid.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::make_rectangular_grid ()
  {
    verbose_cout << "get_materials_table..." << std::flush;
    Table<dim, types::material_id> materials_table;
    materials.get_materials_table(materials_table, assem_per_dim);
    verbose_cout << " Done!" << std::endl;

    // Make the grid
    verbose_cout << "Construction of the mesh..." << std::flush;
    Point<dim> p1 = (
                    (dim == 1) ?
                                 Point<dim>(0.0) :
                                 ((dim == 2) ?
                                               Point<dim>(0.0, 0.0) :
                                               Point<dim>(0.0, 0.0, 0.0)));
    if (dim == 1)
      // In 1D it is not implemented the colorize=true flag.
      GridGenerator::subdivided_hyper_rectangle(tria, assembly_pitch, p1,
        materials_table, false);
    else
      GridGenerator::subdivided_hyper_rectangle(tria, assembly_pitch, p1,
        materials_table, true);

    tria.refine_global(n_refinements);

    if (refine_y)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
                                                                    tria.begin_active();
          cell != tria.end(); ++cell)
      {
        // User_index is the assembly number
        cell->set_refine_flag(cut_axis<dim>(2));
      }

      tria.prepare_coarsening_and_refinement();
      tria.execute_coarsening_and_refinement();
    }

    typename Triangulation<dim>::cell_iterator cell = tria.begin(0);
    for (unsigned int i = 0; cell != tria.end(0); ++cell, ++i)
    {
      // User_index is the assembly number
      cell->recursively_set_user_index(i);
    }
    verbose_cout << " Done!" << std::endl;
  }

/**
 * @brief Get an unstructured grid from a Gmsh file (.msh). Only 2D/3D grids are accepted.
 * The mesh file must contain the materials id as Physical surfaces/volumes id and
 * the boundary id as Physical line/surfaces. Only quadrilateral/hexahedral meshes are accepted.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::get_unstructured_grid (
    const std::string &mesh_file)
  {
    AssertRelease(dim == 2 or dim == 3,
      "Only valid dim 2 and 3 with unstructured grids");

    if (dim == 2)
    {
      Assert(fexists(mesh_file), ExcMessage("Mesh_file doesn't exist"));
      GridIn<dim> gridin2;
      std::ifstream f(mesh_file.c_str());
      gridin2.attach_triangulation(tria);
      gridin2.read_msh(f);
      tria.refine_global(n_refinements);

      // Set the pertinent Rod (user_index) for each cell.
      typename Triangulation<dim>::active_cell_iterator cell =
                                                               tria.begin_active();
      for (; cell != tria.end(); ++cell)
      {
        // User_index is the assembly number
        // material_id should be set in the Gmsh file.
        // User_index is the assembly number
        cell->set_user_index(cell->material_id() - 1);
        verbose_cout << "Cell: " << cell << "  subdomain: "
        << (cell->material_id() - 1)
        << std::endl;

      }
      verbose_cout << "Done!" << std::endl;
    }
    if (dim == 3)
    {
      Assert(fexists(mesh_file), ExcMessage("mesh_file doesn't exist"));
      verbose_cout << "Extruding..." << std::flush;

      Triangulation<2, 2> tria2d;
      GridIn<2> gridin2;
      std::ifstream f(mesh_file.c_str());
      gridin2.attach_triangulation(tria2d);
      gridin2.read_msh(f);

      // Read the z data from somewhere
      double height = sum_vector(assembly_pitch[2]);
      unsigned int n_slices = assem_per_dim[2] + 1;

      extrude_triangulation(tria2d, n_slices, height, tria, true); // TODO Changemat
      tria.refine_global(n_refinements);

      // Set the pertinent Rod (user_index) for each cell.
      Triangulation<3>::active_cell_iterator cell = tria.begin_active();
      for (; cell != tria.end(); ++cell)
      {
        // User_index is the assembly number
        cell->set_user_index(cell->material_id() - 1);
        //verbose_cout << "Cell: " << cell << "  subdomain: " << (cell->material_id() - 1)
        //             << std::endl;
      }
      verbose_cout << "Done!" << std::endl;
    }
  }

/**
 * @brief It uses Petsc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::get_parameters_from_command_line ()
  {
    // Booleans
    get_bool_from_options("-adjoint", adjoint);
    get_bool_from_options("-mk", adjoint);
    get_bool_from_options("-out_flag", out_flag);
    get_bool_from_options("-print_grid", print_grid_flag);
    get_bool_from_options("-show_eps_convergence", show_eps_convergence);
    get_bool_from_options("-residual_norm", residual_norm);
    get_bool_from_options("-spectral_index", spectral_index);
    get_bool_from_options("-static_ksp_tol", static_ksp_tol);
    get_bool_from_options("-p_init", p_init);
    get_bool_from_options("-refine_y", refine_y);

    // Integers
    get_uint_from_options("-n_refinements", n_refinements);
    get_uint_from_options("-n_eigenvalues", n_eigenvalues);
    get_uint_from_options("-n_out_ref", n_out_ref);

    // Reals
    get_double_from_options("-tol_eps", tol_eps);
    get_double_from_options("-tol_ksp", tol_ksp);

    // String
    get_string_from_options("-solver_type", solver_type);
    get_string_from_options("-out_file", out_file);
    get_string_from_options("-init_type", init_type);
    get_string_from_options("-modes", t_modes);

    // Enum
    get_enum_from_options("-matrixfree_type", matrixfree_type);

    lower_case(solver_type);
    lower_case(init_type);

    //Others
    // -eps_ncv 2
    // -pc_factor_levels 1
    // -pc_factor_mat_ordering_type rcm
  }

template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::show_cells ()
  {

    bool show_materials = false;
    bool show_boundary = false;
    unsigned int mat;

    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    if (!show_materials and !show_boundary)
      return;

    verbose_cout << "   Show Cells: " << std::endl;
    std::vector<unsigned int> mat_vec = materials.get_materials_vector();
    print_vector(mat_vec);
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        mat = materials.get_material_id<dim>(cell);
        if (show_materials)
        {
          cout << "Cell " << cell;
          cout << " of id " << cell->user_index();
          cout << " have material " << mat + 1;
          cout << std::endl;
          cout << " Sigma_A " << materials.get_sigma_r(0, mat)
               << std::endl;
        }
        // Take care of albedo Boundary Conditions: Boundary integral
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
          {
            types::boundary_id boundary_id =
                                             cell->face(f)->boundary_id();
            AssertIndexRange(boundary_id, boundary_conditions.size());
            if (show_boundary)
            {
              cout << " Cell " << cell->user_index() + 1 << " face "
              << f
              << " with mat "
              << materials.get_material_id<dim>(cell) + 1
              << " at boundary " << int(boundary_id)
              << " with bc conditions "
              << boundary_conditions[boundary_id]
              << std::endl;
            }
          }
      }
  }

/**
 * Create/Allocate the structures that hold the degrees of freedom.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::make_dofs ()
  {

    // Relative to the DOFS
    dof_handler.distribute_dofs(fe);

    verbose_cout << " dof_handler.distribute_dofs(fe); " << std::endl;
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
      local_dofs_vector[g] = locally_owned_dofs;

    n_local_cells = GridTools::count_cells_with_subdomain_association(tria,
      tria.locally_owned_subdomain());

    // Set the Dirlichet boundary conditions of 0 value.
    constraints.clear();
    for (unsigned int c = 0; c < boundary_conditions.size(); c++)
      if (boundary_conditions[c] == 0)
      {
        DoFTools::make_zero_boundary_constraints(dof_handler, c,
          constraints);
      }
    constraints.close();

    // Make the Sparsity Pattern and allocate the Matrices
    n_dofs = dof_handler.n_dofs();
    n_cells = tria.n_active_cells();

    phi.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      phi[eig].reinit(local_dofs_vector, comm);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::load_static_calculation (
    std::string &file)
  {
    verbose_cout << "   making dofs..." << std::flush;
    make_dofs();
    verbose_cout << "  Done! " << std::endl;

    cout << "   Static file read: " << std::endl;
    AssertRelease(fexists(file), "STA_FILENAME doesn't exist");
    std::ifstream input(file.c_str(), std::ios::in);
    std::string keyword;
    std::vector<double> phi_critic;
    unsigned int n_dofs_total;

    // for every line
    for (std::string line; getline(input, line);)
    {
      std::istringstream iss(line);
      keyword.clear();
      iss >> keyword;

      if (is_commentary(keyword))
        continue;
      // First definition Material and XSecs:
      else if (keyword == "keff")
      {
        get_new_valid_line(input, line);
        std::istringstream iss(line);
        cout << "     keff: " << std::flush;
        iss >> keff;
        AssertRelease(!iss.fail(), "keff must be defined in STA file!");
        cout << keff << std::endl;
      }
      else if (keyword == "phi")
      {
        iss >> n_dofs_total;
        AssertRelease(!iss.fail(),
          "The size of the vector phi must be after phi in STA file!");

        get_new_valid_line(input, line); // Next line
        parse_vector(line, phi_critic, n_dofs_total);

        cout << "     phi: " << phi_critic[0] << " ... "
             << phi_critic.back()
             << std::endl;
      }
    }
    // Make critical
    materials.make_critical(keff);

    // Copy phi_critic to PETScWrappers::MPI::Vector
    phi[0].reinit(n_groups, comm, n_dofs, n_dofs);
    for (unsigned int i = 0; i < n_dofs_total; i++)
      phi[0][i] = phi_critic[i];
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::save_static_calculation (
    std::string &file)
  {
    // Erase the content of output file
    std::ofstream out(file.c_str(), std::ios::out);

    const unsigned int precision = 15;
    print_in_file(eigenvalues[0], out, "keff\n", precision);

    print_vector_in_file(phi[0], out,
      "phi " + num_to_str(n_dofs * n_groups) + "\n", true, precision);

    out.close();
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::assemble_system_lambda ()
  {
    // Allocate and assemble block matrices
    T.reinit(materials, boundary_conditions, albedo_factors, matrixfree_type);
    F.reinit(materials, matrixfree_type);

    // Print matrices if it is needed
    print_matrices();

    memory_consumption = T.memory_consumption() + F.memory_consumption();

    cout << "   Memory consumption of matrix elements "
         << memory_consumption * 1e-6
         << " MB" << std::endl;

  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::assemble_system_gamma ()
  {
    // Allocate and assemble block matrices
    // CHANGE
    L.reinit(materials, boundary_conditions, albedo_factors, matrixfree_type);
    G.reinit(materials, matrixfree_type);

    memory_consumption = L.memory_consumption() + G.memory_consumption();

    cout << "   Memory consumption of matrix elements "
         << memory_consumption * 1e-6
         << " MB" << std::endl;

  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::assemble_system_alpha ()
  {

    // Allocate and assemble block matrices
    A.reinit(materials, boundary_conditions, albedo_factors, matrixfree_type);
    V.reinit(materials, matrixfree_type);

    memory_consumption = A.memory_consumption() + V.memory_consumption();

    cout << "   Memory consumption of matrix elements "
         << memory_consumption * 1e-6
         << " MB" << std::endl;

  }

/**
 * @brief Print matrices in the given file in a matlab way if a file is specified
 *  with '-print_matrices_matlab' command line option.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode StaticDiffusion<dim, n_fe_degree>::print_matrices ()
  {

    //PetscErrorCode ierr;

    std::string print_matrices_matlab;

    get_string_from_options("-print_matrices_matlab", print_matrices_matlab);

    if (!print_matrices_matlab.empty())
    {
      AssertRelease(matrixfree_type == full_allocated,
        "-print_matrices_matlab must be used with -allocate_matrices option");

      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          print_matrix_in_matlab(T.block(g1, g2),
            "T_" + num_to_str(g1 + 1) + num_to_str(g2 + 1), out, 8);
          print_matrix_in_matlab(F.block(g1, g2),
            "F_" + num_to_str(g1 + 1) + num_to_str(g2 + 1), out, 8);
        }

      out << "T= [";
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
      {
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          out << "T_" + num_to_str(g1 + 1) + num_to_str(g2 + 1) << " ";
        }
        out << ";" << std::endl;
      }
      out << "];";

      out << "F= [";
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
      {
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          out << "F_" + num_to_str(g1 + 1) + num_to_str(g2 + 1) << " ";
        }
        out << ";" << std::endl;
      }
      out << "];";

      out.close();
    }

    return 0;
  }

/**
 * @brief Solve the eigenvalue problem
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::solve_eps ()
  {
    EPSSolver<dim, n_fe_degree> *solver = NULL;

    if (t_modes == "lambda")
      solver = new EPSSolver<dim, n_fe_degree>(solver_type, T, F,
        n_eigenvalues, timer, show_eps_convergence,
        verbose_cout.is_active(), tria, dof_handler, fe);
    else if (t_modes == "alpha")
      solver = new EPSSolver<dim, n_fe_degree>(solver_type, A, V,
        n_eigenvalues, timer, show_eps_convergence,
        verbose_cout.is_active(), tria, dof_handler, fe);
    else if (t_modes == "gamma")
      solver = new EPSSolver<dim, n_fe_degree>(solver_type, L, G,
        n_eigenvalues, timer, show_eps_convergence,
        verbose_cout.is_active(), tria, dof_handler, fe);

    // Select some solver options
    if (solver_type == "slepc_2g" or solver_type == "ks"
        or solver_type == "slepc_7g")
      p_init = false;
    solver->p_init = p_init;
    solver->to_init = to_init;
    solver->tol_ksp = tol_ksp;
    solver->tol_eps = tol_eps;
    solver->init_type = "multilevel-sp1";
    solver->input_file = input_file;
    solver->n_groups = n_groups;
    solver->equations = "diffusion";
    solver->out_file = out_file;
    solver->adjoint = adjoint;

    solver->solve(eigenvalues, phi);

    if (t_modes == "alpha")
    {
      for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
        eigenvalues[eig] = 1.0 / eigenvalues[eig];
    }

    if (adjoint)
    {
      phi_adj = solver->phi_adjoint;
    }

    // Print the Eigenvalues calculated
    cout << std::setprecision(6) << "   The eigenvalues are:" << "     Time = "
         << timer.cpu_time()
         << " s." << std::endl;
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      cout << "      K" << i << " : " << eigenvalues[i] << std::endl;

    delete solver;

    return;
  }

/**
 * @brief Normalize the problem to to mean neutron density power equal 1.
 * Also, calculate and print the mean values per assembly.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::postprocess ()
  {

    // Make reactor critical
    //	materials.make_critical(eigenvalues[0]);

    std::vector<double> norm = std::vector<double>(n_eigenvalues);

    phi_serial.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      phi_serial[eig].reinit(n_groups, n_dofs);
      for (unsigned int g = 0; g < n_groups; g++)
        phi_serial[eig].block(g) = phi[eig].block(g);
    }

    if (this_mpi_process == 0)
    {
      // Create Folder if Does not exist
      std::size_t pos = out_file.find_last_of('/');
      std::string folderPath = (pos != std::string::npos) ? out_file.substr(0, pos) : "";
      if (!std::filesystem::exists(folderPath))
        std::filesystem::create_directory(folderPath);
      // Create and erase the content of output file
      std::ofstream out(out_file.c_str(), std::ios::out);

      // Print the eigenvalues in the outFile
      print_logo(out);
      print_vector_in_file(eigenvalues, out, "The Eigenvalues are: \n");
      out << "Problem File: " << input_file << "\n";
      print_in_file(timer.cpu_time(), out, "CPU Time: ", 4);
      out << "Equations: " << "DIFFUSION" << "\n";
      if (geo_type == "Composed")
      {
        print_in_file(n_refinements_radial, out, "Radial Refinements: "); // @suppress("Invalid arguments")
        print_in_file(n_refinements_axial, out, "Axial Refinements: "); // @suppress("Invalid arguments")
      }
      else
        print_in_file(n_refinements, out, "Global Refinements: "); // @suppress("Invalid arguments")

      print_in_file(n_fe_degree, out, "Degree of FE: ", 1);
      print_in_file(dof_handler.n_dofs(), out, "DoFs per Group: ", 1);
      print_in_file(dof_handler.n_dofs() * n_groups, out, "Total DoFs: ", 1);
      print_in_file(tria.n_active_cells(), out, "Number of active cells: ",
        1);
      print_vector_in_file(assem_per_dim, out, "Mesh Size: ", true);
      print_in_file("", out, "\n", 1);
      out.close();

      std::vector<IndexSet> locally_relevant_dofs(n_groups);
      for (unsigned int g = 0; g < n_groups; ++g)
        DoFTools::extract_locally_relevant_dofs(dof_handler,
          locally_relevant_dofs[g]);

      // Initialize all that  is needed to iterate over dofs and cells
      QGauss<dim> quadrature_formula(n_fe_degree + 1);

      FEValues<dim> fe_values(fe, quadrature_formula,
        update_values | update_quadrature_points
        | update_volume_elements
        | update_JxW_values);

      unsigned int n_q_points = quadrature_formula.size();
      const unsigned int n_cells_out = n_assemblies;
      double power_cell = 0;
      double phi_cell;
      double sigma_f;
      double vol;
      unsigned int index;

      // Initialize and resize the vectors where it is stored the solution
      std::vector<double> volume_per_assembly;
      power_per_assembly.resize(n_eigenvalues,
        std::vector<double>(n_cells_out, 0.0));

      phi_per_assembly.resize(n_eigenvalues,
        std::vector<std::vector<double> >(n_groups,
          std::vector<double>(n_cells_out, 0.0)));

      std::vector<double> local_phi(n_q_points);

      // For all defined eigenvalues
      for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      {
        double volume = 0.0;
        norm[eig] = 0.0;

        // Initialize Values
        volume_per_assembly.assign(n_cells_out, 0.0);

        // Iterate over every cell
        typename DoFHandler<dim>::active_cell_iterator cell =
                                                              dof_handler.begin_active(),
            endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
          fe_values.reinit(cell);
          power_cell = 0;
          index = cell->user_index();
          for (unsigned int g = 0; g < n_groups; ++g)
          {
            sigma_f = materials.get_sigma_f(g,
              materials.get_material_id<dim>(cell));
            fe_values.get_function_values(phi_serial[eig].block(g),
              local_phi);

            phi_cell = 0.0;
            for (unsigned int q = 0; q < n_q_points; q++)
              phi_cell += local_phi[q] * fe_values.JxW(q);

            power_cell += sigma_f * phi_cell;

            phi_per_assembly[eig][g][index] += phi_cell;
          }

          vol = cell->measure();
          volume += vol;
          norm[eig] += std::abs(power_cell);
          power_per_assembly[eig][index] += power_cell;
          volume_per_assembly[index] += vol;

        }
        norm[eig] /= volume;

        // Normalize the values of the power and fluxes per cell
        normalize_vector(power_per_assembly[eig], norm[eig]);
        for (unsigned int g = 0; g < n_groups; ++g)
          normalize_vector(phi_per_assembly[eig][g], norm[eig]);

        // Normalize the values of the power and fluxes per assembly
        normalize_vector(power_per_assembly[eig], volume_per_assembly);
        for (unsigned int g = 0; g < n_groups; ++g)
          normalize_vector(phi_per_assembly[eig][g], volume_per_assembly);

        // Normalize the result to make the first value positive:
        bool change_flag = false;
        unsigned int i = 0;
        while (change_flag == false)
        {
          if (power_per_assembly[eig][i] < -0.0001)
          {
            normalize_vector(power_per_assembly[eig], -1);
            for (unsigned int g = 0; g < n_groups; ++g)
            {
              normalize_vector(phi_per_assembly[eig][g], -1);
            }
            norm[eig] *= -1.0;
            change_flag = true;
          }
          else if (power_per_assembly[eig][i] > +0.0001)
          {
            change_flag = true;
          }
          i++;
        }

        //	Normalize also the Fluxes passed to the .vtk
        for (unsigned int g = 0; g < n_groups; ++g)
          normalize_vector(phi_serial[eig].block(g), norm[eig]);

        // ------------------------------------------------------------//
        // ------------------------------------------------------------//
        // Calculate the axial power distribution
        std::vector<double> power_axial;
        std::vector<std::vector<double> > phi_axial;
        std::vector<double> volume_per_plane;
        unsigned int plane;

        if (dim == 3)
        {
          power_axial.resize(assem_per_dim[2], 0.0);
          phi_axial.resize(n_groups,
            std::vector<double>(assem_per_dim[2]));
          volume_per_plane.resize(assem_per_dim[2], 0.0);

          for (unsigned int i = 0; i < power_per_assembly[eig].size();
              i++)
            if (power_per_assembly[eig][i] > 1e-5)
            {
              plane = materials.plane(i);
              power_axial[plane] += power_per_assembly[eig][i]
                                    * volume_per_assembly[i];
              volume_per_plane[plane] += volume_per_assembly[i];

              // for (unsigned int g = 0; g < n_groups; ++g)
              //   phi_axial[g][plane] += phi_per_assembly[eig][g][i] * volume_per_assembly[i];
            }

          // Normalize axial power
          normalize_vector(power_axial, volume_per_plane);
          for (unsigned int g = 0; g < n_groups; ++g)
            normalize_vector(phi_axial[g], volume_per_plane);
        }

        // ------------------------------------------------------------//
        // ------------------------------------------------------------//
        // Print the Powers Distribution and the Fluxes
        const unsigned int presicion = 6;
        std::ofstream out_stream(out_file.c_str(), std::ios::app);
        print_cell_distribution_in_file(dim, power_per_assembly[eig],
          assem_per_dim, out_stream, materials, "Neutron Power\n",
          false, presicion);

        if (dim == 3)
          print_vector_in_file(power_axial, out_stream,
            "Axial power distribution:\n");

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          print_cell_distribution_in_file(dim, phi_per_assembly[eig][g],
            assem_per_dim, out_stream, materials,
            "Group " + num_to_str(g + 1) + " flux\n", false,
            presicion);
          if (dim == 3)
            print_vector_in_file(phi_axial[g], out_stream,
              "Group " + num_to_str(g + 1)
              + "  axial flux distribution:\n");
        }

        out_stream << "\n\n";
        out_stream.close();

        power_axial.clear();
        phi_axial.clear();
        volume_per_plane.clear();
      }
    }

    MPI_Bcast(&norm, n_eigenvalues, MPIU_REAL, 0, comm);

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      for (unsigned int g = 0; g < n_groups; g++)
        phi[eig].block(g) /= norm[eig];

  }

/**
 * @brief Function that creates the output files .vtk.
 * It should be extended to create the power.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::output_results () const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    std::string filename_vtk = out_file + ".vtk";

    // Find the last '/' character
    std::size_t pos = filename_vtk.find_last_of('/');
    // Extract the folder path
    std::string folderPath =
                             (pos != std::string::npos) ?
                                                          filename_vtk.substr(0, pos) :
                                                          "";
    if (!std::filesystem::exists(folderPath))
    {
      std::filesystem::create_directory(folderPath);
      cout << "Folder created: " << folderPath << std::endl;
    }

    std::vector<Vector<double> > power_per_cell(n_eigenvalues,
      Vector<double>(tria.n_active_cells()));
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      power_per_cell[eig].reinit(tria.n_active_cells());

    std::vector<std::vector<Vector<double> > > flux_per_cell(n_eigenvalues,
      std::vector<Vector<double> >(n_groups,
        Vector<double>(tria.n_active_cells())));
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_groups; ++g)
        flux_per_cell[eig][g].reinit(tria.n_active_cells());

    // Power and Flux per assembly
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      for (unsigned int g = 0; g < n_groups; ++g)
        data_out.add_data_vector(phi_serial[eig].block(g),
          "phi_g" + num_to_str(g + 1) + "_eig_"
          + num_to_str(eig + 1));

      unsigned int c = 0;
      typename Triangulation<dim>::active_cell_iterator cell =
                                                               tria.begin_active();
      for (; cell != tria.end(); ++cell, ++c)
      {
        power_per_cell[eig][c] =
                                 power_per_assembly[eig][cell->user_index()];
        for (unsigned int g = 0; g < n_groups; ++g)
          flux_per_cell[eig][g][c] =
                                     phi_per_assembly[eig][g][cell->user_index()];
      }
    }

    // Materials id
    Vector<double> mat_id_cells(tria.n_active_cells());
    unsigned int c = 0;
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active();
    for (; cell != dof_handler.end(); ++cell, ++c)
    {
      mat_id_cells[c] = static_cast<double>(materials.get_material_id<dim>(
                                              cell)
                                            + 1);
    }

    // Attach data to vectors
    data_out.add_data_vector(mat_id_cells, "Material_id");

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      data_out.add_data_vector(power_per_cell[eig],
        "Power_eig" + num_to_str(eig + 1));
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_groups; ++g)
        data_out.add_data_vector(flux_per_cell[eig][g],
          "phi_assembly_g_" + num_to_str(g) + "_eig"
          + num_to_str(eig));

    std::ofstream output(filename_vtk.c_str());
    data_out.build_patches(n_out_ref);
    data_out.write_vtk(output);
  }

/**
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and time-line.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::clear_vectors ()
  {
    constraints.clear();
    dof_handler.clear();
    tria.clear();
  }
/**
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and time-line.
 */
template <int dim, int n_fe_degree>
  void StaticDiffusion<dim, n_fe_degree>::run (ParameterHandler &prm)
  {
    verbose_cout << "   Input files read" << std::endl;
    timer.start();

    if (t_modes != "lambda")
      cout << "   Computation of static " + t_modes + " modes problem..."
           << std::endl;
    PetscLogDouble memory;

    if (print_grid_flag == true)
    {
      std::string mesh_filename = out_file.substr(0, out_file.size() - 4);
      print_grid(tria, mesh_filename + ".eps");
    }

    verbose_cout << "   making dofs..." << std::flush;
    make_dofs();
    verbose_cout << "  Done! " << std::endl;

    if (type_perturbation == "Rods")
    {
      verbose_cout << "   move_bars...  " << std::flush;
      perturbation.move_bars_static();
      verbose_cout << "  Done bars! " << std::endl;
    }
    if (type_perturbation == "Mechanical_Vibration")
    {
      AssertRelease(geo_type == "Rectangular",
        "This perturbation is only implemented for rectangular geometries");
      perturbation.mechanical_vibration_static();
    }

    cout << "   Grid Done." << " Time = " << timer.cpu_time() << " s."
         << std::endl;
    cout << "      Equations:  DIFFUSION" << std::endl;
    cout << "      Number of Energy Groups:  " << n_groups << std::endl;
    cout << "      Matrix-free type: " << enum_to_string(matrixfree_type)
         << std::endl;
    cout << "      Number of active cells:  " << tria.n_active_cells()
    << " (by partition:";
    for (unsigned int p = 0; p < n_mpi_processes; ++p)
      cout << (p == 0 ? ' ' : '+')
      << (GridTools::count_cells_with_subdomain_association(tria, p));
    cout << ")" << std::endl;
    cout << "      Number of DoFs per block: " << dof_handler.n_dofs()
    << " (by partition:";
    for (unsigned int p = 0; p < n_mpi_processes; ++p)
      cout << (p == 0 ? ' ' : '+')
      << (DoFTools::count_dofs_with_subdomain_association(dof_handler,
        p));
    cout << ")" << std::endl;
    cout << "      Number of Total DoFs: " << dof_handler.n_dofs() * n_groups
    << std::endl;

    if (geo_type == "Composed")
    {
      cout << "      Refs_Radial: " << n_refinements_radial << std::endl;
      if (dim > 2)
        cout << "      Refs_Axial: " << n_refinements_axial << std::endl;
    }
    else
      cout << "      Refinements: " << n_refinements << std::endl;
    cout << "      Degree of FE: " << n_fe_degree << std::endl;
    cout << "      Dofs per cell: " << fe.dofs_per_cell << std::endl;
    cout << std::endl;

    show_cells();

    verbose_cout << "   assembling system " + t_modes + " modes..."
                 << std::flush;
    if (t_modes == "alpha")
      assemble_system_alpha();
    else if (t_modes == "gamma")
      assemble_system_gamma();
    else if (t_modes == "lambda")
      assemble_system_lambda();
    verbose_cout << "  Done! " << std::endl;

    cout << "   Matrices assembled and sent to " + solver_type + " solver"
         << " Time = "
         << timer.cpu_time() << " s." << std::endl;

    verbose_cout << "  solve eps problem... " << std::endl;
    solve_eps();
    verbose_cout << "  Done! " << std::endl;
    if (to_init) // we have finished initialization
      return;

    PetscMemoryGetCurrentUsage(&memory);
    cout << "   Current Memory " << memory * 1e-6 << " MB" << std::endl;

    MPI_Barrier(comm);

    // Get transient options
    bool rom_static = prm.get_bool("ROM_Static");
    if (!rom_static)
    {
      verbose_cout << "Make the reactor critical " << std::endl;
      materials.make_critical(eigenvalues[0]);
    }
    verbose_cout << "postprocess..." << std::flush;
    postprocess();
    verbose_cout << "Done!" << std::endl;

    if (out_flag)
    {
      verbose_cout << "output_results..." << std::flush;
      output_results();
      verbose_cout << "Done!" << std::endl;
    }

    MPI_Barrier(comm);

    // Clear memory
    phi_serial.clear();

  }

template class StaticDiffusion<1, 1> ;
template class StaticDiffusion<1, 2> ;
template class StaticDiffusion<1, 3> ;
template class StaticDiffusion<1, 4> ;
template class StaticDiffusion<1, 5> ;

template class StaticDiffusion<2, 1> ;
template class StaticDiffusion<2, 2> ;
template class StaticDiffusion<2, 3> ;
template class StaticDiffusion<2, 4> ;
template class StaticDiffusion<2, 5> ;

template class StaticDiffusion<3, 1> ;
template class StaticDiffusion<3, 2> ;
template class StaticDiffusion<3, 3> ;
template class StaticDiffusion<3, 4> ;
template class StaticDiffusion<3, 5> ;

