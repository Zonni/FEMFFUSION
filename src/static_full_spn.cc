/**
 *
 * @file   static_full_spn.cc
 * @brief  Implementation of the class StaticFullSPN and the main functions of
 *  the FemFusion program.
 */
#include <deal.II/lac/solver_selector.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/fe/fe_tools.h>

#include <set>

#include "../include/static_full_spn.h"
#include "../include/static_diffusion.h"
#include "../include/femffusion.h"
#include "../include/test.h"
#include "../include/prob_geom.h"
#include "../include/input_geom.h"
#include "../include/matrix_operators/matrix_operators_petsc.h"
#include "../include/printing.h"

using namespace dealii;

/**
 * @brief Constructor of the main class StaticFullSPN.
 * Reads the input file and it reads or builds the grid.
 */
template <int dim, int n_fe_degree>
  StaticFullSPN<dim, n_fe_degree>::StaticFullSPN (
    ParameterHandler &prm,
    std::string input_file,
    const bool verbose,
    const bool silent,
    const bool _to_init)
  :
      comm(PETSC_COMM_WORLD),
      verbose_cout(std::cout, verbose),
      cout(std::cout, !silent),
      to_init(_to_init),
      n_groups(prm.get_integer("Energy_Groups")),
      n_moments(prm.get_integer("N_SPN") + 1),
      n_components(n_moments / 2 * (1 + dim)),
      n_blocks(n_groups * n_components),
      fe(QGaussLobatto<1>(n_fe_degree + 1)),
      dof_handler(tria),
      fe_system(FE_Q<dim>(QGaussLobatto<1>(n_fe_degree + 1)), n_blocks),
      dof_handler_system(tria),
      input_file(input_file),
      materials(verbose_cout)
  {
    verbose_cout << "Start of the program " << std::endl;
    AssertRelease(n_fe_degree > 0, "FE cannot be 0");

    // General Options
    n_eigenvalues = prm.get_integer("N_Eigenvalues");
    n_refinements = prm.get_integer("N_Refinements");

    // Output
    out_file = prm.get("Output_Filename");
    out_flag = prm.get_bool("Output_Flag");
    print_grid_flag = prm.get_bool("Print_Grid_Flag");
    n_out_ref = prm.get_integer("Out_Refinements");

    // Solver
    show_eps_convergence = not silent;
    solver_type = prm.get("Solver_Type");
    renumbering = prm.get("Renumbering");
    tol_eps = prm.get_double("EPS_Tolerance");
    tol_ksp = prm.get_double("KSP_Tolerance");
    residual_norm = true;
    static_ksp_tol = prm.get_bool("Static_KSP_Tolerance");
    matrixfree_type = string_to_enum(prm.get("Matrix_Free_Type"));
    std::string refinement_model = prm.get("Refinement_Model");
    p_init = prm.get_bool("P_Init");

    // Albedo Factors
    std::string str = prm.get("Albedo_Factors");
    trim(str);
    if (str.length() > 0)
      parse_vector(str, albedo_factors, 0);

    init_type = "SP1";

    // Get changes in the parameters through the command line
    get_parameters_from_command_line();

    // ---------------------------------------------------------------------------------
    // GEOMETRY SHAPE
    geo_type = prm.get("Geometry_Type");
    verbose_cout << "geo_type  " << geo_type << std::endl;
    listen_to_material_id = (geo_type == "Composed");

    if (geo_type != "Composed" and geo_type != "Coresim")
      get_mesh_shape(prm);

    // ---------------------------------------------------------------------------------
    // XSEC
    xs_file = prm.get("XSECS_Filename");
    std::string xs_type = prm.get("XSEC_Type");

    precursors_file = prm.get("PREC_Filename");
    delta_xs_file = prm.get("DELTAXSECS_Filename");

    verbose_cout << "xs_file  " << xs_file << std::endl;
    verbose_cout << "xs_type  " << xs_type << std::endl;

    // Parse  Materials Cross sections
    std::vector<unsigned int> geo_ps;
    if (xs_type == "Valkin")
    {
      parse_vector(prm.get("Geometry_Points"), geo_ps,
        2 * assem_per_dim[1], default_geometry_points(assem_per_dim));
    }

    if (geo_type != "Composed")
    {
      // Parse  Materials Cross sections
      materials.reinit(xs_file,
        xs_type,
        n_groups,
        assem_per_dim,
        n_assemblies,
        listen_to_material_id,
        geo_ps,
        precursors_file);
    }

    // ---------------------------------------------------------------------------- //
    // GEOMETRY
    // Mesh File
    if (geo_type == "Unstructured")
    {
      std::string mesh_file = prm.get("Mesh_Filename");
      verbose_cout << "mesh_file " << mesh_file << std::endl;
      get_unstructured_grid(mesh_file);

      verbose_cout << "parsing Boundary_Conditions... " << std::flush;
      parse_vector(prm.get("Boundary_Conditions"), boundary_conditions);
      if (verbose_cout.is_active())
        print_vector(boundary_conditions, false);
      verbose_cout << "Done!" << std::endl;

      // Get the number of assemblies (number of not coarse cells)
      n_assemblies = materials.get_n_assemblies();
    }
    else if (geo_type == "Hexagonal")
    {
      std::string mesh_file = prm.get("Mesh_Filename");
      verbose_cout << "mesh_file " << mesh_file << std::endl;

      verbose_cout << "parsing assembly_pitch... " << std::endl;
      // Pin Pitch definition
      assembly_pitch.resize(dim);
      parse_vector(prm.get("Cell_Pitch_x"), assembly_pitch[0], assem_per_dim[0]);
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

      AssertRelease(xs_type == "XML",
        "Composed geo_type needs XML xs_type ");
      std::string geom_file = prm.get("Geometry_Filename");
      verbose_cout << "geom_file " << geom_file << std::endl;

      verbose_cout << "load geometry file... " << std::flush;
      InputGeom input_geometry;
      input_geometry.load(geom_file);
      verbose_cout << "Done!" << std::endl;

      // Construct
      assem_per_dim.resize(3, 1);
      for (unsigned int d = 0; d < dim; ++d)
      {
        if (d < 2)
          assem_per_dim[d] = input_geometry.planes[0].n_nodes[d];
        else
          assem_per_dim[d] = input_geometry.core.n_planes;
      }

      n_assemblies = assem_per_dim[0] * assem_per_dim[1] * assem_per_dim[2];

      materials.reinit(xs_file,
        xs_type,
        n_groups,
        assem_per_dim,
        n_assemblies,
        listen_to_material_id,
        geo_ps,
        precursors_file);

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
          tria.n_active_cells(),
          tria.n_active_cells());

        unsigned int i = 0;
        typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
        for (; cell != tria.end(); ++cell, ++i)
        {
          cell->set_user_index(user_indices[i]);
        }

        verbose_cout << "Done!" << std::endl;
      }
      else
      {
        cout << " N_Refs_Radial " << n_refinements_radial << std::endl;
        AssertRelease(n_refinements == 0,
          "Not valid n_refinements for geo_type == Composed");

        verbose_cout << "make_composed_geometry... " << std::flush;
        make_composed_geometry(input_geometry, tria,
          n_refinements_radial,
          n_refinements_axial,
          refinement_model);
        verbose_cout << "Done!" << std::endl;
      }
      // Set Geometry Matrix
      verbose_cout << "reading geometry_matrix... make" << std::endl;
      std::vector<unsigned int> geo_ps = default_geometry_points(assem_per_dim);
      materials.set_geometry_matrix(assem_per_dim, geo_ps);
      verbose_cout << " Done!" << std::endl;

      // Copy Boundary Conditions
      boundary_conditions = input_geometry.core.boundary[0];
      for (unsigned int d = 1; d < dim; ++d)
        boundary_conditions.insert(boundary_conditions.end(),
          input_geometry.core.boundary[0].begin(),
          input_geometry.core.boundary[0].end());

      AssertRelease(
        prm.get("Boundary_Conditions") == "",
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
        typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
        for (; cell != tria.end(); ++cell, ++i)
        {
          user_indices.push_back(cell->user_index());
        }

        std::ofstream vec_file_stream(vec_file.c_str(), std::ios::out);
        print_vector_in_file(user_indices, vec_file_stream, "User_indices\n", false, 1);
        verbose_cout << " Done!" << std::endl;
      }
    }
    else if (geo_type == "Rectangular")
    {
      verbose_cout << "parsing assembly_pitch..." << std::endl;
      // Pin Pitch definition
      assembly_pitch.resize(dim);
      parse_vector(prm.get("Cell_Pitch_x"), assembly_pitch[0], assem_per_dim[0]);
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
    // Bars
    n_bars = 0;
    // Default Value for reflector material and bar_top_pos
    bars_top_pos = 0.0;

    std::string bar_file = prm.get("Bar_Filename");
    if (bar_file != "no.bars")
    {
      verbose_cout << "parsing bar_file... " << bar_file << std::flush;
      parse_bar_file(bar_file);
      verbose_cout << "Done!" << std::endl;
      materials_no_bars = materials.get_materials_vector();
    }

    run();
  }

/**
 * @brief Get the information about the Mesh Size and and Shape.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::get_mesh_shape (ParameterHandler &prm)
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
      std::vector<std::vector<unsigned int> > geometry_matrix_2d;
      //parse_matrix(prm.get("Geometry_Matrix"), geometry_matrix_2d,
      //  assem_per_dim[1], assem_per_dim[0]);
      materials.set_geometry_matrix(assem_per_dim, prm.get("Geometry_Matrix"));
      AssertRelease(prm.get("Geometry_Points") == "",
        "Geometry_Matrix and Geometry_Points cannot be defined at the same time");
      verbose_cout << "Done!" << std::endl;

      // Calculation n_assemblies (approximately)
      n_assemblies = 0;
      for (unsigned int r = 0; r < assem_per_dim[1]; ++r)
        for (unsigned int c = 0; c < assem_per_dim[0]; ++c)
          if (geometry_matrix_2d[r][c] != 0)
            n_assemblies++;
      n_assemblies *= assem_per_dim[2];
      verbose_cout << "n_assemblies: " << n_assemblies << std::endl;
    }
    else // if geometry point exists and default complete rectangular problem
    {
      verbose_cout << "parsing Geometry_Points... " << std::flush;
      std::vector<unsigned int> geo_ps;
      parse_vector(prm.get("Geometry_Points"), geo_ps,
        2 * assem_per_dim[1], default_geometry_points(assem_per_dim));

      materials.set_geometry_matrix(assem_per_dim, geo_ps);
      verbose_cout << " Done!" << std::endl;

      // Calculation n_assemblies
      n_assemblies = 0;
      for (unsigned int i = 0; i < geo_ps.size(); i += 2)
        n_assemblies += geo_ps[i + 1] - geo_ps[i] + 1;
      n_assemblies *= assem_per_dim[2];
      verbose_cout << "n_assemblies: " << n_assemblies << std::endl;
    }
  }

/**
 * @brief Construct a rectangular grid.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::make_rectangular_grid ()
  {
    verbose_cout << "get_materials_table..." << std::flush;
    Table<dim, types::material_id> materials_table;
    materials.get_materials_table(materials_table, assem_per_dim);
    verbose_cout << " Done!" << std::endl;

    // Make the grid
    verbose_cout << "Construction the mesh..." << std::flush;
    Point<dim> p1 = ((dim == 1) ? Point<dim>(0.0) :
                                  ((dim == 2) ? Point<dim>(0.0, 0.0) :
                                                Point<dim>(0.0, 0.0, 0.0)));
    if (dim == 1)
      // In 1D it is not implement the colorize=true flag.
      GridGenerator::subdivided_hyper_rectangle(tria, assembly_pitch, p1,
        materials_table,
        false);
    else
      GridGenerator::subdivided_hyper_rectangle(tria, assembly_pitch, p1,
        materials_table,
        true);

    tria.refine_global(n_refinements);

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
 * the boundary id as Physical line/surfaces. Only quadrirateral/hexahedral meshes are accepted.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::get_unstructured_grid (const std::string &mesh_file)
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
        cell->set_user_index(cell->material_id() - 1);
      }
      verbose_cout << "Done!" << std::endl;
    }
    if (dim == 3)
    {
      Assert(fexists(mesh_file), ExcMessage("mesh_file doesn't exist"));
      verbose_cout << "Extruding.." << std::flush;

      Triangulation<2, 2> tria2d;
      GridIn<2> gridin2;
      std::ifstream f(mesh_file.c_str());
      gridin2.attach_triangulation(tria2d);
      gridin2.read_msh(f);
      // Read the z data from somewhere
      // TODO extrude_triangulation to allow z non equally distributed cells
      double height = sum_vector(assembly_pitch[2]);
      unsigned int n_slices = assem_per_dim[2] + 1;

      extrude_triangulation(tria2d, n_slices, height, tria, true);
      tria.refine_global(n_refinements);

      // Set the pertinent Rod (user_index) for each cell.
      Triangulation<3>::active_cell_iterator cell = tria.begin_active();
      for (; cell != tria.end(); ++cell)
      {
        // User_index is the assembly number
        cell->set_user_index(cell->material_id() - 1);
      }
      verbose_cout << "Done!" << std::endl;
    }
  }

/**
 * @brief This Functions parses the Bar file completing
 *  bars_position vector
 *  bar_materials
 *  bar_points
 *  bars_top_pos
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::parse_bar_file (std::string BarFile)
  //
  {
    Assert(fexists(BarFile), ExcMessage("Bar file doesn't exist"));
    std::ifstream input(BarFile.c_str(), std::ios::in);
    std::string sub, keyword;

    verbose_cout << "parsing Move_Banks...  " << std::endl;

    // for every line
    for (std::string line; getline(input, line);)
    {
      std::istringstream iss(line);
      keyword.clear();
      iss >> keyword;

      if (keyword == "#" or keyword == "!" or keyword == "//") // Commentary
        continue;

      else if (keyword == "") //Blank line
        continue;

      // First definition Material and XSecs:
      else if (keyword == "Bank_Configuration")
      {
        bars_position.reserve(n_assemblies / (assem_per_dim[2])); // Reserve space for allocation
        parse_multiline_vector(input, assem_per_dim[1], bars_position, true);

        if (verbose_cout.is_active())
          print_vector(bars_position);
      }
      else if (keyword == "Move_Banks")
      {
        if (verbose_cout.is_active())
          cout << "parsing Move_Banks... ";

        unsigned int bar;
        std::string str;
        double time, height;
        unsigned int bar_mat;
        iss >> n_bars;
        Assert(!iss.fail(),
          ExcMessage("It must be defined the number of bars."));
        Assert(n_bars>0,
          ExcMessage("It must be defined the number of bars."));
        Assert(n_bars<1000,
          ExcMessage("It cannot be defined more than 1000 bars."));

        bar_materials.resize(n_bars);
        bar_points.resize(n_bars);

        for (unsigned int j = 0; j < n_bars; j++)
        {
          get_new_valid_line(input, line);
          std::istringstream iss(line);
          iss >> str;

          // Bar number
          bar = Utilities::string_to_int(str) - 1;
          Assert(bar == j,
            ExcMessage("Bars should be defined in order!"));

          // Bar Materials
          iss >> bar_mat;
          Assert(! iss.fail(),
            ExcMessage("There are not enough (well) bars specified"));
          Assert(bar_mat > 0,
            ExcMessage("There are not enough (well) bars specified"));
          bar_materials[bar] = bar_mat - 1;

          iss >> bar_mat;
          Assert(! iss.fail(),
            ExcMessage("There are not enough (well) bars specified"));
          Assert(bar_mat >= 2,
            ExcMessage("There are not enough (well) bars specified"));

          bar_points[bar].resize(bar_mat);

          for (unsigned int p = 0; p < bar_mat; p++)
          {
            // Bar height and times
            iss >> time;
            Assert(! iss.fail(),
              ExcMessage("There are not enough (well) bars specified"));
            iss >> height;
            Assert(! iss.fail(),
              ExcMessage("There are not enough (well) bars specified"));
            bar_points[bar][p] = std::make_pair(time, height);
          }
        }
      }
      else if (keyword == "BarsTopPosition")
      {
        double num;
        iss >> num;
        Assert(! iss.fail(),
          ExcMessage("There is not specified a valid BarsTopPosition"));
        bars_top_pos = num;
      }
      else
        // Error
        Assert(false, ExcMessage("Invalid Header in Bar_file: " + keyword));
    }

    Assert(bars_top_pos > 0.0,
      ExcMessage("There is not specified a valid BarsTopPosition"));
    verbose_cout << "Done!" << std::endl;

    if (verbose_cout.is_active())
      for (unsigned int bar = 0; bar < bar_points.size(); bar++)
        for (unsigned int p = 0; p < bar_points[bar].size(); p++)
          cout << " bar time " << bar_points[bar][p].first
               << " height  "
               << bar_points[bar][p].second
               << std::endl;
  }

/**
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::get_parameters_from_command_line ()
  {
    // Booleans
    get_bool_from_options("-out_flag", out_flag);
    get_bool_from_options("-print_grid", print_grid_flag);
    get_bool_from_options("-show_eps_convergence", show_eps_convergence);
    get_bool_from_options("-static_ksp_tol", static_ksp_tol);
    get_bool_from_options("-residual_norm", residual_norm);
    get_bool_from_options("-p_init", p_init);

    // Integers
    get_uint_from_options("-n_refinements", n_refinements);
    get_uint_from_options("-n_eigenvalues", n_eigenvalues);
    get_uint_from_options("-n_out_ref", n_out_ref);

    // Reals
    get_double_from_options("-tol_eps", tol_eps);
    get_double_from_options("-tol_ksp", tol_ksp);

    // String
    get_string_from_options("-solver_type", solver_type);
    get_string_from_options("-renumbering", renumbering);
    get_string_from_options("-init_type", init_type);

    // Enums
    get_enum_from_options("-matrixfree_type", matrixfree_type);

    lower_case(renumbering);
    lower_case(solver_type);
    lower_case(init_type);

    //Others
    // -eps_ncv 2
    // -pc_factor_levels 1
    // -pc_factor_mat_ordering_type rcm
  }

/**
 * Create/Allocate the structures that hold the degrees of freedom.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::make_dofs ()
  {
    // Relative to the DOFS
    dof_handler.distribute_dofs(fe);

    // Make the Sparsity Pattern and allocate the Matrices
    n_dofs = dof_handler.n_dofs();
    n_cells = tria.n_active_cells();

    //
    // FeSystem
    //
    // Relative to the DOFS
    dof_handler_system.distribute_dofs(fe_system);
    DoFRenumbering::component_wise(dof_handler_system);

    n_dofs_system = dof_handler_system.n_dofs();

    // Make the Sparsity Pattern and allocate the Matrices
    DynamicSparsityPattern csp1(n_dofs_system);
    DynamicSparsityPattern csp2(n_dofs_system);
    SparsityPattern sp1, sp2; // @suppress("Abstract class cannot be instantiated")

    DoFTools::make_sparsity_pattern(dof_handler_system, csp1, constraints_system, true);
    DoFTools::make_sparsity_pattern(dof_handler_system, csp2, constraints_system, true);

    // constraints_system.condense(csp1);
    // constraints_system.condense(csp2);

    sp1.copy_from(csp1);
    sp2.copy_from(csp2);

    A.reinit(sp1);
    B.reinit(sp2);

  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::assemble_system ()
  {

    unsigned int moment_i, moment_j;
    unsigned int group_i, group_j;
    unsigned int mi, mj;
    unsigned int vacuum = 0.0;

    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler_system.get_fe(),
      quadrature_formula,
      update_values | update_gradients |
      update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(
      dof_handler_system.get_fe(),
      face_quadrature_formula,
      update_values | update_quadrature_points |
      update_JxW_values
      | update_normal_vectors);

    const unsigned int dofs_per_cell = fe_system.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_A(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_B(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_bound(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Indicates which components of a vector-valued finite element constitute a single scalar
    std::vector<std::vector<FEValuesExtractors::Scalar>> scalar_flux(n_moments / 2,
      std::vector<FEValuesExtractors::Scalar>(n_groups));
    std::vector<std::vector<FEValuesExtractors::Vector>> current(n_moments / 2,
      std::vector<FEValuesExtractors::Vector>(n_groups));

    unsigned int n_vec_moments = 2;

    for (unsigned int mm = 0; mm < n_moments / 2; mm++)
      for (unsigned int g = 0; g < n_groups; g++)
      {
        scalar_flux[mm][g] = FEValuesExtractors::Scalar(mm * (dim + 1) * n_groups + g);
//        std::cout << "scalar" << mm * (dim + 1) * n_groups + g << std::endl;
        current[mm][g] = FEValuesExtractors::Vector(
          (mm + 1) * n_groups + mm * dim * n_groups + g * dim);
//        std::cout << "vector" << (mm + 1) * n_groups + mm * dim * n_groups + g * dim
//                  << std::endl;
      }

    typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler_system.begin_active(),
        endc = dof_handler_system.end();
    for (; cell != endc; ++cell)
    {

      fe_values.reinit(cell);
      const unsigned int mat = materials.get_material_id<dim>(cell);

      cell_A = 0;
      cell_B = 0;

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {

          const unsigned int block_i = fe_system.system_to_component_index(i).first;
          const unsigned int block_j = fe_system.system_to_component_index(j).first;

          block_to_moment_group(block_i, moment_i, group_i);
          block_to_moment_group(block_j, moment_j, group_j);

          if (moment_i == 0 and moment_j == 0)
          {
            if (group_i == group_j)
            {
              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_A(i, j) += materials.get_sigma_r(group_i, mat)
                                * fe_values[scalar_flux[0][group_i]].value(i, q)
                                * fe_values[scalar_flux[0][group_j]].value(j, q)
                                * fe_values.JxW(q);
            }
            else
            {
              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_A(i, j) += -materials.get_sigma_s(group_j, group_i, mat)
                                * fe_values[scalar_flux[0][group_i]].value(i, q)
                                * fe_values[scalar_flux[0][group_j]].value(j, q)
                                * fe_values.JxW(q);
            }

            for (unsigned int q = 0; q < n_q_points; ++q)
              cell_B(i, j) += materials.get_xi_nu_sigma_f(group_j, group_i, mat)
                              * fe_values[scalar_flux[0][group_i]].value(i, q)
                              * fe_values[scalar_flux[0][group_j]].value(j, q)
                              * fe_values.JxW(q);
          }
          else if (moment_j - moment_i == 1 or moment_i - moment_j == 1)
          {
            if (group_i == group_j)
            {
              if (moment_i % 2 == 0)
              {
                for (unsigned int q = 0; q < n_q_points; ++q)
                  cell_A(i, j) +=
                      -spn_coeff[moment_i][moment_j] // Cambio de signo debido a la discretizacion FEM
                      * fe_values[scalar_flux[moment_i / n_vec_moments][group_i]].gradient(
                        i, q)
                      * fe_values[current[moment_j / n_vec_moments][group_j]].value(j, q)
                      * fe_values.JxW(q);
              }
              else
              {
                for (unsigned int q = 0; q < n_q_points; ++q)
                  cell_A(i, j) +=
                      spn_coeff[moment_i][moment_j]
                      * fe_values[current[moment_i / n_vec_moments][group_i]].value(i, q)
                      * fe_values[scalar_flux[moment_j / n_vec_moments][group_j]].gradient(
                        j, q)
                      * fe_values.JxW(q);
              }
            }
          }

          else if (moment_i == moment_j and moment_i > 0)
          {
            if (moment_i % 2 == 0)
            {
              if (group_i == group_j)
              {
                for (unsigned int q = 0; q < n_q_points; ++q)
                  cell_A(i, j) += materials.get_sigma_t(group_i, mat)
                      * fe_values[scalar_flux[moment_i / n_vec_moments][group_i]].value(i,
                        q)
                      * fe_values[scalar_flux[moment_j / n_vec_moments][group_j]].value(j,
                        q)
                      * fe_values.JxW(q);
              }
            }
            else
            {

              if (group_i == group_j)
              {
                for (unsigned int q = 0; q < n_q_points; ++q)
                  cell_A(i, j) += materials.get_sigma_t(group_i, mat)
                      * fe_values[current[moment_i / n_vec_moments][group_i]].value(i, q)
                      * fe_values[current[moment_j / n_vec_moments][group_j]].value(j, q)
                      * fe_values.JxW(q);
              }
            }
          }
        }

      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
        {
          types::boundary_id boundary_id = cell->face(f)->boundary_id();
          AssertIndexRange(boundary_id, boundary_conditions.size());

          if (boundary_conditions[boundary_id] > 1)
          {
            fe_face_values.reinit(cell, f);
            switch (boundary_conditions[boundary_id])
            {
              case 2: // Vacuum BC
                vacuum = 1.0;
                break;
              default: // Custom Albedo BC
                vacuum = 100.0;
                break;
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const unsigned int block_i = fe_system.system_to_component_index(i).first;
                const unsigned int block_j = fe_system.system_to_component_index(j).first;

                block_to_moment_group(block_i, moment_i, group_i);
                block_to_moment_group(block_j, moment_j, group_j);

                if (group_i == group_j and moment_j % 2 == 0 and moment_i % 2 == 0)
                {

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                  {

                    mi = moment_i / n_vec_moments;
                    mj = moment_j / n_vec_moments;
                    cell_A(i, j) += vacuum *
                                    bound_vacuum_coeff[mi][mj]
                                    * fe_face_values[scalar_flux[mi][group_i]].value(i, q)
                                    * fe_face_values[scalar_flux[mj][group_j]].value(j, q)
                                    * fe_face_values.JxW(q);
                  }
                }
              }
          }

        }

      cell->get_dof_indices(local_dof_indices);

      constraints_system.distribute_local_to_global(
        cell_A, local_dof_indices, A);
      constraints_system.distribute_local_to_global(
        cell_B, local_dof_indices, B);

    }

    //
    A.compress(VectorOperation::add);
    B.compress(VectorOperation::add);

    print_matrices();

  }

template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::block_to_moment_group (unsigned int block,
    unsigned int &moment,
    unsigned int &group)
  {

    if (block < n_groups)
    {
      moment = 0;
      group = block;
    }
    else if (block < n_groups * (1 + dim))
    {
      moment = 1;
      group = (block - n_groups) / dim;
    }
    else if (block < (dim + 1) * n_groups + n_groups)
    {
      moment = 2;
      group = (block - n_groups * (1 + dim));
    }
    else
    {
      moment = 3;
      group = (block - (dim + 2) * n_groups) / dim;
    }

    return;
  }

/**
 * @brief Print matrices in the given file in a Matlab way if a file is specified
 *  with '-print_matrices_matlab' command line option.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode
  StaticFullSPN<dim, n_fe_degree>::print_matrices ()
  {
    PetscErrorCode ierr;
    std::string print_matrices_matlab;

    ierr = get_string_from_options("-print_matrices",
      print_matrices_matlab);

    if (!print_matrices_matlab.empty())
    {

      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);

      print_matrix_in_matlab(A, "A", out, 13);
      print_matrix_in_matlab(B, "B", out, 13);

      out.close();
    }

    return ierr;
  }

/**
 * @brief Solve the eigenvalue problem
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::solve_eps ()
  {

    EPS eps;
    ST st;
    KSP ksp;
    PC pc;
    EPSCreate(comm, &eps); // @suppress("Invalid arguments")
    EPSSetDimensions(eps, n_eigenvalues, PETSC_DEFAULT, PETSC_DEFAULT);

    EPSSetOperators(eps, B, A);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetTolerances(eps, tol_eps, 1e5);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    EPSSetType(eps, EPSGD);
    EPSGetST(eps, &st);
    STGetKSP(st, &ksp);
    KSPGetPC(ksp, &pc);
    KSPSetType(ksp, KSPPREONLY);
    PCSetType(pc, PCILU);
    PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
    PCFactorSetMatOrderingType(pc, MATORDERINGRCM);

    //PCFactorSetLevels(pc, 10);
    //KSPSetUp(ksp);
    EPSSetFromOptions(eps);
    KSPSetFromOptions(ksp);
    EPSSetUp(eps);
    EPSSolve(eps);

    EPSConvergedReason reason;
    EPSGetConvergedReason(eps, &reason);
    AssertRelease(reason > 0,
      "EPS not converged, reason: " + num_to_str(static_cast<int>(reason)));

    phi_sol.resize(n_eigenvalues);
    phi.resize(n_eigenvalues);

    eigenvalues.resize(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      phi_sol[eig].reinit(comm, n_dofs_system, n_dofs_system); // @suppress("Invalid arguments")
      phi[eig].resize(n_components);
      for (unsigned int nc = 0; nc < n_components; nc++)
        phi[eig][nc].reinit(n_groups, comm, n_dofs, n_dofs); // @suppress("Invalid arguments")
    }

    std::vector<std::vector<PETScWrappers::MPI::Vector>> phi_sep(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      phi_sep[eig].resize(n_groups * n_components);
      for (unsigned int nb = 0; nb < n_groups * n_components; nb++)
        phi_sep[eig][nb].reinit(comm, n_dofs, n_dofs); // @suppress("Invalid arguments")
    }

    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      EPSGetEigenpair(eps, eig, &eigenvalues[eig], NULL,
        phi_sol[eig], PETSC_NULL);

      separate_vectors(dof_handler, dof_handler_system, phi_sol[eig], phi_sep[eig]);
    }

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      for (unsigned int nc = 0; nc < n_components; nc++)
        for (unsigned int ng = 0; ng < n_groups; ng++)
          phi[eig][nc].block(ng) = phi_sep[eig][nc * n_groups + ng];

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int nb = 0; nb < n_groups * n_components; nb++)
        phi_sep[eig][nb].clear();

    EPSDestroy(&eps);

    A.clear();
    B.clear();
  }

/**
 * @brief Solve the eigenvalue problem
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::solve_eps_B ()
  {
    EPS eps;
    ST st;
    KSP ksp;
    PC pc;
    EPSCreate(comm, &eps); // @suppress("Invalid arguments")
    EPSSetDimensions(eps, n_eigenvalues, PETSC_DEFAULT, PETSC_DEFAULT);

    EPSSetOperators(eps, B, A);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetTolerances(eps, tol_eps, 1e5);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    EPSSetType(eps, EPSKRYLOVSCHUR);
    EPSGetST(eps, &st);
    STGetKSP(st, &ksp);
    KSPGetPC(ksp, &pc);
    //KSPSetType(ksp, KSPPREONLY);
//    PCSetType(pc, PCILU);
//    PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
//    PCFactorSetMatOrderingType(pc, MATORDERINGRCM);

//    KSPSetFromOptions(ksp);
    EPSSetFromOptions(eps);
    EPSSetUp(eps);
    EPSSolve(eps);

    EPSConvergedReason reason;
    EPSGetConvergedReason(eps, &reason);
    AssertRelease(reason > 0,
      "EPS not converged, reason: " + num_to_str(static_cast<int>(reason)));

    phi_sol.resize(n_eigenvalues);
    phi.resize(n_eigenvalues);

    eigenvalues.resize(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      phi_sol[eig].reinit(comm, n_dofs_system, n_dofs_system); // @suppress("Invalid arguments")
      phi[eig].resize(n_components);
      for (unsigned int nc = 0; nc < n_components; nc++)
        phi[eig][nc].reinit(n_groups, comm, n_dofs, n_dofs); // @suppress("Invalid arguments")
    }

    std::vector<std::vector<PETScWrappers::MPI::Vector>> phi_sep(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      phi_sep[eig].resize(n_groups * n_components);
      for (unsigned int nb = 0; nb < n_groups * n_components; nb++)
        phi_sep[eig][nb].reinit(comm, n_dofs, n_dofs); // @suppress("Invalid arguments")
    }

    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      EPSGetEigenpair(eps, eig, &eigenvalues[eig], NULL,
        phi_sol[eig], PETSC_NULL);

      separate_vectors(dof_handler, dof_handler_system, phi_sol[eig], phi_sep[eig]);
    }

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      for (unsigned int nc = 0; nc < n_components; nc++)
        for (unsigned int ng = 0; ng < n_groups; ng++)
          phi[eig][nc].block(ng) = phi_sep[eig][nc * n_groups + ng];

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int nb = 0; nb < n_groups * n_components; nb++)
        phi_sep[eig][nb].clear();

//    PCDestroy(&pc);
//    KSPDestroy(&ksp);
    EPSDestroy(&eps);

    A.clear();
    B.clear();
  }

/**
 * @brief Normalize the problem to to mean neutron density power equal 1.
 * Also, calculate and print the mean values per assembly.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::postprocess ()
  {
    // Make reactor critical
    materials.make_critical(eigenvalues[0]);

    // Erase the content of output file
    std::ofstream
    out(out_file.c_str(), std::ios::out);

    // Print the eigenvalues in the outFile
    print_logo(out);
    print_vector_in_file(eigenvalues, out, "The Eigenvalues are: \n", false, 6);
    out << "Problem File: " << input_file << "\n";
    print_in_file(timer.cpu_time(), out, "CPU Time: ", 4);
    if (geo_type == "Composed")
    {
      print_in_file(n_refinements_radial, out, "Radial Refinements: "); // @suppress("Invalid arguments")
      print_in_file(n_refinements_axial, out, "Axial Refinements: "); // @suppress("Invalid arguments")
    }
    else
      print_in_file(n_refinements, out, "Global Refinements: "); // @suppress("Invalid arguments")

    print_in_file(n_fe_degree, out, "Degree of FE: ", 1);
    print_in_file(n_groups, out, "Number of Groups: ", 1);
    print_in_file(n_moments / 2 + 1, out, "Full SPN: ", 1);
    print_in_file(dof_handler.n_dofs(), out, "DoFs per Group: ", 1);
    print_in_file(tria.n_active_cells(), out, "Number of active cells: ", 1);
    print_vector_in_file(assem_per_dim, out, "Mesh Size: ", true, 1);
    print_in_file("", out, "\n", 1);
    out.close();

    // Initialize all that  is needed to iterate over dofs and cells
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    FEValues<dim> fe_values(
      fe, quadrature_formula,
      update_values | update_quadrature_points | update_volume_elements
      | update_JxW_values);

    unsigned int n_q_points = quadrature_formula.size();
    unsigned int n_cells_out = n_assemblies;
    double power_cell = 0;
    double phi_cell;
    double sigma_f;

    // Initialize and resize the vectors where it is stored the solution
    power_per_assembly.resize(n_eigenvalues, std::vector<double>(n_cells_out));
    phi_per_assembly.resize(
      n_eigenvalues,
      std::vector<std::vector<double> >(n_groups,
        std::vector<double>(n_cells_out, 0.0)));
    std::vector<double> volume_per_assembly(n_cells_out);
    std::vector<double> local_phi(n_q_points);

    // For all defined eigenvalues
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      double volume = 0.0;
      double norm = 0.0;

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

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          sigma_f = materials.get_sigma_f(g, materials.get_material_id<dim>(cell));

          fe_values.get_function_values(phi[eig][0].block(g), local_phi);

          phi_cell = 0.0;
          for (unsigned int q = 0; q < n_q_points; q++)
            phi_cell += local_phi[q] * fe_values.JxW(q);

          power_cell += sigma_f * phi_cell;

          phi_per_assembly[eig][g][cell->user_index()] = phi_cell;
        }

        volume += cell->measure();
        norm += power_cell;
        power_per_assembly[eig][cell->user_index()] += power_cell;
        volume_per_assembly[cell->user_index()] += cell->measure();
      }

      norm /= volume;
      //cout << "  Volume "   << volume << std::endl;

      // Normalize the values of the power and flows per cell
      normalize_vector(power_per_assembly[eig], norm);
      for (unsigned int g = 0; g < n_groups; ++g)
        normalize_vector(phi_per_assembly[eig][g], norm);

      // Normalize the values of the power and flows per Rod
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
          for (unsigned int c = 0; c < n_components; ++c)
            for (unsigned int g = 0; g < phi[eig][c].n_blocks(); ++g)
            {
              normalize_vector(phi_per_assembly[eig][g], -1);
              normalize_vector(phi[eig][c].block(g), -1);
            }
          change_flag = true;
        }
        else if (power_per_assembly[eig][i] > +0.0001)
        {
          change_flag = true;
        }
        i++;
      }

      // Normalize also the Fluxes passed to the .vtk
      for (unsigned int c = 0; c < n_components; ++c)
        for (unsigned int g = 0; g < n_groups; ++g)
          normalize_vector(phi[eig][c].block(g), norm);

      // Calculate the axial power distribution
      std::vector<double> power_axial;
      std::vector<std::vector<double> > phi_axial;
      std::vector<double> volume_per_plane;
      unsigned int n_assemblies_per_plane = n_assemblies / assem_per_dim[2];
      Assert(n_assemblies_per_plane > 0,
        ExcMessage("n_assemblies cannot be 0"));
      if (dim == 3)
      {
        power_axial.resize(assem_per_dim[2], 0.0);
        phi_axial.resize(n_groups, std::vector<double>(assem_per_dim[2]));
        volume_per_plane.resize(assem_per_dim[2], 0.0);

        for (unsigned int i = 0; i < power_per_assembly[eig].size(); i++)
          if (power_per_assembly[eig][i] > 1e-5)
          {
            power_axial[i / n_assemblies_per_plane] +=
                power_per_assembly[eig][i] * volume_per_assembly[i];

            volume_per_plane[i / n_assemblies_per_plane] +=
                                                            volume_per_assembly[i];

            for (unsigned int g = 0; g < n_groups; ++g)
              phi_axial[g][i / n_assemblies_per_plane] +=
                  phi_per_assembly[eig][g][i] * volume_per_assembly[i];
          }

        // Normalize axial power
        normalize_vector(power_axial, volume_per_plane);
        for (unsigned int g = 0; g < n_groups; ++g)
          normalize_vector(phi_axial[g], volume_per_plane);
      }

      // Print the Powers Distribution and the Fluxes
      const unsigned int presicion = 6;
      std::ofstream out_stream(out_file.c_str(), std::ios::app);
      out_stream << "MODE " + Utilities::int_to_string(eig + 1) << "\n";
      print_cell_distribution_in_file(dim, power_per_assembly[eig],
        assem_per_dim, out_stream, materials, "Neutron Power\n", false, presicion);

      //      if (dim == 3)
      //        print_vector_in_file(power_axial,
      //          out_stream,
      //          "Axial power distribution:\n");

      for (unsigned int g = 0; g < n_groups; ++g)
      {
        print_cell_distribution_in_file(dim, phi_per_assembly[eig][g],
          assem_per_dim, out_stream, materials,
          "Group " + num_to_str(g + 1) + " flux\n", false, presicion);
        //        if (dim == 3)
        //          print_vector_in_file(phi_axial[g], out_stream,
        //            "Group " + num_to_str(g + 1) + "  axial flux distribution:\n");
      }
      // Add Some blank lines
      out_stream << "\n\n";
      out_stream.close();
    }

    std::vector<std::vector<PETScWrappers::MPI::Vector> > phi_sep(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      phi_sep[eig].resize(n_groups * n_components);
      for (unsigned int nb = 0; nb < n_groups * n_components; nb++)
        phi_sep[eig][nb].reinit(comm, n_dofs, n_dofs); // @suppress("Invalid arguments")
    }

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      for (unsigned int nc = 0; nc < n_components; nc++)
        phi[eig][nc].compress(VectorOperation::insert);

    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      for (unsigned int nc = 0; nc < n_components; nc++)
        for (unsigned int ng = 0; ng < n_groups; ng++)
          phi_sep[eig][nc * n_groups + ng] = phi[eig][nc].block(ng);

    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      join_vectors(dof_handler, dof_handler_system, phi_sep[eig], phi_sol[eig]);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int nb = 0; nb < n_groups * n_components; nb++)
        phi_sep[eig][nb].clear();

  }

/**
 * @brief Function that creates the output files .vtk.
 * It should be extended to create the power.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::output_results () const
  {
    DataOut<dim, DoFHandler<dim> > data_out;
    data_out.attach_dof_handler(dof_handler);
    std::string filename_vtk = out_file + ".vtk";

    std::vector<Vector<double> > power_per_cell(n_eigenvalues,
      Vector<double>(tria.n_active_cells()));
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      power_per_cell[eig].reinit(tria.n_active_cells());

    std::vector<std::vector<Vector<double> > > flux_per_cell(
      n_eigenvalues,
      std::vector<Vector<double> >(n_groups,
        Vector<double>(tria.n_active_cells())));
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_groups; ++g)
        flux_per_cell[eig][g].reinit(tria.n_active_cells());

    // Power and Fulx per assembly
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      for (unsigned int g = 0; g < n_groups; ++g)
        data_out.add_data_vector(
          phi[0][0].block(g),
          "phi_g" + num_to_str(g + 1) + "_eig_" + num_to_str(eig + 1));

      unsigned int c = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
                                                                    tria.begin_active();
          cell != tria.end(); ++cell, ++c)
      {
        power_per_cell[eig][c] = power_per_assembly[eig][cell->user_index()];
        for (unsigned int g = 0; g < n_groups; ++g)
          flux_per_cell[eig][g][c] = phi_per_assembly[eig][g][cell->user_index()];
      }
    }

    // Materials id
    Vector<double> mat_id_cells(tria.n_active_cells());
    unsigned int c = 0;
    for (typename DoFHandler<dim>::active_cell_iterator cell =
                                                               dof_handler.begin_active();
        cell != dof_handler.end(); ++cell, ++c)
    {
      mat_id_cells[c] = static_cast<double>(materials.get_material_id<dim>(cell));
    }

    // Attach data to vectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      data_out.add_data_vector(power_per_cell[eig],
        "Power_eig" + num_to_str(eig + 1));
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_groups; ++g)
        data_out.add_data_vector(
          flux_per_cell[eig][g],
          "phi_assembly_g_" + num_to_str(g) + "_eig" + num_to_str(eig));

    data_out.add_data_vector(mat_id_cells, "Material_id");

    std::ofstream output(filename_vtk.c_str());
    data_out.build_patches(n_out_ref);
    data_out.write_vtk(output);
  }

/**
 * @brief Move control bars.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::move_bars ()
  {
    unsigned int bar;
    double bar_pos = 0;
    for (unsigned int plant_pos = 0; plant_pos < bars_position.size(); ++plant_pos)
    {
      if (bars_position[plant_pos] > 0)
      {
        bar = bars_position[plant_pos] - 1;
        for (unsigned int p = 1; p < bar_points[bar].size(); p++)
        {
          bar_pos = bar_points[bar][p - 1].second;
          break;
        }
        Assert(bar_pos < 9e5,
          ExcMessage("Error in time of the bars definition"));

        move_bar_volume_homogenized(plant_pos, bar_pos, bar_materials[bar],
          bar);
      }
    }
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::move_bar_volume_homogenized (
    unsigned int plant_pos,
    double bar_pos,
    unsigned int mat_bar,
    unsigned int bar)
  {
    const unsigned int move_dim = dim - 1;
    unsigned int n_assemblies_per_plane = n_assemblies
                                          / assem_per_dim[move_dim];
    Assert(n_assemblies_per_plane > 0, ExcMessage("n_assemblies cannot be 0"));

    double maxp = 0.0;
    double minp = 0.0;
    unsigned int mat_no_bar;
    double frac;
    std::vector<bool> is_done(n_assemblies, false);
    unsigned int averaged_mat = n_mats + bar;

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc =
               dof_handler.end();
    for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    {
      // We only update each assembly material once
      // However we must iterate over all cells in order to get the z position
      // of the cells

      if (is_done[cell->user_index()] == false)
      {
        is_done[cell->user_index()] = true;
        if (plant_pos == cell->user_index() % n_assemblies_per_plane)
        {

          getMaxMinVertex(cell, move_dim, maxp, minp);
          // The bar_pos is in the middle of the cell.
          // We create a new material at the end of the defined materials
          // that has volume-averaged cross sections
          if ((bar_pos - minp > 1e-8) and (maxp - bar_pos > 1e-8))
          {
            Assert(
              (bars_top_pos - minp > 1e-8 and maxp - bars_top_pos > 1e-8) == false,
              ExcNotImplemented());

            // Calculate the fraction of the cell occupied by the bar
            frac = (maxp - bar_pos) / (maxp - minp);

            mat_no_bar = materials_no_bars[cell->user_index()];
            materials.create_new_mixed_mat(averaged_mat,
              frac,
              mat_bar,
              mat_no_bar,
              cell->user_index());

          }

          // The bar occupy all the cell
          else if (bar_pos - minp < 1e-8 && maxp - bars_top_pos < 1e-8)
          {
            materials.set_materials_id(cell->user_index(), mat_bar);
          }
          // The bar does no occupy any space in the cell
          else
          {
            mat_no_bar = materials_no_bars[cell->user_index()];
            materials.set_materials_id(cell->user_index(), mat_no_bar);
          }
        }
      }
    }
  }

/**
 *
 */
//template <int dim, int n_fe_degree>
//  void
//  StaticFullSPN<dim, n_fe_degree>::output_to_noise ()
//  {
//    std::string filename_sta = out_file + ".sta";
//    // Erase the content of output file
//    std::ofstream out(filename_sta.c_str(), std::ios::out);
//
//    const unsigned int precision = 12;
//    print_in_file(eigenvalues[0], out, "keff\n", precision);
//
//    print_vector_in_file(phi_sol[0], out, "phi " + num_to_str(phi_sol[0].size()) + "\n",
//      true,
//      precision);
//
//    std::cout << "n_dofs * n_components: " << n_dofs * n_components << std::endl;
//    std::cout << "phi_sol[0].size(): " << phi_sol[0].size() << std::endl;
//
//    out.close();
//  }
/**
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and the time-line.
 */
template <int dim, int n_fe_degree>
  void
  StaticFullSPN<dim, n_fe_degree>::run ()
  {
    verbose_cout << "   Input files read" << std::endl;
    timer.start();

    if (print_grid_flag == true)
    {
      std::string mesh_filename = out_file.substr(0,
        out_file.size() - 4);
      print_grid(tria, mesh_filename + ".eps");
    }

    // Get Maximum memory
    PetscLogDouble memory;
    PetscMemorySetGetMaximumUsage();

    verbose_cout << "   making dofs..." << std::flush;
    make_dofs();
    verbose_cout << "  Done! " << std::endl;

    verbose_cout << "   move_bars...  " << std::flush;
    move_bars();
    verbose_cout << "  Done! " << std::endl;

    cout << "   Grid Done." << " Time = " << timer.cpu_time() << " s." << std::endl;
    cout << "      Equations:  P" << n_moments - 1 << std::endl;
    cout << "      Number of active cells:  " << tria.n_active_cells() << std::endl;
    cout << "      Number of DoFs per block: " << dof_handler.n_dofs() << std::endl;
    cout << "      Number of Total DoFs: " << dof_handler_system.n_dofs()
    << std::endl;
    if (geo_type == "Composed")
    {
      cout << "      Refs_Radial: " << n_refinements_radial << std::endl;
      cout << "      Refs_Axial: " << n_refinements_axial << std::endl;
    }
    else
      cout << "      Refinements: " << n_refinements << std::endl;
    cout << "      Degree of FE: " << n_fe_degree << std::endl;
    cout << "      Dofs per cell: " << fe.dofs_per_cell << std::endl;
    cout << std::endl;

    verbose_cout << "   assembling system..." << std::flush;
    assemble_system();
    verbose_cout << "  Done! " << std::endl;

    bool solver_B = false;
    get_bool_from_options("-eps_solver_B", solver_B);
    if (solver_B)
    {
      cout << "   Matrices assembled and sent to solve_eps_B solver"
           << " Time = "
           << timer.cpu_time() << " s." << std::endl;

      verbose_cout << "  solve_eps_B problem... " << std::endl;
      solve_eps_B();
      verbose_cout << "  Done! " << std::endl;
    }
    else
    {
      cout << "   Matrices assembled and sent to solve_eps solver"
           << " Time = "
           << timer.cpu_time() << " s." << std::endl;

      verbose_cout << "  solve eps problem... " << std::endl;
      solve_eps();
      verbose_cout << "  Done! " << std::endl;
    }

    if (to_init) // we have finishes initialization
      return;

    // Some Prints
    PetscMemoryGetMaximumUsage(&memory);
    cout << "   Max Memory " << memory * 1e-6 << " MB." << std::endl;
    cout << "   Memory consumption of matrix elements "
         << memory_consumption * 1e-6
         << " MB." << std::endl;
    cout << std::setprecision(6);
    // Print the Eigenvalues calculated
    cout << "   The eigenvalues are:" << "     Time = " << timer.cpu_time() << " s."
         << std::endl;
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
    {
      cout << "      K" << i << " : " << eigenvalues[i] << std::endl;
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

    cout << "                           Time: " << timer.cpu_time() << " s." << std::endl;

  }

template class StaticFullSPN<1, 1> ;
template class StaticFullSPN<1, 2> ;
template class StaticFullSPN<1, 3> ;
template class StaticFullSPN<1, 4> ;
template class StaticFullSPN<1, 5> ;

template class StaticFullSPN<2, 1> ;
template class StaticFullSPN<2, 2> ;
template class StaticFullSPN<2, 3> ;
template class StaticFullSPN<2, 4> ;
template class StaticFullSPN<2, 5> ;

template class StaticFullSPN<3, 1> ;
template class StaticFullSPN<3, 2> ;
template class StaticFullSPN<3, 3> ;
template class StaticFullSPN<3, 4> ;
template class StaticFullSPN<3, 5> ;

