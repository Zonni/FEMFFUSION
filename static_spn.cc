/**
 *
 * @file   static_spn.cc
 * @brief  Implementation of the class StaticSPN and the main functions of
 *  the FemFusion program.
 */

#include <deal.II/lac/solver_selector.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/fe/fe_tools.h>

#include "static_spn.h"
#include "static_diffusion.h"
#include "femffusion.h"
#include "eps_solver.h"
#include "test.h"
#include "prob_geom.h"
#include "input_geom.h"
#include "matrix_operators/matrix_operators_petsc.h"
#include "printing.h"

using namespace dealii;

/**
 * @brief Constructor of the main class StaticSPN.
 * Reads the input file and it reads or builds the grid.
 */
template <int dim, int n_fe_degree>
  StaticSPN<dim, n_fe_degree>::StaticSPN (ParameterHandler &prm,
    std::string input_file,
    const bool verbose,
    const bool silent,
    const bool _to_init) :
      comm(PETSC_COMM_WORLD),
      n_mpi_processes(Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
      n_local_cells(numbers::invalid_unsigned_int),
      verbose_cout(std::cout, verbose and this_mpi_process == 0),
      cout(std::cout, !silent and this_mpi_process == 0),
      to_init(_to_init),
      n_groups(prm.get_integer("Energy_Groups")),
      fe(QGaussLobatto<1>(n_fe_degree + 1)),
      tria(comm),
      dof_handler(tria),
      L(comm, dof_handler, constraints),
      M(comm, dof_handler, constraints),
      input_file(input_file),
      materials(verbose_cout),
      perturbation(prm, materials, dof_handler, verbose_cout)
  {
    verbose_cout << "Start of the program " << std::endl;
    AssertRelease(n_fe_degree > 0, "FE cannot be 0");

    // General Options
    n_eigenvalues = prm.get_integer("N_Eigenvalues");
    n_refinements = prm.get_integer("N_Refinements");
    n_moments = prm.get_integer("N_SPN") / 2 + 1;

    materials.transient = prm.get_bool("Transient");

    // Output
    out_file = prm.get("Output_Filename");
    out_flag = prm.get_bool("Output_Flag");
    print_grid_flag = prm.get_bool("Print_Grid_Flag");
    n_out_ref = prm.get_integer("Out_Refinements");
    out_to_noise = prm.get_bool("Output_To_Noise");

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
      get_unstructured_grid(mesh_file);

      verbose_cout << "parsing Boundary_Conditions... " << std::flush;
      parse_vector(prm.get("Boundary_Conditions"), boundary_conditions);
      if (verbose_cout.is_active())
        print_vector(boundary_conditions, false);
      verbose_cout << "Done!" << std::endl;
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
        cout << " N_Refs_Radial " << n_refinements_radial << std::endl;
        AssertRelease(n_refinements == 0,
          "Not valid n_refinements for geo_type == Composed");

        verbose_cout << "make_composed_geometry... " << std::flush;
        make_composed_geometry(input_geometry, tria, n_refinements_radial,
          n_refinements_axial, refinement_model);
        verbose_cout << "Done!" << std::endl;
      }
      // Set Geometry Matrix
      verbose_cout << "reading geometry_matrix... make" << std::endl;
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
        print_vector_in_file(user_indices, vec_file_stream, "User_indices\n", false, 1);
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
    // Bars
    type_perturbation = prm.get("Type_Perturbation");

    if (type_perturbation == "Rods")
    {
      std::string bar_file = prm.get("Bar_Filename");
      AssertRelease(bar_file != "no.bar", "It is necessary a bar file");
      verbose_cout << "parsing bar_file... " << bar_file << std::flush;
      perturbation.parse_bar_file(bar_file);
      verbose_cout << "Done!" << std::endl;
    }

    if (to_init == false)
    {
      run();
    }
  }

/**
 * @brief Get the information about the Mesh Size and and Shape.
 */
template <int dim, int n_fe_degree>
  void StaticSPN<dim, n_fe_degree>::get_mesh_shape (ParameterHandler &prm)
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
  void StaticSPN<dim, n_fe_degree>::make_rectangular_grid ()
  {
    verbose_cout << "get_materials_table..." << std::flush;
    Table<dim, types::material_id> materials_table;
    materials.get_materials_table(materials_table, assem_per_dim);
    verbose_cout << " Done!" << std::endl;

    // Make the grid
    verbose_cout << "Construction the mesh..." << std::flush;
    Point<dim> p1 = (
                    (dim == 1) ?
                                 Point<dim>(0.0) :
                                 ((dim == 2) ?
                                               Point<dim>(0.0, 0.0) :
                                               Point<dim>(0.0, 0.0, 0.0)));
    if (dim == 1)
      // In 1D it is not implement the colorize=true flag.
      GridGenerator::subdivided_hyper_rectangle(tria, assembly_pitch, p1,
        materials_table, false);
    else
      GridGenerator::subdivided_hyper_rectangle(tria, assembly_pitch, p1,
        materials_table, true);

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
  void StaticSPN<dim, n_fe_degree>::get_unstructured_grid (
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
      typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
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
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void StaticSPN<dim, n_fe_degree>::get_parameters_from_command_line ()
  {
    // Booleans
    get_bool_from_options("-out_flag", out_flag);
    get_bool_from_options("-out_to_noise", out_to_noise);
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
  void StaticSPN<dim, n_fe_degree>::make_dofs ()
  {
    verbose_cout << "  dof_handler.distribute_dofs(fe);" << std::endl;
    // Relative to the DOFS
    dof_handler.distribute_dofs(fe);

    verbose_cout << "  n_blocks" << std::endl;

    const unsigned int n_blocks = n_groups * n_moments;
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_blocks);
    for (unsigned int b = 0; b < n_blocks; ++b)
      local_dofs_vector[b] = locally_owned_dofs;
    n_local_cells = GridTools::count_cells_with_subdomain_association(tria,
      tria.locally_owned_subdomain());
    n_dofs = dof_handler.n_dofs();
    n_cells = tria.n_active_cells();

    //DoFRenumbering::boost::king_ordering(dof_handler);
    //verbose_cout << "   Renumbering king_ordering. " << std::flush;
    verbose_cout << " boundary_conditions.. " << std::flush;
    // Set the Dirlichet boundary conditions of 0 value.
    constraints.clear();
    for (unsigned int c = 0; c < boundary_conditions.size(); c++)
      if (boundary_conditions[c] == 0)
      {
        AssertRelease(matrixfree_type != full_matrixfree and
                      matrixfree_type != non_diagonal,
          "Zero-Flux BC does not work in with matrix-free matrices,"
            " use Albedo_Factors=1e6 to approximate it.");
        DoFTools::make_zero_boundary_constraints(dof_handler, c, constraints);
      }
    constraints.close();
    verbose_cout << " Done!" << std::endl;

    verbose_cout << " u.." << std::endl;
    u.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
    {
      u[eig].reinit(local_dofs_vector, comm);
      u[eig].compress(VectorOperation::insert);
    }
    verbose_cout << " LALAL" << std::endl;
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void StaticSPN<dim, n_fe_degree>::assemble_system ()
  {
    M.reinit(materials, n_moments, matrixfree_type, listen_to_material_id);
    L.reinit(materials, n_moments, boundary_conditions, albedo_factors,
      matrixfree_type, listen_to_material_id);

    // Print matrices if it is needed
    print_matrices();

    memory_consumption = L.memory_consumption() + M.memory_consumption();
  }

/**
 * @brief Print matrices in the given file in a Matlab way if a file is specified
 *  with '-print_matrices_matlab' command line option.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode StaticSPN<dim, n_fe_degree>::print_matrices ()
  {
    std::string print_matrices_matlab;

    get_string_from_options("-print_matrices_matlab", print_matrices_matlab);

    if (!print_matrices_matlab.empty())
    {
      AssertRelease(matrixfree_type == full_allocated,
        "-print_matrices_matlab must be used with -allocate_matrices option");

      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);

      for (unsigned int g1 = 0; g1 < n_groups * n_moments; g1++)
        for (unsigned int g2 = 0; g2 < n_groups * n_moments; g2++)
        {
          std::string name = "L" + std::to_string(g1 + 1)
                             + std::to_string(g2 + 1);
          print_matrix_in_matlab(L.block(g1, g2), name, out, 12);
        }

      out << "L= [";
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
      {
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          out << "L" + num_to_str(g1 + 1) + num_to_str(g2 + 1) << " ";
        }
        out << ";" << std::endl;
      }
      out << "];";

      for (unsigned int g1 = 0; g1 < n_groups * n_moments; g1++)
        for (unsigned int g2 = 0; g2 < n_groups * n_moments; g2++)
        {
          std::string name = "M" + std::to_string(g1 + 1)
                             + std::to_string(g2 + 1);
          print_matrix_in_matlab(M.block(g1, g2), name, out, 12);
        }

      out << "M= [";
      for (unsigned int g1 = 0; g1 < n_groups; g1++)
      {
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          out << "M" + num_to_str(g1 + 1) + num_to_str(g2 + 1) << " ";
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
  void StaticSPN<dim, n_fe_degree>::solve_eps ()
  {

    EPSSolver<dim, n_fe_degree> solver(solver_type, L, M, n_eigenvalues, timer,
      show_eps_convergence, verbose_cout.is_active(), tria, dof_handler,
      fe);

    AssertRelease(solver_type != "slepc_2g" and solver_type != "slepc_7g",
      "The slepc_2g or slepc_7g is not implemente for the spn equations");
    // Select some solver options
    if (solver_type == "ks")
      p_init = false;
    solver.p_init = p_init;
    solver.to_init = to_init;
    solver.tol_ksp = tol_ksp;
    solver.tol_eps = tol_eps;
    solver.init_type = "multilevel-sp1";
    solver.input_file = input_file;
    solver.n_groups = n_groups;
    solver.equations = "spn";
    solver.out_file = out_file;

    solver.solve(eigenvalues, u);

    // Print the Eigenvalues calculated
    cout << std::setprecision(6) << "   The eigenvalues are:" << "     Time = "
         << timer.cpu_time()
         << " s." << std::endl;
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      cout << "      K" << i << " : " << eigenvalues[i] << std::endl;

//	L.clear();
//	M.clear();

    return;
  }

/**
 * @brief Normalize the problem to to mean neutron density power equal 1.
 * Also, calculate and print the mean values per assembly.
 */
template <int dim, int n_fe_degree>
  void StaticSPN<dim, n_fe_degree>::postprocess ()
  {

    materials.make_critical(eigenvalues[0]);

    // Erase the content of output file
    std::ofstream out(out_file.c_str(), std::ios::out);

    // Get the scalar flux phi_0
    phi.resize(n_moments);
    for (unsigned int m = 0; m < n_moments; ++m)
      phi[m].resize(n_eigenvalues);

    for (unsigned int m1 = 0; m1 < n_moments; ++m1)
    {
      for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      {
        phi[m1][eig].reinit(n_groups, comm, n_dofs,
          dof_handler.n_locally_owned_dofs());
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          for (unsigned int m = 0; m < n_moments; ++m)
          {
            phi[m1][eig].block(g).add(u_to_phi_coeff[m1][m],
              u[eig].block(g * n_moments + m));
          }
        }
      }
    }

    std::vector<std::vector<BlockVector<double> > > flux_serial(n_moments);
    for (unsigned int m1 = 0; m1 < n_moments; ++m1)
    {
      flux_serial[m1].resize(n_eigenvalues);
      for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      {
        flux_serial[m1][eig].reinit(n_groups, n_dofs);
        for (unsigned int g = 0; g < n_groups; g++)
          flux_serial[m1][eig].block(g) = phi[m1][eig].block(g);
      }
    }

    if (this_mpi_process == 0)
    {
      // Print the eigenvalues in the outFile
      print_logo(out);
      print_vector_in_file(eigenvalues, out, "The Eigenvalues are: \n", false, 7);
      out << "Problem File: " << input_file << "\n";
      print_in_file(timer.cpu_time(), out, "CPU Time: ", 4);
      out << "Equations: " << "SP" << n_moments * 2 - 1 << "\n";
      if (geo_type == "Composed")
      {
        print_in_file(n_refinements_radial, out, "Radial Refinements: "); // @suppress("Invalid arguments")
        print_in_file(n_refinements_axial, out, "Axial Refinements: "); // @suppress("Invalid arguments")
      }
      else
        print_in_file(n_refinements, out, "Global Refinements: "); // @suppress("Invalid arguments")

      print_in_file(n_fe_degree, out, "Degree of FE: ", 1);
      print_in_file(dof_handler.n_dofs(), out, "DoFs per Group: ", 1);
      print_in_file(dof_handler.n_dofs() * n_groups * n_moments, out, "Total DoFs: ", 1);
      print_in_file(tria.n_active_cells(), out, "Number of active cells: ", 1);
      print_vector_in_file(assem_per_dim, out, "Mesh Size: ", true, 1);
      print_in_file("", out, "\n", 1);
      out.close();

      // Initialize all that  is needed to iterate over dofs and cells
      QGauss<dim> quadrature_formula(n_fe_degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula,
        update_values | update_quadrature_points
        | update_volume_elements
        | update_JxW_values);

      unsigned int n_q_points = quadrature_formula.size();
      unsigned int n_cells_out = n_assemblies;
      double power_cell = 0;
      double phi_cell;
      double sigma_f;

      // Initialize and resize the vectors where it is stored the solution
      power_per_assembly.resize(n_eigenvalues,
        std::vector<double>(n_cells_out, 0.0));
      phi_per_assembly.resize(n_eigenvalues,
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
        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
        typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
          fe_values.reinit(cell);
          power_cell = 0;
          for (unsigned int g = 0; g < n_groups; ++g)
          {
            sigma_f = materials.get_sigma_f(g,
              materials.get_material_id<dim>(cell));
            fe_values.get_function_values(phi[0][eig].block(g),
              local_phi);

            phi_cell = 0.0;
            for (unsigned int q = 0; q < n_q_points; q++)
              phi_cell += local_phi[q] * fe_values.JxW(q);

            power_cell += sigma_f * phi_cell;

            phi_per_assembly[eig][g][cell->user_index()] += phi_cell;
          }

          volume += cell->measure();
          norm += std::abs(power_cell);
          power_per_assembly[eig][cell->user_index()] += power_cell;
          volume_per_assembly[cell->user_index()] += cell->measure();
        }

        norm /= volume;

        // Normalize the values of the power and fluxes per cell
        normalize_vector(power_per_assembly[eig], norm);
        for (unsigned int g = 0; g < n_groups; ++g)
          normalize_vector(phi_per_assembly[eig][g], norm);

        // Normalize the values of the power and fluxes per Rod
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
            u[eig] *= -1.0;
            for (unsigned int g = 0; g < n_groups; ++g)
            {
              normalize_vector(phi_per_assembly[eig][g], -1);
              for (unsigned int m = 0; m < n_moments; m++)
              {
                flux_serial[m][eig].block(g) *= -1;
              }
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
        for (unsigned m = 0; m < n_moments; m++)
          flux_serial[m][eig] /= norm;

        for (unsigned b = 0; b < n_moments * n_groups; b++)
        {
          u[eig].block(b) /= norm;
        }

        // Calculate the axial power distribution
        std::vector<double> power_axial;
        std::vector<std::vector<double> > phi_axial;
        std::vector<double> volume_per_plane;
        unsigned int n_assemblies_per_plane = n_assemblies
                                              / assem_per_dim[2];
        Assert(n_assemblies_per_plane > 0,
          ExcMessage("n_assemblies cannot be 0"));
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
              power_axial[i / n_assemblies_per_plane] +=
                                                         power_per_assembly[eig][i]
                                                         * volume_per_assembly[i];

              volume_per_plane[i / n_assemblies_per_plane] +=
                                                              volume_per_assembly[i];

              for (unsigned int g = 0; g < n_groups; ++g)
                phi_axial[g][i / n_assemblies_per_plane] +=
                                                            phi_per_assembly[eig][g][i]
                                                            * volume_per_assembly[i];
            }

          // Normalize axial power
          normalize_vector(power_axial, volume_per_plane);
          for (unsigned int g = 0; g < n_groups; ++g)
            normalize_vector(phi_axial[g], volume_per_plane);
        }

        // Print the Powers Distribution and the Flows
        std::ofstream out2(out_file.c_str(), std::ios::app);
        out2 << "MODE " + Utilities::int_to_string(eig + 1) << "\n";
        print_cell_distribution_in_file(dim, power_per_assembly[eig],
          assem_per_dim,
          out2,
          materials,
          "Neutron Power\n");
        if (dim == 3)
          print_vector_in_file(power_axial, out_file, "Axial power distribution:\n",
            false, 6);

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          print_cell_distribution_in_file(dim,
            phi_per_assembly[eig][g],
            assem_per_dim,
            out2,
            materials,
            "Group " + num_to_str(g + 1) + " flux\n");
          if (dim == 3)
            print_vector_in_file(phi_axial[g],
              out_file,
              "Group " + num_to_str(g + 1) + "  axial flux distribution:\n", false, 6);
        }

        out2 << "\n\n";
        out2.close();
      }
    }

    for (unsigned int m1 = 0; m1 < n_moments; ++m1)
    {
      for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      {
        for (unsigned int g = 0; g < n_groups; g++)
          phi[m1][eig].block(g) = flux_serial[m1][eig].block(g);
      }
    }

  }

/**
 * @brief Function that creates the output files .vtk.
 * It should be extended to create the power.
 */
template <int dim, int n_fe_degree>
  void StaticSPN<dim, n_fe_degree>::output_results () const
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
          phi[0][eig].block(g),
          "phi_g" + num_to_str(g + 1) + "_eig_" + num_to_str(eig + 1));

      unsigned int c = 0;
      typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
      for (; cell != tria.end(); ++cell, ++c)
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
        "Power_eig" + num_to_str(eig));
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
 *
 */
template <int dim, int n_fe_degree>
  void
  StaticSPN<dim, n_fe_degree>::output_to_noise ()
  {
    std::string filename_sta = out_file + ".sta";
    // Erase the content of output file
    std::ofstream out(filename_sta.c_str(), std::ios::out);

    int precision = 12;
    print_in_file(eigenvalues[0], out, "keff\n", precision);

    print_vector_in_file(u[0], out, "phi " + num_to_str<unsigned int>(u[0].size()) + "\n",
      true,
      precision);

    std::cout << "phi_sol[0].size(): " << u[0].size() << std::endl;

    out.close();
  }

/**
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and the time-line.
 */
template <int dim, int n_fe_degree>
  void StaticSPN<dim, n_fe_degree>::run ()
  {
    verbose_cout << "   Input files read" << std::endl;
    timer.start();

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
      verbose_cout << "  Done! " << std::endl;
    }
    else if (type_perturbation == "Mechanical_Vibration")
    {
      AssertRelease(geo_type == "Rectangular",
        "This perturbation is only implemented for rectangular geometries");
      perturbation.mechanical_vibration_static();
    }

    cout << "   Grid Done." << " Time = " << timer.cpu_time() << " s."
         << std::endl;
    cout << "      Equations:  SP" << n_moments * 2 - 1 << std::endl;
    cout << "      Number of energy groups: " << n_groups << std::endl;
    cout << "      Matrix-free type: " << enum_to_string(matrixfree_type) << std::endl;
    cout << "      Number of assemblies:" << n_assemblies << std::endl;
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
    cout << "      Number of Total DoFs: "
    << dof_handler.n_dofs() * n_groups * n_moments
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

    cout << "   Matrices assembled and sent to " + solver_type + " solver"
         << " Time = "
         << timer.cpu_time() << " s." << std::endl;

    // Get Maximum memory
    PetscLogDouble memory;
    PetscMemorySetGetMaximumUsage();

    verbose_cout << "   solve the eigenvalue problem..." << std::flush;
    solve_eps();
    verbose_cout << "  Done! " << std::endl;

    if (to_init) // we have finishes initialization
      return;

    PetscMemoryGetMaximumUsage(&memory);
    cout << "   Max Memory " << memory * 1e-6 << std::endl;
    cout << "   Memory consumption of matrix elements "
         << memory_consumption * 1e-6
         << " MB" << std::endl;

    verbose_cout << "postprocess..." << std::flush;
    postprocess();
    verbose_cout << "Done!" << std::endl;

    if (out_flag)
    {
      verbose_cout << "output_results..." << std::flush;
      output_results();
      verbose_cout << "Done!" << std::endl;
    }
    if (out_to_noise)
    {
      verbose_cout << "  output_to_noise..." << std::flush;
      output_to_noise();
      verbose_cout << "  Done!" << std::endl;
    }
  }

template class StaticSPN<1, 1> ;
template class StaticSPN<1, 2> ;
template class StaticSPN<1, 3> ;
template class StaticSPN<1, 4> ;
template class StaticSPN<1, 5> ;

template class StaticSPN<2, 1> ;
template class StaticSPN<2, 2> ;
template class StaticSPN<2, 3> ;
template class StaticSPN<2, 4> ;
template class StaticSPN<2, 5> ;

template class StaticSPN<3, 1> ;
template class StaticSPN<3, 2> ;
template class StaticSPN<3, 3> ;
template class StaticSPN<3, 4> ;
template class StaticSPN<3, 5> ;

