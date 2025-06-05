/**
 *
 * @file  femfussion.cc
 * @brief Main file of FEMFFUSION.
 *
 */

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/petsc_vector.h>

#include "../include/femffusion.h"
#include "../include/static_diffusion.h"
#include "../include/static_spn.h"
#include "../include/static_sdpn.h"
#include "../include/test.h"
#include "../include/performance.h"
#include "../include/noise/noise_diffusion.h"
#include "../include/noise/noise_spn.h"
#include "../include/noise/noise_full_spn.h"
#include "../include/time_computation.h"
#include "../include/time_computation_spn.h"
#include "../include/rom/rom_kinetics.h"
#include "../include/rom/rom_static.h"

/**
 * Print the FEMFFUSION logo.
 */
void print_logo (std::ostream &out)
{
  out << "\n";
  out << "      ___  ___  _   _  ___  ___  _ _  __  _  _  _  _ \n";
  out << "     | __|| __|| \\_/ || __|| __|| | |/ _|| |/ \\| \\| |\n";
  out << "     | _| | _| | \\_/ || _| | _| | U |\\_ \\| ( o ) \\\\ |\n";
  out << "     |_|  |___||_| |_||_|  |_|  |___||__/|_|\\_/|_|\\_|\n";
  out << "\n";
  out << "                  VERSION   1.0               \n";
  out << "\n";
}

/**
 *  Declare all the entries defined in the prm file.
 */
void prm_declare_entries (ParameterHandler &prm)
{
  prm.declare_entry("Dimension", "3", Patterns::Integer(),
    "Dimension of the problem (1, 2 or 3)");
  prm.declare_entry("N_Eigenvalues", "1", Patterns::Integer(1, 50),
    "Number of eigenvalues computed");
  prm.declare_entry("N_Refinements", "0", Patterns::Integer(0, 10),
    "Number cell refinements performed");
  prm.declare_entry("N_Refs_Radial", "0", Patterns::Integer(0, 10),
    "Number radial refinements performed");
  prm.declare_entry("N_Refs_Axial", "0", Patterns::Integer(0, 10),
    "Number axial refinements performed");
  prm.declare_entry("FE_Degree", "3", Patterns::Integer(0, 5),
    "Finite element degree");

  // Other Transport approximations
  prm.declare_entry("Transport_Appr", "Diffusion",
    Patterns::Selection("Diffusion | SPN | Full_SPN | SDPN"),
    "Transport Approximation used (Diffusion | SPN | SDPN)");
  prm.declare_entry("SPDN_Type", "Damian", Patterns::Selection("Damian | Nazari"),
    "SDPN Type (Damian | Nazari)");
  prm.declare_entry("N_SPN", "3", Patterns::Integer(1, 7),
    "Set N of the SPN equations used (1, 3 or 5)");

  // Geometry
  prm.declare_entry("Geometry_Type", "Rectangular",
    Patterns::Selection(
      "Rectangular | Hexagonal | Unstructured | Composed"),
    "Grid Type (Rectangular | Hexagonal | Unstructured | Composed)");
  prm.declare_entry("Mesh_Filename", "no.msh", Patterns::FileName(),
    ".msh File where it is the 2d mesh.");
  prm.declare_entry("Geometry_Filename", "", Patterns::FileName(),
    ".xml File where it is declared the mesh in dream way");
  prm.declare_entry("Mesh_Size", "", Patterns::Anything(),
    "Number of cells per dimension");
  prm.declare_entry("Cell_Pitch_x", " ", Patterns::Anything(),
    "Pitch of the cells in x dimension");
  prm.declare_entry("Cell_Pitch_y", "", Patterns::Anything(),
    "Pitch of the cells in y dimension");
  prm.declare_entry("Cell_Pitch_z", "", Patterns::Anything(),
    "Pitch of the cells in z dimension");
  prm.declare_entry("Geometry_Points", "", Patterns::Anything(),
    "Set the first and last existing cell every row in the reactor plant");
  prm.declare_entry("Geometry_Matrix", "", Patterns::Anything(),
    "Declare the geometry of a the reactor \n"
      " Set a matrix where, 0 means no materials, 1 nuclear fuel and 2 reflector");
  prm.declare_entry("Boundary_Conditions", "", Patterns::Anything(),
    "(LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK ) \n"
      "(0  ZeroFlow) (1  Symmetry) (2 Albedo) (3 Vacuum)");
  prm.declare_entry("Albedo_Factors", "", Patterns::Anything(),
    "Set Albedo Factor, beta, defining the BC as,  vec{n} D nabla phi = beta phi");
  prm.declare_entry("Triangulation_Filename", "", Patterns::Anything(),
    "Filename where the triangulation is saved/loaded the output");
  prm.declare_entry("Refinement_Model", "local", Patterns::Anything(),
    "Refinement model for the Composed geometry case, available options are\n"
      "Local | Uniform");

  // Output
  prm.declare_entry("Output_Filename", "out", Patterns::FileName(),
    "Filename where will be written the output");
  prm.declare_entry("Output_Flag", "false", Patterns::Bool(),
    "True/false - Make a .vtk file with the output ");
  prm.declare_entry("Print_Grid_Flag", "false", Patterns::Bool(),
    "True/false - Print the grid");
  prm.declare_entry("Out_Refinements", "1", Patterns::Integer(),
    "Number of Refinements per cell in the .vtk output");
  prm.declare_entry("Out_Interval", "1", Patterns::Integer(),
    "Every this number of steps we will create a .vtk file");
  prm.declare_entry("Output_To_Noise", "false", Patterns::Bool(),
    "True/false - Output static values to noise");

  // Time Step
  // FIXME Check if Time_Step is defined but the Transient==false is not set
  prm.declare_entry("Transient", "false", Patterns::Bool(),
    "True/false - Compute the time-dependent computation");
  prm.declare_entry("Time_Step", "0.0", Patterns::Double(),
    "Delta time iteration");
  prm.declare_entry("Time_End", "0.0", Patterns::Double(),
    "Final time of computation");

  // Save Static Calculation
  prm.declare_entry("Save_Static", "false", Patterns::Bool(),
    "Save_Static");
  prm.declare_entry("STA_Filename", "none.sta", Patterns::FileName(),
    "Load Steady state calculation from a previous one.");

  // Solver Options
  prm.declare_entry("Renumbering", "Reversed_Cuthill_McKee",
    Patterns::Selection("Cuthill_McKee | Reversed_Cuthill_McKee"
      " | Minimum_degree | King_Ordering | None"),
    "Renumbering Type (Cuthill_McKee | Reversed_Cuthill_McKee"
      " | Minimum_degree | King_Ordering | None)");
  prm.declare_entry("Solver_Type", "bifpam",
    Patterns::Selection(
      "power_it | bifpam | slepc_2g | slepc_7g | ks | gd | newton"),
    "Solver Type (power_it | bifpam | slepc_2g | slepc_7g | ks | gd | newton)");
  prm.declare_entry("EPS_Tolerance", "1e-7", Patterns::Double(),
    "Relative tolerance of the Eigenvalue problem solver");
  prm.declare_entry("KSP_Tolerance", "1e-9", Patterns::Double(),
    "Relative tolerance of the linear system solver");
  prm.declare_entry("Static_KSP_Tolerance", "false", Patterns::Bool(),
    "True/false - Use an static tolerance for the associated linear systems");
  prm.declare_entry("Adjoint", "false", Patterns::Bool(),
    "True/false - Do the adjoint problem");
  prm.declare_entry("Spectral_index", "false", Patterns::Bool(),
    "True/false - Compute the spectral_index ");
  prm.declare_entry("Matrix_Free_Type", "non_diagonal",
    Patterns::Selection(
      "full_allocated | non_diagonal | full_matrixfree"),
    "Type of the matrix-free methodology employed (full_allocated | non_diagonal | full_matrixfree)");
  prm.declare_entry("P_Init", "true", Patterns::Bool(),
    "True/false - Use FE_Degree=1 Initialization");

  // Cross Sections
  prm.declare_entry("XSEC_Type", "XS2G",
    Patterns::Selection("XS2G | XSEC | XML | Valkin"),
    "XS type file defined (XS2G | XSEC | XML | Valkin)");
  prm.declare_entry("XSECS_Filename", "XSECS", Patterns::FileName(),
    "Filename where it is stored the material cross section");
  prm.declare_entry("Energy_Groups", "2", Patterns::Integer(1, 120),
    "Set the number of neutron energy groups defined in the calculation");
  prm.declare_entry("Precursors_Flag", "true", Patterns::Bool(),
    "true/false - Select if it is considered the precursors");
  prm.declare_entry("PREC_Filename", "none", Patterns::FileName(),
    "Filename where it is stored the precursors data");

  // Noise Calculation
  prm.declare_entry("Noise_Calculation", "false", Patterns::Bool(),
    "True/false - Activate Noise Calculation");
  prm.declare_entry("DS_Filename", "", Patterns::FileName(),
    "Filename where it is stored the perturbation data");
  prm.declare_entry("DYN_Filename", "", Patterns::FileName(),
    "Filename where it is stored some dynamic data");
  prm.declare_entry("RESULTS_Filename", "", Patterns::FileName(),
    "Filename where the output results are stored in a Matlab .mat structure");
  prm.declare_entry("Perturbation_Type", "Cell_Wise",
    Patterns::Selection("Cell_Wise | Borders | BordersHex"),
    " Type of the perturbation (Cell_Wise or Borders)");
  prm.declare_entry("PC_Noise", "gauss_seidel", Patterns::Anything(),
    "Preconditioner used for the noise complex linear system solver");
  prm.declare_entry("KSP_Noise_Tolerance", "1e-8", Patterns::Double(),
    "Relative tolerance of the noise complex linear system solver");

  // Time Step
  prm.declare_entry("Time_Delta", "0.01", Patterns::Double(),
    "Delta time of the iteration");
  prm.declare_entry("Time_Delta_Updating", "0.0", Patterns::Double(),
    "Delta time for updating the modes");
  prm.declare_entry("Time_End", "1.0", Patterns::Double(),
    "Final time of the computation");
  prm.declare_entry("TS_Solver", "petsc",
    Patterns::Selection("euler | arkode | dealii | petsc"),
    "TS_Solver defined euler | arkode | dealii | petsc");
  prm.declare_entry("Type_Time_Preconditioner", "fixed",
    Patterns::Selection("fixed | good-broyden | bad-broyden "),
    "Time Preconditioner defined fixed | good-broyden | bad-broyden");
  prm.declare_entry("Initial_Time_Preconditioner", "gs-cgilu",
    Patterns::Selection("gs-cgilu | gs-ilu | diagonal"),
    "Time Preconditioner defined fixed | good-broyden | bad-broyden");
  prm.declare_entry("Print_Time_Dependent_Data", "false", Patterns::Bool(),
    "True/ false - Print Radial Output data each time step");

  // Time Variables
  prm.declare_entry("Frequency", "1.0", Patterns::Double(),
    "Frequency of the vibration in Hz");
  prm.declare_entry("Amplitude", "0.0", Patterns::Double(),
    "Amplitude instability in %");
  prm.declare_entry("Amplitudes", "", Patterns::Anything(),
    "Amplitude instabilities in %");
  prm.declare_entry("Out_Phase", "3.141592653589793", Patterns::Double(),
    "Phase for the out of phase");
  prm.declare_entry("Slope_Up", "0.0", Patterns::Anything(),
    "Slope for the ramp perturbation");
  prm.declare_entry("Slope_Down", "0.0", Patterns::Anything(),
    "Slope for the ramp perturbation");
  prm.declare_entry("Cut_Time", "100.0", Patterns::Anything(),
    "Cut Time for the ramp perturbation");
  prm.declare_entry("Material_Changing", "", Patterns::Anything(),
    "Material where the instability is inserted");
//  prm.declare_entry("Material_Changing_1", "0", Patterns::Integer(),
//    "Second Material where the instability is inserted");
  prm.declare_entry("Group_Changing", "0", Patterns::Integer(),
    "Second Material where the instability is inserted");
  prm.declare_entry("XS_Name", "Sigma_f", Patterns::Anything(),
    "Cross Section where the Sigma_f/Sigma_a...");

  prm.declare_entry("Distributed_Time_Scheme", "Implicit-Exponential",
    Patterns::Selection(
      "Implicit-Exponential | Semi-Implicit-Exponential | Semi-Implicit-Euler"),
    "Time Scheme used to solve the distributed time equation");
  prm.declare_entry("Type_Perturbation", "None",
    Patterns::Selection(
      "Flux_Distributed | Single_Material | Out_Of_Phase | Ramp_Two_Mats | Step_Change_Material "
        "| Rods | AECL | Mechanical_Vibration | Read_XS_File | Read_XML_File | C5G7-TD1.1 | Random_XS| None"),
    "Distribution of the instability: Flux_Distributed or Single_Material");
  prm.declare_entry("Perturbation_Function", "Constant",
    Patterns::Selection("Constant | Ramp | Sinus | Ramp_hex | Noise_7g "),
    "Type of instability Constant, Ramp, Sinus, Ramp_hex or Noise_7g");
  prm.declare_entry("Spatial_Modes", "lambda",
    Patterns::Selection("lambda | alpha | gamma"),
    "Type of modes equation: lambda or alpha");
  prm.declare_entry("Bar_Filename", "no.bar", Patterns::FileName(),
    "Filename where it is defined the movement of the Rods");
  prm.declare_entry("Rod_Cusping_Method", "volhom",
    Patterns::Selection("volhom | fluxwei "),
    "Rod cusping method: volhom or fluxwei");

  prm.declare_entry("Save_Time", "false", Patterns::Bool(),
    "# True/false - Activate Save_Time");
  prm.declare_entry("Load_Time", "false", Patterns::Bool(),
    "# True/false - Activate Load_Time");
  prm.declare_entry("Reinit_File", "nofile", Patterns::FileName(),
    "# Filename where the reinit is saved/loaded");

  // Perturbations
  prm.declare_entry("Vibrating_Material", "0", Patterns::Integer(),
    "Material of the  assembly that is vibrating");
  prm.declare_entry("Static_Position", "", Patterns::Anything(),
    "Static position of the assembly that is vibrating (x_left x_right y_left y_right)");
  prm.declare_entry("Direction", "0", Patterns::Integer(),
    "Direction of the vibration (x=0, y=1, z=2)");
  prm.declare_entry("PseudoStatic", "false", Patterns::Bool(),
    "True/false - Make a pseudostatic calculation");
  prm.declare_entry("Read_XS_Filename", " ",
    Patterns::FileName(Patterns::FileName::FileType::input),
    "Filename where the XS of the read_xs_file perturbation.");
  prm.declare_entry("Read_XML_Filename", " ",
    Patterns::FileName(Patterns::FileName::FileType::input),
    "Filename where the XS of the read_xml_file perturbation.");
  prm.declare_entry("XS_Perturbation_Fraction", "0.0", Patterns::Double(),
    "True/false - Make a pseudostatic calculation");

  // ROM variables
  prm.declare_entry("ROM_Static", "false", Patterns::Bool(),
    "# True/false - Activate ROM Static Calculation");
  prm.declare_entry("ROM_Transient", "false", Patterns::Bool(),
    "# True/false - Activate ROM Calculation");
  prm.declare_entry("ROM_Type_Snapshots", "", Patterns::Anything(),
    "Type of Snapshots used in ROM calculation");
  prm.declare_entry("ROM_Time_Break_Snapshots", "", Patterns::Anything(),
    "Time Break of Snapshots used in ROM calculation");
  prm.declare_entry("N_Snapshots", "2", Patterns::Integer(),
    "Number of snapshots for the ROM method");
  prm.declare_entry("N_Test", "50", Patterns::Integer(),
    "Number of test for the Static ROM Method");
  prm.declare_entry("ROM_Slope_Up", "0.0", Patterns::Anything(),
    "Slope for the ramp perturbation");
  prm.declare_entry("ROM_Slope_Down", "0.0", Patterns::Anything(),
    "Slope for the ramp perturbation");
  prm.declare_entry("ROM_Cut_Time", "100.0", Patterns::Anything(),
    "Cut Time for the ramp perturbation");
  prm.declare_entry("ROM_Group_Wise", "Monolithic",
    Patterns::Selection("Group_Wise | Monolithic"),
    "Type of Snapshots used in ROM calculation");

  // LUPOD
  prm.declare_entry("LUPOD_Type", "POD",
    Patterns::Selection("POD | LUPOD | LUPOD_ext | Random | FEM1"),
    "Activate LUPOD or LUPOD_ext techniques to optimize ROM");
  prm.declare_entry("Epsilon_N", "0.0", Patterns::Double(0, 1.0),
    "Epsilon_M of LUPOD technique");
  prm.declare_entry("Epsilon_M", "0.0", Patterns::Double(0, 1.0),
    "Epsilon_M of LUPOD technique");
  prm.declare_entry("N_LUPOD_Points", "0",
    Patterns::Integer(),
    "Number of LUPOD points retained in the cross products");
}

/**
 * @brief Declare prm entries for dynamic file.
 */
void prm_dyn_entries (ParameterHandler &prm)
{
  prm.declare_entry("Neutron_Velocities", "", Patterns::Anything(),
    "Per Group neutron Velocity");
  prm.declare_entry("Frequency", "0.0", Patterns::Double(0.0, 1e6),
    "Frequency of the perturbation and output");
  prm.declare_entry("Lambda_eff", "0.0", Patterns::Double(0.0, 1000),
    "Lambda  effective of neutron precursors ");
  prm.declare_entry("Beta_eff", "0.0", Patterns::Double(0.0, 0.5),
    "Beta effective of neutron precursors");
}

/**
 *  Main function.
 */
int main (int argc,
  char **argv)
{
  try
  {
    using namespace dealii;
    deallog.depth_console(0);

    bool verbose = false;
    bool silent = false;

    // Default input Filename
    std::string input_file = "InputFile.prm";

    // Parse command line in order to find the input_file
    // and remove this option in order not to raise a SLEPC Warning.
    for (int i = 0; i < argc; i++)
    {
      if (!strcmp("-f", argv[i]) or !strcmp("--file", argv[i]))
      {
        input_file = argv[i + 1];
        argv[i][0] = 0;
      }
      if (!strcmp("-v", argv[i]) or !strcmp("--verbose", argv[i]))
      {
        verbose = true;
        argv[i][0] = 0;
      }
      if (!strcmp("-s", argv[i]) or !strcmp("--silent", argv[i]))
      {
        silent = true;
        argv[i][0] = 0;
      }
      if (!strcmp("-t", argv[i]) or !strcmp("--test", argv[i]))
      {
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
          1);
        {
          run_tests();
          return 0;
        }
      }
      if (!strcmp("-p", argv[i]) or !strcmp("--performance", argv[i]))
      {
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        {
          run_performance_matmult();
          run_performance_solvers();
          return 0;
        }
      }
    }

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    {
      if (not silent and Utilities::MPI::this_mpi_process(PETSC_COMM_WORLD) == 0)
        print_logo(std::cout);

      //Check if the input_file exists also in release mode
      if (fexists(input_file) == false)
      {
        if (Utilities::MPI::this_mpi_process(PETSC_COMM_WORLD) == 0)
        {
          std::cout << "ERROR!: Input File .prm Does NOT exist"
                    << std::endl;
          std::cout << "Use -f InputFilename.rm  option "
                    << std::endl;
          // Exit with error
          MPI_Finalize();
          exit(1);
        }
      }

      // ParameterHandler definition. The get of the parameters is done in the
      ParameterHandler prm;
      prm_declare_entries(prm);
      prm.parse_input(input_file);

      // get Dimension and fe_degree of the problem
      int dim = prm.get_integer("Dimension");
      unsigned int fe_degree = prm.get_integer("FE_Degree");
      std::string transport = prm.get("Transport_Appr");
      get_string_from_options("-transport", transport);
      lower_case(transport);
      get_uint_from_options("-fe_degree", fe_degree);

      // Get transient options
      bool transient = prm.get_bool("Transient");
      get_bool_from_options("-transient", transient);

      // Get Noise Calculation
      bool noise = prm.get_bool("Noise_Calculation");

      // Get ROM Calculation
      bool rom_transient = prm.get_bool("ROM_Transient");

      // Get ROM Calculation
      bool rom_static = prm.get_bool("ROM_Static");

      AssertRelease(dim >= 0 and dim <= 3,
        "Spatial dimension not valid,  dim = " + num_to_str(dim));

      // The template parameters must be a constant:
      // This way, we ensure the compilation of all the cases
      // and only create/run the good one
      /* ------------------------------------------------------ */
      /* --------------------    Diffusion  ------------------- */
      /* ------------------------------------------------------ */
      if (transport == "diffusion")
      {
        if (dim == 1)
        {
          if (fe_degree == 1)
          {
            StaticDiffusion<1, 1> static_prob(prm, input_file, verbose, silent, false);

            if (transient)
              TimeNeutronDiffusion<1, 1> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<1, 1> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<1, 1> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<1, 1> rom_pro(prm, static_prob, verbose, silent);

          }
          else if (fe_degree == 2)
          {
            StaticDiffusion<1, 2> static_prob(prm, input_file, verbose, silent, false);

            if (transient)
              TimeNeutronDiffusion<1, 2> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<1, 2> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<1, 2> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<1, 2> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticDiffusion<1, 3> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<1, 3> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<1, 3> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<1, 3> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<1, 3> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticDiffusion<1, 4> static_prob(prm, input_file, verbose, silent, false);

            if (transient)
              TimeNeutronDiffusion<1, 4> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<1, 4> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<1, 4> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<1, 4> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticDiffusion<1, 5> static_prob(prm, input_file, verbose, silent, false);

            if (transient)
              TimeNeutronDiffusion<1, 5> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<1, 5> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<1, 5> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<1, 5> rom_pro(prm, static_prob, verbose, silent);
          }
        }
        else if (dim == 2)
        {
          if (fe_degree == 1)
          {
            StaticDiffusion<2, 1> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<2, 1> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<2, 1> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<2, 1> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<2, 1> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticDiffusion<2, 2> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<2, 2> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<2, 2> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<2, 2> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<2, 2> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticDiffusion<2, 3> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<2, 3> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<2, 3> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<2, 3> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<2, 3> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticDiffusion<2, 4> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<2, 4> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<2, 4> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<2, 4> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<2, 4> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticDiffusion<2, 5> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<2, 5> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<2, 5> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<2, 5> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<2, 5> rom_pro(prm, static_prob, verbose, silent);
          }
        }
        else if (dim == 3)
        {
          if (fe_degree == 1)
          {
            StaticDiffusion<3, 1> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<3, 1> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<3, 1> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<3, 1> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<3, 1> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {

            StaticDiffusion<3, 2> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<3, 2> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<3, 2> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<3, 2> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<3, 2> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticDiffusion<3, 3> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<3, 3> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<3, 3> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<3, 3> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<3, 3> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticDiffusion<3, 4> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<3, 4> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<3, 4> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<3, 4> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<3, 4> rom_pro(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticDiffusion<3, 5> static_prob(prm, input_file, verbose, silent, false);
            if (transient)
              TimeNeutronDiffusion<3, 5> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseDiffusion<3, 5> noise_prob(prm, static_prob, verbose, silent);
            else if (rom_transient)
              ROMKinetics<3, 5> rom_pro(prm, static_prob, verbose, silent);
            else if (rom_static)
              ROMStatic<3, 5> rom_pro(prm, static_prob, verbose, silent);
          }
        }
      }
      /* ------------------------------------------------------ */
      /* --------------------    SPN  ------------------------- */
      /* ------------------------------------------------------ */
      if (transport == "spn")
      {
        AssertRelease(rom_static == false, "ROM STATIC not implemented with SPN");
        if (dim == 1)
        {
          if (fe_degree == 1)
          {
            StaticSPN<1, 1> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<1, 1> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<1, 1> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticSPN<1, 2> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<1, 2> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<1, 2> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticSPN<1, 3> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<1, 3> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<1, 3> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticSPN<1, 4> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<1, 4> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<1, 4> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticSPN<1, 5> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<1, 5> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<1, 5> noise_prob(prm, static_prob, verbose, silent);
          }
        }
        else if (dim == 2)
        {
          if (fe_degree == 1)
          {
            StaticSPN<2, 1> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<2, 1> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<2, 1> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticSPN<2, 2> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<2, 2> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<2, 2> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticSPN<2, 3> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<2, 3> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<2, 3> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticSPN<2, 4> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<2, 4> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<2, 4> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticSPN<2, 5> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<2, 5> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<2, 5> noise_prob(prm, static_prob, verbose, silent);
          }
        }
        else if (dim == 3)
        {
          if (fe_degree == 1)
          {
            StaticSPN<3, 1> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<3, 1> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<3, 1> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticSPN<3, 2> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<3, 2> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<3, 2> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticSPN<3, 3> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<3, 3> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<3, 3> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticSPN<3, 4> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<3, 4> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<3, 4> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticSPN<3, 5> static_prob(prm, input_file, verbose, silent);
            if (transient)
              TimeNeutronSPN<3, 5> time_pro(prm, static_prob, verbose, silent);
            else if (noise)
              NoiseSPN<3, 5> noise_prob(prm, static_prob, verbose, silent);
          }
        }
      }
      /* ------------------------------------------------------ */
      /* ---------------  Simplified Double PN  --------------- */
      /* ------------------------------------------------------ */
      if (transport == "sdpn")
      {
        if (dim == 1)
        {
          if (fe_degree == 1)
          {
            StaticSDPN<1, 1> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticSDPN<1, 2> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticSDPN<1, 3> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticSDPN<1, 4> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticSDPN<1, 5> static_prob(prm, input_file, verbose, silent);
          }
        }
        else if (dim == 2)
        {
          if (fe_degree == 1)
          {
            StaticSDPN<2, 1> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticSDPN<2, 2> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticSDPN<2, 3> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticSDPN<2, 4> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticSDPN<2, 5> static_prob(prm, input_file, verbose, silent);
          }
        }
        else if (dim == 3)
        {
          if (fe_degree == 1)
          {
            StaticSDPN<3, 1> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticSDPN<3, 2> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticSDPN<3, 3> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticSDPN<3, 4> static_prob(prm, input_file, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticSDPN<3, 5> static_prob(prm, input_file, verbose, silent);
          }
        }
      }
      /* ----------------------------------------------------- */
      /* -----------------   FULL SPN  ----------------------- */
      /* ----------------------------------------------------- */
      if (transport == "full_spn")
      {
        AssertRelease(transient == false, "Transient in Full_SPN not Implemented yet.");
        if (dim == 1)
        {
          if (fe_degree == 1)
          {
            StaticFullSPN<1, 1> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<1, 1> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticFullSPN<1, 2> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<1, 2> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticFullSPN<1, 3> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<1, 3> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticFullSPN<1, 4> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<1, 4> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticFullSPN<1, 5> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<1, 5> noise_prob(prm, static_prob, verbose, silent);
          }
        }
        else if (dim == 2)
        {
          if (fe_degree == 1)
          {
            StaticFullSPN<2, 1> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<2, 1> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticFullSPN<2, 2> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<2, 2> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticFullSPN<2, 3> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<2, 3> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticFullSPN<2, 4> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<2, 4> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticFullSPN<2, 5> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<2, 5> noise_prob(prm, static_prob, verbose, silent);
          }
        }
        else if (dim == 3)
        {
          if (fe_degree == 1)
          {
            StaticFullSPN<3, 1> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<3, 1> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 2)
          {
            StaticFullSPN<3, 2> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<3, 2> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 3)
          {
            StaticFullSPN<3, 3> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<3, 3> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 4)
          {
            StaticFullSPN<3, 4> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<3, 4> noise_prob(prm, static_prob, verbose, silent);
          }
          else if (fe_degree == 5)
          {
            StaticFullSPN<3, 5> static_prob(prm, input_file, verbose, silent);
            if (noise)
              NoiseFullSPN<3, 5> noise_prob(prm, static_prob, verbose, silent);
          }
        }
      }

    }

    // If no exceptions are thrown,
    // then we tell the program to stop monkeying around and exit nicely:
    if (not silent and Utilities::MPI::this_mpi_process(PETSC_COMM_WORLD) == 0)
      std::cout << std::endl << "Job done." << std::endl;
  }

  // All the while, we are watching out if any exceptions should
  // have been generated. If that is so, we panic...
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl << exc.what()
              << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl << "Aborting!"
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
