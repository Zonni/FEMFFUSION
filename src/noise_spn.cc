/**
 * @file   noise_spn.cc
 * @brief Implementation of the neutron noise
 */

#include <deal.II/base/parameter_handler.h>

#include <deal.II/base/numbers.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <complex>
#include <cmath>
#include <math.h>  //  For M_PI and atan2
#include <string>

#include "../include/matrix_operators/matrix_operators_noise_spn.h"
#include "../include/complex_solver.h"
#include "../include/static_diffusion.h"
#include "../include/femffusion.h"
#include "../include/matlab_io.h"
#include "../include/utils.h"
#include "../include/printing.h"
#include "../include/noise_spn.h"

using namespace dealii;
//using namespace std::complex_literals;
typedef std::complex<double> complex;

/**
 * @brief Constructor of the main class Noise.
 *
 */
template <int dim, int n_fe_degree>
  NoiseSPN<dim, n_fe_degree>::NoiseSPN (ParameterHandler &prm,
    StaticSPN<dim, n_fe_degree> &static_problem,
    const bool verbose,
    const bool silent) :
      verbose_cout(std::cout, verbose),
      cout(std::cout, !silent),
      comm(static_problem.comm),
      n_groups(static_problem.n_groups),
      n_moments(static_problem.n_moments),
      n_blocks(2 * static_problem.n_groups * static_problem.n_moments),
      n_assemblies(static_problem.n_assemblies),
      n_dofs(static_problem.n_dofs),
      assem_per_dim(static_problem.assem_per_dim),
      assembly_pitch(static_problem.assembly_pitch),
      // boundary_conditions(static_problem.boundary_conditions),
      materials(static_problem.materials),
      tria(static_problem.tria),
      fe_block(
        static_problem.fe),
      dof_handler_block(static_problem.dof_handler),
      //dof_handler_block(tria),
      boundary_conditions(static_problem.boundary_conditions),
      albedo_factors(static_problem.albedo_factors),
      A(comm, dof_handler_block, constraints),
      phi_per_assembly(static_problem.phi_per_assembly[0]),
      pert(materials, dim, verbose),
      timer(static_problem.timer)
  {
    cout << std::endl;
    cout << "   ---------------------" << std::endl;
    cout << "     NOISE CALCULATION  " << std::endl;
    cout << "   ---------------------" << std::endl << std::endl;

    // Copy some values From Static
    n_out_ref = static_problem.n_out_ref;
    out_file = static_problem.out_file;
    out_flag = static_problem.out_flag;
    tol_ksp_noise = prm.get_double("KSP_Noise_Tolerance");
    pc_complex = prm.get("PC_Noise");

    // ---------------------------------------------------------------------------------
    // Make reactor critical
    keff = static_problem.eigenvalues[0];

    // Copy phi critic and normalize it
    u_crit.reinit(n_blocks, comm, n_dofs, n_dofs);
    for (unsigned int b = 0; b < n_blocks / 2; b++)
      u_crit.block(2 * b) = static_problem.u[0].block(b);

    if (static_problem.geo_type == "Rectangular"
        or static_problem.geo_type == "Hexagonal"
        or static_problem.geo_type == "Composed")
    {
      // -------------------------------------------------------------------
      // Perturbation
      std::string dxs_file = prm.get("DS_Filename");
      get_string_from_options("-dxs_file", dxs_file);
      verbose_cout << "dxs_file " << dxs_file << std::endl;
      std::string xs_type = prm.get("XSEC_Type");
      std::string pert_type = prm.get("Perturbation_Type");
      verbose_cout << " pert.reinit... " << std::flush;
      pert.reinit(dxs_file, xs_type, pert_type);
      verbose_cout << " Done! " << std::endl;

      // -------------------------------------------------------------------
      // Dynamic data
      std::string dyn_file = prm.get("DYN_Filename");
      verbose_cout << "dyn_file " << dyn_file << std::endl;

      ParameterHandler prm_dyn;
      prm_dyn_entries(prm_dyn);
      prm_dyn.parse_input(dyn_file); // @suppress("Ambiguous problem")
      const double freq = prm_dyn.get_double("Frequency");
      AssertRelease(freq > 0.0, "frequency must be >0");
      double _beta, _lambda;
      if (materials.get_n_precursors() == 0)
      {
        _beta = prm_dyn.get_double("Beta_eff");
        AssertRelease(_beta > 0.0, "beta must be >0");
        _lambda = prm_dyn.get_double("Lambda_eff");
        AssertRelease(_lambda > 0.0, "lambda must be >0");

        materials.set_n_precursors(1);
        materials.set_beta_prec(0, _beta);
        materials.set_lambda_prec(0, _lambda);
        materials.set_default_spectra();

        // Neutron Velocities
        std::vector<double> velocities;
        parse_vector(prm_dyn.get("Neutron_Velocities"), velocities,
          n_groups);
        AssertRelease(velocities.size() == n_groups,
          "Not correct number of Neutron_Velocities entries");
        for (unsigned int g = 0; g < n_groups; g++)
        {
          AssertRelease(velocities[g] > 0.0,
            "Each Velocity must be greater than 0");
          materials.set_velocity(g, velocities[g]);
        }
      }

      //omega = 2 * M_PI * freq;
      pert.set_frequency(freq);

    }
    else
      AssertRelease(false, "Invalid geo_type used in noise calculation");

    // Output in Matlab Format
    results_file = prm.get("RESULTS_Filename");
    verbose_cout << "results_file " << results_file << std::endl;

    // Detectors
    detectors_file = prm.get("Detectors_Out_File");
    if (detectors_file != "")
    {
      verbose_cout << "detectors_file " << detectors_file << std::endl;

      parse_vector(prm.get("Detectors_Pos_Idx"), detectors_idx);
      for (unsigned int i = 0; i < detectors_idx.size(); i++)
        detectors_idx[i]--; // Convert to C-style idx
      verbose_cout << "detectors_pos: " << std::flush;
      if (verbose)
        print_vector(detectors_idx);

      // Z position of detectors
      parse_vector(prm.get("Detector_levels"), detector_levels);
      n_detectors = detector_levels.size() / 2 * detectors_idx.size();
      verbose_cout << "n_detectors " << n_detectors << std::endl;
      for (unsigned int i = 0; i < detector_levels.size(); i++)
        detector_levels[i]--; // Convert to C-style idx
      verbose_cout << "detector_levels: " << std::flush;
      if (verbose)
        print_vector(detector_levels);
    }

    // Get changes in the parameters through the command line
    get_parameters_from_command_line();

    this->run();
  }

/**
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::get_parameters_from_command_line ()
  {

    // Booleans
    get_bool_from_options("-out_flag", out_flag);

    // Integers
    //get_uint_from_options("-n_refinements", n_refinements);

    // Reals
    get_double_from_options("-tol_ksp_noise", tol_ksp_noise);
    double freq = 0.0;
    get_double_from_options("-frequency", freq);
    if (freq > 1e-10)
      pert.set_frequency(freq);

    // String
    get_string_from_options("-pc_complex", pc_complex);
    get_string_from_options("-dtc_file", detectors_file);
    //get_string_from_options("-refinement_model", refinement_model);
    // lower_case(refinement_model);
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::assemble_system ()
  {

    system_rhs.reinit(n_blocks, comm, n_dofs, n_dofs);
    delta_u.reinit(n_blocks, comm, n_dofs, n_dofs);

    assemble_rhs();

    MatrixFreeType matrixfree_type_A = non_diagonal;
    //get_int_from_options("-matrixfree_type_A", matrixfree_type_A);
    get_enum_from_options("-matrixfree_type_A", matrixfree_type_A);
    verbose_cout << "A.reinit ..." << std::flush;
    cout << "      Matrixfree_type_A  " << enum_to_string(matrixfree_type_A) << std::endl;
    A.reinit(materials, pert, n_moments, boundary_conditions, albedo_factors,
      matrixfree_type_A);
    verbose_cout << "  Done!" << std::endl;

    memory_consumption = A.memory_consumption();

//    cout << "      Real_block (1,1) norm infty: " << A.block(0,0).linfty_norm() << std::endl;
//    cout << "      Imag block (1,1) norm infty: " << A.block(1,0).linfty_norm() << std::endl;
//
//    cout << "      Real_block (2,2) norm infty: " << A.block(2,2).linfty_norm() << std::endl;
//    cout << "      Imag block (2,2) norm infty: " << A.block(3,2).linfty_norm() << std::endl;
//
//    cout << "      Real_block (1,2) norm infty: " << A.block(0,2).linfty_norm() << std::endl;
//    cout << "      Imag block (1,2) norm infty: " << A.block(1,2).linfty_norm() << std::endl;
//
//    cout << "      Real_block (2,1) norm infty: " << A.block(2,0).linfty_norm() << std::endl;
//    cout << "      Imag block (2,1) norm infty: " << A.block(3,0).linfty_norm() << std::endl;

  }

/**
 * @brief Compute Noise equation rhs: B*\phi_0
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::assemble_rhs ()
  {
    NoiseBMatrixSPN<dim, n_fe_degree> B(comm, dof_handler_block, constraints);
    verbose_cout << "B.reinit ..." << std::flush;
    MatrixFreeType matrixfree_type_B = full_allocated;

    if (pert.pert_type == "Cell_Wise")
      matrixfree_type_B = full_matrixfree;

    B.reinit(materials, pert, n_moments, matrixfree_type_B);
    verbose_cout << "  Done!" << std::endl;

    verbose_cout << "  Matrixfree_type_B " << matrixfree_type_B << std::endl;
    cout << "      u_crit  norm: " << u_crit.l2_norm() << std::endl;
    verbose_cout << " B.vmult ..." << std::flush;
    B.vmult(system_rhs, u_crit);
    verbose_cout << "  Done!" << std::endl;
  }

/**
 * @brief Solve the eigenvalue problem with the power iteration method.
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::solve ()
  {

    Timer timer_solve;
    timer_solve.start();
    bool show_ksp_convergence = false;

    ComplexSolver<dim, n_fe_degree> complex_solver(A, timer_solve, materials, pert,
      show_ksp_convergence);
    complex_solver.tol_ksp = tol_ksp_noise;
    complex_solver.pc_complex = pc_complex;

    cout << "      Complex preconditioner: " << complex_solver.pc_complex << std::endl;

    complex_solver.solve(system_rhs, delta_u);

    cout << "   Complex system solved. It took " << timer_solve.cpu_time()
    << " s and "
    << complex_solver.n_ksp_iterations << " iterations to solve the linear system."
    << std::endl;

    cout << "      system_rhs norm: " << system_rhs.l2_norm() << " " << std::endl;
    cout << "      delta_phi norm: " << delta_u.l2_norm() << std::endl;

  }

/**
 * @brief Normalize the problem to to mean neutron density power equal 1.
 * Also, calculate and print the mean values per assembly.
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::postprocess ()
  {
    // Erase the content of output file
    std::ofstream out(out_file.c_str(), std::ios::app);

    // Initialize all that  is needed to iterate over dofs and cells
    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler_block.get_fe(), quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values);

    unsigned int n_q_points = quadrature_formula.size();
    unsigned int n_cells_out = n_assemblies;

    std::vector<std::vector<std::vector<double> > > delta_cell_u(n_groups,
      std::vector<std::vector<double>>(n_moments,
        std::vector<double>(2)));

    // Initialize and resize the vectors where it is stored the solution
    noise_per_assembly.resize(n_groups, std::vector<complex>(n_cells_out, 0.0));

    std::vector<double> volume_per_assembly(n_cells_out, 0.0);

    std::vector<std::vector<std::vector<std::vector<double>>>> local_delta_u(n_groups,
      std::vector<std::vector<std::vector<double>>>(n_moments,
        std::vector<std::vector<double>>(2,
          std::vector<double>(n_q_points))));

    // Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler_block.begin_active(),
        endc = dof_handler_block.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_moments; m++)
        {
          delta_cell_u[g][m][0] = 0.0;
          delta_cell_u[g][m][1] = 0.0;
        }

      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_moments; m++)
        {
          fe_values.get_function_values(delta_u.block(2 * A.gm_to_b(g, m)),
            local_delta_u[g][m][0]);
          fe_values.get_function_values(delta_u.block(2 * A.gm_to_b(g, m) + 1),
            local_delta_u[g][m][1]);
        }

      for (unsigned int q = 0; q < n_q_points; q++)
      {
        for (unsigned int g = 0; g < n_groups; g++)
          for (unsigned int m = 0; m < n_moments; m++)
          {
            delta_cell_u[g][m][0] += local_delta_u[g][m][0][q]
                                     * fe_values.JxW(q);
            delta_cell_u[g][m][1] += local_delta_u[g][m][1][q]
                                     * fe_values.JxW(q);
          }
      }

      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_moments; m++)
          noise_per_assembly[g][cell->user_index()] +=
                                                       u_to_phi_coeff[0][m]
                                                       * complex(delta_cell_u[g][m][0],
                                                         delta_cell_u[g][m][1]);

      volume_per_assembly[cell->user_index()] += cell->measure();
    }

    // Remove the added volume to noise
    for (unsigned int g = 0; g < n_groups; ++g)
      normalize_vector(noise_per_assembly[g], volume_per_assembly);

    // Print the noise data
    std::ofstream out_stream(out_file.c_str(), std::ios::app);
    const unsigned int presicion = 6;
    out_stream << "NOISE RESULTS\n";
    for (unsigned int gr = 0; gr < n_groups; gr++)
      print_cell_distribution_in_file(dim, noise_per_assembly[gr],
        assem_per_dim, out_stream, materials,
        "Flux Noise Group " + num_to_str(gr + 1) + "\n", false,
        presicion);
    out_stream << "\n\n";
    out_stream.close();

    A.clear();
  }

#ifdef MATIO
  /**
   * @brief Function that creates the output files .vtk.
   * It should be extended to create the power.
   */
  template <int dim, int n_fe_degree>
  void
  NoiseSPN<dim, n_fe_degree>::output_results_matlab () const
  {

    mat_t *matfp;
    matfp = Mat_CreateVer(results_file.c_str(), NULL, MAT_FT_DEFAULT);
    verbose_cout << "results_file: " << results_file << std::flush;
    AssertRelease(matfp != NULL, "Error creating the RESULTS.mat file.");

    if (dim == 1 or dim == 2)
    {
      std::vector<unsigned int> dims =
      { assem_per_dim[0], assem_per_dim[1]};

      const unsigned int n_cells_out = assem_per_dim[1] * assem_per_dim[0];
      Vector<double> phi(n_cells_out);
      Vector<double> dphi_real(n_cells_out);
      Vector<double> dphi_imag(n_cells_out);

      write_matlab_number(matfp, "keff", keff);

      for (unsigned int g = 0; g < n_groups; g++)
      {
        unsigned int idx = 0;
        unsigned int idx2 = 0;
        for (unsigned int k = 0; k < assem_per_dim[2]; k++)
        for (unsigned int j = 0; j < assem_per_dim[1]; j++)
        for (unsigned int i = 0; i < assem_per_dim[0]; i++)
        {

          if (materials.exist(i, j, k))
          {
            phi[idx2] = phi_per_assembly[g][idx];
            dphi_real[idx2] = noise_per_assembly[g][idx].real();
            dphi_imag[idx2] = noise_per_assembly[g][idx].imag();

            idx++;
          }
          else // No material
          {
            dphi_real[idx2] = std::nan("");
            dphi_imag[idx2] = std::nan("");
            phi[idx2] = std::nan("");
          }
          idx2++;
        }
        write_matlab_matrix(matfp, "FLX" + num_to_str(g + 1), dims, phi);
        write_matlab_matrix(matfp, "dFLX" + num_to_str(g + 1), dims,
            dphi_real, dphi_imag);
      }
    }
    else // dim == 3
    {
      std::vector<unsigned int> dims =
      { assem_per_dim[0], assem_per_dim[1],
        assem_per_dim[2]};

      const unsigned int n_cells_out = assem_per_dim[2] * assem_per_dim[1]
      * assem_per_dim[0];
      Vector<double> phi(n_cells_out);
      Vector<double> dphi_real(n_cells_out);
      Vector<double> dphi_imag(n_cells_out);

      write_matlab_number(matfp, "keff", keff);

      for (unsigned int g = 0; g < n_groups; g++)
      {
        unsigned int idx = 0;
        unsigned int idx2 = 0;
        for (unsigned int k = 0; k < assem_per_dim[2]; k++)
        for (unsigned int j = 0; j < assem_per_dim[1]; j++)
        for (unsigned int i = 0; i < assem_per_dim[0]; i++)
        {

          if (materials.exist(i,j,k))
          {
            phi[idx2] = phi_per_assembly[g][idx];
            dphi_real[idx2] = noise_per_assembly[g][idx].real();
            dphi_imag[idx2] = noise_per_assembly[g][idx].imag();

            idx++;
          }
          else // No material
          {
            dphi_real[idx2] = std::nan("");
            dphi_imag[idx2] = std::nan("");
            phi[idx2] = std::nan("");
          }
          idx2++;
        }
        write_matlab_matrix(matfp, "FLX" + num_to_str(g + 1), dims, phi);
        write_matlab_matrix(matfp, "dFLX" + num_to_str(g + 1), dims,
            dphi_real, dphi_imag);
      }
    }

    Mat_Close(matfp);
  }
#endif

/**
 * @brief Function that creates the output files .vtk.
 * It should be extended to create the power.
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::output_results_vtk () const
  {

    // Plot Magnitude and Phase of the noise
    DataOut<dim, DoFHandler<dim> > data_out;
    data_out.attach_dof_handler(dof_handler_block);
    std::string filename_vtk = out_file + ".vtk";

    // Get Magnitude and the phase of the neutron noise
    std::vector<std::vector<double>> dphi(n_groups, std::vector<double>(2));
    std::vector<Vector<double>> delta_phi_mag(n_groups);
    std::vector<Vector<double>> delta_phi_arg(n_groups);
    std::vector<Vector<double>> phi(n_groups);
    std::vector<Vector<double> > rel_noise(n_groups);

    for (unsigned g = 0; g < n_groups; g++)
    {
      delta_phi_mag[g].reinit(n_dofs);
      delta_phi_arg[g].reinit(n_dofs);
      phi[g].reinit(n_dofs);
      rel_noise[g].reinit(n_dofs);
    }

    for (unsigned g = 0; g < n_groups; g++)
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        dphi[g][0] = 0.0;
        dphi[g][1] = 0.0;
        phi[g][i] = 0.0;

        for (unsigned int m = 0; m < n_moments; m++)
        {
          dphi[g][0] += u_to_phi_coeff[0][m]
                        * delta_u(
                          g * n_moments * n_dofs * 2 + 2 * m * n_dofs
                          + i);
          dphi[g][1] += u_to_phi_coeff[0][m]
                        * delta_u(
                          g * n_moments * n_dofs * 2 + 2 * m * n_dofs
                          + n_dofs
                          + i);
          phi[g][i] += u_to_phi_coeff[0][m]
                       * u_crit(
                         g * n_moments * n_dofs * 2 + 2 * m * n_dofs
                         + i);
        }

        delta_phi_mag[g][i] = sqrt(
          dphi[g][0] * dphi[g][0] + dphi[g][1] * dphi[g][1]);
        delta_phi_arg[g][i] = atan2(dphi[g][1], dphi[g][0]) * 180 / M_PI;
        rel_noise[g][i] = delta_phi_mag[g][i] / phi[g][i] * 100;

      }
      data_out.add_data_vector(phi[g], "Static_Flux_g" + num_to_str(g + 1));
      data_out.add_data_vector(delta_phi_mag[g],
        "Noise_g" + num_to_str(g + 1) + "_Magnitude");
      data_out.add_data_vector(delta_phi_arg[g],
        "Noise_g" + num_to_str(g + 1) + "_Phase");
      data_out.add_data_vector(rel_noise[g],
        "RelNoise_g" + num_to_str(g + 1));
    }

    std::ofstream output(filename_vtk.c_str());
    data_out.build_patches(n_out_ref);
    data_out.write_vtk(output);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::output_detectors () const
  {
    // Calculate average Radial Distribution of noise
    // TODO Holes not allowed
    AssertRelease(dim == 3, "output_detectors Only allowed for dim==3");
    const unsigned int n_assemblies_per_plane = n_assemblies / assem_per_dim[2];
    const unsigned int n_axial_dtc = detector_levels.size() / 2;
    std::vector<std::vector<complex> > noise_at_detectors(n_groups,
      std::vector<complex>(n_detectors, complex(0.0, 0.0)));
    std::vector<double> norm(n_detectors, 0.0);

    for (unsigned int gr = 0; gr < n_groups; gr++)
    {
      noise_at_detectors[gr].resize(n_detectors);
      for (unsigned int a = 0; a < n_assemblies; a++)
      {
// Radial Position OK
        for (unsigned int dt_rad = 0; dt_rad < detectors_idx.size();
            dt_rad++)
        {
          if (detectors_idx[dt_rad] == (a % n_assemblies_per_plane))
          {
            // Axial Position OK
            const unsigned int a_axl = (a / n_assemblies_per_plane);
            for (unsigned int dt_axl = 0; dt_axl < n_axial_dtc;
                dt_axl++)
            {
              if ((a_axl >= detector_levels[2 * dt_axl])
                  and (a_axl <= detector_levels[2 * dt_axl + 1]))
              {
                const unsigned int dtc = dt_rad * n_axial_dtc
                                         + dt_axl;
                noise_at_detectors[gr][dtc] +=
                                               noise_per_assembly[gr][a]
                                               * assembly_pitch[2][dt_axl];
                norm[dtc] += assembly_pitch[2][dt_axl];
              }
            }
          }
        }
      }
    }
    // Normalize
    for (unsigned int gr = 0; gr < n_groups; gr++)
      for (unsigned int dtc = 0; dtc < n_detectors; dtc++)
      {
        AssertRelease(norm[dtc] != 0.0,
          num_to_str(dtc) + " detector not defined well");
        noise_at_detectors[gr][dtc] /= norm[dtc];
      }

    // Print
    std::ofstream out_stream(detectors_file.c_str(), std::ios::out);
    out_stream
    << "# det_num  det_radial  det_axial  noise_ieal_part  noise_imag_part\n";
    const unsigned int dtc_group = 1;    // Detectors only get thermal group
    for (unsigned int dtc = 0; dtc < n_detectors; dtc++)
    {
      out_stream << dtc + 1 << "   " << (dtc / n_axial_dtc) + 1 << "  "
                 << (dtc % n_axial_dtc) + 1
                 << "  "
                 << noise_at_detectors[dtc_group][dtc].real()
                 << "   "
                 << noise_at_detectors[dtc_group][dtc].imag()
                 << "\n";
    }
    out_stream.close();
  }

/**
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and the time-line.
 */
template <int dim, int n_fe_degree>
  void NoiseSPN<dim, n_fe_degree>::run ()
  {
    cout << "   Noise input files read." << std::endl;
    // Get Maximum memory
    PetscLogDouble memory;
    PetscMemorySetGetMaximumUsage();

    verbose_cout << "  assemble_system..." << std::flush;
    assemble_system();
    verbose_cout << "  Done! " << std::endl;

    cout << "      Number of DoFs per block: " << n_dofs << std::endl;
    cout << "      Number of Total DoFs: " << n_blocks * n_dofs << std::endl;
    cout << "   Matrices assembled. " << " Time = " << timer.cpu_time() << " s."
         << std::endl;

    verbose_cout << "  solve... " << std::flush;
    solve();
    verbose_cout << "  Done! " << std::endl;

    cout << "   Complex system solved. " << " Time = " << timer.cpu_time()
         << " s."
         << std::endl;
    PetscMemoryGetMaximumUsage(&memory);
    cout << "      Max Memory " << memory * 1e-6 << " MB." << std::endl;
    cout << "      Memory consumption of matrix elements "
         << memory_consumption * 1e-6
         << " MB" << std::endl;

    verbose_cout << "  postprocess... " << std::flush;
    postprocess();
    verbose_cout << "  Done! " << std::endl;

    if (results_file != "")
    {
#ifdef MATIO
      verbose_cout << "   output_results_matlab..." << std::flush;
      output_results_matlab();
      verbose_cout << "  Done!" << std::endl;
#else
      cout
      << "   WARNING! Output results in matlab only valid if matio is installed."
      << std::endl;
#endif
    }

    if (out_flag)
    {
      verbose_cout << "  output_results..." << std::flush;
      output_results_vtk();
      verbose_cout << "  Done!" << std::endl;
    }
    if (detectors_file != "")
    {
      verbose_cout << "  output_detectors..." << std::flush;
      output_detectors();
      verbose_cout << "  Done!" << std::endl;

    }
  }

template class NoiseSPN<1, 1> ;
template class NoiseSPN<1, 2> ;
template class NoiseSPN<1, 3> ;
template class NoiseSPN<1, 4> ;
template class NoiseSPN<1, 5> ;

template class NoiseSPN<2, 1> ;
template class NoiseSPN<2, 2> ;
template class NoiseSPN<2, 3> ;
template class NoiseSPN<2, 4> ;
template class NoiseSPN<2, 5> ;

template class NoiseSPN<3, 1> ;
template class NoiseSPN<3, 2> ;
template class NoiseSPN<3, 3> ;
template class NoiseSPN<3, 4> ;
template class NoiseSPN<3, 5> ;

