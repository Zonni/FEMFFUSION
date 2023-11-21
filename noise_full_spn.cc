/**
 * @file   noise_full_spn.cc
 * @brief Implementation of the Full SPN neutron noise.
 */

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <complex>
#include <cmath>
#include <math.h>  //  For M_PI and atan2
#include <string>

#include "static_diffusion.h"
#include "femffusion.h"
#include "matlab_io.h"
#include "utils.h"
#include "printing.h"
#include "noise_full_spn.h"
#include "static_full_spn.h"

using namespace dealii;
//using namespace std::complex_literals;
typedef std::complex<double> complex;

/**
 * @brief Constructor of the main class Noise.
 *
 */
template <int dim, int n_fe_degree>
  NoiseFullSPN<dim, n_fe_degree>::NoiseFullSPN (
    ParameterHandler &prm,
    StaticFullSPN<dim, n_fe_degree> &static_problem,
    const bool verbose,
    const bool silent)
  :
      verbose_cout(std::cout, verbose),
      cout(std::cout, !silent),
      n_groups(static_problem.n_groups),
      n_moments(static_problem.n_moments),
      n_even_moments(n_moments / 2),
      n_blocks(2 * static_problem.n_blocks),
      n_blocks_per_gr(static_problem.n_components),
      n_assemblies(static_problem.n_assemblies),
      n_dofs(static_problem.n_dofs),
      assem_per_dim(static_problem.assem_per_dim),
      assembly_pitch(static_problem.assembly_pitch),
      // boundary_conditions(static_problem.boundary_conditions),
      materials(static_problem.materials),
      tria(static_problem.tria),
      fe_block(static_problem.fe),
      dof_handler_block(static_problem.dof_handler),
      //constraints_block(static_problem.constraints),
      fe_system(static_problem.fe, n_blocks), // One for the real part and the other for the imaginary
      dof_handler(tria),
      //dof_handler_block(tria),
      boundary_conditions(static_problem.boundary_conditions),
      albedo_factors(static_problem.albedo_factors),
      phi_per_assembly(static_problem.phi_per_assembly[0]),
      pert(materials, dim, verbose),
      timer(static_problem.timer)
  {
    timer.restart();
    cout << std::endl;
    cout << "   ---------------------" << std::endl;
    cout << "     NOISE CALCULATION  " << std::endl;
    cout << "   ---------------------" << std::endl << std::endl;

    // Copy some values From Static
    n_out_ref = static_problem.n_out_ref;
    out_file = static_problem.out_file;
    out_flag = static_problem.out_flag;
    tol_ksp_noise = 1e-12;

    // Copy phi critic and normalize it
    phi_crit.reinit(PETSC_COMM_WORLD, // @suppress("Invalid arguments")
      2 * n_dofs * n_groups * n_blocks_per_gr,
      2 * n_dofs * n_groups * n_blocks_per_gr); // This vector is complex

    // Make reactor critical
    keff = static_problem.eigenvalues[0];
    //    static_problem.materials.make_critical(keff);
    //    std::ofstream out(out_file.c_str(), std::ios::out);
    //    out.close();
    //    static_problem.postprocess(); // Normalize it

    // Initialize
    for (unsigned int i = 0; i < static_problem.phi_sol[0].size(); i++)
    {
      phi_crit[i] = static_problem.phi_sol[0][i];
    }

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

      omega = 2 * M_PI * freq;
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
 * @brief
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::setup_system ()
  {
    dof_handler.distribute_dofs(fe_system);
    DoFRenumbering::component_wise(dof_handler);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sp.copy_from(dsp);
    A.reinit(sp);

    system_rhs.reinit(PETSC_COMM_WORLD, dof_handler.n_dofs(), dof_handler.n_dofs());
    delta_phi.reinit(PETSC_COMM_WORLD, dof_handler.n_dofs(), dof_handler.n_dofs());

  }

/**
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::get_parameters_from_command_line ()
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
      omega = 2 * M_PI * freq;

    // String
    get_string_from_options("-dtc_file", detectors_file);
    //get_string_from_options("-refinement_model", refinement_model);
    // lower_case(refinement_model);
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::assemble_system ()
  {

    if (pert.pert_type == "Cell_Wise")
    {
      verbose_cout << "assemble_rhs_cell_wise();..." << std::flush;
      assemble_rhs_cell_wise();
      verbose_cout << "  Done! " << std::endl;
    }
    else if (pert.pert_type == "Borders")
    {
      verbose_cout << "  assemble_system_borders_2g..." << std::flush;
      assemble_rhs_borders();
      verbose_cout << "  Done! " << std::endl;
    }
    else if (pert.pert_type == "BordersHex")
    {
      verbose_cout << "  assemble_system_bordershex_2g..." << std::flush;
      assemble_rhs_bordershex();
      verbose_cout << "  Done! " << std::endl;
    }
    else
      AssertRelease(false,
        "Not valid pert_type, only valid:  Cell_Wise | Borders | BordersHex");

    verbose_cout << "  assemble_A..." << std::flush;
    assemble_A();
    verbose_cout << "  Done! " << std::endl;

    memory_consumption = A.memory_consumption();

    // Print in Matlab
    std::string print_matrices_matlab;
    get_string_from_options("-print_system_noise_matlab", print_matrices_matlab);
    if (!print_matrices_matlab.empty())
    {
      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);
      print_matrix_in_matlab(A, "A", out, 8);
      print_vector_in_matlab(system_rhs, "b", out, 8);
    }
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void
  NoiseFullSPN<dim, n_fe_degree>::block_to_moment_group (unsigned int block,
    unsigned int &moment,
    unsigned int &group,
    unsigned int &complex)
  {

    complex = block / (n_groups * n_blocks_per_gr);
    unsigned int block_real = block % (n_groups * n_blocks_per_gr);

    if (block_real < n_groups)
    {
      moment = 0;
      group = block_real;
    }
    else if (block_real < n_groups * (1 + dim))
    {
      moment = 1;
      group = (block_real - n_groups) / dim;
    }
    else if (block_real < (dim + 1) * n_groups + n_groups)
    {
      moment = 2;
      group = (block_real - n_groups * (1 + dim));
    }
    else
    {
      moment = 3;
      group = (block_real - (dim + 2) * n_groups) / dim;
    }

  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::assemble_A ()
  {
    double factor = 0.0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    const unsigned int dofs_per_cell = fe_system.dofs_per_cell;

    FEValues<dim> fe_values(fe_system, quadrature_formula,
      update_values | update_gradients | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe_system, face_quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values
      | update_normal_vectors);

    FullMatrix<double> cell_matrix_A(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    unsigned int cc_i, moment_i, group_i;
    unsigned int cc_j, moment_j, group_j;

    std::vector<std::vector<std::vector<std::vector<complex> > > > A_factor(n_groups,
      std::vector<std::vector<std::vector<complex>>>(n_groups,
        std::vector<std::vector<complex>>(n_moments,
          std::vector<complex>(n_moments))));

    verbose_cout << "omega " << omega << std::endl;

    // Indicates which components of a vector-valued finite element constitute a single scalar
    std::vector<std::vector<std::vector<FEValuesExtractors::Scalar>> > scalar_flux(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Scalar> >(n_groups,
        std::vector<FEValuesExtractors::Scalar>(2)));
    std::vector<std::vector<std::vector<FEValuesExtractors::Vector>> > current(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Vector> >(n_groups,
        std::vector<FEValuesExtractors::Vector>(2)));

    for (unsigned int mm = 0; mm < n_even_moments; mm++)
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int cc = 0; cc < 2; cc++)
        {
          scalar_flux[mm][g][cc] = FEValuesExtractors::Scalar(
            (mm * (dim + 1) * n_groups + g)
            + cc * n_blocks_per_gr * n_groups);
          current[mm][g][cc] = FEValuesExtractors::Vector(
            ((mm + 1) * n_groups + mm * dim * n_groups + g * dim) +
            cc * n_blocks_per_gr * n_groups);
        }

    std::vector<complex> prec_factor(n_groups);
    complex divisor;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {

      cell_matrix_A = 0;
      fe_values.reinit(cell);

      // Get the material coefficients:
      const unsigned int mat_id = materials.get_material_id<dim>(cell);

      for (unsigned int g = 0; g < n_groups; g++)
      {
        prec_factor[g] = (1 - materials.get_delayed_fraction_sum(mat_id)) *
                         materials.get_prompt_spectra(mat_id, g);

        for (unsigned int p = 0; p < materials.get_n_precursors(); p++)
        {
          divisor = complex(materials.get_delayed_decay_constant(mat_id, p), omega);
          prec_factor[g] += materials.get_delayed_decay_constant(mat_id, p)
                            * materials.get_delayed_fraction(mat_id, p)
                            / (divisor)
                            * materials.get_delayed_spectra(mat_id, p, g);
        }
      }

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          for (unsigned int m1 = 0; m1 < n_moments; m1++)
            for (unsigned int m2 = 0; m2 < n_moments; m2++)
            {
              if (m1 == 0 and m2 == 0)
              {
                if (g1 == g2)
                {
                  A_factor[g1][g2][m1][m2] = complex(0,
                                               omega / materials.get_velocitiy(mat_id,
                                                 g1))
                                             + materials.get_sigma_r(g1, mat_id);
                }
                else
                {
                  A_factor[g1][g2][m1][m2] = -materials.get_sigma_s(g2, g1, mat_id);

                }

                A_factor[g1][g2][m1][m2] += -prec_factor[g1]
                                            * materials.get_nu_sigma_f(g2, mat_id);
              }
              else if (m1 == m2)
              {
                if (g1 == g2) // Only isotropic scattering
                  A_factor[g1][g2][m1][m2] = complex(0,
                                               omega / materials.get_velocitiy(mat_id,
                                                 g1))
                                             + materials.get_sigma_t(g1, mat_id);

              }
            }

      //      std::cout << "cell " << cell << std::endl;
      //      std::cout << "A_factor_00 " << A_factor[0][0][0][0] << std::endl;
      //      std::cout << "A_factor_01 " << A_factor[0][0][0][1] << std::endl;
      //      std::cout << "A_factor_10 " << A_factor[0][0][1][0] << std::endl;
      //      std::cout << "A_factor_11 " << A_factor[0][0][1][1] << std::endl;
      //
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          const unsigned int block_i = fe_system.system_to_component_index(i).first;
          const unsigned int block_j = fe_system.system_to_component_index(j).first;

          block_to_moment_group(block_i, moment_i, group_i, cc_i);
          block_to_moment_group(block_j, moment_j, group_j, cc_j);

          if (cc_i == cc_j) // Real
          {
            if ((moment_i == 0 and moment_j == 0) or (moment_i == 2 and moment_j == 2))
            {
              const unsigned int moment_type = moment_i / 2;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_matrix_A(i, j) +=
                    A_factor[group_i][group_j][moment_i][moment_j].real()
                    * fe_values[scalar_flux[moment_type][group_i][cc_i]].value(i, q)
                    * fe_values[scalar_flux[moment_type][group_j][cc_j]].value(j, q)
                    * fe_values.JxW(q);

            }
            else if ((moment_i == 1 and moment_j == 1)
                     or (moment_i == 3 and moment_j == 3))
            {
              const unsigned int moment_type = moment_i / 2;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_matrix_A(i, j) +=
                    A_factor[group_i][group_j][moment_i][moment_j].real()
                    * fe_values[current[moment_type][group_i][cc_i]].value(i, q)
                    * fe_values[current[moment_type][group_j][cc_j]].value(j, q)
                    * fe_values.JxW(q);
            }
            else if (moment_j - moment_i == 1 or moment_i - moment_j == 1)
            {
              if (group_i == group_j)
              {
                if (moment_i % 2 == 0)
                {
                  const unsigned int moment_type_i = moment_i / 2;
                  const unsigned int moment_type_j = moment_j / 2;

                  for (unsigned int q = 0; q < n_q_points; ++q)
                    cell_matrix_A(i, j) +=
                        -grad_coeff_full_spn[moment_i][moment_j] // Debido a la discretzacion FEM
                        * fe_values[scalar_flux[moment_type_i][group_i][cc_i]].gradient(i,
                          q)
                        * fe_values[current[moment_type_j][group_j][cc_j]].value(j, q)
                        * fe_values.JxW(q);
                }
                else
                {
                  const unsigned int moment_type_i = moment_i / 2;
                  const unsigned int moment_type_j = moment_j / 2;

                  for (unsigned int q = 0; q < n_q_points; ++q)
                    cell_matrix_A(i, j) +=
                        grad_coeff_full_spn[moment_i][moment_j]
                        * fe_values[current[moment_type_i][group_i][cc_i]].value(i, q)
                        * fe_values[scalar_flux[moment_type_j][group_j][cc_j]].gradient(j,
                          q)
                        * fe_values.JxW(q);
                }
              }
            }
          }
          else // Imag
          {
            if ((moment_i == 0 and moment_j == 0) or (moment_i == 2 and moment_j == 2))
            {
              const unsigned int moment_type = moment_i / 2;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_matrix_A(i, j) += pow(-1, cc_j)
                    * A_factor[group_i][group_j][moment_i][moment_j].imag()
                    * fe_values[scalar_flux[moment_type][group_i][cc_i]].value(i, q)
                    * fe_values[scalar_flux[moment_type][group_j][cc_j]].value(j, q)
                    * fe_values.JxW(q);

            }
            else if ((moment_i == 1 and moment_j == 1)
                     or (moment_i == 3 and moment_j == 3))
            {
              const unsigned int moment_type = moment_i / 2;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_matrix_A(i, j) += pow(-1, cc_j)
                    * A_factor[group_i][group_j][moment_i][moment_j].imag()
                    * fe_values[current[moment_type][group_i][cc_i]].value(i, q)
                    * fe_values[current[moment_type][group_j][cc_j]].value(j, q)
                    * fe_values.JxW(q);
            }
          }

        }

      // Boundary
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
        {
          types::boundary_id boundary_id = cell->face(f)->boundary_id();
          AssertIndexRange(boundary_id, boundary_conditions.size());

          if (boundary_conditions[boundary_id] > 1) // Comprobar
          {
            fe_face_values.reinit(cell, f);
            switch (boundary_conditions[boundary_id])
            {
              case 2: // Vacuum BC
                factor = 1.0;
                break;
              default: // Custom Albedo BC
                AssertRelease(false, "BC Not Implemented");
                factor = 100.0;
                break;
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const unsigned int block_i =
                                             fe_system.system_to_component_index(i).first;
                const unsigned int block_j =
                                             fe_system.system_to_component_index(j).first;

                block_to_moment_group(block_i, moment_i, group_i, cc_i);
                block_to_moment_group(block_j, moment_j, group_j, cc_j);

                if (group_i == group_j and moment_i % 2 == 0 and moment_j % 2 == 0
                    and cc_i == cc_j)
                {
                  const unsigned int mt_i = moment_i / 2;
                  const unsigned int mt_j = moment_j / 2;

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                  {

                    cell_matrix_A(i, j) +=
                        factor * bound_vacuum_coeff[mt_i][mt_j]
                        * fe_face_values[scalar_flux[mt_i][group_i][cc_i]].value(i, q)
                        * fe_face_values[scalar_flux[mt_j][group_j][cc_j]].value(j, q)
                        * fe_face_values.JxW(q);
                  }
                }
              }
          }

        }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix_A, local_dof_indices, A);

    }

    A.compress(VectorOperation::add);
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::assemble_rhs_cell_wise ()
  {
    // Setup B
    PETScWrappers::SparseMatrix B;
    B.reinit(sp);

    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int dofs_per_cell = fe_system.dofs_per_cell;
    unsigned int cc_i, moment_i, group_i;
    unsigned int cc_j, moment_j, group_j;

    FEValues<dim> fe_values(fe_system, quadrature_formula,
      update_values | update_gradients | update_JxW_values);

    FullMatrix<double> cell_matrix_B(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<std::vector<std::vector<complex> > > B_factor(n_groups,
      std::vector<std::vector<complex>>(n_groups,
        std::vector<complex>(n_moments)));

    // Indicates which components of a vector-valued finite element constitute a single scalar
    std::vector<std::vector<std::vector<FEValuesExtractors::Scalar>> > scalar_flux(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Scalar> >(n_groups,
        std::vector<FEValuesExtractors::Scalar>(2)));
    std::vector<std::vector<std::vector<FEValuesExtractors::Vector>> > current(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Vector> >(n_groups,
        std::vector<FEValuesExtractors::Vector>(2)));

    std::vector<complex> prec_factor(n_groups);
    complex divisor;

    for (unsigned int mm = 0; mm < n_even_moments; mm++)
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int cc = 0; cc < 2; cc++)
        {
          scalar_flux[mm][g][cc] = FEValuesExtractors::Scalar(
            (mm * (dim + 1) * n_groups + g)
            + cc * n_blocks_per_gr * n_groups);
          current[mm][g][cc] = FEValuesExtractors::Vector(
            ((mm + 1) * n_groups + mm * dim * n_groups + g * dim) +
            cc * n_blocks_per_gr * n_groups);
        }

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {

      cell_matrix_B = 0;
      fe_values.reinit(cell);

      // Get the material coefficients:
      const unsigned int mat_id = materials.get_material_id<dim>(cell);

      for (unsigned int g = 0; g < n_groups; g++)
      {
        prec_factor[g] = (1 - materials.get_delayed_fraction_sum(mat_id)) *
                         materials.get_prompt_spectra(mat_id, g);

        for (unsigned int p = 0; p < materials.get_n_precursors(); p++)
        {
          divisor = complex(materials.get_delayed_decay_constant(mat_id, p), omega);
          prec_factor[g] += materials.get_delayed_decay_constant(mat_id, p)
                            * materials.get_delayed_fraction(mat_id, p)
                            / (divisor)
                            * materials.get_delayed_spectra(mat_id, p, g);
        }
      }

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          for (unsigned int m = 0; m < n_moments; m++)
          {
            if (m == 0)
            {
              if (g1 == g2)
              {
                B_factor[g1][g2][m] = -pert.get_delta_sigma_r(g1, mat_id);
              }
              else
              {
                B_factor[g1][g2][m] = pert.get_delta_sigma_s(g2, g1, mat_id);
              }

              B_factor[g1][g2][m] += prec_factor[g1]
                                     * pert.get_delta_sigma_f(g2, mat_id)
                                     / keff;
            }
            else
            {
              if (g1 == g2) // Only isotropic scattering
                B_factor[g1][g2][m] = -pert.get_delta_sigma_t(g1, mat_id);

            }
          }

      //
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {

          const unsigned int block_i = fe_system.system_to_component_index(i).first;
          const unsigned int block_j = fe_system.system_to_component_index(j).first;

          // Get moment, group and complex_component from the block number
          block_to_moment_group(block_i, moment_i, group_i, cc_i);
          block_to_moment_group(block_j, moment_j, group_j, cc_j);

          if (cc_i == cc_j) // Real
          {
            if ((moment_i == 0 and moment_j == 0) or (moment_i == 2 and moment_j == 2))
            {
              const unsigned int moment_type = moment_i / 2;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_matrix_B(i, j) += B_factor[group_i][group_j][moment_i].real()
                    * fe_values[scalar_flux[moment_type][group_i][cc_i]].value(i, q)
                    * fe_values[scalar_flux[moment_type][group_j][cc_j]].value(j, q)
                    * fe_values.JxW(q);

            }
            else if ((moment_i == 1 and moment_j == 1)
                     or (moment_i == 3 and moment_j == 3))
            {
              const unsigned int moment_type = moment_i / 2;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_matrix_B(i, j) += B_factor[group_i][group_j][moment_i].real()
                    * fe_values[current[moment_type][group_i][cc_i]].value(i, q)
                    * fe_values[current[moment_type][group_j][cc_j]].value(j, q)
                    * fe_values.JxW(q);
            }
          }
          else // Imag
          if ((moment_i == 0 and moment_j == 0) or (moment_i == 2 and moment_j == 2))
          {
            const unsigned int moment_type = moment_i / 2;

            for (unsigned int q = 0; q < n_q_points; ++q)
              cell_matrix_B(i, j) += pow(-1, cc_j)
                  * B_factor[group_i][group_j][moment_i].imag()
                  * fe_values[scalar_flux[moment_type][group_i][cc_i]].value(i, q)
                  * fe_values[scalar_flux[moment_type][group_j][cc_j]].value(j, q)
                  * fe_values.JxW(q);

          }
          else if ((moment_i == 1 and moment_j == 1)
                   or (moment_i == 3 and moment_j == 3))
          {
            const unsigned int moment_type = moment_i / 2;

            for (unsigned int q = 0; q < n_q_points; ++q)
              cell_matrix_B(i, j) += pow(-1, cc_j)
                  * B_factor[group_i][group_j][moment_i].imag()
                  * fe_values[current[moment_type][group_i][cc_i]].value(i, q)
                  * fe_values[current[moment_type][group_j][cc_j]].value(j, q)
                  * fe_values.JxW(q);
          }
        }

      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_matrix_B, local_dof_indices, B);
    }
    B.compress(VectorOperation::add);
    // Get system rhs
    B.vmult(system_rhs, phi_crit);

    ////////////////////////////////////////
    std::string print_matrices_matlab;
    get_string_from_options("-print_system_noise_matlab_B",
      print_matrices_matlab);

    if (!print_matrices_matlab.empty())
    {
      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);
      print_matrix_in_matlab(B, "B", out, 8);
    }
    B.clear();
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::assemble_rhs_borders ()
  {
    PETScWrappers::SparseMatrix B;
    B.reinit(sp);

    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    const unsigned int dofs_per_cell = fe_system.dofs_per_cell;

    FEFaceValues<dim> fe_face_values(fe_system, face_quadrature_formula,
      update_values | update_JxW_values);

    FullMatrix<double> cell_matrix_B(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // -------------------------------------------------------------------------- //
    unsigned int moment_i, group_i, cc_i;
    unsigned int moment_j, group_j, cc_j;

    std::vector<std::vector<std::vector<complex> > > B_factor(n_groups,
      std::vector<std::vector<complex>>(n_groups,
        std::vector<complex>(n_moments)));

    std::vector<complex> prec_factor(n_groups);
    complex divisor;

    // Indicates which components of a vector-valued finite element constitute a single scalar
    std::vector<std::vector<std::vector<FEValuesExtractors::Scalar>> > scalar_flux(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Scalar> >(n_groups,
        std::vector<FEValuesExtractors::Scalar>(2)));
    std::vector<std::vector<std::vector<FEValuesExtractors::Vector>> > current(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Vector> >(n_groups,
        std::vector<FEValuesExtractors::Vector>(2)));

    for (unsigned int mm = 0; mm < n_even_moments; mm++)
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int cc = 0; cc < 2; cc++)
        {
          scalar_flux[mm][g][cc] = FEValuesExtractors::Scalar(
            (mm * (dim + 1) * n_groups + g)
            + cc * n_blocks_per_gr * n_groups);
          current[mm][g][cc] = FEValuesExtractors::Vector(
            ((mm + 1) * n_groups + mm * dim * n_groups + g * dim) +
            cc * n_blocks_per_gr * n_groups);
        }

    // -------------------------------------------------------------------------- //
    for (auto cell : dof_handler.cell_iterators_on_level(0))
    {
      unsigned int n_active_children = cell->number_of_children();
      auto active_cells = GridTools::get_active_child_cells<DoFHandler<dim> >(cell);
      if (active_cells.size() == 0) /* get_active_child_cells only selects children */
      {
        active_cells.push_back(cell);
      }

      // -------------------------------------------------------------------------- //
      for (unsigned int ch = 0; ch < n_active_children; ch++)
      {
        //  if (cell_child->active())
        typename DoFHandler<dim>::cell_iterator cell_child_active = active_cells[ch];
        const unsigned int mat_id = materials.get_material_id<dim>(cell_child_active);

        for (unsigned int face = 0;
            face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
          // We could have used child_cell_on_face()
          const unsigned int face_pert_id = pert.get_pertubation_face_id(
            mat_id, face, ch, n_active_children);
          if (face_pert_id != static_cast<unsigned int>(-1))
          {
            fe_face_values.reinit(cell_child_active, face);

            for (unsigned int g = 0; g < n_groups; g++)
            {
              prec_factor[g] = (1 - materials.get_delayed_fraction_sum(mat_id)) *
                               materials.get_prompt_spectra(mat_id, g);

              for (unsigned int p = 0; p < materials.get_n_precursors(); p++)
              {
                divisor = complex(materials.get_delayed_decay_constant(mat_id, p), omega);
                prec_factor[g] += materials.get_delayed_decay_constant(mat_id, p)
                                  * materials.get_delayed_fraction(mat_id, p)
                                  / (divisor)
                                  * materials.get_delayed_spectra(mat_id, p, g);
              }
            }

            cell_matrix_B = 0;

            // Get the material coefficients:
            const unsigned int mat_id = materials.get_material_id<dim>(cell);
            for (unsigned int g1 = 0; g1 < n_groups; g1++)
              for (unsigned int g2 = 0; g2 < n_groups; g2++)
                for (unsigned int m = 0; m < n_moments; m++)
                {
                  if (m == 0)
                  {
                    if (g1 == g2)
                    {
                      B_factor[g1][g2][m] = -pert.get_delta_sigma_r(g1, mat_id);
                    }
                    else
                    {
                      B_factor[g1][g2][m] = pert.get_delta_sigma_s(g2, g1, mat_id);
                    }

                    B_factor[g1][g2][m] += prec_factor[g1]
                                           * pert.get_delta_sigma_f(g2, mat_id)
                                           / keff;
                  }
                  else
                  {
                    if (g1 == g2) // Only isotropic scattering
                      B_factor[g1][g2][m] = -pert.get_delta_sigma_t(g1, mat_id);

                  }
                }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {

                const unsigned int block_i = fe_system.system_to_component_index(i).first;
                const unsigned int block_j = fe_system.system_to_component_index(j).first;

                // Get moment, group and complex_component from the block number
                block_to_moment_group(block_i, moment_i, group_i, cc_i);
                block_to_moment_group(block_j, moment_j, group_j, cc_j);

                if (cc_i == cc_j) // Real
                {
                  if ((moment_i == 0 and moment_j == 0) or (moment_i == 2
                                                            and moment_j == 2))
                  {
                    const unsigned int moment_type = moment_i / 2;

                    for (unsigned int q = 0; q < n_face_q_points; ++q)
                      cell_matrix_B(i, j) += B_factor[group_i][group_j][moment_i].real()
                          * fe_face_values[scalar_flux[moment_type][group_i][cc_i]].value(
                            i, q)
                          * fe_face_values[scalar_flux[moment_type][group_j][cc_j]].value(
                            j, q)
                          * fe_face_values.JxW(q);

                  }
                  else if ((moment_i == 1 and moment_j == 1)
                           or (moment_i == 3 and moment_j == 3))
                  {
                    const unsigned int moment_type = moment_i / 2;

                    for (unsigned int q = 0; q < n_face_q_points; ++q)
                      cell_matrix_B(i, j) += B_factor[group_i][group_j][moment_i].real()
                          * fe_face_values[current[moment_type][group_i][cc_i]].value(
                            i, q)
                          * fe_face_values[current[moment_type][group_j][cc_j]].value(
                            j, q)
                          * fe_face_values.JxW(q);
                  }
                }
                else // Imag
                if ((moment_i == 0 and moment_j == 0) or (moment_i == 2 and moment_j == 2))
                {
                  const unsigned int moment_type = moment_i / 2;

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                    cell_matrix_B(i, j) += pow(-1, cc_j)
                        * B_factor[group_i][group_j][moment_i].imag()
                        * fe_face_values[scalar_flux[moment_type][group_i][cc_i]].value(
                          i, q)
                        * fe_face_values[scalar_flux[moment_type][group_j][cc_j]].value(
                          j, q)
                        * fe_face_values.JxW(q);

                }
                else if ((moment_i == 1 and moment_j == 1)
                         or (moment_i == 3 and moment_j == 3))
                {
                  const unsigned int moment_type = moment_i / 2;

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                    cell_matrix_B(i, j) += pow(-1, cc_j)
                        * B_factor[group_i][group_j][moment_i].imag()
                        * fe_face_values[current[moment_type][group_i][cc_i]].value(i, q)
                        * fe_face_values[current[moment_type][group_j][cc_j]].value(j, q)
                        * fe_face_values.JxW(q);
                }
              }

            cell_child_active->get_dof_indices(local_dof_indices);

            constraints.distribute_local_to_global(cell_matrix_B,
              local_dof_indices, B);

          } // End IF there is perturbation
        } // End FOR face
      } // End FOR active children
    } // End FOR cell at level(0)

    // Get system rhs
    B.compress(VectorOperation::add);
    B.vmult(system_rhs, phi_crit);

    ////////////////////////////////////////
    std::string print_matrices_matlab;
    get_string_from_options("-print_system_noise_matlab_B",
      print_matrices_matlab);

    if (!print_matrices_matlab.empty())
    {
      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);
      print_matrix_in_matlab(B, "B", out, 8);
    }
    B.clear();
  }

/**
 * @brief Assemble the matrices or prepare structure in the matrix-free cases.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::assemble_rhs_bordershex ()
  {
    PETScWrappers::SparseMatrix B;
    B.reinit(sp);

    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    const unsigned int dofs_per_cell = fe_system.dofs_per_cell;
    std::vector<unsigned int> quad_in_hex(n_assemblies, 0);

    FEFaceValues<dim> fe_face_values(fe_system, face_quadrature_formula,
      update_values | update_JxW_values);

    FullMatrix<double> cell_matrix_B(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // -------------------------------------------------------------------------- //
    unsigned int moment_i, group_i, cc_i;
    unsigned int moment_j, group_j, cc_j;

    std::vector<std::vector<std::vector<complex> > > B_factor(n_groups,
      std::vector<std::vector<complex>>(n_groups,
        std::vector<complex>(n_moments)));

    // Indicates which components of a vector-valued finite element constitute a single scalar
    std::vector<std::vector<std::vector<FEValuesExtractors::Scalar>> > scalar_flux(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Scalar> >(n_groups,
        std::vector<FEValuesExtractors::Scalar>(2)));
    std::vector<std::vector<std::vector<FEValuesExtractors::Vector>> > current(
      n_even_moments,
      std::vector<std::vector<FEValuesExtractors::Vector> >(n_groups,
        std::vector<FEValuesExtractors::Vector>(2)));

    std::vector<complex> prec_factor(n_groups);
    complex divisor;

    for (unsigned int mm = 0; mm < n_even_moments; mm++)
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int cc = 0; cc < 2; cc++)
        {
          scalar_flux[mm][g][cc] = FEValuesExtractors::Scalar(
            (mm * (dim + 1) * n_groups + g)
            + cc * n_blocks_per_gr * n_groups);
          current[mm][g][cc] = FEValuesExtractors::Vector(
            ((mm + 1) * n_groups + mm * dim * n_groups + g * dim) +
            cc * n_blocks_per_gr * n_groups);
        }

    // Integrate Interiors Cells
    // -------------------------------------------------------------------------- //
    for (auto cell : dof_handler.cell_iterators_on_level(0))
    {
      unsigned int n_active_children = cell->number_of_children();
      auto active_cells = GridTools::get_active_child_cells<DoFHandler<dim> >(
        cell);
      if (active_cells.size() == 0) /* get_active_child_cells only selects children */
      {
        active_cells.push_back(cell);
      }

      // -------------------------------------------------------------------------- //
      for (unsigned int ch = 0; ch < n_active_children; ch++)
      {
        typename DoFHandler<dim>::cell_iterator cell_child_active = active_cells[ch];
        const unsigned quad_id = quad_in_hex[cell_child_active->user_index()]
                                 / n_active_children;
        quad_in_hex[cell_child_active->user_index()]++;
        const unsigned int mat_id = materials.get_material_id<dim>(cell_child_active);
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
          // We could have used child_cell_on_face()
          const unsigned int face_pert_id = pert.get_pertubation_face_id_hex(mat_id,
            quad_id, face, ch, n_active_children);

          if (face_pert_id != static_cast<unsigned int>(-1))
          {
            fe_face_values.reinit(cell_child_active, face);
            cell_matrix_B = 0;

            for (unsigned int g = 0; g < n_groups; g++)
            {
              prec_factor[g] = (1 - materials.get_delayed_fraction_sum(mat_id)) *
                               materials.get_prompt_spectra(mat_id, g);

              for (unsigned int p = 0; p < materials.get_n_precursors(); p++)
              {
                divisor = complex(materials.get_delayed_decay_constant(mat_id, p), omega);
                prec_factor[g] += materials.get_delayed_decay_constant(mat_id, p)
                                  * materials.get_delayed_fraction(mat_id, p)
                                  / (divisor)
                                  * materials.get_delayed_spectra(mat_id, p, g);
              }
            }

            // Get the material coefficients:
            const unsigned int mat_id = materials.get_material_id<dim>(cell);
            for (unsigned int g1 = 0; g1 < n_groups; g1++)
              for (unsigned int g2 = 0; g2 < n_groups; g2++)
                for (unsigned int m = 0; m < n_moments; m++)
                {
                  if (m == 0)
                  {
                    if (g1 == g2)
                    {
                      B_factor[g1][g2][m] = -pert.get_delta_sigma_r(g1, mat_id);
                    }
                    else
                    {
                      B_factor[g1][g2][m] = pert.get_delta_sigma_s(g2, g1, mat_id);
                    }

                    B_factor[g1][g2][m] += prec_factor[g1]
                                           * pert.get_delta_sigma_f(g2, mat_id)
                                           / keff;
                  }
                  else
                  {
                    if (g1 == g2) // Only isotropic scattering
                      B_factor[g1][g2][m] = -pert.get_delta_sigma_t(g1, mat_id);

                  }
                }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {

                const unsigned int block_i = fe_system.system_to_component_index(i).first;
                const unsigned int block_j = fe_system.system_to_component_index(j).first;

                // Get moment, group and complex_component from the block number
                block_to_moment_group(block_i, moment_i, group_i, cc_i);
                block_to_moment_group(block_j, moment_j, group_j, cc_j);

                if (cc_i == cc_j) // Real
                {
                  if ((moment_i == 0 and moment_j == 0) or (moment_i == 2
                                                            and moment_j == 2))
                  {
                    const unsigned int moment_type = moment_i / 2;

                    for (unsigned int q = 0; q < n_face_q_points; ++q)
                      cell_matrix_B(i, j) += B_factor[group_i][group_j][moment_i].real()
                          * fe_face_values[scalar_flux[moment_type][group_i][cc_i]].value(
                            i, q)
                          * fe_face_values[scalar_flux[moment_type][group_j][cc_j]].value(
                            j, q)
                          * fe_face_values.JxW(q);

                  }
                  else if ((moment_i == 1 and moment_j == 1)
                           or (moment_i == 3 and moment_j == 3))
                  {
                    const unsigned int moment_type = moment_i / 2;

                    for (unsigned int q = 0; q < n_face_q_points; ++q)
                      cell_matrix_B(i, j) += B_factor[group_i][group_j][moment_i].real()
                          * fe_face_values[current[moment_type][group_i][cc_i]].value(
                            i, q)
                          * fe_face_values[current[moment_type][group_j][cc_j]].value(
                            j, q)
                          * fe_face_values.JxW(q);
                  }
                }
                else // Imag
                if ((moment_i == 0 and moment_j == 0) or (moment_i == 2 and moment_j == 2))
                {
                  const unsigned int moment_type = moment_i / 2;

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                    cell_matrix_B(i, j) += pow(-1, cc_j)
                        * B_factor[group_i][group_j][moment_i].imag()
                        * fe_face_values[scalar_flux[moment_type][group_i][cc_i]].value(
                          i, q)
                        * fe_face_values[scalar_flux[moment_type][group_j][cc_j]].value(
                          j, q)
                        * fe_face_values.JxW(q);

                }
                else if ((moment_i == 1 and moment_j == 1)
                         or (moment_i == 3 and moment_j == 3))
                {
                  const unsigned int moment_type = moment_i / 2;

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                    cell_matrix_B(i, j) += pow(-1, cc_j)
                        * B_factor[group_i][group_j][moment_i].imag()
                        * fe_face_values[current[moment_type][group_i][cc_i]].value(i, q)
                        * fe_face_values[current[moment_type][group_j][cc_j]].value(j, q)
                        * fe_face_values.JxW(q);
                }
              }

            cell_child_active->get_dof_indices(local_dof_indices);

            constraints.distribute_local_to_global(cell_matrix_B,
              local_dof_indices, B);

          } // End IF there is perturbation
        } // End FOR face
      } // End FOR active children
    } // End FOR cell at level(0)

    // Get system rhs
    B.compress(VectorOperation::add);
    B.vmult(system_rhs, phi_crit);

    ////////////////////////////////////////
    std::string print_matrices_matlab;
    get_string_from_options("-print_system_noise_matlab_B",
      print_matrices_matlab);

    if (!print_matrices_matlab.empty())
    {
      std::ofstream out(print_matrices_matlab.c_str(), std::ios::out);
      print_matrix_in_matlab(B, "B", out, 8);
    }
    B.clear();
  }

/**
 * @brief Solve the eigenvalue problem with the power iteration method.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::solve ()
  {
    const unsigned int max_its = 1e6;
    cout << "      system_rhs norm: " << system_rhs.l2_norm() << " " << std::endl; // @suppress("Invalid overload")

    PC pc;
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp); // @suppress("Function cannot be resolved") // @suppress("Invalid arguments")
    KSPGetPC(ksp, &pc);
    KSPSetTolerances(ksp, tol_ksp_noise, 0.0, PETSC_DEFAULT, max_its);
    KSPSetType(ksp, KSPBICG);
    // KSPSetType(ksp, KSPFGMRES);
    PCSetType(pc, PCILU);
    PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    KSPSolve(ksp, system_rhs, delta_phi);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    AssertRelease(reason > 0,
      "Noise KSP not converged, reason: " + num_to_str(static_cast<int>(reason)));

    cout << "      delta_phi norm: " << delta_phi.l2_norm() << std::endl;
    A.clear();
  }

/**
 * @brief Get the scalar neutron noise
 *
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::postprocess ()
  {
    // Erase the content of output file
    std::ofstream out(out_file.c_str(), std::ios::app);

    // Initialize all that  is needed to iterate over dofs and cells
    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(fe_system, quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values);

    unsigned int n_q_points = quadrature_formula.size();
    unsigned int n_cells_out = n_assemblies;

    std::vector<std::vector<double>> delta_cell_phi(n_groups, std::vector<double>(2));

    // Initialize and resize the vectors where it is stored the solution
    noise_per_assembly.resize(n_groups, std::vector<complex>(n_cells_out, 0.0));

    std::vector<double> volume_per_assembly(n_cells_out, 0.0);

    std::vector<std::vector<std::vector<double> > > local_delta_phi(
      n_groups,
      std::vector<std::vector<double> >(2,
        std::vector<double>(n_q_points)));

    std::vector<std::vector<FEValuesExtractors::Scalar>> delta_scalar_flux(n_groups,
      std::vector<FEValuesExtractors::Scalar>(2));

    for (unsigned int g = 0; g < n_groups; g++)
      for (unsigned int cc = 0; cc < 2; cc++)
      {
        delta_scalar_flux[g][cc] = FEValuesExtractors::Scalar(
          g + cc * n_blocks_per_gr * n_groups);
      }

    // Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      for (unsigned int g = 0; g < n_groups; g++)
      {
        delta_cell_phi[g][0] = 0.0;
        delta_cell_phi[g][1] = 0.0;
      }

      for (unsigned int g = 0; g < n_groups; g++)
      {
        fe_values[delta_scalar_flux[g][0]].get_function_values(delta_phi,
          local_delta_phi[g][0]);
        fe_values[delta_scalar_flux[g][1]].get_function_values(delta_phi,
          local_delta_phi[g][1]);
      }

      for (unsigned int q = 0; q < n_q_points; q++)
      {
        for (unsigned int g = 0; g < n_groups; g++)
        {
          delta_cell_phi[g][0] += local_delta_phi[g][0][q]
                                  * fe_values.JxW(q);
          delta_cell_phi[g][1] += local_delta_phi[g][1][q]
                                  * fe_values.JxW(q);
        }
      }

      for (unsigned int g = 0; g < n_groups; g++)
        noise_per_assembly[g][cell->user_index()] += complex(delta_cell_phi[g][0],
          delta_cell_phi[g][1]);

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
  }

#ifdef MATIO
/**
 * @brief Function that creates the output files .vtk.
 * It should be extended to create the power.
 */
template <int dim, int n_fe_degree>
void
NoiseFullSPN<dim, n_fe_degree>::output_results_matlab () const
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
  void NoiseFullSPN<dim, n_fe_degree>::output_results_vtk () const
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
    std::vector<Vector<double>> rel_noise(n_groups);

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
        for (unsigned int cc = 0; cc < 2; cc++)
          dphi[g][cc] = delta_phi((g + cc * n_blocks_per_gr * n_groups) * n_dofs + i);

        phi[g][i] = phi_crit((g) * n_dofs + i);

        // Magnitud y fase
        delta_phi_mag[g][i] = sqrt(dphi[g][0] * dphi[g][0] + dphi[g][1] * dphi[g][1]);
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
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and the time-line.
 */
template <int dim, int n_fe_degree>
  void NoiseFullSPN<dim, n_fe_degree>::run ()
  {
    cout << "   Noise input files read." << std::endl;
    // Get Maximum memory
    PetscLogDouble memory;
    PetscMemorySetGetMaximumUsage();

    verbose_cout << "  setup_system..." << std::flush;
    setup_system();
    verbose_cout << "  Done! " << std::endl;

    assemble_system();

    cout << "      Number of DoFs per block: " << n_dofs << std::endl;
    cout << "      Number of Total DoFs: " << dof_handler.n_dofs() << std::endl;
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
      AssertRelease(false, "Not implemented");
      verbose_cout << "  output_detectors..." << std::flush;
      //output_detectors();
      verbose_cout << "  Done!" << std::endl;

    }
  }

template class NoiseFullSPN<1, 1> ;
template class NoiseFullSPN<1, 2> ;
template class NoiseFullSPN<1, 3> ;
template class NoiseFullSPN<1, 4> ;
template class NoiseFullSPN<1, 5> ;

template class NoiseFullSPN<2, 1> ;
template class NoiseFullSPN<2, 2> ;
template class NoiseFullSPN<2, 3> ;
template class NoiseFullSPN<2, 4> ;
template class NoiseFullSPN<2, 5> ;

template class NoiseFullSPN<3, 1> ;
template class NoiseFullSPN<3, 2> ;
template class NoiseFullSPN<3, 3> ;
template class NoiseFullSPN<3, 4> ;
template class NoiseFullSPN<3, 5> ;

