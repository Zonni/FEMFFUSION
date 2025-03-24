/*
 * @file   rom_kinetics.cc
 * @brief  Implementation of a (time dependent) Reduced Order Model
 * using the Proper Orthogonal decomposition
 */
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/qr.h>

#include <string>
#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscts.h>
#include <petscis.h>
#include <petscmat.h>
#include <slepcsvd.h>
#include <slepcbv.h>

#include "../../include/rom/rom_static.h"
#include "../../include/rom/ihs.h"
#include "../../include/utils.h"
#include "../../include/io/materials.h"
#include "../../include/io/perturbation.h"
#include "../../include/io/printing.h"
#include "../../include/eps_solvers/eps_solver.h"
#include "../../include/static_diffusion.h"

using namespace dealii;

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  ROMStatic<dim, n_fe_degree>::ROMStatic (
    ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &_static_problem,
    const bool verbose,
    const bool silent) :
      comm(MPI_COMM_WORLD),
      n_mpi_processes(Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
      verbose_cout(std::cout, verbose and this_mpi_process == 0),
      cout(std::cout, !silent and this_mpi_process == 0),
      n_groups(_static_problem.n_groups),
      n_dofs(_static_problem.n_dofs),
      dof_handler(_static_problem.dof_handler),
      constraints(_static_problem.constraints),
      boundary_conditions(_static_problem.boundary_conditions),
      n_assemblies(_static_problem.n_assemblies),
      perturbation(_static_problem.perturbation),
      static_problem(_static_problem),
      T(comm, dof_handler, constraints),
      F(comm, dof_handler, constraints),
      assem_per_dim(_static_problem.materials.assem_per_dim)
  {

    cout << std::endl << std::endl;
    cout << "---------------------------" << std::endl;
    cout << "-------  ROM - POD  -------" << std::endl;
    cout << "---------------------------" << std::endl;
    cout << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
      local_dofs_vector[g] = locally_owned_dofs;

    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Out parameters
    out_file = static_problem.out_file;
    get_bool_from_options("-out_flag", out_flag);
    out_flag = static_problem.out_flag;
    get_string_from_options("-out_file", out_file);
    n_out_ref = static_problem.n_out_ref;

    // Material parameters
    albedo_factors = static_problem.albedo_factors;

    // 	Time parameters
    type_perturbation = prm.get("Type_Perturbation");

    verbose_cout << "Initialize the perturbation class in ROM" << std::endl;

    perturbation.init_transient();

    perturbation_frac = prm.get_double("XS_Perturbation_Fraction");

    // Reinit ROM data
    n_snap = prm.get_integer("N_Snapshots");
    get_uint_from_options("-n_snap", n_snap);

    // Get type of snapshots
    type_snapshots = prm.get("ROM_Type_Snapshots");
    get_string_from_options("-type_snapshots", type_snapshots);

    // Get Group wise type
    rom_group_wise = prm.get("ROM_Group_Wise");
    get_string_from_options("-rom_group_wise", rom_group_wise);

    // Get snapshots
    cout << "   Snapshots type " << type_snapshots << std::endl;
    get_snapshots(static_problem, snapshots, type_snapshots);

    // LUPOD
    LUPOD_flag = prm.get_bool("LUPOD_Flag");
    get_bool_from_options("-LUPOD", LUPOD_flag);
    epsilon_N = prm.get_double("Epsilon_N");
    get_double_from_options("-epsilon_N", epsilon_N);
    epsilon_M = prm.get_double("Epsilon_M");
    get_double_from_options("-epsilon_M", epsilon_M);

    this->run();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::get_snapshots (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    std::string type_snap)
  {

    if (type_snap == "modes")
      get_snapshots_modes(static_problem, _snapshots);
    else if (type_snap == "bank12_ihs")
      get_snapshots_bar_ihs(static_problem, _snapshots);

    else if (type_snap == "IHS_XS")
      get_snapshots_IHS_XS(static_problem, perturbation_frac);
    else
      AssertRelease(false, "Invalid Type ROM of Snapshots");
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::get_snapshots_modes (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    n_snap = static_problem.n_eigenvalues;
    _snapshots.resize(n_snap);

    for (unsigned int ns = 0; ns < n_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    for (unsigned int ns = 0; ns < n_snap; ns++)
    {
      static_problem.phi[ns].compress(VectorOperation::insert);
      _snapshots[ns] = static_problem.phi[ns];
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::get_snapshots_bar (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    unsigned int bar_bank)
  {
    unsigned int min_s_snap = _snapshots.size();
    unsigned int max_s_snap = _snapshots.size() + n_snap;

    _snapshots.resize(max_s_snap);
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    double bars_bottom_pos = 0.0;
    AssertRelease(n_snap > 1, "The number of snapshots must be greater than 1");
    double step_bar = (perturbation.bars_top_pos - bars_bottom_pos)
                      / (n_snap - 1);

    unsigned int bar_material = perturbation.bar_materials[bar_bank - 1];

    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
    {

      unsigned int bar_pos_z = (ns - min_s_snap) * step_bar + bars_bottom_pos;
      for (unsigned int plant_pos = 0;
          plant_pos < perturbation.bars_position.size(); ++plant_pos)
      {
        if (perturbation.bars_position[plant_pos] == bar_bank)
        {
          perturbation.move_bar_volume_homogenized(plant_pos, bar_pos_z,
            bar_material, bar_bank - 1);
        }
      }

      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.assemble_system_lambda();
      static_problem.solve_eps();
      static_problem.phi[0].compress(VectorOperation::insert);
      _snapshots[ns] = static_problem.phi[0];
      cout << "     Step " << ns << ", Bar Position: " << bar_pos_z
      << ", Eigenvalue: "
      << std::setprecision(6)
      << static_problem.eigenvalues[0]
      << std::endl;
    }

    n_snap = max_s_snap;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::get_snapshots_IHS_XS (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    const double &frac_pert)
  {
    AssertRelease(frac_pert != 0.0,
      "If XS_Perturbation_Fraction is not set this does not make sense");
    const unsigned int n_mats = static_problem.materials.get_n_mats();
    // IHS - Parameters
    int duplication = 5;
    int *x;

    // Perturbate 7 XS of each Material: Sigma_r1, r2, s12, f1, f2, tr1 and tr2.
    const unsigned int n_dims_perturbated = n_mats * 7;
    //std::cout << " n_dims_perturbated " << n_dims_perturbated << std::endl;
    XS_org.reinit(n_dims_perturbated);
    Vector<double> pert_min(n_dims_perturbated), pert_max(n_dims_perturbated);
    AssertRelease(n_groups == 2, "Only Implemented for n_groups==2");
    for (unsigned int mat = 0; mat < n_mats; mat++)
    {
      XS_org[mat * 7 + 0] = static_problem.materials.get_sigma_r(0, mat);
      XS_org[mat * 7 + 1] = static_problem.materials.get_sigma_r(1, mat);
      XS_org[mat * 7 + 2] = static_problem.materials.get_sigma_s(0, 1, mat);
      XS_org[mat * 7 + 3] = static_problem.materials.get_nu_sigma_f(0, mat);
      XS_org[mat * 7 + 4] = static_problem.materials.get_nu_sigma_f(1, mat);
      XS_org[mat * 7 + 5] = static_problem.materials.get_sigma_tr(0, mat);
      XS_org[mat * 7 + 6] = static_problem.materials.get_sigma_tr(1, mat);
    }

    pert_min.equ(1 - frac_pert, XS_org);
    pert_max.equ(1 + frac_pert, XS_org);

    //  Get the IHS
    std::vector<Vector<double> > xs_sample(n_snap,
      Vector<double>(n_dims_perturbated));
    for (unsigned int ns = 0; ns < n_snap; ns++)
      xs_sample[ns].reinit(n_dims_perturbated);

    //std::cout << " ihs.. "  << std::endl;
    x = ihs(n_dims_perturbated, n_snap, duplication, seed);

    for (unsigned int ns = 0; ns < n_snap; ns++)
      for (unsigned int d = 0; d < n_dims_perturbated; d++)
      {
        //cout << "ns: " << ns << "  d: " << d << " " <<  pert_min[d] << "  XS: " << flush;
        xs_sample[ns][d] = pert_min[d]
                           + static_cast<double>(x[n_dims_perturbated * ns + d])
                             / n_snap
                             * (pert_max[d] - pert_min[d]);
        //std::cout << XS_org[d] << " -> " << xs_sample[ns][d] << std::endl;
      }

    delete[] x;

    // Resize snapshots Matrix
    snapshots.resize(n_snap);
    for (unsigned int ns = 0; ns < n_snap; ns++)
      snapshots[ns].reinit(local_dofs_vector, comm);

    cout << "  Creating " << n_snap << " snapshots with IHS and " << frac_pert
         << " perturbation fraction."
         << std::endl;
    // for each snapshot
    for (unsigned int ns = 0; ns < n_snap; ns++)
    {
      // Perturbate
      update_xs(static_problem.materials, xs_sample[ns]);

      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.assemble_system_lambda();
      static_problem.solve_eps();
      static_problem.phi[0].compress(VectorOperation::insert);
      // Normalize // Maybe it is necessary to run postprocess without printing anything
      static_problem.phi[0] *= static_problem.phi[0].mean_value();
      static_problem.phi[0] /= static_problem.phi[0].linfty_norm();
      snapshots[ns] = static_problem.phi[0];
      cout << "     Sample " << ns << ", Eigenvalue: " << std::setprecision(6)
      << static_problem.eigenvalues[0]
      << std::endl;
    }

    // Reset XS
    update_xs(static_problem.materials, XS_org);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::get_pertubation_random_XS (
    const unsigned int n_test,
    std::mt19937 &gen,
    std::vector<Vector<double> > &xs_sample)
  {
    AssertRelease(perturbation_frac != 0.0,
      "If XS_Perturbation_Fraction is not set this does not make sense");

    // Perturbate 7 XS of each Material: Sigma_r1, r2, s12, f1, f2, tr1 and tr2.
    const unsigned int n_mats = static_problem.materials.get_n_mats();
    const unsigned int n_dims_perturbated = n_mats * 7;

    xs_sample.resize(n_test, Vector<double>(n_dims_perturbated));
    for (unsigned int t = 0; t < n_test; t++)
      xs_sample[t].reinit(n_dims_perturbated);

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    Vector<double> pert_min(n_dims_perturbated), pert_max(n_dims_perturbated);

    pert_min.equ(1 - perturbation_frac, XS_org);
    pert_max.equ(1 + perturbation_frac, XS_org);

    for (unsigned int t = 0; t < n_test; t++)
      for (unsigned int d = 0; d < n_dims_perturbated; ++d)
        xs_sample[t][d] = pert_min[d] + dist(gen) * (pert_max[d] - pert_min[d]);

    cout << "   ------------------------------------   " << std::endl;
    cout << std::endl;
    cout << "   " << n_test << " Tests with random XS created and "
         << std::setprecision(2)
         << perturbation_frac
         << " perturbation fraction."
         << std::setprecision(6)
         << std::endl;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::get_snapshots_bar_ihs (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {

    int n_bars = 2;
    // IHS
    int duplication = 5;
    int seed = 17;
    int *x;
    std::vector<std::vector<double>> sample_pos(n_bars,
      std::vector<double>(n_snap));

    //
    //  Get the points.
    //
    x = ihs(n_bars, n_snap, duplication, seed);
    for (int d = 0; d < n_bars; d++)
      for (unsigned int i = 0; i < n_snap; i++)
      {
        sample_pos[d][i] = static_cast<double>(x[n_bars * i + d]) / n_snap;
      }

    delete[] x;

    // Compute Snapshots
    unsigned int min_s_snap = _snapshots.size();
    unsigned int max_s_snap = _snapshots.size() + n_snap;

    _snapshots.resize(max_s_snap);
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    // for each snapshot
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
    {
      //for each bar bank
      for (unsigned int nb = 1; nb < static_cast<unsigned int>(n_bars) + 1; nb++)
      {

        double bars_bottom_pos = 0.0;
        unsigned int bar_material = perturbation.bar_materials[nb - 1];

        unsigned int bar_pos_z = sample_pos[nb - 1][ns]
                                 * (perturbation.bars_top_pos - bars_bottom_pos);

        cout << "bank_bar" << nb << "bar_pos: " << bar_pos_z << std::endl;

        for (unsigned int plant_pos = 0;
            plant_pos < perturbation.bars_position.size(); ++plant_pos)
        {
          if (perturbation.bars_position[plant_pos] == nb)
          {
            perturbation.move_bar_volume_homogenized(plant_pos,
              bar_pos_z, bar_material, nb - 1);
          }
        }
      }

      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.assemble_system_lambda();
      static_problem.solve_eps();
      static_problem.phi[0].compress(VectorOperation::insert);
      _snapshots[ns] = static_problem.phi[0];
      cout << "     Step " << ns
      << ", Eigenvalue: "
      << static_problem.eigenvalues[0]
      << std::endl;

    }
    n_snap = max_s_snap;
  }

/*
 * @brief Compute POD basis
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::compute_POD_basis_monolithic (
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {

    cout << "   POD ----- Monolithic" << std::endl;
    // TODO SET AN ENERGY THRSHOLD
    dim_rom = _snapshots.size();
    get_uint_from_options("-dim_rom", dim_rom);

    snap_basis.resize(dim_rom);
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      snap_basis[dr].reinit(local_dofs_vector, comm);

    // Create the matrix with the snapshot to apply the SVD
    Mat Mat_snap;
    PetscInt i_snap, i_sv;
    PetscInt *idm = new PetscInt[n_dofs * n_groups];
    PetscScalar *values_snap = new PetscScalar[n_dofs * n_groups];

    for (unsigned int j = 0; j < n_dofs * n_groups; j++)
      idm[j] = j;

    MatCreateDense(comm, n_dofs * n_groups, dim_rom, n_dofs * n_groups, dim_rom, NULL,
      &Mat_snap);
    for (i_snap = 0; i_snap < static_cast<int>(dim_rom); i_snap++)
    {
      for (unsigned k = 0; k < n_dofs * n_groups; k++)
        values_snap[k] = _snapshots[i_snap][k];
      MatSetValues(Mat_snap, n_dofs * n_groups, idm, 1, &i_snap, values_snap,
        INSERT_VALUES);
    }
    MatAssemblyBegin(Mat_snap, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Mat_snap, MAT_FINAL_ASSEMBLY);

    SVD svd;
    SVDCreate(comm, &svd);
    SVDSetOperators(svd, Mat_snap, NULL);
    SVDSetDimensions(svd, dim_rom, 2 * dim_rom, dim_rom);
    SVDSetProblemType(svd, SVD_STANDARD);
    SVDSetType(svd, SVDTRLANCZOS);

    SVDSetFromOptions(svd);
    SVDSolve(svd);

    PETScWrappers::MPI::Vector u;
    u.reinit(comm, n_dofs * n_groups, n_dofs * n_groups);
    std::vector<double> singular_values(dim_rom);
    for (i_sv = 0; i_sv < static_cast<int>(dim_rom); i_sv++)
    {
      SVDGetSingularTriplet(svd, i_sv, &(singular_values[i_sv]), u, NULL);
      copy_to_BlockVector(snap_basis[i_sv], u);
      cout << "      Singular value " << i_sv << ": " << std::scientific
           << singular_values[i_sv]
           << std::fixed << std::endl;
    }

    MatDestroy(&Mat_snap);
    SVDDestroy(&svd);
    u.clear();
  }

/*
 * @brief Compute POD basis
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::compute_POD_basis_group_wise (
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    cout << "   POD ---- GROUP WISE " << std::endl;
    unsigned int J = n_dofs;
    unsigned int K = _snapshots.size();
    dim_rom = n_groups * _snapshots.size();
    get_uint_from_options("-dim_rom", dim_rom);

    cout << "   K " << K << " n_dofs " << n_dofs << std::endl;

    // Create the matrix with the snapshot to apply the SVD
    std::vector<LAPACKFullMatrix<double> > S(n_groups, LAPACKFullMatrix<double>(J, K));
    for (unsigned int g = 0; g < n_groups; g++)
      S[g].reinit(J, K);

    // Copy to Full Matrix
    for (unsigned int g = 0; g < n_groups; g++)
      for (unsigned int k = 0; k < K; k++)
        for (unsigned int j = 0; j < J; j++)
          S[g](j, k) = _snapshots[k].block(g)[j];

    // Create U
    std::vector<LAPACKFullMatrix<double> > U(n_groups, LAPACKFullMatrix<double>(J, K));
    for (unsigned int g = 0; g < n_groups; g++)
      U[g].reinit(J, K);

    std::vector<Vector<double> > Sigma(n_groups, Vector<double>(K));
    for (unsigned int g = 0; g < n_groups; g++)
      Sigma[g].reinit(K, true);
    for (unsigned int g = 0; g < n_groups; g++)
    {
      S[g].compute_svd();
      // Retrieve values from SVD
      U[g] = S[g].get_svd_u();

      for (unsigned int i = 0; i < K; i++)
      {
        Sigma[g](i) = S[g].singular_value(i);
        cout << "   Singular value g" << g << "  index " << i << ": "
             << std::scientific
             << Sigma[g](i)    //        << std::fixed
             << std::endl;
      }
      cout << std::fixed;
    }

    // Determine the number of modes to retain based on epsilon_M
    //FIXME
    //    double total_energy = Sigma.l2_norm();
    //    unsigned int M = 0;
    //    for (unsigned int m = 0; m < N; ++m)
    //    {
    //      double leftout_energy = 0.0;
    //      for (unsigned int m2 = m + 1; m2 < N; ++m2)
    //        leftout_energy += Sigma(m2) * Sigma(m2);
    //      leftout_energy = sqrt(leftout_energy);
    //
    //      if (leftout_energy / total_energy <= epsilon_M)
    //      {
    //        M = m + 1;
    //        break;
    //      }
    //    }
    //dim_rom = M;
    // ---------------------------------------

    // Copy snap_basis
    snap_basis.resize(dim_rom);
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      snap_basis[dr].reinit(local_dofs_vector, comm);

    for (unsigned int g = 0; g < n_groups; g++)
      for (unsigned int k = 0; k < K; k++)
      {
        for (unsigned int j = 0; j < J; j++)
        {
          snap_basis[k + g * K].block(g)(j) = U[g](j, k);
        }
        snap_basis[k].compress(VectorOperation::insert);
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::compute_LUPOD_basis_monolithic (
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    unsigned int J = n_dofs * n_groups;
    unsigned int K = _snapshots.size();

    cout << "   LUPOD ---- MONOLITHIC " << std::endl;
    //cout << " J " << J << "  K " << K << " SIZE " << _snapshots[0].size() << std::endl;
    // Create the matrix with the snapshot to apply the SVD
    FullMatrix<double> S(J, K);
    // Copy to Full Matrix
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < J; j++)
        S(j, k) = _snapshots[k][j];

    //cout << "S.print_formatted " << std::endl;
    //S.print_formatted(std::cout, 18, true, 0, "0.0");

    //
    // --------------------------------------------------------------------
    // LUPOD Technology
    // Inputs:
    // S - Snapshot matrix (J x K)
    // epsilon_N - threshold for selection of collocation points
    // epsilon_M - threshold for retention of modes in POD
    FullMatrix<double> S_mod(S);

    snaps.reserve(K);
    points.reserve(K);

    cout << "   epsilon_N  " << std::scientific << epsilon_N << std::endl;
    cout << "   epsilon_M  " << std::scientific << epsilon_M << std::endl;
    unsigned int j_max = 0, k_max = 0, it = 0;
    double pivot, max_val, error = 1e6;
    while (error > epsilon_N)
    {
      // Step 1.1: Find the index of the largest absolute value in S_mod
      max_val = 0.0;
      for (unsigned int j = 0; j < J; j++)
        for (unsigned int k = 0; k < K; k++)
          if (std::abs(S_mod(j, k)) > max_val)
          {
            max_val = std::abs(S_mod(j, k));
            j_max = j;
            k_max = k;
          }

      std::vector<double> temp_col(J);
      // Step 1.2: Perform Gaussian elimination to set the j_index-th row to zero
      pivot = S_mod(j_max, k_max);
      for (unsigned int j = 0; j < J; j++)
        S_mod(j, k_max) /= pivot;

      // Update other columns
      // in FullMatrix<double> could be done with add_row
      for (unsigned int k = 0; k < K; ++k)
        if (k != k_max)
        {
          for (unsigned int j = 0; j < J; ++j)
            temp_col[j] = S_mod(j, k) - S_mod(j_max, k) * S_mod(j, k_max);

          for (unsigned int j = 0; j < J; ++j)
            S_mod(j, k) = temp_col[j];
        }

      // Store the indices
      points.push_back(j_max);
      snaps.push_back(k_max);

      // Remove the selected snapshot
      for (unsigned int j = 0; j < J; ++j)
        S_mod(j, k_max) = 0;

      //cout << "j_max: " << j_max << "  k_max: " << k_max << std::endl;

      //cout << "S_mod " << std::endl;
      //S_mod.print_formatted(std::cout, 7, true, 0, "0.0");
      // Step 3: Check for convergence based on the Frobenius norm
      error = S_mod.frobenius_norm() / S.frobenius_norm();

      //cout << "error: " << std::scientific << error << std::fixed << std::endl;
      AssertRelease(it <= J, "Something went wrong in LUPOD");
      it++;
    }

    unsigned int N = snaps.size();  // Number of collocation points

    // Step 4: Perform SVD on the reduced snapshot matrix
    LAPACKFullMatrix<double> S_reduced(points.size(), snaps.size());
    for (unsigned int j = 0; j < N; j++)
      for (unsigned int k = 0; k < N; k++)
        S_reduced(j, k) = S(points[j], snaps[k]);

    //S_reduced.print_formatted(std::cout, 14, true, 0, "0.0", 1.0, 1e-15);
    LAPACKFullMatrix<double> vt;
    S_reduced.compute_svd();
    // Retrieve values from SVD
    LAPACKFullMatrix<double> U_reduced = S_reduced.get_svd_u();
    vt = S_reduced.get_svd_vt();
    Vector<double> Sigma(N);
    LAPACKFullMatrix<double> Sigma_mat_inv(N, N);
    for (unsigned int i = 0; i < N; i++)
    {
      Sigma(i) = S_reduced.singular_value(i);
      cout << "   Singular value " << i << ": " << std::scientific << Sigma(i)
           << std::fixed
           << std::endl;
      Sigma_mat_inv(i, i) = 1 / Sigma(i);
    }

    cout << "   LUPOD Points:  " << std::flush;
    print_vector(points);

    // Determine the number of modes to retain based on epsilon_M
    double total_energy = Sigma.l2_norm();
    unsigned int M = 0;
    for (unsigned int m = 0; m < N; ++m)
    {
      double leftout_energy = 0.0;
      for (unsigned int m2 = m + 1; m2 < N; ++m2)
        leftout_energy += Sigma(m2) * Sigma(m2);
      leftout_energy = sqrt(leftout_energy);

      if (leftout_energy / total_energy <= epsilon_M)
      {
        M = m + 1;
        break;
      }
    }

    // Perform QR decomposition using LAPACK
    LAPACKFullMatrix<double> Q, R, aux(N, N), R_inv;
    FullMatrix<double> P(N, N);
    qr(U_reduced, Q, R);

    // P = V*inv(Sigma_mat)*inv(R);
    R.invert();

    vt.Tmmult(aux, Sigma_mat_inv); // We get V transposed by get_svd_vt();
    aux.mmult(P, R);

    FullMatrix<double> S2(J, M);
    FullMatrix<double> U_full(J, M);

    // S2= S(:, snaps)
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int k = 0; k < M; k++)
        S2(j, k) = S(j, snaps[k]);

    // U_full = S(:, snaps)*P;
    S2.mmult(U_full, P);

    //    cout << "U_full: " << std::endl;
    //    U_full.print_formatted(std::cout, 6, true);
    //    cout << "U_reduced " << std::endl;
    //    U_reduced.print_formatted(std::cout, 6, true);

    // ---------------------------------------
    // Copy snap_basis_red
    dim_rom = M;
    snap_basis_red.resize(M, Vector<double>(N));
    for (unsigned int j = 0; j < M; ++j)
    {
      snap_basis_red[j].reinit(N);
      for (unsigned int i = 0; i < N; ++i)
      {
        snap_basis_red[j](i) = U_reduced(i, j);
      }
    }

    // Copy snap_basis
    snap_basis.resize(dim_rom);
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      snap_basis[dr].reinit(local_dofs_vector, comm);
    for (unsigned int j = 0; j < dim_rom; j++)
    {
      for (unsigned int i = 0; i < J; i++)
      {
        snap_basis[j](i) = U_full(i, j);
      }
      snap_basis[j].compress(VectorOperation::insert);
    }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::compute_LUPOD_basis_group_wise (
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    unsigned int J = n_dofs * n_groups;
    unsigned int K = _snapshots.size();

    cout << "   LUPOD ---- GROUP WISE " << std::endl;
    //cout << " J " << J << "  K " << K << " SIZE " << _snapshots[0].size() << std::endl;
    // Create the matrix with the snapshot to apply the SVD
    FullMatrix<double> S(J, K);
    // Copy to Full Matrix
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < J; j++)
        S(j, k) = _snapshots[k][j];

    //cout << "S.print_formatted " << std::endl;
    //S.print_formatted(std::cout, 18, true, 0, "0.0");

    //
    // --------------------------------------------------------------------
    // LUPOD Technology
    // Inputs:
    // S - Snapshot matrix (J x K)
    // epsilon_N - threshold for selection of collocation points
    // epsilon_M - threshold for retention of modes in POD
    FullMatrix<double> S_mod(S);

    snaps.reserve(K);
    points.reserve(K);

    cout << "   epsilon_N  " << std::scientific << epsilon_N << std::endl;
    cout << "   epsilon_M  " << std::scientific << epsilon_M << std::endl;
    unsigned int j_max = 0, k_max = 0, it = 0;
    double pivot, max_val, error = 1e6;
    while (error > epsilon_N)
    {
      // Step 1.1: Find the index of the largest absolute value in S_mod
      max_val = 0.0;
      for (unsigned int j = 0; j < J; j++)
        for (unsigned int k = 0; k < K; k++)
          if (std::abs(S_mod(j, k)) > max_val)
          {
            max_val = std::abs(S_mod(j, k));
            j_max = j;
            k_max = k;
          }

      std::vector<double> temp_col(J);
      // Step 1.2: Perform Gaussian elimination to set the j_index-th row to zero
      pivot = S_mod(j_max, k_max);
      for (unsigned int j = 0; j < J; j++)
        S_mod(j, k_max) /= pivot;

      // Update other columns
      // in FullMatrix<double> could be done with add_row
      for (unsigned int k = 0; k < K; ++k)
        if (k != k_max)
        {
          for (unsigned int j = 0; j < J; ++j)
            temp_col[j] = S_mod(j, k) - S_mod(j_max, k) * S_mod(j, k_max);

          for (unsigned int j = 0; j < J; ++j)
            S_mod(j, k) = temp_col[j];
        }

      // Store the indices
      points.push_back(j_max);
      snaps.push_back(k_max);

      // Remove the selected snapshot
      for (unsigned int j = 0; j < J; ++j)
        S_mod(j, k_max) = 0;

      //std::cout << "j_max: " << j_max << "  k_max: " << k_max << std::endl;

      //std::cout << "S_mod " << std::endl;
      //S_mod.print_formatted(std::cout, 7, true, 0, "0.0");
      // Step 3: Check for convergence based on the Frobenius norm
      error = S_mod.frobenius_norm() / S.frobenius_norm();

      //std::cout << "error: " << std::scientific << error << std::fixed << std::endl;
      AssertRelease(it <= J, "Something went wrong in LUPOD");
      it++;
    }

    unsigned int N = snaps.size();  // Number of collocation points

    // Step 4: Perform SVD on the reduced snapshot matrix
    LAPACKFullMatrix<double> S_reduced(points.size(), snaps.size());
    for (unsigned int j = 0; j < N; j++)
      for (unsigned int k = 0; k < N; k++)
        S_reduced(j, k) = S(points[j], snaps[k]);

    //S_reduced.print_formatted(std::cout, 14, true, 0, "0.0", 1.0, 1e-15);
    LAPACKFullMatrix<double> vt;
    S_reduced.compute_svd();

    // Retrieve values from SVD
    LAPACKFullMatrix<double> U_reduced = S_reduced.get_svd_u();
    vt = S_reduced.get_svd_vt();
    Vector<double> Sigma(N);
    LAPACKFullMatrix<double> Sigma_mat_inv(N, N);
    for (unsigned int i = 0; i < N; i++)
    {
      Sigma(i) = S_reduced.singular_value(i);
      cout << "   Singular value " << i << ": " << std::scientific << Sigma(i)
           << std::fixed
           << std::endl;
      Sigma_mat_inv(i, i) = 1 / Sigma(i);
    }

    cout << "   LUPOD Points:  " << std::flush;
    print_vector(points);

    // Determine the number of modes to retain based on epsilon_M
    double total_energy = Sigma.l2_norm();
    unsigned int M = 0;
    for (unsigned int m = 0; m < N; ++m)
    {
      double leftout_energy = 0.0;
      for (unsigned int m2 = m + 1; m2 < N; ++m2)
        leftout_energy += Sigma(m2) * Sigma(m2);
      leftout_energy = sqrt(leftout_energy);

      if (leftout_energy / total_energy <= epsilon_M)
      {
        M = m + 1;
        break;
      }
    }

    // Perform QR decomposition using LAPACK
    LAPACKFullMatrix<double> Q, R, aux(N, N), R_inv;
    FullMatrix<double> P(N, N);
    qr(U_reduced, Q, R);

    // P = V*inv(Sigma_mat)*inv(R);
    R.invert();

    vt.Tmmult(aux, Sigma_mat_inv); // We get V transposed by get_svd_vt();
    aux.mmult(P, R);

    FullMatrix<double> S2(J, M);
    FullMatrix<double> U_full(J, M);

    // S2= S(:, snaps)
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int k = 0; k < M; k++)
        S2(j, k) = S(j, snaps[k]);

    // U_full = S(:, snaps)*P;
    S2.mmult(U_full, P);

    //    std::cout << "U_full: " << std::endl;
    //    U_full.print_formatted(std::cout, 6, true);
    //    std::cout << "U_reduced " << std::endl;
    //    U_reduced.print_formatted(std::cout, 6, true);

    // ---------------------------------------
    // Copy snap_basis_red
    dim_rom = M;
    snap_basis_red.resize(M, Vector<double>(N));
    for (unsigned int j = 0; j < M; ++j)
    {
      snap_basis_red[j].reinit(N);
      for (unsigned int i = 0; i < N; ++i)
      {
        snap_basis_red[j](i) = U_reduced(i, j);
      }
    }

    // Copy snap_basis
    snap_basis.resize(dim_rom);
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      snap_basis[dr].reinit(local_dofs_vector, comm);
    for (unsigned int j = 0; j < dim_rom; j++)
    {
      for (unsigned int i = 0; i < J; i++)
      {
        snap_basis[j](i) = U_full(i, j);
      }
      snap_basis[j].compress(VectorOperation::insert);
    }
  }

/*
 * @brief Perform QR decomposition using LAPACK
 * @input A
 * @output Q and R
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::qr (
    const LAPACKFullMatrix<double> &A,
    LAPACKFullMatrix<double> &Q,
    LAPACKFullMatrix<double> &R)
  {
    const int m = A.m();
    const int n = A.n();

    int info;
    std::vector<double> tau(n); // Householder reflectors
    int lwork = -1;
    double query_work;

    // Manually extract matrix data into column-major order
    std::vector<double> A_data(A.begin(), A.end());

    // Query optimal workspace size
    dgeqrf_(&m, &n, A_data.data(), &n, tau.data(), &query_work, &lwork, &info);
    lwork = static_cast<int>(query_work);
    std::vector<double> work(lwork);

    // Compute QR decomposition
    dgeqrf_(&m, &n, A_data.data(), &n, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
      std::cerr << "QR factorization failed with error code: " << info << std::endl;
      exit(-1);
    }

    // Extract R (upper triangular part)
    R.reinit(m, n);
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        R(i, j) = (j >= i and abs(A_data[i + j * n]) > 1e-15) ? A_data[i + j * n] : 0.0; // Store R in upper part

      // Compute Q explicitly
    dorgqr_(&m, &n, &n, A_data.data(), &n, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
      std::cerr << "Q computation failed with error code: " << info << std::endl;
      exit(-1);
    }

    // Store Q in a deal.II matrix
    Q.reinit(m, n);
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        Q(i, j) = A_data[i + j * n];
  }

/*
 * @brief Update XS.
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::update_xs (
    Materials &materials,
    const Vector<double> &xs)
  {
    for (unsigned int mat = 0; mat < materials.get_n_mats(); mat++)
    {
      materials.set_sigma_r(xs[mat * 7 + 0], 0, mat);
      materials.set_sigma_r(xs[mat * 7 + 1], 1, mat);
      materials.set_sigma_s(xs[mat * 7 + 2], 0, 1, mat);
      materials.set_nu_sigma_f(xs[mat * 7 + 3], 0, mat);
      materials.set_nu_sigma_f(xs[mat * 7 + 4], 1, mat);
      materials.set_sigma_tr(xs[mat * 7 + 5], 0, mat);
      materials.set_sigma_tr(xs[mat * 7 + 6], 1, mat);
    }
  }

/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::assemble_ROM_matrices ()
  {

    MatrixFreeType matrixfree_type = full_allocated;
    // Allocate and assemble block matrices
    T.reinit(static_problem.materials, boundary_conditions, albedo_factors,
      matrixfree_type);
    F.reinit(static_problem.materials, matrixfree_type);

    MatCreateSeqDense(PETSC_COMM_SELF, dim_rom, dim_rom, NULL, &romT);
    MatCreateSeqDense(PETSC_COMM_SELF, dim_rom, dim_rom, NULL, &romF);

    if (LUPOD_flag)
    {
      //std::cout << "LUPOD TRUE" << std::endl;
      //assemble_matrices();

      //romT = FullMatrix<double>(dim_rom);
      //romF = FullMatrix<double>(dim_rom);
      unsigned int n_points = points.size();

      Vector<double> auxT(n_points), auxF(n_points);

      for (unsigned int b2 = 0; b2 < dim_rom; b2++)
      {
        for (unsigned int i = 0; i < n_points; i++)
        {
          T.vmult_row(auxT[i], snap_basis[b2], points[i]);
          F.vmult_row(auxF[i], snap_basis[b2], points[i]);
        }

        for (unsigned int b1 = 0; b1 < dim_rom; b1++)
        {
          MatSetValue(romT, b1, b2, auxT * snap_basis_red[b1], INSERT_VALUES);
          MatSetValue(romF, b1, b2, auxF * snap_basis_red[b1], INSERT_VALUES);
        }
      }
    }
    else
    {
      //std::cout << "LUPOD FALSE" << std::endl;

      PETScWrappers::MPI::BlockVector auxT(snap_basis[0]);
      PETScWrappers::MPI::BlockVector auxF(snap_basis[0]);

      for (unsigned int b2 = 0; b2 < dim_rom; b2++)
      {
        T.vmult(auxT, snap_basis[b2]);
        F.vmult(auxF, snap_basis[b2]);

        for (unsigned int b1 = 0; b1 < dim_rom; b1++)
        {
          MatSetValue(romT, b1, b2, auxT * snap_basis[b1], INSERT_VALUES);
          MatSetValue(romF, b1, b2, auxF * snap_basis[b1], INSERT_VALUES);
        }
      }
    }

    // Assemble the matrix (mandatory after setting values)
    MatAssemblyBegin(romT, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(romT, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(romF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(romF, MAT_FINAL_ASSEMBLY);

    F.clear();
    T.clear();
  }

/*
 * @brief Compute the largest eigenvalue with FEMFFUSION
 * @input
 * @output
 */
template <int dim, int n_fe_degree>
  double ROMStatic<dim, n_fe_degree>::solve_eigenvalue (PETScWrappers::MPI::BlockVector &phi)
  {
    EPS eps;
    PETScWrappers::MPI::Vector phi_red(comm, dim_rom, dim_rom);

    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, romF, romT);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    //EPSSetType(eps, EPSLAPACK);
    EPSSolve(eps);

    double eigr, eigi;
    EPSGetEigenpair(eps, 0, &eigr, &eigi, phi_red, NULL);
    AssertRelease(abs(eigi) < 1e-12,
      " The largest  eigenvalue of the reactor is complex!");

    for (unsigned int j = 0; j < n_dofs * n_groups; j++)
      for (unsigned int k = 0; k < dim_rom; k++)
      {
        phi[j] += snap_basis[k][j] * phi_red[k];
      }
    phi.compress(VectorOperation::add);

    EPSDestroy(&eps);
    MatDestroy(&romT);
    MatDestroy(&romF);

    return eigr;
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::run ()
  {
    cout << std::endl;
    cout << "------------ START OF THE ROM TIME LOOP ------------------"
         << std::endl;

    cout << "   Type of perturbation: " << type_perturbation << std::endl
         << std::endl;

    cout << std::fixed
         << "   Compute POD basis...                CPU Time = "
         << timer.cpu_time() << " s." << std::endl;

    if (rom_group_wise == "Monolithic")
    {
      if (LUPOD_flag)
        compute_LUPOD_basis_monolithic(snapshots);
      else
        compute_POD_basis_monolithic(snapshots);
    }
    else if (rom_group_wise == "Group_Wise")
    {
      if (LUPOD_flag)
        compute_LUPOD_basis_group_wise(snapshots);
      else
        compute_POD_basis_group_wise(snapshots);
    }
    else
      AssertRelease(false, "rom_group_wise must be Monolithic or Group_Wise");

    double time_get_snap = timer.cpu_time();
    cout << "   Time to get and process snapshots: " << time_get_snap << " s."
         << std::endl;
    //    -------------------------------------------------------------------------    //
    // Create Tests
    const unsigned int n_test = 10;
    std::vector<Vector<double> > xs_test;
    std::mt19937 gen(seed); // Random number generator for uniform distribution [0.0,1.0]
    get_pertubation_random_XS(n_test, gen, xs_test);

    //    -------------------------------------------------------------------------    //
    Vector<double> eig_fom(n_test);
    phi_fom.resize(n_test, PETScWrappers::MPI::BlockVector(local_dofs_vector, comm));
    for (unsigned int t = 0; t < n_test; t++)
      phi_fom[t].reinit(local_dofs_vector, comm);
    for (unsigned int t = 0; t < n_test; t++)
    {
      verbose_cout << "Update XS... " << std::flush;
      update_xs(static_problem.materials, xs_test[t]);
      verbose_cout << " Done!" << std::endl;

      //static_problem.tol_eps =
      //static_problem.tol_ksp =
      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.assemble_system_lambda();
      static_problem.solve_eps();
      static_problem.phi[0].compress(VectorOperation::insert);
      // Normalize // Maybe it is necessary to run postprocess without printing anything
      static_problem.phi[0] *= static_problem.phi[0].mean_value();
      static_problem.phi[0] /= static_problem.phi[0].linfty_norm();
      phi_fom[t] = static_problem.phi[0];
      eig_fom[t] = static_problem.eigenvalues[0];
      cout << "      Test FOM " << t << ", Eigenvalue: " << std::setprecision(5)
      << static_problem.eigenvalues[0]
      << std::endl;
    }

    double time_fom = timer.cpu_time() - time_get_snap;
    cout << "   Time FOM: " << time_fom << " s." << std::endl;
    cout << "   ------------------------------------   " << std::endl;
    cout << std::endl;
    //    -------------------------------------------------------------------------    //
    // ROM
    Vector<double> eig_rom(n_test);
    phi_rom.resize(n_test, PETScWrappers::MPI::BlockVector(local_dofs_vector, comm));
    for (unsigned int t = 0; t < n_test; t++)
      phi_rom[t].reinit(local_dofs_vector, comm);

    for (unsigned int t = 0; t < n_test; t++)
    {
      verbose_cout << "Update XS... " << std::flush;
      update_xs(static_problem.materials, xs_test[t]);
      verbose_cout << " Done!" << std::endl;

      verbose_cout << "assemble_ROM_matrices... " << std::flush;
      assemble_ROM_matrices();
      verbose_cout << " Done!" << std::endl;
      verbose_cout << "   Solve the ROM system...                 CPU Time = "
                   << timer.cpu_time()
                   << " s." << std::endl;

      eig_rom[t] = solve_eigenvalue(phi_rom[t]);
      // Postprocess
      phi_rom[t].compress(VectorOperation::insert);
      phi_rom[t] *= phi_rom[t].mean_value();
      phi_rom[t] /= phi_rom[t].linfty_norm();

      cout << "      Test ROM " << t << ", Eigenvalue: " << std::setprecision(5)
           << eig_rom[t]
           << std::endl;
      verbose_cout << " Done!" << std::endl;
    }
    double time_rom = timer.cpu_time() - (time_fom + time_get_snap);
    cout << "   Time ROM: " << time_rom << " s." << std::endl;
    cout << "   ------------------------------------   " << std::endl;
    cout << std::endl;
    // -------------------------------------------------------------------------    //
    // Postprocess FOM vs ROM

    Vector<double> delta_keff(n_test);
    Vector<double> rms_phi(n_test);
    for (unsigned int t = 0; t < n_test; t++)
    {
      delta_keff[t] = 1e5 * std::abs(eig_rom(t) - eig_fom(t));
      // RMS phi
      for (unsigned int i = 0; i < n_dofs * n_groups; i++)
        rms_phi[t] += std::pow((phi_rom[t][i] - phi_fom[t][i]), 2);
      rms_phi[t] = 100 * sqrt(rms_phi[t] / n_dofs * n_groups);
    }
    cout << std::scientific;
    cout << "   Mean Delta Keff = " << delta_keff.mean_value() << " pcm." << std::endl;
    cout << "   Max  Delta Keff = " << delta_keff.linfty_norm() << " pcm." << std::endl;
    cout << "   Mean RMS Phi = " << rms_phi.mean_value() << " %." << std::endl;
    cout << "   Max  RMS Phi = " << rms_phi.linfty_norm() << " %." << std::endl;
    cout << std::fixed << std::endl;
    cout << "   Finished in " << timer.cpu_time() << " s." << std::endl;
  }

template class ROMStatic<1, 1> ;
template class ROMStatic<1, 2> ;
template class ROMStatic<1, 3> ;
template class ROMStatic<1, 4> ;
template class ROMStatic<1, 5> ;

template class ROMStatic<2, 1> ;
template class ROMStatic<2, 2> ;
template class ROMStatic<2, 3> ;
template class ROMStatic<2, 4> ;
template class ROMStatic<2, 5> ;
//
template class ROMStatic<3, 1> ;
template class ROMStatic<3, 2> ;
template class ROMStatic<3, 3> ;
template class ROMStatic<3, 4> ;
template class ROMStatic<3, 5> ;

