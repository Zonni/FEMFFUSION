/*
 * @file   rom_static.cc
 * @brief  Implementation of a the steady_state Reduced Order Model
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
#include "../../include/rom/rom_utils.h"
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
      n_assemblies(_static_problem.n_assemblies),
      perturbation(_static_problem.perturbation),
      static_problem(_static_problem),
      T(comm, dof_handler, _static_problem.constraints),
      F(comm, dof_handler, _static_problem.constraints),
      assem_per_dim(_static_problem.materials.assem_per_dim),
      Tred(comm, dof_handler, _static_problem.constraints),
      Fred(comm, dof_handler, _static_problem.constraints)
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

    // 	Perturbation Parameters
    type_perturbation = prm.get("Type_Perturbation");

    verbose_cout << "Initialize the perturbation class in ROM" << std::endl;

    perturbation.init_transient();

    perturbation_frac = prm.get_double("XS_Perturbation_Fraction");
    get_double_from_options("-xs_perturbation_fraction", perturbation_frac);

    // Reinit ROM data
    n_snap = prm.get_integer("N_Snapshots");
    get_uint_from_options("-n_snap", n_snap);

    // Reinit ROM data
    n_test = prm.get_integer("N_Test");
    get_uint_from_options("-n_test", n_test);

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
    LUPOD_type = prm.get("LUPOD_Type");
    get_string_from_options("-lupod_type", LUPOD_type);
    epsilon_N = prm.get_double("Epsilon_N");
    get_double_from_options("-epsilon_N", epsilon_N);
    epsilon_M = prm.get_double("Epsilon_M");
    get_double_from_options("-epsilon_M", epsilon_M);
    M_req = 0;
    get_uint_from_options("-m_req", M_req); //

    n_LUPOD_points = prm.get_integer("N_LUPOD_Points");
    get_uint_from_options("-n_lupod_points", n_LUPOD_points);

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

    cout << "    Creating " << n_snap << " snapshots with IHS and " << frac_pert
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
 * @brief Prepare ROM
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::prepare_ROM ()
  {

    matrixfree_type = full_allocated;
    MatCreateSeqDense(PETSC_COMM_SELF, dim_rom, dim_rom, NULL, &romT);
    MatCreateSeqDense(PETSC_COMM_SELF, dim_rom, dim_rom, NULL, &romF);

    if (LUPOD_type != "POD")
    {
      //std::cout << " points_per_block " << std::endl;
      //print_vector(points_per_block[0]);
      //print_vector(points_per_block[1]);
      n_LUPOD_points = points_per_block[0].size();

      // Assemble block matrices
      Tred.reinit(static_problem.materials, points_per_block, matrixfree_type);
      Fred.reinit(static_problem.materials, points_per_block, matrixfree_type);
    }

//    else // POD - It is done in assemble_ROM_matrices ()
//    {
//      // Allocate and assemble block matrices
//      T.reinit(static_problem.materials, static_problem.boundary_conditions,
//        static_problem.albedo_factors, matrixfree_type);
//      F.reinit(static_problem.materials, matrixfree_type);
//    }
  }

/*
 * @brief Assemble ROM matrices
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::assemble_ROM_matrices ()
  {

    std::cout << "  dim_rom  " << dim_rom << std::endl;
    if (LUPOD_type == "POD")
    {
      // Allocate and assemble block matrices
      T.reinit(static_problem.materials, static_problem.boundary_conditions,
        static_problem.albedo_factors, matrixfree_type);
      F.reinit(static_problem.materials, matrixfree_type);

      //      std::cout << "Matrix Reduced 0 0" << std::endl;
      //      MatView(*(T.matrix_blocks[0][0]), PETSC_VIEWER_STDOUT_SELF);
      //      std::cout << "Matrix Reduced 0 1" << std::endl;
      //      MatView(*(T.matrix_blocks[0][1]), PETSC_VIEWER_STDOUT_SELF);
      //      std::cout << "Matrix Reduced 1 0" << std::endl;
      //      MatView(*(T.matrix_blocks[1][0]), PETSC_VIEWER_STDOUT_SELF);
      //      std::cout << "Matrix Reduced 1 1" << std::endl;
      //      MatView(*(T.matrix_blocks[1][1]), PETSC_VIEWER_STDOUT_SELF);

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

      F.clear();
      T.clear();

      // Assemble the matrix (mandatory after setting values)
      MatAssemblyBegin(romT, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(romT, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(romF, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(romF, MAT_FINAL_ASSEMBLY);

      //      std::cout << "romT " << std::endl;
      //      MatView(romT, PETSC_VIEWER_STDOUT_SELF);
      //      std::cout << "romF " << std::endl;
      //      MatView(romF, PETSC_VIEWER_STDOUT_SELF);
    }
    else
    {
      // Assemble block matrices
      Tred.assemble_full_matrices(static_problem.materials,
        static_problem.boundary_conditions,
        static_problem.albedo_factors, points_per_block);
      Fred.assemble_full_matrices(static_problem.materials, points_per_block);

      //      std::cout << "Matrix Reduced  0 0 " << std::endl;
      //      Tred.matrix_sp_blocks[0][0].print(std::cout);
      //      std::cout << "Matrix Reduced  0 1 " << std::endl;
      //      Tred.matrix_sp_blocks[0][1].print(std::cout);
      //      std::cout << "Matrix Reduced  1 0 " << std::endl;
      //      Tred.matrix_sp_blocks[1][0].print(std::cout);
      //      std::cout << "Matrix Reduced  1 1 " << std::endl;
      //      Tred.matrix_sp_blocks[1][1].print(std::cout);

      //std::cout << "snap_basis_full"<< std::endl;
      //for (unsigned int dr = 0; dr < dim_rom; dr++)
      //  snap_basis_full[dr].print(std::cout);
      //      std::cout << "snap_basis_red" << std::endl;
      //      for (unsigned int dr = 0; dr < dim_rom; dr++)
      //        snap_basis_red[dr].print(std::cout);

      BlockVector<double> auxT(n_groups, n_LUPOD_points);
      BlockVector<double> auxF(n_groups, n_LUPOD_points);
      for (unsigned int b2 = 0; b2 < dim_rom; b2++)
      {
        Tred.vmult(auxT, snap_basis_full[b2]);
        Fred.vmult(auxF, snap_basis_full[b2]);

        //std::cout << "b " << b2 << std::endl;
        //std::cout << "auxT_b0 " << auxT.block(0).norm_sqr() << std::endl;
        //std::cout << "auxT_b1 " << auxT.block(1).norm_sqr() << std::endl;
        //std::cout << "auxF_b0 " << auxF.block(0).norm_sqr() << std::endl;
        //std::cout << "auxF_b1 " << auxF.block(1).norm_sqr() << std::endl;
        //std::cout << std::endl;

        //std::cout << "b " << b2 << std::endl;
        //std::cout << "snap_basis_red_b0 " << snap_basis_red[b2].block(0).norm_sqr() << std::endl;
        //std::cout << "snap_basis_red_b1 " << snap_basis_red[b2].block(1).norm_sqr() << std::endl;
        //std::cout << std::endl;

        for (unsigned int b1 = 0; b1 < dim_rom; b1++)
        {
          MatSetValue(romT, b1, b2, auxT * snap_basis_red[b1], INSERT_VALUES);
          MatSetValue(romF, b1, b2, auxF * snap_basis_red[b1], INSERT_VALUES);
        }
      }

      // Assemble the matrix (mandatory after setting values)
      MatAssemblyBegin(romT, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(romT, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(romF, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(romF, MAT_FINAL_ASSEMBLY);

      //      std::cout << "romT " << std::endl;
      //      MatView(romT, PETSC_VIEWER_STDOUT_SELF);
      //      std::cout << "romF " << std::endl;
      //      MatView(romF, PETSC_VIEWER_STDOUT_SELF);
    }
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
    //EPSSetTolerances(eps, static_problem.tol_eps, 1e5);
    //EPSSetType(eps, EPSARPACK);
    EPSSolve(eps);

    PetscInt nconv;
    EPSGetConverged(eps, &nconv);
    AssertRelease(nconv >= 1, "REDUCED EIGENVALUE PROBLEM NOT CONVERGED");

    double eigr, eigi;
    EPSGetEigenpair(eps, 0, &eigr, &eigi, phi_red, NULL);

    if (eigi > 1e-9)
      std::cerr << "ERROR: The largest  eigenvalue of the reactor is complex:  Keff = "
                << std::scientific
                << eigr << " + " << eigi << " i." << std::endl;

    if (LUPOD_type == "POD")
    {
      for (unsigned int j = 0; j < n_dofs * n_groups; j++)
        for (unsigned int k = 0; k < dim_rom; k++)
        {
          phi[j] += snap_basis[k][j] * phi_red[k];
        }
    }
    else // LUPOD and others
    {
      for (unsigned int j = 0; j < n_dofs * n_groups; j++)
        for (unsigned int k = 0; k < dim_rom; k++)
        {
          phi[j] += snap_basis_full[k][j] * phi_red[k];
        }
    }
    phi.compress(VectorOperation::add);
    EPSDestroy(&eps);

    return eigr;
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void ROMStatic<dim, n_fe_degree>::run ()
  {
    cout << std::endl;
    cout << "----------------------------------------------------------"
         << std::endl;

    cout << "   Type of perturbation: " << type_perturbation << std::endl
         << std::endl;

    cout << std::fixed
         << "   Compute POD basis...                CPU Time = "
         << timer.wall_time() << " s." << std::endl;

    if (rom_group_wise == "Monolithic")
    {
      if (LUPOD_type == "LUPOD")
      {
        cout << "   LUPOD ----- MONOLITHIC" << std::endl;
        compute_LUPOD_basis_monolithic(snapshots, epsilon_M, epsilon_N,
           points_per_block, dim_rom, snap_basis_full, snap_basis_red);
      }
      else if (LUPOD_type == "LUPOD_ext")
      {
        cout << "   LUPOD EXTENDED  ----- MONOLITHIC" << std::endl;
        cout << "   N_LUPOD_POINTS: " << n_LUPOD_points << std::endl;
        compute_LUPODext_basis_monolithic(snapshots, epsilon_M, M_req, n_LUPOD_points,
           points_per_block, dim_rom, snap_basis_full, snap_basis_red);
      }
      else if (LUPOD_type == "POD")
      {
        cout << "   POD ----- MONOLITHIC" << std::endl;
        compute_POD_basis_monolithic(snapshots, epsilon_M, M_req, dim_rom,
          snap_basis);
      }
    }
    else if (rom_group_wise == "Group_Wise")
    {
      if (LUPOD_type == "POD")
      {
        cout << "   POD ---- GROUP WISE " << std::endl;
        compute_POD_basis_group_wise(snapshots, epsilon_M, M_req, dim_rom, snap_basis);
      }
      else if (LUPOD_type == "LUPOD")
      {
        cout << "   LUPOD ----- GROUP WISE" << std::endl;
        compute_LUPOD_basis_group_wise(snapshots, epsilon_M, epsilon_N,
           points_per_block, dim_rom, snap_basis_full, snap_basis_red);
      }
      else if (LUPOD_type == "LUPOD_ext")
      {
        cout << "   LUPOD EXTENDED----- GROUP WISE" << std::endl;
        cout << "   N_LUPOD_POINTS PER BLOCK: " << n_LUPOD_points << std::endl;
        compute_LUPODext_basis_group_wise(snapshots, epsilon_M, M_req, n_LUPOD_points,
           points_per_block, dim_rom, snap_basis_full, snap_basis_red);
      }
      else if (LUPOD_type == "Random")
      {
        cout << "   LUPOD RANDOM ----- GROUP WISE" << std::endl;
        cout << "   N_RANDOM_POINTS: " << n_LUPOD_points << std::endl;
        compute_random_points_group_wise(snapshots, n_LUPOD_points,
           points_per_block, dim_rom, snap_basis_full, snap_basis_red);
      }
      else if (LUPOD_type == "FEM1")
      {
        cout << "   LUPOD FEM 1 ----- GROUP WISE" << std::endl;
        compute_points_FEM1<dim>(dof_handler, snapshots, n_LUPOD_points,
           points_per_block, dim_rom, snap_basis_full, snap_basis_red);
      }
    }
    else
      AssertRelease(false, "rom_group_wise must be Monolithic or Group_Wise");

    double time_get_snap = timer.wall_time();
    cout << "   Time to get and process snapshots: " << time_get_snap << " s."
         << std::endl;

    //    -------------------------------------------------------------------------    //
    // Create Tests
    std::vector<Vector<double> > xs_test;
    std::mt19937 gen(seed); // Random number generator for uniform distribution [0.0,1.0]
    get_pertubation_random_XS(n_test, gen, xs_test);

    //---------------------------------------------------------------------------//
    //------------------------------- ROM ---------------------------------------//
    Vector<double> eig_rom(n_test);
    phi_rom.resize(n_test, PETScWrappers::MPI::BlockVector(local_dofs_vector, comm));
    for (unsigned int t = 0; t < n_test; t++)
      phi_rom[t].reinit(local_dofs_vector, comm);

    prepare_ROM();

    double timer_eig_starts, time_eig = 0.0;
    for (unsigned int t = 0; t < n_test; t++)
    {
      verbose_cout << "Update XS... " << std::flush;
      update_xs(static_problem.materials, xs_test[t]);
      verbose_cout << " Done!" << std::endl;

      verbose_cout << "assemble_ROM_matrices... " << std::flush;
      assemble_ROM_matrices();
      verbose_cout << " Done!" << std::endl;
      verbose_cout << "   Solve the ROM system...                 CPU Time = "
                   << timer.wall_time()
                   << " s." << std::endl;

      timer_eig_starts = timer.wall_time();
      eig_rom[t] = solve_eigenvalue(phi_rom[t]);
      time_eig += timer.wall_time() - timer_eig_starts;
      // Postprocess
      phi_rom[t].compress(VectorOperation::insert);
      phi_rom[t] *= phi_rom[t].mean_value();
      phi_rom[t] /= phi_rom[t].linfty_norm();

      cout << "      Test ROM " << t << ", Eigenvalue: " << std::setprecision(5)
           << eig_rom[t]
           << std::endl;
      verbose_cout << " Done!" << std::endl;
    }
    double time_rom = timer.wall_time() - time_get_snap;
    cout << "   Time ROM: " << time_rom << " s." << std::endl;
    cout << "   ------------------------------------   " << std::endl;
    cout << std::endl;

    MatDestroy(&romT);
    MatDestroy(&romF);

    //---------------------------------------------------------------------------//
    //------------------------------- FOM ---------------------------------------//
    Vector<double> eig_fom(n_test);
    phi_fom.resize(n_test, PETScWrappers::MPI::BlockVector(local_dofs_vector, comm));
    for (unsigned int t = 0; t < n_test; t++)
      phi_fom[t].reinit(local_dofs_vector, comm);

    double timer_fom_eig_start, time_fom_eig = 0.0;
    for (unsigned int t = 0; t < n_test; t++)
    {
      verbose_cout << "Update XS... " << std::flush;
      update_xs(static_problem.materials, xs_test[t]);
      verbose_cout << " Done!" << std::endl;

      //static_problem.tol_eps =
      //static_problem.tol_ksp =
      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.matrixfree_type = full_allocated;
      static_problem.assemble_system_lambda();
      timer_fom_eig_start = timer.wall_time();
      static_problem.solve_eps();
      time_fom_eig += timer.wall_time() - timer_fom_eig_start;

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

    double time_fom = timer.wall_time() - time_get_snap - time_rom;
    cout << "   Time FOM: " << time_fom << " s." << std::endl;
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
        rms_phi[t] += (phi_rom[t][i] - phi_fom[t][i]) * (phi_rom[t][i] - phi_fom[t][i]);
      rms_phi[t] = 100 * sqrt(rms_phi[t] / n_dofs * n_groups);
    }
    cout << std::scientific;
    cout << "   Mean Delta Keff = " << delta_keff.mean_value() << " pcm." << std::endl;
    cout << "   Max  Delta Keff = " << delta_keff.linfty_norm() << " pcm." << std::endl;
    cout << "   Mean RMS Phi = " << rms_phi.mean_value() << " %." << std::endl;
    cout << "   Max  RMS Phi = " << rms_phi.linfty_norm() << " %." << std::endl;
    cout << std::fixed << setprecision(2);
    cout << "   Speed Up = " << time_fom / time_rom << "\n";
    cout << "   Speed Up Solver = " << time_fom_eig / time_eig << "\n";
    cout << std::endl;
    cout << "   Finished in " << timer.wall_time() << " s." << std::endl;

    // WRITE OUTPUT
    // Create and erase the content of output file
    std::ofstream out(out_file.c_str(), std::ios::out);
    // Print the eigenvalues in the outFile
    print_logo(out);
    out << "\n";
    out << " ----- ROM STATIC PROBLEM -----" << "\n";
    out << "\n";
    out << "Problem File: " << static_problem.input_file << "\n";
    out << "ROM_Type_Snapshots: " << type_snapshots << "\n";
    out << "Type_Perturbation: " << type_perturbation << "\n";
    out << "XS_Perturbation_Fraction: " << perturbation_frac << "\n";
    out << "ROM_Group_Wise: " << rom_group_wise << "\n";

    //
    out << "N_Snapshots: " << n_snap << "\n";
    out << "ROM_DIM: " << dim_rom << "\n";
    out << "LUPOD_type: " << LUPOD_type << "\n";
    out << "N_LUPOD_Points: " << n_LUPOD_points << "\n";
    out << std::scientific;
    out << "Epsilon_M: " << epsilon_M << "\n";
    out << "Epsilon_N: " << epsilon_N << "\n";
    out << "\n";
    //
    out << "N_test: " << n_test << "\n";
    out << "Mean Delta Keff (pcm): " << delta_keff.mean_value() << std::endl;
    out << "Max  Delta Keff (pcm): " << delta_keff.linfty_norm() << std::endl;
    out << "Mean RMS Phi (%): " << rms_phi.mean_value() << std::endl;
    out << "Max  RMS Phi (%): " << rms_phi.linfty_norm() << std::endl;

    // Times
    out << "CPU Time Snapshots (s): " << time_get_snap << "\n";
    out << "CPU Time FOM (s): " << time_fom << "\n";
    out << "CPU Time ROM (s): " << time_rom << "\n";
    out << "Speed Up: " << time_fom / time_rom << "\n";
    out << "Speed Up Solver: " << time_fom_eig / time_eig << "\n";
    out << "\n";
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

