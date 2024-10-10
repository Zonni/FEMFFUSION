/*
 * rom_kinetics.cc
 *
 *  Created on: 28 feb 2024
 *      Author: amanda
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

#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscts.h>
#include <petscis.h>
#include <petscmat.h>
#include <slepcsvd.h>
#include "slepcbv.h"



#include "../include/rom_kinetics.h"
#include "../include/utils.h"
#include "../include/materials.h"
#include "../include/perturbation.h"
#include "../include/eps_solvers/eps_solver.h"
#include "../include/static_diffusion.h"
#include "../include/printing.h"
#include "../include/preconditioner.h"
#include "../include/ihs.h"

#include <string>
#include <math.h>

using namespace dealii;

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  ROMKinetics<dim, n_fe_degree>::ROMKinetics (
    ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    const bool verbose,
    const bool silent,
    const bool run) :
      comm(MPI_COMM_WORLD),
      n_mpi_processes(
        Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(
        Utilities::MPI::this_mpi_process(comm)),
      n_local_cells(
        numbers::invalid_unsigned_int),
      verbose_cout(std::cout,
        verbose and this_mpi_process == 0),
      cout(std::cout,
        !silent and this_mpi_process == 0),
      n_groups(static_problem.n_groups),
      n_dofs(static_problem.n_dofs),
	  n_prec(static_problem.materials.get_n_precursors()),
      dof_handler(static_problem.dof_handler),
      constraints(static_problem.constraints),
      boundary_conditions(static_problem.boundary_conditions),
      n_assemblies(static_problem.n_assemblies),
      materials(static_problem.materials),
      perturbation(static_problem.perturbation),
      V(comm, dof_handler, constraints),
      F(comm, dof_handler, constraints),
      L(comm, dof_handler, constraints),
	  XBF(n_prec,SpectraBetaFission<dim, n_fe_degree>(comm, dof_handler, constraints)),
      assem_per_dim(materials.assem_per_dim)
  {


    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
      local_dofs_vector[g] = locally_owned_dofs;

    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Out parameters
    out_file = static_problem.out_file;
    out_flag = static_problem.out_flag;
    out_interval = prm.get_integer("Out_Interval");
    n_out_ref = static_problem.n_out_ref;
    print_rhs = false;


    std::string geo_type = prm.get("Geometry_Type");
    listen_to_material_id = (geo_type == "Composed");

    matrixfree_type_time = non_diagonal;
    matrixfree_type = static_problem.matrixfree_type;

    // Material parameters
    albedo_factors = static_problem.albedo_factors;
    prec_flag = prm.get_bool("Precursors_Flag");

    // 	Time parameters
    type_perturbation = prm.get("Type_Perturbation");
    t_end = prm.get_double("Time_End");
    time_init_upd = 0.0;

    // Solver parameters
    init_delta_t = prm.get_double("Time_Delta");
    tol_time_ksp = static_problem.tol_ksp;


    step = 0;
    print_step = 1;
    sim_time = 0.0;
    power_total = 0.0;
    aux_power = 1.0;
    step_rom=0;

    verbose_cout << "Initialize the perturbation class" << std::endl;
    perturbation.init_transient();


    // Initialize phi //TODO remove
    phi.reinit(local_dofs_vector, comm);
    static_problem.phi[0].compress(VectorOperation::insert);
    phi = static_problem.phi[0];
    phi_critic = static_problem.phi[0]; // This vector is necessary to compute the noise


    // Reinit ROM data
    init_n_snap=prm.get_integer("N_Snapshots");
	// Initialize snapshots
	get_uint_from_options("-n_snap", init_n_snap);
    n_snap=init_n_snap;
    n_sets_snap=1;

    get_parameters_from_command_line();

    // Get type of snapshots
    parse_vector(prm.get("ROM_Type_Snapshots"), type_snapshots);
    n_sets_snap=type_snapshots.size();
    for (unsigned int ns=0; ns<n_sets_snap; ns++)
    	std::cout<<"type: "<<ns<<": "<<type_snapshots[ns]<<std::endl;



    // Get snapshots
    snapshots.resize(n_sets_snap);
    for (unsigned int ns=0; ns<type_snapshots.size(); ns++){
    	get_snapshots(prm, static_problem, snapshots[ns], type_snapshots[ns]);
    	std::cout<<"snapshots size: "<<ns<<": "<<snapshots[ns].size()<<std::endl;
    }

    if (n_sets_snap>1)
    	parse_vector(prm.get("ROM_Time_Break_Snapshots"), time_intervals_snapshots,n_sets_snap-1);

    time_intervals_snapshots.push_back(t_end);

    // Remove precursors
    if (prec_flag == false)
    {
      n_prec = 0;
      materials.remove_precursors();
    }

    // Initialization
    time_vect.reserve(t_end / init_delta_t);
    power_vector.reserve(t_end / init_delta_t);



    if (run)
    this->run();
  }



/**
 *
 *
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::init_time_computation ()
  {

    verbose_cout << "Assemble XBF..." << std::endl;
    for (unsigned int p=0; p<n_prec; p++)
    XBF[p].reinit(materials, full_matrixfree,p, listen_to_material_id);

    VecCreateSeq(PETSC_COMM_SELF, (n_prec + 1) * dim_rom, &coeffs_n);
    PetscScalar *N;
    VecGetArray(coeffs_n, &N);

      // Compute N[0]=(n_1,n_2,...,n_m)
      for (unsigned int nb = 0; nb < dim_rom; nb++){
        N[nb] = snap_basis[nb]* phi/(snap_basis[nb]*snap_basis[nb]);
      }


      // Compute c[0]=(c11,...,cq1,...,c1k,...,cqk)
      std::vector<std::vector<double>> cjk(dim_rom,
        std::vector<double>(n_prec));
      PETScWrappers::MPI::BlockVector aux(snap_basis[0]);


    if (step_rom == 0)
	{
		for (unsigned int nb = 0; nb < dim_rom; nb++)
		{
			for (unsigned int k = 0; k < n_prec; k++)
			{
				AssertRelease(materials.get_delayed_decay_constant(0, k) > 0.0,
						"Delayed decay constant must be greater than 0.");
				cjk[nb][k] = 1.0 / (materials.get_delayed_decay_constant(0, k))
						* XBF[k].vmult_dot(snap_basis[nb], phi);

			}
		}
	}
	else
	{
		std::vector<std::vector<double>> ajm(dim_rom,
		        std::vector<double>(snap_basis_old.size()));

		for (unsigned int j= 0; j < dim_rom; j++){
			for (unsigned int m= 0; m < snap_basis_old.size(); m++){
			ajm[j][m]=(snap_basis_old[m]* snap_basis[j])/(snap_basis_old[m]*snap_basis_old[m]);


			std::cout<<"a_"<<j<<m<<", value:"<<ajm[j][m]<<std::endl;
			}

		}

		// Test basis
		PETScWrappers::MPI::BlockVector aux(snap_basis[0]);
		PETScWrappers::MPI::BlockVector erroraux(snap_basis[0]);

		for (unsigned int j= 0; j < dim_rom; j++){
			aux=0.0;
			erroraux=0.0;
			for (unsigned int m= 0; m < snap_basis_old.size(); m++){
			aux.add(ajm[j][m],snap_basis_old[m]);
			}

			erroraux.add(1.0,aux, -1.,snap_basis[j]);
			std::cout<<"erroraux: "<<erroraux.l2_norm()<<std::endl;

		}


		for (unsigned int nb = 0; nb < dim_rom; nb++)
		{

			for (unsigned int k = 0; k < n_prec; k++)
			{
				AssertRelease(materials.get_delayed_decay_constant(0, k) > 0.0,
						"Delayed decay constant must be greater than 0.");
				cjk[nb][k]*=0.0;
				for (unsigned int m= 0; m < snap_basis_old.size(); m++){
				cjk[nb][k] += ajm[nb][m]*cjk_old[m][k];
				}

			}
		}
	}


      for (unsigned int k = 0; k < n_prec; k++)
        for (unsigned int nb = 0; nb < dim_rom; nb++)
            N[(k+1) * dim_rom + nb] = cjk[nb][k];



    VecRestoreArray(coeffs_n, &N);




    for (unsigned int k = 0; k < n_prec; k++)
      XBF[k].clear();


  }

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots(
		ParameterHandler &prm,
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
		std::string t_snap)
{



	if (t_snap == "modes")
		get_snapshots_modes(static_problem, _snapshots);
	else if (t_snap == "bank1")
	{
		get_snapshots_bar(static_problem, _snapshots, 1);
	}
	else if (t_snap == "bank2")
	{
		get_snapshots_bar(static_problem, _snapshots, 2);
	}
	else if (t_snap == "bank12")
	{
		get_snapshots_bar(static_problem, _snapshots, 1);
		get_snapshots_bar(static_problem, _snapshots, 2);
	}
	else if (t_snap == "bank6")
	{
		get_snapshots_bar(static_problem, _snapshots, 6);
	}
	else if (t_snap == "bank135")
	{
		get_snapshots_bar(static_problem, _snapshots, 1);
		get_snapshots_bar(static_problem, _snapshots, 3);
		get_snapshots_bar(static_problem, _snapshots, 5);
	}
	else if (t_snap == "bank123")
	{
		get_snapshots_bar(static_problem, _snapshots, 1);
		get_snapshots_bar(static_problem, _snapshots, 2);
		get_snapshots_bar(static_problem, _snapshots, 3);
	}
	else if (t_snap == "bank1356")
	{
		get_snapshots_bar(static_problem, _snapshots, 1);
		get_snapshots_bar(static_problem, _snapshots, 3);
		get_snapshots_bar(static_problem, _snapshots, 5);
		get_snapshots_bar(static_problem, _snapshots, 6);
	}
	else if (t_snap == "bank123456")
	{
		get_snapshots_bar(static_problem, _snapshots, 1);
		get_snapshots_bar(static_problem, _snapshots, 2);
		get_snapshots_bar(static_problem, _snapshots, 3);
		get_snapshots_bar(static_problem, _snapshots, 4);
		get_snapshots_bar(static_problem, _snapshots, 5);
		get_snapshots_bar(static_problem, _snapshots, 6);
	}
	else if (t_snap == "bank12_ihs")
	{
		get_snapshots_bar_ihs(static_problem, _snapshots);
	}
	else if (t_snap == "bank12_time")
	{
		get_snapshots_bar_time_variation(static_problem, _snapshots);
	}
	else if (t_snap == "ramp_time")
	{
		get_snapshots_ramp_time(static_problem, _snapshots);
	}
	else if (t_snap == "ramp_time_dyn")
	{
		get_snapshots_ramp_time_dyn(prm,static_problem, _snapshots);
	}
	else if (t_snap == "ramp_time_mat1")
	{
		get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[0],0);
	}
	else if (t_snap == "ramp_time_mat2")
	{
		get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[1],1);
	}
	else if (t_snap == "ramp_time_mat12")
	{
		get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[0],0);
		get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[1],1);
	}
	else if (t_snap == "ramp_time_change_mat1")
	{
		get_snapshots_ramp_change(prm,static_problem, _snapshots, perturbation.mat_changing[0],0);
	}
	else if (t_snap == "ramp_time_change_mat2")
	{
		get_snapshots_ramp_change(prm,static_problem, _snapshots, perturbation.mat_changing[1],1);
	}
	else if (t_snap == "ramp_time_change_mat12")
	{
		get_snapshots_ramp_change(prm,static_problem, _snapshots, perturbation.mat_changing[0],0);
		get_snapshots_ramp_change(prm,static_problem, _snapshots, perturbation.mat_changing[1],1);
	}
	else
	{
		AssertRelease(false,"Invalid type of Snapshots");
	}

}

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_modes(
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
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_bar(
		StaticDiffusion<dim, n_fe_degree> &static_problem,std::vector<PETScWrappers::MPI::BlockVector> &_snapshots, unsigned int bar_bank)
{


	unsigned int min_s_snap= _snapshots.size();
	unsigned int max_s_snap= _snapshots.size()+ init_n_snap;


	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	double bars_bottom_pos = 0.0;
	AssertRelease(init_n_snap > 1, "The number of snapshots must be greater than 1");
	double step_bar = (perturbation.bars_top_pos - bars_bottom_pos)
			/ (init_n_snap - 1);

	unsigned int bar_material = perturbation.bar_materials[bar_bank - 1];

	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
	{

		unsigned int bar_pos_z = (ns-min_s_snap) * step_bar + bars_bottom_pos;
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
		cout << "Step " << ns << ", Bar Position: " << bar_pos_z
				<< ", Eigenvalue: " << static_problem.eigenvalues[0]
				<< std::endl;
	}


	n_snap=max_s_snap;


}

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_bar_ihs(
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
{

	int n_bars = 2;
	// IHS
	int duplication = 5;
	int seed = 17;
	int *x;
	std::vector<std::vector<double>> sample_pos(n_bars,
			std::vector<double>(init_n_snap));

	//
	//  Get the points.
	//
	x = ihs(n_bars, init_n_snap, duplication, seed);
	for (int d = 0; d < n_bars; d++)
		for (unsigned int i = 0; i < init_n_snap; i++){
			sample_pos[d][i] = static_cast<double>(x[n_bars * i + d]) / init_n_snap;
		}

	delete[] x;

	// Compute Snapshots
	unsigned int min_s_snap = _snapshots.size();
	unsigned int max_s_snap = _snapshots.size() + init_n_snap;

	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	// for each snapshot
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
	{
		//for each bar bank
		for (int nb = 1; nb < n_bars + 1; nb++)
		{

			double bars_bottom_pos = 0.0;
			unsigned int bar_material = perturbation.bar_materials[nb - 1];

			unsigned int bar_pos_z = sample_pos[nb - 1][ns]
					* (perturbation.bars_top_pos - bars_bottom_pos);

			std::cout<<"bank_bar"<<nb<<"bar_pos: "<<bar_pos_z<<std::endl;

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
		cout << "Step " << ns
				<< ", Eigenvalue: " << static_problem.eigenvalues[0]
				<< std::endl;

	}
	n_snap = max_s_snap;

}

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_bar_time_variation(
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots){

	AssertRelease(init_n_snap > 1, "The number of snapshots must be greater than 1");
	double sample_time;

	// Compute Snapshots
	unsigned int min_s_snap = _snapshots.size();
	unsigned int max_s_snap = _snapshots.size() + init_n_snap;

	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	// for each snapshot
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
	{
		//for each bar bank
		sample_time=(t_end-0.0)/(init_n_snap-1)*ns;
		perturbation.move_bars(sample_time);

		static_problem.cout.set_condition(false);
		static_problem.show_eps_convergence = false;
		static_problem.assemble_system_lambda();
		static_problem.solve_eps();
		static_problem.phi[0].compress(VectorOperation::insert);
		_snapshots[ns] = static_problem.phi[0];
		cout << "Step " << ns
				<< ", Eigenvalue: " << static_problem.eigenvalues[0]
				<< std::endl;
	}
	n_snap = max_s_snap;


}

/**
 *
 *
 *
 */
//template<int dim, int n_fe_degree>
//void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp(
//		StaticDiffusion<dim, n_fe_degree> &static_problem,
//		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
//		unsigned int mat_chan, double num_mat)
//{
//
//	double sample_time, sample_time_2;
//
//	unsigned int min_s_snap = _snapshots.size();
//	unsigned int max_s_snap = _snapshots.size() + init_n_snap;
//
//	_snapshots.resize(max_s_snap);
//	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
//		_snapshots[ns].reinit(local_dofs_vector, comm);
//
//	std::vector<int> mat_chan_vector;
//
//	std::cout<<"-------------MAT CHAN"<<mat_chan<<std::endl;
//
//	// for each snapshot
//	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
//	{
//		//for each bar bank
//		sample_time = (t_end - 0.0) / (init_n_snap - 1) * (ns-min_s_snap);
//
//
//		std::cout<<"sample_time: "<<sample_time<<std::endl;
//
//
//		if (sample_time < perturbation.cut_time[num_mat]
//													and std::abs(perturbation.slope_up[num_mat]) > 0 )
//		{
//			std::cout<<"slope_up"<<std::endl;
//			mat_chan_vector.clear();
//			mat_chan_vector.push_back(mat_chan);
//			//sample_time_2 = (perturbation.cut_time[num_mat] - 0.0) / (init_n_snap - 1) * ns;
//			perturbation.modify_xsec(sample_time, mat_chan_vector);
//
//		}
//		else if (sample_time > perturbation.cut_time[num_mat]
//														and std::abs(perturbation.slope_down[num_mat])>0)
//		{
//			std::cout<<"--slope_down"<<std::endl;
//			mat_chan_vector.clear();
//			mat_chan_vector.push_back(-1);
//			mat_chan_vector.push_back(mat_chan);
//			//sample_time_2 = (t_end-perturbation.cut_time[num_mat]) / (init_n_snap - 1) * ns+perturbation.cut_time[num_mat];
//			perturbation.modify_xsec(sample_time, mat_chan_vector);
//
//		}
//		else{
//			continue;
//		}
//
//		static_problem.cout.set_condition(false);
//		static_problem.show_eps_convergence = false;
//		static_problem.assemble_system_lambda();
//		static_problem.solve_eps();
//		static_problem.phi[0].compress(VectorOperation::insert);
//		_snapshots[ns] = static_problem.phi[0];
//		cout << "Step " << ns << ", Eigenvalue: "
//				<< static_problem.eigenvalues[0] << std::endl;
//	}
//
//
//	n_snap = max_s_snap;
//
//}

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp(
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
		unsigned int mat_chan, double num_mat)
{

	double sample_time, sample_time_2;

	unsigned int min_s_snap = _snapshots.size();
	unsigned int max_s_snap = _snapshots.size() + init_n_snap;

	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	std::vector<int> mat_chan_vector;

	std::cout<<"-------------MAT CHAN"<<mat_chan<<std::endl;

	// for each snapshot
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
	{


		if (std::abs(perturbation.slope_up[num_mat]) > 0 )
		{
			std::cout<<"slope_up"<<std::endl;
			mat_chan_vector.clear();
			mat_chan_vector.push_back(mat_chan);
			sample_time_2 = (perturbation.cut_time[num_mat] - 0.0) / (init_n_snap - 1) * (ns-min_s_snap);
			perturbation.modify_xsec(sample_time_2, mat_chan_vector);

		}
		else if (std::abs(perturbation.slope_down[num_mat])>0)
		{
			std::cout<<"--slope_down"<<std::endl;
			mat_chan_vector.clear();
			mat_chan_vector.push_back(-1);
			mat_chan_vector.push_back(mat_chan);
			sample_time_2 = (t_end-perturbation.cut_time[num_mat]) / (init_n_snap - 1) * (ns-min_s_snap) + perturbation.cut_time[num_mat];
			perturbation.modify_xsec(sample_time_2, mat_chan_vector);
			//TODO
			mat_chan_vector.clear();
			mat_chan_vector.push_back(1);
			mat_chan_vector.push_back(-1);
			perturbation.modify_xsec(perturbation.cut_time[0], mat_chan_vector);

		}
		else{
			continue;
		}

		static_problem.cout.set_condition(false);
		static_problem.show_eps_convergence = false;
		static_problem.assemble_system_lambda();
		static_problem.solve_eps();
		static_problem.phi[0].compress(VectorOperation::insert);
		_snapshots[ns] = static_problem.phi[0];
		cout << "Step " << ns << ", Eigenvalue: "
				<< static_problem.eigenvalues[0] << std::endl;
	}


	n_snap = max_s_snap;

}

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp_change(ParameterHandler &prm,
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
		unsigned int mat_chan, double num_mat)
{

	double sample_time, sample_time_2;

	unsigned int min_s_snap = _snapshots.size();
	unsigned int max_s_snap = _snapshots.size() + init_n_snap;

	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	std::vector<int> mat_chan_vector;

	std::vector<double> ROM_Slope_Up;
	std::vector<double> ROM_Slope_Down;
	std::vector<double> ROM_Cut_Time;

	parse_vector(prm.get("ROM_Slope_Up"), ROM_Slope_Up);
	parse_vector(prm.get("ROM_Slope_Down"), ROM_Slope_Down);
	parse_vector(prm.get("ROM_Cut_Time"), ROM_Cut_Time);



	std::cout<<"-------------MAT CHAN"<<mat_chan<<std::endl;

	// for each snapshot
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
	{


		if (std::abs(ROM_Slope_Up[num_mat]) > 0 )
		{
			std::cout<<"slope_up"<<std::endl;
			mat_chan_vector.clear();
			mat_chan_vector.push_back(mat_chan);
			sample_time_2 = (ROM_Cut_Time[num_mat] - 0.0) / (init_n_snap - 1) * (ns-min_s_snap);
			perturbation.modify_xsec(sample_time_2, mat_chan_vector);

		}
		else if (std::abs(ROM_Slope_Down[num_mat])>0)
		{
			std::cout<<"--slope_down"<<std::endl;
			mat_chan_vector.clear();
			mat_chan_vector.push_back(-1);
			mat_chan_vector.push_back(mat_chan);
			sample_time_2 = (t_end-ROM_Cut_Time[num_mat]) / (init_n_snap - 1) * (ns-min_s_snap) + ROM_Cut_Time[num_mat];
			perturbation.modify_xsec(sample_time_2, mat_chan_vector);
			//TODO
			mat_chan_vector.clear();
			mat_chan_vector.push_back(1);
			mat_chan_vector.push_back(-1);
			perturbation.modify_xsec(ROM_Cut_Time[0], mat_chan_vector);

		}
		else{
			continue;
		}

		static_problem.cout.set_condition(false);
		static_problem.show_eps_convergence = false;
		static_problem.assemble_system_lambda();
		static_problem.solve_eps();
		static_problem.phi[0].compress(VectorOperation::insert);
		_snapshots[ns] = static_problem.phi[0];
		cout << "Step " << ns << ", Eigenvalue: "
				<< static_problem.eigenvalues[0] << std::endl;
	}


	n_snap = max_s_snap;

}
/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp_time(
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots){

	AssertRelease(init_n_snap > 1, "The number of snapshots must be greater than 1");
	double sample_time;

	// Compute Snapshots
	unsigned int min_s_snap = _snapshots.size();
	unsigned int max_s_snap = _snapshots.size() + init_n_snap;

	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	// for each snapshot
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
	{
		//for each bar bank
		sample_time=(t_end-0.0)/(init_n_snap-1)*ns;
		std::cout<<"sample time: "<<sample_time<<std::endl;
		perturbation.apply_function_to_perturb(sample_time);
		std::cout<<"done!"<<std::endl;
		static_problem.cout.set_condition(false);
		static_problem.show_eps_convergence = false;
		static_problem.assemble_system_lambda();
		static_problem.solve_eps();
		static_problem.phi[0].compress(VectorOperation::insert);
		_snapshots[ns] = static_problem.phi[0];
		cout << "Step " << ns
				<< ", Eigenvalue: " << static_problem.eigenvalues[0]
				<< std::endl;
	}
	n_snap = max_s_snap;


}

/**
 *
 *
 *
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp_time_dyn(
		ParameterHandler &prm,
		StaticDiffusion<dim, n_fe_degree> &static_problem,
		std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
{

	AssertRelease(init_n_snap > 1,
			"The number of snapshots must be greater than 1");
	double sample_time_step;

	// Compute Snapshots
	unsigned int min_s_snap = _snapshots.size();
	unsigned int max_s_snap = _snapshots.size() + init_n_snap;

	_snapshots.resize(max_s_snap);
	for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
		_snapshots[ns].reinit(local_dofs_vector, comm);

	TimeNeutronDiffusion<dim, n_fe_degree> time_pro(prm, static_problem,
			false, true, false);


	//for each bar bank
	sample_time_step = (t_end - 0.0) / (init_n_snap - 1);
	std::cout << "sample time: " << sample_time_step << std::endl;
	time_pro.init_delta_t = sample_time_step;
	time_pro.delta_t.reserve(t_end / time_pro.init_delta_t);
	time_pro.delta_t[0]=init_delta_t;
	time_pro.init_time_computation();

	while (time_pro.t_end - time_pro.sim_time > -1e-12)
	{

		time_pro.update_xsec();
		time_pro.assemble_matrices();
		time_pro.compute_RHS();
		time_pro.solve_LHS();
		time_pro.phi.compress(VectorOperation::insert);
		_snapshots[time_pro.step] = time_pro.phi;
		cout << "Snapshot " << time_pro.step <<" at time: "<<time_pro.sim_time<< std::endl;

		time_pro.step++;
		time_pro.delta_t.push_back(time_pro.init_delta_t);
		time_pro.sim_time += time_pro.delta_t[time_pro.step];

	}

	n_snap = max_s_snap;

}

/**
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_parameters_from_command_line ()
  {

    // Booleans
    get_bool_from_options("-out_flag", out_flag);
    get_bool_from_options("-prec_flag", prec_flag);
    get_bool_from_options("-print_timefile", print_timefile);
    get_bool_from_options("-print_rhs", print_rhs);


    // Integers
    get_uint_from_options("-n_out_ref", n_out_ref);
    get_uint_from_options("-out_interval", out_interval);

    // Reals
    get_double_from_options("-init_delta_t", init_delta_t);
    get_double_from_options("-t_end", t_end);
    get_double_from_options("-tol_time_ksp", tol_time_ksp);

    // String
    get_string_from_options("-out_file", out_file);
    get_enum_from_options("-matrixfree_type_time", matrixfree_type_time);

  }

/*
 *
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::update_xsec ()
  {

    if (type_perturbation == "Flux_Distributed"
        or type_perturbation == "Single_Material"
        or type_perturbation == "Out_Of_Phase"
        or type_perturbation == "Ramp_Two_Mats")
    {
      verbose_cout << "Apply function to perturbed " << std::endl;
      perturbation.apply_function_to_perturb(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Rods")
    {
      verbose_cout << "Moving rods: time" << sim_time << std::endl;
      perturbation.move_bars(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "AECL")
    {
      verbose_cout << "Perturbed the AECL transient: " << std::endl;
      perturbation.move_th(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Step_Change_Material")
    {
      verbose_cout << "Perturbed the Step_Change_Material: " << std::endl;
      perturbation.step_change_material(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Mechanical_Vibration")
    {
      verbose_cout << "   move_vibrating... " << std::flush;
      perturbation.move_vibrating(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "C5G7-TD1.1")
    {
      verbose_cout << "Apply perturbation C5G7-TD1.1: " << std::endl;
      perturbation.apply_c5G7_perturb(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Read_XS_File")
    {
      verbose_cout << "   move_read_xs_file... " << std::flush;
      perturbation.move_read_xs_file(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else
    {
      AssertRelease(false, "Invalid type of perturbation");
    }
  }

/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::compute_pod_basis (std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {

	dim_rom=_snapshots.size();
    get_uint_from_options("-dim_rom", dim_rom);
    std::cout<<"dim_rom: "<<dim_rom<<std::endl;


	snap_basis.resize(dim_rom);
	for (unsigned int dr=0; dr<dim_rom; dr++){
	snap_basis[dr].reinit(local_dofs_vector, comm);
	snap_basis[dr]*=0.0;
	}

	// Create the matrix with the snapshot to apply the SVD
	Mat Mat_snap;
	PetscInt i_snap, i_sv;
	PetscInt *idm = new PetscInt[n_dofs*n_groups];
	PetscScalar *values_snap = new PetscScalar[n_dofs*n_groups];

	for (unsigned int j=0; j<n_dofs*n_groups; j++)
		idm[j]=j;

	MatCreateDense(comm, n_dofs*n_groups, dim_rom, n_dofs*n_groups, dim_rom, NULL, &Mat_snap);
	for (i_snap=0; i_snap<static_cast<int>(dim_rom); i_snap++){
		for (unsigned k=0; k<n_dofs*n_groups; k++)
			values_snap[k]=_snapshots[i_snap][k];
	MatSetValues(Mat_snap, n_dofs*n_groups, idm, 1, &i_snap, values_snap, INSERT_VALUES);
	}
	MatAssemblyBegin(Mat_snap, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Mat_snap, MAT_FINAL_ASSEMBLY);


	SVD svd;
	SVDCreate(comm,&svd);
	SVDSetOperators(svd,Mat_snap,NULL);
	SVDSetDimensions(svd,dim_rom,2*dim_rom,dim_rom);
	SVDSetProblemType(svd,SVD_STANDARD);
	SVDSetType(svd,SVDTRLANCZOS);

	SVDSetFromOptions(svd);
	SVDSolve(svd);

	PETScWrappers::MPI::Vector u;
	u.reinit(comm, n_dofs*n_groups,n_dofs*n_groups);
	std::vector<double> singular_values(dim_rom);
	for (i_sv=0; i_sv<static_cast<int>(dim_rom); i_sv++){
	SVDGetSingularTriplet(svd,i_sv,&(singular_values[i_sv]),u,NULL);
	copy_to_BlockVector(snap_basis[i_sv], u);
   	cout<<"Singular value "<<i_sv<<": "<<singular_values[i_sv]<<std::endl;
	}

	MatDestroy(&Mat_snap);
	SVDDestroy(&svd);
	u.clear();


  }


/*
 * @brief compute_pod_basis This function must be work for SLEPc 3.21.0
 */
//template <int dim, int n_fe_degree>
//  void ROMKinetics<dim, n_fe_degree>::compute_pod_basis ()
//  {
//
//	// Create a BV type from basis (BVInsertVec) and then apply BVSVDAndRank
//	unsigned int n_snap=snapshots.size(), n_sinvalues=3;
//
//	// Create the matrix with the snapshot to apply the SVD
//	BV Snap_basis;
//	BVCreate(comm, &Snap_basis);
//	BVSetSizes(Snap_basis,n_dofs*n_groups,n_dofs*n_groups,n_snap);
//	BVSetType(Snap_basis,BVVECS);
//
//	PetscInt i_snap;
//	Vec v;
//
//	for (i_snap=0; i_snap<static_cast<int>(n_snap); i_snap++){
//		BVGetColumn(Snap_basis,i_snap,&v);
//		copy_to_Vec(v,snapshots[i_snap]);
//		BVRestoreColumn(Snap_basis,i_snap,&v);
//	}
//
//	BVView(Snap_basis,PETSC_VIEWER_STDOUT_WORLD );
//	PetscScalar A[n_sinvalues*n_snap*n_sinvalues*n_snap];
//	PetscReal sigma;
//	PetscInt rank;
//	// This function is not available in SLEPc 3.15
////	BVSVDAndRank(Snap_basis,n_sinvalues,n_snap,1E-09,BV_SVD_METHOD_QR,&A,&sigma,&rank);
//
//  }


/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::assemble_matrices ()
  {

    // Assemble the mass matrix
    // It is necessary to assemble all matrices even if the operators
    // do not change because the number of material (in rods) can be changed
//	verbose_cout << "Assemble V..." << std::endl;
	V.reinit(materials, full_matrixfree, listen_to_material_id);
	verbose_cout << "Assemble L..." << std::endl;
	L.reinit(materials, boundary_conditions, albedo_factors, full_matrixfree,
			listen_to_material_id);
	verbose_cout << "Assemble F..." << std::endl;
	F.reinit(materials, full_matrixfree, listen_to_material_id);
	verbose_cout << "Assemble XBF..." << std::endl;
	for (unsigned int p=0; p<n_prec; p++)
	XBF[p].reinit(materials, full_matrixfree,p, listen_to_material_id);


  }

/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::assemble_ROM_matrices ()
  {

	assemble_matrices ();

	FullMatrix<double> romV(dim_rom);
	rominvV = FullMatrix<double>(dim_rom);
	romL = FullMatrix<double>(dim_rom);
	romF = FullMatrix<double>(dim_rom);
	romXBF.resize(n_prec);
	for (unsigned int p=0; p<n_prec; p++)
		romXBF[p] = FullMatrix<double>(dim_rom);

	for (unsigned int b1 = 0; b1 < dim_rom; b1++)
	      for (unsigned int b2 = 0; b2 < dim_rom; b2++){
	    	  romV(b1, b2) = V.vmult_dot(snap_basis[b1], snap_basis[b2]);
	    	  romL(b1, b2) = L.vmult_dot(snap_basis[b1], snap_basis[b2]);
	    	  romF(b1, b2) = F.vmult_dot(snap_basis[b1], snap_basis[b2]);
	    	  for(unsigned int p=0; p<n_prec; p++)
	    		 romXBF[p](b1,b2)=XBF[p].vmult_dot(snap_basis[b1], snap_basis[b2]);
	      }

//	std::cout<<"romV:"<<std::endl;
//	romV.print(std::cout,10);


	rominvV.invert(romV);

//	std::cout<<"romL:"<<std::endl;
//	romL.print(std::cout,10);
//
//
//	std::cout<<"romF:"<<std::endl;


	V.clear();
	L.clear();
	F.clear();
	for(unsigned int p=0; p<n_prec; p++)
		XBF[p].clear();

  }




template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::solve_system_petsc ()
  {

    double init_time_step = 1e-2;
    unsigned int max_steps = 1e+5;

    // SOLVER PETSC implemented for a linear system
    TS ts;

    if (step == 0 and print_timefile)
    {
      filename_time = out_file;
      filename_time.erase(filename_time.end() - 4, filename_time.end());
      filename_time = filename_time + "_time.out";
      std::ofstream out(filename_time.c_str(), std::ios::out);
    }

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetProblemType(ts, TS_LINEAR);
    TSSetApplicationContext(ts, this);

    //   Tell the timestepper context where to compute solutions
    TSSetSolution(ts, coeffs_n);

	//       Provide the call-back for the nonlinear function we are
	//       evaluating. Thus whenever the timestepping routines need the
	//       function they will call this routine. Note the final argument
	//       is the application context used by the call-back functions
    TSSetRHSFunction(ts, NULL, FormFunctionROM_system<dim, n_fe_degree>,
        this);


    //   Form the initial guess for the problem
    //   (this can be done from a function FormInitialGuess)

    //   This indicates that we are using pseudo timestepping to
    //   find a steady state solution to the nonlinear problem.
    TSSetType(ts, TSBDF);//
    //TSSetType(ts, TSBEULER);

    //   Set the initial time to start at (this is arbitrary for
    //   steady state problems); and the initial timestep given above
    //    TSSetTimeStep(ts, 1e-5);
    TSSetTimeStep(ts, init_time_step);

    //   Set a large number of timesteps and final duration time
    //   to insure convergence to steady state.
    TSSetMaxSteps(ts, max_steps);
    //    if (update_modes == true)
    TSSetMaxTime(ts, delta_t);
    //    else
    //      TSSetMaxTime(ts, delta_t_petsc);
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_INTERPOLATE);
    TSSetPostStep(ts, PostStep<dim, n_fe_degree>);

    TSSetFromOptions(ts);
    TSSetUp(ts);

    //   Perform the solve. This is where the timestepping takes place.
    TSSolve(ts, coeffs_n);

    const PetscScalar *narray;

    VecGetArrayRead(coeffs_n, &narray);

    phi = 0.0;
    for (unsigned int dr = 0; dr < dim_rom; ++dr){
      phi.add(narray[dr], snap_basis[dr]);
    }


    // For updating the cjk
    cjk_old.resize(dim_rom,
      std::vector<double>(n_prec));
    for (unsigned int nb = 0; nb < dim_rom; nb++)
    	for (unsigned int k = 0; k < n_prec; k++)
    		cjk_old[nb][k]=narray[(k+1) * dim_rom + nb];


	if (std::abs(time_vect.back() - (time_init_upd + delta_t)) > 1e-12)
	{
	    postprocess_time_step();
	    cout << std::setprecision(4) << "Time:   " << time_init_upd + delta_t
	         << " ---> Power:  "
	         << power_total << std::endl;
		power_vector.push_back(power_total);
		time_vect.push_back(time_init_upd + delta_t);
	}
    aux_power = power_total;
    if (out_flag)
    output_results();


    //   Get the number of steps
    PetscInt its;
    TSGetStepNumber(ts, &its);
    cout
    << "                                                  its ts_solver: "
    << its
    << std::endl;

//    total_ts_its += its;

    TSDestroy(&ts);

  }


/*
 * @brief FormFunction definition for Petsc Solver
 */
template <int dim, int n_fe_degree>
  PetscErrorCode PostStep (TS ts)
  {


    void *ctx;
    TSGetApplicationContext(ts, &ctx);

    ROMKinetics<dim, n_fe_degree>* TSobject =
                                                (ROMKinetics<dim, n_fe_degree>*) ctx;
    PetscReal ts_time;
    const PetscScalar *n;
    Vec N;

    TSGetTime(ts, &ts_time);
    TSGetSolution(ts, &N);
    VecGetArrayRead(N, &n);

    // for updating
    double real_time = ts_time + TSobject->time_init_upd;
//    double real_time = ts_time + TSobject->t_init_upd;


    TSobject->phi = 0.0;

    for (unsigned int dr = 0; dr < TSobject->dim_rom; ++dr)
      TSobject->phi.add(n[dr], TSobject->snap_basis[dr]);

    TSobject->sim_time=real_time;
    TSobject->update_xsec();
    TSobject->postprocess_time_step();

//    if (real_time < TSobject->t_end_upd)
    if (real_time < TSobject->time_init_upd + TSobject->delta_t)
    {
      // Save data only if the change in the power is larger than 0.1%
      if (std::abs((TSobject->aux_power - TSobject->power_total) / TSobject->aux_power) > 1e-3)
      {
    	TSobject->step++;
        TSobject->power_vector.push_back(TSobject->power_total);
        TSobject->time_vect.push_back(real_time);
        TSobject->aux_power = TSobject->power_total;
        if (TSobject->out_flag)
        TSobject->output_results();
      }

      TSobject->cout << std::setprecision(4) << std::fixed << "Time:   "
      << real_time
      << " ---> Power:  " << TSobject->power_total
      << std::endl;
    }

    VecRestoreArrayRead(N, &n);

//    TSobject->tseg = transf_time;
//    VecCopy(N, TSobject->nseg);

//   VecDestroy(&N);

//    It does not work
//    Vec error_vec;
//    double err_norm;
//    VecDuplicate(N,&error_vec);
//    VecAssemblyBegin(error_vec);
//    VecAssemblyEnd(error_vec);
//
//    TSGetTimeError(ts,0,&error_vec);
//
//    VecNorm(error_vec,NORM_2,&err_norm);
//    VecDestroy(&error_vec);
//    PetscPrintf(PETSC_COMM_WORLD,"Estimated Error = %E.\n",err_norm);

    return 0;

  }



/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::postprocess_time_step ()
  {

    BlockVector<double> phi_serial;

    phi_serial.reinit(n_groups, n_dofs);
    for (unsigned int g = 0; g < n_groups; g++)
      phi_serial.block(g) = phi.block(g);

    if (this_mpi_process == 0)
    {

// Initialize all that  is needed to iterate over dofs and cells
      QGauss<dim> quadrature_formula(n_fe_degree + 1);

      FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
        update_values | update_volume_elements | update_JxW_values);

      unsigned int n_q_points = quadrature_formula.size();
      unsigned int n_cells_out = n_assemblies;
      double power_cell = 0;
      double phi_cell;
      double sigma_f;

      // Initialize and resize the vectors where it is stored the solution
      power_per_assembly.resize(1, std::vector<double>(n_cells_out, 0.0));
      std::vector<double> volume_per_assembly(n_cells_out, 0.0);
      std::vector<double> local_phi(n_q_points);
      std::vector<double> radial_power;

      // Initialize Values
      volume_per_assembly.assign(n_cells_out, 0.0);

      double volume = 0.0;
      double norm = 0.0;
      unsigned int mat_id;

      // Iterate over every cell
      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {

        fe_values.reinit(cell);
        power_cell = 0;
        mat_id = materials.get_material_id<dim>(cell);

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          sigma_f = materials.get_sigma_f(g, mat_id);

          fe_values.get_function_values(phi_serial.block(g), local_phi);

          phi_cell = 0.0;
          for (unsigned int q = 0; q < n_q_points; q++)
            phi_cell += local_phi[q] * fe_values.JxW(q);

          power_cell += sigma_f * phi_cell;
        }

        volume_per_assembly[cell->user_index()] += cell->measure();
        power_per_assembly[0][cell->user_index()] += power_cell;
        volume += cell->measure();
        norm +=  std::abs(power_cell);

      }

      power_total = norm / volume;


      norm /= volume;

	  // Normalize the values of the power and fluxes per cell
	  normalize_vector(power_per_assembly[0], norm);

      if (print_timefile and (step % print_step == 0))
      {

        std::ofstream out(filename_time.c_str(), std::ios::app);
        out.precision(9);
        out << "Time in step: " << " \n" << step << " " << sim_time
            << " \n";
        out.close();

        print_time_vect.push_back(sim_time);

        print_vector_in_file(power_per_assembly[0],
          filename_time,
          "Power per assembly \n",
          true);
      }

      if (dim == 3 and print_timefile)
      {

        // Calculate the axial power distribution
        std::vector<double> power_axial;
        std::vector<std::vector<double> > phi_axial;
        unsigned int n_planes = materials.assem_per_dim[2];

        unsigned int n_assemblies_per_plane = n_assemblies / n_planes;
        Assert(n_assemblies_per_plane > 0,
          ExcMessage("n_assemblies cannot be 0"));
        radial_power.resize(n_assemblies_per_plane);
        std::vector<double> volume_per_axial(n_assemblies_per_plane);

        for (unsigned int i = 0; i < n_assemblies_per_plane; i++)
        {
          double sum = 0.0;
          for (unsigned int pl = 0; pl < n_planes; pl++)
          {
            if (std::abs(power_per_assembly[0][n_assemblies_per_plane * pl + i]) > 1e-5)
            {
              sum += power_per_assembly[0][n_assemblies_per_plane * pl + i]
                     * volume_per_assembly[n_assemblies_per_plane * pl + i];
              volume_per_axial[i] += volume_per_assembly[n_assemblies_per_plane * pl + i];
            }

          }
          if (std::abs(sum) > 0.0)
            radial_power[i] = sum / n_assemblies_per_plane / volume_per_axial[i];
        }

        print_vector_in_file(radial_power,
          filename_time,
          "Radial power " + Utilities::int_to_string(step) + "\n",
          true);

      }

      if (dim == 3)
      {

        // Normalize the values of the power and flows per cell
        //      normalize_vector(power_per_assembly[0], norm);

        // Calculate the axial power distribution
        std::vector<double> power_axial;
        std::vector<std::vector<double> > phi_axial;
        unsigned int n_planes = materials.assem_per_dim[2];

        unsigned int n_assemblies_per_plane = n_assemblies / n_planes;
        Assert(n_assemblies_per_plane > 0,
          ExcMessage("n_assemblies cannot be 0"));
        radial_power.resize(n_assemblies_per_plane);
        std::vector<double> volume_per_axial(n_assemblies_per_plane);

        for (unsigned int i = 0; i < n_assemblies_per_plane; i++)
        {
          double sum = 0.0;
          for (unsigned int plane = 0; plane < n_planes; plane++)
          {
            if (std::abs(
                  power_per_assembly[0][n_assemblies_per_plane * plane
                                        + i])
                > 1e-5)
            {
              sum += power_per_assembly[0][n_assemblies_per_plane
                                           * plane
                                           + i]
                     * volume_per_assembly[n_assemblies_per_plane
                                           * plane
                                           + i];
              volume_per_axial[i] +=
                                     volume_per_assembly[n_assemblies_per_plane * plane
                                                         + i];
            }

          }
          if (std::abs(sum) > 0.0)
            radial_power[i] = sum / n_assemblies_per_plane
                              / volume_per_axial[i];
        }

        if (print_timefile)
          print_vector_in_file(radial_power,
            filename_time,
            "Radial power " + Utilities::int_to_string(step) + "\n",
            true);

      }

      // These vectors are necessary to make the approximate
      // flux weighting method
      if (dim == 3 and perturbation.rod_cusping_treat == "fluxwei")
      {
    	  std::cout<<"power.. "<<std::endl;
    	  std::cout<<"assem_per_dim[2]: "<<assem_per_dim[2]<<std::endl;
        volume_per_plane.resize(assem_per_dim[2], 0.0);
        unsigned int plane;
        power_axial.resize(assem_per_dim[2], 0.0);

        for (unsigned int i = 0; i < power_per_assembly[0].size(); i++)
        {
          plane = materials.plane(i);
          if (power_per_assembly[0][i] > 1e-5)
          {
            power_axial[plane] += power_per_assembly[0][i];
          }
          volume_per_plane[plane] += volume_per_assembly[i];
        }

        perturbation.power_axial = power_axial;
        perturbation.volume_per_plane = volume_per_plane;

      }

    }

    MPI_Bcast(&power_total, 1, MPIU_REAL, 0, comm);

    phi_norm.reinit(phi);

    phi_norm = phi;
    phi_norm *= 1.0 / power_total;

    MPI_Barrier(comm);

  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::postprocess_noise ()
  {
    PETScWrappers::MPI::BlockVector noise = phi;
    noise -= phi_critic;

    BlockVector<double> noise_serial;

    noise_serial.reinit(n_groups, n_dofs);
    for (unsigned int g = 0; g < n_groups; g++)
      noise_serial.block(g) = noise.block(g);

    if (this_mpi_process == 0)
    {

      // Initialize all that  is needed to iterate over dofs and cells
      QGauss<dim> quadrature_formula(n_fe_degree + 1);
      FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
        update_values | update_volume_elements | update_JxW_values);

      const unsigned int n_q_points = quadrature_formula.size();
      const unsigned int n_cells_out = n_assemblies;
      double noise_cell;

      // Initialize and resize the vectors where it is stored the solution
      std::vector<std::vector<double> > noise_per_assembly(n_groups,
        std::vector<double>(n_cells_out, 0.0));
      std::vector<double> local_noise(n_q_points);
      std::vector<double> volume(n_cells_out, 0.0);

      // Iterate over every cell
      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          fe_values.get_function_values(noise_serial.block(g),
            local_noise);

          noise_cell = 0.0;
          for (unsigned int q = 0; q < n_q_points; q++)
            noise_cell += local_noise[q] * fe_values.JxW(q);

          noise_per_assembly[g][cell->user_index()] += noise_cell;
        }
        volume[cell->user_index()] += cell->measure();
      }

      for (unsigned int c = 0; c < n_cells_out; c++)
        for (unsigned int g = 0; g < n_groups; ++g)
          noise_per_assembly[g][c] /= volume[c];

      std::string noi_file = out_file + ".nos";

      for (unsigned int g = 0; g < n_groups; ++g)
      {

        print_cell_distribution_in_file(dim,
          noise_per_assembly[g],
          assem_per_dim,
          noi_file,
          materials,
          "Noise of group " + num_to_str(g + 1) + " time step " + num_to_str(step)
          + "\n");
      }
      // Add Some blank lines
      std::ofstream out3(noi_file.c_str(), std::ios::app);
      out3 << "\n\n";
      out3.close();
    }

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();
  }

/**
 * @brief Output results
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::output_results ()
  {

    // Create folder if needed
//    std::size_t found = out_file.find_last_of("/\\");
//    std::string path = out_file.substr(0, found);
//    mkdir(path.c_str(), 0777);


    PETScWrappers::MPI::BlockVector noise = phi;
    noise -= phi_critic;

    BlockVector<double> noise_serial;
    BlockVector<double> phi_serial;
    BlockVector<double> phi_norm_serial;

    noise_serial.reinit(n_groups, n_dofs);
    phi_serial.reinit(n_groups, n_dofs);
    phi_norm_serial.reinit(n_groups, n_dofs);

    for (unsigned int g = 0; g < n_groups; g++)
    {
      noise_serial.block(g) = noise.block(g);
      phi_serial.block(g) = phi.block(g);
      phi_norm_serial.block(g) = phi_norm.block(g);
    }

    if (this_mpi_process == 0)
    {

      std::vector<DataComponentInterpretation::DataComponentInterpretation> dci;
      dci.push_back(DataComponentInterpretation::component_is_scalar);

      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(phi_serial.block(0), "Fast_Flux",
        DataOut<dim>::type_dof_data, dci);

      if (n_groups > 1)
      {
        data_out.add_data_vector(phi_serial.block(1), "Thermal_Flux",
          DataOut<dim>::type_dof_data, dci);
        data_out.add_data_vector(phi_norm_serial.block(1),
          "Thermal_Flux_norm", DataOut<dim>::type_dof_data, dci);
        data_out.add_data_vector(noise_serial.block(1), "Thermal_Noise",
          DataOut<dim>::type_dof_data, dci);
      }

      data_out.add_data_vector(noise_serial.block(0), "Fast_Noise",
        DataOut<dim>::type_dof_data, dci);

      data_out.build_patches(n_out_ref);

      std::string filename = out_file + Utilities::int_to_string(step) + ".vtk";
      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
    }

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();
  }



/**
 * @brief
 */
template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::run()
{

	cout << "------------ START OF THE ROM TIME LOOP ------------------"
			<< std::endl;

	cout << "Type of perturbation: " << type_perturbation << std::endl
			<< std::endl;
	double time_get_snap = timer.cpu_time();

	for (unsigned int rom = 0; rom < n_sets_snap; rom++)
	{

		step_rom = rom;
		cout << std::fixed
				<< "   Compute POD basis...                CPU Time = "
				<< timer.cpu_time() << " s." << std::endl;

		if (rom>0){
			if (snap_basis_old.size()>1)
				snap_basis_old.clear();
			snap_basis_old.resize(dim_rom);
			for (unsigned int dr=0; dr<dim_rom; dr++){
				snap_basis_old[dr].reinit(local_dofs_vector, comm);
				snap_basis_old[dr]=snap_basis[dr];
			}
		}


		compute_pod_basis(snapshots[rom]);

		cout << std::fixed
				<< "   Init time computation...                CPU Time = "
				<< timer.cpu_time() << " s." << std::endl;
		init_time_computation();

		if(rom==0){
		verbose_cout << "   Post-processing time_step...   " << std::flush;
		postprocess_time_step();

		power_vector.push_back(power_total);
		time_vect.push_back(0.0);
		if (out_flag)
		output_results();
		}

		cout << std::setprecision(4) << std::fixed << "Time:   " << 0.0000
				<< " ---> Power:  " << power_total << std::endl;

		time_init_upd = time_init_upd + delta_t;

		if (rom == 0)
			delta_t = time_intervals_snapshots[rom];
		else
			delta_t = time_intervals_snapshots[rom]
					- time_intervals_snapshots[rom - 1];

		std::cout<<"Delta_t: "<<delta_t<<std::endl;
		std::cout<<"time_init_upd: "<<time_init_upd<<std::endl;

		verbose_cout << "   Solve the ROM system...                 CPU Time = "
				<< timer.cpu_time() << " s." << std::endl;

		solve_system_petsc();

		MPI_Barrier(comm);
		verbose_cout << "   Post-processing time_step...   " << std::flush;
		postprocess_time_step();
		verbose_cout << "         CPU Time = " << timer.cpu_time() << " s"
				<< std::endl;

		verbose_cout << "      postprocess_noise..." << std::flush;
		postprocess_noise();

	}

	if (this_mpi_process == 0)
	{
		// Print Total power evolution in time
		print_vector_in_file(time_vect, out_file, "Time vector\n", true, 10);
		print_vector_in_file(power_vector, out_file, "Total Power vector\n",
				true, 10);

	}

	cout << "            Finished in " << timer.cpu_time() << " s."			<< std::endl;

	cout << "            time_get_snap " << time_get_snap << " s."			<< std::endl;
	cout << "            time_ROM: " <<  timer.cpu_time() -time_get_snap  << " s."			<< std::endl;


}

/*
 * @brief FormFunction definition for Petsc Solver
 */
template <int dim, int n_fe_degree>
  PetscErrorCode FormFunctionROM_system (TS ts,
    PetscReal time,
    Vec N,
    Vec NDOT,
    void *ctx)
  {



    TSGetApplicationContext(ts, &ctx);
    PetscScalar *ndot;
    const PetscScalar *n;

    VecGetArrayRead(N, &n);
    VecGetArray(NDOT, &ndot);

    ROMKinetics<dim, n_fe_degree> *TSobject =
                                                (ROMKinetics<dim, n_fe_degree> *)ctx;


    unsigned int n_prec = TSobject->n_prec;
    unsigned int dim_rom = TSobject->dim_rom;
    double real_time = TSobject->time_init_upd + time ; //0.0 is the initial time
//    std::cout<<"real time:"<<real_time<<std::endl;

	if (std::abs(TSobject->sim_time - real_time) > 1e-10)
	{
		TSobject->sim_time = real_time;
		TSobject->update_xsec();
		TSobject->assemble_ROM_matrices();
	}

    for (unsigned int i = 0; i < (n_prec + 1) * dim_rom; i++)
      ndot[i] = 0.0;


    Vector<double> ndot_b(dim_rom), n_b(dim_rom);
    Vector<double> auxvec(dim_rom);

    //
    // BLOCK 11 of T
    //
		for (unsigned int dr = 0; dr < dim_rom; dr++)
		  n_b[dr] = n[dr];

		// ndot=-L*n
		TSobject->romL.vmult(auxvec, n_b);
		auxvec *= -1.0;
		// auxvec=Fn+auxvec
		TSobject->romF.vmult_add(auxvec, n_b);
		// ndot_b1=LambdaInv*auxvec;
		TSobject->rominvV.vmult(ndot_b, auxvec);

		for (unsigned int dr = 0; dr < dim_rom; dr++)
		  ndot[dr] = ndot_b[dr];


    for (unsigned int k = 0; k < n_prec; k++)
    {
      unsigned int blockprec = (k + 1) * dim_rom;

      //
      // BLOCKS 1,(2...K+1)
      //
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] = n[blockprec + dr];

      // ndot= Lambdainv*lambda^d
      TSobject->rominvV.vmult(ndot_b, n_b);
//		TSobject->LapackLambda.solve(ndot_sub);
      ndot_b *= TSobject->materials.get_delayed_decay_constant(0, k);

      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[dr] += ndot_b[dr];

      //
      // BLOCKS (2...K+1),1
      //
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] = n[dr];

      // ndot= XBF(k)*n
      TSobject->romXBF[k].vmult(ndot_b, n_b);

      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[blockprec + dr] += ndot_b[dr];

      //
      // BLOCKS (2,2),...,(K+1,K+1)
      //
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] = n[blockprec + dr];

      // ndot= beta*(n+n*coeffM*A_M)
      ndot_b.equ(-1.0 * TSobject->materials.get_delayed_decay_constant(0, k), n_b);

      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[blockprec + dr] += ndot_b[dr];

    }


    VecRestoreArrayRead(N, &n);
    VecRestoreArray(NDOT, &ndot);

    return 0;
  }

template class ROMKinetics<1, 1> ;
template class ROMKinetics<1, 2> ;
template class ROMKinetics<1, 3> ;
template class ROMKinetics<1, 4> ;
template class ROMKinetics<1, 5> ;

template class ROMKinetics<2, 1> ;
template class ROMKinetics<2, 2> ;
template class ROMKinetics<2, 3> ;
template class ROMKinetics<2, 4> ;
template class ROMKinetics<2, 5> ;
//
template class ROMKinetics<3, 1> ;
template class ROMKinetics<3, 2> ;
template class ROMKinetics<3, 3> ;
template class ROMKinetics<3, 4> ;
template class ROMKinetics<3, 5> ;






