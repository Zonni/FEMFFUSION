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
    const bool silent) :
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


    // Solver parameters
    init_delta_t = prm.get_double("Time_Delta");
    tol_time_ksp = static_problem.tol_ksp;



    step = 0;
    print_step = 1;
    delta_t.push_back(init_delta_t);
    sim_time = 0.0;
    power_total = 0.0;


    // Initialize phi //TODO remove
    phi.reinit(local_dofs_vector, comm);
    static_problem.phi[0].compress(VectorOperation::insert);
    phi = static_problem.phi[0];
    phi_critic = static_problem.phi[0];


    // Reinit ROM data
    type_snapshots=' '; // modes or dyn_pos
    n_snap=2;

    get_snapshots(static_problem);

    dim_rom=snapshots.size(); //TODO

    get_parameters_from_command_line();

    snap_basis.resize(dim_rom);
    for (unsigned int dr=0; dr<dim_rom; dr++)
    snap_basis[dr].reinit(local_dofs_vector, comm);

    if (prec_flag == false)
    {
      n_prec = 0;
      materials.remove_precursors();
    }

    // Initialization
    delta_t.reserve(t_end / init_delta_t);
    time_vect.reserve(t_end / init_delta_t);
    power_vector.reserve(t_end / init_delta_t);

    verbose_cout << "Initialize the perturbation class" << std::endl;
    perturbation.init_transient();

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
      for (unsigned int nb = 0; nb < dim_rom; nb++)
        N[nb] = snap_basis[nb]* phi_critic;


      // Compute c[0]=(c11,...,cq1,...,c1k,...,cqk)
      std::vector<std::vector<double>> cjk(dim_rom,
        std::vector<double>(n_prec));
      PETScWrappers::MPI::BlockVector aux(snap_basis[0]);


	for (unsigned int nb = 0; nb < dim_rom; nb++) {
		for (unsigned int k = 0; k < n_prec; k++) {
			AssertRelease(materials.get_delayed_decay_constant(0, k) > 0.0,
					"Delayed decay constant must be greater than 0.");
			cjk[nb][k] = 1.0 / (materials.get_delayed_decay_constant(0, k))
					* XBF[k].vmult_dot(snap_basis[nb], phi_critic);

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
  void
  ROMKinetics<dim, n_fe_degree>::get_snapshots (
      StaticDiffusion<dim, n_fe_degree> &static_problem)
  {


    // Initialize snapshots
    get_uint_from_options("-n_snap", n_snap);

    if (type_snapshots == "modes")
      n_snap = static_problem.n_eigenvalues;

    snapshots.resize (n_snap);
    for (unsigned int ns = 0; ns < n_snap; ns++)
      snapshots[ns].reinit (local_dofs_vector, comm);

    if (type_snapshots == "modes")
      {
	for (unsigned int ns = 0; ns < n_snap; ns++)
	  {
	    static_problem.phi[ns].compress (VectorOperation::insert);
	    snapshots[ns] = static_problem.phi[ns];
	  }
	return;
      }

    unsigned int bar_bank = 1;
    double bars_bottom_pos = 0.0;
    AssertRelease(n_snap>1,"The number of snapshots must be greater than 1");
    double step_bar = (perturbation.bars_top_pos - bars_bottom_pos)
	/ (n_snap - 1);

    unsigned int bar_material = perturbation.bar_materials[bar_bank];

    for (unsigned int ns = 0; ns < n_snap; ns++)
      {

	unsigned int bar_pos_z = ns * step_bar + bars_bottom_pos;

	for (unsigned int plant_pos = 0;
	    plant_pos < perturbation.bars_position.size (); ++plant_pos)
	  {

	    if (perturbation.bars_position[plant_pos] == bar_bank){

	      perturbation.move_bar_volume_homogenized (plant_pos, bar_pos_z,
							bar_material, bar_bank);
	    }

	  }

	static_problem.cout.set_condition (false);
	static_problem.show_eps_convergence = false;
	static_problem.assemble_system_lambda ();
	static_problem.solve_eps ();
	static_problem.phi[0].compress (VectorOperation::insert);
	snapshots[ns] = static_problem.phi[0];
	cout << "Step " << ns << ", Bar Position: " << bar_pos_z
	    << ", Eigenvalue: " << static_problem.eigenvalues[0] << std::endl;
      }

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
    get_uint_from_options("-dim_rom", dim_rom);

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
        or type_perturbation == "Out_Of_Phase")
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
  void ROMKinetics<dim, n_fe_degree>::compute_pod_basis ()
  {


	// Create the matrix with the snapshot to apply the SVD
	Mat Mat_snap;
	PetscInt i_snap, i_sv;
	PetscInt *idm = new PetscInt[n_dofs*n_groups];
	PetscScalar *values_snap = new PetscScalar[n_dofs*n_groups];

	for (unsigned int j=0; j<n_dofs*n_groups; j++)
		idm[j]=j;

	MatCreateDense(comm, n_dofs*n_groups, snapshots.size(), n_dofs*n_groups, snapshots.size(), NULL, &Mat_snap);
	for (i_snap=0; i_snap<static_cast<int>(snapshots.size()); i_snap++){
		for (unsigned k=0; k<n_dofs*n_groups; k++)
			values_snap[k]=snapshots[i_snap][k];
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
   	std::cout<<"Singular value "<<i_sv<<": "<<singular_values[i_sv]<<std::endl;
	}

	MatDestroy(&Mat_snap);
	SVDDestroy(&svd);
	u.clear();


  }


/*
 * @brief Assemble Time System
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

	rominvV.invert(romV);

	V.clear();
	L.clear();
	F.clear();
	for(unsigned int p=0; p<n_prec; p++)
		XBF[p].clear();

  }



/*
 * @brief Compute the RHS called E
 */
/*template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::compute_RHS ()
  {
    double exp_factor;

    // E =R*phi_old
    Rfree.vmult(E, phi);
    Rfree.clear();
    E *= 1 / delta_t[step];

    if (time_scheme == "implicit-exponential")
    {
      // Precursors term
      std::vector<PETScWrappers::MPI::BlockVector> XCk(n_prec);
      for (unsigned int np = 0; np < n_prec; np++)
      {
        XCk[np].reinit(local_dofs_vector, comm);

        for (unsigned int ng = 0; ng < n_groups; ng++)
          XCk[np].block(ng).equ(materials.get_delayed_spectra(0, np, ng),
            PCk[np]);
      }

      for (unsigned int p = 0; p < n_prec; ++p)
      {
        exp_factor = materials.get_delayed_decay_constant(0, p)
                     * exp(
                       -materials.get_delayed_decay_constant(0, p)
                       * delta_t[step]);
        E.add(exp_factor, XCk[p]);
      }
    }
    else if (time_scheme == "semi-implicit-exponential"
             or time_scheme == "semi-implicit-euler")
    {

      PETScWrappers::MPI::BlockVector XCp;
      XCp.reinit(local_dofs_vector, comm);

      // Obtain Ck
      if (time_scheme == "semi-implicit-exponential")
        for (unsigned int p = 0; p < n_prec; p++)
          KSPSolve(kspP, PCk[p], Ck[p]);

      // Precursors term
      for (unsigned int p = 0; p < n_prec; ++p)
      {
        XCp *= 0.0;
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          X[p][g]->vmult_add(XCp.block(g), Ck[p]);
        }
        E.add(1.0, XCp);
      }

//		E.print(std::cout);

    }
    else
    {
      AssertRelease(false, "Time scheme not available");
    }

    if (print_rhs)
    {

      std::string print_time_matrices_matlab;
      get_string_from_options("-print_time_matrices_matlab",
        print_time_matrices_matlab);

      if (!print_time_matrices_matlab.empty())
      {

        std::string name_file;
        name_file.append(print_time_matrices_matlab.begin(),
          print_time_matrices_matlab.end() - 2);
        name_file.append("_s" + std::to_string(step) + "_e.m");
        out_matlab.open(name_file.c_str(), std::ios::out);

        print_block_vector_in_matlab(E, name_file, out_matlab, 12);

        out_matlab.close();
      }

      std::string print_time_matrices_python;
      get_string_from_options("-print_time_matrices_python",
        print_time_matrices_python);

      if (!print_time_matrices_python.empty())
      {

        for (unsigned int g1 = 0; g1 < n_groups; g1++)
        {
          std::string name = "E" + std::to_string(g1);
          print_vector_in_python(E.block(g1), name, out_matlab, 12);
        }

        out_matlab << "E=[];" << std::endl;
        out_matlab << "E=np.array(E);" << std::endl;
        for (unsigned int g1 = 0; g1 < n_groups; g1++)
          out_matlab << "E=np.append(E,E" + std::to_string(g1) + ");";

        for (unsigned int g1 = 0; g1 < n_groups; g1++)
          out_matlab << "del E" + std::to_string(g1) + ";";

        out_matlab.close();
      }
    }

  }
*/

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::solve_system_petsc ()
  {

    double init_time_step_sundials = 1e-2;
    unsigned int max_steps = 1e+2;
//    unsigned int maxl;


// SOLVER PETSC implemented for a linear system
    TS ts;

//    vecn.resize(n_eigenvalues);

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
//       is the application context used by the call-back functions.

    TSSetRHSFunction(ts, NULL, FormFunctionROM_system<dim, n_fe_degree>,
        this);


    		//   Form the initial guess for the problem
    //   (this can be done from a function FormInitialGuess)

    //   This indicates that we are using pseudo timestepping to
    //   find a steady state solution to the nonlinear problem.
    TSSetType(ts, TSBDF);

    //   Set the initial time to start at (this is arbitrary for
    //   steady state problems); and the initial timestep given above
//    TSSetTimeStep(ts, 1e-5);
    TSSetTimeStep(ts, init_time_step_sundials);

    //   Set a large number of timesteps and final duration time
    //   to insure convergence to steady state.
    //   TSSetMaxSteps(ts, max_steps);
//    if (update_modes == true)
    TSSetMaxTime(ts, t_end);
//    else
//      TSSetMaxTime(ts, delta_t_petsc);
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_INTERPOLATE);
//    TSSetPostStep(ts, PostStep<dim, n_fe_degree>);
//    TSSetApplicationContext(ts, this);
//    TSSundialsMonitorInternalSteps(ts, PETSC_TRUE);
//    TSSundialsSetType(ts, SUNDIALS_BDF);
//    //TSSundialsSetMaxl(ts, 50);
//    //TSSundialsSetLinearTolerance(ts,1.0);
//    TSSundialsSetMinTimeStep(ts, 1e-12);
//    TSSundialsSetMaxTimeStep(ts, 1e-1);
//    TSSundialsSetTolerance(ts, tol_sundials_rel, tol_sundials_abs);

    //   Set any additional options from the options database. This
    //   includes all options for the nonlinear and linear solvers used
    //   internally the timestepping routines.
    TSSetFromOptions(ts);
    TSSetUp(ts);

    //   Perform the solve. This is where the timestepping takes place.
    TSSolve(ts, coeffs_n);

    const PetscScalar *narray;

    VecGetArrayRead(coeffs_n, &narray);

    phi = 0.0;
    for (unsigned int dr = 0; dr < dim_rom; ++dr)
      phi.add(narray[dr], snap_basis[dr]);

    postprocess_time_step();

    cout << std::setprecision(4) << "Time:   " << t_end
         << " ---> Power:  "
         << power_total << std::endl;

//    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
//      vecn[eig].push_back(narray[eig]);

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


///*
// * @brief solve_LHS
// */
//template <int dim, int n_fe_degree>
//  void ROMKinetics<dim, n_fe_degree>::solve_LHS ()
//  {
//
//    PETScWrappers::MPI::Vector phivec(comm, n_groups * n_dofs,
//      n_groups * locally_owned_dofs.n_elements());
//    PETScWrappers::MPI::Vector Evec(comm, n_groups * n_dofs,
//      n_groups * locally_owned_dofs.n_elements());
//
//    if (type_perturbation != "Step_Change_Material" or step < 2 or reinit_step< 2)
//    {
//
//
//      // Setup the solver
//      KSPCreate(comm, &ksp);
//      KSPGetPC(ksp, &pc);
//      KSPSetTolerances(ksp, tol_ksp, tol_ksp, PETSC_DEFAULT, 3000);
//      KSPSetType(ksp, KSPGMRES);
//      PCSetType(pc, PCSHELL);
//      PCShellSetApply(pc, apply_preconditioner<dim, n_fe_degree>);
//      PCShellSetContext(pc, this);
//      KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
//      KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
//      KSPSetFromOptions(ksp);
//
//      KSPSetOperators(ksp, shell_T, shell_T);
//      KSPSetUp(ksp);
//    }
//
//    copy_to_Vec(Evec, E);
//    copy_to_Vec(phivec, phi);
//
//    KSPSolve(ksp, Evec, phivec);
//    copy_to_BlockVector(phi, phivec);
//
//    old_its = its;
//    its = 0;
//    KSPGetIterationNumber(ksp, &its);
//    solver_its.push_back(its);
//    cpu_time.push_back(timer.cpu_time());
//    totalits += its;
//
//    if (step % out_interval == 0)
//      cout << "   its: " << its << std::endl;
//
//    solve_precursors();
//
//    phivec.clear();
//    Evec.clear();
//
//    if (type_perturbation != "Step_Change_Material" or step < 1 or reinit_step< 1)
//    {
//      Tfree.clear();
//      MatDestroy(&shell_T);
//      KSPDestroy(&ksp);
//    }
//
//  }

/*
 * @brief solve_LHS
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::solve_precursors ()
  {

//    double expFactor, ak;
//    PETScWrappers::MPI::Vector Mxphi(locally_owned_dofs, comm);
//
//    if (time_scheme == "implicit-exponential"
//        or time_scheme == "semi-implicit-exponential")
//    {
//      // Computation of XC_k
//      for (unsigned int k = 0; k < n_prec; ++k)
//      {
//        Mxphi = 0.0;
//        for (unsigned int ng = 0; ng < n_groups; ++ng)
//          BMfree.vmult_add(ng, k, Mxphi, phi.block(ng));
//
//        // In Ck_update it is accumulated the produced precursors
//        expFactor = exp(
//          -materials.get_delayed_decay_constant(0, k)
//          * delta_t[step]);
//        ak = 1.0 / materials.get_delayed_decay_constant(0, k)
//             * (1 - expFactor);
//        //Computes Ck = expFactor *  Ck[k] + ak* Mxphi;
//        PCk[k].sadd(expFactor, ak, Mxphi);
//
//      }
//
//      Mxphi.clear();
//
//    }
//    else if (time_scheme == "semi-implicit-euler")
//    {
//
//      for (unsigned int p = 0; p < n_prec; p++)
//      {
//        PCk[p] *= 0.0;
//        P.vmult(PCk[p], Ck[p]);
//        PCk[p] /= delta_t[step];
//      }
//
//      PETScWrappers::MPI::Vector Mxphi(locally_owned_dofs, comm);
//
//      // Computation of XC_k
//      for (unsigned int np = 0; np < n_prec; ++np)
//      {
//        Mxphi *= 0.0;
//
//        for (unsigned int ng = 0; ng < n_groups; ++ng)
//          BMfree.vmult_add(ng, np, Mxphi, phi.block(ng));
//
//        PCk[np].add(1.0, Mxphi);
//
//        KSPSolve(kspLP[np], PCk[np], Ck[np]);
//
//      }
//
//      Mxphi.clear();
//
//      if (type_perturbation != "Step_Change_Material" or step < 1 or reinit_step< 1)
//      {
//        for (unsigned int np = 0; np < n_prec; np++)
//        {
//          KSPDestroy(&kspLP[np]);
//          LP[np]->clear();
//        }
//      }
//
//    }
//    else
//    {
//      AssertRelease(false, "Invalid type of time scheme.");
//    }
//
//    BMfree.clear();

  }

/*
 * @brief FormFunction definition for Petsc Solver
 */
template <int dim, int n_fe_degree>
  PetscErrorCode PostStep (TS ts)
  {

	/*
    void *ctx;
    TSGetApplicationContext(ts, &ctx);

    ROMKinetics<dim, n_fe_degree>* TSobject =
                                                (ROMKinetics<dim, n_fe_degree>*) ctx;
    PetscReal deltatime;
    const PetscScalar *n;
    Vec N;

    TSGetTime(ts, &deltatime);
    TSGetSolution(ts, &N);
    VecGetArrayRead(N, &n);

    double transf_time = deltatime + TSobject->t_init_upd;

    TSobject->step++;

    for (unsigned int eig = 0; eig < TSobject->n_eigenvalues; ++eig)
    {
//         power_total+=n[eig]*TSobject->power_modes[eig];
      TSobject->vecn[eig].push_back(n[eig]);
    }

    TSobject->phi = 0.0;

    for (unsigned int eig = 0; eig < TSobject->n_eigenvalues; ++eig)
    {
      TSobject->phi.add(n[eig], TSobject->phi_modes[eig]);
    }

    TSobject->postprocess_time_step(transf_time);

    if (transf_time < TSobject->t_end_upd)
    {
      if (std::abs((TSobject->old_power - TSobject->power_total) / TSobject->old_power) < 2.0)
      {

        TSobject->step_time_out.push_back(TSobject->step);
        TSobject->power_vector.push_back(TSobject->power_total);
        TSobject->time_vect.push_back(transf_time);

        TSobject->old_power = TSobject->power_total;
      }

      TSobject->cout << std::setprecision(4) << std::fixed << "Time:   "
      << transf_time
      << " ---> Power:  " << TSobject->power_total
      << std::endl;
    }

    VecRestoreArrayRead(N, &n);

//    TSobject->tseg = transf_time;
//    VecCopy(N, TSobject->nseg);

//    VecDestroy(&N);

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
    */
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
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::run ()
  {
    cout << "------------ START OF THE ROM TIME LOOP ------------------"
         << std::endl;

    cout << "Type of perturbation: " << type_perturbation << std::endl
             << std::endl;


    cout << "IT IS ASSUMED THAT LAMBDA (DELAYED DECAY CONSTANT) ARE CONSTANT!!"<< std::endl
                 << std::endl;

    verbose_cout << std::fixed
                         << "   Compute POD basis...                CPU Time = "
                         << timer.cpu_time() << " s." << std::endl;
    compute_pod_basis();

//    verbose_cout << std::fixed
//                     << "   Assemble matrices...                CPU Time = "
//                     << timer.cpu_time() << " s." << std::endl;
//    assemble_matrices();

    verbose_cout << std::fixed
                 << "   Init time computation...                CPU Time = "
                 << timer.cpu_time() << " s." << std::endl;
    init_time_computation();

    verbose_cout << "   Post-processing time_step...   " << std::flush;
    postprocess_time_step();


        cout <<" ---> Power:  "
             << power_total << std::endl;

    verbose_cout << "   Solve the ROM system...                 CPU Time = "
                         << timer.cpu_time()
                         << " s." << std::endl;
    solve_system_petsc();


//      if (type_perturbation == "Mechanical_Vibration")
//        perturbation.restore_indices();

      MPI_Barrier(comm);
      verbose_cout << "   Post-processing time_step...   " << std::flush;
      postprocess_time_step();
      verbose_cout << "         CPU Time = " << timer.cpu_time() << " s"
                   << std::endl;

      verbose_cout << "      postprocess_noise..." << std::flush;
      postprocess_noise();

      // ------------------------------------------------------------------------
      //if (step % out_interval == 0)
      //{
//      cout << " Step " << step << " at t=" << sim_time << std::endl;
//      cout << "                         Total Power " << power_total
//           << "   Time = "
//           << timer.cpu_time() << std::endl;
//
//      verbose_cout << "   Step Done." << " Time = " << timer.cpu_time()
//                   << " s."
//                   << std::endl;
      //}

      // ------------------------------------------------------------------------

//      if (out_flag and step % out_interval == 0)
//      {
//        MPI_Barrier(comm);
//        verbose_cout << " Done!" << std::endl;
//        output_results();
//        MPI_Barrier(comm);
//      }
//
//      PetscMemorySetGetMaximumUsage();
//      PetscLogDouble memory;
//      PetscMemoryGetMaximumUsage(&memory);
//      cout << "   Max Memory " << memory * 1e-6 << " MB" << std::endl;

//      verbose_cout << "---------------------------------------------------"
//                   << std::endl;

      // Out things in the future will be a function
//      time_vect.push_back(sim_time);
//      power_vector.push_back(power_total);
//      delta_t.push_back(init_delta_t);

//      step++;


//      MPI_Barrier(comm);


//    if (this_mpi_process == 0)
//    {
//      // Print Total power evolution in time
//      print_vector_in_file(time_vect, out_file, "Time vector\n", true, 10);
//      print_vector_in_file(power_vector, out_file, "Total Power vector\n", true, 10);
//      print_vector_in_file(delta_t, out_file, "Delta t\n", true, 10);
//      print_vector_in_file(cpu_time, out_file, "CPU Time\n", true, 10);
//      print_vector_in_file(solver_its, out_file, "Solver Its\n", true, 10);
//
//      if (print_timefile)
//        print_vector_in_file(print_time_vect, filename_time, "Time\n", true, 10);
//    }
//
//    if (out_matlab.is_open())
//      out_matlab.close();

//    cout << "Total its: " << totalits << ", mean by it:"
//         << double(totalits) / step
//         << std::endl;


    cout << "            Finished in " << timer.cpu_time() << " s." << std::endl;

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
    double transf_time = time ; //0.0 is the initial time

    if (std::abs(TSobject->sim_time-time)>1e-6){
    std::cout<<"time: "<<time<<std::endl;
    TSobject->sim_time=transf_time;
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






