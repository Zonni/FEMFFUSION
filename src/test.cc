/**
 * @file   test.cc
 * @brief  Implementation some tests. Execute them with ./femffusion -t
 */
#include <deal.II/base/parameter_handler.h>
#include <math.h>

#include "../include/test.h"
#include "../include/utils.h"
#include "../include/femffusion.h"
#include "../include/static_diffusion.h"
#include "../include/static_spn.h"
#include "../include/static_full_spn.h"
#include "../include/time_computation.h"

#include "../include/noise/noise_diffusion.h"
#include "../include/noise/noise_spn.h"
#include "../include/noise/noise_full_spn.h"

#include "../include/rom/rom_kinetics.h"
#include "../include/rom/rom_utils.h"


using namespace dealii;

/**
 * @brief Run all test asserting the right values.
 */
int run_tests ()
{
  std::string input_file, out_file;


  //
  // ---------------------- TEST ROM ----------------------- //
  //
  //
  // Test
  test_POD_groupwise();
  test_LUPOD_extended ();
  test_LUPOD_extended_group_wise();


//  // input_file = "test/1D_hom_5cells_mov_prec/1D_5cells_f1_rom.prm";
//  run_test_rom_LUPOD_1();
//  // input_file = "test/1D_hom_5cells_mov_prec/1D_5cells_f1_rom.prm";
//  run_test_rom_LUPOD_2();



  //input_file = "test/3D_Langenbuch/3D_Langenbuch_rom_rods.prm";
  //run_test_static_rom(input_file, 3);

#ifdef DEBUG
  std::cout<<" Compiled in DEBUG mode! "<<std::endl;
  std::cout<<" Tests does not work in this mode. "<<std::endl;
  exit(0);
#endif

  //-----------------------------------------------------------------//
  std::cout << std::endl;
  std::cout << "----------------------------" << std::endl;
  std::cout << "    RUNNING TESTS           " << std::endl;
  std::cout << "----------------------------" << std::endl;

  // Run test in utils.cc
  run_tests_utils();

  //-----------------------------------------------------------------//
  // TEST 1D_hom_slab_2cm
  input_file = "test/1D_hom_slab_2cm/1D_hom_zerocurrent.prm";
  out_file = "test/1D_hom_slab_2cm/1D_hom_zerocurrent.out";
  test_keff_problem(input_file, 2.5);

  std::vector<double> power(
      { 1.0 });
  chek_output_vector(out_file, "Neutron Power", power);
  std::vector<double> eig(
      { 2.5 });
  chek_output_vector(out_file, "The Eigenvalues are:", eig);
  std::vector<double> ff(
      { 10.0 });
  chek_output_vector(out_file, "Group 1 flux", ff);
  std::vector<double> tf(
      { 0.0 });
  chek_output_vector(out_file, "Group 2 flux", tf);

  input_file = "test/1D_hom_slab_2cm/1D_hom_zeroflux.prm";
  test_keff_problem(input_file, 0.271012);

  input_file = "test/1D_hom_slab_2cm/1D_hom_vacuum.prm";
  test_keff_problem(input_file, 0.587489);

  //
  //-----------------------------------------------------------------//
  // TEST 1D_hom_slab_10cm
  input_file = "test/1D_hom_slab_10cm/1D_hom_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/1D_hom_slab_10cm/1D_hom_zeroflux.prm";
  test_keff_problem(input_file, 1.88113);

  input_file = "test/1D_hom_slab_10cm/1D_hom_vacuum.prm";
  test_keff_problem(input_file, 1.98917);

  //
  //-----------------------------------------------------------------//
  // TEST 1D_hom_slab_2cm_SP3
  input_file = "test/1D_hom_slab_2cm_SP3/1D_hom_SP3.prm";
  test_keff_problem(input_file, 0.6529562);

  input_file = "test/1D_hom_slab_2cm_SP3/1D_hom_SP5.prm";
  test_keff_problem(input_file, 0.660523);

  input_file = "test/1D_hom_slab_2cm_SP3/1D_hom_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/1D_hom_slab_2cm_SP3/1D_hom_zeroflux.prm";
  test_keff_problem(input_file, 0.364217);

  //
  //-----------------------------------------------------------------//
  // TEST 1D_hom_slab_2g_2cm
  input_file = "test/1D_hom_slab_2g_2cm/1D_hom_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/1D_hom_slab_2g_2cm/1D_hom_zeroflux.prm";
  test_keff_problem(input_file, 0.271012);

  input_file = "test/1D_hom_slab_2g_2cm/1D_hom_vacuum.prm";
  test_keff_problem(input_file, 0.587489);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_hom_slab_2cm
  // pseudo 2D problem -> Must be the same results as 1D_hom_2cm
  input_file = "test/2D_hom_slab_2cm/2D_hom_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/2D_hom_slab_2cm/2D_hom_zeroflux.prm";
  test_keff_problem(input_file, 0.271012);

  input_file = "test/2D_hom_slab_2cm/2D_hom_vacuum.prm";
  test_keff_problem(input_file, 0.587489);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_hom_2cm
  // Real 2D problem with 2g
  input_file = "test/2D_hom/2D_hom_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/2D_hom/2D_hom_zeroflux.prm";
  test_keff_problem(input_file, 0.143272);

  input_file = "test/2D_hom/2D_hom_vacuum.prm";
  test_keff_problem(input_file, 0.332854);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_hom_multigroup
  //
  input_file = "test/2D_hom_multigroup/2D_hom_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/2D_hom_multigroup/2D_hom_zeroflux.prm";
  test_keff_problem(input_file, 0.143272);

  input_file = "test/2D_hom_multigroup/2D_hom_vacuum.prm";
  test_keff_problem(input_file, 0.332854);

  //
  //-----------------------------------------------------------------//
  // TEST   2D_biblis_test
  // Testing degenerated eigenvalues (2 or more equal eigenvalues due to symmetry).
  input_file = "test/2D_biblis_test/biblis_FE3.prm";
  out_file = "test/2D_biblis_test/biblis_FE3.out";
  test_keff_problem(input_file, 1.02533);
  std::vector<double> eigs(
      { 1.02533, 1.01838, 1.01838, 1.00778, 1.00692 });
  chek_output_vector(out_file, "The Eigenvalues are:", eigs, eigs.size(), 1e-5);

  input_file = "test/2D_biblis_test/biblis_SP3.prm";
  test_keff_problem(input_file, 1.02201);

  // Bifpam solver with FE=1 initialization
  input_file = "test/2D_biblis_test/biblis_bipfam.prm";
  test_keff_problem(input_file, 1.02533);
  out_file = "test/2D_biblis_test/biblis_bifpam.out";
  chek_output_vector(out_file, "The Eigenvalues are:", eigs, eigs.size(), 1e-5);

  // Bifpam solver with FE=1 initialization for SP3
  input_file = "test/2D_biblis_test/biblis_SP3_bifpam.prm";
  test_keff_problem(input_file, 1.02585);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_hex_hom
  input_file = "test/2D_hex_hom/2D_hex_zerocurrent.prm";
  test_keff_problem(input_file, 0.626177);

  input_file = "test/2D_hex_hom/2D_hex_vacuum.prm";
  test_keff_problem(input_file, 0.531227);

  input_file = "test/2D_hex_hom/2D_hex_zeroflux.prm";
  test_keff_problem(input_file, 0.516308);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_pin_structured
  input_file = "test/2D_pin_structured/pin_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/2D_pin_structured/pin_zeroflux.prm";
  test_keff_problem(input_file, 0.143272);

  input_file = "test/2D_assembly_structured/assembly_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);
  out_file = "test/2D_assembly_structured/assembly_zerocurrent.out";
  std::vector<double> power2(
      { 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
        1.000000e+00,
        1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
        1.000000e+00,
        1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
        1.000000e+00,
        1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
        1.000000e+00,
        1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00 });
  chek_output_vector(out_file, "Neutron Power", power2, 5);

  input_file = "test/2D_assembly_structured/assembly_vacuum.prm";
  test_keff_problem(input_file, 0.332682);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_pin_unstructured
  input_file = "test/2D_pin_unstructured/pin_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/2D_pin_unstructured/pin_zeroflux.prm";
  test_keff_problem(input_file, 0.143272);

  input_file = "test/2D_pin_unstructured/pin_vacuum.prm";
  test_keff_problem(input_file, 0.332854);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_pin_save_load
  std::remove("test/2D_pin_save_load/mesh.tri.vec");
  std::remove("test/2D_pin_save_load/mesh.tri");

  input_file = "test/2D_pin_save_load/pin_vacuum.prm";
  test_keff_problem(input_file, 0.330654);
  // Make it twice to ensure load mesh is working properly
  input_file = "test/2D_pin_save_load/pin_vacuum.prm";
  test_keff_problem(input_file, 0.330654);

  //
  //-----------------------------------------------------------------//
  // TEST 2D_assembly_unstructured
  input_file = "test/2D_assembly_unstructured/assembly_zerocurrent.prm";
  test_keff_problem(input_file, 2.5);

  input_file = "test/2D_assembly_unstructured/assembly_vacuum.prm";
  test_keff_problem(input_file, 0.332805);

  //
  //-----------------------------------------------------------------//
  // 3D HOLES
  input_file = "test/3D_holes/holes.prm";
  out_file = "test/3D_holes/holes.out";
  test_keff_problem(input_file, 0.225029);
  // FIXME
  /*
   std::vector<double> holes_power(
   { 9.263313e-01, 9.774135e-01, 9.922412e-01,
   7.589613e-01,
   1.377407e+00, 1.155487e+00, 8.344543e-01, 9.777047e-01 });
   std::vector<double> holes_plane1(
   { 9.263313e-01, 9.774135e-01 });
   std::vector<double> holes_plane2(
   { 9.922412e-01, 7.589613e-01, 1.377407e+00, 1.155487e+00 });
   std::vector<double> holes_plane3(
   { 8.344543e-01, 9.777047e-01 });
   //chek_output_vector(out_file, "Axial power distribution:", axial_power, 3, 1e-6);
   chek_output_vector(out_file, "Neutron Power", holes_power, 5, 1e-6);
   chek_output_vector(out_file, "Plane 1", holes_plane1, 2, 1e-6);
   chek_output_vector(out_file, "Plane 2", holes_plane2, 2, 1e-6);
   chek_output_vector(out_file, "Plane 3", holes_plane3, 1, 1e-6);
   */
  //
  //-----------------------------------------------------------------//
  // 3D_IAEA
  input_file = "test/3D_IAEA_test/IAEA_FE1.prm";
  test_keff_problem(input_file, 1.04493);

  //
  //-----------------------------------------------------------------//
  // 3D_Roseton
  // with .bar file
  input_file = "test/3D_roseton_test/roseton_FE1.prm";
  test_keff_problem(input_file, 0.801294);

  //
  //-----------------------------------------------------------------//
  // TEST 3D_pin_structured
  input_file = "test/3D_structured/pin2.prm";
  test_keff_problem(input_file, 0.147951);

  input_file = "test/3D_structured/pin_axialrefs.prm";
  test_keff_problem(input_file, 0.147951);

  input_file = "test/3D_structured/assembly.prm";
  test_keff_problem(input_file, 0.117764);

  //
  // --------------------------------------------------------------- //
  // FULL Matrix free tests
  input_file = "test/3D_Langenbuch/3D_Langenbuch_fe1_alpha.prm";
  test_keff_problem(input_file, 14.5996, 2e-4);
  input_file = "test/3D_Langenbuch/3D_Langenbuch_fe1_gamma.prm";
  test_keff_problem(input_file, 1.00031);

  //
  // --------------------------------------------------------------- //
  // ALPHA AND GAMMA MODES
  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_matrixfree/biblis_FE1_fullmatrixfree.prm";
  test_keff_problem(input_file, 1.02180);
  input_file = "test/2D_biblis_matrixfree/biblis_FE1_fullmatfree_zerocurrent.prm";
  test_keff_problem(input_file, 1.02183);

  // --------------------------------------------------------------- //
  //
  // --------------------------------------------------------------- //
  // NOISE
  // --------------------------------------------------------------- //

  //
  // --------------------------------------------------------------- //
  // 2D BIBLIS Noise
  // --------------------------------------------------------------- //
  input_file = "test/2D_biblis_noise/biblis_FE2.prm";
  test_noise_problem(input_file, 0.816894, 1e-4);

  //
  // --------------------------------------------------------------- //
  // 3D IAEA Noise
  // --------------------------------------------------------------- //
  input_file = "test/3D_IAEA_noise/IAEA_FE1.prm";
  test_noise_problem(input_file, 3.22249);

  //
  // --------------------------------------------------------------- //
  // 2D Rectangular Borders Perturbation
  // --------------------------------------------------------------- //
  input_file = "test/2D_test_vibration/2D_test.prm";
  test_noise_problem(input_file, 0.465861);

  input_file = "test/2D_test_vibration/2D_test_ref.prm";
  test_noise_problem(input_file, 0.389968);

  //
  // --------------------------------------------------------------- //
  // 2D Hexagonal Borders Perturbation
  // --------------------------------------------------------------- //
  input_file = "test/2D_test_vibration_hex/2D_test.prm";
  test_noise_problem(input_file, 0.465927);

  //
  // --------------------------------------------------------------- //
  // 1D SPN Noise
  // --------------------------------------------------------------- //
  // Not converged
  //input_file = "test/1D_noise_SPN_1g/1D_noise_diffusion.prm";
  //test_noise_problem(input_file, 26.0953);
  input_file = "test/1D_noise_SPN_1g/1D_noise_sp1.prm";
  test_noise_problem(input_file, 26.095597);
  input_file = "test/1D_noise_SPN_1g/1D_noise_fsp1.prm";
  test_noise_problem(input_file, 28.2524);

  input_file = "test/1D_noise_SPN/1D_noise_diffusion.prm";
  test_noise_problem(input_file, 0.343722);
  input_file = "test/1D_noise_SPN/1D_noise_sp3.prm";
  test_noise_problem(input_file, 0.3433418244);
  input_file = "test/1D_noise_SPN/1D_noise_fsp3.prm";
  test_noise_problem(input_file, 0.349394);

  //
  // --------------------------------------------------------------- //
  //  TEST MESHES

  // Test rectangular mesh + distributed + Rods:

  double power_double[] =
                            { 9.999999e-01, 1.166752, 1.14599 };
  std::vector<double> power_rect(power_double,
    power_double + sizeof(power_double) / sizeof(double));

  input_file = "test/rectangular/rectangular_mov.prm";
  test_keff_problem(input_file, 0.978199);
  test_power_evolution(input_file, power_rect, 1e-4);

  // Test hexagonal mesh + distributed + Rods:"
  double power_double_hex_ds[] =
                                   { 1.00, 2.20851 };
  std::vector<double> power_double_hex_d(power_double_hex_ds,
    power_double_hex_ds + sizeof(power_double_hex_ds) / sizeof(double));

  input_file = "test/3D_roseton/roseton_ds.prm";
  test_keff_problem(input_file, 0.801286);
  test_power_evolution(input_file, power_double_hex_d, 1e-4);

  ////-----------------------------------------------------------------//

  // Test Type of Perturbations:

  // TEST 3D Langenbuch
  // Test Langenbuch reactor with Rods:"

  std::vector<double> power_one(2, 1.0);
  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds.prm";
  test_keff_problem(input_file, 1.00023);
  test_power_evolution(input_file, power_one, 1e-4);
  power_one.resize(9, 1.0);
  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds_rods.prm";
  test_power_evolution(input_file, power_one, 1e-4);
  double power_langenbuch[] =
                                { 1, 1.06485, 1.15686 };
  std::vector<double> power_langen(power_langenbuch,
    power_langenbuch + sizeof(power_langenbuch) / sizeof(double));

  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds_rods_mov.prm";
  test_power_evolution(input_file, power_langen, 1e-4);

  // "Test Langenbuch reactor with Out-Of-Phase (sigma_f):"
  double power_langenbuch_oop[] =
                                    { 1, 1.06485, 1.02644 };
  std::vector<double> power_langen_oop(power_langenbuch_oop,
    power_langenbuch_oop
    + sizeof(power_langenbuch_oop) / sizeof(double));

  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds_oop.prm";
  test_power_evolution(input_file, power_langen_oop, 1e-4);

  // Test Langenbuch reactor Single material (All xsecs)

  double power_langenbuch_all[] =
                                    { 1, 0.862199 };
  std::vector<double> power_langen_all(power_langenbuch_all,
    power_langenbuch_all
    + sizeof(power_langenbuch_all) / sizeof(double));

  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds_all.prm";
  test_power_evolution(input_file, power_langen_all, 1e-4);

  // Test CROCUS - Mechanical vibration: << std::endl;
  double power_crocus[] =
                            { 1, 0.988385 };
  std::vector<double> power_crocus_v(power_crocus,
    power_crocus + sizeof(power_crocus) / sizeof(double));

  input_file = "test/2D_CROCUS/crocus_ds.prm";
  test_power_evolution(input_file, power_crocus_v, 1e-4);

  //
  // --------------------------------------------------------------- //
  // TEST REINIT
  double power_langenbuch_reinit[] =
                                       { 1, 1.06485, 1.15686 };
  std::vector<double> power_langen_reinit(power_langenbuch_reinit,
    power_langenbuch_reinit + sizeof(power_langenbuch_reinit) / sizeof(double));

  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds_reinit.prm";
  test_power_evolution(input_file, power_langen_reinit, 1e-4);

  double power_langenbuch_load[] =
                                     { 1.25088, 1.37291 };
  std::vector<double> power_langen_load(power_langenbuch_load,
    power_langenbuch_load + sizeof(power_langenbuch_load) / sizeof(double));

  input_file = "test/3D_Langenbuch/3D_Langenbuch_ds_reinit_load.prm";
  test_power_evolution(input_file, power_langen_load, 1e-4);

  //
  // --------------------------------------------------------------- //
  // TEST 1D_C5G7 7 ENERGY GROUPS
  // 7 energy group tests. Distributed TS. 1D_C5G7 Benchmark (3 cells)

  // All precursor data are constant and they are in .xsec
  double power_7_groups_S[] =
                                { 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000,
                                  1.0000000,
                                  1.0000000, 1.0000000, 1.0000000, 1.0000000,
                                  1.0000000 };
  std::vector<double> power_7_groups_s(power_7_groups_S,
    power_7_groups_S + sizeof(power_7_groups_S) / sizeof(double));

  //  Static cases:
  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_xsec_static_ix.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_xsec_static_six.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_xsec_static_sil.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  // Dynamic cases:

  double power_7_groups_d[] =
                                { 1.0000000002, 1.0138348297, 1.0287905141e+00,
                                  1.0442883214,
                                  1.0603234562, 1.0769205086, 1.0941071617,
                                  1.1119130335,
                                  1.1303697698, 1.1495112109, 1.1693735859 };
  std::vector<double> power_7_groups(power_7_groups_d,
    power_7_groups_d + sizeof(power_7_groups_d) / sizeof(double));

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_xsec_din_ix.prm";
  test_power_evolution(input_file, power_7_groups, 1e-8);

  double power_7_groups_Six[] =
                                  { 1.0000000005, 1.0137762653, 1.0286647111,
                                    1.0440912562,
                                    1.0600511877, 1.0765688164, 1.0936714975,
                                    1.1113884904,
                                    1.1297510494, 1.1487925888, 1.1685488691 };
  std::vector<double> power_7_groups_six(power_7_groups_Six,
    power_7_groups_Six + sizeof(power_7_groups_Six) / sizeof(double));

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_xsec_din_six.prm";
  test_power_evolution(input_file, power_7_groups_six, 1e-8);

  double power_7_groups_Sil[] =
                                  { 1.0000000005, 1.0137762654, 1.0286642713,
                                    1.0440898816,
                                    1.0600483613, 1.0765640017, 1.0936641359,
                                    1.1113779978,
                                    1.1297368142, 1.1487739687, 1.1685251872 };
  std::vector<double> power_7_groups_sil(power_7_groups_Sil,
    power_7_groups_Sil + sizeof(power_7_groups_Sil) / sizeof(double));

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_xsec_din_sil.prm";
  test_power_evolution(input_file, power_7_groups_sil, 1e-8);

  // ---------------------------------------------------------------------//
  // All precursor data are constant and they are in .prec
  // Static cases:
  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_prec_static_ix.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_prec_static_six.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_prec_static_sil.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  // Dynamic cases:
  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_prec_din_ix.prm";
  test_power_evolution(input_file, power_7_groups, 1e-8);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_prec_din_six.prm";
  test_power_evolution(input_file, power_7_groups_six, 1e-8);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_all_cte_prec_din_sil.prm";
  test_power_evolution(input_file, power_7_groups_sil, 1e-8);

  // ---------------------------------------------------------------------//

  // "All precursor data are constant and they are in .prec"
  // Static cases

  input_file = "test/1D_C5G7_3cells/1D_C5G7_lambda_cte_prec_static_six.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  input_file = "test/1D_C5G7_3cells/1D_C5G7_lambda_cte_prec_static_sil.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  // Dynamic cases:
  double power_7_groups_Six_lc[] =
                                     { 1.0000000006, 1.0324616869, 1.0703276409,
                                       1.1115297783,
                                       1.1562287762, 1.2048465741, 1.2579011496,
                                       1.3160073514,
                                       1.3798977852, 1.4504517870, 1.5287335746 };
  std::vector<double> power_7_groups_six_lc(power_7_groups_Six_lc,
    power_7_groups_Six_lc
    + sizeof(power_7_groups_Six_lc) / sizeof(double));

  input_file = "test/1D_C5G7_3cells/1D_C5G7_lambda_cte_prec_din_six.prm";
  test_power_evolution(input_file, power_7_groups_six_lc, 1e-8);

  double power_7_groups_Sil_lc[] =
                                     { 1.0000000006, 1.0324616869, 1.0703264191,
                                       1.1115257135,
                                       1.1562199633, 1.2048307899, 1.2578757852,
                                       1.3159693200,
                                       1.3798434057, 1.4503766373, 1.5286322907 };
  std::vector<double> power_7_groups_sil_lc(power_7_groups_Sil_lc,
    power_7_groups_Sil_lc
    + sizeof(power_7_groups_Sil_lc) / sizeof(double));

  input_file = "test/1D_C5G7_3cells/1D_C5G7_lambda_cte_prec_din_sil.prm";
  test_power_evolution(input_file, power_7_groups_sil_lc, 1e-8);

  double power_7_groups_Fi_lc[] =
                                    { 1.0000000006, 1.0326097328, 1.0706703798,
                                      1.1120983123,
                                      1.1570566742, 1.2059735264, 1.2593748343,
                                      1.3178854675,
                                      1.3822506645, 1.4533658287, 1.5323158197 };
  std::vector<double> power_7_groups_fi_lc(power_7_groups_Fi_lc,
    power_7_groups_Fi_lc
    + sizeof(power_7_groups_Fi_lc) / sizeof(double));

  // ---------------------------------------------------------------------//
  // "All precursor data change with the material
  // Static cases:
  input_file = "test/1D_C5G7_3cells/1D_C5G7_static_sil.prm";
  test_power_evolution(input_file, power_7_groups_s, 1e-4);

  // Dynamic cases:
  double power_7_groups_Sil_nc[] =
                                     { 1.0000000004, 1.0324616866, 1.0703249149,
                                       1.1115206784,
                                       1.1562089786, 1.2048109939, 1.2578437785,
                                       1.3159210390,
                                       1.3797739639, 1.4502801178, 1.5285014788 };
  std::vector<double> power_7_groups_sil_nc(power_7_groups_Sil_nc,
    power_7_groups_Sil_nc
    + sizeof(power_7_groups_Sil_nc) / sizeof(double));

  input_file = "test/1D_C5G7_3cells/1D_C5G7_din_sil.prm";
  test_power_evolution(input_file, power_7_groups_sil_nc, 1e-8);

  double power_7_groups_Fi_nc[] =
                                    { 1.0000000001, 1.0326082639, 1.0706654713,
                                      1.1120876210,
                                      1.1570374359, 1.2059424682, 1.2593280577,
                                      1.3178183043,
                                      1.3821574851, 1.4532397975, 1.5321485764 };
  std::vector<double> power_7_groups_fi_nc(power_7_groups_Fi_nc,
    power_7_groups_Fi_nc
    + sizeof(power_7_groups_Fi_nc) / sizeof(double));

  // ---------------------------------------------------------------------//

// TODO

  std::cout << std::endl;
  std::cout << "  ALL TEST PASSED!" << std::endl;
  std::cout << std::endl;
  return 1;
}

/**
 * @brief Run a time dependent problem.
 */
void test_power_evolution (const std::string &input_file,
  const std::vector<double> &ref_power,
  const double tol)
{
  std::cout << "Trans ... " << input_file << " " << std::flush;
  const bool verb = false;
  const bool sil = true;
  ParameterHandler prm;
  prm_declare_entries(prm);
  AssertRelease(fexists(input_file),
    "ERROR!: Input File .prm Does NOT exist\n Use -f input_file.prm  option");
  prm.parse_input(input_file);

  int dim = prm.get_integer("Dimension");
  int fe_degree = prm.get_integer("FE_Degree");
  std::vector<double> test_power;

  std::string distributed_time_scheme = "implicit-exponential";
  distributed_time_scheme = prm.get("Distributed_Time_Scheme");
  lower_case(distributed_time_scheme);
  get_string_from_options("-time_scheme", distributed_time_scheme);

  if (dim == 1)
  {
    if (fe_degree == 1)
    {
      StaticDiffusion<1, 1> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<1, 1> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;

    }
    else if (fe_degree == 2)
    {
      StaticDiffusion<1, 2> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<1, 2> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;

    }
    else if (fe_degree == 3)
    {
      StaticDiffusion<1, 3> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<1, 3> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;

    }
    else if (fe_degree == 4)
    {
      StaticDiffusion<1, 4> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<1, 4> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;

    }
    else if (fe_degree == 5)
    {
      StaticDiffusion<1, 5> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<1, 5> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
  }
  else if (dim == 2)
  {
    if (fe_degree == 1)
    {
      StaticDiffusion<2, 1> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<2, 1> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 2)
    {
      StaticDiffusion<2, 2> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<2, 2> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 3)
    {
      StaticDiffusion<2, 3> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<2, 3> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;

    }
    else if (fe_degree == 4)
    {
      StaticDiffusion<2, 4> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<2, 4> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 5)
    {
      StaticDiffusion<2, 5> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<2, 5> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
  }
  else if (dim == 3)
  {
    if (fe_degree == 1)
    {
      StaticDiffusion<3, 1> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<3, 1> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 2)
    {
      StaticDiffusion<3, 2> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<3, 2> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 3)
    {
      StaticDiffusion<3, 3> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<3, 3> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 4)
    {
      StaticDiffusion<3, 4> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<3, 4> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
    else if (fe_degree == 5)
    {
      StaticDiffusion<3, 5> static_prob(prm, input_file, verb, sil, false);
      TimeNeutronDiffusion<3, 5> time_pro(prm, static_prob, verb, sil);
      test_power = time_pro.power_vector;
    }
  }

  AssertRelease(
    is_similar(test_power[test_power.size() - 1],
      ref_power[ref_power.size() - 1], tol),
    "  Error solving test " + input_file + '\n' + "  Last power:  "
    + num_to_str(test_power[test_power.size() - 1])
    + " and it should be "
    + num_to_str(ref_power[ref_power.size() - 1])
    + '.');
  //  assert_vectors_similar(test_power, ref_power, tol);
  std::cout << " Passed!" << std::endl;
}

/**
 * @brief Run an Static test.
 */
void test_keff_problem (const std::string &input_file,
  const double reference_keff,
  const double tol)
{
  std::cout << "Static... " << input_file << " " << std::flush;

  ParameterHandler prm;
  prm_declare_entries(prm);
  AssertRelease(fexists(input_file),
    "ERROR!: Input File .prm Does NOT exist\n Use -f input_file.rm  option");
  prm.parse_input(input_file);

  int dim = prm.get_integer("Dimension");
  int fe_degree = prm.get_integer("FE_Degree");
  double calculated_keff = 0.0;
  std::string transport = prm.get("Transport_Appr");
  lower_case(transport);

  // The template parameters must be a constant:
  // This way, we ensure the compilation of all the cases
  // and only create/run the good one
  if (transport == "diffusion")
  {
    if (dim == 1)
    {
      if (fe_degree == 1)
      {
        StaticDiffusion<1, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticDiffusion<1, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticDiffusion<1, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticDiffusion<1, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticDiffusion<1, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
    else if (dim == 2)
    {
      if (fe_degree == 1)
      {
        StaticDiffusion<2, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticDiffusion<2, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticDiffusion<2, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticDiffusion<2, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticDiffusion<2, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
    else if (dim == 3)
    {
      if (fe_degree == 1)
      {
        StaticDiffusion<3, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticDiffusion<3, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticDiffusion<3, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticDiffusion<3, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticDiffusion<3, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
  }
  /* ------------------------------------------------------ */
  /* -----------------------    SPN   --------------------- */
  /* ------------------------------------------------------ */
  else if (transport == "spn")
  {
    if (dim == 1)
    {
      if (fe_degree == 1)
      {
        StaticSPN<1, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticSPN<1, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticSPN<1, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticSPN<1, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticSPN<1, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
    else if (dim == 2)
    {
      if (fe_degree == 1)
      {
        StaticSPN<2, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticSPN<2, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticSPN<2, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticSPN<2, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticSPN<2, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
    else if (dim == 3)
    {
      if (fe_degree == 1)
      {
        StaticSPN<3, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticSPN<3, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticSPN<3, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticSPN<3, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];

      }
      else if (fe_degree == 5)
      {
        StaticSPN<3, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
  }
  /* ------------------------------------------------------ */
  /* --------------------- FUll  SPN   -------------------- */
  /* ------------------------------------------------------ */
  else if (transport == "spn")
  {
    if (dim == 1)
    {
      if (fe_degree == 1)
      {
        StaticFullSPN<1, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticFullSPN<1, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticFullSPN<1, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticFullSPN<1, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticFullSPN<1, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
    else if (dim == 2)
    {
      if (fe_degree == 1)
      {
        StaticFullSPN<2, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticFullSPN<2, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticFullSPN<2, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticFullSPN<2, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 5)
      {
        StaticFullSPN<2, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
    else if (dim == 3)
    {
      if (fe_degree == 1)
      {
        StaticFullSPN<3, 1> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 2)
      {
        StaticFullSPN<3, 2> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 3)
      {
        StaticFullSPN<3, 3> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
      else if (fe_degree == 4)
      {
        StaticFullSPN<3, 4> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];

      }
      else if (fe_degree == 5)
      {
        StaticFullSPN<3, 5> problem(prm, input_file, false, true);
        calculated_keff = problem.eigenvalues[0];
      }
    }
  }
  else
  {
    AssertRelease(false, "Transport Mode not found!");
  }

  AssertRelease(is_similar(calculated_keff, reference_keff, tol),
    "  Error solving test " + input_file + '\n' + "  Calculated k_eff "
    + num_to_str(calculated_keff)
    + " and it should be "
    + num_to_str(reference_keff)
    + '.');

  std::cout << " Passed!" << std::endl;
}

/**
 *
 */
void chek_output_vector (const std::string &out_file,
  const std::string &headline,
  const std::vector<double> &vector_correct,
  const unsigned int n_lines,
  const double tol)
{
  std::cout << "Output... " << out_file << "  " << headline << std::flush;
  std::vector<double> vector_test;

  parse_vector_in_file(out_file, headline, vector_test, n_lines,
    vector_correct.size());

  assert_vectors_similar(vector_test, vector_correct, tol);
  std::cout << " Passed!" << std::endl;
}

/**
 *  @brief Test noise problem through the mean and the l2 norm of the delta_phi.
 */
void test_noise_problem (
  const std::string &input_file,
  const double ref_l2norm_delta_phi,
  const double tol)
{
  std::cout << "Noise ... " << input_file << std::flush;

  ParameterHandler prm;
  prm_declare_entries(prm);
  AssertRelease(fexists(input_file),
    "ERROR!: Input File .prm Does NOT exist\n Use -f input_file.rm  option");
  prm.parse_input(input_file);

  int dim = prm.get_integer("Dimension");
  int fe_degree = prm.get_integer("FE_Degree");
  double calculated_l2norm = 0.0;
  std::string transport = prm.get("Transport_Appr");
  lower_case(transport);
  bool noise = prm.get_bool("Noise_Calculation");
  AssertRelease(noise == true, "Must be a Noise calculation");
  const bool verbose = false;
  const bool silent = true;

  // The template parameters must be a constant:
  // This way, we ensure the compilation of the 3 problems
  // and only create/run the good one
  if (transport == "diffusion")
  {
    if (dim == 1)
    {
      if (fe_degree == 1)
      {
        StaticDiffusion<1, 1> static_prob(prm, input_file, false, true);
        NoiseDiffusion<1, 1> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 2)
      {
        StaticDiffusion<1, 2> static_prob(prm, input_file, false, true);
        NoiseDiffusion<1, 2> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 3)
      {
        StaticDiffusion<1, 3> static_prob(prm, input_file, false, true);
        NoiseDiffusion<1, 3> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 4)
      {
        StaticDiffusion<1, 4> static_prob(prm, input_file, false, true);
        NoiseDiffusion<1, 4> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 5)
      {
        StaticDiffusion<1, 5> static_prob(prm, input_file, false, true);
        NoiseDiffusion<1, 5> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
    }
    if (dim == 2)
    {
      if (fe_degree == 1)
      {
        StaticDiffusion<2, 1> static_prob(prm, input_file, false, true);
        NoiseDiffusion<2, 1> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 2)
      {
        StaticDiffusion<2, 2> static_prob(prm, input_file, false, true);
        NoiseDiffusion<2, 2> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 3)
      {
        StaticDiffusion<2, 3> static_prob(prm, input_file, false, true);
        NoiseDiffusion<2, 3> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 4)
      {
        StaticDiffusion<2, 4> static_prob(prm, input_file, false, true);
        NoiseDiffusion<2, 4> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 5)
      {
        StaticDiffusion<2, 5> static_prob(prm, input_file, false, true);
        NoiseDiffusion<2, 5> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
    }
    if (dim == 3)
    {
      if (fe_degree == 1)
      {
        StaticDiffusion<3, 1> static_prob(prm, input_file, false, true);
        NoiseDiffusion<3, 1> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 2)
      {
        StaticDiffusion<3, 2> static_prob(prm, input_file, false, true);
        NoiseDiffusion<3, 2> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 3)
      {
        StaticDiffusion<3, 3> static_prob(prm, input_file, false, true);
        NoiseDiffusion<3, 3> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 4)
      {
        StaticDiffusion<3, 4> static_prob(prm, input_file, false, true);
        NoiseDiffusion<3, 4> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      if (fe_degree == 5)
      {
        StaticDiffusion<3, 5> static_prob(prm, input_file, false, true);
        NoiseDiffusion<3, 5> noise_prob(prm, static_prob, false, true);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
    }
  }
  if (transport == "spn")
  {
    if (dim == 1)
    {
      if (fe_degree == 1)
      {
        StaticSPN<1, 1> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<1, 1> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();

      }
      else if (fe_degree == 2)
      {
        StaticSPN<1, 2> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<1, 2> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 3)
      {
        StaticSPN<1, 3> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<1, 3> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 4)
      {
        StaticSPN<1, 4> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<1, 4> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 5)
      {
        StaticSPN<1, 5> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<1, 5> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
    }
    else if (dim == 2)
    {
      if (fe_degree == 1)
      {
        StaticSPN<2, 1> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<2, 1> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 2)
      {
        StaticSPN<2, 2> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<2, 2> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 3)
      {
        StaticSPN<2, 3> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<2, 3> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 4)
      {
        StaticSPN<2, 4> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<2, 4> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 5)
      {
        StaticSPN<2, 5> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<2, 5> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
    }
    else if (dim == 3)
    {
      if (fe_degree == 1)
      {
        StaticSPN<3, 1> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<3, 1> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 2)
      {
        StaticSPN<3, 2> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<3, 2> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 3)
      {
        StaticSPN<3, 3> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<3, 3> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 4)
      {
        StaticSPN<3, 4> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<3, 4> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
      else if (fe_degree == 5)
      {
        StaticSPN<3, 5> static_prob(prm, input_file, verbose, silent);
        NoiseSPN<3, 5> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_u.l2_norm();
      }
    }
  }
  if (transport == "full_spn")
  {
    if (dim == 1)
    {
      if (fe_degree == 1)
      {
        StaticFullSPN<1, 1> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<1, 1> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();

      }
      else if (fe_degree == 2)
      {
        StaticFullSPN<1, 2> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<1, 2> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();

      }
      else if (fe_degree == 3)
      {
        StaticFullSPN<1, 3> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<1, 3> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 4)
      {
        StaticFullSPN<1, 4> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<1, 4> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 5)
      {
        StaticFullSPN<1, 5> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<1, 5> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
    }
    else if (dim == 2)
    {
      if (fe_degree == 1)
      {
        StaticFullSPN<2, 1> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<2, 1> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 2)
      {
        StaticFullSPN<2, 2> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<2, 2> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();

      }
      else if (fe_degree == 3)
      {
        StaticFullSPN<2, 3> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<2, 3> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 4)
      {
        StaticFullSPN<2, 4> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<2, 4> noise_prob(prm, static_prob, verbose, silent);
      }
      else if (fe_degree == 5)
      {
        StaticFullSPN<2, 5> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<2, 5> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
    }
    else if (dim == 3)
    {
      if (fe_degree == 1)
      {
        StaticFullSPN<3, 1> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<3, 1> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 2)
      {
        StaticFullSPN<3, 2> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<3, 2> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 3)
      {
        StaticFullSPN<3, 3> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<3, 3> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 4)
      {
        StaticFullSPN<3, 4> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<3, 4> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
      else if (fe_degree == 5)
      {
        StaticFullSPN<3, 5> static_prob(prm, input_file, verbose, silent);
        NoiseFullSPN<3, 5> noise_prob(prm, static_prob, verbose, silent);
        calculated_l2norm = noise_prob.delta_phi.l2_norm();
      }
    }
  }

  AssertRelease(calculated_l2norm != 0.0, "Test not done!");

  AssertRelease(is_similar(calculated_l2norm, ref_l2norm_delta_phi, tol),
    "  Error solving test " + input_file + '\n'
    + "  Calculated l2norm "
    + num_to_str(calculated_l2norm)
    + " and it should be "
    + num_to_str(ref_l2norm_delta_phi) + '.');

  std::cout << " Passed!" << std::endl;
}

/**
 * @brief test_sum_vector()
 */
void test_sum_vector ()
{
  std::vector<double> vector(2);
  vector[0] = 1.5;
  vector[1] = 2.5;
  AssertRelease(sum_vector(vector) == 4.0, "Error in sum_vector()");
}

/**
 * @brief test_compute_eigenvalues()
 */
void test_compute_eigenvalues ()
{
  double eigenvalue;
  const int n = 3;
  std::vector<double> eigenvector(n);
  LAPACKFullMatrix<double> matrix(n);
  matrix(0, 0) = 1.0;
  matrix(0, 1) = 1.0;
  matrix(0, 2) = 2.0;
  matrix(1, 0) = 1.0;
  matrix(1, 1) = 1.0;
  matrix(1, 2) = 2.0;
  matrix(2, 0) = 1.0;
  matrix(2, 1) = 1.0;
  matrix(2, 2) = 3.0;
  compute_max_eig(matrix, eigenvalue, eigenvector);
  AssertRelease(is_similar(eigenvalue, 4.56155281280882),
    "Error 1 in compute_eigenvalues");
  AssertRelease(eigenvector.size() == n, "Error 2 in compute_eigenvalues");
  AssertRelease(is_similar(eigenvector[0], -0.524114),
    "Error 3 in compute_eigenvalues");
  AssertRelease(is_similar(eigenvector[1], -0.524114),
    "Error 4 in compute_eigenvalues");
  AssertRelease(is_similar(eigenvector[2], -0.671273),
    "Error 5 in compute_eigenvalues");
}

/**
 * @brief test_lower_case()
 */
void test_lower_case ()
{
  AssertRelease(lower_case("WORD") == "word", "Error in Test 1");
  AssertRelease(lower_case("word") == "word", "Error in Test 2");
  AssertRelease(lower_case("WoRd") == "word", "Error in Test 3");
  AssertRelease(lower_case("123") == "123", "Error in Test 4");
  AssertRelease(lower_case("Num123NUM") == "num123num", "Error in Test 5");
}

/**
 * @brief test_copy_to
 */
void test_copy_to ()
{
  PETScWrappers::MPI::BlockVector vec1;
  const unsigned int n_blocks = 3;
  const unsigned int n_size_per_block = 2;
  const unsigned int n_size = n_blocks * n_size_per_block;
  vec1.reinit(n_blocks, PETSC_COMM_WORLD, n_size_per_block, n_size_per_block);

  Vec vec2, vec3;
  VecCreateMPI(PETSC_COMM_WORLD, n_size, n_size, &vec2);
  VecCreateMPI(PETSC_COMM_WORLD, n_size, n_size, &vec3);
  VecAssemblyBegin(vec2);
  for (int i = 0; i < static_cast<int>(n_size); ++i)
  {
    double val = static_cast<double>(i);
    VecSetValues(vec2, 1, &i, &val, INSERT_VALUES);
  }
  VecAssemblyEnd(vec2);

  copy_to_BlockVector(vec1, vec2);
  copy_to_Vec(vec3, vec1);
}

/**
 * @brief test_round
 */
void test_round ()
{
  double num = 1.23456789;
  AssertRelease(round(num) == 1.0, "Error in Test 1");
  AssertRelease(round(num, 0) == 1.0, "Error in Test 2");
  AssertRelease(round(num, 1) == 1.2, "Error in Test 3");
  AssertRelease(round(num, 2) == 1.23, "Error in Test 4");
}

/**
 *
 */
void test_child_pos ()
{
  std::vector<unsigned int> pos;
  std::vector<unsigned int> ref_pos;

  //  child_pos (const unsigned int child,
  //             const unsigned int n_levels,
  //             const unsigned int dim)
  // dim == 1
  const double tol = 1e-4;
  ref_pos =
              { 3 };
  pos = child_pos(3, 3, 1);
  assert_vectors_similar(pos, ref_pos, tol);

  // dim==2
  ref_pos =
              { 0, 0 };
  pos = child_pos(0, 1, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 1, 1 };
  pos = child_pos(3, 1, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 0, 1 };
  pos = child_pos(2, 1, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 2, 1 };
  pos = child_pos(6, 2, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 1, 3 };
  pos = child_pos(11, 2, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 5, 5 };
  pos = child_pos(51, 3, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 1, 6 };
  pos = child_pos(41, 3, 2);
  assert_vectors_similar(pos, ref_pos, tol);

  // dim==3
  ref_pos =
              { 1, 0, 0 };
  pos = child_pos(1, 1, 3);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 0, 1, 1 };
  pos = child_pos(6, 1, 3);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 3, 1, 0 };
  pos = child_pos(11, 2, 3);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 1, 3, 2 };
  pos = child_pos(51, 2, 3);
  assert_vectors_similar(pos, ref_pos, tol);

  ref_pos =
              { 3, 0, 2 };
  pos = child_pos(41, 2, 3);
  assert_vectors_similar(pos, ref_pos, tol);
}

/**
 *
 */
void test_child_at_face ()
{

  // child_at_face ( face, child, n_children, dim)
  // AssertRelease(child_at_face(face, child, n_children, dim) == false);

  // dim=1
  AssertRelease(child_at_face(1, 0, 1, 1) == true, "Failed Test 0");
  AssertRelease(child_at_face(0, 0, 2, 1) == true, "Failed Test 1");
  AssertRelease(child_at_face(1, 0, 8, 1) == false, "Failed Test 2");

  // dim=2
  AssertRelease(child_at_face(1, 56, 64, 2) == false, "Failed Test 3");
  AssertRelease(child_at_face(0, 0, 4, 2) == true, "Failed Test 4");
  AssertRelease(child_at_face(2, 20, 64, 2) == true, "Failed Test 5");
  AssertRelease(child_at_face(3, 20, 64, 2) == false, "Failed Test 6");
  AssertRelease(child_at_face(3, 46, 64, 2) == true, "Failed Test 7");

  // dim=3
  AssertRelease(child_at_face(0, 0, 8, 3) == true, "Failed Test 8");
  AssertRelease(child_at_face(3, 2, 8, 3) == true, "Failed Test 9");

  AssertRelease(child_at_face(5, 3, 8, 3) == false, "Failed Test 10");
  AssertRelease(child_at_face(4, 20, 64, 3) == false, "Failed Test 11");
  AssertRelease(child_at_face(3, 8, 64, 3) == false, "Failed Test 12");
  AssertRelease(child_at_face(1, 57, 64, 3) == true, "Failed Test 13");
  AssertRelease(child_at_face(5, 55, 64, 3) == true, "Failed Test 14");

}

/**
 * @brief run_tests_utils
 */
void run_tests_utils ()
{
  std::cout << "Unit Testing utils..." << std::flush;

  test_sum_vector();
  test_compute_eigenvalues();
  test_lower_case();
  test_round();
  test_copy_to();
  test_child_pos();
  test_child_at_face();

  std::cout << " Passed!" << std::endl;
}

///**
// * @brief run_test_static_rom
// */
//void run_test_static_rom (
//  const std::string &input_file,
//  unsigned int n_tests)
//{
//  std::cout << "Testing ROM Static.. " << input_file << " " << std::flush;
//
//  ParameterHandler prm;
//  prm_declare_entries(prm);
//  AssertRelease(fexists(input_file),
//    "ERROR!: Input File .prm Does NOT exist\n Use -f input_file.rm  option");
//  prm.parse_input(input_file);
//
//  int dim = prm.get_integer("Dimension");
//  AssertRelease(dim == 3,
//    "This test is only implemented for problems of dimension 3.");
//
//  int fe_degree = prm.get_integer("FE_Degree");
//  AssertRelease(fe_degree == 1,
//    "This test is only implemented for problems of FE Degree 1.");
//
//  std::string transport = prm.get("Transport_Appr");
//  lower_case(transport);
//  AssertRelease(transport == "diffusion",
//    "This test is only implemented for diffusion.");
//
//  std::string type_per = prm.get("Type_Perturbation");
//  AssertRelease(type_per == "Rods",
//    "This test is implemented for Bars Transient");
//
//  bool rom = prm.get_bool("ROM_Transient");
//  AssertRelease(rom, "This test is implemented for ROM Computation");
//
//  double rom_eig;
//  std::vector<double> static_eig(n_tests);
//  std::vector<double> rom_eigs(n_tests);
//
//  StaticDiffusion<3, 1> static_prob(prm, input_file, false, true);
//  ROMKinetics<3, 1> rom_prob(prm, static_prob, false, true, false);
//
//  const double time_step = rom_prob.t_end / n_tests;
//
//  rom_prob.snapshots.resize(1);
//  rom_prob.snapshots[0].resize(n_tests);
//  for (unsigned int nt = 0; nt < n_tests; nt++)
//  {
//    static_prob.perturbation.move_bars(nt * time_step);
//
//    // Compute eigenvalues with static_problem
//    static_prob.assemble_system_lambda();
//    static_prob.solve_eps();
//    static_eig[nt] = static_prob.eigenvalues[0];
//    static_prob.phi[0].compress(VectorOperation::insert);
//    rom_prob.snapshots[0][nt] = static_prob.phi[0];
//  }
//
//  for (unsigned int nt = 0; nt < n_tests; nt++)
//  {
//    LAPACKFullMatrix<double> romT(n_tests);
//    LAPACKFullMatrix<double> romTM(n_tests);
//    LAPACKFullMatrix<double> romM(n_tests);
//
//    static_prob.perturbation.move_bars(nt * time_step);
//
//    // Compute eigenvalues with ROM method
//    rom_prob.compute_pod_basis_LUPOD_monolithic(rom_prob.snapshots[0]);
//
//    // We use the matrices of the static problem because F in ROM is multiplied by (1-beta)
//    static_prob.assemble_system_lambda();
//
//    for (unsigned int b1 = 0; b1 < n_tests; b1++)
//      for (unsigned int b2 = 0; b2 < n_tests; b2++)
//      {
//        romT(b1, b2) = static_prob.T.vmult_dot(rom_prob.snap_basis[b1],
//          rom_prob.snap_basis[b2]);
//        romM(b1, b2) = static_prob.F.vmult_dot(rom_prob.snap_basis[b1],
//          rom_prob.snap_basis[b2]);
//      }
//
//    static_prob.T.clear();
//    static_prob.F.clear();
//
//    romT.invert();
//    romT.mmult(romTM, romM);
//    romTM.compute_eigenvalues(true);
//    for (unsigned int e = 0; e < n_tests; e++)
//      rom_eigs[e] = romTM.eigenvalue(e).real();
//    rom_eig = *std::max_element(rom_eigs.begin(), rom_eigs.end());
//
//    AssertRelease(std::abs(rom_eig - static_eig[nt]) < 1e-5,
//      "  Error solving test " + input_file + '\n'
//      + "  Calculated ROM "
//      + num_to_str(rom_eig)
//      + " and calculated Static "
//      + num_to_str(static_eig[nt])
//      + '.');
//  }
//
//  std::cout << " Passed!" << std::endl;
//}

/**
 * @brief run_test_rom_LUPOD
 * Test LUPOD Technique
 */
//void run_test_rom_LUPOD_1 ()
//{
//  std::string input_file = "test/1D_hom_5cells_mov_prec/1D_5cells_f1_rom_LUPOD.prm";
//  std::cout << "Testing ROM LUPOD 1... " << input_file << " " << std::flush;
//
//  ParameterHandler prm;
//  prm_declare_entries(prm);
//  AssertRelease(fexists(input_file),
//    "ERROR!: Input File .prm Does NOT exist\n Use -f input_file.rm  option");
//  prm.parse_input(input_file);
//
//  std::string transport = prm.get("Transport_Appr");
//  lower_case(transport);
//  AssertRelease(transport == "diffusion",
//    "This test is only implemented for diffusion.");
//
//  int fe_degree = prm.get_integer("FE_Degree");
//  AssertRelease(fe_degree == 1, "This test is implemented for FE_Degree==1");
//
//  bool rom = prm.get_bool("ROM_Transient");
//  AssertRelease(rom, "This test is implemented for ROM Computation");
//
//  StaticDiffusion<1, 1> static_prob(prm, input_file, false, true);
//  ROMKinetics<1, 1> rom_prob(prm, static_prob, false, true, true);
//  const double tol = 1e-5;
//  std::vector<unsigned int> points_reference =
//                                                 { 3, 1, 4 };
//  std::vector<unsigned int> snaps_reference =
//                                                { 3, 1, 4 };
//  //  std::cout << "rom_prob.snaps" << std::endl;
//  //  print_vector(rom_prob.snaps);
//  //  std::cout << " rom_prob.points" << std::endl;
//  //  print_vector(rom_prob.points);
//  //  std::cout << " rom_prob.epsilon_M " << rom_prob.epsilon_M << std::endl;
//  //  std::cout << " rom_prob.epsilon_N " << rom_prob.epsilon_N << std::endl;
//
//  assert_vectors_similar(rom_prob.snaps, snaps_reference, tol);
//  assert_vectors_similar(rom_prob.points, points_reference, tol);
//
//  PETScWrappers::MPI::Vector dst0;
//  std::vector<std::vector<double> >
//  U_red_ref =
//                {
//                    { -7.554602e-01, +4.166607e-01, -5.056419e-01 },
//                    { -4.931731e-01, -8.696971e-01, +2.018046e-02 },
//                    { -4.313468e-01, +2.646145e-01, +8.625074e-01 }
//                };
//
//  // Test Elements
//  for (unsigned int i = 0; i < rom_prob.snap_basis_red.size(); ++i)
//    for (unsigned int j = 0; j < rom_prob.snap_basis_red[0].size(); ++j)
//    {
//      AssertRelease(is_similar(U_red_ref[j][i], rom_prob.snap_basis_red[i][j], tol),
//        "Error in U_red_ref[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
//    }
//
//  // singular_values = { 1.52834, 0.587737, 0.0365, 3.23611e-16}
//
//  // Test U_full
//  std::vector<std::vector<double> >
//  U_full_ref =
//                 {
//                     { +1.775534e-17, -9.415731e-17, -2.753979e-15 },
//                     { -4.931731e-01, -8.696971e-01, +2.018046e-02 },
//                     { -7.366730e-01, -3.752936e-01, +8.965586e-01 },
//                     { -7.554602e-01, +4.166607e-01, -5.056419e-01 },
//                     { -4.313468e-01, +2.646145e-01, +8.625074e-01 },
//                     { +5.697307e-17, -1.716882e-16, -1.704035e-15 },
//                 };
//
//  // Test Elements
//  for (unsigned int i = 0; i < rom_prob.snap_basis.size(); ++i)
//    for (unsigned int j = 0; j < rom_prob.snap_basis[0].size(); ++j)
//    {
//      AssertRelease(is_similar(U_full_ref[j][i], rom_prob.snap_basis[i][j], tol),
//        "Error in U_full[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
//    }
//
//  std::cout << " Passed!" << std::endl;
//}
//
///**
// * @brief run_test_rom_LUPOD
// * Test LUPOD technique
// */
//void run_test_rom_LUPOD_2 ()
//{
//  std::string input_file = "test/1D_hom_5cells_mov_prec/1D_5cells_f1_rom_LUPOD2.prm";
//  std::cout << "Testing ROM LUPOD 2... " << input_file << " " << std::flush;
//
//  ParameterHandler prm;
//  prm_declare_entries(prm);
//  AssertRelease(fexists(input_file),
//    "ERROR!: Input File .prm Does NOT exist\n Use -f input_file.rm  option");
//  prm.parse_input(input_file);
//
//  std::string transport = prm.get("Transport_Appr");
//  lower_case(transport);
//  AssertRelease(transport == "diffusion",
//    "This test is only implemented for diffusion.");
//
//  int fe_degree = prm.get_integer("FE_Degree");
//  AssertRelease(fe_degree == 2, "This test is implemented for FE_Degree==1");
//
//  bool rom = prm.get_bool("ROM_Transient");
//  AssertRelease(rom, "This test is implemented for ROM Computation");
//
//  StaticDiffusion<1, 2> static_prob(prm, input_file, false, true);
//  ROMKinetics<1, 2> rom_prob(prm, static_prob, false, true, true);
//  const double tol = 1e-6;
//  std::vector<unsigned int> snaps_reference =
//                                                { 1, 2 };
//  std::vector<unsigned int> points_reference =
//                                                 { 4, 5 };
//
//  assert_vectors_similar(rom_prob.snaps, snaps_reference, tol);
//  assert_vectors_similar(rom_prob.points, points_reference, tol);
//  std::vector<std::vector<double> >
//  U_red_ref =
//                {
//                    { -0.729106, -0.684401 },
//                    { -0.684401, +0.729106 }
//                };
//
//  // Test Elements
//  for (unsigned int i = 0; i < rom_prob.snap_basis_red.size(); ++i)
//    for (unsigned int j = 0; j < rom_prob.snap_basis_red[0].size(); ++j)
//    {
//      AssertRelease(is_similar(U_red_ref[j][i], rom_prob.snap_basis_red[i][j], tol),
//        "Error in U_red[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
//    }
//
//  // Test U_full
//  std::vector<std::vector<double> > U_full_ref =
//                                                   {
//                                                       { -0.000000001266899,
//                                                         +0.000000001404289 },
//                                                       { -0.588889993842957,
//                                                         -0.762939586299597 },
//                                                       { -0.330700849029543,
//                                                         -0.481982160696203 },
//                                                       { -0.758143532925241,
//                                                         -0.159583583017522 },
//                                                       { -0.729105631228219,
//                                                         -0.684401182429795 },
//                                                       { -0.684401182429795,
//                                                         +0.729105631228219 },
//                                                       { -0.741655760998246,
//                                                         +0.443931888864489 },
//                                                       { -0.375852267581478,
//                                                         +0.439099210411635 },
//                                                       { -0.561536596313876,
//                                                         +0.710273354101464 },
//                                                       { -0.000000002552288,
//                                                         +0.000000002829072 },
//                                                       { -0.180285863103333,
//                                                         +0.163023661494686 },
//                                                   };
//
//  // Test Elements
//  for (unsigned int i = 0; i < rom_prob.snap_basis.size(); ++i)
//    for (unsigned int j = 0; j < rom_prob.snap_basis[0].size(); ++j)
//    {
//      AssertRelease(is_similar(U_full_ref[j][i], rom_prob.snap_basis[i][j], tol),
//        "Error in U_full[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
//    }
//  // singular_values = { 0.740887,: 0.337672}
//
//  // Test vmult_row
//  rom_prob.assemble_matrices();
//  PETScWrappers::MPI::BlockVector src(rom_prob.n_groups, rom_prob.comm, rom_prob.n_dofs,
//    rom_prob.n_dofs);
//  double dst = 0.0;
//  src = 1.0; // Set all vector entries to ones
//  rom_prob.F.vmult_row(dst, src, 3);
//
//  AssertRelease(is_similar(dst, 0.178074, tol), "  Error in vmult_row ");
//
//  std::cout << " Passed!" << std::endl;
//}
//
