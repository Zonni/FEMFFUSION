/*
 * perturbation.cc
 *
 *  Created on: 2 dic. 2020
 *      Author: amanda
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include "perturbation.h"
#include "materials.h"
#include "utils.h"

#include <algorithm>
#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace dealii;

/**
 *
 */
template <int dim>
  Perturbation<dim>::Perturbation (ParameterHandler &prm,
    Materials &materials,
    const DoFHandler<dim> &dh,
    ConditionalOStream &verbose_cout) :
      prm(prm),
      materials(materials),
      dof_handler(dh),
      verbose_cout(
        verbose_cout)
  {

    // 	Time parameters
    type_perturbation = prm.get("Type_Perturbation");
    perturbation_function = prm.get("Perturbation_Function");
    frequency = prm.get_double("Frequency");
    xs_amplitude = prm.get_double("Amplitude");
    out_phase = prm.get_double("Out_Phase");
    slope = prm.get_double("Slope");
    slope_up = prm.get_double("Slope_Up");
    slope_down = prm.get_double("Slope_Down");
    cut_time = prm.get_double("Cut_Time");
    mat_changing = prm.get_integer("Material_Changing") - 1;
    mat_changing_1 = prm.get_integer("Material_Changing_1") - 1;
    group_changing = prm.get_integer("Group_Changing") - 1;
    xs_pert_name = prm.get("XS_Name");
    lower_case(xs_pert_name);
    rod_cusping_treat = prm.get("Rod_Cusping_Method");

    if (type_perturbation == "AECL")
      materials_no_bars = materials.get_materials_vector();

    if ((perturbation_function == "Sinus" and xs_pert_name == "all") or perturbation_function
                                                                        == "Noise_7g")
    {
      parse_vector(prm.get("Amplitudes"), amplitudes);
    }

    if (type_perturbation == "AECL")
      xs_pert_name = "sigma_a";

    if (type_perturbation == "Out_Of_Phase")
      perturbation_function = "Sinus";

    // vibration parameters
    direction = prm.get_integer("Direction");
    AssertRelease(direction < dim, "Invalid Direction");
    if (type_perturbation == "Mechanical_Vibration")
      parse_vector(prm.get("Static_Position"), vib_pos_static, 2 * dim);
    mat_vib = prm.get_integer("Vibrating_Material") - 1;

    n_groups = materials.get_n_groups();
    materials.save_n_mats_init();

    // Bars
    n_bars = 0;
    bars_top_pos = 0.0;
    n_groups = 0;

  }

/*
 * This function initialize the perturbation class with the parameters
 */
template <int dim>
  void Perturbation<dim>::init_transient ()
  {

    n_groups = materials.get_n_groups();
    materials.save_initial_xsec();
    materials.save_n_mats_init();

    // Save the initial materials
    xsec_init.resize(n_groups);
    xsec_init_all.resize(4, std::vector<std::vector<double>>(n_groups));

    for (unsigned int ng = 0; ng < n_groups; ng++)
    {
      xsec_init[ng].resize(materials.get_n_mats());
      xsec_init_all[0][ng].resize(materials.get_n_mats());
      xsec_init_all[1][ng].resize(materials.get_n_mats());
      xsec_init_all[2][ng].resize(materials.get_n_mats());
      xsec_init_all[3][ng].resize(materials.get_n_mats());
      if (xs_pert_name == "sigma_f")
        xsec_init[ng] = materials.get_nu_sigma_f(ng);
      else if (xs_pert_name == "sigma_a")
        xsec_init[ng] = materials.get_sigma_r(ng);
      else if (xs_pert_name == "all")
      {
        xsec_init_all[0][ng] = materials.get_sigma_t(ng);
        xsec_init_all[1][ng] = materials.get_sigma_r(ng);
        xsec_init_all[2][ng] = materials.get_sigma_f(ng);
      }
    }

    if (xs_pert_name == "all")
    {
      AssertRelease(n_groups == 2,
        "This perturbation is only implemented for 2 energy groups");
      xsec_init_all[3][0] = materials.get_sigma_s(0, 1);
      amplitudes[2] /= materials.keff;
      amplitudes[6] /= materials.keff;
    }

    if (type_perturbation == "Read_XS_File")
    {
      verbose_cout << "  get_read_xs_file... " << std::flush;
      std::string read_xs_filename = prm.get("Read_XS_Filename");
      get_read_xs_file(read_xs_filename);
      verbose_cout << "   Done!" << std::endl;
    }

    get_parameters_from_command_line();
  }

template <int dim>
  void Perturbation<dim>::get_parameters_from_command_line ()
  {

    get_double_from_options("-frequency", frequency);
    get_double_from_options("-xs_amplitude", xs_amplitude);
    get_double_from_options("-out_phase", out_phase);
    get_double_from_options("-slope_up", slope_up);

    get_string_from_options("-rod_cusping_treat", rod_cusping_treat);

  }

/**
 * @brief This Functions parses the Bar file completing
 *  bars_position vector
 *  bar_materials
 *  bar_points
 *  bars_top_pos
 */
template <int dim>
  void Perturbation<dim>::parse_bar_file (std::string BarFile)
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
        bars_position.reserve(
          materials.n_assemblies / (materials.assem_per_dim[2])); // Reserve space for allocation
        parse_multiline_vector(input, materials.assem_per_dim[1],
          bars_position, true);

        if (verbose_cout.is_active())
          print_vector(bars_position);
      }
      else if (keyword == "Move_Banks")
      {
        verbose_cout << "parsing Move_Banks... ";

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
          Assert(bar == j, ExcMessage("Bars should be defined in order!"));

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
      else if (keyword == "RadialBarPosition")
      {
        int num, bank_pos, bank_num;
        iss >> num;
        bank_num = num;
        iss >> num;
        bank_pos = num;

        bars_position.resize(
          materials.n_assemblies / (materials.assem_per_dim[2]));

        bars_position[bank_pos - 1] = bank_num;
        if (verbose_cout.is_active())
          print_vector(bars_position);
      }
      else
        // Error
        Assert(false, ExcMessage("Invalid Header in Bar_file: " + keyword));
    }

    Assert(bars_top_pos > 0.0,
      ExcMessage("There is not specified a valid BarsTopPosition"));
    verbose_cout << "Done!" << std::endl;

    for (unsigned int bar = 0; bar < bar_points.size(); bar++)
      for (unsigned int p = 0; p < bar_points[bar].size(); p++)
        verbose_cout << " bar time " << bar_points[bar][p].first
                     << " height  "
                     << bar_points[bar][p].second << std::endl;

    materials_no_bars = materials.get_materials_vector();

  }

/**
 * @brief Move control bars.
 */
template <int dim>
  void Perturbation<dim>::move_bars_static ()
  {
    unsigned int bar;
    double bar_pos = 0;
    for (unsigned int plant_pos = 0; plant_pos < bars_position.size();
        ++plant_pos)
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
        std::cout << "  Bar " << bar + 1 << " in pos " << bar_pos << " cm "
                     << std::endl;

        move_bar_volume_homogenized(plant_pos, bar_pos, bar_materials[bar], bar);
      }
    }
  }

/**
 * @brief Move control bars.
 */
template <int dim>
  void Perturbation<dim>::mechanical_vibration_static ()
  {

    std::vector<unsigned int> indices_changed;
    std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_changed;
    if (frequency > 0.0 and xs_amplitude > 0.0)
      move_volume_homogenized(vib_pos_static, mat_vib, indices_changed,
        cells_changed);

  }

/**
 * @brief Move control bars.
 */
template <int dim>
  void Perturbation<dim>::move_bars (double sim_time)
  {
    double bar_pos_z;
    unsigned int bar;

    std::vector<double> bar_pos_z_vec(n_bars);
    // Compute the new position and direction of movement of all the bars.
    for (unsigned int bar = 0; bar < n_bars; ++bar)
    {
      for (unsigned int p = 1; p < bar_points[bar].size(); p++)
      {
        if (bar_points[bar][p].first - sim_time > -1e-8)
        {
          AssertRelease(
            bar_points[bar][p].first - bar_points[bar][p - 1].first
            > 1e-10,
            "The times that define the bar possition must be different");
          bar_pos_z_vec[bar] = bar_points[bar][p - 1].second
                               + (bar_points[bar][p].second
                                  - bar_points[bar][p - 1].second)
                                 / (bar_points[bar][p].first
                                    - bar_points[bar][p - 1].first)
                                 * (sim_time - bar_points[bar][p - 1].first);

          break;
        }
      }

      // If the simulation continues after the last time definition of the bar,
      // Stop the bar in the last position set.
      if (sim_time > bar_points[bar].back().first)
      {
        bar_pos_z_vec[bar] = bar_points[bar].back().second;
      }

      Assert(bar_points[bar][0].first == 0.0,
        ExcMessage("Error in initial time of the bars definition"));
      Assert(bar_pos_z_vec[bar] < 9e5,
        ExcMessage("Error in time of the bars definition"));
    }

// Move the bars  one by one (by )

    for (unsigned int plant_pos = 0; plant_pos < bars_position.size();
        ++plant_pos)
    {
      if (bars_position[plant_pos] > 0)
      {
        bar = bars_position[plant_pos] - 1;
        bar_pos_z = bar_pos_z_vec[bar];

        verbose_cout << "  Bar " << bar + 1 << " in pos " << bar_pos_z
                     << " cm "
                     << std::endl;

        if (rod_cusping_treat == "volhom")
          move_bar_volume_homogenized(plant_pos, bar_pos_z,
            bar_materials[bar], bar);
        else if (rod_cusping_treat == "fluxwei")
          move_bar_flux_weighting(plant_pos, bar_pos_z,
            bar_materials[bar], bar);
        else
          AssertRelease(false, "Invalid type of rod_cusping method");
      }
    }
  }

/**
 * @brief
 */
template <int dim>
  void Perturbation<dim>::move_bar_volume_homogenized (unsigned int plant_pos,
    double bar_pos,
    unsigned int mat_bar,
    unsigned int)
  {
    const unsigned int move_dim = dim - 1;
    unsigned int n_assemblies_per_plane = materials.n_assemblies
                                          / materials.assem_per_dim[move_dim];
    Assert(n_assemblies_per_plane > 0, ExcMessage("n_assemblies cannot be 0"));

    const double tol = 1e-18;

    double maxp = 0.0;
    double minp = 0.0;
    unsigned int mat_no_bar;
    double frac;
    std::vector<bool> is_done(materials.n_assemblies, false);
    unsigned int averaged_mat = materials.get_n_mats();
    std::cout << "plant_pos " << plant_pos << std::endl;
    std::cout << "mat_bar " << mat_bar << std::endl;

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
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
          if ((bar_pos - minp > tol) and (maxp - bar_pos > tol))
          {
            Assert(
              (bars_top_pos - minp > tol and maxp - bars_top_pos > tol) == false,
              ExcNotImplemented());

            // Calculate the fraction of the cell occupied by the bar
            frac = (maxp - bar_pos) / (maxp - minp);



            mat_no_bar = materials_no_bars[cell->user_index()];
            materials.create_new_mixed_mat(averaged_mat, frac, mat_bar,
              mat_no_bar, cell->user_index());

            std::cout << " frac " << frac << std::endl;
            std::cout << " mat_no_bar " << mat_no_bar << std::endl;
            std::cout << " Averaged R " << materials.get_sigma_r(7, averaged_mat) << std::endl;
            std::cout << " BAR R " << materials.get_sigma_r(7, mat_bar) << std::endl;
            std::cout << " No Bar R " << materials.get_sigma_r(7, mat_no_bar) << std::endl;

          }

          // The bar occupy all the cell
          else if (bar_pos - minp < 1e-8 && maxp - bars_top_pos < 1e-8)
          //          else if (bar_pos - minp < 1e-8 )
          {
//            if (maxp - bars_top_pos>1e-8) {
//              mat_no_bar = materials_no_bars[cell->user_index()];
//              materials.create_new_mixed_mat(averaged_mat,
//                0.0,
//                mat_bar,
//                mat_no_bar,
//                cell->user_index());
//              averaged_mat++;
//            }
//            else
            materials.set_materials_id(cell->user_index(), mat_bar);

          }
          // The bar does no occupy any space in the cell
          else
          {
            mat_no_bar = materials_no_bars[cell->user_index()];
            materials.set_materials_id(cell->user_index(), mat_no_bar); //
          }

        }

      }

    }

    is_done.clear();
  }

/**
 * @brief
 */
template <int dim>
  void Perturbation<dim>::move_bar_flux_weighting (unsigned int plant_pos,
    double bar_pos,
    unsigned int mat_bar,
    unsigned int)
  {

//	AssertRelease(false,
//			"We must obtain the vectors volume_per_plane and power_axial");

    const unsigned int move_dim = dim - 1;
    unsigned int n_assemblies_per_plane = materials.n_assemblies
                                          / materials.assem_per_dim[move_dim];
    Assert(n_assemblies_per_plane > 0, ExcMessage("n_assemblies cannot be 0"));

    const double tol = 1e-15;
    double maxp = 0.0;
    double minp = 0.0;
    unsigned int mat_no_bar;
    double frac;
    std::vector<bool> is_done(materials.n_assemblies, false);
    // Add the new material to the end
    unsigned int averaged_mat = materials.get_n_mats();
    double flux_bar = 1.0;
    double flux_nobar = 1.0;

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
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
//          std::cout<<"bar"<<bar<<"cell :"<<cell->user_index()<<std::endl;
          getMaxMinVertex(cell, move_dim, maxp, minp);
          // The bar_pos is in the middle of the cell.
          // We create a new material at the end of the defined materials
          // that has volume-averaged cross sections
          if ((bar_pos - minp > tol) and (maxp - bar_pos > tol))
          //          if ((bar_pos  > minp) and (maxp > bar_pos ) )
          {
            // Calculate the fraction of the cell occupied by the bar
            frac = (maxp - bar_pos) / (maxp - minp);

            unsigned int plane = cell->user_index()
                                 / n_assemblies_per_plane;

            // Compute flux_nobar
            flux_nobar = (volume_per_plane[plane - 1]
                          * power_axial[plane - 1]
                          + (1 - frac) * volume_per_plane[plane]
                            * power_axial[plane])
                         / (volume_per_plane[plane - 1]
                            + (1 - frac) * volume_per_plane[plane]);

            flux_bar = (volume_per_plane[plane + 1]
                        * power_axial[plane + 1]
                        + frac * volume_per_plane[plane]
                          * power_axial[plane])
                       / (volume_per_plane[plane + 1]
                          + frac * volume_per_plane[plane]);

            mat_no_bar = materials_no_bars[cell->user_index()];

            materials.create_new_mixed_mat_flux(averaged_mat, frac,
              mat_bar, mat_no_bar, flux_bar, flux_nobar,
              cell->user_index());

            averaged_mat++;

          }
          // The bar occupy all the cell
          else if (bar_pos - minp < tol && maxp - bars_top_pos < tol)
          //            else if ((bar_pos <= minp ) )
          {
            materials.set_materials_id(cell->user_index(), mat_bar);
          }
          else
          {
            mat_no_bar = materials_no_bars[cell->user_index()];
            materials.set_materials_id(cell->user_index(), mat_no_bar);

          }

        }

      }
    }

    is_done.clear();
  }

/**
 * @brief
 */
template <int dim>
  void Perturbation<dim>::move_th (double sim_time)
  {

    AssertRelease(materials.get_n_groups() == 2,
      "This is only implemented for 2 energy groups");

    std::vector<double> new_xsec(n_groups);
    for (int nmat = 0; nmat < 26; nmat++)
    {
      new_xsec[0] = xsec_init[0][nmat];
      new_xsec[1] = xsec_init[1][nmat];

      if (nmat == 4 or nmat == 5 or nmat == 9 or nmat == 10 or nmat == 16
          or nmat == 17
          or nmat == 21 or nmat == 22)
      {
        if (sim_time <= 0.4)
          new_xsec[1] += -1e-04 * sim_time;
        else if (sim_time > 0.4)
          new_xsec[1] += -1e-04 * 0.4 - 8.88889e-06 * (sim_time - 0.4);
      }

      materials.modify_xsec(xs_pert_name, nmat, new_xsec);
    }

    if (sim_time > 0.6)
    {

      double maxp = 0.0;
      double minp = 0.0;
      unsigned int mat_no_bar;
      double frac;
      unsigned int averaged_mat = materials.get_n_mats() - 1;
      double cte = 6.150e-04;

      double time_th = sim_time - 0.6;
      double bar_pos = 520 * time_th;

      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();
      for (cell = dof_handler.begin_active(); cell != endc; ++cell)
      {
        // We only update each assembly material once
        // However we must iterate over all cells in order to get the z position
        // of the cells

        unsigned mat_id = materials_no_bars[cell->user_index()];

        // This affect only to this cells
        if (mat_id == 1 or mat_id == 3 or mat_id == 6 or mat_id == 8
            or mat_id == 11
            or mat_id == 13 or mat_id == 15
            or mat_id == 17
            or mat_id == 18 or mat_id == 20
            or mat_id == 21
            or mat_id == 22 or mat_id == 23)
        {

          getMaxMinVertex(cell, 1, maxp, minp);

          // The bar_pos is in the middle of the cell.
          // We create a new material at the end of the defined materials
          // that has volume-averaged cross sections
          if ((bar_pos - minp > 1e-14) and (maxp - bar_pos > 1e-14))
          {
            averaged_mat++;
            // Calculate the fraction of the cell occupied by the bar
            frac = (maxp - bar_pos) / (maxp - minp);

            materials.create_new_added_mat(averaged_mat, frac, cte,
              mat_id, cell->user_index());

          }
          // The bar occupy all the cell
          else if (bar_pos - maxp > -1e-14 && maxp - 780 < 1e-16)
          {
            averaged_mat++;
            materials.create_new_added_mat(averaged_mat, 1.0, cte,
              mat_id, cell->user_index());
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
 * @brief Perturbate a xsec for C5G7-TD benchmmark
 */
template <int dim>
  void Perturbation<dim>::apply_c5G7_perturb (double sim_time)
  {
    materials.modify_xsec_c5g7_td11(sim_time);
  }


/**
 * @brief
 */
template <int dim>
  void Perturbation<dim>::step_change_material (double sim_time)
  {

    if (sim_time > 0.0)
      materials.change_mat_value(mat_changing, mat_changing_1);

  }

/**
 * @brief Perturbate a xsec with a function
 */
template <int dim>
  void Perturbation<dim>::apply_function_to_perturb (double sim_time)
  {

    if (type_perturbation == "Flux_Distributed")
    {
      for (int nmat = 0; nmat < static_cast<int>(materials.get_n_mats());
          nmat++)
        modify_xsec(sim_time, nmat);
    }
    else if (type_perturbation == "Single_Material")
    {
      modify_xsec(sim_time, mat_changing);
    }
    else if (type_perturbation == "Out_Of_Phase")
    {
      modify_xsec(sim_time, mat_changing, mat_changing_1);
    }
  }

/**
 * @brief Perturbate a xsec with a function
 */
template <int dim>
  void Perturbation<dim>::modify_xsec (double sim_time,
    unsigned int nmat,
    unsigned int nmat2)
  {

    std::vector<double> new_xsec(n_groups);
    std::vector<std::vector<double>> new_xsec_all(4,
      std::vector<double>(n_groups));

    // The xs_amplitudes in the previous inputs (12/2020) are in %
    if (perturbation_function == "Sinus")
    {

      if (group_changing > -1)
      {
        AssertRelease(xs_pert_name != "all",
          "This perturbation is only valid for sigmaf or sigmaa");
        for (unsigned int ng = 0; ng < n_groups; ng++)
          new_xsec[ng] = xsec_init[ng][nmat];

        new_xsec[group_changing] += xs_amplitude
                                    * sin(2 * M_PI * frequency * sim_time);

      }
      else
      {
        unsigned int xsec_am = 0;
        for (unsigned int ng = 0; ng < n_groups; ng++)
        {
          if (xs_pert_name == "sigma_f" or xs_pert_name == "sigma_a")
            new_xsec[ng] = xsec_init[ng][nmat]
                           + xs_amplitude
                             * sin(2 * M_PI * frequency * sim_time);
          else if (xs_pert_name == "all")
            for (unsigned xsec = 0; xsec < 4; xsec++)
            {
              new_xsec_all[xsec][ng] = xsec_init_all[xsec][ng][nmat]
                                       + amplitudes[xsec_am]
                                         * sin(2 * M_PI * frequency * sim_time);
              xsec_am++;
            }
        }
      }

    }
    else if (perturbation_function == "Constant")
    {
      for (unsigned int ng = 0; ng < n_groups; ng++)
      {
        if (sim_time > 0.0)
          new_xsec[ng] = xsec_init[ng][nmat] + 1e-4;
        else
          new_xsec[ng] = xsec_init[ng][nmat];
      }
    }
    else if (perturbation_function == "Ramp")
    {
      if (group_changing > -1)
      {
        for (unsigned int ng = 0; ng < n_groups; ng++)
          new_xsec[ng] = xsec_init[ng][nmat];

        if (sim_time < cut_time)
          new_xsec[group_changing] += xsec_init[group_changing][nmat]
                                      * (slope_up * sim_time);
        else
          new_xsec[group_changing] += xsec_init[group_changing][nmat]
                                      * (slope_up * cut_time)
                                      - (xsec_init[group_changing][nmat]
                                         + xsec_init[group_changing][nmat]
                                           * (slope_down * cut_time))
                                        * (slope_down * (sim_time - cut_time));
      }
      else
      {
        for (unsigned int ng = 0; ng < n_groups; ng++)
          if (sim_time < cut_time)
            new_xsec[ng] = xsec_init[ng][nmat]
                           * (1 + (slope_up * sim_time));
          else
            new_xsec[ng] = xsec_init[ng][nmat]
                           + xsec_init[ng][nmat] * (slope_up * cut_time)
                           - (xsec_init[ng][nmat]
                              + xsec_init[ng][nmat]
                                * (slope_down * cut_time))
                             * (slope_down * (sim_time - cut_time));
      }
    }
    else if (perturbation_function == "Ramp_hex")
    {
      new_xsec[0] = xsec_init[0][nmat];
      if (sim_time <= 1.0)
      {
        new_xsec[1] = 0.118870 * (1 - sim_time) + 0.016917 * sim_time;
      }
      else if ((sim_time > 1.0) and sim_time < 2.0)
      {
        new_xsec[1] = 0.118870 * (sim_time - 1) + 0.016917 * (2 - sim_time);
      }
      else
      {
        new_xsec[1] = 0.118870;
      }
    }
    else if (perturbation_function == "Noise_7g")
    {
      std::cout << "n_mat" << nmat << std::endl;
      materials.modify_xsec_7g(xs_pert_name, sim_time, amplitudes, nmat);
    }

    if (xs_pert_name == "sigma_f" or xs_pert_name == "sigma_a")
      materials.modify_xsec(xs_pert_name, nmat, new_xsec);
    else if (xs_pert_name == "all")
      materials.modify_xsec_all(xs_pert_name, nmat, new_xsec_all);

    if (nmat2 > 0
        and (xs_pert_name == "sigma_f" or xs_pert_name == "sigma_a"))
    {
      AssertRelease(perturbation_function == "Sinus",
        "This perturbation is only permited for the sinus instability");
      for (unsigned int ng = 0; ng < n_groups; ng++)
        new_xsec[ng] = xsec_init[ng][nmat2]
                       + xs_amplitude
                         * sin(2 * M_PI * frequency * sim_time + M_PI);
      materials.modify_xsec(xs_pert_name, nmat2, new_xsec);
    }

  }

/*
 * @brief Move Vibrating
 */
template <int dim>
  void Perturbation<dim>::move_vibrating (double sim_time)
  {
    //double vib_pos_min, vib_pos_max;
    std::vector<double> vib_pos = vib_pos_static;

    // Compute the new position and direction of movement of all the bars.
    vib_pos[2 * direction] = vib_pos_static[2 * direction]
                             + xs_amplitude * sin(2 * M_PI * frequency * sim_time);
    vib_pos[2 * direction + 1] = vib_pos_static[2 * direction + 1]
                                 + xs_amplitude * sin(2 * M_PI * frequency * sim_time);
    verbose_cout << "   Vibrating in "
                 << round(xs_amplitude * sin(2 * M_PI * frequency * sim_time) * 10,
                   4)
                 << " mm" << " assembly goes from "
                 << round(vib_pos[2 * direction], 4)
                 << " cm to "
                 << round(vib_pos[2 * direction + 1], 4)
                 << " cm." << std::endl;

    move_volume_homogenized(vib_pos, mat_vib, indices_changed, cells_changed);
    //move_volume_homogenized(vib_pos_min, vib_pos_max, direction, mat_vib);
    //std::cout << "materials_vector" << std::endl;
    //print_vector(materials.get_materials_vector());
  }

/**
 * @brief Move volume homogenized
 * For dim == 1
 */
template <int dim>
  void Perturbation<dim>::move_volume_homogenized (
    const std::vector<double> &vib_pos,
    const unsigned int &mat_bar,
    std::vector<unsigned int> &indices_changed,
    std::vector<typename DoFHandler<dim>::active_cell_iterator> &cells_changed)
  {

    AssertRelease(dim == 1, "This function only works for dim==1");
    double maxp = 0.0;
    double minp = 0.0;
    unsigned int mat_no_bar;
    ConditionalOStream cell_cout(std::cout, false);

    const unsigned int n_mats = materials.n_mats_init;
    const double eps = 1e-8;
    double frac;
    const int move_dim = 0; // vibration in x direction
    const double vib_pos_left = vib_pos[0];
    const double vib_pos_righ = vib_pos[1];
    unsigned int new_mat;
    unsigned int new_user_id = materials.n_assemblies;

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    {
      // plant_bar_pos = cell->user_index() % n_assemblies_per_plane;
      //if (bar_in_position[plant_bar_pos] - 1 == bar)
      //{
      cell_cout << "Cell " << cell << std::flush;

      getMaxMinVertex(cell, move_dim, maxp, minp);

      // The vib_pos_left is in the middle of the cell.
      // We create a new material at the end of the defined materials
      // that has volume-averaged cross sections
      if ((vib_pos_left - minp > eps) and (maxp - vib_pos_left > eps))
      {
        // It is not implemented when
        // vib_pos_left and vib_pos_right are in the middle of the same cell
        Assert(
          (vib_pos_righ - minp > eps and maxp - vib_pos_righ > eps) == false,
          ExcNotImplemented());

        // Calculate the fraction of the cell occupied by the bar
        frac = (maxp - vib_pos_left) / (maxp - minp);

        mat_no_bar = materials.get_original_material_id<dim>(cell);
        new_mat = n_mats;
        indices_changed.push_back(cell->user_index());
        cells_changed.push_back(cell);
        cell->set_user_index(new_user_id);
        materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
          mat_no_bar, cell->user_index());

        cell_cout << " Left in the Middle";
        cell_cout << " User Index: " << cell->user_index();
        cell_cout << " vib_pos_left: " << vib_pos_left;
        cell_cout << " minp: " << minp;
        cell_cout << " maxp: " << maxp;
        cell_cout << " frac " << frac;
        cell_cout << " new_mat " << new_mat;
        cell_cout << " mat_bar " << mat_bar;
        cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
      }
      // The vib_pos_right is in the middle of the cell.
      // We create a new material at the end of the defined materials
      // that has volume-averaged cross sections
      else if ((vib_pos_righ - minp > eps) and (maxp - vib_pos_righ > eps))
      {
        // Calculate the fraction of the cell occupied by the bar
        frac = (vib_pos_righ - minp) / (maxp - minp);

        mat_no_bar = materials.get_original_material_id<dim>(cell);
        new_mat = n_mats + 1;
        indices_changed.push_back(cell->user_index());
        cells_changed.push_back(cell);
        cell->set_user_index(new_user_id + 1);
        materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
          mat_no_bar, cell->user_index());

        cell_cout << " Right in the Middle";
        cell_cout << " User Index: " << cell->user_index();
        cell_cout << " vib_pos_righ: " << vib_pos_righ;
        cell_cout << " minp: " << minp;
        cell_cout << " maxp: " << maxp;
        cell_cout << " frac " << frac;
        cell_cout << " new_mat " << new_mat;
        cell_cout << " mat_bar " << mat_bar;
        cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
      }
      // The bar occupy all the cell
      else if (vib_pos_left - minp < eps and maxp - vib_pos_righ < eps)
      {
        materials.set_materials_id(cell->user_index(), mat_bar);

        cell_cout << " The bar is here ";
        cell_cout << " User Index: " << cell->user_index();
        cell_cout << " mat_bar " << mat_bar << std::endl;

      }
      // The bar does no occupy any space in the cell
      else
      {
        mat_no_bar = materials.get_original_material_id<dim>(cell);
        materials.set_materials_id(cell->user_index(), mat_no_bar);

        cell_cout << " The bar is NOT here ";
        cell_cout << " User Index: " << cell->user_index();
        cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
      }
    }
    cell_cout << " Number of cells changed " << indices_changed.size()
              << std::endl;
  }

/**
 * @brief
 */
template <>
  void Perturbation<2>::move_volume_homogenized (
    const std::vector<double> &vib_pos,
    const unsigned int &mat_bar,
    std::vector<unsigned int> &indices_changed,
    std::vector<typename DoFHandler<2>::active_cell_iterator> &cells_changed)
  {

    std::vector<double> maxp(2);
    std::vector<double> minp(2);
    unsigned int mat_no_bar;
    unsigned int new_mat = materials.n_mats_init;
    const double eps = 1e-8;
    double frac;
    unsigned int new_user_id = materials.n_assemblies;

    ConditionalOStream cell_cout(std::cout, false);
    typename DoFHandler<2>::active_cell_iterator cell =
                                                        dof_handler.begin_active(), endc =
        dof_handler.end();
    for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    {
      // plant_bar_pos = cell->user_index() % n_assemblies_per_plane;
      //if (bar_in_position[plant_bar_pos] - 1 == bar)
      //{
      cell_cout << "Cell " << cell << std::flush;

      for (unsigned int d = 0; d < 2; d++)
        getMaxMinVertex(cell, d, maxp[d], minp[d]);

      // The vib_pos_left is in the middle of the cell.
      // We create a new material at the end of the defined materials
      // that has volume-averaged cross sections
      if ((vib_pos[0] - minp[0] > eps) and (maxp[0] - vib_pos[0] > eps))
      {
        // It is not implemented when
        // vib_pos_left and vib_pos_right are in the middle of the same cell
        Assert(
          (vib_pos[0] - minp[0] > eps and maxp[0] - vib_pos[1] > eps) == false,
          ExcNotImplemented());

        // vib_pos[2] in the middle of the cell
        if ((vib_pos[2] - minp[1] > eps) and (maxp[1] - vib_pos[2] > eps))
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (maxp[0] - vib_pos[0]) / (maxp[0] - minp[0])
                 * (maxp[1] - vib_pos[2])
                 / (maxp[1] - minp[1]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Left in the Middle, Down in the Middle ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // vib_pos[3] in the middle of the cell
        else if ((vib_pos[3] - minp[1] > eps)
                 and (maxp[1] - vib_pos[3] > eps))
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (maxp[0] - vib_pos[0]) / (maxp[0] - minp[0])
                 * (vib_pos[3] - minp[1])
                 / (maxp[1] - minp[1]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Left in the Middle, Up in the Middle ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // The bar occupy all the cell in y
        else if (vib_pos[2] - minp[1] < eps
                 and maxp[1] - vib_pos[3] < eps)
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (maxp[0] - vib_pos[0]) / (maxp[0] - minp[0]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Left in the Middle, bar here";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // The bar does no occupy any space in the cell
        else
        {
          mat_no_bar = materials.get_original_material_id<2>(cell);
          materials.set_materials_id(cell->user_index(), mat_no_bar);

          cell_cout << " The bar is NOT here ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
        }
      }
      // The vib_pos_right is in the middle of the cell.
      // We create a new material at the end of the defined materials
      // that has volume-averaged cross sections
      else if ((vib_pos[1] - minp[0] > eps)
               and (maxp[0] - vib_pos[1] > eps))
      {
        // vib_pos[2] in the middle of the cell
        if ((vib_pos[2] - minp[1] > eps) and (maxp[1] - vib_pos[2] > eps))
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (maxp[0] - vib_pos[1]) / (maxp[0] - minp[0])
                 * (maxp[1] - vib_pos[2])
                 / (maxp[1] - minp[1]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Right in the Middle, Down in the Middle ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        else if ((vib_pos[3] - minp[1] > eps)
                 and (maxp[1] - vib_pos[3] > eps))
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (vib_pos[1] - minp[0]) / (maxp[0] - minp[0])
                 * (vib_pos[3] - minp[1])
                 / (maxp[1] - minp[1]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Right in the Middle, Up in the Middle ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // The bar occupy all the cell in y
        else if (vib_pos[2] - minp[1] < eps
                 and maxp[1] - vib_pos[3] < eps)
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (vib_pos[1] - minp[0]) / (maxp[0] - minp[0]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Right in the Middle, bar here";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // The bar does no occupy any space in the cell
        else
        {
          mat_no_bar = materials.get_original_material_id<2>(cell);
          materials.set_materials_id(cell->user_index(), mat_no_bar);

          cell_cout << " The bar is NOT here ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
        }
      }
      // The bar occupy all the cell in x direction
      else if (vib_pos[0] - minp[0] < eps and maxp[0] - vib_pos[1] < eps)
      {
        // vib_pos[2] in the middle of the cell
        if ((vib_pos[2] - minp[1] > eps) and (maxp[1] - vib_pos[2] > eps))
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (maxp[1] - vib_pos[2]) / (maxp[1] - minp[1]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Right in the Middle, Down in the Middle ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // vib_pos[3] in the middle of the cell
        else if ((vib_pos[3] - minp[1] > eps)
                 and (maxp[1] - vib_pos[3] > eps))
        {
          // Calculate the fraction of the cell occupied by the bar
          frac = (vib_pos[3] - minp[1]) / (maxp[1] - minp[1]);

          mat_no_bar = materials.get_original_material_id<2>(cell);
          indices_changed.push_back(cell->user_index());
          cells_changed.push_back(cell);
          cell->set_user_index(new_user_id);
          materials.create_new_mixed_mat_toni(new_mat, frac, mat_bar,
            mat_no_bar, cell->user_index());

          cell_cout << " Right in the Middle, Up in the Middle ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " frac " << frac;
          cell_cout << " new_mat " << new_mat;
          cell_cout << " mat_bar " << mat_bar;
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
          new_mat++;
          new_user_id++;
        }
        // The bar occupy all the cell in y
        else if (vib_pos[2] - minp[1] < eps
                 and maxp[1] - vib_pos[3] < eps)
        {
          materials.set_materials_id(cell->user_index(), mat_bar);

          cell_cout << " The bar is here ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " mat_bar " << mat_bar << std::endl;

        }
        // The bar does not occupy any space in the cell
        else
        {
          mat_no_bar = materials.get_original_material_id<2>(cell);
          materials.set_materials_id(cell->user_index(), mat_no_bar);

          cell_cout << " The bar is NOT here ";
          cell_cout << " User Index: " << cell->user_index();
          cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
        }
      }
      // The bar does not occupy any space in the cell
      else
      {
        mat_no_bar = materials.get_original_material_id<2>(cell);
        materials.set_materials_id(cell->user_index(), mat_no_bar);

        cell_cout << " The bar is NOT here ";
        cell_cout << " User Index: " << cell->user_index();
        cell_cout << " mat_no_bar " << mat_no_bar << std::endl;
      }
    }
  }

template void Perturbation<1>::move_volume_homogenized (
  const std::vector<double> &vib_pos,
  const unsigned int &mat_bar,
  std::vector<unsigned int> &indices_changed,
  std::vector<typename DoFHandler<1>::active_cell_iterator> &cells_changed);

template void Perturbation<2>::move_volume_homogenized (
  const std::vector<double> &vib_pos,
  const unsigned int &mat_bar,
  std::vector<unsigned int> &indices_changed,
  std::vector<typename DoFHandler<2>::active_cell_iterator> &cells_changed);

template void Perturbation<3>::move_volume_homogenized (
  const std::vector<double> &vib_pos,
  const unsigned int &mat_bar,
  std::vector<unsigned int> &indices_changed,
  std::vector<typename DoFHandler<3>::active_cell_iterator> &cells_changed);

/**
 * @brief
 */
template <int dim>
  void Perturbation<dim>::get_read_xs_file (const std::string &xs_file)
  {
    // perturbed_xs
    Assert(fexists(xs_file), ExcMessage("read_xs_file doesn't exist"));
    std::ifstream input(xs_file.c_str(), std::ios::in);
    std::string str, keyword;
    // unsigned int mat;

    unsigned int n_perturbed_materials, n_time_steps;
    double dob;
    unsigned int uint;
    const unsigned int n_xsecs = 9;

    // for every line
    for (std::string line; getline(input, line);)
    {
      std::istringstream iss(line);
      keyword.clear();
      iss >> keyword;

      if (is_commentary(keyword))
        continue;
      // First definition Material and XSecs:
      else if (keyword == "N_Perturbed_Materials")
      {
        verbose_cout << "    N_Perturbed_Materials " << std::flush;
        iss >> n_perturbed_materials;
        Assert(!iss.fail(),
          ExcMessage("It must be defined the number N_Perturbed_Materials!"));
        perturbed_xs.resize(n_perturbed_materials);
        perturbed_materials.resize(n_perturbed_materials);
        verbose_cout << n_perturbed_materials << "   Done! " << std::endl;
      }
      else if (keyword == "Perturbed_Materials")
      {
        verbose_cout << "    Perturbed_Materials " << std::flush;
        int number;
        for (unsigned int i = 0; i < n_perturbed_materials; i++)
        {
          iss >> number;
          Assert(!iss.fail(),
            ExcMessage("Enough perturbed_materials must be defined!"));
          perturbed_materials[i] = number - 1;
        }
        if (verbose_cout.is_active())
          print_vector(perturbed_materials, false);

        verbose_cout << " Done! " << std::endl;
      }
      else if (keyword == "N_Time_Steps")
      {
        verbose_cout << "    N_Time_Steps " << std::flush;
        iss >> n_time_steps;
        Assert(!iss.fail(),
          ExcMessage("It must be defined the number N_Time_Steps!"));
        verbose_cout << n_time_steps << "   Done! " << std::endl;
        perturbed_times.resize(n_time_steps);
      }
      else if (keyword == "XSECS")
      {
        verbose_cout << "  parsing XSecs..." << std::endl;

        // Resize perturbed_materials
        perturbed_xs.resize(n_time_steps);
        for (unsigned int t = 0; t < n_time_steps; t++)
        {
          perturbed_xs[t].resize(n_perturbed_materials);
          for (unsigned int mat = 0; mat < n_perturbed_materials; mat++)
            perturbed_xs[t][mat].resize(n_xsecs);
        }

        for (unsigned int t = 0; t < n_time_steps; t++)
        {
          str = get_new_valid_line(input, line);
          std::istringstream iss(line);
          iss >> str;
          AssertRelease(str == "Time", "It is expected a Time");

          iss >> dob;
          Assert(!iss.fail(),
            ExcMessage("It must be defined the number of materials defined!"));
          perturbed_times[t] = dob;
          verbose_cout << "      Time " << perturbed_times[t]
                       << std::endl;

          for (unsigned int mat = 0; mat < n_perturbed_materials; mat++)
          {
            str = get_new_valid_line(input, line);

            std::istringstream iss(line);
            iss >> uint;
            AssertRelease((uint - 1) == perturbed_materials[mat],
              "The line must begin by the perturbed material: "
              + num_to_str(uint));

            for (unsigned int xs = 0; xs < n_xsecs; xs++)
            {
              // sigma_tr1
              iss >> dob;
              perturbed_xs[t][mat][xs] = dob;
            }
            verbose_cout << "      ";
            if (verbose_cout.is_active())
              print_vector(perturbed_xs[t][mat]);
          }
        }
      }
    }
  }

/*
 * @brief Change the xs following the generic absorber of variable strength.
 */
template <int dim>
  void Perturbation<dim>::move_read_xs_file (double sim_time)
  {
    double frac;
    unsigned int pert_mat;
    double sigma_tr1, sigma_tr2, sigma_a1, sigma_a2;
    double nu_sigma_f1, nu_sigma_f2, sigma_f1, sigma_f2, sigma_12;
    for (unsigned int st = 0; st < perturbed_times.size(); st++)
      if (sim_time >= perturbed_times[st]
          and sim_time < perturbed_times[st + 1])
      {
        frac = (sim_time - perturbed_times[st])
               / (perturbed_times[st + 1] - perturbed_times[st]);
        verbose_cout << "Update XS...   t= " << sim_time << ", " << frac
                     << std::endl;
        for (unsigned int mat = 0; mat < perturbed_materials.size();
            mat++)
        {
          // Group 1
          pert_mat = perturbed_materials[mat];
          sigma_tr1 = (1.0 - frac) * perturbed_xs[st][mat][0]
                      + frac * perturbed_xs[st + 1][mat][0];
          sigma_a1 = (1.0 - frac) * perturbed_xs[st][mat][1]
                     + frac * perturbed_xs[st + 1][mat][1];
          nu_sigma_f1 = (1.0 - frac) * perturbed_xs[st][mat][2]
                        + frac * perturbed_xs[st + 1][mat][2];
          sigma_f1 = (1.0 - frac) * perturbed_xs[st][mat][3]
                     + frac * perturbed_xs[st + 1][mat][3];
          sigma_12 = (1.0 - frac) * perturbed_xs[st][mat][4]
                     + frac * perturbed_xs[st + 1][mat][4];
          sigma_tr2 = (1.0 - frac) * perturbed_xs[st][mat][5]
                      + frac * perturbed_xs[st + 1][mat][5];
          sigma_a2 = (1.0 - frac) * perturbed_xs[st][mat][6]
                     + frac * perturbed_xs[st + 1][mat][6];
          nu_sigma_f2 = (1.0 - frac) * perturbed_xs[st][mat][7]
                        + frac * perturbed_xs[st + 1][mat][7];
          sigma_f2 = (1.0 - frac) * perturbed_xs[st][mat][8]
                     + frac * perturbed_xs[st + 1][mat][8];

          materials.set_diffusion_coefficient(1.0 / (3 * sigma_tr1), 0,
            pert_mat);
          materials.set_diffusion_coefficient(1.0 / (3 * sigma_tr2), 1,
            pert_mat);
          materials.set_sigma_r(sigma_a1 + sigma_12, 0, pert_mat);
          materials.set_sigma_r(sigma_a2, 1, pert_mat);
          materials.set_nu_sigma_f(nu_sigma_f1 / materials.keff, 0,
            pert_mat);
          materials.set_nu_sigma_f(nu_sigma_f2 / materials.keff, 1,
            pert_mat);
          materials.set_sigma_f(sigma_f1, 0, pert_mat);
          materials.set_sigma_f(sigma_f2, 1, pert_mat);
          materials.set_sigma_s(sigma_12, 0, 1, pert_mat);
        }
      }
  }

/*
 * @brief Restore the user indices.
 */
template <int dim>
  void Perturbation<dim>::restore_indices ()
  {

    for (unsigned int i = 0; i < indices_changed.size(); i++)
    {
      cells_changed[i]->set_user_index(indices_changed[i]);
    }

    cells_changed.clear();
    indices_changed.clear();
  }

template class Perturbation<1> ;
template class Perturbation<2> ;
template class Perturbation<3> ;
