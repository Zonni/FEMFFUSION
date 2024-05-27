/**
 * @file   perturbation.h
 * @brief
 */

#ifndef PERTURBATION_H_
#define PERTURBATION_H_

#include "materials.h"
#include "femffusion.h"

using namespace dealii;

template <int dim>
  class Perturbation
  {
    public:

    Perturbation (ParameterHandler &prm,
      Materials &materials,
      const DoFHandler<dim> &dof_handler,
      ConditionalOStream &verbose_cout);

    ParameterHandler &prm;
    Materials &materials;
    const DoFHandler<dim> &dof_handler;
    ConditionalOStream verbose_cout;

    std::vector<std::vector<double>> xsec_init;
    std::vector<std::vector<std::vector<double>>> xsec_init_all;

    // Instability Function
    std::vector<int> mat_changing;
//	int mat_changing_1;
    double frequency;
    double xs_amplitude;
    std::vector<double> amplitudes;
    double out_phase;
    std::vector<double> slope_up, slope_down, cut_time;
    std::string type_perturbation, tmodes, perturbation_function;
    int group_changing;

    // Data
    unsigned int n_groups;
    std::string xs_pert_name;

    // Geometry data
//  std::vector<unsigned int> assem_per_dim;

// TODO move_bar_flux_weighting
    std::vector<double> power_axial;
    std::vector<double> volume_per_plane;

    // Bars
    std::string rod_cusping_treat;
    unsigned int n_bars;
    double bars_top_pos;

    std::vector<std::vector<std::pair<double, double> > > bar_points;
    std::vector<unsigned int> bar_materials, bars_position;
    std::vector<unsigned int> materials_no_bars;

    // Mechanical vibration
    std::vector<double> vib_pos_static;
    unsigned int direction;
    unsigned int mat_vib;
    std::vector<unsigned int> indices_changed;
    std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_changed;
    std::vector<double> xs_initial;

    // Perturbation read_xs_file
    std::vector<std::vector<std::vector<double> > > perturbed_xs; // xs_perturbed[t][mat][xs]
    std::vector<double> perturbed_times;
    std::vector<unsigned int> perturbed_materials;

    std::vector<std::vector<double>> final_sigma_t, final_sigma_r, final_sigma_f,
        final_nu_sigma_f, final_sigma_tr;
    std::vector<std::vector<std::vector<double>>> final_sigma_s;

    void init_transient ();

    void get_parameters_from_command_line ();

    // Bars related
    void parse_bar_file (std::string BarFile);

    void move_bars_static ();

    void mechanical_vibration_static ();

    void move_read_xml_file (double sim_time,
      double t_end);

    void apply_c5G7_perturb (double sim_time);

    void move_bars (double sim_time);

    void move_bar_volume_homogenized (unsigned int bar_plant_pos,
      double bar_pos,
      unsigned int mat_bar,
      unsigned int bar);

    void move_bar_flux_weighting (unsigned int bar_plant_pos,
      double bar_pos,
      unsigned int mat_bar,
      unsigned int bar);

    void move_th (double sim_time);

    void step_change_material (double sim_time);

    void modify_xsec (double sim_time,
    		std::vector<int> mat_chan);

    void apply_function_to_perturb (double sim_time);

    void move_vibrating (double sim_time);

    void move_volume_homogenized (const std::vector<double> &vib_pos,
      const unsigned int &mat_bar,
      std::vector<unsigned int> &indices_changed,
      std::vector<typename DoFHandler<dim>::active_cell_iterator> &cells_changed);

    void move_read_xs_file (double sim_time);

    void get_read_xs_file (const std::string &xs_file);

    void read_xml_final_file (const std::string &xs_file);


    void restore_indices ();

  };

#endif /* PERTURBATION_H_ */
