/**
 * @file   matrix_operators_sdpn.cc
 * @brief  Implementation of TransportMatrixSDPN and FissionMAtrix classes to handle block matrices.
 */

#include "../../include/matrix_operators/matrix_operators_sdpn.h"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h>

#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

std::vector<std::vector<double> > get_sdpn_u_to_phi_coeff (std::string &sdpn_type,
  unsigned int n_moments)
{

  const unsigned int max_moment = n_moments;
  std::vector<std::vector<double> > sdpn_u_to_phi_coeff(max_moment,
    std::vector<double>(max_moment));

  if (sdpn_type == "nazari")
  {
    if (n_moments == 2)
    {
      sdpn_u_to_phi_coeff =
                              {
                                  { 1, -3. / 4 },
                                  { 0, 3. / 8 },
                              };
    }
    if (n_moments == 3)
    {
      sdpn_u_to_phi_coeff =
                              {
                                  { 1, -5. / 144, +17. / 144 },
                                  { 0, 5. / 288, 17. / 288 },
                                  { 0, 1. / 64, -3. / 64 }
                              };
    }
    if (n_moments == 4)
    {
      sdpn_u_to_phi_coeff =
                              {
                                  { 1.0, -2. / 3, +145. / 5082, -73. / 1694 },
                                  { 0, 0, 0, 0 },
                                  { 0, 0, 0, 0 },
                                  { 0, 0, 0, 0 }
                              };
    }
  }
  else if (sdpn_type == "damian")
  {
    if (n_moments == 2)
    {
      sdpn_u_to_phi_coeff =
                              {
                                  { 1.0, -2. / 3. },
                                  { 0, 1. / 3 },
                              };
    }
    if (n_moments == 3)
    {
      sdpn_u_to_phi_coeff =
                              {
                                  { 1.0, -2. / 3, 8. / 15, },
                                  { 0, 0, 0 },
                                  { 0, 0, 0 }
                              };
    }
    if (n_moments == 4)
    {
      sdpn_u_to_phi_coeff =
                              {
                                  { 1.0, -2. / 3, 8. / 15, -16. / 35 },
                                  { 0, 0, 0, 0 },
                                  { 0, 0, 0, 0 },
                                  { 0, 0, 0, 0 }
                              };
    }
  }
  else
    AssertRelease(false, "Not valid sdpn type: " + sdpn_type);

  return sdpn_u_to_phi_coeff;
}

void define_sdpn_coeffs (std::string &sdpn_type,
  unsigned int n_moments)
{
  const unsigned int max_moment = n_moments;

  sdpn_coeff2.resize(max_moment,
    std::vector<std::vector<double> >(max_moment, std::vector<double>(max_moment)));
  sdpn_marshack_coeff2.resize(max_moment, std::vector<double>(max_moment));
  sdpn_diff_coeff.resize(max_moment);

  if (sdpn_type == "nazari")
  {
    if (n_moments == 2) // SDP1
    {
      sdpn_coeff2 =
                      {
                          {
                              { 1., -3. / 4 },
                              { -3. / 4, 9. / 16 },
                          },
                          {
                              { 0.0, 0.0 },
                              { 0.0, 3. / 4 },
                          }
                      };

      sdpn_marshack_coeff2 =
                               {
                                   { +1. / 2, -1. / 8 },
                                   { -1. / 8, +9. / 32 },
                               };

      sdpn_diff_coeff =
                          { 1. / 3, 1. / 16 };
    }
    if (n_moments == 3) // SDP2
    {
      sdpn_coeff2 =
                      {
                          {
                              { 1., -5. / 144, +17. / 144 },
                              { -5. / 144, +25. / (96. * 216), -(5. * 17) / (96. * 216) },
                              { +17. / 144, -85. / 20736, +289. / 20736 },
                          },
                          {
                              { 0., 0., 0. },
                              { 0., 5. * 25 / (96. * 864), -(5. * 85) / (96. * 864) },
                              { 0., -425. / 82944, 1445. / 82944 },
                          },
                          {
                              { 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 0.0 },
                              { 0.0, 255. / 31744, 765. / 31744 },
                          }
                      };

      sdpn_marshack_coeff2 =
                               {
                                   { +1. / 2, -1. / 96, +1. / 96 },
                                   { -5. / 768, +155. / 110592, -95. / 110592 },
                                   { 17. / 1488, -17. / 214272, +1343. / 214272 },
                               };

      sdpn_diff_coeff =
                          { 1. / 3, 5. / 6144, 17. / 23808 };
    }

    if (n_moments == 4) // SDP3
    {
      sdpn_coeff2 =
            {
                {
                    { 1, -(2. / 3), 145. / 5082, -(73. / 1694) },
                    { -(2. / 3), 4. / 9, -(145. / 7623), 73. / 2541 },
                    { -(3193625. / 76372296), 3193625. / 114558444, -(463075625.
                                                                      / 388124008272),
                      233134625. / 129374669424 },
                    { -(12685519. / 242746608), 12685519. / 364119912, -(2171665.
                                                                         / 1456479648),
                      14977. / 6650592 },
                },
                {
                    { 0, 0, 0, 0 },
                    { 0, 5. / 9, -(725. / 30492), 365. / 10164 },
                    { 0, 15968125. / 458233776, -(2315378125. / 1552496033088),
                      1165673125. / 517498677696 },
                    { 0, 63427595. / 1456479648, -(10858325. / 5825918592), 74885.
                                                                            / 26602368 },
                },
                {
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, -(55569075. / 36490291376), 83928465. / 36490291376 },
                    { 0, 0, -(2171665. / 1232405856), 14977. / 5627424 },

                },
                {
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 33341445. / 111617361856, -(139880775. / 504022149631) },
                    { 0, 0, -(612843863. / 106808507520), 73123. / 13716846 },
                }
            };

      sdpn_marshack_coeff2 =
                               {
                                   { 1. / 2, -1. / 8, +417. / 125000, -2. / 847 },
                                   { -1. / 8, +7. / 24, -1137. / 200000, +7. / 1936 },
                                   { +1. / 16, -41. / 384, +559. / 50000, -235. / 54208 },
                                   { -515. / 411136, 103. / 51392, -169. / 1000000,
                                     -215167. / 696464384 },
                               };

      sdpn_diff_coeff =
                          { 1. / 3, 1. / 7, 740921. / 3317299216, 881. / 15347520 };
    }
  }
  else if (sdpn_type == "damian")
  {
    if (n_moments == 2) //SDP1 Damian
    {
      sdpn_coeff2 =
                      {
                          {
                              { 1., -2. / 3 },
                              { -6. / 7., 4. / 7. },
                          },
                          {
                              { 0.0, 0.0 },
                              { 0.0, 5. / 7 },
                          }
                      };

      sdpn_marshack_coeff2 =
                               {
                                   { +1. / 2, -1. / 8 },
                                   { -2. / 7, +23. / 42 },
                               };

      sdpn_diff_coeff =
                          { 1. / 3, 1. / 7 };
    }
    if (n_moments == 3) // Damian SDP2
    {
      std::cout << "Damian SDP2 " << std::endl;
      sdpn_coeff2 =
                      {
                          {
                              { 1., -(2. / 3), 8. / 15 },
                              { -(2. / 3), 4. / 9, -(16. / 45) },
                              { 32. / 33, -(64. / 99), 256. / 495 },
                          },
                          {
                              { 0., 0., 0. },
                              { 0., 5. / 9, -(4. / 9) },
                              { 0., -(80. / 99), 64. / 99 },
                          },
                          {
                              { 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 36. / 55 },
                          },
                      };

      sdpn_marshack_coeff2 =
                               {
                                   { 1. / 2, -(1. / 8), 1. / 16 },
                                   { -2. / 21, (31. / 126), -59. / 1260 },
                                   { (8. / 33), -38. / 99, (287. / 495) },

                               };

      sdpn_diff_coeff =
                          { 1. / 3, 1. / 7, 1. / 11 };
    }
    if (n_moments == 4) // Damian SDP3
    {
      sdpn_coeff2 =
                      {
                          {
                              { 1., -(2. / 3), 8. / 15, -(16. / 35) },
                              { -(2. / 3), 4. / 9, -(16. / 45), 32. / 105 },
                              { 8. / 15, -(16. / 45), 64. / 225, -(128. / 525) },
                              { -(176. / 125), 352. / 375, -(1408. / 1875), 2816. / 4375 }
                          },
                          {
                              { 0., 0., 0., 0. },
                              { 0., 5. / 9, -(4. / 9), 8. / 21 },
                              { 0., -(4. / 9), 16. / 45, -(32. / 105) },
                              { 0., 88. / 75, -(352. / 375), 704. / 875 },
                          },
                          {
                              { 0.0, 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 9. / 25, -(54. / 175) },
                              { 0.0, 0.0, -(594. / 625), 3564. / 4375 },
                          },
                          {
                              { 0.0, 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 0.0, 143. / 175 },
                          }
                      };

      sdpn_marshack_coeff2 =
            {
                { 1. / 2, -(1. / 8), 1. / 16, -(5. / 128) },
                { -(1. / 8), 7. / 24, -(41. / 384), 1. / 16 },
                { -(2. / 143), 4. / 429, 1087. / 17160, 3611. / 40040 },
                { -(10. / 39), 46. / 117, -(2417. / 4680), 1891. / 2730 },
            };

      sdpn_diff_coeff =
                          { 1. / 3, 1. / 7, 1. / 11, 1. / 15 };
    }

  }
  else
    AssertRelease(false, "Not valid sdpn type: " + sdpn_type);
}

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  TransportMatrixSDPN<dim, n_fe_degree>::TransportMatrixSDPN (const MPI_Comm &comm,
    const DoFHandler<dim> &dof_handler,
    const AffineConstraints<double> &constraints) :
      TransportMatrixBase<dim, n_fe_degree>(comm, dof_handler, constraints),
      dof_handler(
        dof_handler),
      constraints(constraints)
  {
    // Silly initializations
    n_moments = 0;
    n_groups = 0;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixSDPN<dim, n_fe_degree>::reinit (const Materials &materials,
    const unsigned int _n_moments,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const MatrixFreeType &_matrixfree_type,
    bool listen_to_material_id)
  {
    this->matrixfree_type = _matrixfree_type;
    n_groups = materials.get_n_groups();
    n_moments = _n_moments;
    this->n_blocks = n_moments * n_groups;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal)
    {
      reinit_non_diagonal(materials, listen_to_material_id,
        boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_matrixfree)
    {
      reinit_full_matrixfree(materials, listen_to_material_id,
        boundary_conditions, albedo_factors);
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(this->n_blocks,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));
      this->locally_owned_dofs = dof_handler.locally_owned_dofs();

      DynamicSparsityPattern dsp(dof_handler.n_dofs());

      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm, this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int b1 = 0; b1 < this->n_blocks; b1++)
        for (unsigned int b2 = 0; b2 < this->n_blocks; b2++)
        {
          this->matrix_blocks[b1][b2] =
                                        new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[b1][b2]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs, this->sp, this->comm);
        }
      assemble_full_matrices(materials, boundary_conditions, albedo_factors);

    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixSDPN<dim, n_fe_degree>::reinit_non_diagonal (
    const Materials &materials,
    bool listen_to_material_id,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler,
      this->locally_relevant_dofs);

    // Build Transport mass coefficients
    coeffs_L.resize(n_groups,
      std::vector<std::vector<std::vector<Vector<double> > > >(n_groups,
        std::vector<std::vector<Vector<double> > >(n_moments,
          std::vector<Vector<double> >(n_moments,
            Vector<double>(materials.get_n_mats())))));

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        for (unsigned int mi = 0; mi < n_moments; mi++)
          for (unsigned int mj = 0; mj < n_moments; mj++)
          {
            for (unsigned int mat = 0; mat < materials.get_n_mats();
                mat++)
            {
              if (gi == gj)
              {
                coeffs_L[gi][gj][mi][mj][mat] = sdpn_coeff2[0][mi][mj]
                                                * materials.get_sigma_r(gi, mat);
                for (unsigned int m3 = 1; m3 < n_moments; m3++)
                {
                  coeffs_L[gi][gj][mi][mj][mat] +=
                                                   sdpn_coeff2[m3][mi][mj]
                                                   * materials.get_sigma_t(gi, mat);
                }
              }
              else
                coeffs_L[gi][gj][mi][mj][mat] = -sdpn_coeff2[0][mi][mj]
                                                * materials.get_sigma_s(gj, gi, mat);

            }
          }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_JxW_values);
    // FIME Esto Rompe zero flux SP3
    additional_data.mapping_update_flags_boundary_faces = (update_JxW_values);
    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
    unsigned int bi, bj;
    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int mi = 0; mi < n_moments; mi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mj = 0; mj < n_moments; mj++)
          {
            bi = gm_to_b(gi, mi);
            bj = gm_to_b(gj, mj);
            if (bi != bj)
            {
              this->mass_mf_blocks[bi][bj] = new MassOperator<dim,
                  n_fe_degree, double>(matfree_data);
              if (gi == gj)
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  boundary_conditions,
                  sdpn_marshack_coeff2[mi][mj],
                  listen_to_material_id); // With BC
              else
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  listen_to_material_id); // Without BC
            }
          }

    //  --------- Diagonal Blocks  ---------
    this->matrix_blocks.resize(this->n_blocks,
      std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));

    DynamicSparsityPattern dsp(this->locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

    SparsityTools::distribute_sparsity_pattern(dsp,
      this->locally_owned_dofs,
      this->comm,
      this->locally_relevant_dofs);
    this->sp.copy_from(dsp);

    for (unsigned int b = 0; b < this->n_blocks; b++)
    {
      this->matrix_blocks[b][b] = new PETScWrappers::MPI::SparseMatrix;
      this->matrix_blocks[b][b]->reinit(this->locally_owned_dofs,
        this->locally_owned_dofs, this->sp, this->comm);
    }

    assemble_diagonal_matrices(materials, boundary_conditions, albedo_factors);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixSDPN<dim, n_fe_degree>::reinit_full_matrixfree (
    const Materials &materials,
    bool listen_to_material_id,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {

    this->locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler,
      this->locally_relevant_dofs);

    // Build Transport mass coefficients
    coeffs_L.resize(n_groups,
      std::vector<std::vector<std::vector<Vector<double> > > >(n_groups,
        std::vector<std::vector<Vector<double> > >(n_moments,
          std::vector<Vector<double> >(n_moments,
            Vector<double>(materials.get_n_mats())))));

    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int gj = 0; gj < n_groups; gj++)
        for (unsigned int mi = 0; mi < n_moments; mi++)
          for (unsigned int mj = 0; mj < n_moments; mj++)
          {
            for (unsigned int mat = 0; mat < materials.get_n_mats();
                mat++)
            {
              if (gi == gj)
              {
                coeffs_L[gi][gj][mi][mj][mat] = sdpn_coeff2[0][mi][mj]
                                                * materials.get_sigma_r(gi, mat);
                for (unsigned int m3 = 1; m3 < n_moments; m3++)
                {
                  coeffs_L[gi][gj][mi][mj][mat] +=
                                                   sdpn_coeff2[m3][mi][mj]
                                                   * materials.get_sigma_t(gi, mat);
                }
              }
              else
                coeffs_L[gi][gj][mi][mj][mat] = -sdpn_coeff2[0][mi][mj]
                                                * materials.get_sigma_s(gj, gi, mat);

            }
          }

    //  --------- Matrix-Free Blocks  ---------
    //  Initialize Matrix free data
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_gradients
                                            | update_JxW_values);

    bool boundary = false;
    for (unsigned int bc = 0; bc < boundary_conditions.size(); bc++)
      if (boundary_conditions[bc] > 1)
      {
        boundary = true;
        break;
      }

    if (boundary == true)
      additional_data.mapping_update_flags_boundary_faces = (update_values
                                                             | update_JxW_values
                                                             | update_quadrature_points);

    MappingQ1<dim> mapping;
    matfree_data.reinit(mapping, dof_handler, constraints, QGauss<1>(n_fe_degree + 1),
      additional_data);

    this->mass_mf_blocks.resize(this->n_blocks,
      std::vector<MassOperator<dim, n_fe_degree, double>*>(
        this->n_blocks));
    unsigned int bi, bj;
    for (unsigned int gi = 0; gi < n_groups; gi++)
      for (unsigned int mi = 0; mi < n_moments; mi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mj = 0; mj < n_moments; mj++)
          {
            bi = gm_to_b(gi, mi);
            bj = gm_to_b(gj, mj);
            if (bi != bj)
            {
              this->mass_mf_blocks[bi][bj] = new MassOperator<dim,
                  n_fe_degree, double>(matfree_data);
              if (gi == gj)
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  boundary_conditions,
                  sdpn_marshack_coeff2[mi][mj],
                  listen_to_material_id); // With BC
              else
                this->mass_mf_blocks[bi][bj]->reinit(constraints,
                  materials.get_materials_vector(),
                  coeffs_L[gi][gj][mi][mj],
                  listen_to_material_id); // Without BC
            }
            //            std::cout << "L(" << bi << ", " << bj << ") = " << std::endl;
          }

    coeffs_grad.resize(this->n_blocks);
    coeffs_val.resize(this->n_blocks);
    coeffs_bound.resize(this->n_blocks);

    unsigned int b;
    for (unsigned int g = 0; g < n_groups; g++)
      for (unsigned int m = 0; m < n_moments; m++)
      {
        b = gm_to_b(g, m);
        coeffs_grad[b].reinit(materials.get_n_mats());
        coeffs_val[b].reinit(materials.get_n_mats());
        coeffs_bound[b] = sdpn_marshack_coeff2[m][m];
        for (unsigned int mat = 0; mat < materials.get_n_mats(); mat++)
        {
          coeffs_grad[b][mat] = sdpn_diff_coeff[m]
                                / (materials.get_sigma_t(g, mat));
          coeffs_val[b][mat] = coeffs_L[g][g][m][m][mat];
        }
      }

    this->poison_mf_blocks.resize(this->n_blocks);

    // Diagonal blocks
    for (unsigned int b = 0; b < this->n_blocks; b++)
    {
      this->poison_mf_blocks[b] =
                                  new PoissonOperator<dim, n_fe_degree, double>(
                                    matfree_data);

      this->poison_mf_blocks[b]->reinit(b, this->constraints, materials,
        materials.get_materials_vector(), coeffs_val[b],
        coeffs_grad[b], boundary_conditions, albedo_factors, true,
        coeffs_bound[b]);
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  unsigned int TransportMatrixSDPN<dim, n_fe_degree>::gm_to_b (
    const unsigned int group,
    const unsigned int moment) const
  {
    return n_moments * group + moment;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixSDPN<dim, n_fe_degree>::assemble_diagonal_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma, bc_coeff;
    double factor = 0.0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_bound(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    bool show_materials = false;
    bool show_boundary = false;
    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        const unsigned int mat = materials.get_material_id<dim>(cell);

        cell_grad = 0;
        cell_val = 0;
        cell_bound = 0;

        //
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
            {
              val = fe_values.shape_value(i, q_p)
                    * fe_values.shape_value(j, q_p)
                    * fe_values.JxW(q_p);

              cell_val(i, j) += val;

              grad = fe_values.shape_grad(i, q_p)
                     * fe_values.shape_grad(j, q_p)
                     * fe_values.JxW(q_p);

              cell_grad(i, j) += grad;
            }

        // Take care of albedo Boundary Conditions: Boundary integral
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
          {
            types::boundary_id boundary_id =
                                             cell->face(f)->boundary_id();
            AssertIndexRange(boundary_id, boundary_conditions.size());

            if (boundary_conditions[boundary_id] > 1)
            {
              fe_face_values.reinit(cell, f);
              switch (boundary_conditions[boundary_id])
              {
                case 2: // Vacuum BC
                  factor = 1.0;
                  break;
                default: // Custom Albedo BC
                  factor =
                           albedo_factors[boundary_conditions[boundary_id]
                                          - 3];
                  break;
              }

              for (unsigned int q_p = 0; q_p < n_face_q_points; ++q_p)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < dofs_per_cell;
                      ++j)
                  {
                    val = fe_face_values.shape_value(i, q_p)
                          * fe_face_values.shape_value(j, q_p)
                          * fe_face_values.JxW(q_p);

                    cell_bound(i, j) += factor * val;
                  }
            }
          }

        // Distribute in the Sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
          for (unsigned int mi = 0; mi < n_moments; ++mi)
          {
            // Get the material coefficients:
            D = sdpn_diff_coeff[mi] / (materials.get_sigma_tr(gi, mat));
            sigma = coeffs_L[gi][gi][mi][mi][mat];
            bc_coeff = sdpn_marshack_coeff2[mi][mi];
            cell_L.equ(D, cell_grad, sigma, cell_val, bc_coeff,
              cell_bound);
            constraints.distribute_local_to_global(cell_L,
              local_dof_indices,
              *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(gi,
                mi)]));
          }

      }

    // Compress the matrices
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      this->matrix_blocks[bi][bi]->compress(VectorOperation::add);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixSDPN<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors)
  {
    double val, grad;
    double D, sigma, bc_coeff;
    double factor = 0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_bound(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    bool show_materials = false;
    bool show_boundary = false;
    get_bool_from_options("-show_boundary", show_boundary);
    get_bool_from_options("-show_materials", show_materials);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      const unsigned int mat = materials.get_material_id<dim>(cell);

      cell_grad = 0;
      cell_val = 0;
      cell_bound = 0;

      //
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
          {
            val = fe_values.shape_value(i, q_p)
                  * fe_values.shape_value(j, q_p)
                  * fe_values.JxW(q_p);

            cell_val(i, j) += val;

            grad = fe_values.shape_grad(i, q_p)
                   * fe_values.shape_grad(j, q_p)
                   * fe_values.JxW(q_p);

            cell_grad(i, j) += grad;
          }

      // Take care of albedo Boundary Conditions: Boundary integral
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
        {
          types::boundary_id boundary_id = cell->face(f)->boundary_id();
          AssertIndexRange(boundary_id, boundary_conditions.size());

          if (boundary_conditions[boundary_id] > 1)
          {
            fe_face_values.reinit(cell, f);
            switch (boundary_conditions[boundary_id])
            {
              case 2: // Vacuum BC
                factor = 1.0;
                break;
              default: // Custom Albedo BC
                factor = albedo_factors[boundary_conditions[boundary_id]
                                        - 3];
                break;
            }

            for (unsigned int q_p = 0; q_p < n_face_q_points; ++q_p)
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  val = fe_face_values.shape_value(i, q_p)
                        * fe_face_values.shape_value(j, q_p)
                        * fe_face_values.JxW(q_p);

                  cell_bound(i, j) += factor * val;
                }
          }
        }

      // Distribute in the Sparse matrix
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          for (unsigned int mi = 0; mi < n_moments; ++mi)
            for (unsigned int mj = 0; mj < n_moments; ++mj)
            {
              // Get the material coefficients:
              if (gi == gj)
              {
                if (mi == mj)
                {
                  D = sdpn_diff_coeff[mi] / (materials.get_sigma_tr(gi, mat));
                  sigma = sdpn_coeff2[0][mi][mj]
                          * materials.get_sigma_r(gi, mat);
                  for (unsigned int m3 = 1; m3 < n_moments; ++m3)
                    sigma += sdpn_coeff2[m3][mi][mj]
                             * materials.get_sigma_t(gi, mat);
                  bc_coeff = sdpn_marshack_coeff2[mi][mj];
                  cell_L.equ(D, cell_grad, sigma, cell_val,
                    bc_coeff, cell_bound);
                  constraints.distribute_local_to_global(cell_L,
                    local_dof_indices,
                    *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)]));
                }
                else
                {
                  sigma = sdpn_coeff2[0][mi][mj]
                          * materials.get_sigma_r(gi, mat);
                  for (unsigned int m3 = 1; m3 < n_moments; ++m3)
                    sigma += sdpn_coeff2[m3][mi][mj]
                             * materials.get_sigma_t(gi, mat);
                  bc_coeff = sdpn_marshack_coeff2[mi][mj];
                  cell_L.equ(sigma, cell_val, bc_coeff,
                    cell_bound);
                  constraints.distribute_local_to_global(cell_L,
                    local_dof_indices,
                    *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)]));
                }
              }
              else
              {
                {
                  sigma = -sdpn_coeff2[0][mi][mj]
                          * materials.get_sigma_s(gj, gi, mat);
                  cell_L.equ(sigma, cell_val);
                  constraints.distribute_local_to_global(cell_L,
                    local_dof_indices,
                    *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)]));
                }
              }
            }
    }

    //
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks; ++bj)
        this->matrix_blocks[bi][bj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t TransportMatrixSDPN<dim, n_fe_degree>::memory_consumption () const
  {
    return 0;
    /*
     std::size_t memory = 0;
     if (this->matrixfree_type == non_diagonal)
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     {
     if (i == j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     else
     memory += this->mass_mf_blocks[i][j]->memory_consumption();
     }
     else
     {
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     }
     return memory;
     */
  }

//-------------------------------------------------------------------------------------------//
//-----------------------------------    FissionMatrix    ---------------------------------- //
//-------------------------------------------------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  FisionMatrixSDPN<dim, n_fe_degree>::FisionMatrixSDPN (const MPI_Comm &comm,
    const DoFHandler<dim> &_dof_handler,
    const AffineConstraints<double> &_constraints) :
      FisionMatrixBase<dim, n_fe_degree>(comm, _dof_handler, _constraints),
      dof_handler(
        _dof_handler),
      constraints(_constraints)
  {
    n_moments = 0;
  }

/**
 * @brief Star the Operator
 */
template <int dim, int n_fe_degree>
  void FisionMatrixSDPN<dim, n_fe_degree>::reinit (const Materials &materials,
    const unsigned int _n_moments,
    const MatrixFreeType &_matrixfree_type,
    bool listen_to_material_id)
  {
    const unsigned int n_groups = materials.get_n_groups();
    n_moments = _n_moments;
    this->n_blocks = _n_moments * n_groups;
    this->matrixfree_type = _matrixfree_type;
    this->n_dofs_block = dof_handler.n_dofs();

    if (this->matrixfree_type == non_diagonal
        or this->matrixfree_type == full_matrixfree)
    {
      const unsigned int n_mats = materials.get_n_mats();
      this->mass_mf_blocks.resize(this->n_blocks,
        std::vector<MassOperator<dim, n_fe_degree, double>*>(
          this->n_blocks));
      coeffs_M.resize(n_groups,
        std::vector<std::vector<std::vector<Vector<double> > > >(
          n_groups,
          std::vector<std::vector<Vector<double> > >(n_moments,
            std::vector<Vector<double> >(n_moments,
              Vector<double>(n_mats)))));
      //
      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
          for (unsigned int mi = 0; mi < n_moments; mi++)
            for (unsigned int mj = 0; mj < n_moments; mj++)
              for (unsigned int mat = 0; mat < n_mats; mat++)
              {
                coeffs_M[to_g][from_g][mi][mj][mat] =
                                                      sdpn_coeff2[0][mi][mj]
                                                      * materials.get_xi_nu_sigma_f(
                                                        from_g, to_g, mat);
              }

      //  --------- Matrix-Free Blocks  ---------
      //  Initialize Matrix free data
      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
          dealii::MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values);
      MappingQ1<dim> mapping;
      matfree_data.reinit(mapping, dof_handler, constraints,
        QGauss<1>(n_fe_degree + 1), additional_data);

      for (unsigned int gi = 0; gi < n_groups; gi++)
        for (unsigned int gj = 0; gj < n_groups; gj++)
          for (unsigned int mi = 0; mi < n_moments; mi++)
            for (unsigned int mj = 0; mj < n_moments; mj++)
            {
              this->mass_mf_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)] =
                                                                       new MassOperator<
                                                                           dim,
                                                                           n_fe_degree,
                                                                           double>(
                                                                         matfree_data);
              this->mass_mf_blocks[gm_to_b(gi, mi)][gm_to_b(gj, mj)]->reinit(
                constraints, materials.get_materials_vector(),
                coeffs_M[gi][gj][mi][mj],
                listen_to_material_id);
            }
    }
    else if (this->matrixfree_type == full_allocated)
    {
      // Resize matrix_blocks
      this->matrix_blocks.resize(this->n_blocks,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(this->n_blocks));

      this->locally_owned_dofs = dof_handler.locally_owned_dofs();
      DynamicSparsityPattern dsp(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);

      SparsityTools::distribute_sparsity_pattern(dsp,
        this->locally_owned_dofs,
        this->comm, this->locally_relevant_dofs);
      this->sp.copy_from(dsp);

      for (unsigned int bi = 0; bi < this->n_blocks; bi++)
        for (unsigned int bj = 0; bj < this->n_blocks; bj++)
        {
          this->matrix_blocks[bi][bj] =
                                        new PETScWrappers::MPI::SparseMatrix;
          this->matrix_blocks[bi][bj]->reinit(this->locally_owned_dofs,
            this->locally_owned_dofs, this->sp, this->comm);
        }

      assemble_full_matrices(materials);
    }
    else
      AssertRelease(false,
        "Invalid matrixfree_type: " + this->matrixfree_type);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void FisionMatrixSDPN<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials)
  {
    double val;
    double xi_nu_sigma_f;
    unsigned int mat;
    std::vector<types::global_dof_index> local_dof_indices(
      dof_handler.get_fe().dofs_per_cell);
    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_M(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      mat = materials.get_material_id<dim>(cell);

      cell_val = 0;
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          for (unsigned int q_p = 0; q_p < n_q_points; ++q_p)
          {
            val = fe_values.shape_value(i, q_p)
                  * fe_values.shape_value(j, q_p)
                  * fe_values.JxW(q_p);

            cell_val(i, j) += val;
          }

      for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          for (unsigned int mi = 0; mi < n_moments; ++mi)
            for (unsigned int mj = 0; mj < n_moments; ++mj)
            {
              // Get the material coefficients:
              xi_nu_sigma_f = sdpn_coeff2[0][mi][mj]
                              * materials.get_xi_nu_sigma_f(gj, gi, mat);
              cell_M.equ(xi_nu_sigma_f, cell_val);
              constraints.distribute_local_to_global(cell_M,
                local_dof_indices,
                *(this->matrix_blocks[gm_to_b(gi, mi)][gm_to_b(
                  gj, mj)]));
            }
    }

    // Compress
    for (unsigned int bi = 0; bi < this->n_blocks; ++bi)
      for (unsigned int bj = 0; bj < this->n_blocks; ++bj)
        this->matrix_blocks[bi][bj]->compress(VectorOperation::add);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  unsigned int FisionMatrixSDPN<dim, n_fe_degree>::gm_to_b (
    const unsigned int group,
    const unsigned int moment) const
  {
    return n_moments * group + moment;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  std::size_t FisionMatrixSDPN<dim, n_fe_degree>::memory_consumption () const
  {
    return 0;
    /*
     std::size_t memory = 0;
     if (this->matrixfree_type == non_diagonal)
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     {
     if (i == j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     else
     memory += this->mass_mf_blocks[i][j]->memory_consumption();
     }
     else
     {
     for (unsigned int i = 0; i < this->n_blocks; ++i)
     for (unsigned int j = 0; j < this->n_blocks; ++j)
     memory += this->matrix_blocks[i][j]->memory_consumption();
     }
     return memory;
     */
  }

// ----------- Explicit Instantations ----------- //

template class TransportMatrixSDPN<1, 1> ;
template class TransportMatrixSDPN<1, 2> ;
template class TransportMatrixSDPN<1, 3> ;
template class TransportMatrixSDPN<1, 4> ;
template class TransportMatrixSDPN<1, 5> ;

template class TransportMatrixSDPN<2, 1> ;
template class TransportMatrixSDPN<2, 2> ;
template class TransportMatrixSDPN<2, 3> ;
template class TransportMatrixSDPN<2, 4> ;
template class TransportMatrixSDPN<2, 5> ;

template class TransportMatrixSDPN<3, 1> ;
template class TransportMatrixSDPN<3, 2> ;
template class TransportMatrixSDPN<3, 3> ;
template class TransportMatrixSDPN<3, 4> ;
template class TransportMatrixSDPN<3, 5> ;

template class FisionMatrixSDPN<1, 1> ;
template class FisionMatrixSDPN<1, 2> ;
template class FisionMatrixSDPN<1, 3> ;
template class FisionMatrixSDPN<1, 4> ;
template class FisionMatrixSDPN<1, 5> ;

template class FisionMatrixSDPN<2, 1> ;
template class FisionMatrixSDPN<2, 2> ;
template class FisionMatrixSDPN<2, 3> ;
template class FisionMatrixSDPN<2, 4> ;
template class FisionMatrixSDPN<2, 5> ;

template class FisionMatrixSDPN<3, 1> ;
template class FisionMatrixSDPN<3, 2> ;
template class FisionMatrixSDPN<3, 3> ;
template class FisionMatrixSDPN<3, 4> ;
template class FisionMatrixSDPN<3, 5> ;

