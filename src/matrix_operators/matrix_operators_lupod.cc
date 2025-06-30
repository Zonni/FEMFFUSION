/**
 * @file   matrix_operators_lupod.cc
 * @brief
 */

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

#include <deal.II/matrix_free/matrix_free.h>

#include <algorithm>
#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "../../include/matrix_operators/matrix_operators_lupod.h"

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  TransportMatrixReduced<dim, n_fe_degree>::TransportMatrixReduced (
    const MPI_Comm &_comm,
    const DoFHandler<dim> &_dof_handler,
    const AffineConstraints<double> &_constraints) :
      comm(_comm),
      dof_handler(_dof_handler),
      constraints(_constraints)
  {
    n_blocks = 0;
    matrixfree_type = full_allocated;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixReduced<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const std::vector<std::vector<unsigned int> > &points_per_block,
    const MatrixFreeType &_matrixfree_type)
  {

    const unsigned int n_groups = materials.get_n_groups();
    n_blocks = n_groups;
    matrixfree_type = _matrixfree_type;
    AssertRelease(matrixfree_type == full_allocated,
      "matrixfree_type must be full_allocated");

//    if (matrixfree_type == non_diagonal)
//    {
//      //reinit_non_diagonal(materials, boundary_conditions, albedo_factors);
//    }
//    else if (matrixfree_type == full_matrixfree)
//    {
//      //reinit_full_matrixfree(materials, boundary_conditions, albedo_factors);
//    }
    if (matrixfree_type == full_allocated)
    {
      unsigned int n_LUPOD_points = points_per_block[0].size();
      // Resize matrix_sp_blocks
      std::vector<DynamicSparsityPattern> dsp(n_groups);
      matrix_sp_blocks.resize(n_groups,
        std::vector<SparseMatrix<double>>(n_groups));
      for (unsigned int g = 0; g < n_groups; g++)
        dsp[g].reinit(n_LUPOD_points, this->dof_handler.n_dofs());

      // -----------------------------------------------------------------------------------------
      //DoFTools::make_sparsity_pattern(this->dof_handler, dsp, this->constraints, true);
      std::vector<types::global_dof_index> dofs_on_this_cell;
      unsigned int point_index = 0;
      for (const auto &cell : this->dof_handler.active_cell_iterators())
      {
        const unsigned int n_dofs_per_cell = cell->get_fe().n_dofs_per_cell();
        dofs_on_this_cell.resize(n_dofs_per_cell);
        cell->get_dof_indices(dofs_on_this_cell);
        for (unsigned int g = 0; g < n_groups; g++)
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
          {
            point_index = binarySearch(points_per_block[g], dofs_on_this_cell[i]);
            if (point_index != static_cast<unsigned int>(-1))
            {
              for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
              {
                dsp[g].add(point_index, dofs_on_this_cell[j]);
              }
            }
          }
      }

      // Copy dsp to SparsityPattern
      sp.resize(n_groups);
      for (unsigned int g = 0; g < n_groups; g++)
        sp[g].copy_from(dsp[g]);

      for (unsigned int g = 0; g < materials.get_n_groups(); g++)
        sp[g].compress();

      //--------------------------------------------------------------------------

      //assemble_full_matrices(materials, boundary_conditions, albedo_factors, points_per_block);
    }

    else
      AssertRelease(false, "Invalid matrixfree_type: " + matrixfree_type);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixReduced<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<unsigned int> &boundary_conditions,
    const std::vector<double> &albedo_factors,
    const std::vector<std::vector<unsigned int> > &points)
  {

    for (unsigned int g1 = 0; g1 < materials.get_n_groups(); g1++)
      for (unsigned int g2 = 0; g2 < materials.get_n_groups(); g2++)
      {
        matrix_sp_blocks[g1][g2].reinit(this->sp[g1]);
      }

    double val, grad;
    double D, sigma_r, sigma_s;
    double factor = 0;
    QGauss<dim> quadrature_formula(n_fe_degree + 1);
    QGauss<dim - 1> face_quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(this->dof_handler.get_fe(),
      quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);
    FEFaceValues<dim> fe_face_values(this->dof_handler.get_fe(),
      face_quadrature_formula,
      update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = this->dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    unsigned int point_index;

    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_grad(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double> > bound(materials.get_n_groups(),
      FullMatrix<double>(dofs_per_cell, dofs_per_cell));
    for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
      bound[g] = FullMatrix<double>(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<types::global_dof_index> local_row_indices;

    typename DoFHandler<dim>::active_cell_iterator cell =
        this->dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator endc = this->dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        const unsigned int mat = materials.get_material_id<dim>(cell);

        cell->get_dof_indices(local_dof_indices);

        for (unsigned int gi = 0; gi < materials.get_n_groups(); gi++)
        {
          cell_grad = 0;
          cell_val = 0;
          local_row_indices.clear();
          local_row_indices.reserve(dofs_per_cell);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            point_index = binarySearch(points[gi], local_dof_indices[i]);
            if (point_index != static_cast<unsigned int>(-1))
            {
              local_row_indices.push_back(point_index);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                for (unsigned int qp = 0; qp < n_q_points; ++qp)
                {
                  val = fe_values.shape_value(i, qp)
                        * fe_values.shape_value(j, qp)
                        * fe_values.JxW(qp);

                  cell_val(i, j) += val;

                  grad = fe_values.shape_grad(i, qp)
                         * fe_values.shape_grad(j, qp)
                         * fe_values.JxW(qp);

                  cell_grad(i, j) += grad;
                }
            }
            else
            {
              local_row_indices.push_back(0);
            }
          }

          //for (unsigned int g = 0; g < materials.get_n_groups(); ++g)
          bound[gi] = 0;

          // Take care of albedo Boundary Conditions: Boundary integral
          for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
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
                    factor = 0.5;
                    break;
                  default: // Custom Albedo BC
                    factor = albedo_factors[(boundary_conditions[boundary_id] - 3)
                                            * materials.get_n_groups()
                                            + gi];
                    break;
                }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  point_index = binarySearch(points[gi], local_dof_indices[i]);
                  if (point_index != static_cast<unsigned int>(-1))
                  {
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      for (unsigned int qp = 0; qp < n_face_q_points; ++qp)
                      {
                        val = fe_face_values.shape_value(i, qp)
                              * fe_face_values.shape_value(j, qp)
                              * fe_face_values.JxW(qp);
                        bound[gi](i, j) += factor * val;
                      }
                  }
                }
              }

            }
          }

          //std::cout << "Cell " << cell << std::endl;
          //std::cout << "at boundary? " <<  cell->at_boundary() << std::endl;

//          std::cout << " cell_grad " << std::endl;
//          cell_grad.print_formatted(std::cout);
//          std::cout << " cell_val " << std::endl;
//          cell_val.print_formatted(std::cout);
//          std::cout << " bound[gi] "  << std::endl;
//          bound[gi].print_formatted(std::cout);
//          std::cout << " local_row_indices " << std::endl;
//          print_vector(local_row_indices);
          // Distribute in the Sparse matrix
          for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          {
            // Get the material coefficients:
            if (gi == gj)
            {
              D = materials.get_diffusion_coefficient(gi, mat);
              sigma_r = materials.get_sigma_r(gi, mat);
              cell_L.equ(D, cell_grad, sigma_r, cell_val, 1.0, bound[gi]);
              constraints.distribute_local_to_global(cell_L,
                local_row_indices,
                local_dof_indices,
                matrix_sp_blocks[gi][gj]);
            }
            else
            {
              sigma_s = -materials.get_sigma_s(gj, gi, mat);
              cell_L.equ(sigma_s, cell_val);
              constraints.distribute_local_to_global(cell_L,
                local_row_indices,
                local_dof_indices,
                matrix_sp_blocks[gi][gj]);
            }
          }
        }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        this->matrix_sp_blocks[gi][gj].compress(VectorOperation::add);
  }

//-------------------------------------------------//
//  --------------- FissionMatrix    --------------//
//-------------------------------------------------//
/**
 * @brief Constructor of FissionMatrix. Just copy references to DoFHandler and AffineConstraints<double>
 */
template <int dim, int n_fe_degree>
  FisionMatrixReduced<dim, n_fe_degree>::FisionMatrixReduced (
    const MPI_Comm &_comm,
    const DoFHandler<dim> &_dof_handler,
    const AffineConstraints<double> &_constraints) :
      comm(_comm),
      dof_handler(_dof_handler),
      constraints(_constraints)
  {
    n_blocks = 0;
    matrixfree_type = full_matrixfree;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixReduced<dim, n_fe_degree>::reinit (
    const Materials &materials,
    const std::vector<std::vector<unsigned int> > &points_per_block,
    const MatrixFreeType &_matrixfree_type)
  {
    const unsigned int n_groups = materials.get_n_groups();
    n_blocks = n_groups;
    matrixfree_type = _matrixfree_type;

    AssertRelease(matrixfree_type == full_allocated,
      "Other matrixfree_types not implemented");

//    if (matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree)
//    {
//      const unsigned int n_mats = materials.get_n_mats();
//      this->mass_mf_blocks.resize(this->n_blocks,
//        std::vector<MassOperator<dim, n_fe_degree, double>*>(this->n_blocks));
//      this->coeffs.resize(n_groups,
//        std::vector<Vector<double> >(n_groups));
//// Fill coeffs
//      for (unsigned int from_g = 0; from_g < n_groups; from_g++)
//        for (unsigned int to_g = 0; to_g < n_groups; to_g++)
//        {
//          this->coeffs[to_g][from_g].reinit(n_mats);
//          for (unsigned int mat = 0; mat < n_mats; mat++)
//          {
//            this->coeffs[to_g][from_g][mat] = materials.get_xi_nu_sigma_f(from_g, to_g,
//              mat);
//          }
//        }
//
//      //  --------- Matrix-Free Blocks  ---------
//      //  Initialize Matrix free data
//      typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
//      additional_data.tasks_parallel_scheme =
//          dealii::MatrixFree<dim, double>::AdditionalData::none;
//      additional_data.mapping_update_flags = (update_values | update_JxW_values);
//      this->matfree_data.reinit(this->dof_handler,
//        this->constraints,
//        QGauss<1>(n_fe_degree + 1),
//        additional_data);
//
//      for (unsigned int gi = 0; gi < n_groups; gi++)
//        for (unsigned int gj = 0; gj < n_groups; gj++)
//        {
//          this->mass_mf_blocks[gi][gj] =
//                                         new MassOperator<dim, n_fe_degree, double>(
//                                           this->matfree_data);
//          this->mass_mf_blocks[gi][gj]->reinit(
//            this->constraints,
//            materials.get_materials_vector(),
//            this->coeffs[gi][gj],
//            materials.listen_to_material_id);
//        }
//    }
//    else if (this->matrixfree_type == full_allocated)
//    {
    // Resize matrix_sp_blocks

    matrix_sp_blocks.resize(n_groups,
      std::vector<SparseMatrix<double> >(n_groups));

    std::vector<DynamicSparsityPattern> dsp(n_groups);
    for (unsigned int g = 0; g < n_groups; g++)
      dsp[g].reinit(points_per_block[0].size(), this->dof_handler.n_dofs());

    // -----------------------------------------------------------------------------------------
    //DoFTools::make_sparsity_pattern(this->dof_handler, dsp, this->constraints, true);
    std::vector<types::global_dof_index> dofs_on_this_cell;
    dofs_on_this_cell.reserve(
      this->dof_handler.get_fe_collection().max_dofs_per_cell());
    for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
      dofs_on_this_cell.resize(dofs_per_cell);
      cell->get_dof_indices(dofs_on_this_cell);
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          unsigned int point_index = binarySearch(points_per_block[g],
            dofs_on_this_cell[i]);
          if (point_index != static_cast<unsigned int>(-1))
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              dsp[g].add(point_index, dofs_on_this_cell[j]);
            }
          }
        }
    }

    //--------------------------------------------------------------------------
    sp.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; g++)
      sp[g].copy_from(dsp[g]);

    for (unsigned int g1 = 0; g1 < n_groups; g1++)
    {
      this->sp[g1].compress();
      for (unsigned int g2 = 0; g2 < n_groups; g2++)
      {
        matrix_sp_blocks[g1][g2].reinit(this->sp[g1]);
      }
    }

    for (unsigned int g = 0; g < materials.get_n_groups(); g++)
      sp[g].compress();

    //this->assemble_full_matrices(materials, points_per_block);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void TransportMatrixReduced<dim, n_fe_degree>::vmult (
    BlockVector<double> &dst,
    const BlockVector<double> &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    // Set dst to 0
    dst = 0.0;

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixReduced<dim, n_fe_degree>::vmult_add (
    const unsigned int row,
    const unsigned int col,
    Vector<double> &dst,
    const Vector<double> &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    matrix_sp_blocks[row][col].vmult_add(dst, src);

//    if ((matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree) and row
//        != col)
//      mass_mf_blocks[row][col]->vmult_add(dst, src);
//    else if (matrixfree_type == full_matrixfree)
//      poison_mf_blocks[row]->vmult_add(dst, src);
//    else
//      matrix_blocks[row][col]->vmult_add(dst, src);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void TransportMatrixReduced<dim, n_fe_degree>::clear ()
  {

//    if (matrixfree_type == non_diagonal)
//    {
//      for (unsigned int row = 0; row < n_blocks; ++row)
//        for (unsigned int col = 0; col < n_blocks; ++col)
//        {
//          if (row == col)
//          {
//            matrix_blocks[row][col]->clear();
//            delete matrix_blocks[row][col];
//
//          }
//          else
//          {
//            mass_mf_blocks[row][col]->clear();
//            delete mass_mf_blocks[row][col];
//          }
//        }
//    }
//    else if (matrixfree_type == full_matrixfree)
//    {
//      for (unsigned int row = 0; row < n_blocks; ++row)
//        for (unsigned int col = 0; col < n_blocks; ++col)
//          if (row != col)
//          {
//            mass_mf_blocks[row][col]->clear();
//            delete mass_mf_blocks[row][col];
//          }
//          else
//          {
//            poison_mf_blocks[row]->clear();
//            delete poison_mf_blocks[row];
//          }
//    }
//    else
//      for (unsigned int row = 0; row < n_blocks; ++row)
//        for (unsigned int col = 0; col < n_blocks; ++col)
//          matrix_blocks[row][col]->clear();

    for (unsigned int row = 0; row < n_blocks; ++row)
      for (unsigned int col = 0; col < n_blocks; ++col)
        matrix_sp_blocks[row][col].clear();

    n_blocks = 0;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixReduced<dim, n_fe_degree>::assemble_full_matrices (
    const Materials &materials,
    const std::vector<std::vector<unsigned int> > &points)
  {
    for (unsigned int g1 = 0; g1 < materials.get_n_groups(); g1++)
      for (unsigned int g2 = 0; g2 < materials.get_n_groups(); g2++)
      {
        matrix_sp_blocks[g1][g2].reinit(this->sp[g1]);
      }

    //-------------------------------------------------------------------
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
    std::vector<unsigned int> local_row_indices;
    //local_row_indices.reserve(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        mat = materials.get_material_id<dim>(cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
        {
          cell_val = 0;
          local_row_indices.clear();
          local_row_indices.reserve(dofs_per_cell);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            unsigned int point_index = binarySearch(points[gi], local_dof_indices[i]);
            if (point_index != static_cast<unsigned int>(-1))
            {
              local_row_indices.push_back(point_index);

              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                for (unsigned int qp = 0; qp < n_q_points; ++qp)
                {
                  val = fe_values.shape_value(i, qp)
                        * fe_values.shape_value(j, qp)
                        * fe_values.JxW(qp);

                  cell_val(i, j) += val;
                }
              }
            }
            else
            {
              local_row_indices.push_back(0);
            }
          }

          for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
          {
            // Get the material coefficients:
            xi_nu_sigma_f = materials.get_xi_nu_sigma_f(gj, gi, mat);
            cell_M.equ(xi_nu_sigma_f, cell_val);

            constraints.distribute_local_to_global(cell_M,
              local_row_indices,
              local_dof_indices,
              matrix_sp_blocks[gi][gj]);
          }
        }
      }

    for (unsigned int gi = 0; gi < materials.get_n_groups(); ++gi)
      for (unsigned int gj = 0; gj < materials.get_n_groups(); ++gj)
        matrix_sp_blocks[gi][gj].compress(VectorOperation::add);
  }

/**
 * @brief Complete matrix-vector multiplication.
 * dst = TransportMatrixBase * src
 */
template <int dim, int n_fe_degree>
  void FisionMatrixReduced<dim, n_fe_degree>::vmult (
    BlockVector<double> &dst,
    const BlockVector<double> &src) const
  {
    AssertDimension(dst.n_blocks(), n_blocks);

    // Set dst to 0
    dst = 0;

    for (unsigned int i = 0; i < n_blocks; i++)
      for (unsigned int j = 0; j < n_blocks; j++)
      {
        vmult_add(i, j, dst.block(i), src.block(j));
      }
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void FisionMatrixReduced<dim, n_fe_degree>::vmult_add (
    const unsigned int row,
    const unsigned int col,
    Vector<double> &dst,
    const Vector<double> &src) const
  {
    AssertIndexRange(row, n_blocks);
    AssertIndexRange(col, n_blocks);

    matrix_sp_blocks[row][col].vmult_add(dst, src);

//    if ((matrixfree_type == non_diagonal or matrixfree_type == full_matrixfree) and row
//        != col)
//      mass_mf_blocks[row][col]->vmult_add(dst, src);
//    else if (matrixfree_type == full_matrixfree)
//      poison_mf_blocks[row]->vmult_add(dst, src);
//    else
//      matrix_blocks[row][col]->vmult_add(dst, src);
  }

template <int dim, int n_fe_degree>
  void FisionMatrixReduced<dim, n_fe_degree>::clear ()
  {
//    if (matrixfree_type == non_diagonal)
//    {
//      for (unsigned int row = 0; row < n_blocks; ++row)
//        for (unsigned int col = 0; col < n_blocks; ++col)
//        {
//          if (row == col)
//          {
//            matrix_blocks[row][col]->clear();
//            delete matrix_blocks[row][col];
//
//          }
//          else
//          {
//            mass_mf_blocks[row][col]->clear();
//            delete mass_mf_blocks[row][col];
//          }
//        }
//    }
//    else if (matrixfree_type == full_matrixfree)
//    {
//      for (unsigned int row = 0; row < n_blocks; ++row)
//        for (unsigned int col = 0; col < n_blocks; ++col)
//          if (row != col)
//          {
//            mass_mf_blocks[row][col]->clear();
//            delete mass_mf_blocks[row][col];
//          }
//          else
//          {
//            poison_mf_blocks[row]->clear();
//            delete poison_mf_blocks[row];
//          }
//    }
//    else
//      for (unsigned int row = 0; row < n_blocks; ++row)
//        for (unsigned int col = 0; col < n_blocks; ++col)
//          matrix_blocks[row][col]->clear();

    for (unsigned int row = 0; row < n_blocks; ++row)
      for (unsigned int col = 0; col < n_blocks; ++col)
        matrix_sp_blocks[row][col].clear();

    n_blocks = 0;
  }

// ----------- Explicit Instantations ----------- //

template class TransportMatrixReduced<1, 1> ;
template class TransportMatrixReduced<1, 2> ;
template class TransportMatrixReduced<1, 3> ;
template class TransportMatrixReduced<1, 4> ;
template class TransportMatrixReduced<1, 5> ;

template class TransportMatrixReduced<2, 1> ;
template class TransportMatrixReduced<2, 2> ;
template class TransportMatrixReduced<2, 3> ;
template class TransportMatrixReduced<2, 4> ;
template class TransportMatrixReduced<2, 5> ;

template class TransportMatrixReduced<3, 1> ;
template class TransportMatrixReduced<3, 2> ;
template class TransportMatrixReduced<3, 3> ;
template class TransportMatrixReduced<3, 4> ;
template class TransportMatrixReduced<3, 5> ;

template class FisionMatrixReduced<1, 1> ;
template class FisionMatrixReduced<1, 2> ;
template class FisionMatrixReduced<1, 3> ;
template class FisionMatrixReduced<1, 4> ;
template class FisionMatrixReduced<1, 5> ;

template class FisionMatrixReduced<2, 1> ;
template class FisionMatrixReduced<2, 2> ;
template class FisionMatrixReduced<2, 3> ;
template class FisionMatrixReduced<2, 4> ;
template class FisionMatrixReduced<2, 5> ;

template class FisionMatrixReduced<3, 1> ;
template class FisionMatrixReduced<3, 2> ;
template class FisionMatrixReduced<3, 3> ;
template class FisionMatrixReduced<3, 4> ;
template class FisionMatrixReduced<3, 5> ;
