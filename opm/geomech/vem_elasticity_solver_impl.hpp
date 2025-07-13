//==============================================================================
//!
//! \file elasticity_upscale_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class - template implementations
//!
//==============================================================================
#ifndef OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP
#define OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/vem/vem.hpp>
#include <opm/geomech/vem/vemutils.hpp>

#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <tuple>
#include <vector>

namespace Opm::Elasticity
{
template <typename GridType>
void
VemElasticitySolver<GridType>::assignToVoigt(
    Dune::BlockVector<Dune::FieldVector<double, 6>>& voigt_stress,
    const Dune::BlockVector<Dune::FieldVector<double, 1>>& vemstress)
{
    // This is used in flow!!!
    for (std::size_t i = 0; i < voigt_stress.size(); ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            voigt_stress[i][k] = vemstress[i * 6 + k];
        }

        // transforming from xx,yy,zz,xy,xz,yz ordring to voigt i.e. xx,yy,zz,zy,zx,xy
        voigt_stress[i][3] = vemstress[i * 6 + 4];
        voigt_stress[i][5] = vemstress[i * 6 + 3];
        voigt_stress[i][4] = vemstress[i * 6 + 5];
    }
}

template <typename GridType>
void
VemElasticitySolver<GridType>::assignToVoigtSymMat(
    Dune::BlockVector<Dune::FieldVector<double, 6>>& voigt_stress,
    const std::vector<std::array<double, 6>>& vemstress)
{
    for (std::size_t i = 0; i < voigt_stress.size(); ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            voigt_stress[i][k] = vemstress[i][k];
        }

        // transforming from xx,yy,zz,xy,xz,yz ordring to voigt i.e. xx,yy,zz,zy,zx,xy
        voigt_stress[i][3] = vemstress[i][4]; // yz??
        voigt_stress[i][5] = vemstress[i][3]; // xy
        voigt_stress[i][4] = vemstress[i][5]; // xz??
    }
}

template <typename GridType>
void
VemElasticitySolver<GridType>::calculateStressPrecomputed(const Vector& dispalldune)
{
    OPM_TIMEBLOCK(calculateStressPrecomputed);

    // NB stressmat is defined in linear indices not block linear indices
    // Vector dispalldune;
    // dispalldune.resize(3 * grid_.leafGridView().size(3));
    // this->expandSolution(dispalldune, this->u);
    //  Dune::BlockVector< DuneFieldVector<double,1> >
    Vector stress(6 * grid_.leafGridView().size(0));

    stressmat_.mv(dispalldune, stress);
    stress_.resize(num_cells_);

    assignToVoigt(stress_, stress);
}

template <typename GridType>
void
VemElasticitySolver<GridType>::expandDisp(std::vector<double>& dispall, const bool expand)
{
    dispall.resize(3 * grid_.leafGridView().size(3));

    if (expand) {
        Vector dispalldune;
        dispalldune.resize(3 * grid_.leafGridView().size(3));

        this->expandSolution(dispalldune, this->u);

        for (std::size_t i = 0; i < dispall.size(); ++i) {
            dispall[i] = dispalldune[i]; // fieldvector<double,1> can be converted to double
        }
    } else {
        dispall = this->u;
    }
}

template <typename GridType>
void
VemElasticitySolver<GridType>::calculateStress()
{
    bool expand = true;

    OPM_TIMEBLOCK(calculateStressFull);

    // assumes the grid structure is made
    num_cells_ = grid_.leafGridView().size(0); // entities of codim 0

    // assemble the mechanical system
    // const int numdof =
    //  const int tot_num_faces = accumulate(num_cell_faces_, num_cell_faces_ +
    //  num_cells_, 0); const int tot_num_fcorners = accumulate(num_face_corners_,
    //  &num_face_corners_[0] + tot_num_faces, 0); const int tot_num_nodes =
    //  *max_element(face_corners_, face_corners + tot_num_fcorners) + 1;
    stress_.resize(num_cells_);

    std::vector<std::array<double, 6>> stress;
    stress.resize(num_cells_);

    std::vector<double> dispall;
    this->expandDisp(dispall, expand);

    std::vector<std::tuple<int, int, double>> stressmat;
    vem::compute_stress_3D(&coords_[0],
                           num_cells_,
                           &num_cell_faces_[0],
                           &num_face_corners_[0],
                           &face_corners_[0],
                           &ymodule_[0],
                           &pratio_[0],
                           dispall,
                           stress,
                           stressmat,
                           false,
                           true);

    // copy to dune definitions
    stress_.resize(num_cells_);
    assignToVoigtSymMat(stress_, stress);
}

template <typename GridType>
void
VemElasticitySolver<GridType>::calculateStrainPrecomputed(const Vector& dispalldune)
{
    OPM_TIMEBLOCK(calculateStrainPrecomputed);

    // NB stressmat is defined in linear indices not block linear indices
    // Vector dispalldune;
    // dispalldune.resize(3 * grid_.leafGridView().size(3));
    // this->expandSolution(dispalldune, this->u);
    //  Dune::BlockVector< DuneFieldVector<double,1> >
    Vector strain(6 * grid_.leafGridView().size(0));

    strainmat_.mv(dispalldune, strain);
    strain_.resize(num_cells_);

    assignToVoigt(strain_, strain);
}

template <typename GridType>
void
VemElasticitySolver<GridType>::calculateStrain()
{
    bool expand = true; // assumes that matrices is caculated with reduced_boundary = true
    OPM_TIMEBLOCK(calculateStressFull);

    // assumes the grid structure is made
    num_cells_ = grid_.leafGridView().size(0); // entities of codim 0

    // assemble the mechanical system
    // const int numdof =
    //  const int tot_num_faces = accumulate(num_cell_faces_, num_cell_faces_ +
    //  num_cells_, 0); const int tot_num_fcorners = accumulate(num_face_corners_,
    //  &num_face_corners_[0] + tot_num_faces, 0); const int tot_num_nodes =
    //  *max_element(face_corners_, face_corners + tot_num_fcorners) + 1;
    strain_.resize(num_cells_);
    std::vector<std::array<double, 6>> strain;
    strain.resize(num_cells_);

    std::vector<double> dispall;
    this->expandDisp(dispall, expand);

    std::vector<std::tuple<int, int, double>> stressmat;
    vem::compute_stress_3D(&coords_[0],
                           num_cells_,
                           &num_cell_faces_[0],
                           &num_face_corners_[0],
                           &face_corners_[0],
                           &ymodule_[0],
                           &pratio_[0],
                           dispall,
                           strain,
                           stressmat,
                           false,
                           false);

    // copy to dune definitions
    strain_.resize(num_cells_);

    assignToVoigtSymMat(strain_, strain);
}

template <typename GridType>
void
VemElasticitySolver<GridType>::assemble(const Vector& pressure,
                                        const bool do_matrix,
                                        const bool do_vector,
                                        const bool reduce_boundary)
{
    OPM_TIMEBLOCK(assemble);

    A.getLoadVector() = 0.0;

    if (do_matrix) {
        A.getOperator() = 0.0;

        vem::getGridVectors(grid_, coords_, num_cell_faces_, num_face_corners_, face_corners_);

        const auto num_fixed_dofs = std::get<0>(dirichlet_);
        const auto& fixed_dof_ixs = std::get<1>(dirichlet_);
        const auto& fixed_dof_values = std::get<2>(dirichlet_);

        // neumann boundary conditions
        const int num_neumann_faces = 0;
        num_cells_ = grid_.leafGridView().size(0); // entities of codim 0

        // assemble the mechanical system
        std::vector<std::tuple<int, int, double>> A_entries;

        {
            OPM_TIMEBLOCK(assembleVEM);

            const auto stability_choice = vem::StabilityChoice::D_RECIPE;

            vem::assemble_mech_system_3D(&coords_[0],
                                         num_cells_,
                                         &num_cell_faces_[0],
                                         &num_face_corners_[0],
                                         &face_corners_[0],
                                         &ymodule_[0],
                                         &pratio_[0],
                                         &body_force_[0],
                                         num_fixed_dofs,
                                         &fixed_dof_ixs[0],
                                         &fixed_dof_values[0],
                                         num_neumann_faces,
                                         nullptr,
                                         nullptr,
                                         A_entries,
                                         rhs_force_,
                                         stability_choice,
                                         reduce_boundary);
        }

        this->makeDuneSystemMatrix(A_entries);

        {
            OPM_TIMEBLOCK(setUpExtraStructuresForVEM);

            // make indexing for div operator i.e. all nodes to dofs
            std::vector<int> dof_idx(grid_.leafGridView().size(3) * 3);
            std::iota(dof_idx.begin(), dof_idx.end(), 0);

            std::set_difference(dof_idx.begin(),
                                dof_idx.end(),
                                fixed_dof_ixs.begin(),
                                fixed_dof_ixs.end(),
                                std::back_inserter(idx_free_));

            //
            std::vector<double> rhs_tmp(pressure.size(), 0);
            std::vector<double> pressure_tmp(pressure.size(), 0);
            std::vector<std::tuple<int, int, double>> divmat;
            vem::potential_gradient_force_3D(&coords_[0],
                                             num_cells_,
                                             &num_cell_faces_[0],
                                             &num_face_corners_[0],
                                             &face_corners_[0],
                                             &pressure_tmp[0],
                                             rhs_tmp,
                                             divmat,
                                             true);

            // sort(divmat_.begin(),
            //      divmat_.end(),
            //      [](const auto& aa, const auto& bb) { return std::get<1>(aa) <
            //      std::get<1>(bb); });
            std::vector<int> global_to_dof(grid_.leafGridView().size(3) * 3, -1);
            for (std::size_t i = 0; i < idx_free_.size(); ++i) {
                if (reduce_boundary) {
                    // new numbering of dofs
                    global_to_dof[idx_free_[i]] = i;
                } else {
                    global_to_dof[idx_free_[i]] = idx_free_[i];
                }
            }

            // renumber and eliminate dof fized
            std::vector<std::tuple<int, int, double>> divmatdof;
            for (const auto& elem : divmat) {
                const int I = global_to_dof[std::get<0>(elem)];
                if (I > -1) {
                    divmatdof.emplace_back(I, std::get<1>(elem), std::get<2>(elem));
                }
            }

            // finaly make dune matrix
            divmat_.setBuildMode(Matrix::implicit);

            // map from dof=3*nodes at a cell (ca 3*3*3) to cell
            divmat_.setImplicitBuildModeParameters(3 * 3 * 3, 0.4);

            if (reduce_boundary) {
                divmat_.setSize(idx_free_.size(), num_cells_);
            } else {
                divmat_.setSize(grid_.leafGridView().size(3) * 3, num_cells_);
            }

            makeDuneMatrixCompressed(divmatdof, divmat_);

            // also make stress matrix
            std::vector<double> dispall(grid_.leafGridView().size(3) * 3);
            std::vector<std::array<double, 6>> stresstmp(grid_.leafGridView().size(0));
            std::vector<std::tuple<int, int, double>> stressmat;
            {
                OPM_TIMEBLOCK(setUpStressStressMatrix);
                vem::compute_stress_3D(&coords_[0],
                                       num_cells_,
                                       &num_cell_faces_[0],
                                       &num_face_corners_[0],
                                       &face_corners_[0],
                                       &ymodule_[0],
                                       &pratio_[0],
                                       dispall,
                                       stresstmp,
                                       stressmat,
                                       true,
                                       true);
            }

            stressmat_.setBuildMode(Matrix::implicit);
            stressmat_.setImplicitBuildModeParameters(3 * 3 * 3, 0.4);
            stressmat_.setSize(num_cells_ * 6, dispall.size());

            makeDuneMatrixCompressed(stressmat, stressmat_);

            std::vector<std::tuple<int, int, double>> strainmat;
            {
                OPM_TIMEBLOCK(setUpStressStrainMatrix);
                vem::compute_stress_3D(&coords_[0],
                                       num_cells_,
                                       &num_cell_faces_[0],
                                       &num_face_corners_[0],
                                       &face_corners_[0],
                                       &ymodule_[0],
                                       &pratio_[0],
                                       dispall,
                                       stresstmp,
                                       strainmat,
                                       true,
                                       false);
            }

            strainmat_.setBuildMode(Matrix::implicit);
            strainmat_.setImplicitBuildModeParameters(3 * 3 * 3, 0.4);
            strainmat_.setSize(num_cells_ * 6, dispall.size());

            makeDuneMatrixCompressed(strainmat, strainmat_);
        }
    }

    if (do_vector) {
        OPM_TIMEBLOCK(calculateRHS);
        // NB rhs_force_ is calculated by matrix call
        std::vector<double> rhs_pressure;
        std::vector<double> std_pressure(pressure.size(), 0.0);

        for (std::size_t i = 0; i < pressure.size(); ++i) {
            std_pressure[i] = pressure[i][0];
        }

        std::vector<std::tuple<int, int, double>> divmat;
        vem::potential_gradient_force_3D(&coords_[0],
                                         num_cells_,
                                         &num_cell_faces_[0],
                                         &num_face_corners_[0],
                                         &face_corners_[0],
                                         &std_pressure[0],
                                         rhs_pressure,
                                         divmat,
                                         false);

        // Sign is added here  i.e \div \sigma =
        std::vector<double> rhs(rhs_force_);
        if (reduce_boundary) {
            assert(rhs_force_.size() == idx_free_.size());
            for (std::size_t i = 0; i < idx_free_.size(); ++i) {
                rhs[i] += rhs_pressure[idx_free_[i]];
            }
        } else {
            assert(rhs_force_.size() == grid_.leafGridView().size(3) * 3);
            assert(rhs.size() == rhs_pressure.size());
            for (std::size_t i = 0; i < rhs_pressure.size(); ++i) {
                rhs[i] += rhs_pressure[i];
            }
        }

        // end initialization
        auto& b = A.getLoadVector();
        b.resize(rhs.size());

        for (std::size_t i = 0; i < rhs.size(); ++i) {
            b[i] = rhs[i];
        }
    }
}

template <typename GridType>
void
VemElasticitySolver<GridType>::solve()
{
    try {
        Dune::InverseOperatorResult r;
        Vector& rhs = A.getLoadVector();
        u.resize(rhs.size());
        u = 0;
        tsolver_->apply(u, rhs, r);
    } catch (Dune::ISTLError& e) {
        std::cerr << "exception thrown " << e << std::endl;
    }
}

} // namespace Opm::Elasticity

#endif // OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP
