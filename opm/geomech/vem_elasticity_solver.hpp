//==============================================================================
//!
//! \file elasticity_upscale.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class
//!
//==============================================================================

#ifndef VEM_ELASTICITY_SOLVER_HPP_
#define VEM_ELASTICITY_SOLVER_HPP_

#include <dune/common/fmatrix.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/ilu.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <opm/common/TimingMacros.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>
#include <opm/elasticity/shapefunctions.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include <opm/elasticity/asmhandler.hpp>
#include <opm/elasticity/boundarygrid.hh>
#include <opm/elasticity/elasticity.hpp>
#include <opm/elasticity/elasticity_preconditioners.hpp>
#include <opm/elasticity/logutils.hpp>
#include <opm/elasticity/materials.hh>
#include <opm/elasticity/matrixops.hpp>
#include <opm/elasticity/meshcolorizer.hpp>
#include <opm/elasticity/mortar_evaluator.hpp>
#include <opm/elasticity/mortar_schur.hpp>
#include <opm/elasticity/mortar_schur_precond.hpp>
#include <opm/elasticity/mortar_utils.hpp>
#include <opm/elasticity/mpc.hh>
#include <opm/elasticity/uzawa_solver.hpp>

#include <opm/geomech/DuneCommunicationHelpers.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

namespace Opm::Elasticity
{
//! \brief The main driver class
template <class GridType>
class VemElasticitySolver
{
public:
    //! \brief Dimension of our grid
    static const int dim = GridType::dimension;

    //! \brief A basic number
    using ctype = typename GridType::LeafGridView::ctype;

    //! \brief A vectorial node value
    using NodeValue = Dune::FieldVector<double, dim>;

    //! \brief A global coordinate
    using GlobalCoordinate =
        typename GridType::LeafGridView::template Codim<1>::Geometry::GlobalCoordinate;

    //! \brief A set of indices
    using LeafIndexSet = typename GridType::LeafGridView::IndexSet;

    //! \brief An iterator over grid cells
    using LeafIterator = typename GridType::LeafGridView::template Codim<0>::Iterator;

    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;

    //! \brief The linear operator
    ASMHandler<GridType> A; // NB NB names have to change

    // Matrix stiffnessMatrix_;
    //! \brief The solution vectors
    Vector u; // NB NB names have to change

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif

    //! \brief Main constructor
    //! \param[in] gv_ The grid to operate on
    //! \param[in] tol_ The tolerance to use when deciding whether or not a coordinate
    //! falls on a plane/line/point. \sa tol \param[in] Escale_ A scale value for
    //! E-moduluses to avoid numerical issues \param[in] file The eclipse grid file
    //! \param[in] rocklist If not blank, file is a rocklist \param[in] verbose If true,
    //! give verbose output
    VemElasticitySolver(const GridType& grid)
        : A(grid)
        , grid_(grid)
    {
    }

    void setMaterial(const std::vector<std::shared_ptr<Material>>& materials_)
    {
        materials = materials_;
    }

    void setMaterial(const std::vector<double>& ymodule, const std::vector<double>& pratio)
    {
        ymodule_ = ymodule;
        pratio_ = pratio;
    }

    //! \brief Find boundary coordinates
    //! \param[out] min The miminum coordinates of the grid
    //! \param[out] max The maximum coordinates of the grid
    void initForAssembly()
    {
    }

    void fixNodes(const std::vector<std::tuple<std::size_t, MechBCValue>>& bc_nodes)
    {
        // void fixNodesVem(const std::vector<std::size_t>& fixed_nodes){
        std::vector<int> fixed_dof_ixs;
        std::vector<double> fixed_dof_values;

        for (int i = 0; i < int(bc_nodes.size()); ++i) {
            const int node = std::get<0>(bc_nodes[i]);
            const auto& mechbcvalue = std::get<1>(bc_nodes[i]);
            const auto& mask = mechbcvalue.fixeddir;
            const auto& disp = mechbcvalue.disp;

            for (int m = 0; m < 3; ++m) {
                if (mask[m]) {
                    fixed_dof_ixs.push_back(3 * node + m);
                    fixed_dof_values.push_back(disp[m]);
                }
            }
        }

        const int num_fixed_dofs = (int)fixed_dof_ixs.size();

        dirichlet_ = {num_fixed_dofs, fixed_dof_ixs, fixed_dof_values};
    }

    void expandSolution(Vector& result, const Vector& U)
    {
        std::fill(result.begin(), result.end(), std::nan("1"));

        const auto& fixed_dof_values = std::get<2>(dirichlet_);
        const auto& fixed_dof_ixs = std::get<1>(dirichlet_);

        for (int i = 0; i != int(fixed_dof_ixs.size()); ++i) {
            result[fixed_dof_ixs[i]] = fixed_dof_values[i];
        }

        for (int i = 0, cur_ix = 0; i != int(result.size()); ++i) {
            if (std::isnan(result[i][0])) {
                result[i] = U[cur_ix++];
            }
        }
    }

    //! \brief Assemble (optionally) stiffness matrix A and load vector
    //! \param[in] loadcase The strain load case. Set to -1 to skip
    //! \param[in] matrix Whether or not to assemble the matrix
    void assemble(const Vector& pressure, bool matrix, bool vector, bool reduce_system);

    //! \brief Solve Au = b for u
    //! \param[in] loadcase The load case to solve
    void solve();

    const CommunicationType* comm() const
    {
        return comm_.get();
    }

    void setCopyRowsToZero()
    {
        auto& matrix = A.getOperator();

        assert(matrix.N() == matrix.M());
        assert(matrix.N() == dofParallelIndexSet_->size());

        for (const auto& index : *dofParallelIndexSet_) {
            const auto at = index.local().attribute();
            if (at == Dune::OwnerOverlapCopyAttributeSet::copy) {
                const auto lind = index.local().local();

                matrix[lind] = 0.0; // set full row to searo
                matrix[lind][lind] = 1.0; // set diagonal to 1 (a bit slow due to seach)
            }
        }
    }

    // //! \param[in] params The linear solver parameters
    void setupSolver(const Opm::PropertyTree& prm)
    {
        OPM_TIMEBLOCK(setupLinearSolver);

#if HAVE_MPI
        // grid_.comm() is something like collective communication
        // comm_ here is the entity communication

        comm_.reset(
            new CommunicationType(grid_.comm(), Dune::SolverCategory::overlapping, /*free*/ false));
#endif

        if (this->comm_->communicator().size() > 1) {
            vertexParallelIndexSet_.reset(
                new CommunicationType::ParallelIndexSet(Opm::makeEntityEntityCommunication<3>(grid_)));
            assert(grid_.size(3) == vertexParallelIndexSet_->size());
            dofParallelIndexSet_.reset(new CommunicationType::ParallelIndexSet(
                Opm::entityToDofIndexSet(*vertexParallelIndexSet_, 3)));
            // vertexRemoteIndexSet_.new(
            // RemotePallelIndexSet(vertex_parallindex,vertex_parallindex, grid_.comm()));
            // vertexRemoteIndexSet_->rebuid<false>(); auto dofParallelIndexSet =
            // Opm::entityToDofIndexSet(*vertexParallelIndexSet_, 3);
            comm_->indexSet() = *dofParallelIndexSet_;
            assert(grid_.size(3) * 3 == dofParallelIndexSet_->size());
            comm_->remoteIndices().rebuild<false>(); // = *vertexRemoteParallelIndexSet_;

            // OPM_THROW(std::runtime_error,"Parallel for mechanics not implemented");
            const std::function<Vector()> weightsCalculator;
            const std::size_t pressureIndex = 0;

            // NB NB ! need other operator in parallel
            setCopyRowsToZero();
            std::cout << "Make parallel solver" << std::endl;
            comm_->communicator().barrier();

            using ParOperatorType
                = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, CommunicationType>;
            this->pop_ = std::make_unique<ParOperatorType>(A.getOperator(), *comm_);

            using FlexibleSolverType = Dune::FlexibleSolver<ParOperatorType>;
            tsolver_ = std::make_unique<FlexibleSolverType>(
                *this->pop_, *comm_, prm, weightsCalculator, pressureIndex);
        } else {
            const std::size_t pressureIndex = 0; // Dummy
            const std::function<Vector()> weightsCalculator; // Dummy

            using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;

            this->sop_ = std::make_unique<SeqOperatorType>(A.getOperator());
            this->tsolver_ = std::make_unique<FlexibleSolverType>(
                *this->sop_, prm, weightsCalculator, pressureIndex);
        }
    }

    template <int comp>
    void averageStress(Dune::BlockVector<Dune::FieldVector<ctype, comp>>& sigmacells,
                       const Vector& uarg);

    void setBodyForce(double gravity)
    {
        OPM_TIMEBLOCK(setBodyFrorce);

        const int num_cells = grid_.leafGridView().size(0); // entities of codim 0

        // assemble the mechanical system
        body_force_.resize(3 * num_cells, 0.0);
        for (int i = 0; i < num_cells; ++i) {
            body_force_[3 * i + 2] = 2000 * gravity;
        }
    }

    static void makeDuneMatrixCompressed(const std::vector<std::tuple<int, int, double>>& A_entries,
                                         Matrix& mat)
    {
        OPM_TIMEBLOCK(makeDuneMatrixCompressed);
        {
            OPM_TIMEBLOCK(buildMatrixImplicite);

            // build mode need to be implicite and corredect matrix settings has to be
            // done mat = 0;
            for (const auto& matel : A_entries) {
                const int i = std::get<0>(matel);
                const int j = std::get<1>(matel);

                mat.entry(i, j) = 0;
            }

            mat.compress();
        }

        {
            OPM_TIMEBLOCK(setMatrixValues);

            mat = 0;

            for (const auto& matel : A_entries) {
                const int i = std::get<0>(matel);
                const int j = std::get<1>(matel);

                mat[i][j] += std::get<2>(matel);
            }
        }
    }

    void updateRhsWithGrad(const Vector& mechpot)
    {
        OPM_TIMEBLOCK(updateRhsWithGrad);

        Vector& b = A.getLoadVector();
        b = 0;
        b.resize(rhs_force_.size());
        divmat_.mv(mechpot, b);

        // end initialization
        for (std::size_t i = 0; i < rhs_force_.size(); ++i) {
            b[i] += rhs_force_[i];
        }
    }

    void makeDuneSystemMatrix(const std::vector<std::tuple<int, int, double>>& A_entries)
    {
        OPM_TIMEBLOCK(makeDuneSystemMatrix);

        {
            OPM_TIMEBLOCK(buildMatrixImplicite);

            int ncols = 0;
            int nrows = 0;
            for (const auto& matel : A_entries) {
                nrows = std::max(nrows, std::get<0>(matel));
                ncols = std::max(ncols, std::get<1>(matel));
            }

            nrows = nrows + 1;
            ncols = ncols + 1;

            Matrix& MAT = this->A.getOperator();
            MAT = 0;
            MAT.setBuildMode(Matrix::implicit);
            MAT.setImplicitBuildModeParameters(81, 0.4);
            MAT.setSize(nrows, ncols);

            makeDuneMatrixCompressed(A_entries, MAT);
        }
    }

    void calculateStressPrecomputed(const Vector& dispalldune);
    void calculateStrainPrecomputed(const Vector& dispalldune);
    void calculateStrain(); // bool precalculated);
    void calculateStress(); // bool precalculated);

    const Dune::BlockVector<Dune::FieldVector<ctype, 6>>& stress() const
    {
        return stress_;
    }

    const Dune::BlockVector<Dune::FieldVector<ctype, 6>>& strain() const
    {
        return strain_;
    }

private:
    void expandDisp(std::vector<double>& dispall, bool expand);
    void assignToVoigt(Dune::BlockVector<Dune::FieldVector<double, 6>>& voigt_stress,
                       const Dune::BlockVector<Dune::FieldVector<double, 1>>& vemstress);
    void assignToVoigtSymMat(Dune::BlockVector<Dune::FieldVector<double, 6>>& voigt_stress,
                             const std::vector<std::array<double, 6>>& vemstress);

    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;

    using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    std::unique_ptr<SeqOperatorType> sop_;

    using ParOperatorType = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, CommunicationType>;
    std::unique_ptr<ParOperatorType> pop_;

    //! \brief Linear solver
    using SolverPtr = std::unique_ptr<Dune::InverseOperator<Vector, Vector>>;
    SolverPtr tsolver_;

    //! \brief An iterator over grid vertices
    using LeafVertexIterator = typename GridType::LeafGridView::template Codim<dim>::Iterator;

    //! \brief A reference to our grid
    const GridType& grid_;

    int num_cells_ {};
    std::vector<double> coords_ {};
    std::vector<int> num_cell_faces_, num_face_corners_, face_corners_ {};
    std::vector<int> idx_free_ {};

    std::vector<double> rhs_force_ {};

    //! \brief Vector holding material parameters for each active grid cell
    std::vector<std::shared_ptr<Material>> materials {};
    std::vector<double> ymodule_ {};
    std::vector<double> pratio_ {};
    std::vector<double> body_force_ {}; //= set_body_force(num_cells, bfcase)
    std::tuple<int, std::vector<int>, std::vector<double>>
        dirichlet_ {}; // set_dirichlet(coords, dircase)

    Dune::BlockVector<Dune::FieldVector<ctype, 6>> stress_;
    Dune::BlockVector<Dune::FieldVector<ctype, 6>> strain_;
    Matrix stressmat_; // from all dofs (not eliminating bc) to cell
    Matrix strainmat_;
    Matrix divmat_; // from cell pressure to active dofs

    std::shared_ptr<CommunicationType> comm_;
    std::shared_ptr<CommunicationType::ParallelIndexSet> vertexParallelIndexSet_;
    std::shared_ptr<CommunicationType::ParallelIndexSet> dofParallelIndexSet_;
};

} // namespace Opm::Elasticity

#include "vem_elasticity_solver_impl.hpp"

#endif // VEM_ELASTICITY_SOLVER_HPP_
