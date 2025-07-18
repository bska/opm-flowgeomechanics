/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_FRACTURE_HH
#define OPM_FRACTURE_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/grid/CpGrid.hpp>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>

// for linear solve
#include <dune/istl/matrixmarket.hh>

#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/geomech/GeometryHelpers.hpp>
#include <opm/geomech/GridStretcher.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{
template <typename Scalar>
class ConnFracStatistics;
}

namespace Opm::Properties
{
template <class TypeTag, class MyType>
struct FluidSystem;

template <class TypeTag, class MyType>
struct NumPhases;
} // namespace Opm::Properties

namespace Opm
{
struct WellInfo
{
    std::string name {};
    int perf {};
    int well_cell {};
    int global_index {};
    int segment {};
    std::optional<std::pair<double, double>> perf_range {};
};

struct RuntimePerforation;

/// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
class Fracture
{
public:
    using Grid = Dune::FoamGrid<2, 3>;
    using Point3D = Dune::FieldVector<double, 3>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using DynamicMatrix = Dune::DynamicMatrix<double>;

    void init(const std::string& well,
              const int perf,
              const int well_cell,
              const int global_index,
              const int segment,
              const std::optional<std::pair<double, double>>& perf_range,
              const Point3D& origo,
              const Point3D& normal,
              const PropertyTree& prm);

    void grow(int layers, int method);
    std::string name() const;
    void write(int reportStep = -1) const;
    void writemulti(double time) const;
    void updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree);

    // solver related
    void updateReservoirProperties();
    void removeCells();

    template <class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator,
                                   bool init_constant_vals = false,
                                   bool update_filtercake = true);

    void updateFilterCakeProps(const Opm::WellConnections& connections,
                               const Opm::SingleWellState<double>& wellstate);
    void initFracturePressureFromReservoir();
    void initFractureStates();
    void initFractureWidth();
    void solveFractureWidth();
    void solvePressure();

    template <class TypeTag, class Simulator>
    void solve(const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
               const Simulator& simulator);

    void printPressureMatrix() const; // debug purposes
    void printMechMatrix() const; // debug purposes
    void writeFractureSystem() const;
    void writePressureSystem() const;
    void setFractureGrid(std::unique_ptr<Fracture::Grid> gptr
                         = nullptr); // a hack to allow use of another grid
    std::vector<RuntimePerforation> wellIndices() const;

    WellInfo& wellInfo()
    {
        return wellinfo_;
    }

    const WellInfo& wellInfo() const
    {
        return wellinfo_;
    }

    std::vector<double> leakOfRate() const;
    double injectionPressure() const;

    void setPerfPressure(double perfpressure)
    {
        perf_pressure_ = perfpressure;
    }

    Dune::FieldVector<double, 6> stress(const Dune::FieldVector<double, 3>& obs) const;
    Dune::FieldVector<double, 6> strain(const Dune::FieldVector<double, 3>& obs) const;
    Dune::FieldVector<double, 3> disp(const Dune::FieldVector<double, 3>& obs) const;

    template <typename Scalar>
    void assignGeomechWellState(ConnFracStatistics<Scalar>& stats) const;

    void setActive(bool active)
    {
        active_ = active;
    }

    bool isActive() const
    {
        return active_;
    }

private:
    std::vector<double> redistribute_values(const std::vector<double>& values,
                                            const std::vector<std::vector<CellRef>>& map1,
                                            const std::vector<std::vector<CellRef>>& map2,
                                            const int level,
                                            bool point_wise);

    void redistribute_values(Dune::BlockVector<Dune::FieldVector<double, 1>>& values,
                             const std::vector<std::vector<CellRef>>& map1,
                             const std::vector<std::vector<CellRef>>& map2,
                             const int level,
                             bool point_wise);

    using ResVector = ::Opm::Fracture::Vector;
    using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

    // using Krull = Dune::MultiTypeBlockVector<Dune::BlockVector<double>>;
    // Dune::MultiTypeBlockVector<Dune::BlockVector<double, std::allocator<double>>> dill;
    // Krull tull(2);

    using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
    using FMatrix = Dune::DynamicMatrix<double>; // full matrix

    std::size_t numFractureCells() const
    {
        return grid_->leafGridView().size(0);
    }

    std::vector<int>
    identify_closed(const FMatrix& A, const VectorHP& x, const ResVector& rhs, const int nwells);
    template <class TypeTag, class Simulator>
    void initReservoirProperties(const Simulator& simulator)
    {
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        const auto& problem = simulator.problem();

        // NB burde truleg interpolere
        // NB reservoir dist not calculated
        const std::size_t ncf = reservoir_cells_.size();
        reservoir_perm_.resize(ncf);
        reservoir_cstress_.resize(ncf);
        reservoir_density_.resize(ncf);

        // should be calcualted
        double dist = prm_.get<double>("reservoir.dist");
        reservoir_dist_.resize(ncf, dist);
        double numax = -1e99;
        double Emax = -1e99;
        assert(ncf > 0);

        for (std::size_t i = 0; i < ncf; ++i) {
            const int cell = reservoir_cells_[i];
            if (cell < 0) {
                // probably outside reservoir set all to zero
                reservoir_perm_[i] = 0.0;
                reservoir_cstress_[i] = 1e20;
                reservoir_density_[i] = 1000.0; // water density
                reservoir_dist_[i] = dist;
                continue;
            }

            auto normal = this->cell_normals_[i];
            {
                auto permmat = problem.intrinsicPermeability(cell);
                auto cstress = problem.cStress(cell);
                auto np = normal;
                permmat.mv(normal, np);
                double value = np.dot(normal);
                reservoir_perm_[i] = value;
                reservoir_cstress_[i] = cstress;
            }

            Emax = std::max(Emax, problem.yModule(cell));
            numax = std::max(numax, problem.pRatio(cell));
        }

        E_ = Emax;
        nu_ = numax;

        assert(E_ > 0);
        assert(nu_ > 0 && nu_ < 1);
    }

    void resetWriters();
    Dune::BlockVector<Dune::FieldVector<double, 3>> all_slips() const;
    // helpers for growing grid
    void insertLinear(const std::vector<unsigned int>& inner_indices);
    void insertExp(const std::vector<unsigned int>& inner_indices);
    void initFracture(); // create a new fracture grid from scratch
    Point3D surfaceMap(double x, double y);

    std::unique_ptr<GridStretcher> grid_stretcher_; //@@ experimental, for stretching grids
    std::unique_ptr<Opm::RegularTrimesh> trimesh_; // @@ experimental, implicitly defined grids

    std::vector<std::vector<CellRef>>
        grid_mesh_map_; // @@ index mapping from cells in grid_ to trimesh_.
                        // @@ NB: in general many-to-many
    std::unique_ptr<Grid> grid_;
    Point3D origo_;
    std::array<Point3D, 3> axis_;
    WellInfo wellinfo_;
    bool active_ {false}; // is fracture active?
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
    static constexpr auto VTKFormat = Dune::VTK::ascii;
    std::unique_ptr<::Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>> vtkmultiwriter_;
    std::vector<unsigned int> out_indices_;

    std::vector<double> filtercake_thikness_; // properties of filter cake, if any
    double filtercake_perm_ {0.0}; // permeability of filter cake, if any
    double filtercake_poro_ {0.0}; // permeability of filter cake, if any
    bool has_filtercake_ {false}; // if true, filter cake is used in the model

    // should probably not be needed
    int layers_;
    int nlinear_;

    // help function for solving
    void assemblePressure();
    void addSource();
    void initPressureMatrix();
    void setupPressureSolver();
    void updateFractureRHS();
    void updateLeakoff();
    void updateCellNormals();
    void normalFractureTraction(Dune::BlockVector<Dune::FieldVector<double, 1>>& traction,
                                bool resize = true) const;
    double normalFractureTraction(std::size_t ix) const;

    // one nonlinear iteration of fully coupled system.  Returns 'true' if converged
    bool fullSystemIteration(const double tol);

    void assembleFractureMatrix() const;
    std::vector<double> stressIntensityK1() const;
    int numWellEquations() const
    {
        return prm_.get_child("control").get<std::string>("type") == "rate_well" ? 1 : 0;
    }

    // double well_pressure_;// for now using prm object for definition
    std::vector<CellRef> well_source_cellref_; // references to well cells in the fully resolved TriMesh
    std::vector<int> well_source_; // indices to well cells
    // for reservoir
    std::vector<int> reservoir_cells_;
    // std::vector< Dune::FieldMatrix<double, 3, 3> > reservoir_perm_;
    std::vector<double> reservoir_perm_;
    std::vector<double> reservoir_cstress_;
    std::vector<double> reservoir_mobility_;
    std::vector<double> reservoir_density_;
    std::vector<double> reservoir_cell_z_;
    std::vector<double> reservoir_dist_;
    std::vector<double> reservoir_pressure_;
    std::vector<Dune::FieldVector<double, 6>> reservoir_stress_;

    // only for radom access need to be updater after trid change
    Dune::BlockVector<Dune::FieldVector<double, 3>> cell_normals_;

    // solution variables (only to avoid memory allocation, do not trust their state)
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_width_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> rhs_width_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_pressure_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> rhs_pressure_;

    // transmissibilities
    using Htrans = std::tuple<std::size_t, std::size_t, double, double>;
    std::vector<Htrans> htrans_;
    std::vector<std::tuple<int, double>> perfinj_;
    double perf_pressure_;
    std::vector<double> leakof_;
    //
    PropertyTree prmpressure_;
    using PressureOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;
    mutable std::unique_ptr<Matrix> pressure_matrix_;
    mutable std::unique_ptr<PressureOperatorType> pressure_operator_;
    mutable std::unique_ptr<FlexibleSolverType> pressure_solver_;
    mutable std::unique_ptr<Matrix> coupling_matrix_; // will be updated by `fullSystemIteration`

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
    // using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    // using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    // using DynamicMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    mutable std::unique_ptr<DynamicMatrix> fracture_matrix_;

    // function ensuring that the fracture matrix exists, and returning a reference to it
    DynamicMatrix& fractureMatrix() const
    {
        if (fracture_matrix_ == nullptr)
            assembleFractureMatrix();
        return *fracture_matrix_;
    }

    double E_;
    double nu_;
    double min_width_; // minimum width of fracture, used for convergence criterion
    double gravity_ {0.0}; //{9.81}; // gravity acceleration, used for leakoff calculations
    std::vector<double> fracture_dgh_; // gravity contribution to fracture pressure, used
                                       // for leakoff calculations
    PropertyTree prm_;
    double total_WI_well_ {0.0}; // total well index for the well, used for leakoff calculations
};

} // namespace Opm

#include "Fracture_impl.hpp"

#endif // OPM_FRACTURE_HH
