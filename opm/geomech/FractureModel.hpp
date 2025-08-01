/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_FRACTUREMODEL_HPP_INCLUDED
#define OPM_FRACTUREMODEL_HPP_INCLUDED

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>
#include <opm/grid/utility/compressedToCartesian.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/FractureWell.hpp>
#include <opm/geomech/GeometryHelpers.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Opm
{
class RuntimePerforation;

class ScheduleState;

template <typename Scalar>
class WellState;
} // namespace Opm

namespace Opm
{
PropertyTree makeDefaultFractureParam();

class FractureModel
{
public:
    using Point3D = Dune::FieldVector<double, 3>;
    using Segment = std::array<unsigned int, 2>;

    template <class Grid>
    FractureModel(const Grid& grid, const std::vector<Well>& wells, const PropertyTree&);

    /// Initialise fracture objects.
    ///
    /// Initialisation method selected by the property tree's "type" node,
    /// which defaults to the "well_seed" method (WSEED keyword).
    ///
    /// \param[in] sched Dynamic objects in current run, especially
    /// including well fracturing seed points and fracturing plane normal
    /// vectors in addition to all current well objects.
    void addFractures(const ScheduleState& sched);

    void updateFractureReservoirCells()
    {
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cell_search_tree_);
            }
        }
    }

    // void updateFractureReservoirCells(const Dune::CpGrid& cpgrid)
    template <class Grid>
    void updateFractureReservoirCells(const Grid& cpgrid)
    {
        external::cvf::ref<external::cvf::BoundingBoxTree> cellSearchTree;
        external::buildBoundingBoxTree(cellSearchTree, cpgrid);

        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cellSearchTree, cpgrid);
            }
        }
    }

    void addWell(const std::string& name,
                 const std::vector<FractureWell::Connection>& conns,
                 const std::vector<Point3D>& points,
                 const std::vector<std::array<unsigned, 2>>& segments);

    void write(int ReportStep = -1) const;
    void writemulti(double time) const;

    template <class TypeTag, class Simulator>
    void solve(const Simulator& simulator)
    {
        for (auto& fractures : this->well_fractures_) {
            for (auto& fracture : fractures) {
                if (fracture.isActive()) {
                    std::cout << "Solving fracture " << fracture.name() << std::endl;
                    fracture.template solve<TypeTag>(cell_search_tree_, simulator);
                }
            }
        }
    }

    void updateReservoirProperties();
    void initFractureStates();

    template <class TypeTag, class Simulator>
    void initReservoirProperties(const Simulator& simulator)
    {
        for (auto& fractures : this->well_fractures_) {
            for (auto& fracture : fractures) {
                fracture.template updateReservoirProperties<TypeTag>(simulator, true);
            }
        }
    }

    template <class TypeTag, class Simulator>
    void updateReservoirAndWellProperties(const Simulator& simulator)
    {
        this->updateReservoirProperties<TypeTag>(simulator);
        this->updateWellProperties<TypeTag>(simulator);
    }

    template <class TypeTag, class Simulator>
    void updateReservoirWellProperties(const Simulator& simulator)
    {
        for (auto& well : wells_) {
            well.template updateReservoirProperties<TypeTag>(simulator);
        }
    }

    std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname) const;

    template <typename Scalar>
    void assignGeomechWellState(WellState<Scalar>& wellState) const;

    bool addPertsToSchedule()
    {
        return prm_.get<bool>("addperfs_to_schedule");
    }

    // probably this should be collected in one loop since all do full loop over fracture
    // ++ well
    Dune::FieldVector<double, 6> stress(const Dune::FieldVector<double, 3>& obs) const;
    Dune::FieldVector<double, 6> strain(const Dune::FieldVector<double, 3>& obs) const;
    Dune::FieldVector<double, 3> disp(const Dune::FieldVector<double, 3>& obs) const;

    PropertyTree& getParam()
    {
        return prm_;
    }

private:
    bool vtkwritewells_ = false; // write wells to VTK files

    template <class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator)
    {
        for (auto& fractures : this->well_fractures_) {
            for (auto& fracture : fractures) {
                fracture.template updateReservoirProperties<TypeTag>(simulator);
            }
        }
    }

    template <class TypeTag, class Simulator>
    void updateWellProperties(const Simulator& simulator); // update all well related properties

    std::vector<FractureWell> wells_;
    std::vector<std::vector<Fracture>> well_fractures_;
    PropertyTree prm_;
    external::cvf::ref<external::cvf::BoundingBoxTree> cell_search_tree_;

    /// Initialise fractures perpendicularly to each reservoir connection.
    void addFracturesPerpWell();
    void addFracturesTensile();

    /// Initialise fractures in each seed identified in the WSEED keyword.
    ///
    /// \param[in] sched Dynamic objects in current run, especially
    /// including well fracturing seed points and fracturing plane normal
    /// vectors in addition to all current well objects.
    void addFracturesWellSeed(const ScheduleState& sched);
};

} // namespace Opm

#include "FractureModel_impl.hpp"

#endif // OPM_FRACTUREMODEL_HPP_INCLUDED
