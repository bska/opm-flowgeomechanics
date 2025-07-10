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

#ifndef GEOMETRYHELPERS_HPP_INCLUDED
#define GEOMETRYHELPERS_HPP_INCLUDED

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp>

// for resinsight search

#include <external/resinsight/CommonCode/cvfStructGrid.h>
#include <external/resinsight/LibGeometry/cvfBoundingBox.h>
#include <external/resinsight/ReservoirDataModel/RigCellGeometryTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>
#include <external/resinsight/ReservoirDataModel/cvfGeometryTools.h>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>
#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractor.hpp>

#include <array>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace Opm
{
// copy from base vangaurd
class GeometryHelper
{
public:
    using Point3D = Dune::FieldVector<double, 3>;

    template <class Grid>
    explicit GeometryHelper(const Grid& grid)
    {
        this->init(grid);
    }

    int compressedIndex(const int cartesianCellIdx) const
    {
        // using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        auto index_pair = this->cartesianToCompressed_.find(cartesianCellIdx);

        return (index_pair == this->cartesianToCompressed_.end()) ? -1 : index_pair->second;
    }

    const Point3D& centroid(const int cell_index) const
    {
        return centroids_[cell_index];
    }

private:
    std::unordered_map<int, int> cartesianToCompressed_;
    std::vector<Point3D> centroids_;
    std::vector<bool> is_interior_;

    template <class Grid>
    void init(const Grid& grid)
    {
        const auto num_cells = grid.leafGridView().size(0);

        this->is_interior_.resize(num_cells);
        this->centroids_.resize(num_cells);

        const auto elemMapper
            = Dune::MultipleCodimMultipleGeomTypeMapper {grid.leafGridView(), Dune::mcmgElementLayout()};

        const auto cartmapper = Dune::CartesianIndexMapper {grid};

        for (const auto& element : elements(grid.leafGridView())) {
            const auto elemIdx = elemMapper.index(element);

            this->cartesianToCompressed_.insert_or_assign(cartmapper.cartesianIndex(elemIdx), elemIdx);
            this->is_interior_[elemIdx] = element.partitionType() == Dune::InteriorEntity;
            this->centroids_[elemIdx] = element.geometry().center();
        }
    }
};

template <class Grid>
int
findCell(const Grid& grid, const Dune::FieldVector<double, 3>& point)
{
    // Brute force linear search.
    const auto elemMapper
        = Dune::MultipleCodimMultipleGeomTypeMapper {grid.leafGridView(), Dune::mcmgElementLayout()};

    for (const auto& element : elements(grid.leafGridView())) {
        const auto elemIdx = elemMapper.index(element);
        const auto& geom = element.geometry();

        if (Dune::referenceElement(geom).checkInside(geom.local(point))) {
            return elemIdx;
        }
    }

    return -1;
}

} // namespace Opm

namespace external
{
std::vector<std::size_t> findCloseCellIndices(const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree,
                                              cvf::BoundingBox& bb);

// copy from RigEclipseWellLogExtractor
void buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree,
                          const Opm::EclipseGrid& m_grid);
void buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Dune::CpGrid& grid);

template <class Grid>
void
buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Grid& grid)
{
    using GridView = typename Grid::LeafGridView;

    const auto& gv = grid.leafGridView();
    const std::size_t cellCount = gv.size(0);

    std::vector<std::size_t> cellIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> cellBoundingBoxes;

    std::array<double, 3> cornerPointArray;
    cvf::Vec3d cornerPoint;

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    const ElementMapper mapper(gv, Dune::mcmgElementLayout()); // used id sets interally

    for (const auto& element : Dune::elements(gv)) {
        const int index = mapper.index(element);
        const auto& geom = element.geometry();

        cvf::BoundingBox cellBB;
        cvf::Vec3d cornerPoint;

        // NB order should not matter when adding to bounding box: dune ordring and
        // resinsight ordering is different
        //  dune 0 1 2 3 4 5 6 7 is resinsight 0 1 3 2 4 5 7 6 (i think)
        for (std::size_t l = 0; l < geom.corners(); ++l) {
            const auto cornerPointArray = geom.corner(l);
            cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
            cellBB.add(cornerPoint);
        }

        cellIndicesForBoundingBoxes.emplace_back(index);
        cellBoundingBoxes.emplace_back(cellBB);
    }

    m_cellSearchTree = new cvf::BoundingBoxTree;
    m_cellSearchTree->buildTreeFromBoundingBoxes(cellBoundingBoxes, &cellIndicesForBoundingBoxes);
}

} // namespace external

#endif // GEOMETRYHELPERS_HPP_INCLUDED
