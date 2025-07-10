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

#include <config.h>

#include <opm/geomech/GeometryHelpers.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace external
{

std::vector<std::size_t>
findCloseCellIndices(const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, cvf::BoundingBox& bb)
{
    std::vector<std::size_t> closeCells;

    m_cellSearchTree->findIntersections(bb, &closeCells);

    return closeCells;
}

void
buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Dune::CpGrid& grid)
{
    using GridView = Dune::CpGrid::LeafGridView;
    const auto& gv = grid.leafGridView();

    std::vector<std::size_t> cellIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> cellBoundingBoxes;

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    const ElementMapper mapper(gv, Dune::mcmgElementLayout()); // used id sets interally

    for (const auto& element : Dune::elements(gv)) {
        const int index = mapper.index(element);
        const auto& geom = element.geometry();

        assert(geom.corners() == 8);

        cvf::BoundingBox cellBB;
        cvf::Vec3d cornerPoint;

        // NB order should not matter when adding to bounding box: dune ordring and
        // resinsight ordering is different
        //  dune 0 1 2 3 4 5 6 7 is resinsight 0 1 3 2 4 5 7 6 (i think)
        for (std::size_t l = 0; l < 8; ++l) {
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

void
buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Opm::EclipseGrid& m_grid)
{
    const auto nx = m_grid.getNX();
    const auto ny = m_grid.getNY();
    const auto nz = m_grid.getNZ();

    const std::size_t cellCount = nx * ny * nz;

    std::vector<std::size_t> cellIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> cellBoundingBoxes;

    const std::size_t threadCellCount = cellCount;

    std::vector<std::size_t> threadIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> threadBoundingBoxes;

    threadIndicesForBoundingBoxes.reserve(threadCellCount);
    threadBoundingBoxes.reserve(threadCellCount);

    std::array<double, 3> cornerPointArray;
    cvf::Vec3d cornerPoint;

    for (int cIdx = 0; cIdx < (int)cellCount; ++cIdx) {
        const auto& [i, j, k] = m_grid.getIJK(cIdx);

        cvf::BoundingBox cellBB;

        for (std::size_t l = 0; l < 8; ++l) {
            cornerPointArray = m_grid.getCornerPos(i, j, k, l);
            cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
            cellBB.add(cornerPoint);
        }

        if (cellBB.isValid()) {
            threadIndicesForBoundingBoxes.emplace_back(cIdx);
            threadBoundingBoxes.emplace_back(cellBB);
        }
    }

    threadIndicesForBoundingBoxes.shrink_to_fit();
    threadBoundingBoxes.shrink_to_fit();

    cellIndicesForBoundingBoxes.insert(cellIndicesForBoundingBoxes.end(),
                                       threadIndicesForBoundingBoxes.begin(),
                                       threadIndicesForBoundingBoxes.end());

    cellBoundingBoxes.insert(
        cellBoundingBoxes.end(), threadBoundingBoxes.begin(), threadBoundingBoxes.end());

    m_cellSearchTree = new cvf::BoundingBoxTree;
    m_cellSearchTree->buildTreeFromBoundingBoxes(cellBoundingBoxes, &cellIndicesForBoundingBoxes);
}

} // namespace external
