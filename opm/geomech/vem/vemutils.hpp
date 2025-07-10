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

#ifndef VEM_UTILS_HPP
#define VEM_UTILS_HPP

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>

#ifdef HAVE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <opm/geomech/boundaryutils.hh>

#include <array>
#include <vector>

namespace vem
{
using PolyGrid = Dune::PolyhedralGrid<3, 3>;

#ifdef HAVE_ALUGRID
using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
#endif

void getGridVectors(const PolyGrid& grid,
                    std::vector<double>& coords,
                    std::vector<int>& num_cell_faces,
                    std::vector<int>& num_face_corners,
                    std::vector<int>& face_corners);

void getGridVectorsDune(const PolyGrid& grid,
                        std::vector<double>& coords,
                        std::vector<int>& num_cell_faces,
                        std::vector<int>& num_face_corners,
                        std::vector<int>& face_corners);

template <class GridType>
void
getGridVectors(const GridType& grid,
               std::vector<double>& coords,
               std::vector<int>& num_cell_faces,
               std::vector<int>& num_face_corners,
               std::vector<int>& face_corners)
{
    static constexpr int dim = GridType::dimension;
    const auto& gv = grid.leafGridView();

    for (const auto& v : vertices(gv)) {
        const auto c = v.geometry().corner(0);
        coords.insert(coords.end(), c.begin(), c.end());
    }

    for (const auto& cell : elements(gv)) {
        num_cell_faces.push_back(6);

        for (int i = 0; i < 6; ++i) {
            const auto faceDir = Opm::Elasticity::faceToFaceDir(i);
            const auto nodes = Opm::Elasticity::faceDirToNodes(faceDir);
            num_face_corners.push_back(nodes.size());

            for (const auto& nind : nodes) {
                const auto global_ind = gv.indexSet().subIndex(cell, nind, dim);
                face_corners.push_back(global_ind);
            }
        }
    }
}

void getGridVectors(const Dune::CpGrid& grid,
                    std::vector<double>& coords,
                    std::vector<int>& num_cell_faces,
                    std::vector<int>& num_face_corners,
                    std::vector<int>& face_corners);

} // namespace vem

#endif // VEM_UTILS_HPP
