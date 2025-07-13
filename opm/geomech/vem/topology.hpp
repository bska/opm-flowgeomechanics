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

#ifndef VEM_TOPOLOGY_HPP_INCLUDED
#define VEM_TOPOLOGY_HPP_INCLUDED

#include <array>
#include <tuple>
#include <vector>

namespace vem
{
using IndexPair = std::array<int, 2>;

// Identify cell faces with individual cells and their faces.  Each entry in
// the return vector represents a cell face.  Its first element is the index
// of the cell, and the second element is the index of the face.
std::vector<IndexPair> cellfaces_cells_faces(const int num_cells, const int* const num_cell_faces);

// identify matching cell faces
std::tuple<std::vector<IndexPair>, std::vector<int>>
cellfaces_matching_faces(const int num_cells,
                         const int* const num_cell_faces,
                         const int* const num_face_corners,
                         const int* const face_corners);

std::vector<std::array<double, 3>> cellface_centroids(const double* const coords,
                                                      const std::vector<int>& cellface_ixs,
                                                      const int* const num_face_corners,
                                                      const int* const face_corners);

std::vector<std::tuple<double, IndexPair>>
mutual_distances(const std::vector<std::array<double, 3>>& points);

// return vector with indices to the cellfaces considered to be 'top' and 'bottom'
// for each cell.
std::vector<std::array<int, 2>> identify_top_bottom_faces(const double* const coords,
                                                          const int num_cells,
                                                          const int* const num_cell_faces,
                                                          const int* const num_face_corners,
                                                          const int* const face_corners);


} // namespace vem

#endif // VEM_TOPOLOGY_HPP_INCLUDED
