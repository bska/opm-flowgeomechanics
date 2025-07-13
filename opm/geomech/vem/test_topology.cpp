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

#include <opm/geomech/vem/topology.hpp>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

int
main([[maybe_unused]] int varnum, char** vararg)
{
    std::ifstream is(vararg[1]);

    // read number of cells and number of cell faces per cell
    int num_cells {};
    is >> num_cells;
    std::vector<int> num_cell_faces(num_cells);
    for (int i = 0; i < num_cells; ++i) {
        is >> num_cell_faces[i];
    }
    const int total_num_cellfaces = std::accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

    // read nodes
    int num_nodes {};
    is >> num_nodes;
    std::vector<double> nodes(num_nodes * 3);
    for (int i = 0; i < num_nodes * 3; ++i) {
        is >> nodes[i];
    }

    // read number of corners per cell face
    std::vector<int> num_face_corners(total_num_cellfaces);
    for (int i = 0; i < total_num_cellfaces; ++i) {
        is >> num_face_corners[i];
    }

    // read face corners
    const int total_num_face_corners
        = std::accumulate(num_face_corners.begin(), num_face_corners.end(), 0);
    std::vector<int> face_corners(total_num_face_corners);
    for (int i = 0; i < total_num_face_corners; ++i) {
        is >> face_corners[i];
    }

    // mention all cell face indices per default
    std::vector<std::size_t> cellface_ixs(total_num_cellfaces);
    for (auto i = 0 * total_num_cellfaces; i < total_num_cellfaces; ++i) {
        cellface_ixs[i] = static_cast<std::size_t>(i);
    }

    const auto tbfaces = vem::identify_top_bottom_faces(
        &nodes[0], num_cells, &num_cell_faces[0], &num_face_corners[0], &face_corners[0]);

    // print all tbfaces to screen
    std::cout << "----- Top bottom faces: ---\n";
    for (const auto e : tbfaces) {
        std::cout << e[0] << " " << e[1] << '\n';
    }

    std::cout.flush();

    return 0;
}
