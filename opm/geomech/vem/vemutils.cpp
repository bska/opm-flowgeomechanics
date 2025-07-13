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

#include <opm/geomech/vem/vemutils.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace vem
{
using PolyGrid = Dune::PolyhedralGrid<3, 3>;

#ifdef HAVE_ALUGRID
using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
#endif

void
getGridVectors(const PolyGrid& grid,
               std::vector<double>& coords,
               std::vector<int>& num_cell_faces,
               std::vector<int>& num_face_corners,
               std::vector<int>& face_corners)
{
    const UnstructuredGrid& ungrid = grid;
    static constexpr int dim = PolyGrid::dimension;

    coords.resize(ungrid.number_of_nodes * dim);
    for (int i = 0; i < ungrid.number_of_nodes; ++i) {
        for (int j = 0; j < dim; ++j) {
            coords[3 * i + j] = ungrid.node_coordinates[3 * i + j];
        }
    }

    num_cell_faces.resize(ungrid.number_of_cells);
    for (int i = 0; i < ungrid.number_of_cells; ++i) {
        num_cell_faces[i] = ungrid.cell_facepos[i + 1] - ungrid.cell_facepos[i];
    }

    num_face_corners.resize(ungrid.cell_facepos[ungrid.number_of_cells]);

    int tot_num_face_corners = 0;
    for (int cell = 0; cell < ungrid.number_of_cells; ++cell) {
        for (grid_size_t hface = ungrid.cell_facepos[cell]; hface < ungrid.cell_facepos[cell + 1];
             hface++) {
            int face = ungrid.cell_faces[hface]; // ungrid.cell_facepos[i]];
            int num_local_corners = ungrid.face_nodepos[face + 1] - ungrid.face_nodepos[face];
            num_face_corners[hface] = num_local_corners;
            // NB maybe we should order with outwards normal
            if (cell == ungrid.face_cells[2 * face]) {
                for (grid_size_t j = ungrid.face_nodepos[face]; j < ungrid.face_nodepos[face + 1]; ++j) {
                    face_corners.push_back(ungrid.face_nodes[j]);
                }
            } else {
                // flip orientation for hface
                for (grid_size_t j = ungrid.face_nodepos[face + 1] - 1; j >= ungrid.face_nodepos[face];
                     --j) {
                    face_corners.push_back(ungrid.face_nodes[j]);
                }
            }

            tot_num_face_corners += num_local_corners;
        }
    }

    assert(face_corners.size() == static_cast<std::size_t>(tot_num_face_corners));
}

void
getGridVectorsDune(const PolyGrid& grid,
                   std::vector<double>& coords,
                   std::vector<int>& num_cell_faces,
                   std::vector<int>& num_face_corners,
                   std::vector<int>& face_corners)
{
    const auto& gv = grid.leafGridView();

    // start VEM assembly
    // make global point coordinate vector
    for (const auto& v : vertices(gv)) {
        const auto c = v.geometry().corner(0);
        coords.insert(coords.end(), c.begin(), c.end());
    }

    const auto& ixset = gv.indexSet();

    // count cell faces
    for (const auto& c : elements(gv)) {
        num_cell_faces.push_back(Dune::subEntities(c, Dune::Codim<1> {}).size());
    }

    const int tot_num_cfaces = std::accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

    // count face corners
    for (const auto& c : elements(gv)) {
        for (const auto& f : Dune::subEntities(c, Dune::Codim<1> {})) {
            num_face_corners.push_back(f.geometry().corners());
        }
    }

    // establish all face corners
    for (const auto& c : elements(gv)) {
        for (const auto& f : Dune::subEntities(c, Dune::Codim<1> {})) {
            for (int i = 0, f_ix = 0; f_ix != tot_num_cfaces && i != num_face_corners[f_ix];
                 ++i, f_ix += (i == num_face_corners[f_ix])) {
                face_corners.push_back(ixset.subIndex(f, i, 3));
            }
        }
    }

    // correct order of nodes in each face, to ensure they are mentioned in
    // clockwise or counterclockwise order. @@ This is a hack that might only work
    // for 4-faces!

    // body force

    // dirichlet boundary conditions
}



void
getGridVectors(const Dune::CpGrid& grid,
               std::vector<double>& coords,
               std::vector<int>& num_cell_faces,
               std::vector<int>& num_face_corners,
               std::vector<int>& face_corners)
{
    const auto& gv = grid.leafGridView();

    for (const auto& v : vertices(gv)) {
        const auto c = v.geometry().corner(0);
        coords.insert(coords.end(), c.begin(), c.end());
    }

    for (const auto& cell : elements(gv)) {
        const int cellIdx = gv.indexSet().index(cell);
        const int nf = grid.numCellFaces(cellIdx);
        num_cell_faces.push_back(nf);

        for (int f = 0; f < nf; ++f) {
            const auto face = grid.cellFace(cellIdx, f);
            const auto faceSize = grid.numFaceVertices(face);

            num_face_corners.push_back(faceSize);
            const auto out_cell = grid.faceCell(face, 1);
            const auto in_cell = grid.faceCell(face, 0);

            if (out_cell != cellIdx) {
                assert(in_cell == cellIdx);
                for (int v = 0; v < faceSize; ++v) {
                    const int fv = grid.faceVertex(face, v);
                    face_corners.push_back(fv);
                }
            } else {
                for (int v = faceSize - 1; v > -1; --v) {
                    const int fv = grid.faceVertex(face, v);
                    face_corners.push_back(fv);
                }
            }
        }
    }
}

} // namespace vem
