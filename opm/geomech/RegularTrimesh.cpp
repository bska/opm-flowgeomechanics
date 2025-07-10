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

#include <opm/geomech/RegularTrimesh.hpp>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <opm/geomech/GeometryHelpers.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
const bool DEBUG_DUMP_GRIDS = false;
static int DEBUG_CURRENT_GRID_ITERATION_COUNT = 0; // @@
static int DEBUG_GRID_COUNT = 0;

std::array<Opm::EdgeRef, 3>
cell2edges(const Opm::CellRef& cell)
{
    return (cell[2] == 0) ? std::array<Opm::EdgeRef, 3> {Opm::EdgeRef {cell[0], cell[1], 0},
                                                         Opm::EdgeRef {cell[0], cell[1], 1},
                                                         Opm::EdgeRef {cell[0], cell[1], 2}}
                          : std::array<Opm::EdgeRef, 3> {Opm::EdgeRef {cell[0], cell[1] + 1, 0},
                                                         Opm::EdgeRef {cell[0] + 1, cell[1], 1},
                                                         Opm::EdgeRef {cell[0], cell[1], 2}};
}

std::array<Opm::CellRef, 2>
edge2cells(const Opm::EdgeRef& edge)
{
    return (edge[2] == 0) ? std::array<Opm::CellRef, 2> {Opm::CellRef {edge[0], edge[1], 0},
                                                         Opm::CellRef {edge[0], edge[1] - 1, 1}}
        : (edge[2] == 1)  ? std::array<Opm::CellRef, 2> {Opm::CellRef {edge[0], edge[1], 0},
                                                         Opm::CellRef {edge[0] - 1, edge[1], 1}}
                          : std::array<Opm::CellRef, 2> {Opm::CellRef {edge[0], edge[1], 1},
                                                         Opm::CellRef {edge[0], edge[1], 0}};
}

std::array<Opm::CellRef, 3>
cellNeighbors(const Opm::CellRef& cell)
{
    return (cell[2] == 0) ? std::array<Opm::CellRef, 3> {Opm::CellRef {cell[0], cell[1], 1},
                                                         Opm::CellRef {cell[0] - 1, cell[1], 1},
                                                         Opm::CellRef {cell[0], cell[1] - 1, 1}}
                          : std::array<Opm::CellRef, 3> {Opm::CellRef {cell[0], cell[1], 0},
                                                         Opm::CellRef {cell[0] + 1, cell[1], 0},
                                                         Opm::CellRef {cell[0], cell[1] + 1, 0}};
}

// += operator for Opm::Coord3D
Opm::Coord3D&
operator+=(Opm::Coord3D& lhs, const Opm::Coord3D& rhs)
{
    for (int i = 0; i != 3; ++i) {
        lhs[i] += rhs[i];
    }

    return lhs;
}

// *= operator for Opm::Coord3D
Opm::Coord3D&
operator*=(Opm::Coord3D& lhs, const double rhs)
{
    for (int i = 0; i != 3; ++i) {
        lhs[i] *= rhs;
    }

    return lhs;
}

// /= operator for Opm::Coord3D
Opm::Coord3D&
operator/=(Opm::Coord3D& lhs, const double rhs)
{
    for (int i = 0; i != 3; ++i) {
        lhs[i] /= rhs;
    }

    return lhs;
}

// + operator for Opm::Coord3D
Opm::Coord3D
operator+(const Opm::Coord3D& lhs, const Opm::Coord3D& rhs)
{
    Opm::Coord3D result(lhs);
    result += rhs;

    return result;
}

// / operator for Opm::Coord3D
Opm::Coord3D
operator/(const Opm::Coord3D& lhs, const double rhs)
{
    Opm::Coord3D result(lhs);
    result /= rhs;

    return result;
}

// * operator for Opm::Coord3D
Opm::Coord3D
operator*(const Opm::Coord3D& lhs, const double rhs)
{
    Opm::Coord3D result(lhs);
    result *= rhs;

    return result;
}

// * operator for Opm::Coord3D with double on lhs
Opm::Coord3D
operator*(const double lhs, const Opm::Coord3D& rhs)
{
    return rhs * lhs;
}

// == operator for Opm::NodeRef
bool
operator==(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs)
{
    return (lhs[0] == rhs[0]) && (lhs[1] == rhs[1]);
}

// < operator for Opm::NodeRef
bool
operator<(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs)
{
    return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
}

// == operator for Opm::EdgeRef
bool
operator==(const Opm::EdgeRef& lhs, const Opm::EdgeRef& rhs)
{
    return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
}

// < operator for Opm::EdgeRef
bool
operator<(const Opm::EdgeRef& lhs, const Opm::EdgeRef& rhs)
{
    return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1])
        || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]);
}

Opm::NodeRef
operator+(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs)
{
    return Opm::NodeRef {lhs[0] + rhs[0], lhs[1] + rhs[1]};
}

// ----------------------------------------------------------------------------
std::array<Opm::EdgeRef, 3>
half_edges(const Opm::CellRef& cell)
// ----------------------------------------------------------------------------
{
    return (cell[2] == 0) ? std::array<Opm::EdgeRef, 3> {Opm::EdgeRef {cell[0], cell[1], 0},
                                                         Opm::EdgeRef {cell[0], cell[1], 1},
                                                         Opm::EdgeRef {cell[0], cell[1], 2}}
                          : std::array<Opm::EdgeRef, 3> {Opm::EdgeRef {cell[0], cell[1] + 1, 0},
                                                         Opm::EdgeRef {cell[0] + 1, cell[1], 1},
                                                         Opm::EdgeRef {cell[0], cell[1], 2}};
}

// ----------------------------------------------------------------------------
std::array<Opm::CellRef, 3>
neigh_cells(const Opm::CellRef& cell)
// ----------------------------------------------------------------------------
{
    return (cell[2] == 0) ? std::array<Opm::CellRef, 3> {Opm::CellRef {cell[0], cell[1] - 1, 1},
                                                         Opm::CellRef {cell[0] - 1, cell[1], 1},
                                                         Opm::CellRef {cell[0], cell[1], 1}}
                          :

                          std::array<Opm::CellRef, 3> {Opm::CellRef {cell[0], cell[1] + 1, 0},
                                                       Opm::CellRef {cell[0] + 1, cell[1], 0},
                                                       Opm::CellRef {cell[0], cell[1], 0}};
}

// ----------------------------------------------------------------------------
std::array<bool, 3>
identify_boundary(const Opm::CellRef& cell, const Opm::RegularTrimesh& mesh)
// ----------------------------------------------------------------------------
{
    // return 'true' for the triangle edges that lie on the boundary
    const auto ncells = neigh_cells(cell);
    return {!mesh.isActive(ncells[0]), !mesh.isActive(ncells[1]), !mesh.isActive(ncells[2])};
}

// ----------------------------------------------------------------------------
Opm::RegularTrimesh
remove_mesh_corners(const Opm::RegularTrimesh& mesh)
// ----------------------------------------------------------------------------
{
    Opm::RegularTrimesh result(mesh);

    bool redo = true;

    while (redo) {
        redo = false;

        for (const auto& cell : result.boundaryCells()) {
            const auto bnd = identify_boundary(cell, result);
            if (std::count(bnd.begin(), bnd.end(), true) > 1) {
                result.setInactive(cell);
                redo = true;
            }
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<std::array<Opm::NodeRef, 3>>
single_split(const Opm::CellRef& cell, const std::array<bool, 3>& bnd)
// ----------------------------------------------------------------------------
{
    assert(bnd[0] + bnd[1] + bnd[2] == 1);

    const int i = cell[0];
    const int j = cell[1];

    const Opm::NodeRef A
        = (cell[2] == 0) ? Opm::NodeRef {2 * i, 2 * j} : Opm::NodeRef {2 * (i + 1), 2 * (j + 1)};
    const Opm::NodeRef B
        = (cell[2] == 0) ? Opm::NodeRef {2 * (i + 1), 2 * j} : Opm::NodeRef {2 * i, 2 * (j + 1)};
    const Opm::NodeRef C
        = (cell[2] == 0) ? Opm::NodeRef {2 * i, 2 * (j + 1)} : Opm::NodeRef {2 * (i + 1), 2 * j};

    const std::array<Opm::NodeRef, 3> tri = bnd[0] ? std::array<Opm::NodeRef, 3> {A, B, C}
        : bnd[1]                                   ? std::array<Opm::NodeRef, 3> {C, A, B}
                                                   : std::array<Opm::NodeRef, 3> {B, C, A};
    const Opm::NodeRef X = {(tri[0][0] + tri[1][0]) / 2, (tri[0][1] + tri[1][1]) / 2};

    return std::vector<std::array<Opm::NodeRef, 3>> {std::array<Opm::NodeRef, 3> {tri[0], X, tri[2]},
                                                     std::array<Opm::NodeRef, 3> {X, tri[1], tri[2]}};
}

// ----------------------------------------------------------------------------
std::vector<std::array<Opm::NodeRef, 3>>
tesselate_coarsecell(const Opm::CellRef& cell, const Opm::RegularTrimesh& mesh, const bool skip = false)
// ----------------------------------------------------------------------------
{
    if (skip) {
        const auto& nodes = mesh.cellNodes(cell);
        return std::vector<std::array<Opm::NodeRef, 3>> {
            std::array<Opm::NodeRef, 3> {Opm::RegularTrimesh::coarse_to_fine(nodes[0], 1),
                                         Opm::RegularTrimesh::coarse_to_fine(nodes[1], 1),
                                         Opm::RegularTrimesh::coarse_to_fine(nodes[2], 1)}};
    }

    // tessellate this cell as if its boundary intersects with a refined triangulation
    const std::array<bool, 3> bnd = identify_boundary(cell, mesh);

    const int sum_bnd = bnd[0] + bnd[1] + bnd[2];
    assert(sum_bnd <= 1); // before calling this function , the grid should have been
                          // 'rounded' so that there are no boundary cells with more than
                          // one boundary edge.

    if (sum_bnd == 1)
        // triangle split in the middle
        return single_split(cell, bnd);

    // if function is called on interior cell, just return the triangle associated
    // with that cell
    const auto& nodes = mesh.cellNodes(cell);
    return std::vector<std::array<Opm::NodeRef, 3>> {
        std::array<Opm::NodeRef, 3> {Opm::RegularTrimesh::coarse_to_fine(nodes[0]),
                                     Opm::RegularTrimesh::coarse_to_fine(nodes[1]),
                                     Opm::RegularTrimesh::coarse_to_fine(nodes[2])}};
}

// ----------------------------------------------------------------------------
bool
is_boundary_cell(const Opm::CellRef& cell, const Opm::RegularTrimesh& mesh)
// ----------------------------------------------------------------------------
{
    const auto bnd = identify_boundary(cell, mesh);
    return std::count(bnd.begin(), bnd.end(), true) > 0;
}

} // namespace

namespace Opm
{
// ----------------------------------------------------------------------------
RegularTrimesh::RegularTrimesh(const int layers,
                               const std::array<double, 3>& origin,
                               const std::array<double, 3>& axis1,
                               const std::array<double, 3>& axis2,
                               const std::array<double, 2>& edgelen)
    // ----------------------------------------------------------------------------
    : origin_(origin)
    , axis1_(RegularTrimesh::normalize(axis1))
    , axis2_(RegularTrimesh::normalize(axis2))
    , edgelen_(edgelen)
{
    cellinfo_[{0, 0, 0}] = CellAttributes {}; // set a single seed cell

    for (int i = 0; i != layers; ++i) {
        expandGrid();
    }
}

// ----------------------------------------------------------------------------
RegularTrimesh::RegularTrimesh(const double radius,
                               const std::array<double, 3>& origin,
                               const std::array<double, 3>& axis1,
                               const std::array<double, 3>& axis2,
                               const std::array<double, 2>& edgelen)
    // ----------------------------------------------------------------------------
    : origin_(origin)
    , axis1_(RegularTrimesh::normalize(axis1))
    , axis2_(RegularTrimesh::normalize(axis2))
    , edgelen_(edgelen)
{
    // NB to be checked
    const double scale = std::min(edgelen_[0], edgelen_[1]);
    const double R2 = (radius * radius) / (scale * scale);

    const double denom2
        = (5.0 / 4.0) - (axis1_[0] * axis2_[0] + axis1_[1] * axis2_[1] + axis1_[2] * axis2_[2]);

    const int c = std::ceil(std::max(std::sqrt(R2 / denom2), radius));

    for (int i = (-c - 1); i <= c; ++i) {
        for (int j = (-c - 1); j <= c; ++j) {
            for (int k = 0; k != 2; ++k) {
                if (dist(origin_, cellCentroid({i, j, k})) < radius) {
                    setActive(CellRef {i, j, k});
                }
            }
        }
    }

    auto tmp = remove_mesh_corners(*this);

    this->swap(tmp);
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::swap(RegularTrimesh& other)
// ----------------------------------------------------------------------------
{
    std::swap(cellinfo_, other.cellinfo_);
    std::swap(origin_, other.origin_);
    std::swap(axis1_, other.axis1_);
    std::swap(axis2_, other.axis2_);
    std::swap(edgelen_, other.edgelen_);
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::cellIndices() const
// ----------------------------------------------------------------------------
{
    std::vector<CellRef> indices;

    for (const auto& cell : cellinfo_) {
        indices.push_back(cell.first);
    }

    std::sort(indices.begin(), indices.end());

    return indices;
}

// ----------------------------------------------------------------------------
std::vector<EdgeRef>
RegularTrimesh::edgeIndices() const
// ----------------------------------------------------------------------------
{
    std::vector<EdgeRef> all_edges = all_half_edges_();

    // keep only unique elements
    const auto last = std::unique(all_edges.begin(), all_edges.end());

    // remove duplicates from all_edges
    all_edges.erase(last, all_edges.end());

    return all_edges;
}

// ----------------------------------------------------------------------------
std::vector<NodeRef>
RegularTrimesh::nodeIndices() const
// ----------------------------------------------------------------------------
{
    std::vector<NodeRef> indices;

    for (const auto& entry : cellinfo_) {
        const auto& cell = entry.first;
        for (int i = 0; i != 3; ++i) {
            indices.push_back(
                i == 0 ? (cell[2] == 0 ? NodeRef {cell[0], cell[1]} : NodeRef {cell[0] + 1, cell[1] + 1})
                    : i == 1 ? NodeRef {cell[0] + 1, cell[1]}
                             : NodeRef {cell[0], cell[1] + 1});
        }
    }

    std::sort(indices.begin(), indices.end());

    return {indices.begin(), std::unique(indices.begin(), indices.end())};
}

// ----------------------------------------------------------------------------
std::vector<EdgeRef>
RegularTrimesh::boundaryEdges() const
// ----------------------------------------------------------------------------
{
    // make a vector of all edges.  Count internal edges twice
    const std::vector<EdgeRef> all_edges = all_half_edges_();

    // boundary edges are those that are not duplicated
    std::vector<EdgeRef> result;
    for (auto it = all_edges.begin(); it != all_edges.end(); ++it) {
        if (it == all_edges.end() - 1 || *it != *(it + 1)) {
            result.push_back(*it);
        } else {
            ++it; // skip the next one, which we already know is duplicated
        }
    }

    // sort entries
    std::sort(result.begin(), result.end());

    return result;
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::boundaryCells() const
// ----------------------------------------------------------------------------
{
    std::vector<CellRef> result;
    for (const auto& edge : boundaryEdges()) {
        for (const auto& cell : edge2cells(edge)) {
            if (isActive(cell)) {
                result.push_back(cell);
            }
        }
    }

    // remove any duplicates
    std::sort(result.begin(), result.end());

    const auto last = std::unique(result.begin(), result.end());

    // shrink result to remove duplicates
    result.erase(last, result.end());

    return result;
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::interiorCells() const
// ----------------------------------------------------------------------------
{
    std::vector<CellRef> result;

    const std::vector<CellRef> bcells = boundaryCells();
    const std::vector<CellRef> allcells = cellIndices();

    std::set_difference(
        allcells.begin(), allcells.end(), bcells.begin(), bcells.end(), std::back_inserter(result));

    return result;
}

// ----------------------------------------------------------------------------
std::vector<std::pair<std::array<unsigned int, 3>, std::array<CellRef, 2>>>
RegularTrimesh::boundary_smoothing_triangles_() const
// ----------------------------------------------------------------------------
{
    // identify all 'internal' edges within boundary cells
    const auto bcells = boundaryCells();
    const auto bedges = boundaryEdges();

    // std::set<EdgeRef> bedges_set(bedges.begin(), bedges.end()); // better for search?
    std::array<std::set<EdgeRef>, 3> internal_edges;
    for (const auto& cell : bcells) {
        for (const auto& e : cell2edges(cell)) {
            // if (bedges_set.find(e) == bedges_set.end())
            internal_edges[e[2]].insert(e);
        }
    }

    // identify internal edges that 'line up' along one of the three cardinal grid
    // directions
    const std::array<int, 3> ioffsets {-1, 2, 1}, joffsets {2, -1, 1};
    std::array<std::vector<EdgeRef>, 3>
        candidate_sites; // candidates for where to place a smoothing triangle

    for (int i = 0; i != 3; ++i) {
        for (const auto& iedge : internal_edges[i]) {
            if (internal_edges[i].find(EdgeRef {iedge[0] + ioffsets[i], iedge[1] + joffsets[i], i})
                != internal_edges[i].end()) {
                candidate_sites[i].push_back(iedge);
            }
        }
    }

    // determine smoothing triangles
    std::vector<NodeRef> corners;
    std::vector<CellRef> neigh_cells;

    // neighbors that should be inactive in order to create a smoothing triangle
    // (template will be rotated depending on the direction considered)
    const std::array<CellRef, 6> check_template {
        {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}, {0, -1, 1}, {1, -1, 0}, {1, -2, 1}}};

    // template of smoothing triangle to add (will be rotated depending on the direction
    // considered)
    const std::array<NodeRef, 3> smooth_template {{{0, 0}, {1, -1}, {0, 1}}};
    const std::array<std::array<int, 2>, 3> ocell {
        {{-1, 1}, {1, -1}, {1, 1}}}; // location of 'opposing' cell
    const std::array<int, 3> dir_rot {
        0, 2, 1}; // how many times to rotate 60 degrees to align with corresp. triangle edge

    const auto cell_position
        = [](const EdgeRef& edge, const CellRef& template_cell, const int rotation) -> CellRef {
        CellRef result(template_cell);
        for (int i = 0; i != rotation; ++i) {
            rotate60(result);
        }

        for (int i = 0; i != 2; ++i) {
            result[i] += edge[i];
        }

        if (edge[2] != 0) {
            result[0] += 1;
        }

        if (edge[2] != 1) {
            result[1] += 1;
        }

        return result;
    };

    const auto node_rotate_n = [](NodeRef n, const int times) -> NodeRef {
        for (int i = 0; i != times; ++i) {
            rotate60(n);
        }

        return n;
    };

    // function to check that the cells from check_template are inactive
    const auto check_active = [&](const EdgeRef& edge, const int rotation) -> bool {
        for (const auto& template_cell : check_template) {
            if (isActive(cell_position(edge, template_cell, rotation % 6))) {
                return false;
            }
        }

        return true;
    };

    for (int dir = 0; dir != 3; ++dir) { // three cardinal directions
        for (const auto& edge : candidate_sites[dir]) {
            for (int side = 0; side != 2; ++side) { // right and left
                if (check_active(edge, dir_rot[dir] + 3 * side)) {
                    const CellRef n1 {
                        node2cell(edge2node(edge), dir == 2)}; // cell associated with first edge
                    const CellRef n2 {
                        node2cell(edge2node(edge) + ocell[dir], dir != 2)}; // opposing cell

                    if (!isActive(n1) || !isActive(n2)) {
                        continue;
                    }

                    // make smoothing triangle
                    for (int i = 0; i != 3; ++i) {
                        corners.push_back(
                            {edge2node(edge) + NodeRef {dir != 0, dir != 1}
                             + node_rotate_n(smooth_template[i], dir_rot[dir] + 3 * side)});
                    }

                    // keep track of cell neighbors to the new smoothing cell
                    neigh_cells.push_back(n1);
                    neigh_cells.push_back(n2);
                }
            }
        }
    }

    // prepare result by changing NodeRefs to indices
    const std::vector<unsigned int> node_ixs = noderefs_to_indices_(corners);

    assert(node_ixs.size()
           == 3 * (neigh_cells.size() / 2)); // each traingle has 3 nodes and 2 cell neighbors

    const int N = node_ixs.size() / 3;
    std::vector<std::pair<std::array<unsigned int, 3>, std::array<CellRef, 2>>> result;
    for (int i = 0; i != N; ++i) {
        result.push_back({{node_ixs[3 * i], // triangle corner 1
                           node_ixs[3 * i + 1], // triangle corner 2
                           node_ixs[3 * i + 2]}, // triangle corner 3
                          {neigh_cells[2 * i], // cell neighbor 1
                           neigh_cells[2 * i + 1]}}); // cell neighbor 2
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<EdgeRef>
RegularTrimesh::all_half_edges_() const
// ----------------------------------------------------------------------------
{
    // make a vector of all edges.  Count internal edges twice.  Sort the result
    std::vector<EdgeRef> all_edges;

    for (const auto& entry : cellinfo_) {
        for (const auto& edge : half_edges(entry.first)) {
            all_edges.push_back(edge);
        }
    }

    std::sort(all_edges.begin(), all_edges.end());

    return all_edges;
}

// ----------------------------------------------------------------------------
Coord3D
RegularTrimesh::nodeCoord(const NodeRef& node) const
// ----------------------------------------------------------------------------
{
    Coord3D result(origin_);

    result += node[0] * edgelen_[0] * axis1_;
    result += node[1] * edgelen_[1] * axis2_;

    return result;
}

// ----------------------------------------------------------------------------
Coord3D
RegularTrimesh::cellCentroid(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
    Coord3D result
        = cell[2] == 0 ? nodeCoord({cell[0], cell[1]}) : nodeCoord({cell[0] + 1, cell[1] + 1});

    result += nodeCoord({cell[0] + 1, cell[1]});
    result += nodeCoord({cell[0], cell[1] + 1});

    result /= 3;

    return result;
}

// ----------------------------------------------------------------------------
Coord3D
RegularTrimesh::edgeCentroid(const EdgeRef& edge) const
// ----------------------------------------------------------------------------
{
    return (edge[2] == 0) ? (nodeCoord({edge[0], edge[1]}) + nodeCoord({edge[0] + 1, edge[1]})) / 2
        : (edge[2] == 1)  ? (nodeCoord({edge[0], edge[1]}) + nodeCoord({edge[0], edge[1] + 1})) / 2
                          : (nodeCoord({edge[0], edge[1] + 1}) + nodeCoord({edge[0] + 1, edge[1]})) / 2;
}

// ----------------------------------------------------------------------------
std::vector<Coord3D>
RegularTrimesh::cellCentroids() const
// ----------------------------------------------------------------------------
{
    std::vector<Coord3D> result;

    for (const auto& entry : cellinfo_) {
        result.push_back(cellCentroid(entry.first));
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<Coord3D>
RegularTrimesh::edgeCentroids() const
// ----------------------------------------------------------------------------
{
    std::vector<Coord3D> result;

    for (const auto& edge : edgeIndices()) {
        result.push_back(edgeCentroid(edge));
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<Coord3D>
RegularTrimesh::nodeCoords() const
// ----------------------------------------------------------------------------
{
    std::vector<Coord3D> result;

    for (const auto& node : nodeIndices()) {
        result.push_back(nodeCoord(node));
    }

    return result;
}

//------------------------------------------------------------------------------
std::pair<std::vector<std::array<NodeRef, 3>>, std::vector<std::vector<CellRef>>>
RegularTrimesh::getTriangles() const
//------------------------------------------------------------------------------
{
    std::vector<std::array<NodeRef, 3>> result;
    std::vector<std::vector<CellRef>> trivialmap;

    const auto cells = cellIndices();
    for (const auto& cref : cells) {
        result.push_back(cellNodes(cref));
        trivialmap.push_back(std::vector<CellRef> {cref});
    }

    return std::make_pair(result, trivialmap);
}

//------------------------------------------------------------------------------
std::pair<std::vector<std::array<NodeRef, 3>>, std::vector<std::vector<CellRef>>>
RegularTrimesh::getMultiresTriangles(const std::vector<CellRef>& fixed_cells,
                                     const int max_levels,
                                     const int cellnum_threshold) const
//------------------------------------------------------------------------------
{
    RegularTrimesh mesh(*this);

    const int CELLNUM_THRESHOLD = cellnum_threshold;

    std::vector<std::array<NodeRef, 3>> result_triangles;
    std::vector<std::vector<CellRef>> result_cells;

    int level = 0;
    while (mesh.numCells() > CELLNUM_THRESHOLD && level < max_levels) {
        // create new coarsened level
        auto coarsened = coarsen_mesh(mesh, level == 0 ? fixed_cells : std::vector<CellRef>());
        auto& new_mesh = coarsened.first;
        auto& uncoarsened_cells = coarsened.second;

        // add uncoarsened cells into result
        for (const auto& cell : uncoarsened_cells) {
            for (const auto& tri : tesselate_coarsecell(cell, mesh, level == 0)) {
                result_triangles.push_back({coarse_to_fine(tri[0], level - 1),
                                            coarse_to_fine(tri[1], level - 1),
                                            coarse_to_fine(tri[2], level - 1)});
                result_cells.push_back(coarse_to_fine(cell, level));
            }
        }

        // increment level and swap meshes
        ++level;
        mesh.swap(new_mesh);
    }

    // add remaining cells into result
    for (const auto& cell : mesh.cellIndices()) {
        for (const auto& tri : tesselate_coarsecell(cell, mesh, level == 0)) {
            result_triangles.push_back({coarse_to_fine(tri[0], level - 1),
                                        coarse_to_fine(tri[1], level - 1),
                                        coarse_to_fine(tri[2], level - 1)});

            result_cells.push_back(coarse_to_fine(cell, level));
        }
    }

    return std::make_pair(result_triangles, result_cells);
}

//------------------------------------------------------------------------------
std::pair<RegularTrimesh, std::vector<CellRef>>
RegularTrimesh::coarsen_mesh(const RegularTrimesh& mesh, const std::vector<CellRef>& fixed_cells)
//------------------------------------------------------------------------------
{
    RegularTrimesh tmp_mesh(mesh);

    for (const auto& cell : fixed_cells) {
        tmp_mesh.setInactive(cell);
    }

    tmp_mesh.contractGrid(); // exclude boundary cells from coarsening

    const RegularTrimesh coarsened = remove_mesh_corners(tmp_mesh.coarsen(true));

    // identify which fine-scale cells have been covered by the coarse grid
    std::vector<CellRef> covered_finecells; // will be gradually filled below
    for (const auto& cell : coarsened.cellIndices()) {
        for (const auto& cref : coarse_to_fine(cell)) {
            covered_finecells.push_back(cref);
        }
    }

    std::vector<CellRef> all_finecells = mesh.cellIndices(); // these are sorted
    std::vector<CellRef> uncovered_finecells;

    std::sort(covered_finecells.begin(),
              covered_finecells.end()); // required by set_difference

    std::set_difference(all_finecells.begin(),
                        all_finecells.end(),
                        covered_finecells.begin(),
                        covered_finecells.end(),
                        std::back_inserter(uncovered_finecells));

    return std::make_pair(coarsened, uncovered_finecells);
}

//------------------------------------------------------------------------------
std::vector<std::tuple<unsigned int, unsigned int, double>>
RegularTrimesh::createGridToGridMap(const std::vector<std::vector<CellRef>>& map1,
                                    const std::vector<std::vector<CellRef>>& map2,
                                    const int level)
//------------------------------------------------------------------------------
{
    std::vector<std::tuple<unsigned int, unsigned int, double>> result;

    // create the inverse (multi)map of map2
    std::map<CellRef, std::vector<unsigned int>> map2_inv;
    for (std::size_t i = 0; i != map2.size(); ++i) {
        for (const auto& cell : map2[i]) {
            map2_inv[cell].push_back(i);
        }
    }

    for (unsigned int i = 0; i != map1.size(); ++i) {
        const auto& mv1 = map1[i];
        const double w1 = 1.0 / mv1.size(); // weight

        for (const auto& cell : mv1) {
            // create a CellRef for the coarse cell
            const CellRef cref = fine_to_coarse(cell, level);
            const auto& target_cells = map2_inv[cref];
            const double w2 = 1.0 / target_cells.size();

            for (const auto j : target_cells) {
                result.emplace_back(i, j, w1 * w2);
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------
std::vector<unsigned int>
RegularTrimesh::noderefs_to_indices_(const std::vector<NodeRef>& noderefs) const
//------------------------------------------------------------------------------
{
    std::vector<unsigned int> result;

    const auto nix = nodeIndices();

    std::map<NodeRef, unsigned int> nodemap;
    for (std::size_t i = 0; i != nix.size(); ++i) {
        nodemap.insert_or_assign(nix[i], i);
    }

    for (const auto& node : noderefs) {
        result.push_back(nodemap[node]);
    }

    return result;
}

//------------------------------------------------------------------------------
std::vector<unsigned int>
RegularTrimesh::cellrefs_to_indices_(const std::vector<CellRef>& cellrefs) const
//------------------------------------------------------------------------------
{
    std::vector<unsigned int> result;
    result.reserve(cellrefs.size());

    const auto cix = cellIndices();

    std::map<CellRef, unsigned int> cellmap;
    for (std::size_t i = 0; i != cix.size(); ++i) {
        cellmap.insert_or_assign(cix[i], i);
    }

    for (const auto& cell : cellrefs) {
        result.push_back(cellmap[cell]);
    }

    return result;
}

//------------------------------------------------------------------------------
std::vector<std::array<unsigned int, 3>>
RegularTrimesh::ref2index_(const std::vector<std::array<NodeRef, 3>>& vec) const
//------------------------------------------------------------------------------
{
    //  collect all NodeRefs in one vector
    std::vector<NodeRef> all_nodes;
    for (const auto& tri : vec) {
        all_nodes.insert(all_nodes.end(), tri.begin(), tri.end());
    }

    // compute linear indices
    const std::vector<unsigned int> node_indices = noderefs_to_indices_(all_nodes);

    // return vector with indices rather than NodeRef
    std::vector<std::array<unsigned int, 3>> result;
    result.reserve(vec.size());

    for (std::size_t i = 0; i != vec.size(); ++i) {
        result.push_back({node_indices[3 * i], node_indices[3 * i + 1], node_indices[3 * i + 2]});
    }

    return result;
}

//------------------------------------------------------------------------------
std::tuple<std::unique_ptr<Grid>, std::vector<std::vector<CellRef>>, std::map<CellRef, int>>
RegularTrimesh::createDuneGrid(const int coarsen_levels,
                               const std::vector<CellRef>& fixed_cells,
                               const bool add_smoothing_triangles,
                               const int cellnum_threshold) const
//------------------------------------------------------------------------------
{
    Dune::GridFactory<Grid> factory;

    // define points
    for (const auto& node : nodeCoords()) {
        factory.insertVertex(Dune::FieldVector<double, 3> {node[0], node[1], node[2]});
    }

    // define triangles
    auto tmp = coarsen_levels > 0 ? getMultiresTriangles(fixed_cells, coarsen_levels, cellnum_threshold)
                                  : getTriangles();

    const auto& triangles = ref2index_(tmp.first); // std::vector<std::array<uint, 3>>

    auto& fsmap = tmp.second;

    for (const auto& tri : triangles) {
        factory.insertElement(Dune::GeometryTypes::simplex(2),
                              std::vector<unsigned int> {tri[0], tri[1], tri[2]});
    }

    // create map from trimesh boundary cells to Dune grid cells
    std::map<CellRef, int> boundary_map {};
    std::map<CellRef, int> fsmap_inv {};
    for (auto i = 0 * fsmap.size(); i != fsmap.size(); ++i) {
        auto& map_i = fsmap[i];

        if (std::size(map_i) == 1) { // there's a one-to-one mapping -> already a fine-scale cell
            fsmap_inv.insert_or_assign(map_i.front(),
                                       static_cast<int>(i)); // mapping from Trimesh to Dune grid cells
        }
    }

    for (const auto& cell : boundaryCells()) {
        boundary_map[cell] = fsmap_inv[cell]; // mapping from Trimesh cells to Dune grid cells
    }

    if (add_smoothing_triangles) {
        const auto smoothing_triangles = boundary_smoothing_triangles_();
        for (const auto& tri : smoothing_triangles) {
            factory.insertElement(Dune::GeometryTypes::simplex(2),
                                  std::vector<unsigned int> {tri.first[0], tri.first[1], tri.first[2]});

            // make the associated boundary cells in the TriMesh map to the
            // smoothing triangle in the Dune grid
            boundary_map[tri.second[0]] = fsmap.size();
            boundary_map[tri.second[1]] = fsmap.size();

            fsmap.push_back(std::vector<CellRef> {tri.second[0], tri.second[1]});
        }
    }

    return std::make_tuple(factory.createGrid(), fsmap, boundary_map);
}

// ----------------------------------------------------------------------------
std::pair<NodeRef, NodeRef>
RegularTrimesh::edgeNodes(const EdgeRef& e) const
// ----------------------------------------------------------------------------
{
    return e[2] == 0 ? std::make_pair(NodeRef {e[0], e[1]}, NodeRef {e[0] + 1, e[1]})
        : e[2] == 1  ? std::make_pair(NodeRef {e[0], e[1]}, NodeRef {e[0], e[1] + 1})
                     : std::make_pair(NodeRef {e[0], e[1] + 1}, NodeRef {e[0] + 1, e[1]});
}

// ----------------------------------------------------------------------------
std::vector<std::pair<std::size_t, std::size_t>>
RegularTrimesh::edgeNodeIndices(bool only_boundary) const
// ----------------------------------------------------------------------------
{
    const auto nodeindices = nodeIndices();
    const auto edgeindices = only_boundary ? boundaryEdges() : edgeIndices();

    // make mapping from node to index
    std::map<NodeRef, std::size_t> node2index;
    for (std::size_t i = 0; i != nodeindices.size(); ++i) {
        node2index.insert_or_assign(nodeindices[i], i);
    }

    // map back to indices, for all edges in the mesh
    std::vector<std::pair<std::size_t, std::size_t>> result;
    for (const auto& edge : edgeindices) {
        const auto nodes = edgeNodes(edge);
        result.emplace_back(node2index[nodes.first], node2index[nodes.second]);
    }

    return result;
}

// ----------------------------------------------------------------------------
std::pair<Coord3D, Coord3D>
RegularTrimesh::edgeNodeCoords(const EdgeRef& edge) const
// ----------------------------------------------------------------------------
{
    return std::make_pair(nodeCoord(edgeNodes(edge).first), nodeCoord(edgeNodes(edge).second));
}

// ----------------------------------------------------------------------------
std::vector<std::pair<Coord3D, Coord3D>>
RegularTrimesh::edgeNodeCoords() const
// ----------------------------------------------------------------------------
{
    std::vector<std::pair<Coord3D, Coord3D>> result;

    const auto edgeindices = edgeIndices();

    result.reserve(edgeindices.size());

    for (const auto& edge : edgeindices) {
        result.push_back(edgeNodeCoords(edge));
    }

    return result;
}

// ----------------------------------------------------------------------------
std::array<NodeRef, 3>
RegularTrimesh::cellNodes(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    return cell[2] == 0 ? std::array<NodeRef, 3> {NodeRef {cell[0], cell[1]},
                                                  NodeRef {cell[0] + 1, cell[1]},
                                                  NodeRef {cell[0], cell[1] + 1}}
                        : std::array<NodeRef, 3> {NodeRef {cell[0] + 1, cell[1]},
                                                  NodeRef {cell[0], cell[1] + 1},
                                                  NodeRef {cell[0] + 1, cell[1] + 1}};
}

// ----------------------------------------------------------------------------
std::vector<std::array<std::size_t, 3>>
RegularTrimesh::cellNodesLinear() const
// ----------------------------------------------------------------------------
{
    std::vector<std::array<std::size_t, 3>> result;

    const auto nodeindices = nodeIndices();
    const auto cellindices = cellIndices();

    const auto findnode = [&nodeindices](const NodeRef& node) {
        auto nodePos = std::find(nodeindices.begin(), nodeindices.end(), node);

        return static_cast<std::size_t>(std::distance(nodeindices.begin(), nodePos));
    };

    for (const auto& cell : cellindices) {
        const auto noderefs = cellNodes(cell);
        result.push_back({findnode(noderefs[0]), findnode(noderefs[1]), findnode(noderefs[2])});
    }

    return result;
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::writeMatlabTriangulation(std::ostream& out) const
// ----------------------------------------------------------------------------
{
    const auto nodeindices = nodeIndices();
    const auto nodecoords = nodeCoords();

    // define the vertices of the triangulation
    out << "vertices = [";
    for (const auto& node : nodecoords) {
        out << node[0] << " " << node[1] << " " << node[2] << "; ";
    }
    out << "];\n";

    // define the triangles of the triangulation
    const auto cellnodes = cellNodesLinear();
    out << "triangles = [";
    for (const auto& cell : cellnodes) {
        out << cell[0] + 1 << " " << cell[1] + 1 << " " << cell[2] + 1 << "; ";
    }
    out << "];\n";

    // create a triangulation object
    out << "tri = triangulation(triangles, vertices);\n";

    // plot the triangulation, using triplot
    out << "figure;\n";
    out << "triplot(tri);\n";

    // set axis equal
    out << "axis equal;\n";
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::inner_ring_cells() // static function
// ----------------------------------------------------------------------------
{
    // utility function to quickly construct a vector with reference to all
    // cells having the origin as a corner
    return std::vector<CellRef> {{0, 0, 0}, {-1, 0, 1}, {-1, 0, 0}, {-1, -1, 1}, {0, -1, 0}, {0, -1, 1}};
}

// ----------------------------------------------------------------------------
void
writeMeshToVTK(const RegularTrimesh& mesh,
               const char* const filename,
               const int coarsen_levels,
               const std::vector<CellRef>& fixed_cells,
               const bool add_smoothing_triangles)
// ----------------------------------------------------------------------------
{
    const auto grid = mesh.createDuneGrid(coarsen_levels, fixed_cells, add_smoothing_triangles);

    const auto& g = std::get<0>(grid);
    const auto& mesh_fsmap = std::get<1>(grid);
    const auto& boundary_map = std::get<2>(grid);

    // write grid to file
    auto vtkwriter = Dune::VTKWriter {std::get<0>(grid)->leafGridView(), Dune::VTK::nonconforming};
    // write flag to file
    if (coarsen_levels == 0) {
        std::vector<int> flags = mesh.getCellFlags();

        if (add_smoothing_triangles) {
            flags.resize(std::get<0>(grid)->size(0), -1);
        }

        vtkwriter.addCellData(flags, "flag");
        vtkwriter.write(filename);
    } else {
        vtkwriter.write(filename);
    }
}

// ----------------------------------------------------------------------------
void
writeMeshToVTKDebug(const RegularTrimesh& mesh,
                    const char* const filename,
                    const int coarsen_levels,
                    const bool add_smoothing_triangles)
// ----------------------------------------------------------------------------
{
    writeMeshToVTK(mesh, filename, coarsen_levels, std::vector<CellRef>(), add_smoothing_triangles);
}

// ----------------------------------------------------------------------------
void
writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const char* const filename)
// ----------------------------------------------------------------------------
{
    std::ofstream file(filename);

    const std::vector<Coord3D> points = mesh.nodeCoords();
    const std::vector<std::pair<std::size_t, std::size_t>> edges = mesh.edgeNodeIndices(true);

    // header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << edges.size()
         << "\">\n";
    // points
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";
    for (const auto& point : points) {
        file << "          " << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // cells (edges)
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& edge : edges) {
        file << "          " << edge.first << " " << edge.second << "\n";
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= edges.size(); ++i) {
        file << "          " << i * 2 << "\n";
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < edges.size(); ++i) {
        file << "          3\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // Footer
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
}

// ----------------------------------------------------------------------------
std::size_t
RegularTrimesh::numActive() const
// ----------------------------------------------------------------------------
{
    return cellinfo_.size();
}

// ----------------------------------------------------------------------------
bool
RegularTrimesh::isActive(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
    return cellinfo_.find(cell) != cellinfo_.end();
}

// ----------------------------------------------------------------------------
bool
RegularTrimesh::setInactive(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    if (!isActive(cell)) {
        return false;
    }

    // remove this cell from cellinfo
    cellinfo_.erase(cell);

    return true;
}

// ----------------------------------------------------------------------------
bool
RegularTrimesh::setActive(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    if (isActive(cell)) {
        return false;
    }

    cellinfo_.insert({cell, CellAttributes()});

    assert(cellinfo_[cell].flag == 0);

    return true;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::expandGrid(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    int result = 0;

    const auto edges = cell2edges(cell);
    for (const auto& e : edges) {
        for (const auto& c : edge2cells(e)) {
            result += setActive(c);
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::contractGrid()
// ----------------------------------------------------------------------------
{
    // identify boundary cells and remove them
    const std::vector<CellRef> bcells = boundaryCells();

    // identify all cells that are not boundary cells, using set_difference
    for (const auto& bc : bcells) {
        cellinfo_.erase(bc);
    }

    return static_cast<int>(bcells.size());
}


// ----------------------------------------------------------------------------
int
RegularTrimesh::expandGrid(const std::vector<CellRef>& cells)
// ----------------------------------------------------------------------------
{
    int result = 0;

    for (const auto& c : cells) {
        result += expandGrid(c);
    }

    return result;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::expandGrid()
// ----------------------------------------------------------------------------
{
    // regular expansion in all directions
    return expandGrid(boundaryCells());
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::removeSawtooths()
// ----------------------------------------------------------------------------
{
    std::vector<CellRef> candidates;

    const auto bedges = boundaryEdges();
    for (const auto& edge : bedges) {
        for (const auto& cell : edge2cells(edge)) {
            if (!isActive(cell)) {
                candidates.push_back(cell);
            }
        }
    }

    std::sort(candidates.begin(), candidates.end());

    // inactive cells adjacent to more than one boundary edge should be activated
    for (auto it = candidates.begin(); it != candidates.end(); ++it) {
        const auto range = std::equal_range(it, candidates.end(), *it);
        if (range.second - range.first > 1) {
            setActive(*it);
            it = range.second - 1;
        }
    }
}

// ----------------------------------------------------------------------------
std::size_t
RegularTrimesh::linearCellIndex(const CellRef& cell) const
{
    // cellinfo is a map from CellRef to CellAttributes.
    assert(cellinfo_.find(cell) != cellinfo_.end());
    return std::distance(cellinfo_.begin(), cellinfo_.find(cell));
}

// ----------------------------------------------------------------------------
CellRef
RegularTrimesh::cellIndex(const std::size_t index) const
{
    auto it = cellinfo_.begin();

    std::advance(it, index);

    return it->first;
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::setCellFlag(const CellRef& cell, const int value)
// ----------------------------------------------------------------------------
{
    cellinfo_[cell].flag = value;
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::setCellFlags(const std::vector<CellRef>& cells, const int value)
// ----------------------------------------------------------------------------
{
    for (const auto& cell : cells) {
        setCellFlag(cell, value);
    }
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::setAllFlags(const int value)
// ----------------------------------------------------------------------------
{
    for (auto& it : cellinfo_) {
        it.second.flag = value;
    }
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::getCellFlag(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
    const auto it = cellinfo_.find(cell);
    return it->second.flag;
}

// ----------------------------------------------------------------------------
std::vector<int>
RegularTrimesh::getCellFlags() const
// ----------------------------------------------------------------------------
{
    std::vector<int> result;

    for (const auto& el : cellinfo_) {
        result.push_back(el.second.flag);
    }

    return result;
}

// ----------------------------------------------------------------------------
RegularTrimesh
RegularTrimesh::refine() const
// ----------------------------------------------------------------------------
{
    std::map<CellRef, CellAttributes> new_cells;

    for (const auto& e : cellinfo_) {
        for (const auto& c : coarse_to_fine(e.first)) {
            new_cells[c] = e.second;
        }
    }

    return RegularTrimesh {new_cells, origin_, axis1_, axis2_, {edgelen_[0] / 2, edgelen_[1] / 2}};
}

// ----------------------------------------------------------------------------
RegularTrimesh
RegularTrimesh::coarsen(bool strict) const
// ----------------------------------------------------------------------------
{
    std::map<CellRef, CellAttributes> new_cells;

    for (const auto& e : cellinfo_) {
        const CellRef& c = e.first;
        if (!strict || !is_boundary_cell(c, *this)) {
            if ((std::abs(c[0]) % 2 == 1 && std::abs(c[1]) % 2 == 1 && c[2] == 0)
                || (std::abs(c[0]) % 2 == 0 && std::abs(c[1]) % 2 == 0 && c[2] == 1)) {
                new_cells[fine_to_coarse(c)] = e.second;
            }
        }
    }

    return RegularTrimesh {new_cells, origin_, axis1_, axis2_, {edgelen_[0] * 2, edgelen_[1] * 2}};
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::coarse_to_fine(const CellRef& c, const int level)
// ----------------------------------------------------------------------------
{
    assert(level >= 0); // can't have negative refinement

    if (level == 0) {
        return std::vector<CellRef> {c}; // no refinement, return the cell itself
    }

    const std::vector<CellRef> fine_cells = (c[2] == 0)
        ? std::vector<CellRef> {CellRef {2 * c[0], 2 * c[1], 0},
                                CellRef {2 * c[0] + 1, 2 * c[1], 0},
                                CellRef {2 * c[0], 2 * c[1] + 1, 0},
                                CellRef {2 * c[0], 2 * c[1], 1}}
        : std::vector<CellRef> {CellRef {2 * c[0] + 1, 2 * c[1] + 1, 0},
                                CellRef {2 * c[0] + 1, 2 * c[1], 1},
                                CellRef {2 * c[0], 2 * c[1] + 1, 1},
                                CellRef {2 * c[0] + 1, 2 * c[1] + 1, 1}};

    if (level == 1) {
        return fine_cells;
    }

    // call recursively if level is higher than 1
    std::vector<CellRef> result;
    for (const auto& cell : fine_cells) {
        const std::vector<CellRef> finer_cells = coarse_to_fine(cell, level - 1);
        result.insert(result.end(), finer_cells.begin(), finer_cells.end());
    }

    return result;
}

// ----------------------------------------------------------------------------
CellRef
RegularTrimesh::fine_to_coarse(const CellRef& cell, const int levels)
// ----------------------------------------------------------------------------
{
    int i = cell[0], j = cell[1], k = cell[2];

    for (int l = 0; l < levels; ++l) {
        const int impairs = std::abs(i % 2) + std::abs(j % 2) + k;

        if (i < 0) {
            --i;
        }
        if (j < 0) {
            --j;
        }

        i /= 2;
        j /= 2;

        k = (impairs > 1) ? 1 : 0;
    }

    return {i, j, k};
}

// ----------------------------------------------------------------------------
NodeRef
RegularTrimesh::coarse_to_fine(const NodeRef& node, const int levels)
// ----------------------------------------------------------------------------
{
    if (levels < 0) {
        return {node[0] / (1 << abs(levels)), node[1] / (1 << abs(levels))};
    } else {
        return {node[0] * (1 << levels), node[1] * (1 << levels)};
    }
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::interior_coarsegrid_() const
// ----------------------------------------------------------------------------
{
    std::map<CellRef, int> cell_count;
    for (const auto& e : cellinfo_) {
        const auto& [countPos, inserted] = cell_count.try_emplace(e.first, 1);

        if (!inserted) {
            ++countPos->second;
        }
    }

    std::vector<CellRef> result;
    for (const auto& e : cell_count) {
        if (e.second == 4) { // this cell is fully covered by all its fine cells
            result.push_back(e.first);
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<CellRef>
RegularTrimesh::activeNeighborCells(const std::vector<CellRef>& cells) const
// ----------------------------------------------------------------------------
{
    std::set<CellRef> result;
    for (const auto& cell : cells) {
        for (const auto& edge : cell2edges(cell)) {
            for (const auto& c : edge2cells(edge)) {
                if (isActive(c)) {
                    result.insert(c);
                }
            }
        }
    }

    for (const auto& c : cells) {
        result.erase(c); // remove the original cells
    }

    return {result.begin(), result.end()};
}

// ----------------------------------------------------------------------------
std::tuple<RegularTrimesh, int>
expand_to_criterion(
    const RegularTrimesh& mesh,
    std::function<std::vector<double>(const RegularTrimesh&, const int level)> score_function,
    const double threshold,
    const std::vector<CellRef>& fixed_cells,
    const int target_cellcount,
    const int cellcount_threshold)
{
    RegularTrimesh working_mesh = mesh; // make a working copy of the mesh;
    std::vector<RegularTrimesh> last_meshes; // keep track of meshes at each level before coarsening
    int cur_level = 0; // start expanding mesh at finest level
    int iter_count = 0; // keep track of iterations on current level
    const int max_iter = 5; // 5; // maximum number of iterations on current level
    int roof = std::numeric_limits<int>::max();

    ++DEBUG_GRID_COUNT; // @@ keeping track of grids to output for debugging/monitoring
                        // purposes
    DEBUG_CURRENT_GRID_ITERATION_COUNT = 0; //@@ same

    // determine starting level
    // const int target_cellcount = 50; // target number of cells in the final mesh
    // const int cellcount_threshold = 4*target_cellcount; // target number of cells in
    // the initial mesh
    const int max_cellcount = 2000; // maximum number of cells in the final mesh

    auto fixed_on_level = [&fixed_cells](const int level) -> std::vector<CellRef> {
        if (level == 0) {
            for (const auto& cell : fixed_cells) {
                std::cout << "{" << cell[0] << ", " << cell[1] << ", " << cell[2] << "} ";
            }

            std::cout << std::endl;
            return fixed_cells;
        } else {
            std::vector<CellRef> result;
            for (const auto& cell : fixed_cells) {
                result.push_back(RegularTrimesh::fine_to_coarse(cell, level));
            }

            // dump vector result to std::cout
            for (const auto& cell : result) {
                std::cout << "{" << cell[0] << ", " << cell[1] << ", " << cell[2] << "} ";
            }

            std::cout << std::endl;

            return result;
        }
    };

    while (working_mesh.numCells() > cellcount_threshold) {
        last_meshes.push_back(working_mesh);
        working_mesh = working_mesh.coarsen(true);
        working_mesh.setCellFlags(fixed_on_level(cur_level),
                                  1); // set fixed cells at this level
        working_mesh.removeSawtooths();
        ++cur_level;
    }

    std::cout << "---------- Starting propagation at level: " << cur_level << " --------" << std::endl;
    while (true) { // keep looping as long as grid need expansion
        if (DEBUG_DUMP_GRIDS) {
            const std::string filename = "current_grid_" + std::to_string(DEBUG_GRID_COUNT) + "_"
                + std::to_string(DEBUG_CURRENT_GRID_ITERATION_COUNT++);

            writeMeshToVTKDebug(working_mesh, filename.c_str(), 0, 1);
        }

        const std::vector<double> bnd_scores = score_function(working_mesh, cur_level);
        const std::vector<CellRef> bnd_cells = mesh.boundaryCells();
        assert(bnd_scores.size() == bnd_cells.size());

        std::vector<CellRef> expand_cells;

        for (std::size_t i = 0; i != bnd_scores.size(); ++i) {
            if (bnd_scores[i] > threshold) {
                expand_cells.push_back(bnd_cells[i]);
            }
        }

        if (expand_cells.size() == 0) {
            if (cur_level == 0 || working_mesh.numCells() >= target_cellcount) {
                break;
            }

            working_mesh.contractGrid();
            working_mesh = working_mesh.refine();

            // ensure we did not lose cells that were already inherited from
            // finer level when the coarse mesh was created
            const auto prev_cells = last_meshes.back().cellIndices();
            for (const auto& cell : prev_cells) {
                working_mesh.setActive(cell);
            }

            last_meshes.pop_back();
            working_mesh.removeSawtooths();

            roof = cur_level--;
            std::cout << "** -------- Refining to level -------- " << cur_level << std::endl;
            iter_count = 0;
        } else if (iter_count >= max_iter && cur_level < roof - 1) {
            // expansion is going too slowly, move to coarser level
            last_meshes.push_back(working_mesh);

            working_mesh = working_mesh.coarsen(true);
            working_mesh.setCellFlags(fixed_on_level(cur_level),
                                      1); // set fixed cells at this level
            working_mesh.removeSawtooths();
            ++cur_level;

            std::cout << "** -------- Coarsening to level ------- " << cur_level << std::endl;
            iter_count = 0;
        } else {
            // expanding grid at current level
            working_mesh.expandGrid(expand_cells);
            working_mesh.removeSawtooths();
            iter_count++;
        }
    }

    std::cout << " ** ---------- CONVERGED BOUNDARY MESH ---------- **" << std::endl;

    return {working_mesh, cur_level};
}

} // namespace Opm
