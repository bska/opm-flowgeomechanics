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

#include <opm/geomech/GridStretcher.hpp>

#include <dune/common/fvector.hh> // FieldVector

#include <dune/grid/common/mcmgmapper.hh> // for element mapper
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <opm/geomech/convex_boundary.hpp>
#include <opm/geomech/param_interior.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <vector>

// ============================================================================
namespace // anonymous
// ============================================================================
{

using GridView = Opm::GridStretcher::Grid::LeafGridView;
using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
using CoordType = Opm::GridStretcher::CoordType;

// ----------------------------------------------------------------------------
double
adjust_disp_to_convex(const double* const startpt,
                      const double* const endpt,
                      const double* const pt,
                      const double* const dir)
// ----------------------------------------------------------------------------
{
    // solve for alpha in the equation:
    // alpha * dir + t * (endpt-startpt) = enpt - pt
    // This will ensure that pt + dir * alpha lies on the edge between startpt and
    // endpt.
    const double vx = endpt[0] - startpt[0];
    const double vy = endpt[1] - startpt[1];
    const double dx = endpt[0] - pt[0];
    const double dy = endpt[1] - pt[1];

    return (vy * dx - vx * dy) / (dir[0] * vy - dir[1] * vx);
}

// ----------------------------------------------------------------------------
template <typename T>
std::vector<T>
extract_elements(const std::vector<T>& vec, const std::vector<std::size_t>& ixs)
// ----------------------------------------------------------------------------
{
    auto result = std::vector<T> {};
    result.reserve(ixs.size());

    for (const auto& ix : ixs) {
        result.push_back(vec[ix]);
    }

    return result;
}

// ----------------------------------------------------------------------------
double
dist3D(const Dune::FieldVector<double, 3>& p1, const Dune::FieldVector<double, 3>& p2)
// ----------------------------------------------------------------------------
{
    return std::hypot(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
}

// ----------------------------------------------------------------------------
std::vector<double>
pick_nodes(const std::vector<std::size_t>& ixs, const std::vector<double>& coords2D)
// ----------------------------------------------------------------------------
{
    std::vector<double> result;

    for (const auto& i : ixs) {
        result.push_back(coords2D[2 * i]);
        result.push_back(coords2D[2 * i + 1]);
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<double>
pickCoords2D(const std::vector<CoordType>& coords3D)
// ----------------------------------------------------------------------------
{
    // pick the two coordinate axes with the largest spread
    const std::size_t N = coords3D.size(); // total number of nodes (interior and boundary)

    auto low = std::array {coords3D[0][0], coords3D[0][1], coords3D[0][2]};
    auto high = low;

    for (std::size_t i = 0; i != N; ++i) {
        for (int d = 0; d != 3; ++d) {
            low[d] = std::min(low[d], coords3D[i][d]);
            high[d] = std::max(high[d], coords3D[i][d]);
        }
    }

    const std::array<double, 3> span {high[0] - low[0], high[1] - low[1], high[2] - low[2]};

    const std::size_t min_ix = std::distance(span.begin(), std::min_element(span.begin(), span.end()));
    const std::size_t ix1 = (min_ix + 1) % 3;
    const std::size_t ix2 = (min_ix + 2) % 3;

    std::vector<double> coords2D;
    for (std::size_t i = 0; i != N; ++i) {
        coords2D.push_back(coords3D[i][ix1]);
        coords2D.push_back(coords3D[i][ix2]);
    }

    return coords2D;
}

// ----------------------------------------------------------------------------
std::size_t
find_coord_in(const CoordType& c, const std::vector<CoordType>& cvec)
// ----------------------------------------------------------------------------
{
    const double TOL = 1e-3;

    return std::find_if(
               cvec.begin(), cvec.end(), [&](const CoordType& el) { return dist3D(el, c) < TOL; })
        - cvec.begin();
}

// ----------------------------------------------------------------------------
CoordType
normalize(const CoordType& vec)
{
    return vec / vec.two_norm();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
CoordType
cross(const CoordType& v1, const CoordType& v2)
// ----------------------------------------------------------------------------
{
    return CoordType {
        v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]};
}

// ----------------------------------------------------------------------------
std::vector<CoordType>
node_coordinates(const Opm::GridStretcher::Grid& g)
// ----------------------------------------------------------------------------
{
    std::vector<CoordType> result;

    for (const auto& v : vertices(g.leafGridView())) {
        result.push_back(v.geometry().corner(0));
    }

    return result;
}

std::vector<std::size_t>
boundary_node_indices(const Opm::GridStretcher::Grid& grid)
{
    std::vector<std::size_t> bix;

    {
        auto gv = grid.leafGridView();
        for (const auto& el : Dune::elements(gv)) {
            if (!el.hasBoundaryIntersections()) {
                continue;
            }

            const auto& refEl = Dune::referenceElement(el);

            for (const auto& is : Dune::intersections(gv, el)) {
                if (!is.boundary()) {
                    continue;
                }

                const auto inside = is.indexInInside();
                const auto faceSize = refEl.size(inside, /*face*/ 1, /*node*/ 2);

                for (int i = 0; i < faceSize; ++i) {
                    // this is a 2 grid where faces=cells=codim 1 nodes is of codim 2
                    const auto corner = refEl.subEntity(
                        /*facenr*/ inside, /*face*/ 1, /*nodenum*/ i, /*node*/ 2);

                    const auto cornerIndex = gv.indexSet().subIndex(el, corner, 2);

                    bix.push_back(cornerIndex);
                }
            }
        }
    }

    // make unique
    std::sort(bix.begin(), bix.end());
    bix.erase(std::unique(bix.begin(), bix.end()), bix.end());

    // ensure correct order (we use convex hull algorithm for that)
    const auto ncoords = node_coordinates(grid);
    const auto plane = Opm::best2Dplane(ncoords);

    std::vector<CoordType> bcoords;
    bcoords.reserve(bix.size());
    for (const auto& i : bix) {
        bcoords.push_back(ncoords[i]);
    }

    const auto bpts2D = Opm::projectPointsTo2DPlane(bcoords, plane);

    auto bix_ordered = Opm::convex_hull(bpts2D);
    assert(bix_ordered.size() == bix.size()); // should be the case if boundary is convex

    for (auto& i : bix_ordered) {
        i = bix[i];
    }

    return bix_ordered;
}

// ----------------------------------------------------------------------------
std::vector<std::size_t>
complement_of(const std::vector<std::size_t>& vec, const std::size_t N)
// ----------------------------------------------------------------------------
{
    // indices not mentioned in 'vec' are stored in 'result'
    std::vector<std::size_t> result {};

    // indicator vector to identify the indices mentioned in 'vec'
    std::vector<int> tmp(N, 0);
    for (const auto& v : vec) {
        tmp[v] = 1;
    }

    for (std::size_t ix = 0; ix != N; ++ix) {
        if (tmp[ix] == 0) {
            result.push_back(ix);
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<double>
interior_parametrization(const std::vector<std::size_t>& bindices,
                         const std::vector<std::size_t>& iindices,
                         const std::vector<CoordType>& coords3D)
// ----------------------------------------------------------------------------
{
    // pick two out of the three coordinate axes to use for computing the 2D
    // parametrization
    const std::vector<double> coords2D(pickCoords2D(coords3D));

    // sort node coordinates into boundary and interior
    const std::vector<double> bcoords(pick_nodes(bindices, coords2D));
    const std::vector<double> icoords(pick_nodes(iindices, coords2D));

    // call parametrization function
    std::vector<double> result;
    Opm::parametrize_interior_2D(bcoords, icoords, result);

    return result;
}

// ----------------------------------------------------------------------------
Opm::GridStretcher::CellBnodeMap
compute_cell_2_bindices_mapping(const Opm::GridStretcher::Grid& grid,
                                const std::vector<std::size_t>& bindices)
// ----------------------------------------------------------------------------
{
    auto result = Opm::GridStretcher::CellBnodeMap {};

    const auto gview = grid.leafGridView();

    // coordinates of boundary indices @@ We really ought to be able to do this
    // without resorting to geometric calculations!
    // register all node coordinates
    const auto ncoords = node_coordinates(grid);
    const auto bcoords = extract_elements(ncoords, bindices);

    const ElementMapper mapper(gview, Dune::mcmgElementLayout());
    for (const auto& elem : elements(gview)) {
        for (const auto& is : intersections(gview, elem)) {
            if (is.boundary()) {
                // this is a boundary cell
                result[mapper.index(elem)] = {elem.seed(),
                                              find_coord_in(is.geometry().corner(0), bcoords),
                                              find_coord_in(is.geometry().corner(1), bcoords)};
            }
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<double>
bcentroid_param_mat(const Opm::GridStretcher::Grid& grid,
                    const std::vector<std::size_t>& bnindices,
                    const std::vector<std::size_t>& iindices,
                    const std::vector<double>& iparam,
                    const Opm::GridStretcher::CellBnodeMap& c2bix)
// ----------------------------------------------------------------------------
{
    const auto view = grid.leafGridView();

    // helper function to identify the internal node and two boundary nodes of a
    // given boundary triangle
    auto corner_nodes = [&](const Opm::GridStretcher::CellBnodeMap::value_type& pair) {
        const std::size_t bn1 = bnindices[std::get<1>(pair.second)]; // boundary node 1
        const std::size_t bn2 = bnindices[std::get<2>(pair.second)]; // boundary node 2
        const auto elem = grid.entity(std::get<0>(pair.second));

        assert(elem.subEntities(2) == 3); // should have three corners

        std::size_t in(0);
        for (int i = 0; i != 3; ++i) {
            in = view.indexSet().index(elem.subEntity<2>(i));

            if (in != bn1 && in != bn2) {
                break;
            }
        }

        return std::array<std::size_t, 3> {in, bn1, bn2};
    };

    // computing the bcentroid parameter matrix
    std::vector<double> result;
    const std::size_t Nb = bnindices.size();

    for (const auto& pair : c2bix) {
        const auto cnodes = corner_nodes(pair);
        const std::size_t iix
            = std::find(iindices.begin(), iindices.end(), cnodes[0]) - iindices.begin();

        // write in one-third of the parameterization of the internal point
        std::transform(&iparam[iix * Nb],
                       &iparam[(iix + 1) * Nb],
                       std::back_inserter(result),
                       [](const double p) { return p / 3; });

        // add one-third of each of the two boundary points
        for (int i = 1; i != 3; ++i) {
            result[result.size() - Nb + cnodes[i]] += 1.0 / 3.0;
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
Opm::BoundaryNormals
boundary_normals(const Opm::GridStretcher::Grid& grid,
                 const Opm::GridStretcher::CellBnodeMap& c2bix,
                 const std::vector<std::size_t>& bcindices,
                 const std::vector<std::size_t>& bnindices,
                 const std::vector<CoordType>& nodecoords)
// ----------------------------------------------------------------------------
{
    const std::size_t N = size(bcindices);
    // map amounts to displacement vectors for each boundary point

    auto result = Opm::BoundaryNormals {
        std::vector<CoordType>(N, {0.0, 0.0, 0.0}), // nodenormals
        std::vector<CoordType>(N, {0.0, 0.0, 0.0}), // cell "normals"
        std::vector<CoordType>(N, {0.0, 0.0, 0.0}) // edge normals
    };

    int pos = 0;
    for (const auto& bix : bcindices) {
        // compute directional vector from cell centroid to edge centroid
        const auto entry = c2bix.find(bix)->second;
        const auto elem = grid.entity(std::get<0>(entry));
        const auto ccenter = elem.geometry().center();
        const auto ecenter
            = (nodecoords[bnindices[std::get<1>(entry)]] + nodecoords[bnindices[std::get<2>(entry)]])
            / 2;

        const auto dvec = normalize(ecenter - ccenter); // cell centroid to edge centroid

        result.bcell_normals[pos] = dvec;

        // compute outward edge normal
        const auto tangent
            = nodecoords[bnindices[std::get<1>(entry)]] - nodecoords[bnindices[std::get<2>(entry)]];

        const auto enormal = normalize(cross(cross(tangent, dvec), tangent));
        result.bedge_normals[pos] = enormal;

        // accumulate result for average node normal
        result.bnode_normals[std::get<1>(entry)] += enormal;
        result.bnode_normals[std::get<2>(entry)] += enormal;

        ++pos;
    }

    // each entry is the sum of two normals, so we need to normalize again
    std::transform(result.bnode_normals.begin(),
                   result.bnode_normals.end(),
                   result.bnode_normals.begin(),
                   [](const CoordType& c) { return normalize(c); });

    return result;
}

template <typename Key, typename Value>
std::vector<Key>
keyvec(const std::map<Key, Value>& map)
{
    std::vector<Key> result;
    result.reserve(map.size());

    for (const auto& kv : map) {
        result.push_back(kv.first);
    }

    return result;
}

} // end anonymous namespace

// ============================================================================
namespace Opm
// ============================================================================
{
GridStretcher::GridStretcher(Grid& grid)
    : grid_ {grid}
    , nodecoords_ {node_coordinates(grid_)}
    , bnindices_ {boundary_node_indices(grid_)}
    , iindices_ {complement_of(bnindices_, grid_.leafGridView().size(2))}
    , iparam_ {interior_parametrization(bnindices_, iindices_, nodecoords_)}
    , c2bix_ {compute_cell_2_bindices_mapping(grid_, bnindices_)}
    , bcindices_ {keyvec(c2bix_)}
    , bcentroid_param_ {bcentroid_param_mat(grid_, bnindices_, iindices_, iparam_, c2bix_)}
    , boundary_normals_ {boundary_normals(grid_, c2bix_, bcindices_, bnindices_, nodecoords_)}
{
}

// ----------------------------------------------------------------------------
double
GridStretcher::objective(const std::vector<double>& bndisp,
                         const std::vector<double>& dtarget,
                         std::vector<double>& grad,
                         const bool fixed_cell_centroids)
// ----------------------------------------------------------------------------
{
    const std::vector<CoordType>& ncoords = nodecoords();
    const std::vector<CoordType>& normals = bnodenormals();
    const std::size_t Nb = bnindices_.size(); // number of boundary nodes (and cells)

    assert(bndisp.size() == Nb);

    // compute boundary node positions
    std::vector<CoordType> bnodes, bnodes0;

    {
        auto disp_iter = bndisp.begin();

        for (const auto& bix : bnindices_) {
            bnodes0.push_back(ncoords[bix]);
            bnodes.push_back(ncoords[bix] + normals[bix] * *disp_iter++);
        }
    }

    const std::vector<CoordType>& bnodes_for_centroids = fixed_cell_centroids ? bnodes0 : bnodes;

    // compute cell and face centroid positions, and distance vector
    std::vector<CoordType> cell_centroids;
    std::vector<CoordType> face_centroids;
    std::vector<CoordType> distance;

    auto cpar_iter = bcentroid_param_.begin();
    for (const auto& pair : c2bix_) {
        const std::size_t bn1 = std::get<1>(pair.second); // boundary node 1
        const std::size_t bn2 = std::get<2>(pair.second); // boundary node 2

        // compute cell centroid as a linear combination of all boundary nodes
        CoordType cc {0, 0, 0};
        for (std::size_t i = 0; i != Nb; ++i) {
            cc += bnodes_for_centroids[i] * *cpar_iter++;
        }
        cell_centroids.push_back(cc);

        // compute face centroid
        face_centroids.push_back(0.5 * (bnodes[bn1] + bnodes[bn2]));

        // compute distance vector between face and cell centroid
        distance.push_back(face_centroids.back() - cell_centroids.back());
    }

    // compute objective value
    double objval = 0;
    for (std::size_t i = 0; i != distance.size(); ++i) {
        // o = o + (d^2 - dtarget^2) ^ 2
        objval += std::pow(distance[i].two_norm() - dtarget[i], 2);
    }

    objval /= 2;

    // compute gradient
    grad.resize(Nb);

    std::fill(grad.begin(), grad.end(), 0.0);

    auto c2bix_iter = c2bix_.begin();
    for (std::size_t i = 0; i != Nb; ++i, ++c2bix_iter) { // loop over cells
        const double dfac = (distance[i].two_norm() - dtarget[i]) / distance[i].two_norm();

        for (std::size_t j = 0; j != Nb; ++j) { // loop over boundary nodes
            // efac is zero unless this boundary node is part of the current cell
            const double efac
                = (j == std::get<1>(c2bix_iter->second) || j == std::get<2>(c2bix_iter->second)) ? 0.5
                                                                                                 : 0.0;

            // derivative of centroid position with respect to boundary node 'j'
            const double cpar
                = fixed_cell_centroids ? 0 : bcentroid_param_[i * Nb + j]; // d(d_i)/d(b_j);

            const double m = (efac - cpar);

            for (std::size_t d = 0; d != 3; ++d) {
                // loop over dimensions
                grad[j] += dfac * distance[i][d] * m * normals[j][d];
            }
        }
    }

    return objval;
}

// ----------------------------------------------------------------------------
std::vector<double>
GridStretcher::computeBoundaryNodeDisplacements(const std::vector<double>& amounts,
                                                const std::vector<CoordType>& bnodenormals) const
// ----------------------------------------------------------------------------
{
    const std::size_t N = bcindices_.size();

    const auto& bnnorm = bnodenormals.empty() ? boundary_normals_.bnode_normals : bnodenormals;

    assert(bnnorm.size() == N);

    std::vector<double> node_amounts(N, -1 * std::numeric_limits<double>::infinity());

    std::size_t cur_bcell_ix = 0;
    for (const auto& bix : bcindices_) {
        const auto entry = c2bix_.find(bix)->second;
        const std::size_t bnodeix1 = std::get<1>(entry);
        const std::size_t bnodeix2 = std::get<2>(entry);

        const auto& enorm = boundary_normals_.bedge_normals[cur_bcell_ix];
        const auto& cnorm = boundary_normals_.bcell_normals[cur_bcell_ix];
        const auto& bn1 = bnnorm[bnodeix1];
        const auto& bn2 = bnnorm[bnodeix2];

        const double l1 = amounts[cur_bcell_ix] * enorm.dot(cnorm) / enorm.dot(bn1);
        const double l2 = amounts[cur_bcell_ix] * enorm.dot(cnorm) / enorm.dot(bn2);

        node_amounts[bnodeix1] = std::max(node_amounts[bnodeix1], l1); // @@ correct for negative values?
        node_amounts[bnodeix2] = std::max(node_amounts[bnodeix2], l2);

        ++cur_bcell_ix;
    }

    return node_amounts;
}

// ----------------------------------------------------------------------------
void
GridStretcher::expandBoundaryCells(const std::vector<double>& amounts)
// ----------------------------------------------------------------------------
{
    const auto distance = computeBoundaryNodeDisplacements(amounts, boundary_normals_.bnode_normals);

    std::vector<CoordType> displacements(amounts.size());
    for (std::size_t i = 0; i != amounts.size(); ++i) {
        displacements[i] = distance[i] * boundary_normals_.bnode_normals[i];
    }

    applyBoundaryNodeDisplacements(displacements);
}

// ----------------------------------------------------------------------------
void
GridStretcher::rebalanceBoundary()
// ----------------------------------------------------------------------------
{
    // prepare vector with boundary node coordinates @@ are these always in right order?
    const auto ncoords = node_coordinates(grid_);

    std::vector<double> bnodes3D(bnindices_.size() * 3);
    for (std::size_t i = 0; i != bnindices_.size(); ++i) {
        for (std::size_t d = 0; d != 3; ++d) {
            bnodes3D[i * 3 + d] = ncoords[bnindices_[i]][d];
        }
    }

    // project coordinates to suitable 2D plane
    std::vector<double> bnodes2D; // will be filled in the call to project_to_2D
    Axis3D ax = project_to_2D(bnodes3D, bnodes2D);

    std::vector<double> bnodes2D_redist;
    redistribute_2D(bnodes2D, bnodes2D_redist);

    // project back to 3D
    lift_to_3D(bnodes2D_redist, ax, bnodes3D);

    // reformat data to vector of CoordType
    std::vector<CoordType> new_bnode_pos;
    for (std::size_t i = 0; i != bnindices_.size(); ++i) {
        new_bnode_pos.push_back(CoordType {bnodes3D[i * 3], bnodes3D[i * 3 + 1], bnodes3D[i * 3 + 2]});
    }

    // update all other information after boundary node change
    updateAfterBoundaryChange(new_bnode_pos);
}

// ----------------------------------------------------------------------------
void
GridStretcher::updateAfterBoundaryChange(const std::vector<CoordType>& new_bcoords)
// ----------------------------------------------------------------------------
{
    auto ncoords = nodecoords();

    // write in new boundary node positions
    for (std::size_t i = 0; i != bnindices_.size(); ++i) {
        ncoords[bnindices_[i]] = new_bcoords[i];
    }

    // reset all internal nodes to zero
    for (const auto& iix : iindices_) {
        ncoords[iix] = 0;
    }

    // compute new coordinates for internal nodes
    {
        auto ipiter = iparam_.begin();

        for (const auto& iix : iindices_) {
            for (const auto& bnindex : this->bnindices_) {
                ncoords[iix] += ncoords[bnindex] * (*ipiter++);
            }
        }
    }

    // write all new node positions to mesh
    {
        std::size_t vcount = 0;

        for (const auto& vertex : vertices(grid_.leafGridView())) {
            grid_.setPosition(vertex, ncoords[vcount++]);
        }
    }

    // recompute stored geometric information
    nodecoords_ = node_coordinates(grid_);
    boundary_normals_ = boundary_normals(grid_, c2bix_, bcindices_, bnindices_, nodecoords_);
}

// ----------------------------------------------------------------------------
void
GridStretcher::applyBoundaryNodeDisplacements(const std::vector<CoordType>& disp)
// ----------------------------------------------------------------------------
{
    assert(disp.size() == bnindices_.size());

    // compute new boundary node coordinates
    std::vector<CoordType> new_bnode_pos;
    for (std::size_t i = 0; i != bnindices_.size(); ++i) {
        new_bnode_pos.push_back(nodecoords()[bnindices_[i]] + disp[i]);
    }

    updateAfterBoundaryChange(new_bnode_pos);
}

// ----------------------------------------------------------------------------
std::vector<double>
GridStretcher::centroidEdgeDist() const
// ----------------------------------------------------------------------------
{
    const std::size_t N = size(bcindices_);
    std::vector<double> result(N);

    for (std::size_t bc = 0; bc != N; ++bc) {
        const auto entry = c2bix_.find(bcindices_[bc])->second;
        const auto elem = grid_.entity(std::get<0>(entry));
        const auto ccenter = elem.geometry().center();
        const auto ecenter
            = (nodecoords_[bnindices_[std::get<1>(entry)]] + nodecoords_[bnindices_[std::get<2>(entry)]])
            / 2;

        result[bc] = (ccenter - ecenter).two_norm();
    }
    return result;
}

// ----------------------------------------------------------------------------
double
GridStretcher::maxBoxLength() const
// ----------------------------------------------------------------------------
{
    assert(!nodecoords_.empty());

    std::array<double, 3> low {nodecoords_[0][0], nodecoords_[0][1], nodecoords_[0][2]};
    std::array<double, 3> high {low};

    for (const auto& n : nodecoords_) {
        for (int d = 0; d != 3; ++d) {
            low[d] = std::min(low[d], n[d]);
            high[d] = std::max(high[d], n[d]);
        }
    }

    for (int d = 0; d != 3; ++d) {
        high[d] = high[d] - low[d];
    }

    return *std::max_element(high.begin(), high.end());
}

// ----------------------------------------------------------------------------
void
GridStretcher::adjustToConvex(std::vector<double>& disp,
                              std::vector<double>& total_disp,
                              const std::vector<CoordType>& dirs) const
// ----------------------------------------------------------------------------
{
    const std::size_t N = bnindices_.size();
    assert(dirs.size() == N);
    assert(disp.size() == N);
    assert(total_disp.size() == N);

    // compute new boundary node coordinates, before convexity is enforced
    std::vector<CoordType> old_bcoords(dirs.size()), new_bcoords(dirs.size());
    for (std::size_t i = 0; i != N; ++i) {
        old_bcoords[i] = nodecoords_[bnindices_[i]];
        new_bcoords[i] = old_bcoords[i] + dirs[i] * disp[i];
    }

    // reduce geometry to 2D, to allow convex hull algorithm to work
    const std::array<int, 2> plane = best2Dplane(new_bcoords);
    const std::vector<double> bpts2D = projectPointsTo2DPlane(new_bcoords, plane);
    const std::vector<double> bpts2D_old = projectPointsTo2DPlane(old_bcoords, plane);
    const std::vector<double> dirs2D = projectPointsTo2DPlane(dirs, plane);

    // compute convex hull of boundary points
    const std::vector<std::size_t> cvhull_pts = convex_hull(bpts2D);
    assert(cvhull_pts.size() > 2);

    // adjust positions of points that are inside the convex hull so that they
    // reach the convex hull boundary, ensuring that the fracture outline will be
    // convex
    for (std::size_t i = 0; i != cvhull_pts.size(); ++i) {
        const std::size_t cstart = cvhull_pts[i];
        const std::size_t cend = cvhull_pts[(i + 1) % cvhull_pts.size()];

        // adjust boundary points in the interval between bp_start_ix and bp_end_ix
        for (std::size_t j = (cstart + 1) % N; j != cend; j = (j + 1) % N) {
            const double alpha = adjust_disp_to_convex(
                &bpts2D[2 * cstart], &bpts2D[2 * cend], &bpts2D_old[2 * j], &dirs2D[2 * j]);

            const double diff = alpha - disp[j];

            disp[j] = alpha;
            total_disp[j] += diff;
        }
    }
}

// ----------------------------------------------------------------------------
void
GridStretcher::dumpToVTK(const char* filename, const std::vector<std::vector<double>>& data) const
// ----------------------------------------------------------------------------
{
    auto vtkwriter
        = Dune::VTKWriter<Grid::LeafGridView> {grid_.leafGridView(), Dune::VTK::nonconforming};

    for (int i = 0; i != data.size(); ++i) {
        vtkwriter.addCellData(data[i], "data" + std::to_string(i));
    }

    vtkwriter.write(filename);
}

} // end namespace Opm
