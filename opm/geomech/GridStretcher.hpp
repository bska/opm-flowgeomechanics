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

#ifndef OPM_GRIDSTRETCHER_HPP_INCLUDED
#define OPM_GRIDSTRETCHER_HPP_INCLUDED

#include <dune/common/fvector.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/entityseed.hh>

#include <algorithm>
#include <cstddef>
#include <map>
#include <tuple>
#include <vector>

namespace Opm
{

struct BoundaryNormals
{
    std::vector<Dune::FieldVector<double, 3>>
        bnode_normals; // outward, unitary normals of boundary nodes
    std::vector<Dune::FieldVector<double, 3>>
        bcell_normals; // direction from cell centroid to boundary centroid
    std::vector<Dune::FieldVector<double, 3>> bedge_normals; // normal to boundary edge
};

class GridStretcher
{
public:
    using Grid = Dune::FoamGrid<2, 3>;
    using CellSeed = Grid::Codim<0>::EntitySeed;
    using CellBnodeMap = std::map<std::size_t, std::tuple<CellSeed, std::size_t, std::size_t>>;
    using CoordType = Dune::FieldVector<double, 3>;

    explicit GridStretcher(Grid& grid);

    const std::vector<std::size_t>& boundaryNodeIndices() const
    {
        return bnindices_;
    }

    const std::vector<std::size_t>& boundaryCellIndices() const
    {
        return bcindices_;
    }

    // the vector should have one entry per boundary cell, expressing the distance
    // the boundary of that cell should be expanded outwards
    void expandBoundaryCells(const std::vector<double>& amounts); // will modify grid

    // the vector should have two enties per boundary node, specifying its displacement
    // in the x and y direction
    void applyBoundaryNodeDisplacements(const std::vector<CoordType>& disp);

    // make boundary nodes equidistant
    void rebalanceBoundary(); // will modify grid

    std::vector<double> computeBoundaryNodeDisplacements(const std::vector<double>& amounts,
                                                         const std::vector<CoordType>& bnodenormals
                                                         = std::vector<CoordType>()) const;

    std::vector<double> centroidEdgeDist() const;

    const std::vector<CoordType>& nodecoords() const
    {
        return nodecoords_;
    }

    const std::vector<CoordType>& bnodenormals() const
    {
        return boundary_normals_.bnode_normals;
    }

    const std::vector<CoordType>& bcellnormals() const
    {
        return boundary_normals_.bcell_normals;
    }

    const std::vector<CoordType>& bedgenormals() const
    {
        return boundary_normals_.bedge_normals;
    }

    double maxBoxLength() const;

    // objective function and derivatives, when trying to stretch grid to a particular
    // target (in terms of distances between cell and edge centroids for boundary cells)
    double objective(const std::vector<double>& bndisp,
                     const std::vector<double>& dtarget,
                     std::vector<double>& grad,
                     bool fixed_cell_centroids = false);

    void adjustToConvex(std::vector<double>& disp,
                        std::vector<double>& total_disp,
                        const std::vector<CoordType>& dirs) const;

    void dumpToVTK(const char* filename,
                   const std::vector<std::vector<double>>& = std::vector<std::vector<double>>()) const;

    void updateAfterBoundaryChange(const std::vector<CoordType>& new_bcoords);

private:
    // ------------------------------- internal data -------------------------------

    Grid& grid_; // NB: mutable reference to externally owned grid!
    std::vector<CoordType> nodecoords_; // NB! should be updated when grid is updated.

    const std::vector<std::size_t> bnindices_; // indices of boundary gridnodes
    const std::vector<std::size_t> iindices_; // indices of internal gridnodes
    const std::vector<double> iparam_; // parametrization of internal nodes in
                                       // terms of boundary nodes
    const CellBnodeMap c2bix_; // map cell -> bindices (local indexing for boundary nodes)
    const std::vector<std::size_t> bcindices_; // indices of boundary cells
    const std::vector<double> bcentroid_param_; // parametrization of boundary cell centroids

    BoundaryNormals boundary_normals_; // NB! should be updated when grid is updated.
};

} // end namespace Opm

#endif // OPM_GRIDSTRETCHER_HPP_INCLUDED
