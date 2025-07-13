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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/geomech/vem/vem.hpp>

#include <dune/grid/yaspgrid.hh>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <iterator> // for ostream_iterator
#include <vector>

namespace
{

using GridView = Dune::YaspGrid<3>::LeafGridView;
using Cell = Dune::Entity<0, 3, const Dune::YaspGrid<3>, Dune::YaspEntity>;
using Face = Dune::Entity<1, 3, const Dune::YaspGrid<3>, Dune::YaspEntity>;
using IndexSet = Dune::IndexSet<const Dune::YaspGrid<3>,
                                Dune::YaspIndexSet<const Dune::YaspGrid<3>, true>,
                                unsigned int,
                                std::vector<Dune::GeometryType>>;

// ----------------------------------------------------------------------------
std::vector<int>
make_nodemap(const IndexSet& ixset, const Cell& cell)
// ----------------------------------------------------------------------------
{
    std::vector<int> result;

    const int N = Dune::subEntities(cell, Dune::Codim<3> {}).size();
    result.reserve(N);

    for (int i = 0; i != N; ++i) {
        result.push_back(ixset.subIndex(cell, i, 3));
    }

    return result;
}

// ----------------------------------------------------------------------------
std::vector<double>
corner_coords(const Cell& cell)
// ----------------------------------------------------------------------------
{
    std::vector<double> result;

    const int N = Dune::subEntities(cell, Dune::Codim<3> {}).size();

    const auto cellgeo = cell.geometry();

    for (int i = 0; i != N; ++i) {
        const auto cor = cellgeo.corner(i);
        result.insert(result.end(), cor.begin(), cor.end());
    }

    return result;
}

} // Anonymous namespace

// ----------------------------------------------------------------------------
int
main(int argc, char** argv)
// ----------------------------------------------------------------------------
{
    // to prevent compiler from complaining about MPI ??
    Dune::MPIHelper::instance(argc, argv);

    using Grid = Dune::YaspGrid<3>;
    const Grid grid({4.0, 4.0, 4.0}, // physical dimensions
                    {2, 2, 2}); // resolution

    const auto gv = grid.leafGridView();
    const auto& ixset = gv.indexSet();

    const double young = 1;
    const double poisson = 0.25;

    for (const auto& c : elements(gv)) {
        // make local-to-global map of cell corners
        const auto nodemap = make_nodemap(ixset, c);

        // define lookup function from local face corner indexing to local cell corner
        // indexing
        const auto nodemap_lookup = [&nodemap, &ixset](const Face& f, const int i) {
            auto pos = std::find(nodemap.begin(), nodemap.end(), ixset.subIndex(f, i, 3));
            return std::distance(nodemap.begin(), pos);
        };

        // create 'faces' vector
        std::vector<int> faces {};
        for (const auto& f : Dune::subEntities(c, Dune::Codim<1> {})) {
            for (int i = 0; i != f.geometry().corners(); ++i) {
                faces.push_back(nodemap_lookup(f, i));
            }
        }

        // create 'num_face_corners' vector
        std::vector<int> num_face_corners {};
        for (const auto& f : Dune::subEntities(c, Dune::Codim<1> {})) {
            num_face_corners.push_back(f.geometry().corners());
        }

        // correct order of nodes in each face, to ensure they are mentioned in
        // clockwise or counterclockwise order. @@ This is a hack that might only work
        // for 4-faces!
        for (int i = 0, j = 0; i != static_cast<int>(num_face_corners.size());
             j += num_face_corners[i++]) {
            std::swap(faces[j], faces[j + 1]);
        }

        // create local vector of point coordinates
        const std::vector<double> coords = corner_coords(c);

        // assemble element stiffness matrix
        std::array<double, 3> cell_centroid {};
        std::vector<int> indexing {};
        std::vector<double> target(24 * 24);

        const auto stability_choice = vem::StabilityChoice::D_RECIPE;

        vem::assemble_stiffness_matrix_3D(&coords[0],
                                          &faces[0],
                                          &num_face_corners[0],
                                          static_cast<int>(num_face_corners.size()),
                                          young,
                                          poisson,
                                          stability_choice,
                                          cell_centroid,
                                          indexing,
                                          target);

        std::cout << "Centroid: " << cell_centroid[0] << " " << cell_centroid[1] << " "
                  << cell_centroid[2] << '\n'
                  << "Indexing: " << '\n';
        std::copy(indexing.begin(), indexing.end(), std::ostream_iterator<int>(std::cout, " "));

        std::cout << "\nMatrix:\n";
        vem::matprint(&target[0], indexing.size() * 3, indexing.size() * 3, false, 1.0e-9);
        std::cout << '\n';

        std::cout << "Faces are:\n";
        std::copy(faces.begin(), faces.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << '\n';

        std::cout << "Num face corners are:\n";
        std::copy(num_face_corners.begin(),
                  num_face_corners.end(),
                  std::ostream_iterator<int>(std::cout, " "));
        std::cout << '\n';

        std::cout << "Coordinates are:\n";
        for (std::size_t i = 0; i < coords.size(); i += 3) {
            std::cout << coords[i] << " " << coords[i + 1] << " " << coords[i + 2] << '\n';
        }
    }
}
