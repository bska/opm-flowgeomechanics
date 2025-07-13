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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace
{

int
test_grid_refinement()
{
    // make mesh
    const auto cells = std::map<Opm::CellRef, Opm::CellAttributes> {
        std::pair {Opm::CellRef {0, 0, 0}, Opm::CellAttributes {}},
        std::pair {Opm::CellRef {0, 0, 1}, Opm::CellAttributes {}},
        std::pair {Opm::CellRef {1, 0, 0}, Opm::CellAttributes {}},
        std::pair {Opm::CellRef {1, 1, 0}, Opm::CellAttributes {}},
        std::pair {Opm::CellRef {1, 1, 1}, Opm::CellAttributes {}},
    };

    const auto mesh = Opm::RegularTrimesh {cells};

    // make refined mesh
    const auto mesh2 = mesh.refine();

    // make more refined grid
    const auto mesh3 = mesh2.refine().refine().refine();

    // make coarse grid
    const auto mesh4 = mesh2.coarsen();

    // export both meshes to vtk
    writeMeshToVTK(mesh, "initial_grid");
    writeMeshToVTK(mesh2, "refined_grid");
    writeMeshToVTK(mesh3, "more_refined_grid", 1); //  coarsen interior
    writeMeshToVTK(mesh4, "coarse_grid");

    return 0;
}

int
test_circular(const double radius, const int levels, const bool center_well, const bool smoothing)
{
    const auto mesh = Opm::RegularTrimesh {radius};

    writeMeshToVTK(mesh,
                   "circular_grid",
                   levels,
                   center_well ? Opm::RegularTrimesh::inner_ring_cells() : std::vector<Opm::CellRef> {},
                   smoothing);

    std::cout << "Number of cells: " << mesh.numActive() << std::endl;

    return 0;
}

int
expand_grid_test(const int turns)
{
    auto mesh = Opm::RegularTrimesh {};

    for (int i = 0; i != turns; ++i) {
        const auto cells = mesh.cellIndices();
        mesh.expandGrid(cells);
        mesh.removeSawtooths();
    }

    writeMeshToVTK(mesh, "expandedgrid");

    std::cout << "Number of cells: " << mesh.numActive() << std::endl;

    return 0;
}

int
irregular_grid_test()
{
    const auto cellspec = std::array {
        std::array {0, 0, 0},
        std::array {1, 0, 0},
        std::array {0, 0, 1},
        std::array {1, 1, 1},
        std::array {1, 1, 0},
    };

    const auto mesh = Opm::RegularTrimesh {cellspec.begin(), cellspec.end()};

    {
        const auto cells = mesh.cellIndices();

        std::cout << "Cells:\n";
        for (const auto& cell : cells) {
            std::cout << cell[0] << ' ' << cell[1] << ' ' << cell[2] << '\n';
        }
        std::cout << std::endl;
    }

    {
        const auto edges = mesh.edgeIndices();

        std::cout << "Edges:\n";
        for (const auto& edge : edges) {
            std::cout << edge[0] << ' ' << edge[1] << ' ' << edge[2] << '\n';
        }
        std::cout << std::endl;
    }

    {
        const auto nodes = mesh.nodeIndices();

        std::cout << "Nodes:\n";
        for (const auto& node : nodes) {
            std::cout << node[0] << ' ' << node[1] << '\n';
        }
        std::cout << std::endl;
    }

    const auto boundary_edges = mesh.boundaryEdges();
    {
        std::cout << "Boundary edges:\n";
        for (const auto& edge : boundary_edges) {
            std::cout << edge[0] << ' ' << edge[1] << ' ' << edge[2] << '\n';
        }
        std::cout << std::endl;
    }

    {
        const auto cell_centroids = mesh.cellCentroids();

        std::cout << "Cell centroids:\n";
        for (const auto& centroid : cell_centroids) {
            std::cout << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2] << '\n';
        }
        std::cout << std::endl;
    }

    {
        const auto edge_centroids = mesh.edgeCentroids();

        std::cout << "Edge centroids:\n";
        for (const auto& centroid : edge_centroids) {
            std::cout << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2] << '\n';
        }
        std::cout << std::endl;
    }

    {
        const auto node_coords = mesh.nodeCoords();

        std::cout << "Node coords:\n";
        for (const auto& coord : node_coords) {
            std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << '\n';
        }
        std::cout << std::endl;
    }

    {
        // write to file as a matlab triangulation (open file as fstream)
        std::ofstream file("mesh.m");
        mesh.writeMatlabTriangulation(file);
    }

    // write grid to file
    {
        const auto& [grid, fsmap, bcells] = mesh.createDuneGrid();

        auto vtkwriter = Dune::VTKWriter {grid->leafGridView(), Dune::VTK::nonconforming};
        vtkwriter.write("mesh");
    }

    Opm::writeMeshBoundaryToVTK(mesh, "boundary.vtu");

    return 0;
}

} // Anonymous namespace

// write a test progam
int
main(int varnum, char** vararg)
{
    if (varnum == 1) {
        std::cout << "Options are:\n"
                  << "1 - create irregular 5-cell grid\n"
                  << "2 - test grid expansion <n turns>\n"
                  << "3 - test grid refinement\n"
                  << "4 - test circular grid construction <radius> <# of levels> <1/0 "
                     "(presence of center well) <1/0> "
                     "(smoothing)\n"
                  << std::endl;
    } else if (std::atoi(vararg[1]) == 1) {
        return irregular_grid_test();
    } else if (std::atoi(vararg[1]) == 2) {
        return expand_grid_test(atoi(vararg[2]));
    } else if (std::atoi(vararg[1]) == 3) {
        return test_grid_refinement();
    } else if (std::atoi(vararg[1]) == 4) {
        return test_circular(
            std::atof(vararg[2]), std::atoi(vararg[3]), std::atoi(vararg[4]), std::atoi(vararg[5]));
    } else {
        std::cout << "Invalid option given" << std::endl;
    }

    return EXIT_FAILURE;
}
