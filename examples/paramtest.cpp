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

#include <opm/geomech/param_interior.hpp>

#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/matrixmarket.hh>

#include <opm/geomech/GridStretcher.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

using Grid = Dune::FoamGrid<2, 3>;

namespace
{

// ----------------------------------------------------------------------------
int
simpletest(const std::string& fname)
// ----------------------------------------------------------------------------
{
    std::ifstream is(fname);

    int num_p {};
    is >> num_p;

    std::vector<double> points(num_p * 2);
    for (int i = 0; i != num_p * 2; ++i) {
        is >> points[i];
    }

    int num_q {};
    is >> num_q;

    std::vector<double> q(num_q * 2);
    for (int i = 0; i != num_q * 2; ++i) {
        is >> q[i];
    }

    // -------------------------- compute parametrization --------------------------
    std::vector<double> result {};

    Opm::parametrize_interior_2D(points, q, result);

    return 0;
}

double
dist2D(const Dune::FieldVector<double, 3>& p1, const Dune::FieldVector<double, 3>& p2)
{
    return std::hypot(p1[0] - p2[0], p1[1] - p2[1]);
}

// ----------------------------------------------------------------------------
std::vector<int>
boundary_node_indices(const Grid& g)
// ----------------------------------------------------------------------------
{
    // set<int> bnodes;
    auto view = g.leafGridView();

    // using FVec = Dune::FieldVector<double, 3>;
    using CoordType = Dune::FieldVector<double, 3>;

    // register all node coordinates
    std::vector<CoordType> vcoords {};
    for (const auto& v : vertices(view)) {
        vcoords.push_back(v.geometry().corner(0));
    }

    // determine boundary edges
    std::vector<int> count(view.size(1), 0); // number of edges (codim 1)

    for (const auto& elem : elements(view)) {
        const auto nSub = elem.subEntities(1);

        for (auto i = 0 * nSub; i != nSub; ++i) {
            ++count[view.indexSet().index(elem.subEntity<1>(i))];
        }
    }

    // @@ This is silly - there must be a better way to identify boundary nodes
    // than by geometric comparison.
    const double TOL = 1e-3;
    int ix = 0;
    std::set<int> bix {};
    for (auto ep = view.begin<1>(); ep != view.end<1>(); ++ep, ++ix) {
        if (count[ix] >= 2) {
            continue;
        }

        for (int i = 0; i != 2; ++i) {
            const auto pt = ep->geometry().corner(i);

            bool found = false;
            for (auto j = 0 * vcoords.size(); j != vcoords.size() && !found; ++j) {
                if (dist2D(pt, vcoords[j]) < TOL) {
                    found = true;
                    bix.insert(j);
                }
            }
        }
    }

    return {bix.begin(), bix.end()};
}

// ----------------------------------------------------------------------------
int
meshtest(const std::string& fname)
// ----------------------------------------------------------------------------
{
    auto grid = Dune::GmshReader<Grid>::read(fname);

    // extract coordinates

    auto view = grid->leafGridView();

    for (int i = 0; i != 3; ++i) {
        std::cout << view.size(i) << '\n';
    }

    // extract all node coordinates (2D only)
    std::vector<double> coords {};
    for (const auto& vertex : vertices(view)) {
        coords.push_back(vertex.geometry().corner(0)[0]);
        coords.push_back(vertex.geometry().corner(0)[1]);
    }

    // identify boundary nodes
    const auto bindices = boundary_node_indices(*grid);

    // identifying interior nodes
    std::vector<int> tmp(coords.size() / 2, 0);
    for (const auto& b : bindices) {
        tmp[b] = 1;
    }

    std::vector<int> iindices {};
    for (auto ix = 0 * tmp.size(); ix != tmp.size(); ++ix) {
        if (tmp[ix] == 0) {
            iindices.push_back(ix);
        }
    }

    // sort node coordinates into boundary and interior
    std::vector<double> bcoords {};
    std::vector<double> icoords {};
    for (auto i = 0 * tmp.size(); i != tmp.size(); ++i) {
        if (tmp[i] == 0) {
            // interior
            icoords.push_back(coords[2 * i]);
            icoords.push_back(coords[2 * i + 1]);
        } else {
            // boundary
            bcoords.push_back(coords[2 * i]);
            bcoords.push_back(coords[2 * i + 1]);
        }
    }

    // call parametrization function
    std::vector<double> result {};
    Opm::parametrize_interior_2D(bcoords, icoords, result);

    std::cout << "Num boundary vertices: " << bcoords.size() / 2 << '\n'
              << "Num interior vertices: " << icoords.size() / 2 << '\n'
              << "Num params: " << result.size() << '\n';

    for (const auto& xi : result) {
        std::cout << xi << '\n';
    }

    // write result to disk
    auto vtkwriter = Dune::VTKWriter {grid->leafGridView(), Dune::VTK::nonconforming};
    vtkwriter.write("undeformed");

    // modify grid coordinates
    const int Nb = bcoords.size() / 2; // num boundary points
    for (int i = 0; i != Nb; ++i) {
        bcoords[2 * i] += 3; // translate in x direction
        if (i < Nb / 3) {
            bcoords[2 * i + 1] *= 2; // stretch in y direction
        }
    }

    const int Ni = static_cast<int>(icoords.size()) / 2; // num interior points
    std::fill(icoords.begin(), icoords.end(), 0); // reset icoords

    for (int i = 0; i != Ni; ++i) {
        for (int b = 0; b != Nb; ++b) {
            for (int c = 0; c != 2; ++c) {
                icoords[2 * i + c] += bcoords[2 * b + c] * result[Nb * i + b];
            }
        }
    }

    const int N = coords.size() / 2;
    int bc = 0;
    int ic = 0;
    for (int i = 0; i != N; ++i) {
        for (int c = 0; c != 2; ++c) {
            coords[2 * i + c] = (tmp[i] == 0) ? icoords[2 * ic + c] : bcoords[2 * bc + c];
        }

        if (tmp[i] == 0) {
            ++ic;
        } else {
            ++bc;
        }
    }

    int vcount = 0;
    for (const auto& vertex : vertices(view)) {
        grid->setPosition(
            vertex, Dune::FieldVector<double, 3> {coords[2 * vcount], coords[2 * vcount++ + 1], 0});

        std::cout << vertex.geometry().corner(0)[0] << ", " << vertex.geometry().corner(0)[1] << '\n';
    }

    vtkwriter.write("deformed");

    // repeat deformation, this time using a grid stretcher
    auto grid2 = Dune::GmshReader<Grid>::read(fname);

    Opm::GridStretcher gs {*grid2};

    std::vector<double> coords2 {};
    for (const auto& vertex : vertices(grid2->leafGridView())) {
        coords2.push_back(vertex.geometry().corner(0)[0]);
        coords2.push_back(vertex.geometry().corner(0)[1]);
    }

    std::vector<Opm::GridStretcher::CoordType> disp(Nb, {0, 0, 0});
    for (int i = 0; i != Nb; ++i) {
        disp[i][0] += 6; // translate in x direction
        if (i < Nb / 3) {
            disp[i][1] = coords2[2 * i + 1]; // bcoords[2*i+1]/2; // stretch in y direction
        }
    }

    gs.applyBoundaryNodeDisplacements(disp);

    auto vtkwriter2 = Dune::VTKWriter {grid2->leafGridView(), Dune::VTK::nonconforming};
    vtkwriter2.write("deformed2");

    return 0;
}

// ----------------------------------------------------------------------------
int
projectiontest(const std::string& fname, const std::string& ftarget)
// ----------------------------------------------------------------------------
{
    // open file and read 3D points to a std::vector<double>
    std::ifstream is(fname);

    std::vector<double> points;
    while (!is.eof()) {
        double x, y, z;
        is >> x >> y >> z;
        if (is.eof()) {
            break;
        }

        points.push_back(x);
        points.push_back(y);
        points.push_back(z);
    }

    // call projection function
    const auto N = points.size() / 3;
    std::vector<double> result(N * 2);
    const auto axis = Opm::project_to_2D(points, result);

    // write axis, then result to file
    std::ofstream os(ftarget);
    os.precision(16);
    for (int i = 0; i != 3; ++i) {
        for (int dim = 0; dim != 3; ++dim) {
            os << axis[i][dim] << ' ';
        }

        os << '\n';
    }

    for (auto i = 0 * N; i != N; ++i) {
        os << result[2 * i] << ' ' << result[2 * i + 1] << '\n';
    }

    return 0;
}

//----------------------------------------------------------------------------
int
lifttest(const std::string& fname, const std::string& target)
//----------------------------------------------------------------------------
{
    Opm::Axis3D axis {};
    std::vector<double> points {};

    // reading axis and 2D points from file
    {
        std::ifstream is(fname);

        for (int i = 0; i != 3; ++i) {
            for (int dim = 0; dim != 3; ++dim) {
                is >> axis[i][dim];
            }
        }

        double x {};
        double y {};

        while (is >> x >> y) {
            points.push_back(x);
            points.push_back(y);
        }
    }

    // lift the points to 3D, using the axis
    std::vector<double> result {};
    Opm::lift_to_3D(points, axis, result);

    // write result to file, using high precision
    std::ofstream os(target);
    os.precision(16);
    for (auto i = 0 * result.size(); i != result.size() / 3; ++i) {
        os << result[3 * i] << ' ' << result[3 * i + 1] << ' ' << result[3 * i + 2] << '\n';
    }

    return 0;
}

//----------------------------------------------------------------------------
int
loop_repar_test(const std::string& fname, const std::string& target)
// ----------------------------------------------------------------------------
{
    std::vector<double> points {};

    {
        // reading 2D points from dfile
        std::ifstream is(fname);

        double x {};
        double y {};

        while (is >> x >> y) {
            points.push_back(x);
            points.push_back(y);
        }
    }

    // redistribute points
    std::vector<double> result {};
    Opm::redistribute_2D(points, result);

    // write result to file
    std::ofstream os(target);
    for (auto i = 0 * result.size(); i != result.size() / 2; ++i) {
        os << result[2 * i] << ' ' << result[2 * i + 1] << '\n';
    }

    return 0;
}

} // namespace

// ============================================================================
int
main([[maybe_unused]] int varnum, char** vararg)
// ============================================================================
{
    const auto fname = std::string {vararg[1]};
    const auto ftype = fname.substr(fname.size() - 3);

    if (ftype == "txt") {
        return simpletest(fname);
    } else if (ftype == "msh") {
        return meshtest(fname);
    } else if (ftype == "p3d") {
        return projectiontest(fname, vararg[2]);
    } else if (ftype == "p2d") {
        return lifttest(fname, vararg[2]);
    } else if (ftype == "lop") {
        return loop_repar_test(fname, vararg[2]);
    }
}
