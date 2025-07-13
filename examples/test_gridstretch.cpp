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

#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/istl/matrixmarket.hh>

#include <opm/geomech/param_interior.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using Grid = Dune::FoamGrid<2, 3>;
using CoordType = Opm::GridStretcher::CoordType;

namespace
{

void
translate(std::vector<CoordType>& disp, const CoordType& c)
{
    for (auto& e : disp) {
        e += c;
    }
}

void
uniform_scale(std::vector<CoordType>& disp, const std::vector<CoordType>& bcoords, const double fac)
{
    for (std::size_t i = 0; i != disp.size(); ++i) {
        disp[i] += bcoords[i] * (fac - 1);
    }
}

} // end anonymous namespace

int
main([[maybe_unused]] int varnum, char** vararg)
{
    assert(varnum == 3);

    const auto fname = std::string {vararg[1]};
    const auto choice = std::stoi(vararg[2]);
    auto grid = Dune::GmshReader<Grid>::read(fname); // unique_ptr

    std::vector<CoordType> coords;
    for (const auto& vertex : vertices(grid->leafGridView())) {
        coords.push_back({vertex.geometry().corner(0)[0],
                          vertex.geometry().corner(0)[1],
                          vertex.geometry().corner(0)[2]});
    }

    Opm::GridStretcher gs {*grid};

    const auto& bindices = gs.boundaryNodeIndices();

    std::vector<CoordType> bcoords {};
    for (const auto& b : bindices) {
        bcoords.push_back(coords[b]);
    }

    switch (choice) {
    case 1:
        // test boundary node displacements
        {
            std::cout << "Displacement test\n";

            std::vector<CoordType> disp(bindices.size(), {0, 0, 0});

            translate(disp, {4, 0, 0});
            uniform_scale(disp, bcoords, 1.5);

            gs.applyBoundaryNodeDisplacements(disp);
        }
        break;

    case 2:
        // test expansion
        {
            std::cout << "Expanding test\n";

            std::vector<double> amounts(bindices.size(), 0);
            amounts[0] = 0.05;
            amounts[1] = 0.1;
            amounts[3] = 0.1;
            amounts[4] = 0.12;
            amounts[5] = 0.1;
            amounts[11] = 0.05;
            amounts[17] = 0.1;
            amounts[22] = 0.2;

            gs.expandBoundaryCells(amounts);
        }
        break;

    case 3:
        // test gradient
        {
            std::vector<double> grad {};
            std::vector<double> disp(gs.boundaryNodeIndices().size(), 0);

            auto target = gs.centroidEdgeDist();
            std::transform(
                target.begin(), target.end(), target.begin(), [](const double ti) { return ti * 2; });

            const double obj = gs.objective(disp, target, grad, false);

            std::cout << "Objective value: " << obj << '\n';
            std::cout << "Analytical gradient: ";
            std::copy(grad.begin(), grad.end(), std::ostream_iterator<double>(std::cout, " "));
            std::cout << '\n';

            // computing numerical gradient
            const double delta = 1e-6;

            std::vector<double> grad_dummy {};
            for (auto i = 0 * disp.size(); i != disp.size(); ++i) {
                auto dispDelta = disp;
                dispDelta[i] += delta;

                const double obj2 = gs.objective(dispDelta, target, grad_dummy, false);
                std::cout << (obj2 - obj) / delta << " " << grad[i] << '\n';
            }

            std::cout << std::endl;

            return 0;
        }

    case 4: // test optimization
    {
        bool fixed_centroids = true;

        const auto normals = gs.bnodenormals();
        auto target = gs.centroidEdgeDist();

        const double dispfac = 3;
        target[0] *= dispfac;
        target[7] *= dispfac;
        target[17] *= dispfac;
        target[11] *= dispfac;

        std::vector<double> disp(normals.size(), 0);
        std::vector<double> grad;

        double delta = 1;
        const double fac = 2;

        auto objprev = gs.objective(disp, target, grad, fixed_centroids);

        const int MAX_ITER = 100;
        const double tol = 1e-5;
        for (auto i = 0 * MAX_ITER; i != MAX_ITER; ++i) {
            const auto disp_backup = disp;
            const auto grad_backup = grad;

            // compute norm of gradient
            auto gradnorm
                = std::accumulate(grad.begin(), grad.end(), 0.0, [](const double acc, const double gi) {
                      return acc + (gi * gi);
                  });
            gradnorm = std::sqrt(gradnorm);

            // compute steplength
            const double alpha = delta / gradnorm;

            // modify disp

            for (auto j = 0 * disp.size(); j != disp.size(); ++j) {
                disp[j] -= alpha * grad[j];
                disp[j] = std::max(disp[j], 0.0);
            }

            auto obj = gs.objective(disp, target, grad, fixed_centroids);

            if (std::abs(obj - objprev) / obj < tol) {
                break;
            }

            if (obj < objprev) {
                delta *= fac;
                objprev = obj;
                std::cout << i << " : " << obj << '\n';
            } else {
                // reset disp and grad, retry with smaller delta
                disp = disp_backup;
                grad = grad_backup;
                delta /= fac;
                obj = objprev;
                --i;
            }
        }

        // compute the deformed grid geometry
        std::vector<CoordType> mov {};
        mov.reserve(disp.size());
        for (auto i = 0 * disp.size(); i != disp.size(); ++i) {
            mov.push_back(normals[i] * disp[i]);
        }

        gs.applyBoundaryNodeDisplacements(mov);

        break;
    }

    default:
        std::cout << "Choice was not provided.  Must be 1 (displacement) or 2 (expansion).\n";
        return 1;
    }

    auto vtkwriter = Dune::VTKWriter {grid->leafGridView(), Dune::VTK::nonconforming};
    vtkwriter.write("transformed");

    return 0;
}
