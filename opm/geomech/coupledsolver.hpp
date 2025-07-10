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

#ifndef COUPLEDSOLVER_HPP_INCLUDED
#define COUPLEDSOLVER_HPP_INCLUDED

#include <dune/common/fmatrix.hh> // Dune::FieldMatrix

#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector

#include <cstddef>
#include <tuple>
#include <vector>

namespace // anonymous namespace
{
using Htrans = std::tuple<std::size_t, std::size_t, double, double>;
using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
using FMatrix = Dune::DynamicMatrix<double>; // full matrix
} // namespace

namespace Opm
{
void solve_fully_coupled(ResVector& pvec, // output: fracture pressure
                         ResVector& hvec, // output: aperture
                         const SMatrix& pmat, // pressure matrix (sparse)
                         const FMatrix& amat, // aperture matrix (full)
                         const std::vector<Htrans>& htrans,
                         const double rate,
                         const std::vector<std::size_t>& ratecells,
                         const double bhp,
                         const std::vector<std::size_t>& bhpcells,
                         const std::vector<double>& leakoff_fac,
                         const int max_nonlin_iter = 100,
                         const double conv_tol = 1e-5);

} // end namespace Opm

#endif // COUPLEDSOLVER_HPP_INCLUDED
