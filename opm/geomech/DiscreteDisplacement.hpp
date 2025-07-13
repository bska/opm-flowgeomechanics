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

#ifndef DISCRETE_DISPLACEMENT_HPP_INCLUDED
#define DISCRETE_DISPLACEMENT_HPP_INCLUDED

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/bvector.hh>

#include <opm/geomech/CutDe.hpp>
#include <opm/geomech/Math.hpp>

#include <array>

namespace ddm
{
template <class Element>
Dune::FieldVector<double, 3>
normalOfElement(const Element& elem)
{
    // dim = Dune::Codim<Grid::dimension>{}
    // int codim = 2;//vertices in this cases
    std::array<Dune::FieldVector<double, 3>, 3> tri;
    int i = 0;
    for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<2> {})) {
        // assume order is ok and orientation is fine;
        tri[i] = vertex.geometry().center();
        ++i;
    }

    // NB sign should not matter as long as tensile fracture
    Dune::FieldVector<double, 3> e01 = tri[1] - tri[0];
    Dune::FieldVector<double, 3> e02 = tri[2] - tri[0];
    Dune::FieldVector<double, 3> normal = Opm::crossProduct(e01, e02);

    normal /= normal.two_norm();

    return normal;
}

double fractureK1(double dist, double width, double E, double nu);

// double fractureK2(double dist,double edist, double E, double nu);
// double fractureK3(double dist,double edist, double E, double nu);

template <class Element>
std::array<Real3, 3>
getTri(const Element& elem)
{
    // dim = Dune::Codim<Grid::dimension>{}
    // int codim = 2;//vertices in this cases
    std::array<Real3, 3> tri;

    int i = 0;
    for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<2> {})) {
        // assume order is ok and orientation is fine;
        const auto& corner = vertex.geometry().center();
        tri[i] = make3(corner[0], corner[1], corner[2]);
        ++i;
    }

    return tri;
}

template <class Element>
Dune::FieldVector<double, 3>
TDDispFS(const Dune::FieldVector<double, 3>& obs,
         const Element& elem,
         const Dune::FieldVector<double, 3>& slip,
         double nu)
{
    Dune::FieldVector<double, 3> disp;

    const std::array<Real3, 3> tri = getTri(elem);

    const Real3 obs_tmp = make3(obs[0], obs[1], obs[2]);
    const Real3 slip_tmp = make3(slip[0], slip[1], slip[2]);
    const Real3 disp_tmp = disp_fs(obs_tmp, tri, slip_tmp, nu);

    disp[0] = disp_tmp.x;
    disp[1] = disp_tmp.y;
    disp[2] = disp_tmp.z;

    return disp;
}

template <class Element>
const Dune::FieldVector<double, 6>
TDStrainFS(const Dune::FieldVector<double, 3>& obs,
           const Element& elem,
           const Dune::FieldVector<double, 3>& slip,
           double nu)
{
    Dune::FieldVector<double, 6> strain;

    const std::array<Real3, 3> tri = getTri(elem);
    const Real3 obs_tmp = make3(obs[0], obs[1], obs[2]);
    const Real3 slip_tmp = make3(slip[0], slip[1], slip[2]);
    const Real6 strain_tmp = strain_fs(obs_tmp, tri, slip_tmp, nu);

    strain[0] = strain_tmp.x;
    strain[1] = strain_tmp.y;
    strain[2] = strain_tmp.z;
    strain[3] = strain_tmp.a;
    strain[4] = strain_tmp.b;
    strain[5] = strain_tmp.c;

    return strain;
}

double traceSymTensor(const Dune::FieldVector<double, 6>& symtensor);

Dune::FieldMatrix<double, 3, 3> symTensor2Matrix(const Dune::FieldVector<double, 6>& symtensor);

// compute the stress tensor (6 components) from the strain tensor (6 components)
Dune::FieldVector<double, 6>
strainToStress(const double E, const double nu, const Dune::FieldVector<double, 6>& strain);

double tractionSymTensor(const Dune::FieldVector<double, 6>& symtensor,
                         const Dune::FieldVector<double, 3>& normal);

// assembleMatrix(Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>& matrix, const double
// E, const double nu, const Dune::FoamGrid<2, 3>& grid)
void assembleMatrix(Dune::DynamicMatrix<double>& matrix,
                    const double E,
                    const double nu,
                    const Dune::FoamGrid<2, 3>& grid);

Dune::FieldVector<double, 6> strain(const Dune::FieldVector<double, 3>& obs,
                                    const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
                                    const Dune::FoamGrid<2, 3>& grid,
                                    const double E,
                                    const double nu);

Dune::FieldVector<double, 3> disp(const Dune::FieldVector<double, 3>& obs,
                                  const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
                                  const Dune::FoamGrid<2, 3>& grid,
                                  const double E,
                                  const double nu);
} // namespace ddm

#endif // DISCRETE_DISPLACEMENT_HPP_INCLUDED
