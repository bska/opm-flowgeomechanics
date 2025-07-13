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

#include <opm/geomech/DiscreteDisplacement.hpp>

#include <cmath>

namespace ddm
{
double
fractureK1(const double dist, const double width, const double E, const double nu)
{
    const auto mu = E / (2 * (1 + nu)); //??

    auto K1 = (mu * std::sqrt(M_PI) / (2 * std::sqrt(dist) * (1.0 - nu))) * width;
    K1 /= 1.834; // factor found numerically //@@ Odd: changed this from '*' to '/'

    return K1;
}

double
traceSymTensor(const Dune::FieldVector<double, 6>& symtensor)
{
    double trace = 0;

    for (int i = 0; i < 3; ++i) {
        trace += symtensor[i];
    }

    return trace;
}

Dune::FieldVector<double, 6>
strainToStress(const double E, const double nu, const Dune::FieldVector<double, 6>& strain)
{
    Dune::FieldVector<double, 6> stress;

    const double volume_strain = traceSymTensor(strain);
    const double mu = E / (2 * (1 + nu)); //??
    const double lambda = 2 * mu * nu / (1 - 2 * nu);

    for (int i = 0; i < 3; ++i) {
        stress[i] = 2 * mu * strain[i] + (lambda * volume_strain);
        stress[i + 3] += 2 * mu * strain[i + 3]; // [xy,xz,yz]
    }

    return stress;
}

Dune::FieldMatrix<double, 3, 3>
symTensor2Matrix(const Dune::FieldVector<double, 6>& symtensor)
{
    Dune::FieldMatrix<double, 3, 3> mat;

    for (int i = 0; i < 3; ++i) {
        mat[i][i] = symtensor[i];
    }

    mat[2][1] = mat[1][2] = symtensor[3];
    mat[2][0] = mat[0][2] = symtensor[4];
    mat[1][0] = mat[0][1] = symtensor[5];

    return mat;
}

double
tractionSymTensor(const Dune::FieldVector<double, 6>& symtensor,
                  const Dune::FieldVector<double, 3>& normal)
{
    assert(std::abs(normal.two_norm() - 1) < 1e-13);

    double traction = 0.0;
    for (int i = 0; i < 3; ++i) {
        traction += symtensor[i] * normal[i] * normal[i];
    }

    traction += 2 * symtensor[3] * normal[0] * normal[1]; // xy*nx*ny;
    traction += 2 * symtensor[4] * normal[0] * normal[2]; // xz*nx*nz
    traction += 2 * symtensor[5] * normal[1] * normal[2]; // yz*ny*nz

    return traction;
}

void
assembleMatrix(Dune::DynamicMatrix<double>& matrix,
               const double E,
               const double nu,
               const Dune::FoamGrid<2, 3>& grid)
{
    using Grid = Dune::FoamGrid<2, 3>;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    for (const auto& elem1 : elements(grid.leafGridView())) {
        const int idx1 = mapper.index(elem1);
        const auto center = elem1.geometry().center();
        const auto normal = normalOfElement(elem1);

        for (const auto& elem2 : elements(grid.leafGridView())) {
            const int idx2 = mapper.index(elem2);
            // check if this is defined in relative coordinates
            Dune::FieldVector<double, 3> slip; // = make3(1.0,0.0, 0.0);
            slip[0] = 1;
            slip[1] = 0;
            slip[2] = 0;

            // symmetric stress voit notation
            const Dune::FieldVector<double, 6> strain = TDStrainFS(center, elem2, slip, nu);
            const Dune::FieldVector<double, 6> stress = strainToStress(E, nu, strain);

            // matrix relate to pure traction not area weighted
            matrix[idx1][idx2] = tractionSymTensor(stress, normal);
        }
    }
}

Dune::FieldVector<double, 6>
strain(const Dune::FieldVector<double, 3>& obs,
       const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
       const Dune::FoamGrid<2, 3>& grid,
       const double /*E*/,
       const double nu)
{
    Dune::FieldVector<double, 6> strain;
    strain = 0.0;

    using Grid = Dune::FoamGrid<2, 3>;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    for (const auto& elem1 : elements(grid.leafGridView())) {
        const int idx1 = mapper.index(elem1);
        strain += TDStrainFS(obs, elem1, slips[idx1], nu);
    }

    return strain;
}

Dune::FieldVector<double, 3>
disp(const Dune::FieldVector<double, 3>& obs,
     const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
     const Dune::FoamGrid<2, 3>& grid,
     const double /*E*/,
     const double nu)
{
    Dune::FieldVector<double, 3> disp;
    disp = 0.0;

    using Grid = Dune::FoamGrid<2, 3>;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    for (const auto& elem1 : elements(grid.leafGridView())) {
        const int idx1 = mapper.index(elem1);
        disp += TDDispFS(obs, elem1, slips[idx1], nu);
    }

    return disp;
}

} // namespace ddm
