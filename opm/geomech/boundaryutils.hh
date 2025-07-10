// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#ifndef BOUNDARYUTILS_HH
#define BOUNDARYUTILS_HH

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <array>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace Opm::Elasticity
{
template <int dimension>
unsigned
cartesianIndex(const std::array<int, dimension>& coords,
               const std::array<int, dimension>& cartesianDimensions)
{
    unsigned cartIndex = coords[0];
    int factor = cartesianDimensions[0];
    for (unsigned i = 1; i < dimension; ++i) {
        cartIndex += coords[i] * factor;
        factor *= cartesianDimensions[i];
    }

    return cartIndex;
}

inline FaceDir::DirEnum
faceToFaceDir(int insideFaceIdx)
{
    switch (insideFaceIdx) {
    case 0:
        return FaceDir::XMinus;
    case 1:
        return FaceDir::XPlus;
    case 2:
        return FaceDir::YMinus;
    case 3:
        return FaceDir::YPlus;
    case 4:
        return FaceDir::ZMinus;
    case 5:
        return FaceDir::ZPlus;
    default:
        OPM_THROW(std::logic_error, "Internal error in initialization of aquifer.");
    }
}

inline int
faceDirToFace(Opm::FaceDir::DirEnum faceDirection)
{
    switch (faceDirection) {
    case FaceDir::XMinus:
        return 0;
    case FaceDir::XPlus:
        return 1;
    case FaceDir::YMinus:
        return 2;
    case FaceDir::YPlus:
        return 3;
    case FaceDir::ZMinus:
        return 4;
    case FaceDir::ZPlus:
        return 5;
    default:
        OPM_THROW(std::logic_error, "Internal error in initialization of aquifer.");
    }
}

inline std::array<int, 4>
faceDirToNodes(Opm::FaceDir::DirEnum faceDirection)
{
    switch (faceDirection) {
    case FaceDir::XMinus:
        return {0, 4, 6, 2};
    case FaceDir::XPlus:
        return {1, 3, 7, 5};
    case FaceDir::YMinus:
        return {0, 1, 5, 4};
    case FaceDir::YPlus:
        return {3, 2, 6, 7};
    case FaceDir::ZMinus:
        return {0, 2, 3, 1};
    case FaceDir::ZPlus:
        return {4, 5, 7, 6};
    default:
        OPM_THROW(std::logic_error, "Internal error in initialization of aquifer.");
    }
}


template <class BCConfig, class GvType, class CartMapperType>
void
nodesAtBoundary(std::vector<std::tuple<std::size_t, MechBCValue>>& bc_nodes,
                const BCConfig& bcconfigs,
                const BCProp& bcprops,
                const GvType& gv,
                const CartMapperType& cartesianIndexMapper)
{
    constexpr auto dim = 3;

    if (bcprops.size() > 0) {
        // nonTrivialBoundaryConditions_ = true;
        const std::size_t numCartDof = cartesianIndexMapper.cartesianSize();
        const unsigned numElems = gv.size(/*codim=*/0);

        std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            cartesianToCompressedElemIdx[cartesianIndexMapper.cartesianIndex(elemIdx)] = elemIdx;
        }

        std::array<int, 3> cartdim = cartesianIndexMapper.cartesianDimensions();
        for (const auto& bcconfig : bcconfigs) {
            for (const auto& bcprop : bcprops) {
                if (bcprop.index == bcconfig.index) {
                    // double search since structure is strange
                    const auto& bcface = bcconfig;

                    if ((bcface.i1 < 0) || (bcface.j1 < 0) || (bcface.k1 < 0)) {
                        throw std::logic_error("Lower range of BC wrong");
                    }

                    if ((bcface.i2 > cartdim[0]) || (bcface.j2 > cartdim[1])
                        || (bcface.k2 > cartdim[2])) {
                        throw std::logic_error("Upper range of BC wrong");
                    }

                    const auto& type = bcprop.bcmechtype;
                    if (type == Opm::BCMECHType::FREE) {
                        // do nothing
                    } else if (type == Opm::BCMECHType::FIXED) {
                        std::set<std::size_t> effected_cells;
                        for (int i = bcface.i1; i <= bcface.i2; ++i) {
                            for (int j = bcface.j1; j <= bcface.j2; ++j) {
                                for (int k = bcface.k1; k <= bcface.k2; ++k) {
                                    std::array<int, 3> tmp = {i, j, k};
                                    int cartindex = cartesianIndex<3>(tmp, cartdim);
                                    auto elemIdx = cartesianToCompressedElemIdx[cartindex];
                                    if (elemIdx > -1) {
                                        effected_cells.insert(elemIdx);
                                    }
                                }
                            }
                        }

                        MechBCValue bcval = *bcprop.mechbcvalue;
                        for (const auto& cell : elements(gv)) {
                            auto it = effected_cells.find(gv.indexSet().index(cell));
                            if (it != effected_cells.end()) {
                                // fix all noted for now
                                std::array<int, 4> nodes = faceDirToNodes(bcface.dir);
                                for (const auto& nind : nodes) {
                                    auto global_ind = gv.indexSet().subIndex(cell, nind, dim);
                                    bc_nodes.emplace_back(global_ind, bcval);
                                }
                            }
                        }
                    } else {
                        throw std::logic_error("invalid type for BC. Use FREE or RATE");
                    }
                }
            }
        }
    }

    auto compare = [](std::tuple<std::size_t, MechBCValue> const& t1,
                      std::tuple<std::size_t, MechBCValue> const& t2) {
        return std::get<0>(t1) < std::get<0>(t2);
    };

    auto isequal = [](std::tuple<std::size_t, MechBCValue> const& t1,
                      std::tuple<std::size_t, MechBCValue> const& t2) {
        return std::get<0>(t1) == std::get<0>(t2);
    };

    std::sort(bc_nodes.begin(), bc_nodes.end(), compare); // {1 1 2 3 4 4 5}
    auto last = std::unique(bc_nodes.begin(), bc_nodes.end(), isequal);

    // v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
    bc_nodes.erase(last, bc_nodes.end());
}

} // namespace Opm::Elasticity

#endif // BOUNDARYUTILS_HH
