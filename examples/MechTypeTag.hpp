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

#ifndef OPM_MECHTYPETAG_HPP_INCLUDED
#define OPM_MECHTYPETAG_HPP_INCLUDED

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <opm/geomech/BlackoilGeomechWellModel.hpp>
#include <opm/geomech/BlackoilModelGeomech.hpp>
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>

namespace Opm::Properties
{
namespace TTag
{
    struct EclFlowProblemMech
    {
        using InheritsFrom = std::tuple<FlowProblem>; //?? should it be blackoil
    };
} // namespace TTag

// Set the problem class
template <class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemMech>
{
    using type = EclProblemGeoMech<TypeTag>;
};

template <class TypeTag>
struct NonlinearSystem<TypeTag, TTag::EclFlowProblemMech>
{
    using type = BlackoilModelGeomech<TypeTag>;
};

template <class TypeTag>
struct EnableMech<TypeTag, TTag::EclFlowProblemMech>
{
    static constexpr bool value = true;
};

template <class TypeTag>
struct EnableEnergy<TypeTag, TTag::EclFlowProblemMech>
{
    static constexpr bool value = true;
};

template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemMech>
{
    static constexpr bool value = false;
};

template <class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemMech>
{
    static constexpr bool value = false;
};


template <class TypeTag>
struct WellModel<TypeTag, TTag::EclFlowProblemMech>
{
    using type = BlackoilGeomechWellModel<TypeTag>;
    // using type = BlackoilWellModel<TypeTag>;
};

template <class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemMech>
{
    using type = Opm::Simulator<TypeTag>;
};

} // namespace Opm::Properties

namespace Opm::Parameters
{
} // namespace Opm::Parameters

#endif // OPM_MECHTYPETAG_HPP_INCLUDED
