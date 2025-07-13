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

#ifndef OPM_GEOMECH_BLACKOILGEOMECHWELLMODEL_HPP_INCLUDED
#define OPM_GEOMECH_BLACKOILGEOMECHWELLMODEL_HPP_INCLUDED

#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <stdexcept>
#include <vector>

namespace Opm
{
template <typename TypeTag>
class BlackoilGeomechWellModel : public BlackoilWellModel<TypeTag>
{
    using Parent = BlackoilWellModel<TypeTag>;
    using Simulator = typename Parent::Simulator;

public:
    BlackoilGeomechWellModel(Simulator& simulator)
        : Parent(simulator)
    {
    }

    using NeighborSet = typename Parent::NeighborSet;

    void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        if (!this->param_.matrix_add_well_contributions_) {
            return;
        }

        OPM_THROW(std::runtime_error, "Not implemented");
    }

    void createWellContainer(const int reportStepIdx)
    {
        Parent::createWellContainer(reportStepIdx);

        // only add effect of fracture after one report step
        // NB everything is not explicit and ministeps are not considered
        if (reportStepIdx > 0) {
            const auto& problem = this->simulator_.problem();
            const auto& geomechmodel = problem.geomechModel();

            if (problem.hasFractures() && geomechmodel.fractureModelActive()) {
                for (auto& wellPtr : this->well_container_) {
                    const auto& fracturemodel = geomechmodel.fractureModel();

                    auto wellIndices = fracturemodel.getExtraWellIndices(wellPtr->name());
                    wellPtr->addPerforations(wellIndices);
                }
            }
        }
    }
};

} // namespace Opm

#endif // OPM_GEOMECH_BLACKOILGEOMECHWELLMODEL_HPP_INCLUDED
