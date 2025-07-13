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

#ifndef OPM_BLACKOILMODELGEOMECH_HPP_INCLUDED
#define OPM_BLACKOILMODELGEOMECH_HPP_INCLUDED

#include <opm/simulators/flow/BlackoilModel.hpp>

#include <cassert>
#include <iostream>
#include <string>

namespace Opm
{
template <class TypeTag>
class BlackoilModelGeomech : public BlackoilModel<TypeTag>
{
public:
    using Parent = BlackoilModel<TypeTag>;
    using Simulator = typename Parent::Simulator;
    using Scalar = typename Parent::Scalar;
    using ModelParameters = typename Parent::ModelParameters;

    BlackoilModelGeomech(Simulator& simulator,
                         const ModelParameters& param,
                         BlackoilWellModel<TypeTag>& well_model,
                         const bool terminal_output)
        : Parent {simulator, param, well_model, terminal_output}
    {
    }

    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIteration(const int iteration,
                                             const SimulatorTimerInterface& timer,
                                             NonlinearSolverType& nonlinear_solver)
    {
        const auto method
            = this->simulator_.problem().getGeomechParam().template get<std::string>("solver.method");

        if (method == "PostSolve") {
            return Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
        } else if (method == "SeqMechFrac") {
            return this->nonlinearIterationSeqMechFrac(iteration, timer, nonlinear_solver);
        } else if (method == "SeqMech") {
            return this->nonlinearIterationSeqMech(iteration, timer, nonlinear_solver);
        } else if (method == "FullyImplicitMech") {
            assert(false);
        } else {
            assert(false);

            std::cout << "Geomech nonlinearIterationNewton with mechanical solve:" << iteration
                      << std::endl;

            Parent::nonlinearIterationNewton(iteration, timer, nonlinear_solver);
        }

        return {};
    }

    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationSeqMechFrac(const int iteration,
                                                        const SimulatorTimerInterface& timer,
                                                        NonlinearSolverType& nonlinear_solver)
    {
        SimulatorReportSingle report {};

        const auto& prm = this->simulator_.problem().getGeomechParam();

        if (prm.template get<bool>("solver.implicit_flow")) {
            assert(false);
        } else {
            report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
        }

        if (iteration < prm.template get<int>("solver.max_mech_it")) {
            this->simulator_.problem().geomechModel().solveGeomechanics();

            if (this->simulator_.problem().hasFractures()) {
                this->simulator_.problem().geomechModel().solveFractures();
            }

            std::cout << "Geomech nonlinearIteration with mechanical and fracture solve:" << iteration
                      << std::endl;

            if (prm.template get<bool>("fractureparam.addconnections")) {
                std::cout << "Add connections in iterations" << std::endl;
                this->simulator_.problem().addConnectionsToSchedual();
                this->simulator_.problem().wellModel().beginTimeStep();
                this->simulator_.problem().addConnectionsToWell();
            }

            // TODO check convergence properly
            report.converged = false;
        }

        return report;
    }

    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationSeqMech(const int iteration,
                                                    const SimulatorTimerInterface& timer,
                                                    NonlinearSolverType& nonlinear_solver)
    {
        SimulatorReportSingle report {};

        const auto& prm = this->simulator_.problem().getGeomechParam();

        const auto implicit_flow = prm.template get<bool>("method.implicit_flow");

        if (implicit_flow) {
            assert(false);
        } else {
            report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
        }

        if (iteration < prm.template get<int>("method.max_mech_it")) {
            this->simulator_.problem().geomechModel().solveGeomechanics();

            std::cout << "Geomech nonlinearIteration with mechanical solve:" << iteration << std::endl;

            // TODO check convergence properly
            report.converged = false;
        }

        return report;
    }
};

} // namespace Opm

#endif // OPM_BLACKOILMODELGEOMECH_HPP_INCLUDED
