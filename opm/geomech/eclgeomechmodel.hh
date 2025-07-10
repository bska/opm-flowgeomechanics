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

#ifndef OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#define OPM_ECLPROBLEM_GEOMECH_MODEL_HH

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/input/eclipse/Schedule/BCProp.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>

#include <opm/geomech/FlowGeomechLinearSolverParameters.hpp>
#include <opm/geomech/FractureModel.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/geomech/vem_elasticity_solver.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Opm
{
template <typename TypeTag>
class EclGeoMechModel : public BaseAuxiliaryModule<TypeTag>
{
    using Parent = BaseAuxiliaryModule<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using NeighborSet = typename BaseAuxiliaryModule<TypeTag>::NeighborSet;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    using Toolbox = MathToolbox<Evaluation>;
    using SymTensor = Dune::FieldVector<double, 6>;

public:
    EclGeoMechModel(Simulator& simulator)
        : first_solve_(true)
        , write_system_(false)
        , reduce_boundary_(false)
        , simulator_(simulator)
        , elacticitysolver_(simulator.vanguard().grid())
    {
    }

    // ax model things
    void postSolve(GlobalEqVector&)
    {
        std::cout << "Geomech dummy PostSolve Aux" << std::endl;
    }

    void addNeighbors(std::vector<NeighborSet>&) const
    {
        std::cout << "Geomech add neigbors" << std::endl;
    }

    void applyInitial()
    {
        std::cout << "Geomech applyInitial" << std::endl;
    }

    unsigned numDofs() const
    {
        return 0;
    }

    void linearize(SparseMatrixAdapter&, GlobalEqVector&)
    {
        std::cout << "Geomech Dummy Linearize" << std::endl;
    }

    // model things
    void beginIteration()
    {
        // Parent::beginIteration();
        std::cout << "Geomech begin iteration" << std::endl;
    }

    void endIteration()
    {
        // Parent::endIteration();
        std::cout << "Geomech end iteration" << std::endl;
    }

    void beginTimeStep()
    {
        // Parent::beginIteration();
        std::cout << "Geomech begin iteration" << std::endl;
    }

    void endTimeStep()
    {
        // always do post solve
        std::cout << "Geomech model endstimeStep" << std::endl;
        this->solveGeomechAndFracture();
    }

    void solveGeomechAndFracture()
    {
        // Parent::endIteration();
        const auto& problem = simulator_.problem();

        this->solveGeomechanics();

        if (problem.hasFractures()) {
            this->solveFractures();
        }
    }

    void solveFractures()
    {
        OPM_TIMEBLOCK(solveFractures);

        const int reportStepIdx = simulator_.episodeIndex();
        const auto& schedule = this->simulator_.vanguard().schedule();
        const int end_step = schedule.size() - 1;
        const bool no_seeds = schedule[end_step].wseed().empty();

        if (!no_seeds) {
            std::cout << "Fracture seeds found, on this step " << std::endl;
        } else {
            std::cout << "No fracture seeds found, on this step " << std::endl;
        }

        if (fracturemodel_) {
            std::cout << "Fracture model already initialized, solving fractures using "
                         "previous fractures"
                      << std::endl;
        }

        if (!no_seeds && !fracturemodel_) {
            std::cout << "Fracture model not initialized, initializing now. report step" << reportStepIdx
                      << std::endl;

            const auto& problem = simulator_.problem();

            // NB could probably be moved to some initialization
            //  let fracture contain all wells
            auto param = problem.getFractureParam();

            this->include_fracture_contributions_
                = param.template get<bool>("include_fracture_contributions");

            const auto& wells = problem.wellModel().getLocalWells(end_step);
            const auto& grid = simulator_.vanguard().grid();

            param.put("outputdir", Parameters::Get<Parameters::OutputDir>());
            param.put("casename", this->simulator_.vanguard().caseName());

            fracturemodel_ = std::make_unique<FractureModel>(grid, wells, param);

            // not to get the reservoir properties along the well before initialising the
            // well most important stress
            fracturemodel_->updateReservoirWellProperties<TypeTag>(simulator_);

            // add fractures along the wells
            // fracturemodel_->addFractures(schedule[reportStepIdx]);
            fracturemodel_->addFractures(schedule[end_step]);

            fracturemodel_->updateFractureReservoirCells();
            fracturemodel_->initReservoirProperties<TypeTag>(simulator_);
            fracturemodel_->updateReservoirAndWellProperties<TypeTag>(simulator_);
            fracturemodel_->initFractureStates();
        }

        // get reservoir properties on fractures
        // simulator need
        if (fracturemodel_) {
            std::cout << "Frac modelfound, updating reservoir properties and solving fractures"
                      << std::endl;
            fracturemodel_->updateReservoirAndWellProperties<TypeTag>(simulator_);
            fracturemodel_->solve<TypeTag>(simulator_);
        } else {
            std::cout << "Fracture model not initialized, not solving fractures" << std::endl;
        }
    }

    void writeFractureSolution()
    {
        const auto& problem = simulator_.problem();
        if (problem.hasFractures() && fracturemodel_) {
            // write first solution in standard format
            // this may ad some extra output of static variables
            int reportStepIdx = simulator_.episodeIndex();
            if (reportStepIdx == 1) {
                // fracturemodel_->write(reportStepIdx);
                //  hack to get correct number of fracture output
                fracturemodel_->writemulti(0.0);
            }
            double time = simulator_.time();
            fracturemodel_->writemulti(time);
        }
    }

    std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname)
    {
        if (fracturemodel_) {
            return fracturemodel_->getExtraWellIndices(wellname);
        } else {
            return std::vector<RuntimePerforation>();
        }
    }

    void updatePotentialForces()
    {
        std::cout << "Update Forces" << std::endl;
        const std::size_t numDof = simulator_.model().numGridDof();
        const auto& problem = simulator_.problem();

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            const auto& iq = simulator_.model().intensiveQuantities(dofIdx, 0);

            // pressure part
            const auto& fs = iq.fluidState();
            const auto& press = fs.pressure(waterPhaseIdx);

            // const auto& biotcoef = problem.biotCoef(dofIdx); //NB not used
            // thermal part
            // Properties::EnableTemperature
            const auto& poelCoef = problem.poelCoef(dofIdx);
            const double pressure = Toolbox::value(press);
            const double diffpress = (pressure - problem.initPressure(dofIdx));
            const auto& pratio = problem.pRatio(dofIdx);
            const double fac = (1 - pratio) / (1 - 2 * pratio);
            const double pcoeff = poelCoef * fac;

            assert(pcoeff <= 1.0);

            // assert pcoeff == biot
            pressure_[dofIdx] = pressure;
            mechPotentialForce_[dofIdx] = diffpress * pcoeff;
            mechPotentialPressForce_[dofIdx] = diffpress * pcoeff;
            mechPotentialPressForceFracture_[dofIdx] = diffpress * (1.0 - pcoeff);

            const bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();
            if (thermal_expansion) {
                OPM_TIMEBLOCK(addTermalParametersToMech);

                const auto& temp = fs.temperature(waterPhaseIdx); // NB all phases have equal temperature
                const auto& thelcoef = problem.thelCoef(dofIdx);

                // const auto& termExpr = problem.termExpr(dofIdx); //NB not used
                // tcoeff = (youngs*tempExp/(1-pratio))*fac;
                const double tcoeff = thelcoef * fac;

                // assume difftemp = 0 for non termal runs
                const double difftemp = Toolbox::value(temp) - problem.initTemperature(dofIdx);

                mechPotentialForce_[dofIdx] += difftemp * tcoeff;
                mechPotentialTempForce_[dofIdx] = difftemp * tcoeff;
            }

            // NB check sign !!
            mechPotentialForce_[dofIdx] *= 1.0;
        }
    }

    void setupMechSolver()
    {
        OPM_TIMEBLOCK(SetupMechSolver);

        const auto& problem = simulator_.problem();
        const auto& param = problem.getFractureParam();

        reduce_boundary_ = param.template get<bool>("reduce_boundary");
        const bool do_matrix = true; // assemble matrix
        const bool do_vector = true; // assemble matrix

        // set boundary
        elacticitysolver_.setBodyForce(0.0);
        elacticitysolver_.fixNodes(problem.bcNodes());
        elacticitysolver_.initForAssembly();
        elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector, reduce_boundary_);

        FlowLinearSolverParametersGeoMech p;
        p.init<TypeTag>();

        auto prm = setupPropertyTree(p, true, true);
        if (p.linear_solver_print_json_definition_ && (simulator_.gridView().comm().rank() == 0)) {
            std::ostringstream os;
            os << "Property tree for mech linear solver:\n";

            prm.write_json(os, true);

            OpmLog::note(os.str());
        }

        elacticitysolver_.setupSolver(prm);
        elacticitysolver_.comm()->communicator().barrier();

        first_solve_ = false;
        write_system_ = prm.get<int>("verbosity") > 10;
    }

    void writeMechSystem()
    {
        OPM_TIMEBLOCK(WriteMechSystem);

        const auto& problem = simulator_.problem();
        Helper::writeMechSystem(simulator_,
                                elacticitysolver_.A.getOperator(),
                                elacticitysolver_.A.getLoadVector(),
                                elacticitysolver_.comm());

        const int num_points = simulator_.vanguard().grid().size(3);

        Dune::BlockVector<Dune::FieldVector<double, 1>> fixed(3 * num_points);
        fixed = 0.0;

        const auto& bcnodes = problem.bcNodes();
        for (const auto& [node_idx, bcnode] : bcnodes) {
            for (int i = 0; i < 3; ++i) {
                fixed[3 * node_idx + i][0] = bcnode.fixeddir[i];
            }
        }

        Helper::writeVector(simulator_, fixed, "fixed_values_", elacticitysolver_.comm());
    }

    void calculateOutputQuantitiesMech()
    {
        OPM_TIMEBLOCK(CalculateOutputQuantitesMech);

        const auto& grid = simulator_.vanguard().grid();
        const auto& gv = grid.leafGridView();
        static constexpr int dim = Grid::dimension;

        Elasticity::Vector field;
        field.resize(grid.size(dim) * dim);

        if (reduce_boundary_) {
            elacticitysolver_.expandSolution(field, elacticitysolver_.u);
        } else {
            assert(field.size() == elacticitysolver_.u.size());
            field = elacticitysolver_.u;
        }

        this->makeDisplacement(field);

        // update variables used for output to resinsight
        // NB TO DO
        {
            OPM_TIMEBLOCK(calculateStress);
            elacticitysolver_.calculateStressPrecomputed(field);
            elacticitysolver_.calculateStrainPrecomputed(field);
        }

        const auto& linstress = elacticitysolver_.stress();
        const auto& linstrain = elacticitysolver_.strain();

        for (const auto& cell : elements(gv)) {
            const auto cellindex = simulator_.problem().elementMapper().index(cell);

            assert(cellindex == gv.indexSet().index(cell));

            strain_[cellindex] = linstrain[cellindex];
            linstress_[cellindex] = linstress[cellindex];
        }

        const bool verbose = false;
        if (verbose) {
            OPM_TIMEBLOCK(WriteMatrixMarket);

            // debug output to matrixmaket format
            Dune::storeMatrixMarket(elacticitysolver_.A.getOperator(), "A.mtx");
            Dune::storeMatrixMarket(elacticitysolver_.A.getLoadVector(), "b.mtx");
            Dune::storeMatrixMarket(elacticitysolver_.u, "u.mtx");
            Dune::storeMatrixMarket(field, "field.mtx");
            Dune::storeMatrixMarket(mechPotentialForce_, "pressforce.mtx");
        }
    }

    void setupAndUpdateGemechanics()
    {
        OPM_TIMEBLOCK(endTimeStepMech);

        this->updatePotentialForces();

        // for now assemble and set up solver her
        const auto& problem = simulator_.problem();

        if (first_solve_) {
            this->setupMechSolver();
        }

        {
            // reset the rhs even in the first iteration maybe bug in rhs
            // for reduce_boundary=false;
            OPM_TIMEBLOCK(AssembleRhs);

            elacticitysolver_.updateRhsWithGrad(mechPotentialForce_);
        }
    }

    void solveGeomechanics()
    {
        setupAndUpdateGemechanics();

        {
            OPM_TIMEBLOCK(SolveMechanicalSystem);
            elacticitysolver_.solve();

            if (write_system_) {
                this->writeMechSystem();
            }
        }

        this->calculateOutputQuantitiesMech(); // and properties used for fracturing
    }

    template <class Serializer>
    void serializeOp(Serializer&)
    {
        // serializer(tracerConcentration_);
        // serializer(wellTracerRate_);
    }

    // used in eclproblemgeomech
    void init(bool /*restart*/)
    {
        std::cout << "Geomech init" << std::endl;
        const std::size_t numDof = simulator_.model().numGridDof();

        pressure_.resize(numDof);
        mechPotentialForce_.resize(numDof);
        mechPotentialTempForce_.resize(numDof);
        mechPotentialPressForce_.resize(numDof);
        mechPotentialPressForceFracture_.resize(numDof);

        // hopefully temperature and pressure initilized
        celldisplacement_.resize(numDof);

        std::fill(celldisplacement_.begin(), celldisplacement_.end(), 0.0);

        // stress_.resize(numDof);
        linstress_.resize(numDof);
        std::fill(linstress_.begin(), linstress_.end(), 0.0);

        strain_.resize(numDof);
        std::fill(strain_.begin(), strain_.end(), 0.0);

        const auto& gv = simulator_.vanguard().grid().leafGridView();
        displacement_.resize(gv.indexSet().size(3));
    }

    void setMaterial(const std::vector<std::shared_ptr<Opm::Elasticity::Material>>& materials)
    {
        elacticitysolver_.setMaterial(materials);
    }

    void setMaterial(const std::vector<double>& ymodule, const std::vector<double>& pratio)
    {
        elacticitysolver_.setMaterial(ymodule, pratio);
    }

    const Dune::FieldVector<double, 3>& displacement(const std::size_t vertexIndex) const
    {
        return displacement_[vertexIndex];
    }

    double mechPotentialForce(unsigned globalDofIdx) const
    {
        return mechPotentialForce_[globalDofIdx];
    }

    double pressure(unsigned globalDofIdx) const
    {
        return pressure_[globalDofIdx];
    }

    double mechPotentialTempForce(unsigned globalDofIdx) const
    {
        return mechPotentialTempForce_[globalDofIdx];
    }

    double mechPotentialPressForce(unsigned globalDofIdx) const
    {
        return mechPotentialPressForce_[globalDofIdx];
    }

    Dune::FieldVector<double, 3> disp(const std::size_t globalIdx,
                                      const bool with_fracture = false) const
    {
        auto disp = celldisplacement_[globalIdx];

        if (include_fracture_contributions_ && with_fracture && (fracturemodel_ != nullptr)) {
            for (const auto& elem : Dune::elements(simulator_.vanguard().grid().leafGridView())) {
                const auto& center = elem.geometry().center();

                const Dune::FieldVector<double, 3> obs {center[0], center[1], center[2]};

                // check if this is correct stress
                disp += fracturemodel_->disp(obs);
            }
        }

        return disp;
    }

    SymTensor delstress(const std::size_t globalIdx) const
    {
        auto delStress = this->effstress(globalIdx);

        const double effPress = this->mechPotentialForce(globalIdx);
        for (int i = 0; i < 3; ++i) {
            delStress[i] += effPress;
        }

        return delStress;
    }

    const SymTensor& linstress(const std::size_t globalIdx) const
    {
        return linstress_[globalIdx];
    }

    SymTensor effstress(const std::size_t globalIdx) const
    {
        // make stress in with positive with compression
        return -1.0 * linstress_[globalIdx];
    }

    SymTensor strain(std::size_t globalIdx, bool with_fracture = false) const
    {
        auto strain = strain_[globalIdx];

        if (include_fracture_contributions_ && with_fracture && (fracturemodel_ != nullptr)) {
            for (const auto& elem : Dune::elements(simulator_.vanguard().grid().leafGridView())) {
                const auto center = elem.geometry().center();
                const Dune::FieldVector<double, 3> obs {center[0], center[1], center[2]};

                // check if this is correct stress
                strain += fracturemodel_->strain(obs);
            }
        }

        return strain_[globalIdx];
    }

    SymTensor stress(const std::size_t globalIdx, const bool with_fracture = false) const
    {
        auto effStress = this->effstress(globalIdx);
        effStress += simulator_.problem().initStress(globalIdx);

        const double effPress = this->mechPotentialForce(globalIdx);
        for (int i = 0; i < 3; ++i) {
            effStress[i] += effPress;
        }

        if (include_fracture_contributions_ && with_fracture && (fracturemodel_ != nullptr)) {
            for (const auto& elem : Dune::elements(simulator_.vanguard().grid().leafGridView())) {
                const auto center = elem.geometry().center();
                const Dune::FieldVector<double, 3> obs {center[0], center[1], center[2]};

                // check if this is correct stress
                effStress += fracturemodel_->stress(obs);
            }
        }

        return effStress;
    }

    SymTensor fractureStress(const std::size_t globalIdx) const
    {
        const auto& iq = simulator_.model().intensiveQuantities(globalIdx, 0);
        const auto& fs = iq.fluidState();
        const auto& press = Toolbox::value(fs.pressure(waterPhaseIdx));

        auto fracStress = this->stress(globalIdx);

        for (int i = 0; i < 3; ++i) {
            fracStress[i] -= press;
        }

        return fracStress;
    }

    // NB used in output should be eliminated

    double pressureDiff(const unsigned dofIx) const
    {
        return mechPotentialForce_[dofIx];
    }

    void makeDisplacement(const Opm::Elasticity::Vector& field)
    {
        // make displacement on all nodes used for output to vtk
        const auto& grid = simulator_.vanguard().grid();
        const auto& gv = grid.leafGridView();

        const int dim = 3;
        for (const auto& vertex : Dune::vertices(gv)) {
            const auto index = gv.indexSet().index(vertex);
            for (int k = 0; k < dim; ++k) {
                displacement_[index][k] = field[index * dim + k];
            }
        }

        for (const auto& cell : elements(gv)) {
            const auto cellindex = simulator_.problem().elementMapper().index(cell);
            assert(cellindex == gv.indexSet().index(cell));

            celldisplacement_[cellindex] = 0.0;

            const auto& vertices = Dune::subEntities(cell, Dune::Codim<Grid::dimension> {});
            for (const auto& vertex : vertices) {
                const auto nodeidex = gv.indexSet().index(vertex);

                for (int k = 0; k < dim; ++k) {
                    celldisplacement_[cellindex][k] += displacement_[nodeidex][k];
                }
            }

            celldisplacement_[cellindex] /= vertices.size();
        }
    }

    bool fractureModelActive() const
    {
        return this->fracturemodel_ != nullptr;
    }

    const FractureModel& fractureModel() const
    {
        if (!fracturemodel_) {
            std::cout << "Fracture model not initialized, returning nullptr" << std::endl;
            throw std::runtime_error("Fracture model not initialized");
        }

        return *fracturemodel_;
    }

private:
    bool first_solve_ {true};
    bool write_system_ {false};
    bool reduce_boundary_ {false};
    bool include_fracture_contributions_ {false};
    Simulator& simulator_;

    Dune::BlockVector<Dune::FieldVector<double, 1>> pressure_;
    Dune::BlockVector<Dune::FieldVector<double, 1>> mechPotentialForce_;
    Dune::BlockVector<Dune::FieldVector<double, 1>> mechPotentialPressForce_;
    Dune::BlockVector<Dune::FieldVector<double, 1>> mechPotentialPressForceFracture_;
    Dune::BlockVector<Dune::FieldVector<double, 1>> mechPotentialTempForce_;
    Dune::BlockVector<Dune::FieldVector<double, 3>> celldisplacement_;
    Dune::BlockVector<Dune::FieldVector<double, 3>> displacement_;
    Dune::BlockVector<Dune::FieldVector<double, 6>> linstress_; // NB is also stored in esolver
    Dune::BlockVector<Dune::FieldVector<double, 6>> strain_;
    Opm::Elasticity::VemElasticitySolver<Grid> elacticitysolver_;

    std::unique_ptr<FractureModel> fracturemodel_;
};

} // namespace Opm

#endif // OPM_ECLPROBLEM_GEOMECH_MODEL_HH
