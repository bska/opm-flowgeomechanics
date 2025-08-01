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

#ifndef OPM_ECLPROBLEM_GEOMECH_HH
#define OPM_ECLPROBLEM_GEOMECH_HH

#include <opm/common/ErrorMacros.hpp>

#include <opm/elasticity/material.hh>
#include <opm/elasticity/materials.hh>

#include <opm/geomech/FlowGeomechLinearSolverParameters.hpp>
#include <opm/geomech/boundaryutils.hh>
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/vtkgeomechmodule.hh>

#include <opm/elasticity/material.hh>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Parameters
{
struct FractureParamFile
{
    inline static std::string value {"notafile"};
};
} // namespace Opm::Parameters

namespace Opm
{
template <typename TypeTag>
class EclProblemGeoMech : public FlowProblemBlackoil<TypeTag>
{
public:
    using Parent = FlowProblemBlackoil<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using TimeStepper = AdaptiveTimeStepping<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using Toolbox = MathToolbox<Evaluation>;
    using SymTensor = Dune::FieldVector<double, 6>;
    using GeomechModel = EclGeoMechModel<TypeTag>;

    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    EclProblemGeoMech(Simulator& simulator)
        : FlowProblemBlackoil<TypeTag>(simulator)
        , geomechModel_(simulator)
    {
        try {
            // set seed values
            fracture_param_ = PropertyTree {Parameters::Get<Parameters::FractureParamFile>()};
        } catch (...) {
            std::stringstream ss;
            ss << "No fracture parameter file: " << Parameters::Get<Parameters::FractureParamFile>()
               << " : no fractures added ";

            OpmLog::warning(ss.str());

            fracture_param_ = makeDefaultFractureParam();
            fracture_param_.put("hasfractures", false);
        }

        fracture_param_.write_json(std::cout, true);

        hasFractures_ = this->simulator()
                            .vanguard()
                            .eclState()
                            .runspec()
                            .frac(); // fracture_param_.get<bool>("hasfractures");

        addPerfsToSchedule_ = fracture_param_.get<bool>("add_perfs_to_schedule");

        if (this->simulator().vanguard().eclState().runspec().mech()) {
            this->model().addOutputModule(std::make_unique<VtkGeoMechModule<TypeTag>>(simulator));
        }
    }

    static void registerParameters()
    {
        Parent::registerParameters();
        VtkGeoMechModule<TypeTag>::registerParameters();
        FlowLinearSolverParametersGeoMech::registerParameters<TypeTag>();
        Parameters::Register<Parameters::FractureParamFile>("json file defining fracture setting");

        Parameters::SetDefault<Opm::Parameters::EnableOpmRstFile>(true);
        Parameters::SetDefault<Opm::Parameters::EnableVtkOutput>(true);
        Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
        Parameters::SetDefault<Opm::Parameters::EnableAsyncVtkOutput>(false);
        Parameters::SetDefault<Opm::Parameters::EnableAsyncEclOutput>(false);
    }

    void finishInit()
    {
        OPM_TIMEBLOCK(finishInit);

        Parent::finishInit();

        const auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();

        if (eclState.runspec().mech()) {
            const auto& initconfig = eclState.getInitConfig();
            geomechModel_.init(initconfig.restartRequested());
            const auto& fp = eclState.fieldProps();

            for (const auto* key : {"YMODULE", "PRATIO"}) {
                if (!fp.has_double(key)) {
                    std::stringstream ss;
                    ss << "Missing keyword '" << key << '\'';

                    OPM_THROW(std::runtime_error, ss.str());
                }
            }

            ymodule_ = fp.get_double("YMODULE");
            pratio_ = fp.get_double("PRATIO");

            if (fp.has_double("BIOTCOEF")) {
                biotcoef_ = fp.get_double("BIOTCOEF");
                poelcoef_.resize(ymodule_.size());

                for (int i = 0; i < ymodule_.size(); ++i) {
                    poelcoef_[i] = (1 - 2 * pratio_[i]) / (1 - pratio_[i]) * biotcoef_[i];
                }
            } else {
                if (!fp.has_double("POELCOEF")) {
                    OPM_THROW(std::runtime_error, "Missing keyword BIOTCOEF or POELCOEF");
                }

                poelcoef_ = fp.get_double("POELCOEF");
                biotcoef_.resize(ymodule_.size());

                for (int i = 0; i < ymodule_.size(); ++i) {
                    biotcoef_[i] = poelcoef_[i] * (1 - pratio_[i]) / (1 - 2 * pratio_[i]);
                }
            }

            // thermal related
            if (getPropValue<TypeTag, Properties::EnableEnergy>()) {
                if (fp.has_double("THELCOEF")) {
                    thelcoef_ = fp.get_double("THELCOEF");
                    thermexr_.resize(ymodule_.size());

                    for (int i = 0; i < ymodule_.size(); ++i) {
                        thermexr_[i] = thelcoef_[i] * (1 - pratio_[i]) / ymodule_[i];
                    }
                } else {
                    if (!fp.has_double("THERMEXR")) {
                        OPM_THROW(std::runtime_error, "Missing keyword THELCOEF or THERMEXR");
                    }

                    thermexr_ = fp.get_double("THERMEXR");
                    thelcoef_.resize(ymodule_.size());

                    for (int i = 0; i < ymodule_.size(); ++i) {
                        thelcoef_[i] = thermexr_[i] * ymodule_[i] / (1 - pratio_[i]);
                    }
                }
            }

            for (std::size_t i = 0; i < ymodule_.size(); ++i) {
                using IsoMat = Elasticity::Isotropic;

                if (pratio_[i] > 0.5 || pratio_[i] < 0.0) {
                    OPM_THROW(std::runtime_error, "Pratio not valid");
                }

                if (biotcoef_[i] > 1.0 || biotcoef_[i] < 0.0) {
                    OPM_THROW(std::runtime_error, "BIOTCOEF not valid");
                }

                elasticparams_.push_back(std::make_shared<IsoMat>(i, ymodule_[i], pratio_[i]));
            }

            if (fp.has_double("CSTRESS")) {
                cstress_ = fp.get_double("CSTRESS");
            }

            // read mechanical boundary conditions
            // const auto& simulator = this->simulator();
            const auto& vanguard = simulator.vanguard();
            const auto& bcconfigs = vanguard.eclState().getSimulationConfig().bcconfig();
            const auto& bcprops = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;

            // using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
            const auto& gv = this->gridView();

            // const auto& grid = simulator.grid();
            const auto& cartesianIndexMapper = vanguard.cartesianIndexMapper();

            // CartesianIndexMapper cartesianIndexMapper(grid);
            Elasticity::nodesAtBoundary(bc_nodes_, bcconfigs, bcprops, gv, cartesianIndexMapper);

            // using Opm::ParserKeywords::;
            if (initconfig.hasStressEquil()) {
                const std::size_t numCartDof = cartesianIndexMapper.cartesianSize();
                const unsigned numElems = gv.size(/*codim=*/0);

                std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);
                for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
                    cartesianToCompressedElemIdx[cartesianIndexMapper.cartesianIndex(elemIdx)] = elemIdx;
                }

                const auto& stressequil = initconfig.getStressEquil();
                const auto& equilRegionData = fp.get_int("STRESSEQUILNUM");

                // make lambda functions for each regaion
                std::vector<std::function<std::array<double, 6>()>> functors;
                int recnum = 1;
                initstress_.resize(gv.size(0));
                for (const auto& record : stressequil) {
                    const auto datum_depth = record.datumDepth();
                    const auto STRESSXX = record.stressXX();
                    const auto STRESSXXGRAD = record.stressXX_grad();
                    const auto STRESSYY = record.stressYY();
                    const auto STRESSYYGRAD = record.stressYY_grad();
                    const auto STRESSZZ = record.stressZZ();
                    const auto STRESSZZGRAD = record.stressZZ_grad();
                    const auto STRESSXY = record.stressXY();
                    const auto STRESSXZ = record.stressXZ();
                    const auto STRESSYZ = record.stressYZ();

                    const auto stressXYGRAD = record.stressXY_grad();
                    const auto stressXZGRAD = record.stressXZ_grad();
                    const auto stressYZGRAD = record.stressYZ_grad();
                    for (const auto& cell : elements(gv)) {
                        const auto& center = cell.geometry().center();
                        const auto& cellIdx = gv.indexSet().index(cell);
                        assert(cellIdx < equilRegionData.size());
                        const auto& region
                            = equilRegionData[cellIdx]; // cartesianIndexMapper.cartesianIndex(cellIdx)];
                        assert(region <= stressequil.size());
                        if (region == recnum) {
                            Dune::FieldVector<double, 6> initstress;
                            initstress[0] = STRESSXX + STRESSXXGRAD * (center[2] - datum_depth);
                            initstress[1] = STRESSYY + STRESSYYGRAD * (center[2] - datum_depth);
                            initstress[2] = STRESSZZ + STRESSZZGRAD * (center[2] - datum_depth);
                            initstress[3] = STRESSYZ + stressYZGRAD * (center[2] - datum_depth);
                            initstress[4] = STRESSXZ + stressXZGRAD * (center[2] - datum_depth);
                            initstress[5] = STRESSXY + stressXYGRAD * (center[2] - datum_depth);

                            initstress_[cellIdx] = initstress;
                        }
                    }

                    recnum += 1;
                }
            } else {
                OPM_THROW(std::runtime_error, "Missing stress initialization keywords");
            }
        }
    }

    void initialSolutionApplied()
    {
        OPM_TIMEBLOCK(initialSolutionApplied);

        Parent::initialSolutionApplied();

        const auto& simulator = this->simulator();
        const std::size_t numDof = simulator.model().numGridDof();

        initpressure_.resize(numDof);
        inittemperature_.resize(numDof);

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            const auto& iq = this->model().intensiveQuantities(dofIdx, 0);
            const auto& fs = iq.fluidState();
            initpressure_[dofIdx] = Toolbox::value(fs.pressure(waterPhaseIdx));
            inittemperature_[dofIdx] = Toolbox::value(fs.temperature(waterPhaseIdx));
        }

        initstress_.resize(numDof);

        // for now make a copy
        if (simulator.vanguard().eclState().runspec().mech()) {
            this->geomechModel_.setMaterial(ymodule_, pratio_);
            this->geomechModel_.updatePotentialForces();
        }
    }

    void timeIntegration()
    {
        if (this->gridView().comm().rank() == 0) {
            std::cout << "----------------------Start TimeIntegration-------------------\n"
                      << std::flush;
        }

        Parent::timeIntegration();
    }

    void beginTimeStep()
    {
        if (this->gridView().comm().rank() == 0) {
            std::cout << "----------------------Start beginTimeStep-------------------\n" << std::flush;
        }

        Parent::beginTimeStep();

        if (this->simulator().vanguard().eclState().runspec().mech()) {
            if (this->hasFractures()) {
                if (!(cstress_.size() == this->gridView().size(0))) {
                    OPM_THROW(std::runtime_error, "CSTRESS not set but fractures exists");
                }
            }

            geomechModel_.beginTimeStep();

            if (this->hasFractures()) {
                this->wellModel().beginTimeStep(); // just to be sure well conteiner is reinitialized
                this->addConnectionsToWell(); // modify wells WI wiht fracture well
            }
        }
    }

    void endTimeStep()
    {
        if (this->gridView().comm().rank() == 0) {
            std::cout << "----------------------Start endTimeStep-------------------\n" << std::flush;
        }

        // Parent::FlowProblemType::endTimeStep();
        if (this->simulator().vanguard().eclState().runspec().mech()) {
            geomechModel_.endTimeStep();
            if (this->hasFractures() && this->geomechModel().fractureModelActive()) {
                // method for handling extra connections from fractures
                // it is options for not including them in fractures i.e. addconnections
                if (addPerfsToSchedule_) {
                    this->addConnectionsToSchedual();
                } else {
                    // not not working ... more work...
                    // will only work if structure is ok
                    assert(false);
                    this->addConnectionsToWell();
                }

                this->geomechModel_.fractureModel().assignGeomechWellState(this->wellModel_.wellState());
            }
        }

        Parent::endTimeStep();
    }

    void addConnectionsToWell()
    {
        auto& wellcontainer = this->wellModel().localNonshutWells();

        for (auto& wellPtr : wellcontainer) {
            auto wellName = wellPtr->name();
            const auto& wellcons = geomechModel_.getExtraWellIndices(wellName);
            wellPtr->addPerforations(wellcons);
        }
    }

    void addConnectionsToSchedual()
    {
        auto& simulator = this->simulator();
        auto& schedule = simulator.vanguard().schedule();
        const int reportStep = this->episodeIndex();

        std::map<std::string, std::vector<Connection>> extra_perfs;

        for (const auto& wellName : schedule.wellNames(reportStep)) {
            const auto wellcons = geomechModel_.getExtraWellIndices(wellName);

            if (wellcons.empty()) {
                // No extra connections for this well.
                continue;
            }

            const auto& origConns = this->schedule_[reportStep].wells(wellName).getConnections();

            auto extra = std::vector<Connection> {};

            for (const auto& wellconn : wellcons) {
                // simple calculated with upscaling

                // map to cartesian
                const auto cartesianIdx = simulator.vanguard().cartesianIndex(wellconn.cell);

                if (origConns.hasGlobalIndex(cartesianIdx)) {
                    std::cout << "Connection already exists for cell: " << wellconn.cell << std::endl;
                    continue;
                }

                // get ijk
                std::array<int, 3> ijk {};
                simulator.vanguard().cartesianCoordinate(wellconn.cell, ijk);

                // Making preliminary connection to be added in schedule
                // with correct numbering
                auto& connection = extra.emplace_back(ijk[0],
                                                      ijk[1],
                                                      ijk[2],
                                                      cartesianIdx,
                                                      /*complnum*/ -1,
                                                      Connection::State::OPEN,
                                                      Connection::Direction::Z,
                                                      Connection::CTFKind::DynamicFracturing,
                                                      /* satTableId */ -1,
                                                      wellconn.depth,
                                                      Connection::CTFProperties {},
                                                      /* sort_value */ -1,
                                                      /* defaut sattable */ true);

                // only add zero value
                //  connection need to be modified later.
                connection.setCF(wellconn.ctf * 0.0);

                if (wellconn.perf_range.has_value()) {
                    const auto compseg_insert_index = std::numeric_limits<std::size_t>::max();

                    connection.updateSegment(
                        wellconn.segment, wellconn.depth, compseg_insert_index, wellconn.perf_range);
                }
            }

            if (!extra.empty()) {
                const auto* pl = (extra.size() != 1) ? "s" : "";

                std::cout << "Adding " << extra.size() << "extra connection" << pl
                          << " for well: " << wellName << std::endl;

                extra_perfs.insert_or_assign(wellName, std::move(extra));
            }
        }

        if (extra_perfs.empty()) {
            return;
        } else {
            // add to schedule
            // structure will be changed erase matrix, maybe only rebuilding of linear
            // solver is neede
            this->simulator().model().linearizer().eraseMatrix();
            if (this->gridView().comm().rank() == 0) {
                std::cout << "Adding extra connections to "
                             "schedule for report step: "
                          << reportStep << std::endl;
            }
        }

        bool commit_wellstate = false;
        auto sim_update = schedule.modifyCompletions(reportStep, extra_perfs);

        // should not be used
        auto updateTrans = [](const bool) {};

        // alwas rebuild wells
        sim_update.well_structure_changed = true;

        this->actionHandler_.applySimulatorUpdate(reportStep, sim_update, updateTrans, commit_wellstate);
        if (commit_wellstate) {
            this->wellModel().commitWGState();
        }
    }

    void endEpisode()
    {
        Parent::endEpisode();
        geomechModel_.writeFractureSolution();
    }

    const EclGeoMechModel<TypeTag>& geoMechModel() const
    {
        return geomechModel_;
    }

    EclGeoMechModel<TypeTag>& geoMechModel()
    {
        return geomechModel_;
    }

    double initPressure(const unsigned dofIdx) const
    {
        return initpressure_[dofIdx];
    }

    double initTemperature(const unsigned dofIdx) const
    {
        return inittemperature_[dofIdx];
    }

    double initStress(const unsigned dofIdx, const int comp) const
    {
        return initstress_[dofIdx][comp];
    }

    const SymTensor& initStress(const unsigned dofIdx) const
    {
        return initstress_[dofIdx];
    }

    double biotCoef(const unsigned globalIdx) const
    {
        return biotcoef_[globalIdx];
    }

    double thelCoef(const unsigned globalIdx) const
    {
        return thelcoef_[globalIdx];
    }

    double thermExr(const unsigned globalIdx) const
    {
        return thermexr_[globalIdx];
    }

    double poelCoef(const unsigned globalIdx) const
    {
        return poelcoef_[globalIdx];
    }

    const std::vector<std::tuple<std::size_t, MechBCValue>>& bcNodes() const
    {
        return bc_nodes_;
    }

    Dune::FieldVector<double, 6> stress(const std::size_t globalIdx) const
    {
        return geomechModel_.stress(globalIdx);
    }

    bool hasFractures() const
    {
        return hasFractures_;
    }

    PropertyTree getFractureParam() const
    {
        return fracture_param_.get_child("fractureparam");
    }

    PropertyTree getGeomechParam() const
    {
        return fracture_param_;
    }

    GeomechModel& geomechModel()
    {
        return geomechModel_;
    }

    const GeomechModel& geomechModel() const
    {
        return geomechModel_;
    }

    // used for fracture model
    double yModule(std::size_t idx) const
    {
        return ymodule_[idx];
    }

    double pRatio(std::size_t idx) const
    {
        return pratio_[idx];
    }

    double cStress(std::size_t idx) const
    {
        return cstress_[idx];
    }

private:
    GeomechModel geomechModel_;

    std::vector<double> ymodule_;
    std::vector<double> pratio_;
    std::vector<double> biotcoef_;
    std::vector<double> poelcoef_;
    std::vector<double> thermexr_;
    std::vector<double> thelcoef_;
    std::vector<double> cstress_;

    std::vector<double> initpressure_;
    std::vector<double> inittemperature_;
    std::vector<std::tuple<std::size_t, MechBCValue>> bc_nodes_;
    Dune::BlockVector<SymTensor> initstress_;
    std::vector<std::shared_ptr<Opm::Elasticity::Material>> elasticparams_;

    // for fracture calculation
    bool hasFractures_;
    bool addPerfsToSchedule_;
    PropertyTree fracture_param_;
};

} // namespace Opm

#endif // OPM_ECLPROBLEM_GEOMECH_HH
