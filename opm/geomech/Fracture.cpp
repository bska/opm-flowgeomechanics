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

#include <opm/geomech/Fracture.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/common/fmatrixev.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/istl/io.hh> // needed for printSparseMatrix??

#include <opm/grid/polyhedralgrid.hh>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/wells/ConnFracStatistics.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>

#include <opm/geomech/DiscreteDisplacement.hpp>
#include <opm/geomech/GridStretcher.hpp>
#include <opm/geomech/Math.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{
void
Fracture::init(const std::string& well,
               const int perf,
               const int well_cell,
               int global_index,
               const int segment,
               const std::optional<std::pair<double, double>>& perf_range,
               const Point3D& origo,
               const Point3D& normal,
               const Opm::PropertyTree& prm)
{
    OPM_TIMEFUNCTION();
    prm_ = prm;
    min_width_ = prm_.get<double>("config.min_width", 1e-3);
    wellinfo_ = WellInfo({well, perf, well_cell, global_index, segment, perf_range});

    origo_ = origo;
    axis_[2] = normal;
    axis_[0] = Point3D({std::copysign(normal[2], normal[0]),
                        std::copysign(normal[2], normal[1]),
                        -std::copysign(std::abs(normal[0]) + std::abs(normal[1]), normal[2])});
    axis_[1] = crossProduct(axis_[2], axis_[0]);
    double init_scale = prm_.get<double>("config.axis_scale");
    for (int i = 0; i < 3; ++i) {
        axis_[i] /= axis_[i].two_norm();
        axis_[i] *= init_scale;
    }

    std::cout << "axis: {" << axis_[0][0] << ',' << axis_[0][1] << ',' << axis_[0][2] << "}, { "
              << axis_[1][0] << ',' << axis_[1][1] << ',' << axis_[1][2] << "}, { " << axis_[2][0] << ','
              << axis_[2][1] << ',' << axis_[2][2] << '}' << std::endl;

    layers_ = 0;
    nlinear_ = 0;

    const auto method = prm_.get<std::string>("solver.method");

    if (method == "if_propagate_trimesh") {
        // const int trimeshlayers = 4;
        // const double init_scale = prm_.get<double>("config.axis_scale");
        double trires = prm_.get<double>("config.trires");
        const double edgelen = init_scale / trires; // 1;
        const double radius = init_scale;
        const double fac = std::sqrt(3) / 2;
        const std::array<double, 3> ax1 {axis_[0][0], axis_[0][1], axis_[0][2]};
        const std::array<double, 3> ax2 {0.5 * ax1[0] + fac * axis_[1][0],
                                         0.5 * ax1[1] + fac * axis_[1][1],
                                         0.5 * ax1[2] + fac * axis_[1][2]};

        std::cout << "Creating trimesh with radius: " << radius << ", edgelen: " << edgelen
                  << ", ax1: " << ax1[0] << "," << ax1[1] << "," << ax1[2] << ", ax2: " << ax2[0] << ","
                  << ax2[1] << "," << ax2[2] << std::endl;

        trimesh_ = std::make_unique<RegularTrimesh>(radius, // trimeshlayers,
                                                    std::array {origo_[0], origo_[1], origo_[2]},
                                                    ax1,
                                                    ax2,
                                                    std::array {edgelen, edgelen});

        trimesh_->removeSawtooths();

        // identify well cells (since this is not done in setFractureGrid when providing
        // a user-defined grid)
        // std::vector<CellRef> wellcells { {0, 0, 0}, {0, -1, 1}, {0, -1, 0}, {-1, -1,
        // 1},
        //                                  {-1, 0, 0}, {-1, 0, 1} };
        well_source_cellref_ = RegularTrimesh::inner_ring_cells();

        for (const auto& cell : well_source_cellref_) {
            well_source_.push_back(trimesh_->linearCellIndex(cell));
        }

        auto [grid, fsmap, bmap] = trimesh_->createDuneGrid(1, well_source_cellref_);
        grid_mesh_map_ = fsmap;

        setFractureGrid(std::move(grid)); // create the physical grid from trimesh
    } else {
        setFractureGrid();
    }

    setPerfPressure(0.0); // This can be changed by subsequently calling this
                          // function when the fracture is connected to the
                          // reservoir

    // NB: The Fracture object is still not fully initialized, since there is
    // not yet any connection with a surrounding reservoir. The following
    // needs to be done before the fracture can be used:
    // 1) Call `updateReservoirCells` to establish the mapping between frature grid
    //    cells and the cells in the reservoir grid (sets `reservoir_cells_`)
    // 2) Then, call `updateReservoirProperties` to import relevant reservoir properties
    //    to the fracture (`reservoir_XXX_` vectors, as well as `E_` and `nu_`)
}

void
Fracture::resetWriters()
{
    // nead to be reseat if grid is changed ??
    vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(),
                                                                       Dune::VTK::nonconforming);

    std::string outputdir = prm_.get<std::string>("outputdir");
    std::string simName = prm_.get<std::string>("casename") + this->name();
    std::string multiFileName = "";
    if (!vtkmultiwriter_) {
        vtkmultiwriter_ = std::make_unique<Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>>(
            /*async*/ false, grid_->leafGridView(), outputdir, simName, multiFileName);
    } else {
        vtkmultiwriter_->gridViewChanged(grid_->leafGridView()); // need to be called if grid is changed
    }
}

void
Fracture::setupPressureSolver()
{
    Opm::FlowLinearSolverParameters p;
    p.linsolver_ = prm_.get<std::string>("pressuresolver");
    prmpressure_ = Opm::setupPropertyTree(p, true, true);

    const std::size_t pressureIndex = 0; // Dummy
    const std::function<Vector()> weightsCalculator; // Dummy

    pressure_operator_ = std::make_unique<PressureOperatorType>(*pressure_matrix_);

    pressure_solver_ = std::make_unique<FlexibleSolverType>(
        *pressure_operator_, prmpressure_, weightsCalculator, pressureIndex);
}

/**
 * @brief Removes cells from the grid if they are out side reservoir.
 *
 * This function performs the following steps:
 * 1. Copies the current fracture width data to a persistent container.
 * 2. Iterates over all elements in the grid's leaf view and removes elements where the
 * reservoir cell index is negative.
 * 3. Grows the grid and performs post-growth operations.
 * 4. Resizes the fracture width array to match the new grid size.
 * 5. Copies the fracture width data back from the persistent container to the resized
 * array.
 * 6. Resets the writers associated with the Fracture object.
 */
void
Fracture::removeCells()
{
    // copy all to presistent container
    const ElementMapper mapper(grid_->leafGridView(),
                               Dune::mcmgElementLayout()); // used id sets interally

    Dune::PersistentContainer<Grid, double> fracture_width(*grid_, 0);
    fracture_width.resize();
    for (const auto& elem : elements(grid_->leafGridView())) {
        std::size_t eIdx = mapper.index(elem);
        fracture_width[elem] = fracture_width_[eIdx];
    }

    // const auto& indexSet = foamGridLeafView.indexSet();// for indices
    // const auto& indexSet = grid.localIdSet();// presitent numbering
    for (const auto& elem : elements(grid_->leafGridView())) {
        std::size_t eIdx = mapper.index(elem);
        if (reservoir_cells_[eIdx] < 0) {
            grid_->removeElement(elem);
        }
    }

    grid_->grow();
    grid_->postGrow();
    // resize array

    fracture_width_.resize(numFractureCells());
    // copy back from presistent contatiner
    for (const auto& elem : elements(grid_->leafGridView())) {
        std::size_t eIdx = mapper.index(elem);
        fracture_width_[eIdx] = fracture_width[elem];
    }
    this->resetWriters();
}

void
Fracture::updateFilterCakeProps(const Opm::WellConnections& connections,
                                const Opm::SingleWellState<double>& wellstate)
{
    OPM_TIMEFUNCTION();

    // potentially remap filtercake thikness if grid has changed
    assert(filtercake_thikness_.size() == reservoir_cells_.size());

    const auto& perfdata = wellstate.perf_data;
    std::map<int, double> WI_fluxes;
    int water_index = 0;
    int np = 3;
    for (int i = 0; i < perfdata.cell_index.size(); ++i) {
        int cell_index = perfdata.cell_index[i];
        // similar as in WellFilterCake.cpp:148
        const auto& connection_rates = perfdata.phase_rates;
        const double water_rate = std::max(0.0, connection_rates[i * np + water_index]);
        WI_fluxes[cell_index] = water_rate;
    }

    const auto& connection = connections.getFromGlobalIndex(wellinfo_.global_index); // probably wrong
    auto& filter_cake = connection.getFilterCake();
    has_filtercake_ = connection.filterCakeActive();
    if (has_filtercake_) {
        filtercake_poro_ = filter_cake.poro;
        filtercake_perm_ = filter_cake.perm;

        std::map<int, double> reservoir_areas;
        std::map<int, double> reservoir_flux;

        if (!leakof_.empty()) { // should be first step with seed fracture is clean could
                                // have used rates from solve
            assert(leakof_.size() == fracture_pressure_.size());
            const ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
            for (const auto& element : Dune::elements(grid_->leafGridView())) {
                const int eIdx = mapper.index(element);
                const int res_cell = reservoir_cells_[eIdx];
                const double area = element.geometry().volume(); // is the area of this face

                reservoir_areas[res_cell] += area;

                const double dp = fracture_pressure_[eIdx] - reservoir_pressure_[eIdx];
                double flux = leakof_[eIdx] * dp;

                reservoir_flux[res_cell] += flux;

                if (flux < 0) {
                    std::cout << "Negative flux " << flux << " for element index " << eIdx
                              << " with reservoir cell " << res_cell << std::endl;
                    flux = 0.0;
                }

                const double dh = flux / (area * (1 - filtercake_poro_));

                assert(dh >= 0);
                filtercake_thikness_[eIdx] += dh;
            }
        }

        /// for setting or controling the fluxes between multiphase and single phase
        const ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
        for (const auto& [res_cell, flux] : reservoir_flux) {
            if (reservoir_areas.find(res_cell) == reservoir_areas.end()) {
                std::cout << "Reservoir area not found for element index " << res_cell << std::endl;
                continue;
            }

            const double frac_flux = reservoir_flux[res_cell];
            if (res_cell != wellinfo_.well_cell) {
                if (std::abs(flux - WI_fluxes[res_cell]) > 0) {
                    std::cout << "Fracture flux differs from flow flux " << res_cell << '\n'
                              << "Flux: frac " << frac_flux << " vs res " << WI_fluxes[res_cell]
                              << std::endl;
                }
            } else {
                std::cout << "Total WI flux " << WI_fluxes[res_cell] << " for cell " << res_cell
                          << " matches fracture flux " << frac_flux
                          << std::endl; //" for element index " << eIdx << std::endl;
            }
        }
    }
}

Dune::BlockVector<Dune::FieldVector<double, 3>>
Fracture::all_slips() const
{
    Dune::BlockVector<Dune::FieldVector<double, 3>> slips(grid_->leafGridView().size(0));
    slips = 0;

    const ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());

    for (const auto& elem : Dune::elements(grid_->leafGridView())) {
        std::size_t eIdx = mapper.index(elem);
        // only normal slip for now
        slips[eIdx][0] = fracture_width_[eIdx];
    }

    return slips;
}

Dune::FieldVector<double, 3>
Fracture::disp(const Dune::FieldVector<double, 3>& obs) const
{
    return ddm::disp(obs, this->all_slips(), *grid_, E_, nu_);
}

Dune::FieldVector<double, 6>
Fracture::strain(const Dune::FieldVector<double, 3>& obs) const
{
    // for now use full slip in interface even we only calculate normal slip
    return ddm::strain(obs, this->all_slips(), *grid_, E_, nu_);
}

Dune::FieldVector<double, 6>
Fracture::stress(const Dune::FieldVector<double, 3>& obs) const
{
    // for now use full slip in interface even we only calculate normal slip
    return ddm::strainToStress(E_, nu_, this->strain(obs));
}

std::string
Fracture::name() const
{
    std::string name
        = "Fracure_on_" + wellinfo_.name + "_perf_" + std::to_string(wellinfo_.perf) + "_nr";
    return name;
}

void
Fracture::updateCellNormals()
{
    cell_normals_.resize(numFractureCells());

    const ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : elements(grid_->leafGridView())) {
        cell_normals_[elemMapper.index(element)] = ddm::normalOfElement(element);
    }
}

void
Fracture::setFractureGrid(std::unique_ptr<Fracture::Grid> gptr)
{
    if (gptr == nullptr) {
        // create a new grid
        initFracture();

        const int num_exp = prm_.get<int>("config.num_exp");
        if (num_exp > 0) {
            grow(num_exp, 0);
        }

        nlinear_ = layers_;

        const int num_lin = prm_.get<int>("config.num_lin");
        if (num_lin > 0) {
            grow(num_lin, 1);
        }

        grid_->grow();
        grid_->postGrow();
    } else {
        // use grid provided by user
        grid_ = std::move(gptr);
    }

    // compute the cell normals
    updateCellNormals();

    // set object to allow stretching of the grid // @@ do not use with trimesh
    if (trimesh_ == nullptr) {
        grid_stretcher_ = std::unique_ptr<GridStretcher>(new GridStretcher(*grid_));
    }

    // set sparsity of pressure matrix, but do not compute its entries
    initPressureMatrix();

    // Since the grid has been created/updated/changed, any previous mapping to
    // reservoir cells has been invalidated, and the fracture matrix (for
    // mechanics) is obsolete.
    reservoir_cells_.clear();
    fracture_matrix_ = nullptr;

    this->resetWriters();
}

void
Fracture::initFracture()
{
    Dune::GridFactory<Grid> factory; // Dune::FoamGrid<2,3>> factory;

    const std::size_t N = 6;
    const double radius = 1;

    std::vector<unsigned int> inner_indices;
    std::vector<Dune::FieldVector<double, 3>> vertices;

    vertices.push_back(surfaceMap(0.0, 0.0));

    for (std::size_t i = 0; i < N; ++i) {
        inner_indices.push_back(i + 1);

        const double theta = (i * 2 * M_PI) / N;
        const double x = radius * cos(theta);
        const double y = radius * sin(theta);

        vertices.push_back(surfaceMap(x, y));

        // assume the first 6 cells has source
        well_source_.push_back(i);
    }

    for (std::size_t i = 0; i < vertices.size(); i++) {
        factory.insertVertex(vertices[i]);
    }

    std::vector<std::vector<unsigned int>> cornerIDs;
    for (std::size_t i = 0; i < N; ++i) {
        unsigned int next = inner_indices[(i + 1) % N];
        std::vector<unsigned int> cornerID({unsigned(0), unsigned(inner_indices[i]), next});
        cornerIDs.push_back(cornerID);
    }

    for (std::size_t i = 0; i < N; ++i) {
        factory.insertElement(Dune::GeometryTypes::simplex(2), cornerIDs[i]);
    }

    out_indices_ = inner_indices;
    grid_ = factory.createGrid();
}

std::vector<double>
Fracture::stressIntensityK1() const
{
    const std::size_t nc = numFractureCells();

    std::vector<double> stressIntensityK1(nc, std::nan("0"));

    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& elem : elements(grid_->leafGridView())) {
        const auto& elCenter = elem.geometry().center();

        for (const auto& is : Dune::intersections(grid_->leafGridView(), elem)) {
            if (!is.boundary()) {
                continue;
            }

            const auto distC = (is.geometry().center() - elCenter).two_norm();
            const int nIdx = mapper.index(elem);

            stressIntensityK1[nIdx] = ddm::fractureK1(distC, fracture_width_[nIdx], this->E_, this->nu_);
        }
    }

    return stressIntensityK1;
}

void
Fracture::write(int reportStep) const
{
    std::vector<double> K1; // need to be in scope until written

    auto tmp = this->fracture_pressure_; //@@
    tmp.resize(numFractureCells());

    if (!reservoir_cells_.empty()) {
        vtkwriter_->addCellData(reservoir_cells_, "ReservoirCell");
    }

    if (!reservoir_perm_.empty()) {
        vtkwriter_->addCellData(reservoir_perm_, "ReservoirPerm");
    }

    if (!reservoir_cstress_.empty()) {
        vtkwriter_->addCellData(reservoir_cstress_, "ReservoirCStress");
    }

    if (!reservoir_dist_.empty()) {
        vtkwriter_->addCellData(reservoir_dist_, "ReservoirDist");
    }

    //  Note: size() > 0 since BlockVector<> does not have an "empty()" predicate.
    if (fracture_pressure_.size() > 0) {
        // since fracture_pressure_ may have an additional entry for well
        // pressure, we need to ensure we pass the correct size (@@ ideally this
        // should be fixed by moving the well pressure value into perf_pressure_
        // once it has been computed, and shrink fracture_pressure_ accordingly.

        vtkwriter_->addCellData(tmp, "FracturePressure");
    }

    if (!reservoir_pressure_.empty()) {
        vtkwriter_->addCellData(reservoir_pressure_, "ReservoirPressure");
    }

    if (!reservoir_mobility_.empty()) {
        vtkwriter_->addCellData(reservoir_mobility_, "ReservoirMobility");
    }

    //  Note: size() > 0 since BlockVector<> does not have an "empty()" predicate.
    if (fracture_width_.size() > 0) {
        vtkwriter_->addCellData(fracture_width_, "FractureWidth");
        K1 = this->stressIntensityK1();
        vtkwriter_->addCellData(K1, "stressIntensityK1");
    }

    std::string filename = prm_.get<std::string>("outputdir") + "/Static_"
        + prm_.get<std::string>("casename") + this->name();
    if (reportStep > 0) {
        filename = filename + "_step_" + std::to_string(reportStep);
    }

    vtkwriter_->write(filename.c_str());
}

void
Fracture::writemulti(double time) const
{
    std::cout << "Writing fracture data to VTK files at time: " << time << "grid_size"
              << numFractureCells() << std::endl;

    //  need to have copies in case of async outout (and interface to functions)
    std::vector<double> K1 = this->stressIntensityK1();
    std::vector<double> fracture_pressure(numFractureCells(), 0.0);
    std::vector<double> filtercake_thikness = filtercake_thikness_; //.size(),0.0);
    std::vector<double> reservoir_pressure = reservoir_pressure_; //.size(),0.0);
    std::vector<double> reservoir_dist = reservoir_dist_; //.size(),0.0);
    std::vector<double> reservoir_perm = reservoir_perm_; //.size(),0.0);
    std::vector<double> reservoir_cstress = reservoir_cstress_; //.size(),0.0);
    std::vector<double> reservoir_mobility = reservoir_mobility_; //.size(),0.0);
    std::vector<double> reservoir_traction(reservoir_stress_.size(), 0);
    std::vector<double> fracture_force(reservoir_stress_.size(), 0);
    std::vector<double> reservoir_cells(reservoir_cells_.size(), 0.0);
    std::vector<double> fracture_width(fracture_width_.size(), 0.0);
    std::vector<double> flow_width(fracture_width_.size(), 0.0);
    std::vector<double> rhs_width(rhs_width_.size(), 0.0);
    std::vector<double> well_index(numFractureCells(), 0.0);
    std::vector<double> leakofrate = leakOfRate();

    // make map to do it easy
    std::map<int, double> wellIndMap;
    for (const auto& wind : this->wellIndices()) {
        wellIndMap.insert_or_assign(wind.cell, wind.ctf);
    }

    // loop for need things only with fracture
    for (std::size_t i = 0; i < rhs_width_.size(); ++i) {
        rhs_width[i] = rhs_width_[i][0];
    }

    for (std::size_t i = 0; i < fracture_width_.size(); ++i) {
        fracture_width[i] = fracture_width_[i][0];
        flow_width[i] = fracture_width_[i][0] + min_width_;
        fracture_pressure[i] = fracture_pressure_[i][0];
        reservoir_cells[i] = reservoir_cells_[i]; // only converts to double
        reservoir_traction[i] = ddm::tractionSymTensor(reservoir_stress_[i], cell_normals_[i]);
        fracture_force[i] = reservoir_traction[i] - fracture_pressure[i];
        assert(filtercake_thikness[i] >= 0.0);
    }

    for (std::size_t i = 0; i < numFractureCells(); ++i) {
        // maybe slow
        well_index[i] = wellIndMap[reservoir_cells_[i]];
    }

    vtkmultiwriter_->beginWrite(time);

    if (!reservoir_cells.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_cells, "ReservoirCell");
    }

    if (!reservoir_perm.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_perm, "ReservoirPerm");
    }

    if (!reservoir_dist.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_dist, "ReservoirDist");
    }

    if (!reservoir_cstress.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_cstress, "ReservoirCStres");
    }

    if (!fracture_pressure.empty()) {
        vtkmultiwriter_->attachScalarElementData(fracture_pressure, "FracturePressure");
    }

    if (!reservoir_pressure.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_pressure, "ReservoirPressure");
    }

    if (!filtercake_thikness.empty()) {
        vtkmultiwriter_->attachScalarElementData(filtercake_thikness, "FilterCakeThickness");
    }

    if (!reservoir_traction.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_traction, "ReservoirTraction");
    }

    if (!fracture_force.empty()) {
        vtkmultiwriter_->attachScalarElementData(fracture_force, "FractureForce");
    }

    if (!reservoir_mobility.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_mobility, "ReservoirMobility");
    }

    if (!rhs_width.empty()) {
        vtkmultiwriter_->attachScalarElementData(rhs_width, "ForceOnFracture");
    }

    if (!leakofrate.empty()) {
        vtkmultiwriter_->attachScalarElementData(leakofrate, "LeakOfRate");
    }

    if (!well_index.empty()) {
        vtkmultiwriter_->attachScalarElementData(well_index, "WellIndex");
    }

    if (!fracture_width.empty()) {
        vtkmultiwriter_->attachScalarElementData(fracture_width, "FractureWidth");
        vtkmultiwriter_->attachScalarElementData(flow_width, "FlowWidth");
        vtkmultiwriter_->attachScalarElementData(K1, "stressIntensityK1");
    }

    vtkmultiwriter_->endWrite(false);
}

void
Fracture::grow(int layers, int method)
{
    while (layers_ < layers) {
        std::vector<unsigned int> inner_indices = out_indices_;
        out_indices_.resize(0);

        if (method == 0) {
            this->insertExp(inner_indices);
        } else {
            this->insertLinear(inner_indices);
            nlinear_ += 1;
        }

        ++layers_;

        inner_indices = out_indices_;
    }
}

void
Fracture::insertLinear(const std::vector<unsigned int>& inner_indices)
{
    const double radius = 1.0;
    const std::size_t N = inner_indices.size();

    for (std::size_t i = 0; i < N; ++i) {
        double theta = (i * 2 * M_PI) / (N);
        theta += (layers_ - nlinear_) * 0.5 * (2 * M_PI / N);

        const double out_radius = radius + layers_ + 1;
        const double x = out_radius * cos(theta);
        const double y = out_radius * sin(theta);
        const int new_ind = grid_->insertVertex(surfaceMap(x, y));

        out_indices_.push_back(new_ind);
    }

    for (std::size_t i = 0; i < N; ++i) {
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i % N], out_indices_[(i) % (N)], inner_indices[(i + 1) % (N)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }

        {
            std::vector<unsigned int> cornerID(
                {inner_indices[(i + 1) % N], out_indices_[(i) % (N)], out_indices_[(i + 1) % (N)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
    }
}

void
Fracture::insertExp(const std::vector<unsigned int>& inner_indices)
{
    const std::size_t N = inner_indices.size();
    const double radius = 1.0;

    for (std::size_t i = 0; i < N * 2; ++i) {
        const double theta = (i * 2 * M_PI) / (N * 2);
        const double out_radius = radius + layers_ + 1;
        const double x = out_radius * cos(theta);
        const double y = out_radius * sin(theta);
        const int new_ind = grid_->insertVertex(surfaceMap(x, y));

        out_indices_.push_back(new_ind);
    }

    for (std::size_t i = 0; i < N; ++i) {
        {
            std::vector<unsigned int> cornerID({inner_indices[i],
                                                out_indices_[(2 * i) % (N * 2)],
                                                out_indices_[(2 * i + 1) % (N * 2)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }

        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i], out_indices_[(2 * i + 1) % (N * 2)], inner_indices[(i + 1) % N]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }

        {
            std::vector<unsigned int> cornerID({inner_indices[(i + 1) % N],
                                                out_indices_[(2 * i + 1) % (N * 2)],
                                                out_indices_[(2 * i + 2) % (N * 2)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
    }
}

Dune::FieldVector<double, 3>
Fracture::surfaceMap(const double x, const double y)
{
    Point3D vec(0.0);

    vec += x * axis_[0];
    vec += y * axis_[1];
    vec += origo_;

    return vec;
}

void
Fracture::updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree)
{
    reservoir_cells_.resize(numFractureCells());

    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    const ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());

    int tri_divide = 0;
    int tri_outside = 0;
    for (auto& element : elements(grid_->leafGridView())) {
        const auto elemIdx = elemMapper.index(element);
        const auto geom = element.geometry();
        external::cvf::BoundingBox bb;

        using Vec3d = external::cvf::Vec3d;
        {
            // only lock at centroid
            const auto& vertex = geom.center();
            Vec3d point(vertex[0], vertex[1], vertex[2]);
            bb.add(point);
        }

        std::vector<std::size_t> cells = external::findCloseCellIndices(cellSearchTree, bb);
        if (cells.empty()) {
            reservoir_cells_[elemIdx] = cells[0];
            if (cells.size() > 1) {
                tri_divide += 1;
            }
        } else {
            tri_outside += 1;
            reservoir_cells_[elemIdx] = -1;
        }
    }

    std::cout << "For Fracture : " << this->name() << " : " << tri_divide
              << " triangles should be devided" << std::endl;
    std::cout << "For Fracture : " << this->name() << " : " << tri_outside << " triangles outside"
              << std::endl;
    std::cout << "Total triangles: " << numFractureCells() << std::endl;

    auto it = std::find(reservoir_cells_.begin(), reservoir_cells_.end(), -1);

    const auto extended_fractures = prm_.get<bool>("extended_fractures");
    if ((it != reservoir_cells_.end()) && !extended_fractures) {
        std::cout << "Remove fracture outside of model" << std::endl;
        // remove fracture outside of model
        this->removeCells();
        this->updateReservoirCells(cellSearchTree);
    }
}

void
Fracture::updateReservoirProperties()
{
    // updater for standalone test
    const double perm = prm_.get<double>("reservoir.perm");
    const double dist = prm_.get<double>("reservoir.dist");
    const double cstress = prm_.get<double>("KMax");
    const std::size_t nc = numFractureCells();

    reservoir_perm_.resize(nc, perm);
    reservoir_dist_.resize(nc, dist);
    reservoir_mobility_.resize(nc, 1000);
    reservoir_pressure_.resize(nc, 100.0e5);
    reservoir_stress_.resize(nc);
    reservoir_cstress_.resize(nc, cstress);

    for (std::size_t i = 0; i != nc; ++i) {
        reservoir_stress_[i] = Dune::FieldVector<double, 6> {0, 0, 0, 0, 0, 0};
    }

    nu_ = 0.25;
    E_ = 1e9;

    this->initFractureWidth();
}

void
Fracture::addSource()
{
    if (rhs_pressure_.size() == 0) {
        std::size_t nc = numFractureCells();
        rhs_pressure_.resize(nc + numWellEquations());
    }

    assert(rhs_pressure_.size() == numFractureCells() + numWellEquations());
    rhs_pressure_ = 0;
    // add contributions for fracture fracture gravity contributions
    int nc = grid_->leafGridView().size(0);
    fracture_dgh_.resize(nc, 0.0);
    for (auto element : Dune::elements(grid_->leafGridView())) {
        // int i = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView,
        // Dune::mcmgElementLayout>::index(element);
        int i = grid_->leafGridView().indexSet().index(element);
        double z = element.geometry().center()[2];
        fracture_dgh_[i] = gravity_ * reservoir_density_[i] * z;
    }

    // could be put into the assemble loop
    for (const auto& [i, j, t1, t2] : htrans_) {
        const double h1 = fracture_width_[i] + min_width_;
        const double h2 = fracture_width_[j] + min_width_;

        // harmonic mean of surface flow
        double value = 12. / (h1 * h1 * h1 * t1) + 12. / (h2 * h2 * h2 * t2);

        const double mobility = 0.5 * (reservoir_mobility_[i] + reservoir_mobility_[j]);

        value = 1 / value;
        value *= mobility;

        const double dh = (fracture_dgh_[i] - fracture_dgh_[j]);
        rhs_pressure_[i] -= value * dh;
        rhs_pressure_[j] += value * dh;
    }

    for (std::size_t i = 0; i < reservoir_pressure_.size(); ++i) {
        rhs_pressure_[i] += leakof_[i] * reservoir_pressure_[i];
    }

    // gravity contribution from fracture to reservoir
    for (auto element : Dune::elements(grid_->leafGridView())) {
        const int i = grid_->leafGridView().indexSet().index(element);
        const double z = element.geometry().center()[2];

        rhs_pressure_[i] += leakof_[i] * (z - reservoir_cell_z_[i]) * gravity_ * reservoir_density_[i];
    }

    // gravity contributions between fracture cells

    const auto control = prm_.get_child("control");
    const std::string control_type = control.get<std::string>("type");

    if (control_type == "rate") {
        const double scale = well_source_.size();
        const double rate = control.get<double>("rate");
        for (const auto& cell : well_source_) {
            rhs_pressure_[cell] += rate / scale;
        }
    } else if (control_type == "pressure") {
        const double pressure = control.get<double>("pressure");
        for (const auto& perfinj : perfinj_) {
            const int cell = std::get<0>(perfinj);
            const double value = std::get<1>(perfinj);
            const double dh_perf = origo_[2] * gravity_ * reservoir_density_[cell];
            const double dh_cell = fracture_dgh_[cell];
            const double dh = dh_perf - dh_cell;
            rhs_pressure_[cell] += value * (pressure - dh);
        }
    } else if (control_type == "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            const int cell = std::get<0>(perfinj);
            const double value = std::get<1>(perfinj);
            const double dh_perf = origo_[2] * gravity_ * reservoir_density_[cell];
            const double dh_cell = fracture_dgh_[cell];
            const double dh = dh_perf - dh_cell;
            rhs_pressure_[cell] += value * (perf_pressure_ - dh);
        }
    } else if (control_type == "rate_well") {
        assert(numWellEquations() == 1); // @@ for now, we assume there is just one well equation
        const double r = control.get<double>("rate") / 24 / 60 / 60; // convert to m3/sec
        const double WI = control.get<double>("WI");
        const int cell = std::get<0>(perfinj_[0]); // @@ will this be the correct index?
        const double pres = reservoir_pressure_[cell];
        const double lambda = reservoir_mobility_[0]; // @@ only correct if mobility is constant!
        rhs_pressure_[rhs_pressure_.size() - 1] = r + WI * lambda * pres; // well source term
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }
}

double
Fracture::injectionPressure() const
{
    const std::string control_type = prm_.get<std::string>("control.type");

    if (control_type == "rate") {
        double bhp = 0.0;
        double scale = well_source_.size();

        // could have corrected for WI
        for (const auto& cell : well_source_) {
            bhp += fracture_pressure_[cell] / scale;
        }

        return bhp;
    } else if (control_type == "pressure") {
        return prm_.get<double>("control.pressure");
    } else if (control_type == "perf_pressure") {
        return perf_pressure_;
    } else if (control_type == "rate_well") {
        // @@ We should use perf_pressure_ here too, but ensure it is updated
        return fracture_pressure_[fracture_pressure_.size() - 1][0];
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }

    return 0.0;
}

std::vector<double>
Fracture::leakOfRate() const
{
    if (leakof_.empty()) {
        return {};
    }

    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());

    std::vector<double> leakofrate(numFractureCells(), 0);
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const double dp = (fracture_pressure_[eIdx] - reservoir_pressure_[eIdx]);
        const double q = leakof_[eIdx] * dp;

        leakofrate[eIdx] = q / element.geometry().volume();
    }

    return leakofrate;
}

std::vector<RuntimePerforation>
Fracture::wellIndices() const
{
    // find unique reservoir cells
    if (leakof_.empty()) {
        // if pressure is not assembled return empty
        return {};
    }

    std::vector<int> res_cells = reservoir_cells_;
    std::sort(res_cells.begin(), res_cells.end());

    auto last = std::unique(res_cells.begin(), res_cells.end());
    res_cells.erase(last, res_cells.end());

    std::vector<double> q_cells(res_cells.size(), 0.0);
    std::vector<double> p_cells(res_cells.size(), 0.0);

    double q_prev = 0;
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const double dp = (fracture_pressure_[eIdx] - reservoir_pressure_[eIdx]);
        const double q = leakof_[eIdx] * dp;
        const int res_cell = reservoir_cells_[eIdx];

        // just to search
        auto it = std::find(res_cells.begin(), res_cells.end(), res_cell);
        const int ind_wellIdx = it - res_cells.begin();

        assert(it != res_cells.end());
        if (q_prev * q < 0) {
            OPM_THROW(std::runtime_error, "Cross flow in fracture ??");
        }

        q_cells[ind_wellIdx] += q;
        p_cells[ind_wellIdx] = reservoir_pressure_[eIdx]; // is set multiple times
    }

    std::vector<RuntimePerforation> wellIndices(res_cells.size());
    double inj_press = injectionPressure();
    for (std::size_t i = 0; i < res_cells.size(); ++i) {
        auto& perf = wellIndices[i];

        perf.cell = res_cells[i];

        const double dh_res = reservoir_cell_z_[i] * gravity_ * reservoir_density_[i];
        const double perf_density = reservoir_density_[i];
        const double dh_perf = gravity_ * perf_density * origo_[2];

        double WI = q_cells[i] / ((inj_press - dh_perf) - (p_cells[i] - dh_res));

        if (WI < 0.0) {
            std::cout << "Negative WI: " << WI << " for cell: " << res_cells[i] << std::endl;
            WI = 0.0;
        }

        perf.ctf = WI;
        perf.depth = this->origo_[2];
        perf.segment = this->wellinfo_.segment;
        perf.perf_range = this->wellinfo_.perf_range;
    }

    return wellIndices;
}

template <typename Scalar>
void
Fracture::assignGeomechWellState(ConnFracStatistics<Scalar>& stats) const
{
    using Quantity = typename ConnFracStatistics<Scalar>::Quantity;

    constexpr auto pressIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::Pressure);
    constexpr auto rateIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::FlowRate);
    constexpr auto widthIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::Width);

    const auto nCells = this->reservoir_cells_.size();

    stats.reset();

    for (auto cellIx = 0 * nCells; cellIx < nCells; ++cellIx) {
        auto samplePoint = typename ConnFracStatistics<Scalar>::SamplePoint {};

        samplePoint[pressIx] = this->fracture_pressure_[cellIx][0];

        samplePoint[rateIx] = this->leakof_[cellIx]
            * (this->fracture_pressure_[cellIx][0] - this->reservoir_pressure_[cellIx]);

        samplePoint[widthIx] = this->fracture_width_[cellIx][0];

        stats.addSamplePoint(samplePoint);
    }
}

void
Fracture::writePressureSystem() const
{
    if (prm_.get<bool>("write_pressure_system")) {
        Dune::storeMatrixMarket(*pressure_matrix_, "pressure_matrix");
        Dune::storeMatrixMarket(rhs_pressure_, "pressure_rhs");
    }
}

void
Fracture::writeFractureSystem() const
{
    if (prm_.get<bool>("write_fracture_system")) {
        // Dune::storeMatrixMarket(*fracture_matrix_, "fracture_matrix");
        Dune::storeMatrixMarket(rhs_width_, "rhs_width");
    }
}

void
Fracture::solvePressure()
{
    OPM_TIMEFUNCTION();

    const std::size_t nc = numFractureCells();

    assert(numWellEquations() == 0); // @@ not implemented/tested for "rate_well" control
    fracture_pressure_.resize(nc);
    fracture_pressure_ = 1e5;
    assert(pressure_matrix_); // should always be constructed at this pointn

    this->assemblePressure();
    this->addSource(); // probably include reservoir pressure
    this->writePressureSystem();

    try {
        if (!pressure_solver_) {
            this->setupPressureSolver();
        }

        fracture_pressure_.resize(rhs_pressure_.size());
        fracture_pressure_ = 0;

        Dune::InverseOperatorResult r {};
        pressure_solver_->apply(fracture_pressure_, rhs_pressure_, r);
    } catch (Dune::ISTLError& e) {
        std::cerr << "exception thrown " << e << std::endl;
    }
}

// based on the (given) values of fracture_pressure__, compute rhs_width_ and
// fracture_width_
void
Fracture::solveFractureWidth()
{
    fractureMatrix().solve(fracture_width_, rhs_width_);

    const double max_width = prm_.get<double>("solver.max_width");
    const double min_width = prm_.get<double>("solver.min_width");

    for (auto& width : this->fracture_width_) {
        assert(std::isfinite(width));

        if (width > max_width) {
            std::cout << "Limit Fracture width" << std::endl;
            width = max_width;
        }

        if (width < min_width) {
            std::cout << "Remove small Fracture width" << std::endl;
            width = min_width;
        }

        assert(std::isfinite(width));
    }
}

void
Fracture::initFractureStates()
{
    this->initFractureWidth();
    this->initFracturePressureFromReservoir();
}

void
Fracture::initFractureWidth()
{
    std::size_t nc = numFractureCells();
    fracture_width_.resize(nc);
    fracture_width_ = prm_.get<double>("config.initial_fracture_width");
}

void
Fracture::initFracturePressureFromReservoir()
{
    std::size_t nc = reservoir_cells_.size();
    fracture_pressure_.resize(nc + numWellEquations());
    fracture_pressure_ = 0;
    for (std::size_t i = 0; i < nc; ++i) {
        fracture_pressure_[i] = reservoir_pressure_[i];
    }
}

void
Fracture::updateLeakoff()
{
    const std::size_t nc = numFractureCells();
    leakof_.resize(nc, 0.0);

    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());

    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto geom = element.geometry();
        const double area = geom.volume();
        const double res_mob = reservoir_mobility_[eIdx];

        leakof_[eIdx] = res_mob * reservoir_perm_[eIdx] * area / reservoir_dist_[eIdx];

        if (has_filtercake_) {
            double invtrans = 1 / leakof_[eIdx];
            assert(filtercake_thikness_[eIdx] >= 0.0);
            if (filtercake_thikness_[eIdx] > 0.0) {
                double fitercaketrans = res_mob * filtercake_perm_ * area / filtercake_thikness_[eIdx];
                invtrans += 1 / fitercaketrans;
            }

            leakof_[eIdx] = 1 / invtrans;
        }
    }
}

void
Fracture::redistribute_values(Dune::BlockVector<Dune::FieldVector<double, 1>>& values,
                              const std::vector<std::vector<CellRef>>& map1,
                              const std::vector<std::vector<CellRef>>& map2,
                              const int level,
                              const bool point_wise)
{
    std::vector<double> tmp_values(values.size());

    for (std::size_t i = 0; i < values.size(); ++i) {
        tmp_values[i] = values[i][0]; // assuming values is a BlockVector with one component
    }

    tmp_values = redistribute_values(tmp_values, map1, map2, level, point_wise);

    values.resize(tmp_values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i][0] = tmp_values[i]; // assuming values is a BlockVector with one component
    }
}

std::vector<double>
Fracture::redistribute_values(const std::vector<double>& values,
                              const std::vector<std::vector<CellRef>>& map1,
                              const std::vector<std::vector<CellRef>>& map2,
                              const int level,
                              const bool point_wise)
// ----------------------------------------------------------------------------
{
    const auto g2gmap = RegularTrimesh::createGridToGridMap(map1, map2, level);

    std::vector<double> redistributed_values(map2.size(), 0.0);
    std::vector<double> weight(map2.size(), 0.0);

    for (const auto& e : g2gmap) {
        redistributed_values[std::get<1>(e)] += values[std::get<0>(e)] * std::get<2>(e);
        weight[std::get<1>(e)] += std::get<2>(e);
    }

    if (point_wise) {
        for (std::size_t i = 0; i < redistributed_values.size(); ++i) {
            assert(weight[i] >= 0.0);
            if (weight[i] > 0.0) {
                redistributed_values[i] /= weight[i];
            } else {
                redistributed_values[i] = 0.0; // or some other default value
            }
        }
    }

    return redistributed_values;
}

void
Fracture::initPressureMatrix()
{
    // Flow from wells to fracture cells
    const double fWI = prm_.get<double>("fractureWI");

    perfinj_.clear();
    htrans_.clear();

    for (const auto& cell : well_source_) {
        perfinj_.emplace_back(cell, fWI);
    }

    // flow between fracture cells
    const std::size_t nc = numFractureCells() + numWellEquations();

    // leakof_.resize(nc,0.0);
    const ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());

    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto geom = element.geometry();
        const auto eCenter = geom.center();

        // iterator over all intersections
        for (const auto& is : Dune::intersections(grid_->leafGridView(), element)) {
            if (is.boundary()) {
                continue;
            }

            const int nIdx = mapper.index(is.outside());
            if (!(eIdx < nIdx)) {
                continue;
            }

            // calculate distance between the midpoints
            const auto nCenter = is.outside().geometry().center();
            const auto isCenter = is.geometry().center();
            const auto d_inside = eCenter - isCenter;
            const auto d_outside = nCenter - isCenter;

            // probably should use projected distance
            const auto igeom = is.geometry();
            const double area = igeom.volume();
            const double h1 = area / d_inside.two_norm();
            const double h2 = area / d_outside.two_norm();

            {
                Htrans matel(nIdx, eIdx, h1, h2);
                htrans_.push_back(matel);
            }
        }
    }

    pressure_matrix_ = std::make_unique<Matrix>(nc, nc, 4, 0.4, Matrix::implicit);

    auto& matrix = *pressure_matrix_;

    for (const auto& matel : htrans_) {
        const std::size_t i = std::get<0>(matel);
        const std::size_t j = std::get<1>(matel);
        const double zero_entry = 0.0; // 1e-11;

        matrix.entry(i, j) = zero_entry; // 0;
        matrix.entry(j, i) = zero_entry; // 0;
        matrix.entry(j, j) = zero_entry; // 0;
        matrix.entry(i, i) = zero_entry; // 0;
    }

    if (numWellEquations() > 0) {
        assert(numWellEquations() == 1); // @@ for now, we assume there is just one well equation
        // add elements for last row and column
        const int weqix = nc - 1; // index of well equation
        for (int i : well_source_) {
            matrix.entry(i, weqix) = 0.0;
            matrix.entry(weqix, i) = 0.0;
        }

        matrix.entry(weqix, weqix) = 0.0;
    }

    matrix.compress();
}

void
Fracture::assemblePressure()
{
    updateLeakoff();

    auto& matrix = *pressure_matrix_;
    matrix = 0.0;

    // double mobility=1e4; //1e4; // @@ 1.0
    //  get head in all fracture cells
    for (const auto& [i, j, t1, t2] : htrans_) {
        const double h1 = fracture_width_[i] + min_width_;
        const double h2 = fracture_width_[j] + min_width_;

        // harmonic mean of surface flow
        const double mobility = 0.5 * (reservoir_mobility_[i] + reservoir_mobility_[j]);

        double value = 12. / (h1 * h1 * h1 * t1) + 12. / (h2 * h2 * h2 * t2);
        value = mobility / value;

        matrix[i][j] -= value;
        matrix[j][i] -= value;
        matrix[i][i] += value;
        matrix[j][j] += value;
    }

    const auto control = prm_.get_child("control");
    const std::string control_type = control.get<std::string>("type");
    for (std::size_t i = 0; i < leakof_.size(); ++i) {
        // matrix.entry(i, i) += leakof_[i];
        matrix[i][i] += leakof_[i];
    }

    if (control_type == "rate") {
        // no extra tings in matrix
    } else if (control_type == "pressure" || control_type == "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            matrix[cell][cell] += value;
        }
    } else if (control_type == "rate_well") {
        assert(numWellEquations() == 1); // @@ for now, we assume there is just one well equation
        const int nc = numFractureCells() + numWellEquations();
        const double lambda = reservoir_mobility_[0]; // @@ If not constant, this might be wrong
        const double WI_lambda = control.get<double>("WI") * lambda;
        matrix[nc - 1][nc - 1] = WI_lambda;
        // NB: well_source_[i] is assumed to be the same as get<0>(perfinj_[i])
        for (const auto& pi : perfinj_) {
            const int i = std::get<0>(pi);
            const double value = std::get<1>(pi) * lambda;
            matrix[nc - 1][i] = -value; // well equation
            matrix[nc - 1][nc - 1] += value; // well equation
            matrix[i][nc - 1] = -value;
            matrix[i][i] += value;
        }
    } else {
        OPM_THROW(std::runtime_error, "Unknown control of injection into Fracture");
    }
}

double
Fracture::normalFractureTraction(std::size_t eIdx) const
{
    return ddm::tractionSymTensor(reservoir_stress_[eIdx], cell_normals_[eIdx]);
}

void
Fracture::normalFractureTraction(Dune::BlockVector<Dune::FieldVector<double, 1>>& traction,
                                 const bool resize) const
{
    OPM_TIMEFUNCTION();

    const std::size_t nc = numFractureCells();
    if (resize) {
        traction.resize(nc + numWellEquations());
    }

    for (std::size_t eIdx = 0; eIdx < nc; ++eIdx) {
        traction[eIdx] = normalFractureTraction(eIdx);
    }
}

void
Fracture::updateFractureRHS()
{
    assert(numWellEquations() == 0); // @@ not implemented/tested for rate-controlled systems
    rhs_width_ = fracture_pressure_;

    for (std::size_t i = 0; i < rhs_width_.size(); ++i) {
        std::cout << i << " " << rhs_width_[i] << std::endl;
        rhs_width_[i] = rhs_width_[i] - normalFractureTraction(i);
        if (rhs_width_[i] < 0.0) {
            rhs_width_[i] = 0.0; // @@ not entirely accurate, but will avoid
                                 // unphysical negative normal displacements
        }
    }
}

void
Fracture::assembleFractureMatrix() const
{
    OPM_TIMEFUNCTION();

    std::size_t nc = numFractureCells();
    if (!fracture_matrix_) {
        fracture_matrix_ = std::make_unique<DynamicMatrix>();
    }

    fracture_matrix_->resize(nc, nc);
    *fracture_matrix_ = 0.0;

    ddm::assembleMatrix(*fracture_matrix_, E_, nu_, *grid_);
}

void
Fracture::printPressureMatrix() const // debug purposes
{
    Dune::printSparseMatrix(std::cout, *pressure_matrix_, "matname", "linameo");
}

void
Fracture::printMechMatrix() const // debug purposes
{
    Dune::printmatrix(std::cout, fractureMatrix(), "matname", "linameo");
}

template void Fracture::assignGeomechWellState(ConnFracStatistics<float>&) const;
template void Fracture::assignGeomechWellState(ConnFracStatistics<double>&) const;

} // namespace Opm
