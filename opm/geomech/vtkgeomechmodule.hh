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
/*!
 * \file
 * \copydoc Opm::VtkEclTracerModule
 */
#ifndef VTK_GEOMECH_MODULE_HH
#define VTK_GEOMECH_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters
{

// create new type tag for the VTK tracer output
// create the property tags needed for the tracer model
// set default values for what quantities to output

struct VtkWriteGeoMech
{
    static constexpr bool value = true;
};

} // namespace Opm::Parameters

namespace Opm
{
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the tracer model's parameters.
 */
template <class TypeTag>
class VtkGeoMechModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;


    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using VectorBuffer = typename ParentType::VectorBuffer;
    using TensorBuffer = typename ParentType::TensorBuffer;
    using Tensor = Dune::DynamicMatrix<double>;
    using SymTensor = Dune::FieldVector<double, 6>;

public:
    VtkGeoMechModule(const Simulator& simulator)
        : ParentType(simulator)
    {
    }

    /*!
     * \brief Register all run-time parameters for the tracer VTK output
     * module.
     */
    static void registerParameters()
    {
        Parameters::Register<Parameters::VtkWriteGeoMech>("Include geomech quentities "
                                                          "in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (!geoMechOutput_()) {
            return;
        }

        this->resizeScalarBuffer_(pressDiff_, ParentType::BufferType::Element);
        this->resizeVectorBuffer_(disp_, ParentType::BufferType::Vertex);
        this->resizeTensorBuffer_(stress_, ParentType::BufferType::Element);
        this->resizeTensorBuffer_(delstress_, ParentType::BufferType::Element);
        this->resizeTensorBuffer_(strain_, ParentType::BufferType::Element);

        this->resizeScalarBuffer_(pressure_, ParentType::BufferType::Element);
        this->resizeScalarBuffer_(biot_, ParentType::BufferType::Element);
        this->resizeScalarBuffer_(ymodule_, ParentType::BufferType::Element);
        this->resizeScalarBuffer_(pratio_, ParentType::BufferType::Element);
        this->resizeScalarBuffer_(poro_, ParentType::BufferType::Element);
    }

    void processElement(const ElementContext& elemCtx)
    {
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        const auto& problem = elemCtx.problem();
        const auto& geoMechModel = problem.geoMechModel();

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            // get the global index of the dof
            const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            pressDiff_[globalDofIdx] = geoMechModel.pressureDiff(globalDofIdx);
            poro_[globalDofIdx] = problem.porosity(globalDofIdx, /*timeIdx=*/0);
            biot_[globalDofIdx] = problem.biotCoef(globalDofIdx);
            pratio_[globalDofIdx] = problem.pRatio(globalDofIdx);
            ymodule_[globalDofIdx] = problem.yModule(globalDofIdx);

            // give the pressure for the geomech model
            pressure_[globalDofIdx] = geoMechModel.pressure(globalDofIdx);

            {
                const SymTensor& symtensor = geoMechModel.stress(globalDofIdx);
                setTensor(this->stress_[globalDofIdx], symtensor);
            }

            {
                const SymTensor& symtensor = geoMechModel.delstress(globalDofIdx);
                setTensor(this->delstress_[globalDofIdx], symtensor);
            }

            {
                const SymTensor& symtensor = geoMechModel.strain(globalDofIdx);
                setTensor(this->strain_[globalDofIdx], symtensor);
            }
        }

        // all vertices proably do it to many times for now
        auto gv = elemCtx.gridView();
        auto elem = elemCtx.element();

        static constexpr int Dim = 3;
        for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<Dim> {})) {
            auto index = gv.indexSet().index(vertex);
            disp_[index] = geoMechModel.displacement(index);
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        if (!geoMechOutput_()) {
            return;
        }

        auto* const vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (vtkWriter == nullptr) {
            return;
        }

        this->commitScalarBuffer_(
            baseWriter, "pressureDiff", pressDiff_, ParentType::BufferType::Element);

        this->commitVectorBuffer_(baseWriter, "disp", disp_, ParentType::BufferType::Vertex);
        this->commitTensorBuffer_(baseWriter, "stress", stress_, ParentType::BufferType::Element);
        this->commitTensorBuffer_(baseWriter, "delstress", delstress_, ParentType::BufferType::Element);
        this->commitTensorBuffer_(baseWriter, "strain", strain_, ParentType::BufferType::Element);
    }

private:

    /*!
     * \brief Modify the internal buffers according to the intensive
     * quantities relevant for an element
     */
    static void setTensor(Tensor& tensor, const SymTensor& symtensor)
    {
        for (int i = 0; i < 3; ++i) {
            tensor[i][i] = symtensor[i];
        }

        // Voigt notation conversion
        tensor[0][1] = symtensor[5]; // xy
        tensor[0][2] = symtensor[4]; // xz
        tensor[1][2] = symtensor[3]; // yz

        // fix symmetry of tensor
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                tensor[i][j] = tensor[j][i];
            }
        }
    }

    static bool geoMechOutput_()
    {
        static bool val = Parameters::Get<Parameters::VtkWriteGeoMech>();
        return val;
    }

    ScalarBuffer pressDiff_;
    VectorBuffer disp_;
    VectorBuffer symstress_;
    TensorBuffer stress_;
    TensorBuffer delstress_;
    TensorBuffer strain_;

    ScalarBuffer pressure_;
    ScalarBuffer biot_;
    ScalarBuffer ymodule_;
    ScalarBuffer pratio_;
    ScalarBuffer poro_;
};

} // namespace Opm

#endif // VTK_GEOMECH_MODULE_HH
