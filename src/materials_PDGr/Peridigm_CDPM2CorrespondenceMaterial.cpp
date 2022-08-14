/*! \file Peridigm_CDPM2CorrespondenceMaterial.cpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "Peridigm_CDPM2CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "cdpm2_correspondence.h"
#include <Teuchos_Assert.hpp>
#include "correspondence.h" // to use the MatrixMultiply and TransposeMatrix functions

using namespace std;

PeridigmNS::CDPM2CorrespondenceMaterial::CDPM2CorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_unitFlag(0.0), m_youngsModulus(0.0), m_poissonsRatio(0.0), m_fc(0.0), m_ft(0.0), m_fracEnergy(0.0),
    m_hard1(0.0), m_hard2(0.0), m_hard3(0.0), m_hard4(0.0), m_hardqh2(0.0), m_soften(0.0), m_length(0.0),
    m_unrotatedRateOfDeformationFieldId(-1), 
    m_unrotatedCauchyStressFieldId(-1), 
    m_leftStretchTensorFieldId(-1), 
    m_rotationTensorFieldId(-1),
    m_tempDeformationGradientFieldId(-1),
    m_GreenLagrangeStrainFieldId(-1), 
    m_microState1FieldId(-1),
    m_microState2FieldId(-1),
    m_microState3FieldId(-1),
    m_microState4FieldId(-1),
    m_microState5FieldId(-1),
    m_microState6FieldId(-1),
    m_microState7FieldId(-1),
    m_microState8FieldId(-1),
    m_microState9FieldId(-1),
    m_microState10FieldId(-1),
    m_microState11FieldId(-1),
    m_microState12FieldId(-1),
    m_microState13FieldId(-1),
    m_microState14FieldId(-1),
    m_microState15FieldId(-1),
    m_microState16FieldId(-1),
    m_microState17FieldId(-1),
    m_microState18FieldId(-1),
    m_microState19FieldId(-1),
    m_microState20FieldId(-1),
    m_microState21FieldId(-1),
    m_microState22FieldId(-1),
    m_microState23FieldId(-1),
    m_microState24FieldId(-1),
    m_microState25FieldId(-1),
    m_microState26FieldId(-1),
    m_microState27FieldId(-1),
    m_microState28FieldId(-1),
    m_microState29FieldId(-1),
    m_microState30FieldId(-1),
    m_microState31FieldId(-1),
    m_microState32FieldId(-1),
    m_microState33FieldId(-1),
    m_microState34FieldId(-1),
    m_microState35FieldId(-1),
    m_microState36FieldId(-1),
    m_microState37FieldId(-1),
    m_microState38FieldId(-1),
    m_microState39FieldId(-1),
    m_microState40FieldId(-1),
    m_microState41FieldId(-1),
    m_microState42FieldId(-1),
    m_microState43FieldId(-1),
    m_microState44FieldId(-1),
    m_microState45FieldId(-1),
    m_microState46FieldId(-1),
    m_microState47FieldId(-1),
    m_microState48FieldId(-1),
    m_microState49FieldId(-1),
    m_internalEnergyFieldId(-1), 
    m_inelasticEnergyFieldId(-1)
{
  m_youngsModulus = 9.0 * m_bulkModulus * m_shearModulus / (3.0 * m_bulkModulus + m_shearModulus);
  m_poissonsRatio = m_youngsModulus / (2.0 * m_shearModulus) - 1.0;
  m_unitFlag = params.get<double>("Unit_Flag");
  m_fc = params.get<double>("Compressive_Strength");
  m_ft = params.get<double>("Tensile_Strength");
  m_fracEnergy = params.get<double>("Fracture_Energy");
  m_hard1 = params.get<double>("Hardening_Modulus_1");
  m_hard2 = params.get<double>("Hardening_Modulus_2");
  m_hard3 = params.get<double>("Hardening_Modulus_3");
  m_hard4 = params.get<double>("Hardening_Modulus_4");
  m_hardqh2 = params.get<double>("Hardening_Function_qh2");
  m_soften = params.get<double>("Softening");
  m_length = params.get<double>("Character_Length");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_leftStretchTensorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_tempDeformationGradientFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "temp_Deformation_Gradient");
  m_GreenLagrangeStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Green_Lagrange_Strain");
  m_microState1FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_1");
  m_microState2FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_2");
  m_microState3FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_3");
  m_microState4FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_4");
  m_microState5FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_5");
  m_microState6FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_6");
  m_microState7FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_7");
  m_microState8FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_8");
  m_microState9FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_9");
  m_microState10FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_10");
  m_microState11FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_11");
  m_microState12FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_12");
  m_microState13FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_13");
  m_microState14FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_14");
  m_microState15FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_15");
  m_microState16FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_16");
  m_microState17FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_17");
  m_microState18FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_18");
  m_microState19FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_19");
  m_microState20FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_20");
  m_microState21FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_21");
  m_microState22FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_22");
  m_microState23FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_23");
  m_microState24FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_24");
  m_microState25FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_25");
  m_microState26FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_26");
  m_microState27FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_27");
  m_microState28FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_28");
  m_microState29FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_29");
  m_microState30FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_30");
  m_microState31FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_31");
  m_microState32FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_32");
  m_microState33FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_33");
  m_microState34FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_34");
  m_microState35FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_35");
  m_microState36FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_36");
  m_microState37FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_37");
  m_microState38FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_38");
  m_microState39FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_39");
  m_microState40FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_40");
  m_microState41FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_41");
  m_microState42FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_42");
  m_microState43FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_43");
  m_microState44FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_44");
  m_microState45FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_45");
  m_microState46FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_46");
  m_microState47FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_47");
  m_microState48FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_48");
  m_microState49FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_49");
  m_internalEnergyFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Internal_Energy");
  m_inelasticEnergyFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Inelastic_Energy");
  
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_tempDeformationGradientFieldId);
  m_fieldIds.push_back(m_GreenLagrangeStrainFieldId);
  m_fieldIds.push_back(m_microState1FieldId);
  m_fieldIds.push_back(m_microState2FieldId);
  m_fieldIds.push_back(m_microState3FieldId);
  m_fieldIds.push_back(m_microState4FieldId);
  m_fieldIds.push_back(m_microState5FieldId);
  m_fieldIds.push_back(m_microState6FieldId);
  m_fieldIds.push_back(m_microState7FieldId);
  m_fieldIds.push_back(m_microState8FieldId);
  m_fieldIds.push_back(m_microState9FieldId);
  m_fieldIds.push_back(m_microState10FieldId);
  m_fieldIds.push_back(m_microState11FieldId);
  m_fieldIds.push_back(m_microState12FieldId);
  m_fieldIds.push_back(m_microState13FieldId);
  m_fieldIds.push_back(m_microState14FieldId);
  m_fieldIds.push_back(m_microState15FieldId);
  m_fieldIds.push_back(m_microState16FieldId);
  m_fieldIds.push_back(m_microState17FieldId);
  m_fieldIds.push_back(m_microState18FieldId);
  m_fieldIds.push_back(m_microState19FieldId);
  m_fieldIds.push_back(m_microState20FieldId);
  m_fieldIds.push_back(m_microState21FieldId);
  m_fieldIds.push_back(m_microState22FieldId);
  m_fieldIds.push_back(m_microState23FieldId);
  m_fieldIds.push_back(m_microState24FieldId);
  m_fieldIds.push_back(m_microState25FieldId);
  m_fieldIds.push_back(m_microState26FieldId);
  m_fieldIds.push_back(m_microState27FieldId);
  m_fieldIds.push_back(m_microState28FieldId);
  m_fieldIds.push_back(m_microState29FieldId);
  m_fieldIds.push_back(m_microState30FieldId);
  m_fieldIds.push_back(m_microState31FieldId);
  m_fieldIds.push_back(m_microState32FieldId);
  m_fieldIds.push_back(m_microState33FieldId);
  m_fieldIds.push_back(m_microState34FieldId);
  m_fieldIds.push_back(m_microState35FieldId);
  m_fieldIds.push_back(m_microState36FieldId);
  m_fieldIds.push_back(m_microState37FieldId);
  m_fieldIds.push_back(m_microState38FieldId);
  m_fieldIds.push_back(m_microState39FieldId);
  m_fieldIds.push_back(m_microState40FieldId);
  m_fieldIds.push_back(m_microState41FieldId);
  m_fieldIds.push_back(m_microState42FieldId);
  m_fieldIds.push_back(m_microState43FieldId);
  m_fieldIds.push_back(m_microState44FieldId);
  m_fieldIds.push_back(m_microState45FieldId);
  m_fieldIds.push_back(m_microState46FieldId);
  m_fieldIds.push_back(m_microState47FieldId);
  m_fieldIds.push_back(m_microState48FieldId);
  m_fieldIds.push_back(m_microState49FieldId);
  m_fieldIds.push_back(m_internalEnergyFieldId);
  m_fieldIds.push_back(m_inelasticEnergyFieldId);

}

PeridigmNS::CDPM2CorrespondenceMaterial::~CDPM2CorrespondenceMaterial()
{
}

void
PeridigmNS::CDPM2CorrespondenceMaterial::initialize(    const double dt,
                                                        const int numOwnedPoints,
                                                        const int* ownedIDs,
                                                        const int* neighborhoodList,
                                                        PeridigmNS::DataManager& dataManager)
{

  PeridigmNS::CorrespondenceMaterial::initialize(   dt,
                                                    numOwnedPoints,
                                                    ownedIDs,
                                                    neighborhoodList,
                                                    dataManager);

  CORRESPONDENCE::initializeCDPM2Model( m_unitFlag,
                                        m_youngsModulus,
                                        m_poissonsRatio,
                                        m_fc,
                                        m_ft,
                                        m_fracEnergy,
                                        m_hard1,
                                        m_hard2,
                                        m_hard3,
                                        m_hard4,
                                        m_hardqh2,
                                        m_soften,
                                        m_length);


  dataManager.getData(m_tempDeformationGradientFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_GreenLagrangeStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState1FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState2FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState3FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState4FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState5FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState6FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState7FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState8FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState9FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState10FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState11FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState12FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState13FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState14FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState15FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState16FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState17FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState18FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState19FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState20FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState21FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState22FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState23FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState24FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState25FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState26FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState27FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState28FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState29FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState30FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState31FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState32FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState33FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState34FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState35FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState36FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState37FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState38FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState39FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState40FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState41FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState42FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState43FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState44FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState45FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState46FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState47FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState48FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState49FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);


  dataManager.getData(m_tempDeformationGradientFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_GreenLagrangeStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState1FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState2FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState3FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState4FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState5FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState6FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState7FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState8FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState9FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState10FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState11FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState12FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState13FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState14FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState15FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState16FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState17FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState18FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState19FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState20FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState21FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState22FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState23FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState24FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState25FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState26FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState27FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState28FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState29FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState30FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState31FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState32FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState33FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState34FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState35FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState36FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState37FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState38FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState39FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState40FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState41FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState42FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState43FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState44FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState45FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState46FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState47FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState48FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState49FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

}

void
PeridigmNS::CDPM2CorrespondenceMaterial::computeCauchyStress(const double dt,
                                                            const int numOwnedPoints,
                                                            PeridigmNS::DataManager& dataManager) const
{
  // Compute the node-level stress values
  // This is only done for output (visualization) purposes
  double *unrotatedRateOfDeformation;
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  double *leftStretchN;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchN);
  
  double *rotationTensorN;
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  
  double* defGradN;
  dataManager.getData(m_tempDeformationGradientFieldId, PeridigmField::STEP_N)->ExtractView(&defGradN);
  
  double* GreenLagrangeN;
  dataManager.getData(m_GreenLagrangeStrainFieldId, PeridigmField::STEP_N)->ExtractView(&GreenLagrangeN);
  
  CORRESPONDENCE::setIdentityFullTensor(defGradN, numOwnedPoints);
  CORRESPONDENCE::setIdentityFullTensor(GreenLagrangeN, numOwnedPoints);

  for(int iID=0 ; iID<numOwnedPoints ; ++iID, defGradN+=9, GreenLagrangeN+=9, leftStretchN+=9, rotationTensorN+=9){
        CORRESPONDENCE::MatrixMultiply(false, false, 1.0, leftStretchN, rotationTensorN, defGradN);
        *(GreenLagrangeN) = 0.5 * ( *(defGradN) * *(defGradN) + *(defGradN+3) * *(defGradN+3) + *(defGradN+6) * *(defGradN+6) - 1.0 );
        *(GreenLagrangeN+1) = 0.5 * ( *(defGradN) * *(defGradN+1) + *(defGradN+3) * *(defGradN+4) + *(defGradN+6) * *(defGradN+7) );
        *(GreenLagrangeN+2) = 0.5 * ( *(defGradN) * *(defGradN+2) + *(defGradN+3) * *(defGradN+5) + *(defGradN+6) * *(defGradN+8) );
        *(GreenLagrangeN+3) = 0.5 * ( *(defGradN+1) * *(defGradN) + *(defGradN+4) * *(defGradN+3) + *(defGradN+7) * *(defGradN+6) );
        *(GreenLagrangeN+4) = 0.5 * ( *(defGradN+1) * *(defGradN+1) + *(defGradN+4) * *(defGradN+4) + *(defGradN+7) * *(defGradN+7) - 1.0 );
        *(GreenLagrangeN+5) = 0.5 * ( *(defGradN+1) * *(defGradN+2) + *(defGradN+4) * *(defGradN+5) + *(defGradN+7) * *(defGradN+8) );
        *(GreenLagrangeN+6) = 0.5 * ( *(defGradN+2) * *(defGradN) + *(defGradN+5) * *(defGradN+3) + *(defGradN+8) * *(defGradN+6) );
        *(GreenLagrangeN+7) = 0.5 * ( *(defGradN+2) * *(defGradN+1) + *(defGradN+5) * *(defGradN+4) + *(defGradN+8) * *(defGradN+7) );
        *(GreenLagrangeN+8) = 0.5 * ( *(defGradN+2) * *(defGradN+2) + *(defGradN+5) * *(defGradN+5) + *(defGradN+8) * *(defGradN+8) - 1.0 );
  }
  
  double *unrotatedCauchyStressN;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);

  double *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);

  double *state1N;
  double *state2N;
  double *state3N;
  double *state4N;
  double *state5N;
  double *state6N;
  double *state7N;
  double *state8N;
  double *state9N;
  double *state10N;
  double *state11N;
  double *state12N;
  double *state13N;
  double *state14N;
  double *state15N;
  double *state16N;
  double *state17N;
  double *state18N;
  double *state19N;
  double *state20N;
  double *state21N;
  double *state22N;
  double *state23N;
  double *state24N;
  double *state25N;
  double *state26N;
  double *state27N;
  double *state28N;
  double *state29N;
  double *state30N;
  double *state31N;
  double *state32N;
  double *state33N;
  double *state34N;
  double *state35N;
  double *state36N;
  double *state37N;
  double *state38N;
  double *state39N;
  double *state40N;
  double *state41N;
  double *state42N;
  double *state43N;
  double *state44N;
  double *state45N;
  double *state46N;
  double *state47N;
  double *state48N;
  double *state49N;
  dataManager.getData(m_microState1FieldId, PeridigmField::STEP_N)->ExtractView(&state1N);
  dataManager.getData(m_microState2FieldId, PeridigmField::STEP_N)->ExtractView(&state2N);
  dataManager.getData(m_microState3FieldId, PeridigmField::STEP_N)->ExtractView(&state3N);
  dataManager.getData(m_microState4FieldId, PeridigmField::STEP_N)->ExtractView(&state4N);
  dataManager.getData(m_microState5FieldId, PeridigmField::STEP_N)->ExtractView(&state5N);
  dataManager.getData(m_microState6FieldId, PeridigmField::STEP_N)->ExtractView(&state6N);
  dataManager.getData(m_microState7FieldId, PeridigmField::STEP_N)->ExtractView(&state7N);
  dataManager.getData(m_microState8FieldId, PeridigmField::STEP_N)->ExtractView(&state8N);
  dataManager.getData(m_microState9FieldId, PeridigmField::STEP_N)->ExtractView(&state9N);
  dataManager.getData(m_microState10FieldId, PeridigmField::STEP_N)->ExtractView(&state10N);
  dataManager.getData(m_microState11FieldId, PeridigmField::STEP_N)->ExtractView(&state11N);
  dataManager.getData(m_microState12FieldId, PeridigmField::STEP_N)->ExtractView(&state12N);
  dataManager.getData(m_microState13FieldId, PeridigmField::STEP_N)->ExtractView(&state13N);
  dataManager.getData(m_microState14FieldId, PeridigmField::STEP_N)->ExtractView(&state14N);
  dataManager.getData(m_microState15FieldId, PeridigmField::STEP_N)->ExtractView(&state15N);
  dataManager.getData(m_microState16FieldId, PeridigmField::STEP_N)->ExtractView(&state16N);
  dataManager.getData(m_microState17FieldId, PeridigmField::STEP_N)->ExtractView(&state17N);
  dataManager.getData(m_microState18FieldId, PeridigmField::STEP_N)->ExtractView(&state18N);
  dataManager.getData(m_microState19FieldId, PeridigmField::STEP_N)->ExtractView(&state19N);
  dataManager.getData(m_microState20FieldId, PeridigmField::STEP_N)->ExtractView(&state20N);
  dataManager.getData(m_microState21FieldId, PeridigmField::STEP_N)->ExtractView(&state21N);
  dataManager.getData(m_microState22FieldId, PeridigmField::STEP_N)->ExtractView(&state22N);
  dataManager.getData(m_microState23FieldId, PeridigmField::STEP_N)->ExtractView(&state23N);
  dataManager.getData(m_microState24FieldId, PeridigmField::STEP_N)->ExtractView(&state24N);
  dataManager.getData(m_microState25FieldId, PeridigmField::STEP_N)->ExtractView(&state25N);
  dataManager.getData(m_microState26FieldId, PeridigmField::STEP_N)->ExtractView(&state26N);
  dataManager.getData(m_microState27FieldId, PeridigmField::STEP_N)->ExtractView(&state27N);
  dataManager.getData(m_microState28FieldId, PeridigmField::STEP_N)->ExtractView(&state28N);
  dataManager.getData(m_microState29FieldId, PeridigmField::STEP_N)->ExtractView(&state29N);
  dataManager.getData(m_microState30FieldId, PeridigmField::STEP_N)->ExtractView(&state30N);
  dataManager.getData(m_microState31FieldId, PeridigmField::STEP_N)->ExtractView(&state31N);
  dataManager.getData(m_microState32FieldId, PeridigmField::STEP_N)->ExtractView(&state32N);
  dataManager.getData(m_microState33FieldId, PeridigmField::STEP_N)->ExtractView(&state33N);
  dataManager.getData(m_microState34FieldId, PeridigmField::STEP_N)->ExtractView(&state34N);
  dataManager.getData(m_microState35FieldId, PeridigmField::STEP_N)->ExtractView(&state35N);
  dataManager.getData(m_microState36FieldId, PeridigmField::STEP_N)->ExtractView(&state36N);
  dataManager.getData(m_microState37FieldId, PeridigmField::STEP_N)->ExtractView(&state37N);
  dataManager.getData(m_microState38FieldId, PeridigmField::STEP_N)->ExtractView(&state38N);
  dataManager.getData(m_microState39FieldId, PeridigmField::STEP_N)->ExtractView(&state39N);
  dataManager.getData(m_microState40FieldId, PeridigmField::STEP_N)->ExtractView(&state40N);
  dataManager.getData(m_microState41FieldId, PeridigmField::STEP_N)->ExtractView(&state41N);
  dataManager.getData(m_microState42FieldId, PeridigmField::STEP_N)->ExtractView(&state42N);
  dataManager.getData(m_microState43FieldId, PeridigmField::STEP_N)->ExtractView(&state43N);
  dataManager.getData(m_microState44FieldId, PeridigmField::STEP_N)->ExtractView(&state44N);
  dataManager.getData(m_microState45FieldId, PeridigmField::STEP_N)->ExtractView(&state45N);
  dataManager.getData(m_microState46FieldId, PeridigmField::STEP_N)->ExtractView(&state46N);
  dataManager.getData(m_microState47FieldId, PeridigmField::STEP_N)->ExtractView(&state47N);
  dataManager.getData(m_microState48FieldId, PeridigmField::STEP_N)->ExtractView(&state48N);
  dataManager.getData(m_microState49FieldId, PeridigmField::STEP_N)->ExtractView(&state49N);

  double *state1NP1;
  double *state2NP1;
  double *state3NP1;
  double *state4NP1;
  double *state5NP1;
  double *state6NP1;
  double *state7NP1;
  double *state8NP1;
  double *state9NP1;
  double *state10NP1;
  double *state11NP1;
  double *state12NP1;
  double *state13NP1;
  double *state14NP1;
  double *state15NP1;
  double *state16NP1;
  double *state17NP1;
  double *state18NP1;
  double *state19NP1;
  double *state20NP1;
  double *state21NP1;
  double *state22NP1;
  double *state23NP1;
  double *state24NP1;
  double *state25NP1;
  double *state26NP1;
  double *state27NP1;
  double *state28NP1;
  double *state29NP1;
  double *state30NP1;
  double *state31NP1;
  double *state32NP1;
  double *state33NP1;
  double *state34NP1;
  double *state35NP1;
  double *state36NP1;
  double *state37NP1;
  double *state38NP1;
  double *state39NP1;
  double *state40NP1;
  double *state41NP1;
  double *state42NP1;
  double *state43NP1;
  double *state44NP1;
  double *state45NP1;
  double *state46NP1;
  double *state47NP1;
  double *state48NP1;
  double *state49NP1;
  dataManager.getData(m_microState1FieldId, PeridigmField::STEP_NP1)->ExtractView(&state1NP1);
  dataManager.getData(m_microState2FieldId, PeridigmField::STEP_NP1)->ExtractView(&state2NP1);
  dataManager.getData(m_microState3FieldId, PeridigmField::STEP_NP1)->ExtractView(&state3NP1);
  dataManager.getData(m_microState4FieldId, PeridigmField::STEP_NP1)->ExtractView(&state4NP1);
  dataManager.getData(m_microState5FieldId, PeridigmField::STEP_NP1)->ExtractView(&state5NP1);
  dataManager.getData(m_microState6FieldId, PeridigmField::STEP_NP1)->ExtractView(&state6NP1);
  dataManager.getData(m_microState7FieldId, PeridigmField::STEP_NP1)->ExtractView(&state7NP1);
  dataManager.getData(m_microState8FieldId, PeridigmField::STEP_NP1)->ExtractView(&state8NP1);
  dataManager.getData(m_microState9FieldId, PeridigmField::STEP_NP1)->ExtractView(&state9NP1);
  dataManager.getData(m_microState10FieldId, PeridigmField::STEP_NP1)->ExtractView(&state10NP1);
  dataManager.getData(m_microState11FieldId, PeridigmField::STEP_NP1)->ExtractView(&state11NP1);
  dataManager.getData(m_microState12FieldId, PeridigmField::STEP_NP1)->ExtractView(&state12NP1);
  dataManager.getData(m_microState13FieldId, PeridigmField::STEP_NP1)->ExtractView(&state13NP1);
  dataManager.getData(m_microState14FieldId, PeridigmField::STEP_NP1)->ExtractView(&state14NP1);
  dataManager.getData(m_microState15FieldId, PeridigmField::STEP_NP1)->ExtractView(&state15NP1);
  dataManager.getData(m_microState16FieldId, PeridigmField::STEP_NP1)->ExtractView(&state16NP1);
  dataManager.getData(m_microState17FieldId, PeridigmField::STEP_NP1)->ExtractView(&state17NP1);
  dataManager.getData(m_microState18FieldId, PeridigmField::STEP_NP1)->ExtractView(&state18NP1);
  dataManager.getData(m_microState19FieldId, PeridigmField::STEP_NP1)->ExtractView(&state19NP1);
  dataManager.getData(m_microState20FieldId, PeridigmField::STEP_NP1)->ExtractView(&state20NP1);
  dataManager.getData(m_microState21FieldId, PeridigmField::STEP_NP1)->ExtractView(&state21NP1);
  dataManager.getData(m_microState22FieldId, PeridigmField::STEP_NP1)->ExtractView(&state22NP1);
  dataManager.getData(m_microState23FieldId, PeridigmField::STEP_NP1)->ExtractView(&state23NP1);
  dataManager.getData(m_microState24FieldId, PeridigmField::STEP_NP1)->ExtractView(&state24NP1);
  dataManager.getData(m_microState25FieldId, PeridigmField::STEP_NP1)->ExtractView(&state25NP1);
  dataManager.getData(m_microState26FieldId, PeridigmField::STEP_NP1)->ExtractView(&state26NP1);
  dataManager.getData(m_microState27FieldId, PeridigmField::STEP_NP1)->ExtractView(&state27NP1);
  dataManager.getData(m_microState28FieldId, PeridigmField::STEP_NP1)->ExtractView(&state28NP1);
  dataManager.getData(m_microState29FieldId, PeridigmField::STEP_NP1)->ExtractView(&state29NP1);
  dataManager.getData(m_microState30FieldId, PeridigmField::STEP_NP1)->ExtractView(&state30NP1);
  dataManager.getData(m_microState31FieldId, PeridigmField::STEP_NP1)->ExtractView(&state31NP1);
  dataManager.getData(m_microState32FieldId, PeridigmField::STEP_NP1)->ExtractView(&state32NP1);
  dataManager.getData(m_microState33FieldId, PeridigmField::STEP_NP1)->ExtractView(&state33NP1);
  dataManager.getData(m_microState34FieldId, PeridigmField::STEP_NP1)->ExtractView(&state34NP1);
  dataManager.getData(m_microState35FieldId, PeridigmField::STEP_NP1)->ExtractView(&state35NP1);
  dataManager.getData(m_microState36FieldId, PeridigmField::STEP_NP1)->ExtractView(&state36NP1);
  dataManager.getData(m_microState37FieldId, PeridigmField::STEP_NP1)->ExtractView(&state37NP1);
  dataManager.getData(m_microState38FieldId, PeridigmField::STEP_NP1)->ExtractView(&state38NP1);
  dataManager.getData(m_microState39FieldId, PeridigmField::STEP_NP1)->ExtractView(&state39NP1);
  dataManager.getData(m_microState40FieldId, PeridigmField::STEP_NP1)->ExtractView(&state40NP1);
  dataManager.getData(m_microState41FieldId, PeridigmField::STEP_NP1)->ExtractView(&state41NP1);
  dataManager.getData(m_microState42FieldId, PeridigmField::STEP_NP1)->ExtractView(&state42NP1);
  dataManager.getData(m_microState43FieldId, PeridigmField::STEP_NP1)->ExtractView(&state43NP1);
  dataManager.getData(m_microState44FieldId, PeridigmField::STEP_NP1)->ExtractView(&state44NP1);
  dataManager.getData(m_microState45FieldId, PeridigmField::STEP_NP1)->ExtractView(&state45NP1);
  dataManager.getData(m_microState46FieldId, PeridigmField::STEP_NP1)->ExtractView(&state46NP1);
  dataManager.getData(m_microState47FieldId, PeridigmField::STEP_NP1)->ExtractView(&state47NP1);
  dataManager.getData(m_microState48FieldId, PeridigmField::STEP_NP1)->ExtractView(&state48NP1);
  dataManager.getData(m_microState49FieldId, PeridigmField::STEP_NP1)->ExtractView(&state49NP1);

  double *internalEnergyN, *internalEnergyNP1;
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&internalEnergyN);
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&internalEnergyNP1);

  double *inelasticEnergyN, *inelasticEnergyNP1;
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&inelasticEnergyN);
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&inelasticEnergyNP1);

  CORRESPONDENCE::updateCDPM2Stress(unrotatedRateOfDeformation, 
                                   GreenLagrangeN, 
                                   unrotatedCauchyStressN, 
                                   state1N,
                                   state2N,
                                   state3N,
                                   state4N,
                                   state5N,
                                   state6N,
                                   state7N,
                                   state8N,
                                   state9N,
                                   state10N,
                                   state11N,
                                   state12N,
                                   state13N,
                                   state14N,
                                   state15N,
                                   state16N,
                                   state17N,
                                   state18N,
                                   state19N,
                                   state20N,
                                   state21N,
                                   state22N,
                                   state23N,
                                   state24N,
                                   state25N,
                                   state26N,
                                   state27N,
                                   state28N,
                                   state29N,
                                   state30N,
                                   state31N,
                                   state32N,
                                   state33N,
                                   state34N,
                                   state35N,
                                   state36N,
                                   state37N,
                                   state38N,
                                   state39N,
                                   state40N,
                                   state41N,
                                   state42N,
                                   state43N,
                                   state44N,
                                   state45N,
                                   state46N,
                                   state47N,
                                   state48N,
                                   state49N,
                                   internalEnergyN,
                                   inelasticEnergyN,
                                   unrotatedCauchyStressNP1,
                                   state1NP1,
                                   state2NP1,
                                   state3NP1,
                                   state4NP1,
                                   state5NP1,
                                   state6NP1,
                                   state7NP1,
                                   state8NP1,
                                   state9NP1,
                                   state10NP1,
                                   state11NP1,
                                   state12NP1,
                                   state13NP1,
                                   state14NP1,
                                   state15NP1,
                                   state16NP1,
                                   state17NP1,
                                   state18NP1,
                                   state19NP1,
                                   state20NP1,
                                   state21NP1,
                                   state22NP1,
                                   state23NP1,
                                   state24NP1,
                                   state25NP1,
                                   state26NP1,
                                   state27NP1,
                                   state28NP1,
                                   state29NP1,
                                   state30NP1,
                                   state31NP1,
                                   state32NP1,
                                   state33NP1,
                                   state34NP1,
                                   state35NP1,
                                   state36NP1,
                                   state37NP1,
                                   state38NP1,
                                   state39NP1,
                                   state40NP1,
                                   state41NP1,
                                   state42NP1,
                                   state43NP1,
                                   state44NP1,
                                   state45NP1,
                                   state46NP1,
                                   state47NP1,
                                   state48NP1,
                                   state49NP1,
                                   internalEnergyNP1,
                                   inelasticEnergyNP1,
                                   numOwnedPoints,
                                   dt);

  // No bond level stress
}

