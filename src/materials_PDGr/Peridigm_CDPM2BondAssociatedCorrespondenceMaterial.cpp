/*! \file Peridigm_CDPM2BondAssociatedCorrespondenceMaterial.cpp */

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

#include "Peridigm_CDPM2BondAssociatedCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "cdpm2_correspondence.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::CDPM2BondAssociatedCorrespondenceMaterial::CDPM2BondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : BondAssociatedCorrespondenceMaterial(params),
    m_unitFlag(0.0), m_youngsModulus(0.0), m_poissonsRatio(0.0), m_fc(0.0), m_ft(0.0), m_fracEnergy(0.0),
    m_hard1(0.0), m_hard2(0.0), m_hard3(0.0), m_hard4(0.0), m_hardqh2(0.0), m_soften(0.0), m_length(0.0),
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
    m_inelasticEnergyFieldId(-1), 
    m_bondLevelMicroState1FieldId(-1),
    m_bondLevelMicroState2FieldId(-1),
    m_bondLevelMicroState3FieldId(-1),
    m_bondLevelMicroState4FieldId(-1),
    m_bondLevelMicroState5FieldId(-1),
    m_bondLevelMicroState6FieldId(-1),
    m_bondLevelMicroState7FieldId(-1),
    m_bondLevelMicroState8FieldId(-1),
    m_bondLevelMicroState9FieldId(-1),
    m_bondLevelMicroState10FieldId(-1),
    m_bondLevelMicroState11FieldId(-1),
    m_bondLevelMicroState12FieldId(-1),
    m_bondLevelMicroState13FieldId(-1),
    m_bondLevelMicroState14FieldId(-1),
    m_bondLevelMicroState15FieldId(-1),
    m_bondLevelMicroState16FieldId(-1),
    m_bondLevelMicroState17FieldId(-1),
    m_bondLevelMicroState18FieldId(-1),
    m_bondLevelMicroState19FieldId(-1),
    m_bondLevelMicroState20FieldId(-1),
    m_bondLevelMicroState21FieldId(-1),
    m_bondLevelMicroState22FieldId(-1),
    m_bondLevelMicroState23FieldId(-1),
    m_bondLevelMicroState24FieldId(-1),
    m_bondLevelMicroState25FieldId(-1),
    m_bondLevelMicroState26FieldId(-1),
    m_bondLevelMicroState27FieldId(-1),
    m_bondLevelMicroState28FieldId(-1),
    m_bondLevelMicroState29FieldId(-1),
    m_bondLevelMicroState30FieldId(-1),
    m_bondLevelMicroState31FieldId(-1),
    m_bondLevelMicroState32FieldId(-1),
    m_bondLevelMicroState33FieldId(-1),
    m_bondLevelMicroState34FieldId(-1),
    m_bondLevelMicroState35FieldId(-1),
    m_bondLevelMicroState36FieldId(-1),
    m_bondLevelMicroState37FieldId(-1),
    m_bondLevelMicroState38FieldId(-1),
    m_bondLevelMicroState39FieldId(-1),
    m_bondLevelMicroState40FieldId(-1),
    m_bondLevelMicroState41FieldId(-1),
    m_bondLevelMicroState42FieldId(-1),
    m_bondLevelMicroState43FieldId(-1),
    m_bondLevelMicroState44FieldId(-1),
    m_bondLevelMicroState45FieldId(-1),
    m_bondLevelMicroState46FieldId(-1),
    m_bondLevelMicroState47FieldId(-1),
    m_bondLevelMicroState48FieldId(-1),
    m_bondLevelMicroState49FieldId(-1),
    m_bondLevelInternalEnergyFieldId(-1),
    m_bondLevelInelasticEnergyFieldId(-1)
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

  m_bondLevelMicroState1FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_1");
  m_bondLevelMicroState2FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_2");
  m_bondLevelMicroState3FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_3");
  m_bondLevelMicroState4FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_4");
  m_bondLevelMicroState5FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_5");
  m_bondLevelMicroState6FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_6");
  m_bondLevelMicroState7FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_7");
  m_bondLevelMicroState8FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_8");
  m_bondLevelMicroState9FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_9");
  m_bondLevelMicroState10FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_10");
  m_bondLevelMicroState11FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_11");
  m_bondLevelMicroState12FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_12");
  m_bondLevelMicroState13FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_13");
  m_bondLevelMicroState14FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_14");
  m_bondLevelMicroState15FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_15");
  m_bondLevelMicroState16FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_16");
  m_bondLevelMicroState17FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_17");
  m_bondLevelMicroState18FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_18");
  m_bondLevelMicroState19FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_19");
  m_bondLevelMicroState20FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_20");
  m_bondLevelMicroState21FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_21");
  m_bondLevelMicroState22FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_22");
  m_bondLevelMicroState23FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_23");
  m_bondLevelMicroState24FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_24");
  m_bondLevelMicroState25FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_25");
  m_bondLevelMicroState26FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_26");
  m_bondLevelMicroState27FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_27");
  m_bondLevelMicroState28FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_28");
  m_bondLevelMicroState29FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_29");
  m_bondLevelMicroState30FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_30");
  m_bondLevelMicroState31FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_31");
  m_bondLevelMicroState32FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_32");
  m_bondLevelMicroState33FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_33");
  m_bondLevelMicroState34FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_34");
  m_bondLevelMicroState35FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_35");
  m_bondLevelMicroState36FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_36");
  m_bondLevelMicroState37FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_37");
  m_bondLevelMicroState38FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_38");
  m_bondLevelMicroState39FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_39");
  m_bondLevelMicroState40FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_40");
  m_bondLevelMicroState41FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_41");
  m_bondLevelMicroState42FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_42");
  m_bondLevelMicroState43FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_43");
  m_bondLevelMicroState44FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_44");
  m_bondLevelMicroState45FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_45");
  m_bondLevelMicroState46FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_46");
  m_bondLevelMicroState47FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_47");
  m_bondLevelMicroState48FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_48");
  m_bondLevelMicroState49FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_49");
  m_bondLevelInternalEnergyFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Internal_Energy");
  m_bondLevelInelasticEnergyFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Inelastic_Energy");

  m_fieldIds.push_back(m_bondLevelMicroState1FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState2FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState3FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState4FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState5FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState6FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState7FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState8FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState9FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState10FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState11FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState12FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState13FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState14FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState15FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState16FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState17FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState18FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState19FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState20FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState21FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState22FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState23FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState24FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState25FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState26FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState27FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState28FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState29FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState30FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState31FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState32FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState33FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState34FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState35FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState36FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState37FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState38FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState39FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState40FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState41FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState42FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState43FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState44FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState45FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState46FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState47FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState48FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState49FieldId);
  m_fieldIds.push_back(m_bondLevelInternalEnergyFieldId);
  m_fieldIds.push_back(m_bondLevelInelasticEnergyFieldId);
}

PeridigmNS::CDPM2BondAssociatedCorrespondenceMaterial::~CDPM2BondAssociatedCorrespondenceMaterial()
{
}

void
PeridigmNS::CDPM2BondAssociatedCorrespondenceMaterial::initialize(const double dt,
                                                                       const int numOwnedPoints,
                                                                       const int* ownedIDs,
                                                                       const int* neighborhoodList,
                                                                       PeridigmNS::DataManager& dataManager)
{

  PeridigmNS::BondAssociatedCorrespondenceMaterial::initialize(dt,
                                                               numOwnedPoints,
                                                               ownedIDs,
                                                               neighborhoodList,
                                                               dataManager);

  CORRESPONDENCE::initializeCDPM2Model(  m_unitFlag,
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

  dataManager.getData(m_bondLevelMicroState1FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState2FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState3FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState4FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState5FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState6FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState7FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState8FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState9FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState10FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState11FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState12FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState13FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState14FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState15FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState16FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState17FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState18FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState19FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState20FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState21FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState22FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState23FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState24FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState25FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState26FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState27FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState28FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState29FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState30FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState31FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState32FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState33FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState34FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState35FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState36FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState37FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState38FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState39FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState40FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState41FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState42FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState43FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState44FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState45FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState46FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState47FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState48FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState49FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);

  dataManager.getData(m_bondLevelMicroState1FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState2FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState3FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState4FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState5FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState6FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState7FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState8FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState9FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState10FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState11FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState12FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState13FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState14FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState15FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState16FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState17FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState18FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState19FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState20FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState21FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState22FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState23FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState24FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState25FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState26FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState27FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState28FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState29FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState30FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState31FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState32FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState33FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState34FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState35FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState36FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState37FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState38FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState39FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState40FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState41FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState42FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState43FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState44FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState45FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState46FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState47FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState48FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState49FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}

void
PeridigmNS::CDPM2BondAssociatedCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                                                const int numOwnedPoints,
                                                                                const int* neighborhoodList,
                                                                                PeridigmNS::DataManager& dataManager) const
{

}

void
PeridigmNS::CDPM2BondAssociatedCorrespondenceMaterial::computePK2Stress(const double dt,
                                                                             const int numOwnedPoints,
                                                                             const int* neighborhoodList,
                                                                             PeridigmNS::DataManager& dataManager) const
{

  // Compute the node-level stress values
  // This is only done for output (visualization) purposes
  double *strainRate;
  dataManager.getData(m_strainRateFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainRate);

  double *strainN;
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_N)->ExtractView(&strainN);

  double *PK2StressN;
  dataManager.getData(m_PK2StressFieldId, PeridigmField::STEP_N)->ExtractView(&PK2StressN);

  double *PK2StressNP1;
  dataManager.getData(m_PK2StressFieldId, PeridigmField::STEP_NP1)->ExtractView(&PK2StressNP1);

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

  CORRESPONDENCE::updateCDPM2Stress(strainRate, 
                                           strainN, 
                                           PK2StressN, 
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
                                           PK2StressNP1,
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

  // Compute the bond-level stress values
  double *bondLevelStrainRateXX, *bondLevelStrainRateXY, *bondLevelStrainRateXZ;
  double *bondLevelStrainRateYX, *bondLevelStrainRateYY, *bondLevelStrainRateYZ;
  double *bondLevelStrainRateZX, *bondLevelStrainRateZY, *bondLevelStrainRateZZ;
  dataManager.getData(m_bondLevelStrainRateXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateXX);
  dataManager.getData(m_bondLevelStrainRateXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateXY);
  dataManager.getData(m_bondLevelStrainRateXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateXZ);
  dataManager.getData(m_bondLevelStrainRateYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateYX);
  dataManager.getData(m_bondLevelStrainRateYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateYY);
  dataManager.getData(m_bondLevelStrainRateYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateYZ);
  dataManager.getData(m_bondLevelStrainRateZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateZX);
  dataManager.getData(m_bondLevelStrainRateZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateZY);
  dataManager.getData(m_bondLevelStrainRateZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStrainRateZZ);

  double *bondLevelStrainXXN, *bondLevelStrainXYN, *bondLevelStrainXZN;
  double *bondLevelStrainYXN, *bondLevelStrainYYN, *bondLevelStrainYZN;
  double *bondLevelStrainZXN, *bondLevelStrainZYN, *bondLevelStrainZZN;
  dataManager.getData(m_bondLevelStrainXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainXXN);
  dataManager.getData(m_bondLevelStrainXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainXYN);
  dataManager.getData(m_bondLevelStrainXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainXZN);
  dataManager.getData(m_bondLevelStrainYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainYXN);
  dataManager.getData(m_bondLevelStrainYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainYYN);
  dataManager.getData(m_bondLevelStrainYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainYZN);
  dataManager.getData(m_bondLevelStrainZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainZXN);
  dataManager.getData(m_bondLevelStrainZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainZYN);
  dataManager.getData(m_bondLevelStrainZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelStrainZZN);

  double *bondLevelPK2StressXXN, *bondLevelPK2StressXYN, *bondLevelPK2StressXZN;
  double *bondLevelPK2StressYXN, *bondLevelPK2StressYYN, *bondLevelPK2StressYZN;
  double *bondLevelPK2StressZXN, *bondLevelPK2StressZYN, *bondLevelPK2StressZZN;
  dataManager.getData(m_bondLevelPK2StressXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressXXN);
  dataManager.getData(m_bondLevelPK2StressXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressXYN);
  dataManager.getData(m_bondLevelPK2StressXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressXZN);
  dataManager.getData(m_bondLevelPK2StressYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressYXN);
  dataManager.getData(m_bondLevelPK2StressYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressYYN);
  dataManager.getData(m_bondLevelPK2StressYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressYZN);
  dataManager.getData(m_bondLevelPK2StressZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressZXN);
  dataManager.getData(m_bondLevelPK2StressZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressZYN);
  dataManager.getData(m_bondLevelPK2StressZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelPK2StressZZN);

  double *bondLevelPK2StressXXNP1, *bondLevelPK2StressXYNP1, *bondLevelPK2StressXZNP1;
  double *bondLevelPK2StressYXNP1, *bondLevelPK2StressYYNP1, *bondLevelPK2StressYZNP1;
  double *bondLevelPK2StressZXNP1, *bondLevelPK2StressZYNP1, *bondLevelPK2StressZZNP1;
  dataManager.getData(m_bondLevelPK2StressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressXXNP1);
  dataManager.getData(m_bondLevelPK2StressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressXYNP1);
  dataManager.getData(m_bondLevelPK2StressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressXZNP1);
  dataManager.getData(m_bondLevelPK2StressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressYXNP1);
  dataManager.getData(m_bondLevelPK2StressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressYYNP1);
  dataManager.getData(m_bondLevelPK2StressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressYZNP1);
  dataManager.getData(m_bondLevelPK2StressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressZXNP1);
  dataManager.getData(m_bondLevelPK2StressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressZYNP1);
  dataManager.getData(m_bondLevelPK2StressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelPK2StressZZNP1);

  double *bondLevelState1N;
  double *bondLevelState2N;
  double *bondLevelState3N;
  double *bondLevelState4N;
  double *bondLevelState5N;
  double *bondLevelState6N;
  double *bondLevelState7N;
  double *bondLevelState8N;
  double *bondLevelState9N;
  double *bondLevelState10N;
  double *bondLevelState11N;
  double *bondLevelState12N;
  double *bondLevelState13N;
  double *bondLevelState14N;
  double *bondLevelState15N;
  double *bondLevelState16N;
  double *bondLevelState17N;
  double *bondLevelState18N;
  double *bondLevelState19N;
  double *bondLevelState20N;
  double *bondLevelState21N;
  double *bondLevelState22N;
  double *bondLevelState23N;
  double *bondLevelState24N;
  double *bondLevelState25N;
  double *bondLevelState26N;
  double *bondLevelState27N;
  double *bondLevelState28N;
  double *bondLevelState29N;
  double *bondLevelState30N;
  double *bondLevelState31N;
  double *bondLevelState32N;
  double *bondLevelState33N;
  double *bondLevelState34N;
  double *bondLevelState35N;
  double *bondLevelState36N;
  double *bondLevelState37N;
  double *bondLevelState38N;
  double *bondLevelState39N;
  double *bondLevelState40N;
  double *bondLevelState41N;
  double *bondLevelState42N;
  double *bondLevelState43N;
  double *bondLevelState44N;
  double *bondLevelState45N;
  double *bondLevelState46N;
  double *bondLevelState47N;
  double *bondLevelState48N;
  double *bondLevelState49N;
  dataManager.getData(m_bondLevelMicroState1FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState1N);
  dataManager.getData(m_bondLevelMicroState2FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState2N);
  dataManager.getData(m_bondLevelMicroState3FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState3N);
  dataManager.getData(m_bondLevelMicroState4FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState4N);
  dataManager.getData(m_bondLevelMicroState5FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState5N);
  dataManager.getData(m_bondLevelMicroState6FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState6N);
  dataManager.getData(m_bondLevelMicroState7FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState7N);
  dataManager.getData(m_bondLevelMicroState8FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState8N);
  dataManager.getData(m_bondLevelMicroState9FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState9N);
  dataManager.getData(m_bondLevelMicroState10FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState10N);
  dataManager.getData(m_bondLevelMicroState11FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState11N);
  dataManager.getData(m_bondLevelMicroState12FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState12N);
  dataManager.getData(m_bondLevelMicroState13FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState13N);
  dataManager.getData(m_bondLevelMicroState14FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState14N);
  dataManager.getData(m_bondLevelMicroState15FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState15N);
  dataManager.getData(m_bondLevelMicroState16FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState16N);
  dataManager.getData(m_bondLevelMicroState17FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState17N);
  dataManager.getData(m_bondLevelMicroState18FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState18N);
  dataManager.getData(m_bondLevelMicroState19FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState19N);
  dataManager.getData(m_bondLevelMicroState20FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState20N);
  dataManager.getData(m_bondLevelMicroState21FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState21N);
  dataManager.getData(m_bondLevelMicroState22FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState22N);
  dataManager.getData(m_bondLevelMicroState23FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState23N);
  dataManager.getData(m_bondLevelMicroState24FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState24N);
  dataManager.getData(m_bondLevelMicroState25FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState25N);
  dataManager.getData(m_bondLevelMicroState26FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState26N);
  dataManager.getData(m_bondLevelMicroState27FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState27N);
  dataManager.getData(m_bondLevelMicroState28FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState28N);
  dataManager.getData(m_bondLevelMicroState29FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState29N);
  dataManager.getData(m_bondLevelMicroState30FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState30N);
  dataManager.getData(m_bondLevelMicroState31FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState31N);
  dataManager.getData(m_bondLevelMicroState32FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState32N);
  dataManager.getData(m_bondLevelMicroState33FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState33N);
  dataManager.getData(m_bondLevelMicroState34FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState34N);
  dataManager.getData(m_bondLevelMicroState35FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState35N);
  dataManager.getData(m_bondLevelMicroState36FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState36N);
  dataManager.getData(m_bondLevelMicroState37FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState37N);
  dataManager.getData(m_bondLevelMicroState38FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState38N);
  dataManager.getData(m_bondLevelMicroState39FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState39N);
  dataManager.getData(m_bondLevelMicroState40FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState40N);
  dataManager.getData(m_bondLevelMicroState41FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState41N);
  dataManager.getData(m_bondLevelMicroState42FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState42N);
  dataManager.getData(m_bondLevelMicroState43FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState43N);
  dataManager.getData(m_bondLevelMicroState44FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState44N);
  dataManager.getData(m_bondLevelMicroState45FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState45N);
  dataManager.getData(m_bondLevelMicroState46FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState46N);
  dataManager.getData(m_bondLevelMicroState47FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState47N);
  dataManager.getData(m_bondLevelMicroState48FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState48N);
  dataManager.getData(m_bondLevelMicroState49FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState49N);

  double *bondLevelState1NP1;
  double *bondLevelState2NP1;
  double *bondLevelState3NP1;
  double *bondLevelState4NP1;
  double *bondLevelState5NP1;
  double *bondLevelState6NP1;
  double *bondLevelState7NP1;
  double *bondLevelState8NP1;
  double *bondLevelState9NP1;
  double *bondLevelState10NP1;
  double *bondLevelState11NP1;
  double *bondLevelState12NP1;
  double *bondLevelState13NP1;
  double *bondLevelState14NP1;
  double *bondLevelState15NP1;
  double *bondLevelState16NP1;
  double *bondLevelState17NP1;
  double *bondLevelState18NP1;
  double *bondLevelState19NP1;
  double *bondLevelState20NP1;
  double *bondLevelState21NP1;
  double *bondLevelState22NP1;
  double *bondLevelState23NP1;
  double *bondLevelState24NP1;
  double *bondLevelState25NP1;
  double *bondLevelState26NP1;
  double *bondLevelState27NP1;
  double *bondLevelState28NP1;
  double *bondLevelState29NP1;
  double *bondLevelState30NP1;
  double *bondLevelState31NP1;
  double *bondLevelState32NP1;
  double *bondLevelState33NP1;
  double *bondLevelState34NP1;
  double *bondLevelState35NP1;
  double *bondLevelState36NP1;
  double *bondLevelState37NP1;
  double *bondLevelState38NP1;
  double *bondLevelState39NP1;
  double *bondLevelState40NP1;
  double *bondLevelState41NP1;
  double *bondLevelState42NP1;
  double *bondLevelState43NP1;
  double *bondLevelState44NP1;
  double *bondLevelState45NP1;
  double *bondLevelState46NP1;
  double *bondLevelState47NP1;
  double *bondLevelState48NP1;
  double *bondLevelState49NP1;
  dataManager.getData(m_bondLevelMicroState1FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState1NP1);
  dataManager.getData(m_bondLevelMicroState2FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState2NP1);
  dataManager.getData(m_bondLevelMicroState3FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState3NP1);
  dataManager.getData(m_bondLevelMicroState4FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState4NP1);
  dataManager.getData(m_bondLevelMicroState5FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState5NP1);
  dataManager.getData(m_bondLevelMicroState6FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState6NP1);
  dataManager.getData(m_bondLevelMicroState7FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState7NP1);
  dataManager.getData(m_bondLevelMicroState8FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState8NP1);
  dataManager.getData(m_bondLevelMicroState9FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState9NP1);
  dataManager.getData(m_bondLevelMicroState10FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState10NP1);
  dataManager.getData(m_bondLevelMicroState11FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState11NP1);
  dataManager.getData(m_bondLevelMicroState12FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState12NP1);
  dataManager.getData(m_bondLevelMicroState13FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState13NP1);
  dataManager.getData(m_bondLevelMicroState14FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState14NP1);
  dataManager.getData(m_bondLevelMicroState15FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState15NP1);
  dataManager.getData(m_bondLevelMicroState16FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState16NP1);
  dataManager.getData(m_bondLevelMicroState17FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState17NP1);
  dataManager.getData(m_bondLevelMicroState18FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState18NP1);
  dataManager.getData(m_bondLevelMicroState19FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState19NP1);
  dataManager.getData(m_bondLevelMicroState20FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState20NP1);
  dataManager.getData(m_bondLevelMicroState21FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState21NP1);
  dataManager.getData(m_bondLevelMicroState22FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState22NP1);
  dataManager.getData(m_bondLevelMicroState23FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState23NP1);
  dataManager.getData(m_bondLevelMicroState24FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState24NP1);
  dataManager.getData(m_bondLevelMicroState25FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState25NP1);
  dataManager.getData(m_bondLevelMicroState26FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState26NP1);
  dataManager.getData(m_bondLevelMicroState27FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState27NP1);
  dataManager.getData(m_bondLevelMicroState28FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState28NP1);
  dataManager.getData(m_bondLevelMicroState29FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState29NP1);
  dataManager.getData(m_bondLevelMicroState30FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState30NP1);
  dataManager.getData(m_bondLevelMicroState31FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState31NP1);
  dataManager.getData(m_bondLevelMicroState32FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState32NP1);
  dataManager.getData(m_bondLevelMicroState33FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState33NP1);
  dataManager.getData(m_bondLevelMicroState34FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState34NP1);
  dataManager.getData(m_bondLevelMicroState35FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState35NP1);
  dataManager.getData(m_bondLevelMicroState36FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState36NP1);
  dataManager.getData(m_bondLevelMicroState37FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState37NP1);
  dataManager.getData(m_bondLevelMicroState38FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState38NP1);
  dataManager.getData(m_bondLevelMicroState39FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState39NP1);
  dataManager.getData(m_bondLevelMicroState40FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState40NP1);
  dataManager.getData(m_bondLevelMicroState41FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState41NP1);
  dataManager.getData(m_bondLevelMicroState42FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState42NP1);
  dataManager.getData(m_bondLevelMicroState43FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState43NP1);
  dataManager.getData(m_bondLevelMicroState44FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState44NP1);
  dataManager.getData(m_bondLevelMicroState45FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState45NP1);
  dataManager.getData(m_bondLevelMicroState46FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState46NP1);
  dataManager.getData(m_bondLevelMicroState47FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState47NP1);
  dataManager.getData(m_bondLevelMicroState48FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState48NP1);
  dataManager.getData(m_bondLevelMicroState49FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState49NP1);

  double *bondLevelInternalEnergyN, *bondLevelInternalEnergyNP1;
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelInternalEnergyN);
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelInternalEnergyNP1);

  double *bondLevelInelasticEnergyN, *bondLevelInelasticEnergyNP1;
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelInelasticEnergyN);
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelInelasticEnergyNP1);

  CORRESPONDENCE::updateBondLevelCDPM2Stress(bondLevelStrainRateXX,
                                                    bondLevelStrainRateXY,
                                                    bondLevelStrainRateXZ,
                                                    bondLevelStrainRateYX,
                                                    bondLevelStrainRateYY,
                                                    bondLevelStrainRateYZ,
                                                    bondLevelStrainRateZX,
                                                    bondLevelStrainRateZY,
                                                    bondLevelStrainRateZZ,
                                                    bondLevelStrainXXN,
                                                    bondLevelStrainXYN,
                                                    bondLevelStrainXZN,
                                                    bondLevelStrainYXN,
                                                    bondLevelStrainYYN,
                                                    bondLevelStrainYZN,
                                                    bondLevelStrainZXN,
                                                    bondLevelStrainZYN,
                                                    bondLevelStrainZZN,
                                                    bondLevelPK2StressXXN, 
                                                    bondLevelPK2StressXYN, 
                                                    bondLevelPK2StressXZN, 
                                                    bondLevelPK2StressYXN, 
                                                    bondLevelPK2StressYYN, 
                                                    bondLevelPK2StressYZN, 
                                                    bondLevelPK2StressZXN, 
                                                    bondLevelPK2StressZYN, 
                                                    bondLevelPK2StressZZN, 
                                                    bondLevelState1N,
                                                    bondLevelState2N,
                                                    bondLevelState3N,
                                                    bondLevelState4N,
                                                    bondLevelState5N,
                                                    bondLevelState6N,
                                                    bondLevelState7N,
                                                    bondLevelState8N,
                                                    bondLevelState9N,
                                                    bondLevelState10N,
                                                    bondLevelState11N,
                                                    bondLevelState12N,
                                                    bondLevelState13N,
                                                    bondLevelState14N,
                                                    bondLevelState15N,
                                                    bondLevelState16N,
                                                    bondLevelState17N,
                                                    bondLevelState18N,
                                                    bondLevelState19N,
                                                    bondLevelState20N,
                                                    bondLevelState21N,
                                                    bondLevelState22N,
                                                    bondLevelState23N,
                                                    bondLevelState24N,
                                                    bondLevelState25N,
                                                    bondLevelState26N,
                                                    bondLevelState27N,
                                                    bondLevelState28N,
                                                    bondLevelState29N,
                                                    bondLevelState30N,
                                                    bondLevelState31N,
                                                    bondLevelState32N,
                                                    bondLevelState33N,
                                                    bondLevelState34N,
                                                    bondLevelState35N,
                                                    bondLevelState36N,
                                                    bondLevelState37N,
                                                    bondLevelState38N,
                                                    bondLevelState39N,
                                                    bondLevelState40N,
                                                    bondLevelState41N,
                                                    bondLevelState42N,
                                                    bondLevelState43N,
                                                    bondLevelState44N,
                                                    bondLevelState45N,
                                                    bondLevelState46N,
                                                    bondLevelState47N,
                                                    bondLevelState48N,
                                                    bondLevelState49N,
                                                    bondLevelInternalEnergyN,
                                                    bondLevelInelasticEnergyN,
                                                    bondLevelPK2StressXXNP1, 
                                                    bondLevelPK2StressXYNP1, 
                                                    bondLevelPK2StressXZNP1, 
                                                    bondLevelPK2StressYXNP1, 
                                                    bondLevelPK2StressYYNP1, 
                                                    bondLevelPK2StressYZNP1, 
                                                    bondLevelPK2StressZXNP1, 
                                                    bondLevelPK2StressZYNP1, 
                                                    bondLevelPK2StressZZNP1, 
                                                    bondLevelState1NP1,
                                                    bondLevelState2NP1,
                                                    bondLevelState3NP1,
                                                    bondLevelState4NP1,
                                                    bondLevelState5NP1,
                                                    bondLevelState6NP1,
                                                    bondLevelState7NP1,
                                                    bondLevelState8NP1,
                                                    bondLevelState9NP1,
                                                    bondLevelState10NP1,
                                                    bondLevelState11NP1,
                                                    bondLevelState12NP1,
                                                    bondLevelState13NP1,
                                                    bondLevelState14NP1,
                                                    bondLevelState15NP1,
                                                    bondLevelState16NP1,
                                                    bondLevelState17NP1,
                                                    bondLevelState18NP1,
                                                    bondLevelState19NP1,
                                                    bondLevelState20NP1,
                                                    bondLevelState21NP1,
                                                    bondLevelState22NP1,
                                                    bondLevelState23NP1,
                                                    bondLevelState24NP1,
                                                    bondLevelState25NP1,
                                                    bondLevelState26NP1,
                                                    bondLevelState27NP1,
                                                    bondLevelState28NP1,
                                                    bondLevelState29NP1,
                                                    bondLevelState30NP1,
                                                    bondLevelState31NP1,
                                                    bondLevelState32NP1,
                                                    bondLevelState33NP1,
                                                    bondLevelState34NP1,
                                                    bondLevelState35NP1,
                                                    bondLevelState36NP1,
                                                    bondLevelState37NP1,
                                                    bondLevelState38NP1,
                                                    bondLevelState39NP1,
                                                    bondLevelState40NP1,
                                                    bondLevelState41NP1,
                                                    bondLevelState42NP1,
                                                    bondLevelState43NP1,
                                                    bondLevelState44NP1,
                                                    bondLevelState45NP1,
                                                    bondLevelState46NP1,
                                                    bondLevelState47NP1,
                                                    bondLevelState48NP1,
                                                    bondLevelState49NP1,
                                                    bondLevelInternalEnergyNP1,
                                                    bondLevelInelasticEnergyNP1,
                                                    neighborhoodList,
                                                    numOwnedPoints,
                                                    dt);
}
