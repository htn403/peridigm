/*! \file Peridigm_MicroplaneBondAssociatedCorrespondenceMaterial.cpp */

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

#include "Peridigm_MicroplaneBondAssociatedCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "microplane_correspondence.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::MicroplaneBondAssociatedCorrespondenceMaterial::MicroplaneBondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : BondAssociatedCorrespondenceMaterial(params),
    m_youngsModulus(0.0), m_poissonsRatio(0.0),
    m_k1(0.0), m_k2(0.0), m_k3(0.0), m_k4(0.0), 
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
    m_microState50FieldId(-1),
    m_microState51FieldId(-1),
    m_microState52FieldId(-1),
    m_microState53FieldId(-1),
    m_microState54FieldId(-1),
    m_microState55FieldId(-1),
    m_microState56FieldId(-1),
    m_microState57FieldId(-1),
    m_microState58FieldId(-1),
    m_microState59FieldId(-1),
    m_microState60FieldId(-1),
    m_microState61FieldId(-1),
    m_microState62FieldId(-1),
    m_microState63FieldId(-1),
    m_microState64FieldId(-1),
    m_microState65FieldId(-1),
    m_microState66FieldId(-1),
    m_microState67FieldId(-1),
    m_microState68FieldId(-1),
    m_microState69FieldId(-1),
    m_microState70FieldId(-1),
    m_microState71FieldId(-1),
    m_microState72FieldId(-1),
    m_microState73FieldId(-1),
    m_microState74FieldId(-1),
    m_microState75FieldId(-1),
    m_microState76FieldId(-1),
    m_microState77FieldId(-1),
    m_microState78FieldId(-1),
    m_microState79FieldId(-1),
    m_microState80FieldId(-1),
    m_microState81FieldId(-1),
    m_microState82FieldId(-1),
    m_microState83FieldId(-1),
    m_microState84FieldId(-1),
    m_microState85FieldId(-1),
    m_microState86FieldId(-1),
    m_microState87FieldId(-1),
    m_microState88FieldId(-1),
    m_microState89FieldId(-1),
    m_microState90FieldId(-1),
    m_microState91FieldId(-1),
    m_microState92FieldId(-1),
    m_microState93FieldId(-1),
    m_microState94FieldId(-1),
    m_microState95FieldId(-1),
    m_microState96FieldId(-1),
    m_microState97FieldId(-1),
    m_microState98FieldId(-1),
    m_microState99FieldId(-1),
    m_microState100FieldId(-1),
    m_microState101FieldId(-1),
    m_microState102FieldId(-1),
    m_microState103FieldId(-1),
    m_microState104FieldId(-1),
    m_microState105FieldId(-1),
    m_microState106FieldId(-1),
    m_microState107FieldId(-1),
    m_microState108FieldId(-1),
    m_microState109FieldId(-1),
    m_microState110FieldId(-1),
    m_microState111FieldId(-1),
    m_microState112FieldId(-1),
    m_microState113FieldId(-1),
    m_microState114FieldId(-1),
    m_microState115FieldId(-1),
    m_microState116FieldId(-1),
    m_microState117FieldId(-1),
    m_microState118FieldId(-1),
    m_microState119FieldId(-1),
    m_microState120FieldId(-1),
    m_microState121FieldId(-1),
    m_microState122FieldId(-1),
    m_microState123FieldId(-1),
    m_microState124FieldId(-1),
    m_microState125FieldId(-1),
    m_microState126FieldId(-1),
    m_microState127FieldId(-1),
    m_microState128FieldId(-1),
    m_microState129FieldId(-1),
    m_microState130FieldId(-1),
    m_microState131FieldId(-1),
    m_microState132FieldId(-1),
    m_microState133FieldId(-1),
    m_microState134FieldId(-1),
    m_microState135FieldId(-1),
    m_microState136FieldId(-1),
    m_microState137FieldId(-1),
    m_microState138FieldId(-1),
    m_microState139FieldId(-1),
    m_microState140FieldId(-1),
    m_microState141FieldId(-1),
    m_microState142FieldId(-1),
    m_microState143FieldId(-1),
    m_microState144FieldId(-1),
    m_microState145FieldId(-1),
    m_microState146FieldId(-1),
    m_microState147FieldId(-1),
    m_microState148FieldId(-1),
    m_microState149FieldId(-1),
    m_microState150FieldId(-1),
    m_microState151FieldId(-1),
    m_microState152FieldId(-1),
    m_microState153FieldId(-1),
    m_microState154FieldId(-1),
    m_microState155FieldId(-1),
    m_microState156FieldId(-1),
    m_microState157FieldId(-1),
    m_microState158FieldId(-1),
    m_microState159FieldId(-1),
    m_microState160FieldId(-1),
    m_microState161FieldId(-1),
    m_microState162FieldId(-1),
    m_microState163FieldId(-1),
    m_microState164FieldId(-1),
    m_microState165FieldId(-1),
    m_microState166FieldId(-1),
    m_microState167FieldId(-1),
    m_microState168FieldId(-1),
    m_microState169FieldId(-1),
    m_microState170FieldId(-1),
    m_microState171FieldId(-1),
    m_microState172FieldId(-1),
    m_microState173FieldId(-1),
    m_microState174FieldId(-1),
    m_microState175FieldId(-1),
    m_microState176FieldId(-1),
    m_microState177FieldId(-1),
    m_microState178FieldId(-1),
    m_microState179FieldId(-1),
    m_microState180FieldId(-1),
    m_microState181FieldId(-1),
    m_microState182FieldId(-1),
    m_microState183FieldId(-1),
    m_microState184FieldId(-1),
    m_microState185FieldId(-1),
    m_microState186FieldId(-1),
    m_microState187FieldId(-1),
    m_microState188FieldId(-1),
    m_microState189FieldId(-1),
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
    m_bondLevelMicroState50FieldId(-1),
    m_bondLevelMicroState51FieldId(-1),
    m_bondLevelMicroState52FieldId(-1),
    m_bondLevelMicroState53FieldId(-1),
    m_bondLevelMicroState54FieldId(-1),
    m_bondLevelMicroState55FieldId(-1),
    m_bondLevelMicroState56FieldId(-1),
    m_bondLevelMicroState57FieldId(-1),
    m_bondLevelMicroState58FieldId(-1),
    m_bondLevelMicroState59FieldId(-1),
    m_bondLevelMicroState60FieldId(-1),
    m_bondLevelMicroState61FieldId(-1),
    m_bondLevelMicroState62FieldId(-1),
    m_bondLevelMicroState63FieldId(-1),
    m_bondLevelMicroState64FieldId(-1),
    m_bondLevelMicroState65FieldId(-1),
    m_bondLevelMicroState66FieldId(-1),
    m_bondLevelMicroState67FieldId(-1),
    m_bondLevelMicroState68FieldId(-1),
    m_bondLevelMicroState69FieldId(-1),
    m_bondLevelMicroState70FieldId(-1),
    m_bondLevelMicroState71FieldId(-1),
    m_bondLevelMicroState72FieldId(-1),
    m_bondLevelMicroState73FieldId(-1),
    m_bondLevelMicroState74FieldId(-1),
    m_bondLevelMicroState75FieldId(-1),
    m_bondLevelMicroState76FieldId(-1),
    m_bondLevelMicroState77FieldId(-1),
    m_bondLevelMicroState78FieldId(-1),
    m_bondLevelMicroState79FieldId(-1),
    m_bondLevelMicroState80FieldId(-1),
    m_bondLevelMicroState81FieldId(-1),
    m_bondLevelMicroState82FieldId(-1),
    m_bondLevelMicroState83FieldId(-1),
    m_bondLevelMicroState84FieldId(-1),
    m_bondLevelMicroState85FieldId(-1),
    m_bondLevelMicroState86FieldId(-1),
    m_bondLevelMicroState87FieldId(-1),
    m_bondLevelMicroState88FieldId(-1),
    m_bondLevelMicroState89FieldId(-1),
    m_bondLevelMicroState90FieldId(-1),
    m_bondLevelMicroState91FieldId(-1),
    m_bondLevelMicroState92FieldId(-1),
    m_bondLevelMicroState93FieldId(-1),
    m_bondLevelMicroState94FieldId(-1),
    m_bondLevelMicroState95FieldId(-1),
    m_bondLevelMicroState96FieldId(-1),
    m_bondLevelMicroState97FieldId(-1),
    m_bondLevelMicroState98FieldId(-1),
    m_bondLevelMicroState99FieldId(-1),
    m_bondLevelMicroState100FieldId(-1),
    m_bondLevelMicroState101FieldId(-1),
    m_bondLevelMicroState102FieldId(-1),
    m_bondLevelMicroState103FieldId(-1),
    m_bondLevelMicroState104FieldId(-1),
    m_bondLevelMicroState105FieldId(-1),
    m_bondLevelMicroState106FieldId(-1),
    m_bondLevelMicroState107FieldId(-1),
    m_bondLevelMicroState108FieldId(-1),
    m_bondLevelMicroState109FieldId(-1),
    m_bondLevelMicroState110FieldId(-1),
    m_bondLevelMicroState111FieldId(-1),
    m_bondLevelMicroState112FieldId(-1),
    m_bondLevelMicroState113FieldId(-1),
    m_bondLevelMicroState114FieldId(-1),
    m_bondLevelMicroState115FieldId(-1),
    m_bondLevelMicroState116FieldId(-1),
    m_bondLevelMicroState117FieldId(-1),
    m_bondLevelMicroState118FieldId(-1),
    m_bondLevelMicroState119FieldId(-1),
    m_bondLevelMicroState120FieldId(-1),
    m_bondLevelMicroState121FieldId(-1),
    m_bondLevelMicroState122FieldId(-1),
    m_bondLevelMicroState123FieldId(-1),
    m_bondLevelMicroState124FieldId(-1),
    m_bondLevelMicroState125FieldId(-1),
    m_bondLevelMicroState126FieldId(-1),
    m_bondLevelMicroState127FieldId(-1),
    m_bondLevelMicroState128FieldId(-1),
    m_bondLevelMicroState129FieldId(-1),
    m_bondLevelMicroState130FieldId(-1),
    m_bondLevelMicroState131FieldId(-1),
    m_bondLevelMicroState132FieldId(-1),
    m_bondLevelMicroState133FieldId(-1),
    m_bondLevelMicroState134FieldId(-1),
    m_bondLevelMicroState135FieldId(-1),
    m_bondLevelMicroState136FieldId(-1),
    m_bondLevelMicroState137FieldId(-1),
    m_bondLevelMicroState138FieldId(-1),
    m_bondLevelMicroState139FieldId(-1),
    m_bondLevelMicroState140FieldId(-1),
    m_bondLevelMicroState141FieldId(-1),
    m_bondLevelMicroState142FieldId(-1),
    m_bondLevelMicroState143FieldId(-1),
    m_bondLevelMicroState144FieldId(-1),
    m_bondLevelMicroState145FieldId(-1),
    m_bondLevelMicroState146FieldId(-1),
    m_bondLevelMicroState147FieldId(-1),
    m_bondLevelMicroState148FieldId(-1),
    m_bondLevelMicroState149FieldId(-1),
    m_bondLevelMicroState150FieldId(-1),
    m_bondLevelMicroState151FieldId(-1),
    m_bondLevelMicroState152FieldId(-1),
    m_bondLevelMicroState153FieldId(-1),
    m_bondLevelMicroState154FieldId(-1),
    m_bondLevelMicroState155FieldId(-1),
    m_bondLevelMicroState156FieldId(-1),
    m_bondLevelMicroState157FieldId(-1),
    m_bondLevelMicroState158FieldId(-1),
    m_bondLevelMicroState159FieldId(-1),
    m_bondLevelMicroState160FieldId(-1),
    m_bondLevelMicroState161FieldId(-1),
    m_bondLevelMicroState162FieldId(-1),
    m_bondLevelMicroState163FieldId(-1),
    m_bondLevelMicroState164FieldId(-1),
    m_bondLevelMicroState165FieldId(-1),
    m_bondLevelMicroState166FieldId(-1),
    m_bondLevelMicroState167FieldId(-1),
    m_bondLevelMicroState168FieldId(-1),
    m_bondLevelMicroState169FieldId(-1),
    m_bondLevelMicroState170FieldId(-1),
    m_bondLevelMicroState171FieldId(-1),
    m_bondLevelMicroState172FieldId(-1),
    m_bondLevelMicroState173FieldId(-1),
    m_bondLevelMicroState174FieldId(-1),
    m_bondLevelMicroState175FieldId(-1),
    m_bondLevelMicroState176FieldId(-1),
    m_bondLevelMicroState177FieldId(-1),
    m_bondLevelMicroState178FieldId(-1),
    m_bondLevelMicroState179FieldId(-1),
    m_bondLevelMicroState180FieldId(-1),
    m_bondLevelMicroState181FieldId(-1),
    m_bondLevelMicroState182FieldId(-1),
    m_bondLevelMicroState183FieldId(-1),
    m_bondLevelMicroState184FieldId(-1),
    m_bondLevelMicroState185FieldId(-1),
    m_bondLevelMicroState186FieldId(-1),
    m_bondLevelMicroState187FieldId(-1),
    m_bondLevelMicroState188FieldId(-1),
    m_bondLevelMicroState189FieldId(-1),
    m_bondLevelInternalEnergyFieldId(-1),
    m_bondLevelInelasticEnergyFieldId(-1)
{
  m_youngsModulus = 9.0 * m_bulkModulus * m_shearModulus / (3.0 * m_bulkModulus + m_shearModulus);
  m_poissonsRatio = m_youngsModulus / (2.0 * m_shearModulus) - 1.0;
  if(m_poissonsRatio > 0.25){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  M7 Microplane model can only reproduce Poisson's ratio <= 0.25>.\n");
  }

  m_k1 = params.get<double>("K1");
  m_k2 = params.get<double>("K2");
  m_k3 = params.get<double>("K3");
  m_k4 = params.get<double>("K4");

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
  m_microState50FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_50");
  m_microState51FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_51");
  m_microState52FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_52");
  m_microState53FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_53");
  m_microState54FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_54");
  m_microState55FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_55");
  m_microState56FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_56");
  m_microState57FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_57");
  m_microState58FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_58");
  m_microState59FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_59");
  m_microState60FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_60");
  m_microState61FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_61");
  m_microState62FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_62");
  m_microState63FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_63");
  m_microState64FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_64");
  m_microState65FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_65");
  m_microState66FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_66");
  m_microState67FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_67");
  m_microState68FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_68");
  m_microState69FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_69");
  m_microState70FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_70");
  m_microState71FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_71");
  m_microState72FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_72");
  m_microState73FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_73");
  m_microState74FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_74");
  m_microState75FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_75");
  m_microState76FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_76");
  m_microState77FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_77");
  m_microState78FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_78");
  m_microState79FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_79");
  m_microState80FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_80");
  m_microState81FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_81");
  m_microState82FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_82");
  m_microState83FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_83");
  m_microState84FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_84");
  m_microState85FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_85");
  m_microState86FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_86");
  m_microState87FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_87");
  m_microState88FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_88");
  m_microState89FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_89");
  m_microState90FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_90");
  m_microState91FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_91");
  m_microState92FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_92");
  m_microState93FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_93");
  m_microState94FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_94");
  m_microState95FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_95");
  m_microState96FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_96");
  m_microState97FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_97");
  m_microState98FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_98");
  m_microState99FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_99");
  m_microState100FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_100");
  m_microState101FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_101");
  m_microState102FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_102");
  m_microState103FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_103");
  m_microState104FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_104");
  m_microState105FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_105");
  m_microState106FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_106");
  m_microState107FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_107");
  m_microState108FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_108");
  m_microState109FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_109");
  m_microState110FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_110");
  m_microState111FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_111");
  m_microState112FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_112");
  m_microState113FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_113");
  m_microState114FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_114");
  m_microState115FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_115");
  m_microState116FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_116");
  m_microState117FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_117");
  m_microState118FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_118");
  m_microState119FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_119");
  m_microState120FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_120");
  m_microState121FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_121");
  m_microState122FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_122");
  m_microState123FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_123");
  m_microState124FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_124");
  m_microState125FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_125");
  m_microState126FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_126");
  m_microState127FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_127");
  m_microState128FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_128");
  m_microState129FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_129");
  m_microState130FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_130");
  m_microState131FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_131");
  m_microState132FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_132");
  m_microState133FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_133");
  m_microState134FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_134");
  m_microState135FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_135");
  m_microState136FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_136");
  m_microState137FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_137");
  m_microState138FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_138");
  m_microState139FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_139");
  m_microState140FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_140");
  m_microState141FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_141");
  m_microState142FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_142");
  m_microState143FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_143");
  m_microState144FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_144");
  m_microState145FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_145");
  m_microState146FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_146");
  m_microState147FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_147");
  m_microState148FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_148");
  m_microState149FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_149");
  m_microState150FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_150");
  m_microState151FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_151");
  m_microState152FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_152");
  m_microState153FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_153");
  m_microState154FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_154");
  m_microState155FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_155");
  m_microState156FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_156");
  m_microState157FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_157");
  m_microState158FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_158");
  m_microState159FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_159");
  m_microState160FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_160");
  m_microState161FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_161");
  m_microState162FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_162");
  m_microState163FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_163");
  m_microState164FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_164");
  m_microState165FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_165");
  m_microState166FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_166");
  m_microState167FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_167");
  m_microState168FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_168");
  m_microState169FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_169");
  m_microState170FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_170");
  m_microState171FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_171");
  m_microState172FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_172");
  m_microState173FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_173");
  m_microState174FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_174");
  m_microState175FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_175");
  m_microState176FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_176");
  m_microState177FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_177");
  m_microState178FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_178");
  m_microState179FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_179");
  m_microState180FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_180");
  m_microState181FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_181");
  m_microState182FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_182");
  m_microState183FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_183");
  m_microState184FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_184");
  m_microState185FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_185");
  m_microState186FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_186");
  m_microState187FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_187");
  m_microState188FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_188");
  m_microState189FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro_State_189");
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
  m_fieldIds.push_back(m_microState50FieldId);
  m_fieldIds.push_back(m_microState51FieldId);
  m_fieldIds.push_back(m_microState52FieldId);
  m_fieldIds.push_back(m_microState53FieldId);
  m_fieldIds.push_back(m_microState54FieldId);
  m_fieldIds.push_back(m_microState55FieldId);
  m_fieldIds.push_back(m_microState56FieldId);
  m_fieldIds.push_back(m_microState57FieldId);
  m_fieldIds.push_back(m_microState58FieldId);
  m_fieldIds.push_back(m_microState59FieldId);
  m_fieldIds.push_back(m_microState60FieldId);
  m_fieldIds.push_back(m_microState61FieldId);
  m_fieldIds.push_back(m_microState62FieldId);
  m_fieldIds.push_back(m_microState63FieldId);
  m_fieldIds.push_back(m_microState64FieldId);
  m_fieldIds.push_back(m_microState65FieldId);
  m_fieldIds.push_back(m_microState66FieldId);
  m_fieldIds.push_back(m_microState67FieldId);
  m_fieldIds.push_back(m_microState68FieldId);
  m_fieldIds.push_back(m_microState69FieldId);
  m_fieldIds.push_back(m_microState70FieldId);
  m_fieldIds.push_back(m_microState71FieldId);
  m_fieldIds.push_back(m_microState72FieldId);
  m_fieldIds.push_back(m_microState73FieldId);
  m_fieldIds.push_back(m_microState74FieldId);
  m_fieldIds.push_back(m_microState75FieldId);
  m_fieldIds.push_back(m_microState76FieldId);
  m_fieldIds.push_back(m_microState77FieldId);
  m_fieldIds.push_back(m_microState78FieldId);
  m_fieldIds.push_back(m_microState79FieldId);
  m_fieldIds.push_back(m_microState80FieldId);
  m_fieldIds.push_back(m_microState81FieldId);
  m_fieldIds.push_back(m_microState82FieldId);
  m_fieldIds.push_back(m_microState83FieldId);
  m_fieldIds.push_back(m_microState84FieldId);
  m_fieldIds.push_back(m_microState85FieldId);
  m_fieldIds.push_back(m_microState86FieldId);
  m_fieldIds.push_back(m_microState87FieldId);
  m_fieldIds.push_back(m_microState88FieldId);
  m_fieldIds.push_back(m_microState89FieldId);
  m_fieldIds.push_back(m_microState90FieldId);
  m_fieldIds.push_back(m_microState91FieldId);
  m_fieldIds.push_back(m_microState92FieldId);
  m_fieldIds.push_back(m_microState93FieldId);
  m_fieldIds.push_back(m_microState94FieldId);
  m_fieldIds.push_back(m_microState95FieldId);
  m_fieldIds.push_back(m_microState96FieldId);
  m_fieldIds.push_back(m_microState97FieldId);
  m_fieldIds.push_back(m_microState98FieldId);
  m_fieldIds.push_back(m_microState99FieldId);
  m_fieldIds.push_back(m_microState100FieldId);
  m_fieldIds.push_back(m_microState101FieldId);
  m_fieldIds.push_back(m_microState102FieldId);
  m_fieldIds.push_back(m_microState103FieldId);
  m_fieldIds.push_back(m_microState104FieldId);
  m_fieldIds.push_back(m_microState105FieldId);
  m_fieldIds.push_back(m_microState106FieldId);
  m_fieldIds.push_back(m_microState107FieldId);
  m_fieldIds.push_back(m_microState108FieldId);
  m_fieldIds.push_back(m_microState109FieldId);
  m_fieldIds.push_back(m_microState110FieldId);
  m_fieldIds.push_back(m_microState111FieldId);
  m_fieldIds.push_back(m_microState112FieldId);
  m_fieldIds.push_back(m_microState113FieldId);
  m_fieldIds.push_back(m_microState114FieldId);
  m_fieldIds.push_back(m_microState115FieldId);
  m_fieldIds.push_back(m_microState116FieldId);
  m_fieldIds.push_back(m_microState117FieldId);
  m_fieldIds.push_back(m_microState118FieldId);
  m_fieldIds.push_back(m_microState119FieldId);
  m_fieldIds.push_back(m_microState120FieldId);
  m_fieldIds.push_back(m_microState121FieldId);
  m_fieldIds.push_back(m_microState122FieldId);
  m_fieldIds.push_back(m_microState123FieldId);
  m_fieldIds.push_back(m_microState124FieldId);
  m_fieldIds.push_back(m_microState125FieldId);
  m_fieldIds.push_back(m_microState126FieldId);
  m_fieldIds.push_back(m_microState127FieldId);
  m_fieldIds.push_back(m_microState128FieldId);
  m_fieldIds.push_back(m_microState129FieldId);
  m_fieldIds.push_back(m_microState130FieldId);
  m_fieldIds.push_back(m_microState131FieldId);
  m_fieldIds.push_back(m_microState132FieldId);
  m_fieldIds.push_back(m_microState133FieldId);
  m_fieldIds.push_back(m_microState134FieldId);
  m_fieldIds.push_back(m_microState135FieldId);
  m_fieldIds.push_back(m_microState136FieldId);
  m_fieldIds.push_back(m_microState137FieldId);
  m_fieldIds.push_back(m_microState138FieldId);
  m_fieldIds.push_back(m_microState139FieldId);
  m_fieldIds.push_back(m_microState140FieldId);
  m_fieldIds.push_back(m_microState141FieldId);
  m_fieldIds.push_back(m_microState142FieldId);
  m_fieldIds.push_back(m_microState143FieldId);
  m_fieldIds.push_back(m_microState144FieldId);
  m_fieldIds.push_back(m_microState145FieldId);
  m_fieldIds.push_back(m_microState146FieldId);
  m_fieldIds.push_back(m_microState147FieldId);
  m_fieldIds.push_back(m_microState148FieldId);
  m_fieldIds.push_back(m_microState149FieldId);
  m_fieldIds.push_back(m_microState150FieldId);
  m_fieldIds.push_back(m_microState151FieldId);
  m_fieldIds.push_back(m_microState152FieldId);
  m_fieldIds.push_back(m_microState153FieldId);
  m_fieldIds.push_back(m_microState154FieldId);
  m_fieldIds.push_back(m_microState155FieldId);
  m_fieldIds.push_back(m_microState156FieldId);
  m_fieldIds.push_back(m_microState157FieldId);
  m_fieldIds.push_back(m_microState158FieldId);
  m_fieldIds.push_back(m_microState159FieldId);
  m_fieldIds.push_back(m_microState160FieldId);
  m_fieldIds.push_back(m_microState161FieldId);
  m_fieldIds.push_back(m_microState162FieldId);
  m_fieldIds.push_back(m_microState163FieldId);
  m_fieldIds.push_back(m_microState164FieldId);
  m_fieldIds.push_back(m_microState165FieldId);
  m_fieldIds.push_back(m_microState166FieldId);
  m_fieldIds.push_back(m_microState167FieldId);
  m_fieldIds.push_back(m_microState168FieldId);
  m_fieldIds.push_back(m_microState169FieldId);
  m_fieldIds.push_back(m_microState170FieldId);
  m_fieldIds.push_back(m_microState171FieldId);
  m_fieldIds.push_back(m_microState172FieldId);
  m_fieldIds.push_back(m_microState173FieldId);
  m_fieldIds.push_back(m_microState174FieldId);
  m_fieldIds.push_back(m_microState175FieldId);
  m_fieldIds.push_back(m_microState176FieldId);
  m_fieldIds.push_back(m_microState177FieldId);
  m_fieldIds.push_back(m_microState178FieldId);
  m_fieldIds.push_back(m_microState179FieldId);
  m_fieldIds.push_back(m_microState180FieldId);
  m_fieldIds.push_back(m_microState181FieldId);
  m_fieldIds.push_back(m_microState182FieldId);
  m_fieldIds.push_back(m_microState183FieldId);
  m_fieldIds.push_back(m_microState184FieldId);
  m_fieldIds.push_back(m_microState185FieldId);
  m_fieldIds.push_back(m_microState186FieldId);
  m_fieldIds.push_back(m_microState187FieldId);
  m_fieldIds.push_back(m_microState188FieldId);
  m_fieldIds.push_back(m_microState189FieldId);
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
  m_bondLevelMicroState50FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_50");
  m_bondLevelMicroState51FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_51");
  m_bondLevelMicroState52FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_52");
  m_bondLevelMicroState53FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_53");
  m_bondLevelMicroState54FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_54");
  m_bondLevelMicroState55FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_55");
  m_bondLevelMicroState56FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_56");
  m_bondLevelMicroState57FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_57");
  m_bondLevelMicroState58FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_58");
  m_bondLevelMicroState59FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_59");
  m_bondLevelMicroState60FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_60");
  m_bondLevelMicroState61FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_61");
  m_bondLevelMicroState62FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_62");
  m_bondLevelMicroState63FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_63");
  m_bondLevelMicroState64FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_64");
  m_bondLevelMicroState65FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_65");
  m_bondLevelMicroState66FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_66");
  m_bondLevelMicroState67FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_67");
  m_bondLevelMicroState68FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_68");
  m_bondLevelMicroState69FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_69");
  m_bondLevelMicroState70FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_70");
  m_bondLevelMicroState71FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_71");
  m_bondLevelMicroState72FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_72");
  m_bondLevelMicroState73FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_73");
  m_bondLevelMicroState74FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_74");
  m_bondLevelMicroState75FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_75");
  m_bondLevelMicroState76FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_76");
  m_bondLevelMicroState77FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_77");
  m_bondLevelMicroState78FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_78");
  m_bondLevelMicroState79FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_79");
  m_bondLevelMicroState80FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_80");
  m_bondLevelMicroState81FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_81");
  m_bondLevelMicroState82FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_82");
  m_bondLevelMicroState83FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_83");
  m_bondLevelMicroState84FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_84");
  m_bondLevelMicroState85FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_85");
  m_bondLevelMicroState86FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_86");
  m_bondLevelMicroState87FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_87");
  m_bondLevelMicroState88FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_88");
  m_bondLevelMicroState89FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_89");
  m_bondLevelMicroState90FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_90");
  m_bondLevelMicroState91FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_91");
  m_bondLevelMicroState92FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_92");
  m_bondLevelMicroState93FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_93");
  m_bondLevelMicroState94FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_94");
  m_bondLevelMicroState95FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_95");
  m_bondLevelMicroState96FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_96");
  m_bondLevelMicroState97FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_97");
  m_bondLevelMicroState98FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_98");
  m_bondLevelMicroState99FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_99");
  m_bondLevelMicroState100FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_100");
  m_bondLevelMicroState101FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_101");
  m_bondLevelMicroState102FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_102");
  m_bondLevelMicroState103FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_103");
  m_bondLevelMicroState104FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_104");
  m_bondLevelMicroState105FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_105");
  m_bondLevelMicroState106FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_106");
  m_bondLevelMicroState107FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_107");
  m_bondLevelMicroState108FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_108");
  m_bondLevelMicroState109FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_109");
  m_bondLevelMicroState110FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_110");
  m_bondLevelMicroState111FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_111");
  m_bondLevelMicroState112FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_112");
  m_bondLevelMicroState113FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_113");
  m_bondLevelMicroState114FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_114");
  m_bondLevelMicroState115FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_115");
  m_bondLevelMicroState116FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_116");
  m_bondLevelMicroState117FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_117");
  m_bondLevelMicroState118FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_118");
  m_bondLevelMicroState119FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_119");
  m_bondLevelMicroState120FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_120");
  m_bondLevelMicroState121FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_121");
  m_bondLevelMicroState122FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_122");
  m_bondLevelMicroState123FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_123");
  m_bondLevelMicroState124FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_124");
  m_bondLevelMicroState125FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_125");
  m_bondLevelMicroState126FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_126");
  m_bondLevelMicroState127FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_127");
  m_bondLevelMicroState128FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_128");
  m_bondLevelMicroState129FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_129");
  m_bondLevelMicroState130FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_130");
  m_bondLevelMicroState131FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_131");
  m_bondLevelMicroState132FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_132");
  m_bondLevelMicroState133FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_133");
  m_bondLevelMicroState134FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_134");
  m_bondLevelMicroState135FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_135");
  m_bondLevelMicroState136FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_136");
  m_bondLevelMicroState137FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_137");
  m_bondLevelMicroState138FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_138");
  m_bondLevelMicroState139FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_139");
  m_bondLevelMicroState140FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_140");
  m_bondLevelMicroState141FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_141");
  m_bondLevelMicroState142FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_142");
  m_bondLevelMicroState143FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_143");
  m_bondLevelMicroState144FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_144");
  m_bondLevelMicroState145FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_145");
  m_bondLevelMicroState146FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_146");
  m_bondLevelMicroState147FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_147");
  m_bondLevelMicroState148FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_148");
  m_bondLevelMicroState149FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_149");
  m_bondLevelMicroState150FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_150");
  m_bondLevelMicroState151FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_151");
  m_bondLevelMicroState152FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_152");
  m_bondLevelMicroState153FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_153");
  m_bondLevelMicroState154FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_154");
  m_bondLevelMicroState155FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_155");
  m_bondLevelMicroState156FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_156");
  m_bondLevelMicroState157FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_157");
  m_bondLevelMicroState158FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_158");
  m_bondLevelMicroState159FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_159");
  m_bondLevelMicroState160FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_160");
  m_bondLevelMicroState161FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_161");
  m_bondLevelMicroState162FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_162");
  m_bondLevelMicroState163FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_163");
  m_bondLevelMicroState164FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_164");
  m_bondLevelMicroState165FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_165");
  m_bondLevelMicroState166FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_166");
  m_bondLevelMicroState167FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_167");
  m_bondLevelMicroState168FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_168");
  m_bondLevelMicroState169FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_169");
  m_bondLevelMicroState170FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_170");
  m_bondLevelMicroState171FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_171");
  m_bondLevelMicroState172FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_172");
  m_bondLevelMicroState173FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_173");
  m_bondLevelMicroState174FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_174");
  m_bondLevelMicroState175FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_175");
  m_bondLevelMicroState176FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_176");
  m_bondLevelMicroState177FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_177");
  m_bondLevelMicroState178FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_178");
  m_bondLevelMicroState179FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_179");
  m_bondLevelMicroState180FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_180");
  m_bondLevelMicroState181FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_181");
  m_bondLevelMicroState182FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_182");
  m_bondLevelMicroState183FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_183");
  m_bondLevelMicroState184FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_184");
  m_bondLevelMicroState185FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_185");
  m_bondLevelMicroState186FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_186");
  m_bondLevelMicroState187FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_187");
  m_bondLevelMicroState188FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_188");
  m_bondLevelMicroState189FieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Micro_State_189");
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
  m_fieldIds.push_back(m_bondLevelMicroState50FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState51FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState52FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState53FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState54FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState55FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState56FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState57FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState58FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState59FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState60FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState61FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState62FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState63FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState64FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState65FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState66FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState67FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState68FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState69FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState70FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState71FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState72FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState73FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState74FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState75FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState76FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState77FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState78FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState79FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState80FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState81FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState82FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState83FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState84FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState85FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState86FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState87FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState88FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState89FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState90FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState91FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState92FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState93FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState94FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState95FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState96FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState97FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState98FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState99FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState100FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState101FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState102FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState103FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState104FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState105FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState106FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState107FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState108FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState109FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState110FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState111FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState112FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState113FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState114FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState115FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState116FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState117FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState118FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState119FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState120FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState121FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState122FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState123FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState124FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState125FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState126FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState127FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState128FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState129FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState130FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState131FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState132FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState133FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState134FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState135FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState136FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState137FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState138FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState139FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState140FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState141FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState142FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState143FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState144FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState145FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState146FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState147FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState148FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState149FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState150FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState151FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState152FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState153FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState154FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState155FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState156FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState157FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState158FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState159FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState160FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState161FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState162FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState163FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState164FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState165FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState166FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState167FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState168FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState169FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState170FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState171FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState172FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState173FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState174FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState175FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState176FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState177FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState178FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState179FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState180FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState181FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState182FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState183FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState184FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState185FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState186FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState187FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState188FieldId);
  m_fieldIds.push_back(m_bondLevelMicroState189FieldId);
  m_fieldIds.push_back(m_bondLevelInternalEnergyFieldId);
  m_fieldIds.push_back(m_bondLevelInelasticEnergyFieldId);
}

PeridigmNS::MicroplaneBondAssociatedCorrespondenceMaterial::~MicroplaneBondAssociatedCorrespondenceMaterial()
{
}

void
PeridigmNS::MicroplaneBondAssociatedCorrespondenceMaterial::initialize(const double dt,
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

  CORRESPONDENCE::initializeMicroplaneM7Model(m_youngsModulus,
                                              m_poissonsRatio,
                                              m_k1,
                                              m_k2,
                                              m_k3,
                                              m_k4);

  dataManager.getData(m_microState1FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState2FieldId, PeridigmField::STEP_N)->PutScalar(1.0);
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
  dataManager.getData(m_microState50FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState51FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState52FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState53FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState54FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState55FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState56FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState57FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState58FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState59FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState60FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState61FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState62FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState63FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState64FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState65FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState66FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState67FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState68FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState69FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState70FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState71FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState72FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState73FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState74FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState75FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState76FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState77FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState78FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState79FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState80FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState81FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState82FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState83FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState84FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState85FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState86FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState87FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState88FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState89FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState90FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState91FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState92FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState93FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState94FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState95FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState96FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState97FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState98FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState99FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState100FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState101FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState102FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState103FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState104FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState105FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState106FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState107FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState108FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState109FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState110FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState111FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState112FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState113FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState114FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState115FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState116FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState117FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState118FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState119FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState120FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState121FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState122FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState123FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState124FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState125FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState126FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState127FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState128FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState129FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState130FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState131FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState132FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState133FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState134FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState135FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState136FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState137FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState138FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState139FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState140FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState141FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState142FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState143FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState144FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState145FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState146FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState147FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState148FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState149FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState150FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState151FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState152FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState153FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState154FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState155FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState156FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState157FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState158FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState159FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState160FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState161FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState162FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState163FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState164FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState165FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState166FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState167FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState168FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState169FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState170FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState171FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState172FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState173FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState174FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState175FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState176FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState177FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState178FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState179FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState180FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState181FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState182FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState183FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState184FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState185FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState186FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState187FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState188FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_microState189FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);

  dataManager.getData(m_microState1FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState2FieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);
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
  dataManager.getData(m_microState50FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState51FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState52FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState53FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState54FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState55FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState56FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState57FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState58FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState59FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState60FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState61FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState62FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState63FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState64FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState65FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState66FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState67FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState68FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState69FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState70FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState71FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState72FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState73FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState74FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState75FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState76FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState77FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState78FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState79FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState80FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState81FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState82FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState83FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState84FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState85FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState86FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState87FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState88FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState89FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState90FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState91FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState92FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState93FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState94FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState95FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState96FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState97FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState98FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState99FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState100FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState101FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState102FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState103FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState104FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState105FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState106FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState107FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState108FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState109FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState110FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState111FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState112FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState113FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState114FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState115FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState116FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState117FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState118FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState119FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState120FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState121FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState122FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState123FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState124FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState125FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState126FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState127FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState128FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState129FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState130FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState131FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState132FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState133FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState134FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState135FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState136FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState137FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState138FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState139FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState140FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState141FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState142FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState143FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState144FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState145FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState146FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState147FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState148FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState149FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState150FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState151FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState152FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState153FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState154FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState155FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState156FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState157FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState158FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState159FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState160FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState161FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState162FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState163FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState164FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState165FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState166FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState167FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState168FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState169FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState170FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState171FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState172FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState173FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState174FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState175FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState176FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState177FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState178FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState179FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState180FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState181FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState182FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState183FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState184FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState185FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState186FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState187FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState188FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_microState189FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(m_bondLevelMicroState1FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState2FieldId, PeridigmField::STEP_N)->PutScalar(1.0);
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
  dataManager.getData(m_bondLevelMicroState50FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState51FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState52FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState53FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState54FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState55FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState56FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState57FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState58FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState59FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState60FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState61FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState62FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState63FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState64FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState65FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState66FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState67FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState68FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState69FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState70FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState71FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState72FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState73FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState74FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState75FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState76FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState77FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState78FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState79FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState80FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState81FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState82FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState83FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState84FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState85FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState86FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState87FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState88FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState89FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState90FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState91FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState92FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState93FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState94FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState95FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState96FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState97FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState98FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState99FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState100FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState101FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState102FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState103FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState104FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState105FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState106FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState107FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState108FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState109FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState110FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState111FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState112FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState113FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState114FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState115FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState116FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState117FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState118FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState119FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState120FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState121FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState122FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState123FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState124FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState125FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState126FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState127FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState128FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState129FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState130FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState131FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState132FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState133FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState134FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState135FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState136FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState137FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState138FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState139FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState140FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState141FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState142FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState143FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState144FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState145FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState146FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState147FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState148FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState149FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState150FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState151FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState152FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState153FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState154FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState155FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState156FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState157FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState158FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState159FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState160FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState161FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState162FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState163FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState164FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState165FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState166FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState167FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState168FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState169FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState170FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState171FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState172FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState173FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState174FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState175FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState176FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState177FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState178FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState179FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState180FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState181FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState182FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState183FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState184FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState185FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState186FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState187FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState188FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState189FieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);

  dataManager.getData(m_bondLevelMicroState1FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState2FieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);
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
  dataManager.getData(m_bondLevelMicroState50FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState51FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState52FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState53FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState54FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState55FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState56FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState57FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState58FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState59FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState60FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState61FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState62FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState63FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState64FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState65FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState66FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState67FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState68FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState69FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState70FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState71FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState72FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState73FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState74FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState75FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState76FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState77FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState78FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState79FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState80FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState81FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState82FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState83FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState84FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState85FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState86FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState87FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState88FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState89FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState90FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState91FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState92FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState93FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState94FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState95FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState96FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState97FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState98FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState99FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState100FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState101FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState102FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState103FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState104FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState105FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState106FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState107FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState108FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState109FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState110FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState111FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState112FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState113FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState114FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState115FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState116FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState117FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState118FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState119FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState120FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState121FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState122FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState123FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState124FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState125FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState126FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState127FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState128FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState129FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState130FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState131FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState132FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState133FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState134FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState135FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState136FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState137FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState138FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState139FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState140FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState141FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState142FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState143FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState144FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState145FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState146FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState147FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState148FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState149FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState150FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState151FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState152FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState153FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState154FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState155FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState156FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState157FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState158FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState159FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState160FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState161FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState162FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState163FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState164FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState165FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState166FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState167FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState168FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState169FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState170FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState171FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState172FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState173FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState174FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState175FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState176FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState177FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState178FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState179FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState180FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState181FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState182FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState183FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState184FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState185FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState186FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState187FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState188FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelMicroState189FieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}

void
PeridigmNS::MicroplaneBondAssociatedCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                                                const int numOwnedPoints,
                                                                                const int* neighborhoodList,
                                                                                PeridigmNS::DataManager& dataManager) const
{

}

void
PeridigmNS::MicroplaneBondAssociatedCorrespondenceMaterial::computePK2Stress(const double dt,
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
  double *state50N;
  double *state51N;
  double *state52N;
  double *state53N;
  double *state54N;
  double *state55N;
  double *state56N;
  double *state57N;
  double *state58N;
  double *state59N;
  double *state60N;
  double *state61N;
  double *state62N;
  double *state63N;
  double *state64N;
  double *state65N;
  double *state66N;
  double *state67N;
  double *state68N;
  double *state69N;
  double *state70N;
  double *state71N;
  double *state72N;
  double *state73N;
  double *state74N;
  double *state75N;
  double *state76N;
  double *state77N;
  double *state78N;
  double *state79N;
  double *state80N;
  double *state81N;
  double *state82N;
  double *state83N;
  double *state84N;
  double *state85N;
  double *state86N;
  double *state87N;
  double *state88N;
  double *state89N;
  double *state90N;
  double *state91N;
  double *state92N;
  double *state93N;
  double *state94N;
  double *state95N;
  double *state96N;
  double *state97N;
  double *state98N;
  double *state99N;
  double *state100N;
  double *state101N;
  double *state102N;
  double *state103N;
  double *state104N;
  double *state105N;
  double *state106N;
  double *state107N;
  double *state108N;
  double *state109N;
  double *state110N;
  double *state111N;
  double *state112N;
  double *state113N;
  double *state114N;
  double *state115N;
  double *state116N;
  double *state117N;
  double *state118N;
  double *state119N;
  double *state120N;
  double *state121N;
  double *state122N;
  double *state123N;
  double *state124N;
  double *state125N;
  double *state126N;
  double *state127N;
  double *state128N;
  double *state129N;
  double *state130N;
  double *state131N;
  double *state132N;
  double *state133N;
  double *state134N;
  double *state135N;
  double *state136N;
  double *state137N;
  double *state138N;
  double *state139N;
  double *state140N;
  double *state141N;
  double *state142N;
  double *state143N;
  double *state144N;
  double *state145N;
  double *state146N;
  double *state147N;
  double *state148N;
  double *state149N;
  double *state150N;
  double *state151N;
  double *state152N;
  double *state153N;
  double *state154N;
  double *state155N;
  double *state156N;
  double *state157N;
  double *state158N;
  double *state159N;
  double *state160N;
  double *state161N;
  double *state162N;
  double *state163N;
  double *state164N;
  double *state165N;
  double *state166N;
  double *state167N;
  double *state168N;
  double *state169N;
  double *state170N;
  double *state171N;
  double *state172N;
  double *state173N;
  double *state174N;
  double *state175N;
  double *state176N;
  double *state177N;
  double *state178N;
  double *state179N;
  double *state180N;
  double *state181N;
  double *state182N;
  double *state183N;
  double *state184N;
  double *state185N;
  double *state186N;
  double *state187N;
  double *state188N;
  double *state189N;
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
  dataManager.getData(m_microState50FieldId, PeridigmField::STEP_N)->ExtractView(&state50N);
  dataManager.getData(m_microState51FieldId, PeridigmField::STEP_N)->ExtractView(&state51N);
  dataManager.getData(m_microState52FieldId, PeridigmField::STEP_N)->ExtractView(&state52N);
  dataManager.getData(m_microState53FieldId, PeridigmField::STEP_N)->ExtractView(&state53N);
  dataManager.getData(m_microState54FieldId, PeridigmField::STEP_N)->ExtractView(&state54N);
  dataManager.getData(m_microState55FieldId, PeridigmField::STEP_N)->ExtractView(&state55N);
  dataManager.getData(m_microState56FieldId, PeridigmField::STEP_N)->ExtractView(&state56N);
  dataManager.getData(m_microState57FieldId, PeridigmField::STEP_N)->ExtractView(&state57N);
  dataManager.getData(m_microState58FieldId, PeridigmField::STEP_N)->ExtractView(&state58N);
  dataManager.getData(m_microState59FieldId, PeridigmField::STEP_N)->ExtractView(&state59N);
  dataManager.getData(m_microState60FieldId, PeridigmField::STEP_N)->ExtractView(&state60N);
  dataManager.getData(m_microState61FieldId, PeridigmField::STEP_N)->ExtractView(&state61N);
  dataManager.getData(m_microState62FieldId, PeridigmField::STEP_N)->ExtractView(&state62N);
  dataManager.getData(m_microState63FieldId, PeridigmField::STEP_N)->ExtractView(&state63N);
  dataManager.getData(m_microState64FieldId, PeridigmField::STEP_N)->ExtractView(&state64N);
  dataManager.getData(m_microState65FieldId, PeridigmField::STEP_N)->ExtractView(&state65N);
  dataManager.getData(m_microState66FieldId, PeridigmField::STEP_N)->ExtractView(&state66N);
  dataManager.getData(m_microState67FieldId, PeridigmField::STEP_N)->ExtractView(&state67N);
  dataManager.getData(m_microState68FieldId, PeridigmField::STEP_N)->ExtractView(&state68N);
  dataManager.getData(m_microState69FieldId, PeridigmField::STEP_N)->ExtractView(&state69N);
  dataManager.getData(m_microState70FieldId, PeridigmField::STEP_N)->ExtractView(&state70N);
  dataManager.getData(m_microState71FieldId, PeridigmField::STEP_N)->ExtractView(&state71N);
  dataManager.getData(m_microState72FieldId, PeridigmField::STEP_N)->ExtractView(&state72N);
  dataManager.getData(m_microState73FieldId, PeridigmField::STEP_N)->ExtractView(&state73N);
  dataManager.getData(m_microState74FieldId, PeridigmField::STEP_N)->ExtractView(&state74N);
  dataManager.getData(m_microState75FieldId, PeridigmField::STEP_N)->ExtractView(&state75N);
  dataManager.getData(m_microState76FieldId, PeridigmField::STEP_N)->ExtractView(&state76N);
  dataManager.getData(m_microState77FieldId, PeridigmField::STEP_N)->ExtractView(&state77N);
  dataManager.getData(m_microState78FieldId, PeridigmField::STEP_N)->ExtractView(&state78N);
  dataManager.getData(m_microState79FieldId, PeridigmField::STEP_N)->ExtractView(&state79N);
  dataManager.getData(m_microState80FieldId, PeridigmField::STEP_N)->ExtractView(&state80N);
  dataManager.getData(m_microState81FieldId, PeridigmField::STEP_N)->ExtractView(&state81N);
  dataManager.getData(m_microState82FieldId, PeridigmField::STEP_N)->ExtractView(&state82N);
  dataManager.getData(m_microState83FieldId, PeridigmField::STEP_N)->ExtractView(&state83N);
  dataManager.getData(m_microState84FieldId, PeridigmField::STEP_N)->ExtractView(&state84N);
  dataManager.getData(m_microState85FieldId, PeridigmField::STEP_N)->ExtractView(&state85N);
  dataManager.getData(m_microState86FieldId, PeridigmField::STEP_N)->ExtractView(&state86N);
  dataManager.getData(m_microState87FieldId, PeridigmField::STEP_N)->ExtractView(&state87N);
  dataManager.getData(m_microState88FieldId, PeridigmField::STEP_N)->ExtractView(&state88N);
  dataManager.getData(m_microState89FieldId, PeridigmField::STEP_N)->ExtractView(&state89N);
  dataManager.getData(m_microState90FieldId, PeridigmField::STEP_N)->ExtractView(&state90N);
  dataManager.getData(m_microState91FieldId, PeridigmField::STEP_N)->ExtractView(&state91N);
  dataManager.getData(m_microState92FieldId, PeridigmField::STEP_N)->ExtractView(&state92N);
  dataManager.getData(m_microState93FieldId, PeridigmField::STEP_N)->ExtractView(&state93N);
  dataManager.getData(m_microState94FieldId, PeridigmField::STEP_N)->ExtractView(&state94N);
  dataManager.getData(m_microState95FieldId, PeridigmField::STEP_N)->ExtractView(&state95N);
  dataManager.getData(m_microState96FieldId, PeridigmField::STEP_N)->ExtractView(&state96N);
  dataManager.getData(m_microState97FieldId, PeridigmField::STEP_N)->ExtractView(&state97N);
  dataManager.getData(m_microState98FieldId, PeridigmField::STEP_N)->ExtractView(&state98N);
  dataManager.getData(m_microState99FieldId, PeridigmField::STEP_N)->ExtractView(&state99N);
  dataManager.getData(m_microState100FieldId, PeridigmField::STEP_N)->ExtractView(&state100N);
  dataManager.getData(m_microState101FieldId, PeridigmField::STEP_N)->ExtractView(&state101N);
  dataManager.getData(m_microState102FieldId, PeridigmField::STEP_N)->ExtractView(&state102N);
  dataManager.getData(m_microState103FieldId, PeridigmField::STEP_N)->ExtractView(&state103N);
  dataManager.getData(m_microState104FieldId, PeridigmField::STEP_N)->ExtractView(&state104N);
  dataManager.getData(m_microState105FieldId, PeridigmField::STEP_N)->ExtractView(&state105N);
  dataManager.getData(m_microState106FieldId, PeridigmField::STEP_N)->ExtractView(&state106N);
  dataManager.getData(m_microState107FieldId, PeridigmField::STEP_N)->ExtractView(&state107N);
  dataManager.getData(m_microState108FieldId, PeridigmField::STEP_N)->ExtractView(&state108N);
  dataManager.getData(m_microState109FieldId, PeridigmField::STEP_N)->ExtractView(&state109N);
  dataManager.getData(m_microState110FieldId, PeridigmField::STEP_N)->ExtractView(&state110N);
  dataManager.getData(m_microState111FieldId, PeridigmField::STEP_N)->ExtractView(&state111N);
  dataManager.getData(m_microState112FieldId, PeridigmField::STEP_N)->ExtractView(&state112N);
  dataManager.getData(m_microState113FieldId, PeridigmField::STEP_N)->ExtractView(&state113N);
  dataManager.getData(m_microState114FieldId, PeridigmField::STEP_N)->ExtractView(&state114N);
  dataManager.getData(m_microState115FieldId, PeridigmField::STEP_N)->ExtractView(&state115N);
  dataManager.getData(m_microState116FieldId, PeridigmField::STEP_N)->ExtractView(&state116N);
  dataManager.getData(m_microState117FieldId, PeridigmField::STEP_N)->ExtractView(&state117N);
  dataManager.getData(m_microState118FieldId, PeridigmField::STEP_N)->ExtractView(&state118N);
  dataManager.getData(m_microState119FieldId, PeridigmField::STEP_N)->ExtractView(&state119N);
  dataManager.getData(m_microState120FieldId, PeridigmField::STEP_N)->ExtractView(&state120N);
  dataManager.getData(m_microState121FieldId, PeridigmField::STEP_N)->ExtractView(&state121N);
  dataManager.getData(m_microState122FieldId, PeridigmField::STEP_N)->ExtractView(&state122N);
  dataManager.getData(m_microState123FieldId, PeridigmField::STEP_N)->ExtractView(&state123N);
  dataManager.getData(m_microState124FieldId, PeridigmField::STEP_N)->ExtractView(&state124N);
  dataManager.getData(m_microState125FieldId, PeridigmField::STEP_N)->ExtractView(&state125N);
  dataManager.getData(m_microState126FieldId, PeridigmField::STEP_N)->ExtractView(&state126N);
  dataManager.getData(m_microState127FieldId, PeridigmField::STEP_N)->ExtractView(&state127N);
  dataManager.getData(m_microState128FieldId, PeridigmField::STEP_N)->ExtractView(&state128N);
  dataManager.getData(m_microState129FieldId, PeridigmField::STEP_N)->ExtractView(&state129N);
  dataManager.getData(m_microState130FieldId, PeridigmField::STEP_N)->ExtractView(&state130N);
  dataManager.getData(m_microState131FieldId, PeridigmField::STEP_N)->ExtractView(&state131N);
  dataManager.getData(m_microState132FieldId, PeridigmField::STEP_N)->ExtractView(&state132N);
  dataManager.getData(m_microState133FieldId, PeridigmField::STEP_N)->ExtractView(&state133N);
  dataManager.getData(m_microState134FieldId, PeridigmField::STEP_N)->ExtractView(&state134N);
  dataManager.getData(m_microState135FieldId, PeridigmField::STEP_N)->ExtractView(&state135N);
  dataManager.getData(m_microState136FieldId, PeridigmField::STEP_N)->ExtractView(&state136N);
  dataManager.getData(m_microState137FieldId, PeridigmField::STEP_N)->ExtractView(&state137N);
  dataManager.getData(m_microState138FieldId, PeridigmField::STEP_N)->ExtractView(&state138N);
  dataManager.getData(m_microState139FieldId, PeridigmField::STEP_N)->ExtractView(&state139N);
  dataManager.getData(m_microState140FieldId, PeridigmField::STEP_N)->ExtractView(&state140N);
  dataManager.getData(m_microState141FieldId, PeridigmField::STEP_N)->ExtractView(&state141N);
  dataManager.getData(m_microState142FieldId, PeridigmField::STEP_N)->ExtractView(&state142N);
  dataManager.getData(m_microState143FieldId, PeridigmField::STEP_N)->ExtractView(&state143N);
  dataManager.getData(m_microState144FieldId, PeridigmField::STEP_N)->ExtractView(&state144N);
  dataManager.getData(m_microState145FieldId, PeridigmField::STEP_N)->ExtractView(&state145N);
  dataManager.getData(m_microState146FieldId, PeridigmField::STEP_N)->ExtractView(&state146N);
  dataManager.getData(m_microState147FieldId, PeridigmField::STEP_N)->ExtractView(&state147N);
  dataManager.getData(m_microState148FieldId, PeridigmField::STEP_N)->ExtractView(&state148N);
  dataManager.getData(m_microState149FieldId, PeridigmField::STEP_N)->ExtractView(&state149N);
  dataManager.getData(m_microState150FieldId, PeridigmField::STEP_N)->ExtractView(&state150N);
  dataManager.getData(m_microState151FieldId, PeridigmField::STEP_N)->ExtractView(&state151N);
  dataManager.getData(m_microState152FieldId, PeridigmField::STEP_N)->ExtractView(&state152N);
  dataManager.getData(m_microState153FieldId, PeridigmField::STEP_N)->ExtractView(&state153N);
  dataManager.getData(m_microState154FieldId, PeridigmField::STEP_N)->ExtractView(&state154N);
  dataManager.getData(m_microState155FieldId, PeridigmField::STEP_N)->ExtractView(&state155N);
  dataManager.getData(m_microState156FieldId, PeridigmField::STEP_N)->ExtractView(&state156N);
  dataManager.getData(m_microState157FieldId, PeridigmField::STEP_N)->ExtractView(&state157N);
  dataManager.getData(m_microState158FieldId, PeridigmField::STEP_N)->ExtractView(&state158N);
  dataManager.getData(m_microState159FieldId, PeridigmField::STEP_N)->ExtractView(&state159N);
  dataManager.getData(m_microState160FieldId, PeridigmField::STEP_N)->ExtractView(&state160N);
  dataManager.getData(m_microState161FieldId, PeridigmField::STEP_N)->ExtractView(&state161N);
  dataManager.getData(m_microState162FieldId, PeridigmField::STEP_N)->ExtractView(&state162N);
  dataManager.getData(m_microState163FieldId, PeridigmField::STEP_N)->ExtractView(&state163N);
  dataManager.getData(m_microState164FieldId, PeridigmField::STEP_N)->ExtractView(&state164N);
  dataManager.getData(m_microState165FieldId, PeridigmField::STEP_N)->ExtractView(&state165N);
  dataManager.getData(m_microState166FieldId, PeridigmField::STEP_N)->ExtractView(&state166N);
  dataManager.getData(m_microState167FieldId, PeridigmField::STEP_N)->ExtractView(&state167N);
  dataManager.getData(m_microState168FieldId, PeridigmField::STEP_N)->ExtractView(&state168N);
  dataManager.getData(m_microState169FieldId, PeridigmField::STEP_N)->ExtractView(&state169N);
  dataManager.getData(m_microState170FieldId, PeridigmField::STEP_N)->ExtractView(&state170N);
  dataManager.getData(m_microState171FieldId, PeridigmField::STEP_N)->ExtractView(&state171N);
  dataManager.getData(m_microState172FieldId, PeridigmField::STEP_N)->ExtractView(&state172N);
  dataManager.getData(m_microState173FieldId, PeridigmField::STEP_N)->ExtractView(&state173N);
  dataManager.getData(m_microState174FieldId, PeridigmField::STEP_N)->ExtractView(&state174N);
  dataManager.getData(m_microState175FieldId, PeridigmField::STEP_N)->ExtractView(&state175N);
  dataManager.getData(m_microState176FieldId, PeridigmField::STEP_N)->ExtractView(&state176N);
  dataManager.getData(m_microState177FieldId, PeridigmField::STEP_N)->ExtractView(&state177N);
  dataManager.getData(m_microState178FieldId, PeridigmField::STEP_N)->ExtractView(&state178N);
  dataManager.getData(m_microState179FieldId, PeridigmField::STEP_N)->ExtractView(&state179N);
  dataManager.getData(m_microState180FieldId, PeridigmField::STEP_N)->ExtractView(&state180N);
  dataManager.getData(m_microState181FieldId, PeridigmField::STEP_N)->ExtractView(&state181N);
  dataManager.getData(m_microState182FieldId, PeridigmField::STEP_N)->ExtractView(&state182N);
  dataManager.getData(m_microState183FieldId, PeridigmField::STEP_N)->ExtractView(&state183N);
  dataManager.getData(m_microState184FieldId, PeridigmField::STEP_N)->ExtractView(&state184N);
  dataManager.getData(m_microState185FieldId, PeridigmField::STEP_N)->ExtractView(&state185N);
  dataManager.getData(m_microState186FieldId, PeridigmField::STEP_N)->ExtractView(&state186N);
  dataManager.getData(m_microState187FieldId, PeridigmField::STEP_N)->ExtractView(&state187N);
  dataManager.getData(m_microState188FieldId, PeridigmField::STEP_N)->ExtractView(&state188N);
  dataManager.getData(m_microState189FieldId, PeridigmField::STEP_N)->ExtractView(&state189N);

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
  double *state50NP1;
  double *state51NP1;
  double *state52NP1;
  double *state53NP1;
  double *state54NP1;
  double *state55NP1;
  double *state56NP1;
  double *state57NP1;
  double *state58NP1;
  double *state59NP1;
  double *state60NP1;
  double *state61NP1;
  double *state62NP1;
  double *state63NP1;
  double *state64NP1;
  double *state65NP1;
  double *state66NP1;
  double *state67NP1;
  double *state68NP1;
  double *state69NP1;
  double *state70NP1;
  double *state71NP1;
  double *state72NP1;
  double *state73NP1;
  double *state74NP1;
  double *state75NP1;
  double *state76NP1;
  double *state77NP1;
  double *state78NP1;
  double *state79NP1;
  double *state80NP1;
  double *state81NP1;
  double *state82NP1;
  double *state83NP1;
  double *state84NP1;
  double *state85NP1;
  double *state86NP1;
  double *state87NP1;
  double *state88NP1;
  double *state89NP1;
  double *state90NP1;
  double *state91NP1;
  double *state92NP1;
  double *state93NP1;
  double *state94NP1;
  double *state95NP1;
  double *state96NP1;
  double *state97NP1;
  double *state98NP1;
  double *state99NP1;
  double *state100NP1;
  double *state101NP1;
  double *state102NP1;
  double *state103NP1;
  double *state104NP1;
  double *state105NP1;
  double *state106NP1;
  double *state107NP1;
  double *state108NP1;
  double *state109NP1;
  double *state110NP1;
  double *state111NP1;
  double *state112NP1;
  double *state113NP1;
  double *state114NP1;
  double *state115NP1;
  double *state116NP1;
  double *state117NP1;
  double *state118NP1;
  double *state119NP1;
  double *state120NP1;
  double *state121NP1;
  double *state122NP1;
  double *state123NP1;
  double *state124NP1;
  double *state125NP1;
  double *state126NP1;
  double *state127NP1;
  double *state128NP1;
  double *state129NP1;
  double *state130NP1;
  double *state131NP1;
  double *state132NP1;
  double *state133NP1;
  double *state134NP1;
  double *state135NP1;
  double *state136NP1;
  double *state137NP1;
  double *state138NP1;
  double *state139NP1;
  double *state140NP1;
  double *state141NP1;
  double *state142NP1;
  double *state143NP1;
  double *state144NP1;
  double *state145NP1;
  double *state146NP1;
  double *state147NP1;
  double *state148NP1;
  double *state149NP1;
  double *state150NP1;
  double *state151NP1;
  double *state152NP1;
  double *state153NP1;
  double *state154NP1;
  double *state155NP1;
  double *state156NP1;
  double *state157NP1;
  double *state158NP1;
  double *state159NP1;
  double *state160NP1;
  double *state161NP1;
  double *state162NP1;
  double *state163NP1;
  double *state164NP1;
  double *state165NP1;
  double *state166NP1;
  double *state167NP1;
  double *state168NP1;
  double *state169NP1;
  double *state170NP1;
  double *state171NP1;
  double *state172NP1;
  double *state173NP1;
  double *state174NP1;
  double *state175NP1;
  double *state176NP1;
  double *state177NP1;
  double *state178NP1;
  double *state179NP1;
  double *state180NP1;
  double *state181NP1;
  double *state182NP1;
  double *state183NP1;
  double *state184NP1;
  double *state185NP1;
  double *state186NP1;
  double *state187NP1;
  double *state188NP1;
  double *state189NP1;
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
  dataManager.getData(m_microState50FieldId, PeridigmField::STEP_NP1)->ExtractView(&state50NP1);
  dataManager.getData(m_microState51FieldId, PeridigmField::STEP_NP1)->ExtractView(&state51NP1);
  dataManager.getData(m_microState52FieldId, PeridigmField::STEP_NP1)->ExtractView(&state52NP1);
  dataManager.getData(m_microState53FieldId, PeridigmField::STEP_NP1)->ExtractView(&state53NP1);
  dataManager.getData(m_microState54FieldId, PeridigmField::STEP_NP1)->ExtractView(&state54NP1);
  dataManager.getData(m_microState55FieldId, PeridigmField::STEP_NP1)->ExtractView(&state55NP1);
  dataManager.getData(m_microState56FieldId, PeridigmField::STEP_NP1)->ExtractView(&state56NP1);
  dataManager.getData(m_microState57FieldId, PeridigmField::STEP_NP1)->ExtractView(&state57NP1);
  dataManager.getData(m_microState58FieldId, PeridigmField::STEP_NP1)->ExtractView(&state58NP1);
  dataManager.getData(m_microState59FieldId, PeridigmField::STEP_NP1)->ExtractView(&state59NP1);
  dataManager.getData(m_microState60FieldId, PeridigmField::STEP_NP1)->ExtractView(&state60NP1);
  dataManager.getData(m_microState61FieldId, PeridigmField::STEP_NP1)->ExtractView(&state61NP1);
  dataManager.getData(m_microState62FieldId, PeridigmField::STEP_NP1)->ExtractView(&state62NP1);
  dataManager.getData(m_microState63FieldId, PeridigmField::STEP_NP1)->ExtractView(&state63NP1);
  dataManager.getData(m_microState64FieldId, PeridigmField::STEP_NP1)->ExtractView(&state64NP1);
  dataManager.getData(m_microState65FieldId, PeridigmField::STEP_NP1)->ExtractView(&state65NP1);
  dataManager.getData(m_microState66FieldId, PeridigmField::STEP_NP1)->ExtractView(&state66NP1);
  dataManager.getData(m_microState67FieldId, PeridigmField::STEP_NP1)->ExtractView(&state67NP1);
  dataManager.getData(m_microState68FieldId, PeridigmField::STEP_NP1)->ExtractView(&state68NP1);
  dataManager.getData(m_microState69FieldId, PeridigmField::STEP_NP1)->ExtractView(&state69NP1);
  dataManager.getData(m_microState70FieldId, PeridigmField::STEP_NP1)->ExtractView(&state70NP1);
  dataManager.getData(m_microState71FieldId, PeridigmField::STEP_NP1)->ExtractView(&state71NP1);
  dataManager.getData(m_microState72FieldId, PeridigmField::STEP_NP1)->ExtractView(&state72NP1);
  dataManager.getData(m_microState73FieldId, PeridigmField::STEP_NP1)->ExtractView(&state73NP1);
  dataManager.getData(m_microState74FieldId, PeridigmField::STEP_NP1)->ExtractView(&state74NP1);
  dataManager.getData(m_microState75FieldId, PeridigmField::STEP_NP1)->ExtractView(&state75NP1);
  dataManager.getData(m_microState76FieldId, PeridigmField::STEP_NP1)->ExtractView(&state76NP1);
  dataManager.getData(m_microState77FieldId, PeridigmField::STEP_NP1)->ExtractView(&state77NP1);
  dataManager.getData(m_microState78FieldId, PeridigmField::STEP_NP1)->ExtractView(&state78NP1);
  dataManager.getData(m_microState79FieldId, PeridigmField::STEP_NP1)->ExtractView(&state79NP1);
  dataManager.getData(m_microState80FieldId, PeridigmField::STEP_NP1)->ExtractView(&state80NP1);
  dataManager.getData(m_microState81FieldId, PeridigmField::STEP_NP1)->ExtractView(&state81NP1);
  dataManager.getData(m_microState82FieldId, PeridigmField::STEP_NP1)->ExtractView(&state82NP1);
  dataManager.getData(m_microState83FieldId, PeridigmField::STEP_NP1)->ExtractView(&state83NP1);
  dataManager.getData(m_microState84FieldId, PeridigmField::STEP_NP1)->ExtractView(&state84NP1);
  dataManager.getData(m_microState85FieldId, PeridigmField::STEP_NP1)->ExtractView(&state85NP1);
  dataManager.getData(m_microState86FieldId, PeridigmField::STEP_NP1)->ExtractView(&state86NP1);
  dataManager.getData(m_microState87FieldId, PeridigmField::STEP_NP1)->ExtractView(&state87NP1);
  dataManager.getData(m_microState88FieldId, PeridigmField::STEP_NP1)->ExtractView(&state88NP1);
  dataManager.getData(m_microState89FieldId, PeridigmField::STEP_NP1)->ExtractView(&state89NP1);
  dataManager.getData(m_microState90FieldId, PeridigmField::STEP_NP1)->ExtractView(&state90NP1);
  dataManager.getData(m_microState91FieldId, PeridigmField::STEP_NP1)->ExtractView(&state91NP1);
  dataManager.getData(m_microState92FieldId, PeridigmField::STEP_NP1)->ExtractView(&state92NP1);
  dataManager.getData(m_microState93FieldId, PeridigmField::STEP_NP1)->ExtractView(&state93NP1);
  dataManager.getData(m_microState94FieldId, PeridigmField::STEP_NP1)->ExtractView(&state94NP1);
  dataManager.getData(m_microState95FieldId, PeridigmField::STEP_NP1)->ExtractView(&state95NP1);
  dataManager.getData(m_microState96FieldId, PeridigmField::STEP_NP1)->ExtractView(&state96NP1);
  dataManager.getData(m_microState97FieldId, PeridigmField::STEP_NP1)->ExtractView(&state97NP1);
  dataManager.getData(m_microState98FieldId, PeridigmField::STEP_NP1)->ExtractView(&state98NP1);
  dataManager.getData(m_microState99FieldId, PeridigmField::STEP_NP1)->ExtractView(&state99NP1);
  dataManager.getData(m_microState100FieldId, PeridigmField::STEP_NP1)->ExtractView(&state100NP1);
  dataManager.getData(m_microState101FieldId, PeridigmField::STEP_NP1)->ExtractView(&state101NP1);
  dataManager.getData(m_microState102FieldId, PeridigmField::STEP_NP1)->ExtractView(&state102NP1);
  dataManager.getData(m_microState103FieldId, PeridigmField::STEP_NP1)->ExtractView(&state103NP1);
  dataManager.getData(m_microState104FieldId, PeridigmField::STEP_NP1)->ExtractView(&state104NP1);
  dataManager.getData(m_microState105FieldId, PeridigmField::STEP_NP1)->ExtractView(&state105NP1);
  dataManager.getData(m_microState106FieldId, PeridigmField::STEP_NP1)->ExtractView(&state106NP1);
  dataManager.getData(m_microState107FieldId, PeridigmField::STEP_NP1)->ExtractView(&state107NP1);
  dataManager.getData(m_microState108FieldId, PeridigmField::STEP_NP1)->ExtractView(&state108NP1);
  dataManager.getData(m_microState109FieldId, PeridigmField::STEP_NP1)->ExtractView(&state109NP1);
  dataManager.getData(m_microState110FieldId, PeridigmField::STEP_NP1)->ExtractView(&state110NP1);
  dataManager.getData(m_microState111FieldId, PeridigmField::STEP_NP1)->ExtractView(&state111NP1);
  dataManager.getData(m_microState112FieldId, PeridigmField::STEP_NP1)->ExtractView(&state112NP1);
  dataManager.getData(m_microState113FieldId, PeridigmField::STEP_NP1)->ExtractView(&state113NP1);
  dataManager.getData(m_microState114FieldId, PeridigmField::STEP_NP1)->ExtractView(&state114NP1);
  dataManager.getData(m_microState115FieldId, PeridigmField::STEP_NP1)->ExtractView(&state115NP1);
  dataManager.getData(m_microState116FieldId, PeridigmField::STEP_NP1)->ExtractView(&state116NP1);
  dataManager.getData(m_microState117FieldId, PeridigmField::STEP_NP1)->ExtractView(&state117NP1);
  dataManager.getData(m_microState118FieldId, PeridigmField::STEP_NP1)->ExtractView(&state118NP1);
  dataManager.getData(m_microState119FieldId, PeridigmField::STEP_NP1)->ExtractView(&state119NP1);
  dataManager.getData(m_microState120FieldId, PeridigmField::STEP_NP1)->ExtractView(&state120NP1);
  dataManager.getData(m_microState121FieldId, PeridigmField::STEP_NP1)->ExtractView(&state121NP1);
  dataManager.getData(m_microState122FieldId, PeridigmField::STEP_NP1)->ExtractView(&state122NP1);
  dataManager.getData(m_microState123FieldId, PeridigmField::STEP_NP1)->ExtractView(&state123NP1);
  dataManager.getData(m_microState124FieldId, PeridigmField::STEP_NP1)->ExtractView(&state124NP1);
  dataManager.getData(m_microState125FieldId, PeridigmField::STEP_NP1)->ExtractView(&state125NP1);
  dataManager.getData(m_microState126FieldId, PeridigmField::STEP_NP1)->ExtractView(&state126NP1);
  dataManager.getData(m_microState127FieldId, PeridigmField::STEP_NP1)->ExtractView(&state127NP1);
  dataManager.getData(m_microState128FieldId, PeridigmField::STEP_NP1)->ExtractView(&state128NP1);
  dataManager.getData(m_microState129FieldId, PeridigmField::STEP_NP1)->ExtractView(&state129NP1);
  dataManager.getData(m_microState130FieldId, PeridigmField::STEP_NP1)->ExtractView(&state130NP1);
  dataManager.getData(m_microState131FieldId, PeridigmField::STEP_NP1)->ExtractView(&state131NP1);
  dataManager.getData(m_microState132FieldId, PeridigmField::STEP_NP1)->ExtractView(&state132NP1);
  dataManager.getData(m_microState133FieldId, PeridigmField::STEP_NP1)->ExtractView(&state133NP1);
  dataManager.getData(m_microState134FieldId, PeridigmField::STEP_NP1)->ExtractView(&state134NP1);
  dataManager.getData(m_microState135FieldId, PeridigmField::STEP_NP1)->ExtractView(&state135NP1);
  dataManager.getData(m_microState136FieldId, PeridigmField::STEP_NP1)->ExtractView(&state136NP1);
  dataManager.getData(m_microState137FieldId, PeridigmField::STEP_NP1)->ExtractView(&state137NP1);
  dataManager.getData(m_microState138FieldId, PeridigmField::STEP_NP1)->ExtractView(&state138NP1);
  dataManager.getData(m_microState139FieldId, PeridigmField::STEP_NP1)->ExtractView(&state139NP1);
  dataManager.getData(m_microState140FieldId, PeridigmField::STEP_NP1)->ExtractView(&state140NP1);
  dataManager.getData(m_microState141FieldId, PeridigmField::STEP_NP1)->ExtractView(&state141NP1);
  dataManager.getData(m_microState142FieldId, PeridigmField::STEP_NP1)->ExtractView(&state142NP1);
  dataManager.getData(m_microState143FieldId, PeridigmField::STEP_NP1)->ExtractView(&state143NP1);
  dataManager.getData(m_microState144FieldId, PeridigmField::STEP_NP1)->ExtractView(&state144NP1);
  dataManager.getData(m_microState145FieldId, PeridigmField::STEP_NP1)->ExtractView(&state145NP1);
  dataManager.getData(m_microState146FieldId, PeridigmField::STEP_NP1)->ExtractView(&state146NP1);
  dataManager.getData(m_microState147FieldId, PeridigmField::STEP_NP1)->ExtractView(&state147NP1);
  dataManager.getData(m_microState148FieldId, PeridigmField::STEP_NP1)->ExtractView(&state148NP1);
  dataManager.getData(m_microState149FieldId, PeridigmField::STEP_NP1)->ExtractView(&state149NP1);
  dataManager.getData(m_microState150FieldId, PeridigmField::STEP_NP1)->ExtractView(&state150NP1);
  dataManager.getData(m_microState151FieldId, PeridigmField::STEP_NP1)->ExtractView(&state151NP1);
  dataManager.getData(m_microState152FieldId, PeridigmField::STEP_NP1)->ExtractView(&state152NP1);
  dataManager.getData(m_microState153FieldId, PeridigmField::STEP_NP1)->ExtractView(&state153NP1);
  dataManager.getData(m_microState154FieldId, PeridigmField::STEP_NP1)->ExtractView(&state154NP1);
  dataManager.getData(m_microState155FieldId, PeridigmField::STEP_NP1)->ExtractView(&state155NP1);
  dataManager.getData(m_microState156FieldId, PeridigmField::STEP_NP1)->ExtractView(&state156NP1);
  dataManager.getData(m_microState157FieldId, PeridigmField::STEP_NP1)->ExtractView(&state157NP1);
  dataManager.getData(m_microState158FieldId, PeridigmField::STEP_NP1)->ExtractView(&state158NP1);
  dataManager.getData(m_microState159FieldId, PeridigmField::STEP_NP1)->ExtractView(&state159NP1);
  dataManager.getData(m_microState160FieldId, PeridigmField::STEP_NP1)->ExtractView(&state160NP1);
  dataManager.getData(m_microState161FieldId, PeridigmField::STEP_NP1)->ExtractView(&state161NP1);
  dataManager.getData(m_microState162FieldId, PeridigmField::STEP_NP1)->ExtractView(&state162NP1);
  dataManager.getData(m_microState163FieldId, PeridigmField::STEP_NP1)->ExtractView(&state163NP1);
  dataManager.getData(m_microState164FieldId, PeridigmField::STEP_NP1)->ExtractView(&state164NP1);
  dataManager.getData(m_microState165FieldId, PeridigmField::STEP_NP1)->ExtractView(&state165NP1);
  dataManager.getData(m_microState166FieldId, PeridigmField::STEP_NP1)->ExtractView(&state166NP1);
  dataManager.getData(m_microState167FieldId, PeridigmField::STEP_NP1)->ExtractView(&state167NP1);
  dataManager.getData(m_microState168FieldId, PeridigmField::STEP_NP1)->ExtractView(&state168NP1);
  dataManager.getData(m_microState169FieldId, PeridigmField::STEP_NP1)->ExtractView(&state169NP1);
  dataManager.getData(m_microState170FieldId, PeridigmField::STEP_NP1)->ExtractView(&state170NP1);
  dataManager.getData(m_microState171FieldId, PeridigmField::STEP_NP1)->ExtractView(&state171NP1);
  dataManager.getData(m_microState172FieldId, PeridigmField::STEP_NP1)->ExtractView(&state172NP1);
  dataManager.getData(m_microState173FieldId, PeridigmField::STEP_NP1)->ExtractView(&state173NP1);
  dataManager.getData(m_microState174FieldId, PeridigmField::STEP_NP1)->ExtractView(&state174NP1);
  dataManager.getData(m_microState175FieldId, PeridigmField::STEP_NP1)->ExtractView(&state175NP1);
  dataManager.getData(m_microState176FieldId, PeridigmField::STEP_NP1)->ExtractView(&state176NP1);
  dataManager.getData(m_microState177FieldId, PeridigmField::STEP_NP1)->ExtractView(&state177NP1);
  dataManager.getData(m_microState178FieldId, PeridigmField::STEP_NP1)->ExtractView(&state178NP1);
  dataManager.getData(m_microState179FieldId, PeridigmField::STEP_NP1)->ExtractView(&state179NP1);
  dataManager.getData(m_microState180FieldId, PeridigmField::STEP_NP1)->ExtractView(&state180NP1);
  dataManager.getData(m_microState181FieldId, PeridigmField::STEP_NP1)->ExtractView(&state181NP1);
  dataManager.getData(m_microState182FieldId, PeridigmField::STEP_NP1)->ExtractView(&state182NP1);
  dataManager.getData(m_microState183FieldId, PeridigmField::STEP_NP1)->ExtractView(&state183NP1);
  dataManager.getData(m_microState184FieldId, PeridigmField::STEP_NP1)->ExtractView(&state184NP1);
  dataManager.getData(m_microState185FieldId, PeridigmField::STEP_NP1)->ExtractView(&state185NP1);
  dataManager.getData(m_microState186FieldId, PeridigmField::STEP_NP1)->ExtractView(&state186NP1);
  dataManager.getData(m_microState187FieldId, PeridigmField::STEP_NP1)->ExtractView(&state187NP1);
  dataManager.getData(m_microState188FieldId, PeridigmField::STEP_NP1)->ExtractView(&state188NP1);
  dataManager.getData(m_microState189FieldId, PeridigmField::STEP_NP1)->ExtractView(&state189NP1);

  double *internalEnergyN, *internalEnergyNP1;
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&internalEnergyN);
  dataManager.getData(m_internalEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&internalEnergyNP1);

  double *inelasticEnergyN, *inelasticEnergyNP1;
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&inelasticEnergyN);
  dataManager.getData(m_inelasticEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&inelasticEnergyNP1);

  CORRESPONDENCE::updateMicroplaneM7Stress(strainRate, 
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
                                           state50N,
                                           state51N,
                                           state52N,
                                           state53N,
                                           state54N,
                                           state55N,
                                           state56N,
                                           state57N,
                                           state58N,
                                           state59N,
                                           state60N,
                                           state61N,
                                           state62N,
                                           state63N,
                                           state64N,
                                           state65N,
                                           state66N,
                                           state67N,
                                           state68N,
                                           state69N,
                                           state70N,
                                           state71N,
                                           state72N,
                                           state73N,
                                           state74N,
                                           state75N,
                                           state76N,
                                           state77N,
                                           state78N,
                                           state79N,
                                           state80N,
                                           state81N,
                                           state82N,
                                           state83N,
                                           state84N,
                                           state85N,
                                           state86N,
                                           state87N,
                                           state88N,
                                           state89N,
                                           state90N,
                                           state91N,
                                           state92N,
                                           state93N,
                                           state94N,
                                           state95N,
                                           state96N,
                                           state97N,
                                           state98N,
                                           state99N,
                                           state100N,
                                           state101N,
                                           state102N,
                                           state103N,
                                           state104N,
                                           state105N,
                                           state106N,
                                           state107N,
                                           state108N,
                                           state109N,
                                           state110N,
                                           state111N,
                                           state112N,
                                           state113N,
                                           state114N,
                                           state115N,
                                           state116N,
                                           state117N,
                                           state118N,
                                           state119N,
                                           state120N,
                                           state121N,
                                           state122N,
                                           state123N,
                                           state124N,
                                           state125N,
                                           state126N,
                                           state127N,
                                           state128N,
                                           state129N,
                                           state130N,
                                           state131N,
                                           state132N,
                                           state133N,
                                           state134N,
                                           state135N,
                                           state136N,
                                           state137N,
                                           state138N,
                                           state139N,
                                           state140N,
                                           state141N,
                                           state142N,
                                           state143N,
                                           state144N,
                                           state145N,
                                           state146N,
                                           state147N,
                                           state148N,
                                           state149N,
                                           state150N,
                                           state151N,
                                           state152N,
                                           state153N,
                                           state154N,
                                           state155N,
                                           state156N,
                                           state157N,
                                           state158N,
                                           state159N,
                                           state160N,
                                           state161N,
                                           state162N,
                                           state163N,
                                           state164N,
                                           state165N,
                                           state166N,
                                           state167N,
                                           state168N,
                                           state169N,
                                           state170N,
                                           state171N,
                                           state172N,
                                           state173N,
                                           state174N,
                                           state175N,
                                           state176N,
                                           state177N,
                                           state178N,
                                           state179N,
                                           state180N,
                                           state181N,
                                           state182N,
                                           state183N,
                                           state184N,
                                           state185N,
                                           state186N,
                                           state187N,
                                           state188N,
                                           state189N,
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
                                           state50NP1,
                                           state51NP1,
                                           state52NP1,
                                           state53NP1,
                                           state54NP1,
                                           state55NP1,
                                           state56NP1,
                                           state57NP1,
                                           state58NP1,
                                           state59NP1,
                                           state60NP1,
                                           state61NP1,
                                           state62NP1,
                                           state63NP1,
                                           state64NP1,
                                           state65NP1,
                                           state66NP1,
                                           state67NP1,
                                           state68NP1,
                                           state69NP1,
                                           state70NP1,
                                           state71NP1,
                                           state72NP1,
                                           state73NP1,
                                           state74NP1,
                                           state75NP1,
                                           state76NP1,
                                           state77NP1,
                                           state78NP1,
                                           state79NP1,
                                           state80NP1,
                                           state81NP1,
                                           state82NP1,
                                           state83NP1,
                                           state84NP1,
                                           state85NP1,
                                           state86NP1,
                                           state87NP1,
                                           state88NP1,
                                           state89NP1,
                                           state90NP1,
                                           state91NP1,
                                           state92NP1,
                                           state93NP1,
                                           state94NP1,
                                           state95NP1,
                                           state96NP1,
                                           state97NP1,
                                           state98NP1,
                                           state99NP1,
                                           state100NP1,
                                           state101NP1,
                                           state102NP1,
                                           state103NP1,
                                           state104NP1,
                                           state105NP1,
                                           state106NP1,
                                           state107NP1,
                                           state108NP1,
                                           state109NP1,
                                           state110NP1,
                                           state111NP1,
                                           state112NP1,
                                           state113NP1,
                                           state114NP1,
                                           state115NP1,
                                           state116NP1,
                                           state117NP1,
                                           state118NP1,
                                           state119NP1,
                                           state120NP1,
                                           state121NP1,
                                           state122NP1,
                                           state123NP1,
                                           state124NP1,
                                           state125NP1,
                                           state126NP1,
                                           state127NP1,
                                           state128NP1,
                                           state129NP1,
                                           state130NP1,
                                           state131NP1,
                                           state132NP1,
                                           state133NP1,
                                           state134NP1,
                                           state135NP1,
                                           state136NP1,
                                           state137NP1,
                                           state138NP1,
                                           state139NP1,
                                           state140NP1,
                                           state141NP1,
                                           state142NP1,
                                           state143NP1,
                                           state144NP1,
                                           state145NP1,
                                           state146NP1,
                                           state147NP1,
                                           state148NP1,
                                           state149NP1,
                                           state150NP1,
                                           state151NP1,
                                           state152NP1,
                                           state153NP1,
                                           state154NP1,
                                           state155NP1,
                                           state156NP1,
                                           state157NP1,
                                           state158NP1,
                                           state159NP1,
                                           state160NP1,
                                           state161NP1,
                                           state162NP1,
                                           state163NP1,
                                           state164NP1,
                                           state165NP1,
                                           state166NP1,
                                           state167NP1,
                                           state168NP1,
                                           state169NP1,
                                           state170NP1,
                                           state171NP1,
                                           state172NP1,
                                           state173NP1,
                                           state174NP1,
                                           state175NP1,
                                           state176NP1,
                                           state177NP1,
                                           state178NP1,
                                           state179NP1,
                                           state180NP1,
                                           state181NP1,
                                           state182NP1,
                                           state183NP1,
                                           state184NP1,
                                           state185NP1,
                                           state186NP1,
                                           state187NP1,
                                           state188NP1,
                                           state189NP1,
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
  double *bondLevelState50N;
  double *bondLevelState51N;
  double *bondLevelState52N;
  double *bondLevelState53N;
  double *bondLevelState54N;
  double *bondLevelState55N;
  double *bondLevelState56N;
  double *bondLevelState57N;
  double *bondLevelState58N;
  double *bondLevelState59N;
  double *bondLevelState60N;
  double *bondLevelState61N;
  double *bondLevelState62N;
  double *bondLevelState63N;
  double *bondLevelState64N;
  double *bondLevelState65N;
  double *bondLevelState66N;
  double *bondLevelState67N;
  double *bondLevelState68N;
  double *bondLevelState69N;
  double *bondLevelState70N;
  double *bondLevelState71N;
  double *bondLevelState72N;
  double *bondLevelState73N;
  double *bondLevelState74N;
  double *bondLevelState75N;
  double *bondLevelState76N;
  double *bondLevelState77N;
  double *bondLevelState78N;
  double *bondLevelState79N;
  double *bondLevelState80N;
  double *bondLevelState81N;
  double *bondLevelState82N;
  double *bondLevelState83N;
  double *bondLevelState84N;
  double *bondLevelState85N;
  double *bondLevelState86N;
  double *bondLevelState87N;
  double *bondLevelState88N;
  double *bondLevelState89N;
  double *bondLevelState90N;
  double *bondLevelState91N;
  double *bondLevelState92N;
  double *bondLevelState93N;
  double *bondLevelState94N;
  double *bondLevelState95N;
  double *bondLevelState96N;
  double *bondLevelState97N;
  double *bondLevelState98N;
  double *bondLevelState99N;
  double *bondLevelState100N;
  double *bondLevelState101N;
  double *bondLevelState102N;
  double *bondLevelState103N;
  double *bondLevelState104N;
  double *bondLevelState105N;
  double *bondLevelState106N;
  double *bondLevelState107N;
  double *bondLevelState108N;
  double *bondLevelState109N;
  double *bondLevelState110N;
  double *bondLevelState111N;
  double *bondLevelState112N;
  double *bondLevelState113N;
  double *bondLevelState114N;
  double *bondLevelState115N;
  double *bondLevelState116N;
  double *bondLevelState117N;
  double *bondLevelState118N;
  double *bondLevelState119N;
  double *bondLevelState120N;
  double *bondLevelState121N;
  double *bondLevelState122N;
  double *bondLevelState123N;
  double *bondLevelState124N;
  double *bondLevelState125N;
  double *bondLevelState126N;
  double *bondLevelState127N;
  double *bondLevelState128N;
  double *bondLevelState129N;
  double *bondLevelState130N;
  double *bondLevelState131N;
  double *bondLevelState132N;
  double *bondLevelState133N;
  double *bondLevelState134N;
  double *bondLevelState135N;
  double *bondLevelState136N;
  double *bondLevelState137N;
  double *bondLevelState138N;
  double *bondLevelState139N;
  double *bondLevelState140N;
  double *bondLevelState141N;
  double *bondLevelState142N;
  double *bondLevelState143N;
  double *bondLevelState144N;
  double *bondLevelState145N;
  double *bondLevelState146N;
  double *bondLevelState147N;
  double *bondLevelState148N;
  double *bondLevelState149N;
  double *bondLevelState150N;
  double *bondLevelState151N;
  double *bondLevelState152N;
  double *bondLevelState153N;
  double *bondLevelState154N;
  double *bondLevelState155N;
  double *bondLevelState156N;
  double *bondLevelState157N;
  double *bondLevelState158N;
  double *bondLevelState159N;
  double *bondLevelState160N;
  double *bondLevelState161N;
  double *bondLevelState162N;
  double *bondLevelState163N;
  double *bondLevelState164N;
  double *bondLevelState165N;
  double *bondLevelState166N;
  double *bondLevelState167N;
  double *bondLevelState168N;
  double *bondLevelState169N;
  double *bondLevelState170N;
  double *bondLevelState171N;
  double *bondLevelState172N;
  double *bondLevelState173N;
  double *bondLevelState174N;
  double *bondLevelState175N;
  double *bondLevelState176N;
  double *bondLevelState177N;
  double *bondLevelState178N;
  double *bondLevelState179N;
  double *bondLevelState180N;
  double *bondLevelState181N;
  double *bondLevelState182N;
  double *bondLevelState183N;
  double *bondLevelState184N;
  double *bondLevelState185N;
  double *bondLevelState186N;
  double *bondLevelState187N;
  double *bondLevelState188N;
  double *bondLevelState189N;
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
  dataManager.getData(m_bondLevelMicroState50FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState50N);
  dataManager.getData(m_bondLevelMicroState51FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState51N);
  dataManager.getData(m_bondLevelMicroState52FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState52N);
  dataManager.getData(m_bondLevelMicroState53FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState53N);
  dataManager.getData(m_bondLevelMicroState54FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState54N);
  dataManager.getData(m_bondLevelMicroState55FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState55N);
  dataManager.getData(m_bondLevelMicroState56FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState56N);
  dataManager.getData(m_bondLevelMicroState57FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState57N);
  dataManager.getData(m_bondLevelMicroState58FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState58N);
  dataManager.getData(m_bondLevelMicroState59FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState59N);
  dataManager.getData(m_bondLevelMicroState60FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState60N);
  dataManager.getData(m_bondLevelMicroState61FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState61N);
  dataManager.getData(m_bondLevelMicroState62FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState62N);
  dataManager.getData(m_bondLevelMicroState63FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState63N);
  dataManager.getData(m_bondLevelMicroState64FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState64N);
  dataManager.getData(m_bondLevelMicroState65FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState65N);
  dataManager.getData(m_bondLevelMicroState66FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState66N);
  dataManager.getData(m_bondLevelMicroState67FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState67N);
  dataManager.getData(m_bondLevelMicroState68FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState68N);
  dataManager.getData(m_bondLevelMicroState69FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState69N);
  dataManager.getData(m_bondLevelMicroState70FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState70N);
  dataManager.getData(m_bondLevelMicroState71FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState71N);
  dataManager.getData(m_bondLevelMicroState72FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState72N);
  dataManager.getData(m_bondLevelMicroState73FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState73N);
  dataManager.getData(m_bondLevelMicroState74FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState74N);
  dataManager.getData(m_bondLevelMicroState75FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState75N);
  dataManager.getData(m_bondLevelMicroState76FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState76N);
  dataManager.getData(m_bondLevelMicroState77FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState77N);
  dataManager.getData(m_bondLevelMicroState78FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState78N);
  dataManager.getData(m_bondLevelMicroState79FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState79N);
  dataManager.getData(m_bondLevelMicroState80FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState80N);
  dataManager.getData(m_bondLevelMicroState81FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState81N);
  dataManager.getData(m_bondLevelMicroState82FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState82N);
  dataManager.getData(m_bondLevelMicroState83FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState83N);
  dataManager.getData(m_bondLevelMicroState84FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState84N);
  dataManager.getData(m_bondLevelMicroState85FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState85N);
  dataManager.getData(m_bondLevelMicroState86FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState86N);
  dataManager.getData(m_bondLevelMicroState87FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState87N);
  dataManager.getData(m_bondLevelMicroState88FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState88N);
  dataManager.getData(m_bondLevelMicroState89FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState89N);
  dataManager.getData(m_bondLevelMicroState90FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState90N);
  dataManager.getData(m_bondLevelMicroState91FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState91N);
  dataManager.getData(m_bondLevelMicroState92FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState92N);
  dataManager.getData(m_bondLevelMicroState93FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState93N);
  dataManager.getData(m_bondLevelMicroState94FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState94N);
  dataManager.getData(m_bondLevelMicroState95FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState95N);
  dataManager.getData(m_bondLevelMicroState96FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState96N);
  dataManager.getData(m_bondLevelMicroState97FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState97N);
  dataManager.getData(m_bondLevelMicroState98FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState98N);
  dataManager.getData(m_bondLevelMicroState99FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState99N);
  dataManager.getData(m_bondLevelMicroState100FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState100N);
  dataManager.getData(m_bondLevelMicroState101FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState101N);
  dataManager.getData(m_bondLevelMicroState102FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState102N);
  dataManager.getData(m_bondLevelMicroState103FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState103N);
  dataManager.getData(m_bondLevelMicroState104FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState104N);
  dataManager.getData(m_bondLevelMicroState105FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState105N);
  dataManager.getData(m_bondLevelMicroState106FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState106N);
  dataManager.getData(m_bondLevelMicroState107FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState107N);
  dataManager.getData(m_bondLevelMicroState108FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState108N);
  dataManager.getData(m_bondLevelMicroState109FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState109N);
  dataManager.getData(m_bondLevelMicroState110FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState110N);
  dataManager.getData(m_bondLevelMicroState111FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState111N);
  dataManager.getData(m_bondLevelMicroState112FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState112N);
  dataManager.getData(m_bondLevelMicroState113FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState113N);
  dataManager.getData(m_bondLevelMicroState114FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState114N);
  dataManager.getData(m_bondLevelMicroState115FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState115N);
  dataManager.getData(m_bondLevelMicroState116FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState116N);
  dataManager.getData(m_bondLevelMicroState117FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState117N);
  dataManager.getData(m_bondLevelMicroState118FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState118N);
  dataManager.getData(m_bondLevelMicroState119FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState119N);
  dataManager.getData(m_bondLevelMicroState120FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState120N);
  dataManager.getData(m_bondLevelMicroState121FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState121N);
  dataManager.getData(m_bondLevelMicroState122FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState122N);
  dataManager.getData(m_bondLevelMicroState123FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState123N);
  dataManager.getData(m_bondLevelMicroState124FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState124N);
  dataManager.getData(m_bondLevelMicroState125FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState125N);
  dataManager.getData(m_bondLevelMicroState126FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState126N);
  dataManager.getData(m_bondLevelMicroState127FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState127N);
  dataManager.getData(m_bondLevelMicroState128FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState128N);
  dataManager.getData(m_bondLevelMicroState129FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState129N);
  dataManager.getData(m_bondLevelMicroState130FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState130N);
  dataManager.getData(m_bondLevelMicroState131FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState131N);
  dataManager.getData(m_bondLevelMicroState132FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState132N);
  dataManager.getData(m_bondLevelMicroState133FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState133N);
  dataManager.getData(m_bondLevelMicroState134FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState134N);
  dataManager.getData(m_bondLevelMicroState135FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState135N);
  dataManager.getData(m_bondLevelMicroState136FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState136N);
  dataManager.getData(m_bondLevelMicroState137FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState137N);
  dataManager.getData(m_bondLevelMicroState138FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState138N);
  dataManager.getData(m_bondLevelMicroState139FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState139N);
  dataManager.getData(m_bondLevelMicroState140FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState140N);
  dataManager.getData(m_bondLevelMicroState141FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState141N);
  dataManager.getData(m_bondLevelMicroState142FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState142N);
  dataManager.getData(m_bondLevelMicroState143FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState143N);
  dataManager.getData(m_bondLevelMicroState144FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState144N);
  dataManager.getData(m_bondLevelMicroState145FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState145N);
  dataManager.getData(m_bondLevelMicroState146FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState146N);
  dataManager.getData(m_bondLevelMicroState147FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState147N);
  dataManager.getData(m_bondLevelMicroState148FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState148N);
  dataManager.getData(m_bondLevelMicroState149FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState149N);
  dataManager.getData(m_bondLevelMicroState150FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState150N);
  dataManager.getData(m_bondLevelMicroState151FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState151N);
  dataManager.getData(m_bondLevelMicroState152FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState152N);
  dataManager.getData(m_bondLevelMicroState153FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState153N);
  dataManager.getData(m_bondLevelMicroState154FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState154N);
  dataManager.getData(m_bondLevelMicroState155FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState155N);
  dataManager.getData(m_bondLevelMicroState156FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState156N);
  dataManager.getData(m_bondLevelMicroState157FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState157N);
  dataManager.getData(m_bondLevelMicroState158FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState158N);
  dataManager.getData(m_bondLevelMicroState159FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState159N);
  dataManager.getData(m_bondLevelMicroState160FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState160N);
  dataManager.getData(m_bondLevelMicroState161FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState161N);
  dataManager.getData(m_bondLevelMicroState162FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState162N);
  dataManager.getData(m_bondLevelMicroState163FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState163N);
  dataManager.getData(m_bondLevelMicroState164FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState164N);
  dataManager.getData(m_bondLevelMicroState165FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState165N);
  dataManager.getData(m_bondLevelMicroState166FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState166N);
  dataManager.getData(m_bondLevelMicroState167FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState167N);
  dataManager.getData(m_bondLevelMicroState168FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState168N);
  dataManager.getData(m_bondLevelMicroState169FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState169N);
  dataManager.getData(m_bondLevelMicroState170FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState170N);
  dataManager.getData(m_bondLevelMicroState171FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState171N);
  dataManager.getData(m_bondLevelMicroState172FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState172N);
  dataManager.getData(m_bondLevelMicroState173FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState173N);
  dataManager.getData(m_bondLevelMicroState174FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState174N);
  dataManager.getData(m_bondLevelMicroState175FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState175N);
  dataManager.getData(m_bondLevelMicroState176FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState176N);
  dataManager.getData(m_bondLevelMicroState177FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState177N);
  dataManager.getData(m_bondLevelMicroState178FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState178N);
  dataManager.getData(m_bondLevelMicroState179FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState179N);
  dataManager.getData(m_bondLevelMicroState180FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState180N);
  dataManager.getData(m_bondLevelMicroState181FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState181N);
  dataManager.getData(m_bondLevelMicroState182FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState182N);
  dataManager.getData(m_bondLevelMicroState183FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState183N);
  dataManager.getData(m_bondLevelMicroState184FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState184N);
  dataManager.getData(m_bondLevelMicroState185FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState185N);
  dataManager.getData(m_bondLevelMicroState186FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState186N);
  dataManager.getData(m_bondLevelMicroState187FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState187N);
  dataManager.getData(m_bondLevelMicroState188FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState188N);
  dataManager.getData(m_bondLevelMicroState189FieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelState189N);

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
  double *bondLevelState50NP1;
  double *bondLevelState51NP1;
  double *bondLevelState52NP1;
  double *bondLevelState53NP1;
  double *bondLevelState54NP1;
  double *bondLevelState55NP1;
  double *bondLevelState56NP1;
  double *bondLevelState57NP1;
  double *bondLevelState58NP1;
  double *bondLevelState59NP1;
  double *bondLevelState60NP1;
  double *bondLevelState61NP1;
  double *bondLevelState62NP1;
  double *bondLevelState63NP1;
  double *bondLevelState64NP1;
  double *bondLevelState65NP1;
  double *bondLevelState66NP1;
  double *bondLevelState67NP1;
  double *bondLevelState68NP1;
  double *bondLevelState69NP1;
  double *bondLevelState70NP1;
  double *bondLevelState71NP1;
  double *bondLevelState72NP1;
  double *bondLevelState73NP1;
  double *bondLevelState74NP1;
  double *bondLevelState75NP1;
  double *bondLevelState76NP1;
  double *bondLevelState77NP1;
  double *bondLevelState78NP1;
  double *bondLevelState79NP1;
  double *bondLevelState80NP1;
  double *bondLevelState81NP1;
  double *bondLevelState82NP1;
  double *bondLevelState83NP1;
  double *bondLevelState84NP1;
  double *bondLevelState85NP1;
  double *bondLevelState86NP1;
  double *bondLevelState87NP1;
  double *bondLevelState88NP1;
  double *bondLevelState89NP1;
  double *bondLevelState90NP1;
  double *bondLevelState91NP1;
  double *bondLevelState92NP1;
  double *bondLevelState93NP1;
  double *bondLevelState94NP1;
  double *bondLevelState95NP1;
  double *bondLevelState96NP1;
  double *bondLevelState97NP1;
  double *bondLevelState98NP1;
  double *bondLevelState99NP1;
  double *bondLevelState100NP1;
  double *bondLevelState101NP1;
  double *bondLevelState102NP1;
  double *bondLevelState103NP1;
  double *bondLevelState104NP1;
  double *bondLevelState105NP1;
  double *bondLevelState106NP1;
  double *bondLevelState107NP1;
  double *bondLevelState108NP1;
  double *bondLevelState109NP1;
  double *bondLevelState110NP1;
  double *bondLevelState111NP1;
  double *bondLevelState112NP1;
  double *bondLevelState113NP1;
  double *bondLevelState114NP1;
  double *bondLevelState115NP1;
  double *bondLevelState116NP1;
  double *bondLevelState117NP1;
  double *bondLevelState118NP1;
  double *bondLevelState119NP1;
  double *bondLevelState120NP1;
  double *bondLevelState121NP1;
  double *bondLevelState122NP1;
  double *bondLevelState123NP1;
  double *bondLevelState124NP1;
  double *bondLevelState125NP1;
  double *bondLevelState126NP1;
  double *bondLevelState127NP1;
  double *bondLevelState128NP1;
  double *bondLevelState129NP1;
  double *bondLevelState130NP1;
  double *bondLevelState131NP1;
  double *bondLevelState132NP1;
  double *bondLevelState133NP1;
  double *bondLevelState134NP1;
  double *bondLevelState135NP1;
  double *bondLevelState136NP1;
  double *bondLevelState137NP1;
  double *bondLevelState138NP1;
  double *bondLevelState139NP1;
  double *bondLevelState140NP1;
  double *bondLevelState141NP1;
  double *bondLevelState142NP1;
  double *bondLevelState143NP1;
  double *bondLevelState144NP1;
  double *bondLevelState145NP1;
  double *bondLevelState146NP1;
  double *bondLevelState147NP1;
  double *bondLevelState148NP1;
  double *bondLevelState149NP1;
  double *bondLevelState150NP1;
  double *bondLevelState151NP1;
  double *bondLevelState152NP1;
  double *bondLevelState153NP1;
  double *bondLevelState154NP1;
  double *bondLevelState155NP1;
  double *bondLevelState156NP1;
  double *bondLevelState157NP1;
  double *bondLevelState158NP1;
  double *bondLevelState159NP1;
  double *bondLevelState160NP1;
  double *bondLevelState161NP1;
  double *bondLevelState162NP1;
  double *bondLevelState163NP1;
  double *bondLevelState164NP1;
  double *bondLevelState165NP1;
  double *bondLevelState166NP1;
  double *bondLevelState167NP1;
  double *bondLevelState168NP1;
  double *bondLevelState169NP1;
  double *bondLevelState170NP1;
  double *bondLevelState171NP1;
  double *bondLevelState172NP1;
  double *bondLevelState173NP1;
  double *bondLevelState174NP1;
  double *bondLevelState175NP1;
  double *bondLevelState176NP1;
  double *bondLevelState177NP1;
  double *bondLevelState178NP1;
  double *bondLevelState179NP1;
  double *bondLevelState180NP1;
  double *bondLevelState181NP1;
  double *bondLevelState182NP1;
  double *bondLevelState183NP1;
  double *bondLevelState184NP1;
  double *bondLevelState185NP1;
  double *bondLevelState186NP1;
  double *bondLevelState187NP1;
  double *bondLevelState188NP1;
  double *bondLevelState189NP1;
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
  dataManager.getData(m_bondLevelMicroState50FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState50NP1);
  dataManager.getData(m_bondLevelMicroState51FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState51NP1);
  dataManager.getData(m_bondLevelMicroState52FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState52NP1);
  dataManager.getData(m_bondLevelMicroState53FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState53NP1);
  dataManager.getData(m_bondLevelMicroState54FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState54NP1);
  dataManager.getData(m_bondLevelMicroState55FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState55NP1);
  dataManager.getData(m_bondLevelMicroState56FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState56NP1);
  dataManager.getData(m_bondLevelMicroState57FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState57NP1);
  dataManager.getData(m_bondLevelMicroState58FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState58NP1);
  dataManager.getData(m_bondLevelMicroState59FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState59NP1);
  dataManager.getData(m_bondLevelMicroState60FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState60NP1);
  dataManager.getData(m_bondLevelMicroState61FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState61NP1);
  dataManager.getData(m_bondLevelMicroState62FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState62NP1);
  dataManager.getData(m_bondLevelMicroState63FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState63NP1);
  dataManager.getData(m_bondLevelMicroState64FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState64NP1);
  dataManager.getData(m_bondLevelMicroState65FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState65NP1);
  dataManager.getData(m_bondLevelMicroState66FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState66NP1);
  dataManager.getData(m_bondLevelMicroState67FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState67NP1);
  dataManager.getData(m_bondLevelMicroState68FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState68NP1);
  dataManager.getData(m_bondLevelMicroState69FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState69NP1);
  dataManager.getData(m_bondLevelMicroState70FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState70NP1);
  dataManager.getData(m_bondLevelMicroState71FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState71NP1);
  dataManager.getData(m_bondLevelMicroState72FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState72NP1);
  dataManager.getData(m_bondLevelMicroState73FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState73NP1);
  dataManager.getData(m_bondLevelMicroState74FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState74NP1);
  dataManager.getData(m_bondLevelMicroState75FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState75NP1);
  dataManager.getData(m_bondLevelMicroState76FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState76NP1);
  dataManager.getData(m_bondLevelMicroState77FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState77NP1);
  dataManager.getData(m_bondLevelMicroState78FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState78NP1);
  dataManager.getData(m_bondLevelMicroState79FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState79NP1);
  dataManager.getData(m_bondLevelMicroState80FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState80NP1);
  dataManager.getData(m_bondLevelMicroState81FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState81NP1);
  dataManager.getData(m_bondLevelMicroState82FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState82NP1);
  dataManager.getData(m_bondLevelMicroState83FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState83NP1);
  dataManager.getData(m_bondLevelMicroState84FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState84NP1);
  dataManager.getData(m_bondLevelMicroState85FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState85NP1);
  dataManager.getData(m_bondLevelMicroState86FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState86NP1);
  dataManager.getData(m_bondLevelMicroState87FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState87NP1);
  dataManager.getData(m_bondLevelMicroState88FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState88NP1);
  dataManager.getData(m_bondLevelMicroState89FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState89NP1);
  dataManager.getData(m_bondLevelMicroState90FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState90NP1);
  dataManager.getData(m_bondLevelMicroState91FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState91NP1);
  dataManager.getData(m_bondLevelMicroState92FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState92NP1);
  dataManager.getData(m_bondLevelMicroState93FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState93NP1);
  dataManager.getData(m_bondLevelMicroState94FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState94NP1);
  dataManager.getData(m_bondLevelMicroState95FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState95NP1);
  dataManager.getData(m_bondLevelMicroState96FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState96NP1);
  dataManager.getData(m_bondLevelMicroState97FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState97NP1);
  dataManager.getData(m_bondLevelMicroState98FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState98NP1);
  dataManager.getData(m_bondLevelMicroState99FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState99NP1);
  dataManager.getData(m_bondLevelMicroState100FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState100NP1);
  dataManager.getData(m_bondLevelMicroState101FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState101NP1);
  dataManager.getData(m_bondLevelMicroState102FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState102NP1);
  dataManager.getData(m_bondLevelMicroState103FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState103NP1);
  dataManager.getData(m_bondLevelMicroState104FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState104NP1);
  dataManager.getData(m_bondLevelMicroState105FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState105NP1);
  dataManager.getData(m_bondLevelMicroState106FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState106NP1);
  dataManager.getData(m_bondLevelMicroState107FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState107NP1);
  dataManager.getData(m_bondLevelMicroState108FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState108NP1);
  dataManager.getData(m_bondLevelMicroState109FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState109NP1);
  dataManager.getData(m_bondLevelMicroState110FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState110NP1);
  dataManager.getData(m_bondLevelMicroState111FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState111NP1);
  dataManager.getData(m_bondLevelMicroState112FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState112NP1);
  dataManager.getData(m_bondLevelMicroState113FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState113NP1);
  dataManager.getData(m_bondLevelMicroState114FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState114NP1);
  dataManager.getData(m_bondLevelMicroState115FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState115NP1);
  dataManager.getData(m_bondLevelMicroState116FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState116NP1);
  dataManager.getData(m_bondLevelMicroState117FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState117NP1);
  dataManager.getData(m_bondLevelMicroState118FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState118NP1);
  dataManager.getData(m_bondLevelMicroState119FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState119NP1);
  dataManager.getData(m_bondLevelMicroState120FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState120NP1);
  dataManager.getData(m_bondLevelMicroState121FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState121NP1);
  dataManager.getData(m_bondLevelMicroState122FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState122NP1);
  dataManager.getData(m_bondLevelMicroState123FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState123NP1);
  dataManager.getData(m_bondLevelMicroState124FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState124NP1);
  dataManager.getData(m_bondLevelMicroState125FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState125NP1);
  dataManager.getData(m_bondLevelMicroState126FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState126NP1);
  dataManager.getData(m_bondLevelMicroState127FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState127NP1);
  dataManager.getData(m_bondLevelMicroState128FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState128NP1);
  dataManager.getData(m_bondLevelMicroState129FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState129NP1);
  dataManager.getData(m_bondLevelMicroState130FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState130NP1);
  dataManager.getData(m_bondLevelMicroState131FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState131NP1);
  dataManager.getData(m_bondLevelMicroState132FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState132NP1);
  dataManager.getData(m_bondLevelMicroState133FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState133NP1);
  dataManager.getData(m_bondLevelMicroState134FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState134NP1);
  dataManager.getData(m_bondLevelMicroState135FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState135NP1);
  dataManager.getData(m_bondLevelMicroState136FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState136NP1);
  dataManager.getData(m_bondLevelMicroState137FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState137NP1);
  dataManager.getData(m_bondLevelMicroState138FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState138NP1);
  dataManager.getData(m_bondLevelMicroState139FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState139NP1);
  dataManager.getData(m_bondLevelMicroState140FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState140NP1);
  dataManager.getData(m_bondLevelMicroState141FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState141NP1);
  dataManager.getData(m_bondLevelMicroState142FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState142NP1);
  dataManager.getData(m_bondLevelMicroState143FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState143NP1);
  dataManager.getData(m_bondLevelMicroState144FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState144NP1);
  dataManager.getData(m_bondLevelMicroState145FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState145NP1);
  dataManager.getData(m_bondLevelMicroState146FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState146NP1);
  dataManager.getData(m_bondLevelMicroState147FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState147NP1);
  dataManager.getData(m_bondLevelMicroState148FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState148NP1);
  dataManager.getData(m_bondLevelMicroState149FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState149NP1);
  dataManager.getData(m_bondLevelMicroState150FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState150NP1);
  dataManager.getData(m_bondLevelMicroState151FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState151NP1);
  dataManager.getData(m_bondLevelMicroState152FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState152NP1);
  dataManager.getData(m_bondLevelMicroState153FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState153NP1);
  dataManager.getData(m_bondLevelMicroState154FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState154NP1);
  dataManager.getData(m_bondLevelMicroState155FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState155NP1);
  dataManager.getData(m_bondLevelMicroState156FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState156NP1);
  dataManager.getData(m_bondLevelMicroState157FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState157NP1);
  dataManager.getData(m_bondLevelMicroState158FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState158NP1);
  dataManager.getData(m_bondLevelMicroState159FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState159NP1);
  dataManager.getData(m_bondLevelMicroState160FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState160NP1);
  dataManager.getData(m_bondLevelMicroState161FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState161NP1);
  dataManager.getData(m_bondLevelMicroState162FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState162NP1);
  dataManager.getData(m_bondLevelMicroState163FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState163NP1);
  dataManager.getData(m_bondLevelMicroState164FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState164NP1);
  dataManager.getData(m_bondLevelMicroState165FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState165NP1);
  dataManager.getData(m_bondLevelMicroState166FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState166NP1);
  dataManager.getData(m_bondLevelMicroState167FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState167NP1);
  dataManager.getData(m_bondLevelMicroState168FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState168NP1);
  dataManager.getData(m_bondLevelMicroState169FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState169NP1);
  dataManager.getData(m_bondLevelMicroState170FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState170NP1);
  dataManager.getData(m_bondLevelMicroState171FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState171NP1);
  dataManager.getData(m_bondLevelMicroState172FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState172NP1);
  dataManager.getData(m_bondLevelMicroState173FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState173NP1);
  dataManager.getData(m_bondLevelMicroState174FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState174NP1);
  dataManager.getData(m_bondLevelMicroState175FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState175NP1);
  dataManager.getData(m_bondLevelMicroState176FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState176NP1);
  dataManager.getData(m_bondLevelMicroState177FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState177NP1);
  dataManager.getData(m_bondLevelMicroState178FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState178NP1);
  dataManager.getData(m_bondLevelMicroState179FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState179NP1);
  dataManager.getData(m_bondLevelMicroState180FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState180NP1);
  dataManager.getData(m_bondLevelMicroState181FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState181NP1);
  dataManager.getData(m_bondLevelMicroState182FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState182NP1);
  dataManager.getData(m_bondLevelMicroState183FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState183NP1);
  dataManager.getData(m_bondLevelMicroState184FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState184NP1);
  dataManager.getData(m_bondLevelMicroState185FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState185NP1);
  dataManager.getData(m_bondLevelMicroState186FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState186NP1);
  dataManager.getData(m_bondLevelMicroState187FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState187NP1);
  dataManager.getData(m_bondLevelMicroState188FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState188NP1);
  dataManager.getData(m_bondLevelMicroState189FieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelState189NP1);

  double *bondLevelInternalEnergyN, *bondLevelInternalEnergyNP1;
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelInternalEnergyN);
  dataManager.getData(m_bondLevelInternalEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelInternalEnergyNP1);

  double *bondLevelInelasticEnergyN, *bondLevelInelasticEnergyNP1;
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelInelasticEnergyN);
  dataManager.getData(m_bondLevelInelasticEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelInelasticEnergyNP1);

  CORRESPONDENCE::updateBondLevelMicroplaneM7Stress(bondLevelStrainRateXX,
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
                                                    bondLevelState50N,
                                                    bondLevelState51N,
                                                    bondLevelState52N,
                                                    bondLevelState53N,
                                                    bondLevelState54N,
                                                    bondLevelState55N,
                                                    bondLevelState56N,
                                                    bondLevelState57N,
                                                    bondLevelState58N,
                                                    bondLevelState59N,
                                                    bondLevelState60N,
                                                    bondLevelState61N,
                                                    bondLevelState62N,
                                                    bondLevelState63N,
                                                    bondLevelState64N,
                                                    bondLevelState65N,
                                                    bondLevelState66N,
                                                    bondLevelState67N,
                                                    bondLevelState68N,
                                                    bondLevelState69N,
                                                    bondLevelState70N,
                                                    bondLevelState71N,
                                                    bondLevelState72N,
                                                    bondLevelState73N,
                                                    bondLevelState74N,
                                                    bondLevelState75N,
                                                    bondLevelState76N,
                                                    bondLevelState77N,
                                                    bondLevelState78N,
                                                    bondLevelState79N,
                                                    bondLevelState80N,
                                                    bondLevelState81N,
                                                    bondLevelState82N,
                                                    bondLevelState83N,
                                                    bondLevelState84N,
                                                    bondLevelState85N,
                                                    bondLevelState86N,
                                                    bondLevelState87N,
                                                    bondLevelState88N,
                                                    bondLevelState89N,
                                                    bondLevelState90N,
                                                    bondLevelState91N,
                                                    bondLevelState92N,
                                                    bondLevelState93N,
                                                    bondLevelState94N,
                                                    bondLevelState95N,
                                                    bondLevelState96N,
                                                    bondLevelState97N,
                                                    bondLevelState98N,
                                                    bondLevelState99N,
                                                    bondLevelState100N,
                                                    bondLevelState101N,
                                                    bondLevelState102N,
                                                    bondLevelState103N,
                                                    bondLevelState104N,
                                                    bondLevelState105N,
                                                    bondLevelState106N,
                                                    bondLevelState107N,
                                                    bondLevelState108N,
                                                    bondLevelState109N,
                                                    bondLevelState110N,
                                                    bondLevelState111N,
                                                    bondLevelState112N,
                                                    bondLevelState113N,
                                                    bondLevelState114N,
                                                    bondLevelState115N,
                                                    bondLevelState116N,
                                                    bondLevelState117N,
                                                    bondLevelState118N,
                                                    bondLevelState119N,
                                                    bondLevelState120N,
                                                    bondLevelState121N,
                                                    bondLevelState122N,
                                                    bondLevelState123N,
                                                    bondLevelState124N,
                                                    bondLevelState125N,
                                                    bondLevelState126N,
                                                    bondLevelState127N,
                                                    bondLevelState128N,
                                                    bondLevelState129N,
                                                    bondLevelState130N,
                                                    bondLevelState131N,
                                                    bondLevelState132N,
                                                    bondLevelState133N,
                                                    bondLevelState134N,
                                                    bondLevelState135N,
                                                    bondLevelState136N,
                                                    bondLevelState137N,
                                                    bondLevelState138N,
                                                    bondLevelState139N,
                                                    bondLevelState140N,
                                                    bondLevelState141N,
                                                    bondLevelState142N,
                                                    bondLevelState143N,
                                                    bondLevelState144N,
                                                    bondLevelState145N,
                                                    bondLevelState146N,
                                                    bondLevelState147N,
                                                    bondLevelState148N,
                                                    bondLevelState149N,
                                                    bondLevelState150N,
                                                    bondLevelState151N,
                                                    bondLevelState152N,
                                                    bondLevelState153N,
                                                    bondLevelState154N,
                                                    bondLevelState155N,
                                                    bondLevelState156N,
                                                    bondLevelState157N,
                                                    bondLevelState158N,
                                                    bondLevelState159N,
                                                    bondLevelState160N,
                                                    bondLevelState161N,
                                                    bondLevelState162N,
                                                    bondLevelState163N,
                                                    bondLevelState164N,
                                                    bondLevelState165N,
                                                    bondLevelState166N,
                                                    bondLevelState167N,
                                                    bondLevelState168N,
                                                    bondLevelState169N,
                                                    bondLevelState170N,
                                                    bondLevelState171N,
                                                    bondLevelState172N,
                                                    bondLevelState173N,
                                                    bondLevelState174N,
                                                    bondLevelState175N,
                                                    bondLevelState176N,
                                                    bondLevelState177N,
                                                    bondLevelState178N,
                                                    bondLevelState179N,
                                                    bondLevelState180N,
                                                    bondLevelState181N,
                                                    bondLevelState182N,
                                                    bondLevelState183N,
                                                    bondLevelState184N,
                                                    bondLevelState185N,
                                                    bondLevelState186N,
                                                    bondLevelState187N,
                                                    bondLevelState188N,
                                                    bondLevelState189N,
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
                                                    bondLevelState50NP1,
                                                    bondLevelState51NP1,
                                                    bondLevelState52NP1,
                                                    bondLevelState53NP1,
                                                    bondLevelState54NP1,
                                                    bondLevelState55NP1,
                                                    bondLevelState56NP1,
                                                    bondLevelState57NP1,
                                                    bondLevelState58NP1,
                                                    bondLevelState59NP1,
                                                    bondLevelState60NP1,
                                                    bondLevelState61NP1,
                                                    bondLevelState62NP1,
                                                    bondLevelState63NP1,
                                                    bondLevelState64NP1,
                                                    bondLevelState65NP1,
                                                    bondLevelState66NP1,
                                                    bondLevelState67NP1,
                                                    bondLevelState68NP1,
                                                    bondLevelState69NP1,
                                                    bondLevelState70NP1,
                                                    bondLevelState71NP1,
                                                    bondLevelState72NP1,
                                                    bondLevelState73NP1,
                                                    bondLevelState74NP1,
                                                    bondLevelState75NP1,
                                                    bondLevelState76NP1,
                                                    bondLevelState77NP1,
                                                    bondLevelState78NP1,
                                                    bondLevelState79NP1,
                                                    bondLevelState80NP1,
                                                    bondLevelState81NP1,
                                                    bondLevelState82NP1,
                                                    bondLevelState83NP1,
                                                    bondLevelState84NP1,
                                                    bondLevelState85NP1,
                                                    bondLevelState86NP1,
                                                    bondLevelState87NP1,
                                                    bondLevelState88NP1,
                                                    bondLevelState89NP1,
                                                    bondLevelState90NP1,
                                                    bondLevelState91NP1,
                                                    bondLevelState92NP1,
                                                    bondLevelState93NP1,
                                                    bondLevelState94NP1,
                                                    bondLevelState95NP1,
                                                    bondLevelState96NP1,
                                                    bondLevelState97NP1,
                                                    bondLevelState98NP1,
                                                    bondLevelState99NP1,
                                                    bondLevelState100NP1,
                                                    bondLevelState101NP1,
                                                    bondLevelState102NP1,
                                                    bondLevelState103NP1,
                                                    bondLevelState104NP1,
                                                    bondLevelState105NP1,
                                                    bondLevelState106NP1,
                                                    bondLevelState107NP1,
                                                    bondLevelState108NP1,
                                                    bondLevelState109NP1,
                                                    bondLevelState110NP1,
                                                    bondLevelState111NP1,
                                                    bondLevelState112NP1,
                                                    bondLevelState113NP1,
                                                    bondLevelState114NP1,
                                                    bondLevelState115NP1,
                                                    bondLevelState116NP1,
                                                    bondLevelState117NP1,
                                                    bondLevelState118NP1,
                                                    bondLevelState119NP1,
                                                    bondLevelState120NP1,
                                                    bondLevelState121NP1,
                                                    bondLevelState122NP1,
                                                    bondLevelState123NP1,
                                                    bondLevelState124NP1,
                                                    bondLevelState125NP1,
                                                    bondLevelState126NP1,
                                                    bondLevelState127NP1,
                                                    bondLevelState128NP1,
                                                    bondLevelState129NP1,
                                                    bondLevelState130NP1,
                                                    bondLevelState131NP1,
                                                    bondLevelState132NP1,
                                                    bondLevelState133NP1,
                                                    bondLevelState134NP1,
                                                    bondLevelState135NP1,
                                                    bondLevelState136NP1,
                                                    bondLevelState137NP1,
                                                    bondLevelState138NP1,
                                                    bondLevelState139NP1,
                                                    bondLevelState140NP1,
                                                    bondLevelState141NP1,
                                                    bondLevelState142NP1,
                                                    bondLevelState143NP1,
                                                    bondLevelState144NP1,
                                                    bondLevelState145NP1,
                                                    bondLevelState146NP1,
                                                    bondLevelState147NP1,
                                                    bondLevelState148NP1,
                                                    bondLevelState149NP1,
                                                    bondLevelState150NP1,
                                                    bondLevelState151NP1,
                                                    bondLevelState152NP1,
                                                    bondLevelState153NP1,
                                                    bondLevelState154NP1,
                                                    bondLevelState155NP1,
                                                    bondLevelState156NP1,
                                                    bondLevelState157NP1,
                                                    bondLevelState158NP1,
                                                    bondLevelState159NP1,
                                                    bondLevelState160NP1,
                                                    bondLevelState161NP1,
                                                    bondLevelState162NP1,
                                                    bondLevelState163NP1,
                                                    bondLevelState164NP1,
                                                    bondLevelState165NP1,
                                                    bondLevelState166NP1,
                                                    bondLevelState167NP1,
                                                    bondLevelState168NP1,
                                                    bondLevelState169NP1,
                                                    bondLevelState170NP1,
                                                    bondLevelState171NP1,
                                                    bondLevelState172NP1,
                                                    bondLevelState173NP1,
                                                    bondLevelState174NP1,
                                                    bondLevelState175NP1,
                                                    bondLevelState176NP1,
                                                    bondLevelState177NP1,
                                                    bondLevelState178NP1,
                                                    bondLevelState179NP1,
                                                    bondLevelState180NP1,
                                                    bondLevelState181NP1,
                                                    bondLevelState182NP1,
                                                    bondLevelState183NP1,
                                                    bondLevelState184NP1,
                                                    bondLevelState185NP1,
                                                    bondLevelState186NP1,
                                                    bondLevelState187NP1,
                                                    bondLevelState188NP1,
                                                    bondLevelState189NP1,
                                                    bondLevelInternalEnergyNP1,
                                                    bondLevelInelasticEnergyNP1,
                                                    neighborhoodList,
                                                    numOwnedPoints,
                                                    dt);
}
