//! \file Peridigm_MicroplaneBondAssociatedCorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_MICROPLANEBONDASSOCIATEDCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_MICROPLANEBONDASSOCIATEDCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_BondAssociatedCorrespondenceMaterial.hpp"

namespace PeridigmNS {

  class MicroplaneBondAssociatedCorrespondenceMaterial : public BondAssociatedCorrespondenceMaterial{
  public:

    //! Constructor.
    MicroplaneBondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~MicroplaneBondAssociatedCorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Microplane Bond Associated Correspondence"); }

    //! Initialize the derived class
    virtual void initialize(const double dt, 
                            const int numOwnedPoints, 
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager);

    //! Evaluate the Cauchy stress.
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     const int* neighborhoodList,
                                     PeridigmNS::DataManager& dataManager) const;

    //! Evaluate the second Piola-Kirchhoff stress
    virtual void computePK2Stress(const double dt,
                                  const int numOwnedPoints,
                                  const int* neighborhoodList,
                                  PeridigmNS::DataManager& dataManager) const;

    //! Returns the requested material property
    //! A dummy method here.
    virtual double lookupMaterialProperty(const std::string keyname) const {return 0.0;}


  protected:

    double m_youngsModulus;
    double m_poissonsRatio;
    double m_k1;
    double m_k2;
    double m_k3;
    double m_k4;

    // field spec ids for all relevant data
    int m_microState1FieldId;
    int m_microState2FieldId;
    int m_microState3FieldId;
    int m_microState4FieldId;
    int m_microState5FieldId;
    int m_microState6FieldId;
    int m_microState7FieldId;
    int m_microState8FieldId;
    int m_microState9FieldId;
    int m_microState10FieldId;
    int m_microState11FieldId;
    int m_microState12FieldId;
    int m_microState13FieldId;
    int m_microState14FieldId;
    int m_microState15FieldId;
    int m_microState16FieldId;
    int m_microState17FieldId;
    int m_microState18FieldId;
    int m_microState19FieldId;
    int m_microState20FieldId;
    int m_microState21FieldId;
    int m_microState22FieldId;
    int m_microState23FieldId;
    int m_microState24FieldId;
    int m_microState25FieldId;
    int m_microState26FieldId;
    int m_microState27FieldId;
    int m_microState28FieldId;
    int m_microState29FieldId;
    int m_microState30FieldId;
    int m_microState31FieldId;
    int m_microState32FieldId;
    int m_microState33FieldId;
    int m_microState34FieldId;
    int m_microState35FieldId;
    int m_microState36FieldId;
    int m_microState37FieldId;
    int m_microState38FieldId;
    int m_microState39FieldId;
    int m_microState40FieldId;
    int m_microState41FieldId;
    int m_microState42FieldId;
    int m_microState43FieldId;
    int m_microState44FieldId;
    int m_microState45FieldId;
    int m_microState46FieldId;
    int m_microState47FieldId;
    int m_microState48FieldId;
    int m_microState49FieldId;
    int m_microState50FieldId;
    int m_microState51FieldId;
    int m_microState52FieldId;
    int m_microState53FieldId;
    int m_microState54FieldId;
    int m_microState55FieldId;
    int m_microState56FieldId;
    int m_microState57FieldId;
    int m_microState58FieldId;
    int m_microState59FieldId;
    int m_microState60FieldId;
    int m_microState61FieldId;
    int m_microState62FieldId;
    int m_microState63FieldId;
    int m_microState64FieldId;
    int m_microState65FieldId;
    int m_microState66FieldId;
    int m_microState67FieldId;
    int m_microState68FieldId;
    int m_microState69FieldId;
    int m_microState70FieldId;
    int m_microState71FieldId;
    int m_microState72FieldId;
    int m_microState73FieldId;
    int m_microState74FieldId;
    int m_microState75FieldId;
    int m_microState76FieldId;
    int m_microState77FieldId;
    int m_microState78FieldId;
    int m_microState79FieldId;
    int m_microState80FieldId;
    int m_microState81FieldId;
    int m_microState82FieldId;
    int m_microState83FieldId;
    int m_microState84FieldId;
    int m_microState85FieldId;
    int m_microState86FieldId;
    int m_microState87FieldId;
    int m_microState88FieldId;
    int m_microState89FieldId;
    int m_microState90FieldId;
    int m_microState91FieldId;
    int m_microState92FieldId;
    int m_microState93FieldId;
    int m_microState94FieldId;
    int m_microState95FieldId;
    int m_microState96FieldId;
    int m_microState97FieldId;
    int m_microState98FieldId;
    int m_microState99FieldId;
    int m_microState100FieldId;
    int m_microState101FieldId;
    int m_microState102FieldId;
    int m_microState103FieldId;
    int m_microState104FieldId;
    int m_microState105FieldId;
    int m_microState106FieldId;
    int m_microState107FieldId;
    int m_microState108FieldId;
    int m_microState109FieldId;
    int m_microState110FieldId;
    int m_microState111FieldId;
    int m_microState112FieldId;
    int m_microState113FieldId;
    int m_microState114FieldId;
    int m_microState115FieldId;
    int m_microState116FieldId;
    int m_microState117FieldId;
    int m_microState118FieldId;
    int m_microState119FieldId;
    int m_microState120FieldId;
    int m_microState121FieldId;
    int m_microState122FieldId;
    int m_microState123FieldId;
    int m_microState124FieldId;
    int m_microState125FieldId;
    int m_microState126FieldId;
    int m_microState127FieldId;
    int m_microState128FieldId;
    int m_microState129FieldId;
    int m_microState130FieldId;
    int m_microState131FieldId;
    int m_microState132FieldId;
    int m_microState133FieldId;
    int m_microState134FieldId;
    int m_microState135FieldId;
    int m_microState136FieldId;
    int m_microState137FieldId;
    int m_microState138FieldId;
    int m_microState139FieldId;
    int m_microState140FieldId;
    int m_microState141FieldId;
    int m_microState142FieldId;
    int m_microState143FieldId;
    int m_microState144FieldId;
    int m_microState145FieldId;
    int m_microState146FieldId;
    int m_microState147FieldId;
    int m_microState148FieldId;
    int m_microState149FieldId;
    int m_microState150FieldId;
    int m_microState151FieldId;
    int m_microState152FieldId;
    int m_microState153FieldId;
    int m_microState154FieldId;
    int m_microState155FieldId;
    int m_microState156FieldId;
    int m_microState157FieldId;
    int m_microState158FieldId;
    int m_microState159FieldId;
    int m_microState160FieldId;
    int m_microState161FieldId;
    int m_microState162FieldId;
    int m_microState163FieldId;
    int m_microState164FieldId;
    int m_microState165FieldId;
    int m_microState166FieldId;
    int m_microState167FieldId;
    int m_microState168FieldId;
    int m_microState169FieldId;
    int m_microState170FieldId;
    int m_microState171FieldId;
    int m_microState172FieldId;
    int m_microState173FieldId;
    int m_microState174FieldId;
    int m_microState175FieldId;
    int m_microState176FieldId;
    int m_microState177FieldId;
    int m_microState178FieldId;
    int m_microState179FieldId;
    int m_microState180FieldId;
    int m_microState181FieldId;
    int m_microState182FieldId;
    int m_microState183FieldId;
    int m_microState184FieldId;
    int m_microState185FieldId;
    int m_microState186FieldId;
    int m_microState187FieldId;
    int m_microState188FieldId;
    int m_microState189FieldId;
    int m_internalEnergyFieldId;
    int m_inelasticEnergyFieldId;
    
    int m_bondLevelMicroState1FieldId;
    int m_bondLevelMicroState2FieldId;
    int m_bondLevelMicroState3FieldId;
    int m_bondLevelMicroState4FieldId;
    int m_bondLevelMicroState5FieldId;
    int m_bondLevelMicroState6FieldId;
    int m_bondLevelMicroState7FieldId;
    int m_bondLevelMicroState8FieldId;
    int m_bondLevelMicroState9FieldId;
    int m_bondLevelMicroState10FieldId;
    int m_bondLevelMicroState11FieldId;
    int m_bondLevelMicroState12FieldId;
    int m_bondLevelMicroState13FieldId;
    int m_bondLevelMicroState14FieldId;
    int m_bondLevelMicroState15FieldId;
    int m_bondLevelMicroState16FieldId;
    int m_bondLevelMicroState17FieldId;
    int m_bondLevelMicroState18FieldId;
    int m_bondLevelMicroState19FieldId;
    int m_bondLevelMicroState20FieldId;
    int m_bondLevelMicroState21FieldId;
    int m_bondLevelMicroState22FieldId;
    int m_bondLevelMicroState23FieldId;
    int m_bondLevelMicroState24FieldId;
    int m_bondLevelMicroState25FieldId;
    int m_bondLevelMicroState26FieldId;
    int m_bondLevelMicroState27FieldId;
    int m_bondLevelMicroState28FieldId;
    int m_bondLevelMicroState29FieldId;
    int m_bondLevelMicroState30FieldId;
    int m_bondLevelMicroState31FieldId;
    int m_bondLevelMicroState32FieldId;
    int m_bondLevelMicroState33FieldId;
    int m_bondLevelMicroState34FieldId;
    int m_bondLevelMicroState35FieldId;
    int m_bondLevelMicroState36FieldId;
    int m_bondLevelMicroState37FieldId;
    int m_bondLevelMicroState38FieldId;
    int m_bondLevelMicroState39FieldId;
    int m_bondLevelMicroState40FieldId;
    int m_bondLevelMicroState41FieldId;
    int m_bondLevelMicroState42FieldId;
    int m_bondLevelMicroState43FieldId;
    int m_bondLevelMicroState44FieldId;
    int m_bondLevelMicroState45FieldId;
    int m_bondLevelMicroState46FieldId;
    int m_bondLevelMicroState47FieldId;
    int m_bondLevelMicroState48FieldId;
    int m_bondLevelMicroState49FieldId;
    int m_bondLevelMicroState50FieldId;
    int m_bondLevelMicroState51FieldId;
    int m_bondLevelMicroState52FieldId;
    int m_bondLevelMicroState53FieldId;
    int m_bondLevelMicroState54FieldId;
    int m_bondLevelMicroState55FieldId;
    int m_bondLevelMicroState56FieldId;
    int m_bondLevelMicroState57FieldId;
    int m_bondLevelMicroState58FieldId;
    int m_bondLevelMicroState59FieldId;
    int m_bondLevelMicroState60FieldId;
    int m_bondLevelMicroState61FieldId;
    int m_bondLevelMicroState62FieldId;
    int m_bondLevelMicroState63FieldId;
    int m_bondLevelMicroState64FieldId;
    int m_bondLevelMicroState65FieldId;
    int m_bondLevelMicroState66FieldId;
    int m_bondLevelMicroState67FieldId;
    int m_bondLevelMicroState68FieldId;
    int m_bondLevelMicroState69FieldId;
    int m_bondLevelMicroState70FieldId;
    int m_bondLevelMicroState71FieldId;
    int m_bondLevelMicroState72FieldId;
    int m_bondLevelMicroState73FieldId;
    int m_bondLevelMicroState74FieldId;
    int m_bondLevelMicroState75FieldId;
    int m_bondLevelMicroState76FieldId;
    int m_bondLevelMicroState77FieldId;
    int m_bondLevelMicroState78FieldId;
    int m_bondLevelMicroState79FieldId;
    int m_bondLevelMicroState80FieldId;
    int m_bondLevelMicroState81FieldId;
    int m_bondLevelMicroState82FieldId;
    int m_bondLevelMicroState83FieldId;
    int m_bondLevelMicroState84FieldId;
    int m_bondLevelMicroState85FieldId;
    int m_bondLevelMicroState86FieldId;
    int m_bondLevelMicroState87FieldId;
    int m_bondLevelMicroState88FieldId;
    int m_bondLevelMicroState89FieldId;
    int m_bondLevelMicroState90FieldId;
    int m_bondLevelMicroState91FieldId;
    int m_bondLevelMicroState92FieldId;
    int m_bondLevelMicroState93FieldId;
    int m_bondLevelMicroState94FieldId;
    int m_bondLevelMicroState95FieldId;
    int m_bondLevelMicroState96FieldId;
    int m_bondLevelMicroState97FieldId;
    int m_bondLevelMicroState98FieldId;
    int m_bondLevelMicroState99FieldId;
    int m_bondLevelMicroState100FieldId;
    int m_bondLevelMicroState101FieldId;
    int m_bondLevelMicroState102FieldId;
    int m_bondLevelMicroState103FieldId;
    int m_bondLevelMicroState104FieldId;
    int m_bondLevelMicroState105FieldId;
    int m_bondLevelMicroState106FieldId;
    int m_bondLevelMicroState107FieldId;
    int m_bondLevelMicroState108FieldId;
    int m_bondLevelMicroState109FieldId;
    int m_bondLevelMicroState110FieldId;
    int m_bondLevelMicroState111FieldId;
    int m_bondLevelMicroState112FieldId;
    int m_bondLevelMicroState113FieldId;
    int m_bondLevelMicroState114FieldId;
    int m_bondLevelMicroState115FieldId;
    int m_bondLevelMicroState116FieldId;
    int m_bondLevelMicroState117FieldId;
    int m_bondLevelMicroState118FieldId;
    int m_bondLevelMicroState119FieldId;
    int m_bondLevelMicroState120FieldId;
    int m_bondLevelMicroState121FieldId;
    int m_bondLevelMicroState122FieldId;
    int m_bondLevelMicroState123FieldId;
    int m_bondLevelMicroState124FieldId;
    int m_bondLevelMicroState125FieldId;
    int m_bondLevelMicroState126FieldId;
    int m_bondLevelMicroState127FieldId;
    int m_bondLevelMicroState128FieldId;
    int m_bondLevelMicroState129FieldId;
    int m_bondLevelMicroState130FieldId;
    int m_bondLevelMicroState131FieldId;
    int m_bondLevelMicroState132FieldId;
    int m_bondLevelMicroState133FieldId;
    int m_bondLevelMicroState134FieldId;
    int m_bondLevelMicroState135FieldId;
    int m_bondLevelMicroState136FieldId;
    int m_bondLevelMicroState137FieldId;
    int m_bondLevelMicroState138FieldId;
    int m_bondLevelMicroState139FieldId;
    int m_bondLevelMicroState140FieldId;
    int m_bondLevelMicroState141FieldId;
    int m_bondLevelMicroState142FieldId;
    int m_bondLevelMicroState143FieldId;
    int m_bondLevelMicroState144FieldId;
    int m_bondLevelMicroState145FieldId;
    int m_bondLevelMicroState146FieldId;
    int m_bondLevelMicroState147FieldId;
    int m_bondLevelMicroState148FieldId;
    int m_bondLevelMicroState149FieldId;
    int m_bondLevelMicroState150FieldId;
    int m_bondLevelMicroState151FieldId;
    int m_bondLevelMicroState152FieldId;
    int m_bondLevelMicroState153FieldId;
    int m_bondLevelMicroState154FieldId;
    int m_bondLevelMicroState155FieldId;
    int m_bondLevelMicroState156FieldId;
    int m_bondLevelMicroState157FieldId;
    int m_bondLevelMicroState158FieldId;
    int m_bondLevelMicroState159FieldId;
    int m_bondLevelMicroState160FieldId;
    int m_bondLevelMicroState161FieldId;
    int m_bondLevelMicroState162FieldId;
    int m_bondLevelMicroState163FieldId;
    int m_bondLevelMicroState164FieldId;
    int m_bondLevelMicroState165FieldId;
    int m_bondLevelMicroState166FieldId;
    int m_bondLevelMicroState167FieldId;
    int m_bondLevelMicroState168FieldId;
    int m_bondLevelMicroState169FieldId;
    int m_bondLevelMicroState170FieldId;
    int m_bondLevelMicroState171FieldId;
    int m_bondLevelMicroState172FieldId;
    int m_bondLevelMicroState173FieldId;
    int m_bondLevelMicroState174FieldId;
    int m_bondLevelMicroState175FieldId;
    int m_bondLevelMicroState176FieldId;
    int m_bondLevelMicroState177FieldId;
    int m_bondLevelMicroState178FieldId;
    int m_bondLevelMicroState179FieldId;
    int m_bondLevelMicroState180FieldId;
    int m_bondLevelMicroState181FieldId;
    int m_bondLevelMicroState182FieldId;
    int m_bondLevelMicroState183FieldId;
    int m_bondLevelMicroState184FieldId;
    int m_bondLevelMicroState185FieldId;
    int m_bondLevelMicroState186FieldId;
    int m_bondLevelMicroState187FieldId;
    int m_bondLevelMicroState188FieldId;
    int m_bondLevelMicroState189FieldId;
    int m_bondLevelInternalEnergyFieldId;
    int m_bondLevelInelasticEnergyFieldId;
  };
}

#endif // PERIDIGM_MICROPLANEBONDASSOCIATEDCORRESPONDENCEMATERIAL_HPP
