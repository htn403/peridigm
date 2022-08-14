//! \file Peridigm_CDPM2CorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_CDPM2CORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_CDPM2CORRESPONDENCEMATERIAL_HPP

#include "Peridigm_CorrespondenceMaterial.hpp"

namespace PeridigmNS {

  class CDPM2CorrespondenceMaterial : public CorrespondenceMaterial{
  public:

    //! Constructor.
    CDPM2CorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~CDPM2CorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("CDPM2 Correspondence"); }

    //! Initialize the derived class
    virtual void initialize(const double dt, 
                            const int numOwnedPoints, 
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager);

    //! Evaluate the Cauchy stress.
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     PeridigmNS::DataManager& dataManager) const;


    //! Returns the requested material property
    //! A dummy method here.
    virtual double lookupMaterialProperty(const std::string keyname) const {return 0.0;}


    protected:
  
    double m_unitFlag;
    double m_youngsModulus;
    double m_poissonsRatio;
    double m_fc;
    double m_ft;
    double m_fracEnergy;
    double m_hard1;
    double m_hard2;
    double m_hard3;
    double m_hard4;
    double m_hardqh2;
    double m_soften;
    double m_length;

    // field spec ids for all relevant data
    int m_unrotatedRateOfDeformationFieldId;
    int m_unrotatedCauchyStressFieldId; 
    int m_leftStretchTensorFieldId; 
    int m_rotationTensorFieldId;
    int m_tempDeformationGradientFieldId;
    int m_GreenLagrangeStrainFieldId; 
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
    int m_internalEnergyFieldId;
    int m_inelasticEnergyFieldId;
    
  };
}

#endif // PERIDIGM_CDPM2CORRESPONDENCEMATERIAL_HPP
