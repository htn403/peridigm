/*! \file Peridigm_ElasticBondAssociatedCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticBondAssociatedCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticBondAssociatedCorrespondenceMaterial::ElasticBondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : BondAssociatedCorrespondenceMaterial(params),
    m_vonMisesStressFieldId(-1), 
    m_bondLevelVonMisesStressFieldId(-1)
{

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_vonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Von_Mises_Stress");

  m_fieldIds.push_back(m_vonMisesStressFieldId);

  m_bondLevelVonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Von_Mises_Stress");

  m_fieldIds.push_back(m_bondLevelVonMisesStressFieldId);
}

PeridigmNS::ElasticBondAssociatedCorrespondenceMaterial::~ElasticBondAssociatedCorrespondenceMaterial()
{
}

void
PeridigmNS::ElasticBondAssociatedCorrespondenceMaterial::initialize(const double dt,
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

  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
}

void
PeridigmNS::ElasticBondAssociatedCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                                             const int numOwnedPoints,
                                                                             const int* neighborhoodList,
                                                                             PeridigmNS::DataManager& dataManager) const
{

  // Compute the node-level stress values
  // This is only done for output (visualization) purposes
  double *unrotatedCauchyStressN;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);

  double *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);

  double *unrotatedRateOfDeformation;
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  double *vonMisesStress;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);

  CORRESPONDENCE::updateElasticCauchyStress(unrotatedRateOfDeformation, 
                                            unrotatedCauchyStressN, 
                                            unrotatedCauchyStressNP1,
                                            vonMisesStress,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_shearModulus,
                                            dt);

  // Compute the bond-level stress values
  double *bondLevelUnrotatedCauchyStressXXN, *bondLevelUnrotatedCauchyStressXYN, *bondLevelUnrotatedCauchyStressXZN;
  double *bondLevelUnrotatedCauchyStressYXN, *bondLevelUnrotatedCauchyStressYYN, *bondLevelUnrotatedCauchyStressYZN;
  double *bondLevelUnrotatedCauchyStressZXN, *bondLevelUnrotatedCauchyStressZYN, *bondLevelUnrotatedCauchyStressZZN;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXZN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYZN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZZN);

  double *bondLevelUnrotatedCauchyStressXXNP1, *bondLevelUnrotatedCauchyStressXYNP1, *bondLevelUnrotatedCauchyStressXZNP1;
  double *bondLevelUnrotatedCauchyStressYXNP1, *bondLevelUnrotatedCauchyStressYYNP1, *bondLevelUnrotatedCauchyStressYZNP1;
  double *bondLevelUnrotatedCauchyStressZXNP1, *bondLevelUnrotatedCauchyStressZYNP1, *bondLevelUnrotatedCauchyStressZZNP1;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZZNP1);

  double *bondLevelUnrotatedRateOfDeformationXX, *bondLevelUnrotatedRateOfDeformationXY, *bondLevelUnrotatedRateOfDeformationXZ;
  double *bondLevelUnrotatedRateOfDeformationYX, *bondLevelUnrotatedRateOfDeformationYY, *bondLevelUnrotatedRateOfDeformationYZ;
  double *bondLevelUnrotatedRateOfDeformationZX, *bondLevelUnrotatedRateOfDeformationZY, *bondLevelUnrotatedRateOfDeformationZZ;
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZZ);

  double *bondLevelVonMisesStress;
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVonMisesStress);

  CORRESPONDENCE::updateBondLevelElasticCauchyStress(bondLevelUnrotatedRateOfDeformationXX,
                                                     bondLevelUnrotatedRateOfDeformationXY,
                                                     bondLevelUnrotatedRateOfDeformationXZ,
                                                     bondLevelUnrotatedRateOfDeformationYX,
                                                     bondLevelUnrotatedRateOfDeformationYY,
                                                     bondLevelUnrotatedRateOfDeformationYZ,
                                                     bondLevelUnrotatedRateOfDeformationZX,
                                                     bondLevelUnrotatedRateOfDeformationZY,
                                                     bondLevelUnrotatedRateOfDeformationZZ,
                                                     bondLevelUnrotatedCauchyStressXXN, 
                                                     bondLevelUnrotatedCauchyStressXYN, 
                                                     bondLevelUnrotatedCauchyStressXZN, 
                                                     bondLevelUnrotatedCauchyStressYXN, 
                                                     bondLevelUnrotatedCauchyStressYYN, 
                                                     bondLevelUnrotatedCauchyStressYZN, 
                                                     bondLevelUnrotatedCauchyStressZXN, 
                                                     bondLevelUnrotatedCauchyStressZYN, 
                                                     bondLevelUnrotatedCauchyStressZZN, 
                                                     bondLevelUnrotatedCauchyStressXXNP1, 
                                                     bondLevelUnrotatedCauchyStressXYNP1, 
                                                     bondLevelUnrotatedCauchyStressXZNP1, 
                                                     bondLevelUnrotatedCauchyStressYXNP1, 
                                                     bondLevelUnrotatedCauchyStressYYNP1, 
                                                     bondLevelUnrotatedCauchyStressYZNP1, 
                                                     bondLevelUnrotatedCauchyStressZXNP1, 
                                                     bondLevelUnrotatedCauchyStressZYNP1, 
                                                     bondLevelUnrotatedCauchyStressZZNP1, 
                                                     bondLevelVonMisesStress,
                                                     neighborhoodList,
                                                     numOwnedPoints,
                                                     m_bulkModulus,
                                                     m_shearModulus,
                                                     dt);
}

void
PeridigmNS::ElasticBondAssociatedCorrespondenceMaterial::computePK2Stress(const double dt,
                                                                          const int numOwnedPoints,
                                                                          const int* neighborhoodList,
                                                                          PeridigmNS::DataManager& dataManager) const
{

  // Compute the node-level stress values
  // This is only done for output (visualization) purposes
  double *PK2StressN;
  dataManager.getData(m_PK2StressFieldId, PeridigmField::STEP_N)->ExtractView(&PK2StressN);

  double *PK2StressNP1;
  dataManager.getData(m_PK2StressFieldId, PeridigmField::STEP_NP1)->ExtractView(&PK2StressNP1);

  double *strainRate;
  dataManager.getData(m_strainRateFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainRate);

  double *vonMisesStress;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);

  CORRESPONDENCE::updateElasticCauchyStress(strainRate, 
                                            PK2StressN, 
                                            PK2StressNP1,
                                            vonMisesStress,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_shearModulus,
                                            dt);

  // Compute the bond-level stress values
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

  double *bondLevelVonMisesStress;
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVonMisesStress);

  CORRESPONDENCE::updateBondLevelElasticCauchyStress(bondLevelStrainRateXX,
                                                     bondLevelStrainRateXY,
                                                     bondLevelStrainRateXZ,
                                                     bondLevelStrainRateYX,
                                                     bondLevelStrainRateYY,
                                                     bondLevelStrainRateYZ,
                                                     bondLevelStrainRateZX,
                                                     bondLevelStrainRateZY,
                                                     bondLevelStrainRateZZ,
                                                     bondLevelPK2StressXXN, 
                                                     bondLevelPK2StressXYN, 
                                                     bondLevelPK2StressXZN, 
                                                     bondLevelPK2StressYXN, 
                                                     bondLevelPK2StressYYN, 
                                                     bondLevelPK2StressYZN, 
                                                     bondLevelPK2StressZXN, 
                                                     bondLevelPK2StressZYN, 
                                                     bondLevelPK2StressZZN, 
                                                     bondLevelPK2StressXXNP1, 
                                                     bondLevelPK2StressXYNP1, 
                                                     bondLevelPK2StressXZNP1, 
                                                     bondLevelPK2StressYXNP1, 
                                                     bondLevelPK2StressYYNP1, 
                                                     bondLevelPK2StressYZNP1, 
                                                     bondLevelPK2StressZXNP1, 
                                                     bondLevelPK2StressZYNP1, 
                                                     bondLevelPK2StressZZNP1, 
                                                     bondLevelVonMisesStress,
                                                     neighborhoodList,
                                                     numOwnedPoints,
                                                     m_bulkModulus,
                                                     m_shearModulus,
                                                     dt);
}