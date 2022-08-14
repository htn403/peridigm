//! \file cdpm2_correspondence.cxx

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

#include "cdpm2_correspondence.h"
#include "correspondence.h"
//#include "material_utilities.h" // to use element volume
#include <Sacado.hpp>
#include <math.h> // for sqrt

extern"C" {
void cdpm2inputparams_( const double* uFlag, 
                        const double* E, 
                        const double* nu,
                        const double* comStrength, 
                        const double* tenStrength, 
                        const double* fracEnergy, 
                        const double* hard1,
                        const double* hard2,
                        const double* hard3,
                        const double* hard4,
                        const double* hardqh2,
                        const double* soften,
                        const double* length);
}


extern"C" {
void cdpm2material_(const double *dt,
                    double *strainInc,
                    double *oldStrain,
                    double *stressOld, 
                    double *stateOld, 
                    const double *enerInternOld, 
                    const double *enerInelasOld,
                    double *stressNew, 
                    double *stateNew, 
                    double *enerInternNew, 
                    double *enerInelasNew);
}

namespace CORRESPONDENCE {

void initializeCDPM2Model
(
    const double uFlag, 
    const double E, 
    const double nu,
    const double comStrength, 
    const double tenStrength, 
    const double fracEnergy, 
    const double hard1,
    const double hard2,
    const double hard3,
    const double hard4,
    const double hardqh2,
    const double soften,
    const double length
)
{
  cdpm2inputparams_(&uFlag, 
                    &E, 
                    &nu, 
                    &comStrength, 
                    &tenStrength, 
                    &fracEnergy,
                    &hard1, 
                    &hard2, 
                    &hard3, 
                    &hard4, 
                    &hardqh2,
                    &soften,
                    &length);
}


template<typename ScalarT>
void updateCDPM2Stress
(   
    const ScalarT* strainRate,
    const ScalarT* strain, 
    const ScalarT* stressN, 
    const ScalarT* state1N,
    const ScalarT* state2N,
    const ScalarT* state3N,
    const ScalarT* state4N,
    const ScalarT* state5N,
    const ScalarT* state6N,
    const ScalarT* state7N,
    const ScalarT* state8N,
    const ScalarT* state9N,
    const ScalarT* state10N,
    const ScalarT* state11N,
    const ScalarT* state12N,
    const ScalarT* state13N,
    const ScalarT* state14N,
    const ScalarT* state15N,
    const ScalarT* state16N,
    const ScalarT* state17N,
    const ScalarT* state18N,
    const ScalarT* state19N,
    const ScalarT* state20N,
    const ScalarT* state21N,
    const ScalarT* state22N,
    const ScalarT* state23N,
    const ScalarT* state24N,
    const ScalarT* state25N,
    const ScalarT* state26N,
    const ScalarT* state27N,
    const ScalarT* state28N,
    const ScalarT* state29N,
    const ScalarT* state30N,
    const ScalarT* state31N,
    const ScalarT* state32N,
    const ScalarT* state33N,
    const ScalarT* state34N,
    const ScalarT* state35N,
    const ScalarT* state36N,
    const ScalarT* state37N,
    const ScalarT* state38N,
    const ScalarT* state39N,
    const ScalarT* state40N,
    const ScalarT* state41N,
    const ScalarT* state42N,
    const ScalarT* state43N,
    const ScalarT* state44N,
    const ScalarT* state45N,
    const ScalarT* state46N,
    const ScalarT* state47N,
    const ScalarT* state48N,
    const ScalarT* state49N,
    const ScalarT* internalEnergyN,
    const ScalarT* inelasticEnergyN,
    ScalarT* stressNP1, 
    ScalarT* state1NP1,
    ScalarT* state2NP1,
    ScalarT* state3NP1,
    ScalarT* state4NP1,
    ScalarT* state5NP1,
    ScalarT* state6NP1,
    ScalarT* state7NP1,
    ScalarT* state8NP1,
    ScalarT* state9NP1,
    ScalarT* state10NP1,
    ScalarT* state11NP1,
    ScalarT* state12NP1,
    ScalarT* state13NP1,
    ScalarT* state14NP1,
    ScalarT* state15NP1,
    ScalarT* state16NP1,
    ScalarT* state17NP1,
    ScalarT* state18NP1,
    ScalarT* state19NP1,
    ScalarT* state20NP1,
    ScalarT* state21NP1,
    ScalarT* state22NP1,
    ScalarT* state23NP1,
    ScalarT* state24NP1,
    ScalarT* state25NP1,
    ScalarT* state26NP1,
    ScalarT* state27NP1,
    ScalarT* state28NP1,
    ScalarT* state29NP1,
    ScalarT* state30NP1,
    ScalarT* state31NP1,
    ScalarT* state32NP1,
    ScalarT* state33NP1,
    ScalarT* state34NP1,
    ScalarT* state35NP1,
    ScalarT* state36NP1,
    ScalarT* state37NP1,
    ScalarT* state38NP1,
    ScalarT* state39NP1,
    ScalarT* state40NP1,
    ScalarT* state41NP1,
    ScalarT* state42NP1,
    ScalarT* state43NP1,
    ScalarT* state44NP1,
    ScalarT* state45NP1,
    ScalarT* state46NP1,
    ScalarT* state47NP1,
    ScalarT* state48NP1,
    ScalarT* state49NP1,
    ScalarT* internalEnergyNP1,
    ScalarT* inelasticEnergyNP1,
    const int numPoints, 
    const double dt
)
{
  // Hooke's law
  const ScalarT* Edot = strainRate;
  const ScalarT* eps = strain;

  const ScalarT* stressPtrN = stressN;
  const ScalarT* state1PtrN = state1N;
  const ScalarT* state2PtrN = state2N;
  const ScalarT* state3PtrN = state3N;
  const ScalarT* state4PtrN = state4N;
  const ScalarT* state5PtrN = state5N;
  const ScalarT* state6PtrN = state6N;
  const ScalarT* state7PtrN = state7N;
  const ScalarT* state8PtrN = state8N;
  const ScalarT* state9PtrN = state9N;
  const ScalarT* state10PtrN = state10N;
  const ScalarT* state11PtrN = state11N;
  const ScalarT* state12PtrN = state12N;
  const ScalarT* state13PtrN = state13N;
  const ScalarT* state14PtrN = state14N;
  const ScalarT* state15PtrN = state15N;
  const ScalarT* state16PtrN = state16N;
  const ScalarT* state17PtrN = state17N;
  const ScalarT* state18PtrN = state18N;
  const ScalarT* state19PtrN = state19N;
  const ScalarT* state20PtrN = state20N;
  const ScalarT* state21PtrN = state21N;
  const ScalarT* state22PtrN = state22N;
  const ScalarT* state23PtrN = state23N;
  const ScalarT* state24PtrN = state24N;
  const ScalarT* state25PtrN = state25N;
  const ScalarT* state26PtrN = state26N;
  const ScalarT* state27PtrN = state27N;
  const ScalarT* state28PtrN = state28N;
  const ScalarT* state29PtrN = state29N;
  const ScalarT* state30PtrN = state30N;
  const ScalarT* state31PtrN = state31N;
  const ScalarT* state32PtrN = state32N;
  const ScalarT* state33PtrN = state33N;
  const ScalarT* state34PtrN = state34N;
  const ScalarT* state35PtrN = state35N;
  const ScalarT* state36PtrN = state36N;
  const ScalarT* state37PtrN = state37N;
  const ScalarT* state38PtrN = state38N;
  const ScalarT* state39PtrN = state39N;
  const ScalarT* state40PtrN = state40N;
  const ScalarT* state41PtrN = state41N;
  const ScalarT* state42PtrN = state42N;
  const ScalarT* state43PtrN = state43N;
  const ScalarT* state44PtrN = state44N;
  const ScalarT* state45PtrN = state45N;
  const ScalarT* state46PtrN = state46N;
  const ScalarT* state47PtrN = state47N;
  const ScalarT* state48PtrN = state48N;
  const ScalarT* state49PtrN = state49N;
  const ScalarT* internEnergyN = internalEnergyN;
  const ScalarT* inelasEnergyN = inelasticEnergyN;

  ScalarT* stressPtrNP1 = stressNP1;
  ScalarT* state1PtrNP1 = state1NP1;
  ScalarT* state2PtrNP1 = state2NP1;
  ScalarT* state3PtrNP1 = state3NP1;
  ScalarT* state4PtrNP1 = state4NP1;
  ScalarT* state5PtrNP1 = state5NP1;
  ScalarT* state6PtrNP1 = state6NP1;
  ScalarT* state7PtrNP1 = state7NP1;
  ScalarT* state8PtrNP1 = state8NP1;
  ScalarT* state9PtrNP1 = state9NP1;
  ScalarT* state10PtrNP1 = state10NP1;
  ScalarT* state11PtrNP1 = state11NP1;
  ScalarT* state12PtrNP1 = state12NP1;
  ScalarT* state13PtrNP1 = state13NP1;
  ScalarT* state14PtrNP1 = state14NP1;
  ScalarT* state15PtrNP1 = state15NP1;
  ScalarT* state16PtrNP1 = state16NP1;
  ScalarT* state17PtrNP1 = state17NP1;
  ScalarT* state18PtrNP1 = state18NP1;
  ScalarT* state19PtrNP1 = state19NP1;
  ScalarT* state20PtrNP1 = state20NP1;
  ScalarT* state21PtrNP1 = state21NP1;
  ScalarT* state22PtrNP1 = state22NP1;
  ScalarT* state23PtrNP1 = state23NP1;
  ScalarT* state24PtrNP1 = state24NP1;
  ScalarT* state25PtrNP1 = state25NP1;
  ScalarT* state26PtrNP1 = state26NP1;
  ScalarT* state27PtrNP1 = state27NP1;
  ScalarT* state28PtrNP1 = state28NP1;
  ScalarT* state29PtrNP1 = state29NP1;
  ScalarT* state30PtrNP1 = state30NP1;
  ScalarT* state31PtrNP1 = state31NP1;
  ScalarT* state32PtrNP1 = state32NP1;
  ScalarT* state33PtrNP1 = state33NP1;
  ScalarT* state34PtrNP1 = state34NP1;
  ScalarT* state35PtrNP1 = state35NP1;
  ScalarT* state36PtrNP1 = state36NP1;
  ScalarT* state37PtrNP1 = state37NP1;
  ScalarT* state38PtrNP1 = state38NP1;
  ScalarT* state39PtrNP1 = state39NP1;
  ScalarT* state40PtrNP1 = state40NP1;
  ScalarT* state41PtrNP1 = state41NP1;
  ScalarT* state42PtrNP1 = state42NP1;
  ScalarT* state43PtrNP1 = state43NP1;
  ScalarT* state44PtrNP1 = state44NP1;
  ScalarT* state45PtrNP1 = state45NP1;
  ScalarT* state46PtrNP1 = state46NP1;
  ScalarT* state47PtrNP1 = state47NP1;
  ScalarT* state48PtrNP1 = state48NP1;
  ScalarT* state49PtrNP1 = state49NP1;
  ScalarT* internEnergyNP1 = internalEnergyNP1;
  ScalarT* inelasEnergyNP1 = inelasticEnergyNP1;

  ScalarT deps[6];
  ScalarT epsN[6];

  ScalarT sigmaN[6];
  ScalarT sigmaNP1[6];

  ScalarT stateN[49];
  ScalarT stateNP1[49];
  for(int iID=0 ; iID<numPoints ; ++iID, 
      Edot+=9, eps+=9, 
      stressPtrN+=9, ++internEnergyN, ++inelasEnergyN,
      ++state1PtrN, ++state2PtrN, ++state3PtrN, ++state4PtrN, ++state5PtrN,
      ++state6PtrN, ++state7PtrN, ++state8PtrN, ++state9PtrN, ++state10PtrN,
      ++state11PtrN, ++state12PtrN, ++state13PtrN, ++state14PtrN, ++state15PtrN,
      ++state16PtrN, ++state17PtrN, ++state18PtrN, ++state19PtrN, ++state20PtrN,
      ++state21PtrN, ++state22PtrN, ++state23PtrN, ++state24PtrN, ++state25PtrN,
      ++state26PtrN, ++state27PtrN, ++state28PtrN, ++state29PtrN, ++state30PtrN,
      ++state31PtrN, ++state32PtrN, ++state33PtrN, ++state34PtrN, ++state35PtrN,
      ++state36PtrN, ++state37PtrN, ++state38PtrN, ++state39PtrN, ++state40PtrN,
      ++state41PtrN, ++state42PtrN, ++state43PtrN, ++state44PtrN, ++state45PtrN,
      ++state46PtrN, ++state47PtrN, ++state48PtrN, ++state49PtrN, 
      stressPtrNP1+=9, ++internEnergyNP1, ++inelasEnergyNP1,
      ++state1PtrNP1, ++state2PtrNP1, ++state3PtrNP1, ++state4PtrNP1, ++state5PtrNP1,
      ++state6PtrNP1, ++state7PtrNP1, ++state8PtrNP1, ++state9PtrNP1, ++state10PtrNP1,
      ++state11PtrNP1, ++state12PtrNP1, ++state13PtrNP1, ++state14PtrNP1, ++state15PtrNP1,
      ++state16PtrNP1, ++state17PtrNP1, ++state18PtrNP1, ++state19PtrNP1, ++state20PtrNP1,
      ++state21PtrNP1, ++state22PtrNP1, ++state23PtrNP1, ++state24PtrNP1, ++state25PtrNP1,
      ++state26PtrNP1, ++state27PtrNP1, ++state28PtrNP1, ++state29PtrNP1, ++state30PtrNP1,
      ++state31PtrNP1, ++state32PtrNP1, ++state33PtrNP1, ++state34PtrNP1, ++state35PtrNP1,
      ++state36PtrNP1, ++state37PtrNP1, ++state38PtrNP1, ++state39PtrNP1, ++state40PtrNP1,
      ++state41PtrNP1, ++state42PtrNP1, ++state43PtrNP1, ++state44PtrNP1, ++state45PtrNP1,
      ++state46PtrNP1, ++state47PtrNP1, ++state48PtrNP1, ++state49PtrNP1
      ){

          *(stressPtrNP1+0) = sigmaNP1[0];
        *(stressPtrNP1+4) = sigmaNP1[1];
        *(stressPtrNP1+8) = sigmaNP1[2];
        *(stressPtrNP1+5) = sigmaNP1[3];
        *(stressPtrNP1+7) = sigmaNP1[3];
        *(stressPtrNP1+2) = sigmaNP1[4];
        *(stressPtrNP1+6) = sigmaNP1[4];
        *(stressPtrNP1+1) = sigmaNP1[5];
        *(stressPtrNP1+3) = sigmaNP1[5];
          
        // voigt notation
        epsN[0] = *(eps+0);
        epsN[1] = *(eps+4);
        epsN[2] = *(eps+8);
        epsN[3] = *(eps+5);
        epsN[4] = *(eps+2);
        epsN[5] = *(eps+1);

        // voigt notation
        deps[0] = *(Edot+0) * dt;
        deps[1] = *(Edot+4) * dt;
        deps[2] = *(Edot+8) * dt;
        deps[3] = 0.5*(*(Edot+5) + *(Edot+7)) * dt;
        deps[4] = 0.5*(*(Edot+2) + *(Edot+6)) * dt;
        deps[5] = 0.5*(*(Edot+1) + *(Edot+3)) * dt;

        // voigt notation
        sigmaN[0] = *(stressPtrN+0);
        sigmaN[1] = *(stressPtrN+4);
        sigmaN[2] = *(stressPtrN+8);
        sigmaN[3] = *(stressPtrN+5);
        sigmaN[4] = *(stressPtrN+2);
        sigmaN[5] = *(stressPtrN+1);

        stateN[0] = *state1PtrN;
        stateN[1] = *state2PtrN;
        stateN[2] = *state3PtrN;
        stateN[3] = *state4PtrN;
        stateN[4] = *state5PtrN;
        stateN[5] = *state6PtrN;
        stateN[6] = *state7PtrN;
        stateN[7] = *state8PtrN;
        stateN[8] = *state9PtrN;
        stateN[9] = *state10PtrN;
        stateN[10] = *state11PtrN;
        stateN[11] = *state12PtrN;
        stateN[12] = *state13PtrN;
        stateN[13] = *state14PtrN;
        stateN[14] = *state15PtrN;
        stateN[15] = *state16PtrN;
        stateN[16] = *state17PtrN;
        stateN[17] = *state18PtrN;
        stateN[18] = *state19PtrN;
        stateN[19] = *state20PtrN;
        stateN[20] = *state21PtrN;
        stateN[21] = *state22PtrN;
        stateN[22] = *state23PtrN;
        stateN[23] = *state24PtrN;
        stateN[24] = *state25PtrN;
        stateN[25] = *state26PtrN;
        stateN[26] = *state27PtrN;
        stateN[27] = *state28PtrN;
        stateN[28] = *state29PtrN;
        stateN[29] = *state30PtrN;
        stateN[30] = *state31PtrN;
        stateN[31] = *state32PtrN;
        stateN[32] = *state33PtrN;
        stateN[33] = *state34PtrN;
        stateN[34] = *state35PtrN;
        stateN[35] = *state36PtrN;
        stateN[36] = *state37PtrN;
        stateN[37] = *state38PtrN;
        stateN[38] = *state39PtrN;
        stateN[39] = *state40PtrN;
        stateN[40] = *state41PtrN;
        stateN[41] = *state42PtrN;
        stateN[42] = *state43PtrN;
        stateN[43] = *state44PtrN;
        stateN[44] = *state45PtrN;
        stateN[45] = *state46PtrN;
        stateN[46] = *state47PtrN;
        stateN[47] = *state48PtrN;
        stateN[48] = *state49PtrN;

        // call the cdpm2 function
        cdpm2material_( &dt,
                        deps, 
                        epsN,
                        sigmaN, 
                        stateN, 
                        internEnergyN, 
                        inelasEnergyN,
                        sigmaNP1, 
                        stateNP1, 
                        internEnergyNP1, 
                        inelasEnergyNP1);

        // recast the computed values back to the corresponding format
        // voigt notation
        *(stressPtrNP1+0) = sigmaNP1[0];
        *(stressPtrNP1+4) = sigmaNP1[1];
        *(stressPtrNP1+8) = sigmaNP1[2];
        *(stressPtrNP1+5) = sigmaNP1[3];
        *(stressPtrNP1+7) = sigmaNP1[3];
        *(stressPtrNP1+2) = sigmaNP1[4];
        *(stressPtrNP1+6) = sigmaNP1[4];
        *(stressPtrNP1+1) = sigmaNP1[5];
        *(stressPtrNP1+3) = sigmaNP1[5];

        *state1PtrNP1 = stateNP1[0];
        *state2PtrNP1 = stateNP1[1];
        *state3PtrNP1 = stateNP1[2];
        *state4PtrNP1 = stateNP1[3];
        *state5PtrNP1 = stateNP1[4];
        *state6PtrNP1 = stateNP1[5];
        *state7PtrNP1 = stateNP1[6];
        *state8PtrNP1 = stateNP1[7];
        *state9PtrNP1 = stateNP1[8];
        *state10PtrNP1 = stateNP1[9];
        *state11PtrNP1 = stateNP1[10];
        *state12PtrNP1 = stateNP1[11];
        *state13PtrNP1 = stateNP1[12];
        *state14PtrNP1 = stateNP1[13];
        *state15PtrNP1 = stateNP1[14];
        *state16PtrNP1 = stateNP1[15];
        *state17PtrNP1 = stateNP1[16];
        *state18PtrNP1 = stateNP1[17];
        *state19PtrNP1 = stateNP1[18];
        *state20PtrNP1 = stateNP1[19];
        *state21PtrNP1 = stateNP1[20];
        *state22PtrNP1 = stateNP1[21];
        *state23PtrNP1 = stateNP1[22];
        *state24PtrNP1 = stateNP1[23];
        *state25PtrNP1 = stateNP1[24];
        *state26PtrNP1 = stateNP1[25];
        *state27PtrNP1 = stateNP1[26];
        *state28PtrNP1 = stateNP1[27];
        *state29PtrNP1 = stateNP1[28];
        *state30PtrNP1 = stateNP1[29];
        *state31PtrNP1 = stateNP1[30];
        *state32PtrNP1 = stateNP1[31];
        *state33PtrNP1 = stateNP1[32];
        *state34PtrNP1 = stateNP1[33];
        *state35PtrNP1 = stateNP1[34];
        *state36PtrNP1 = stateNP1[35];
        *state37PtrNP1 = stateNP1[36];
        *state38PtrNP1 = stateNP1[37];
        *state39PtrNP1 = stateNP1[38];
        *state40PtrNP1 = stateNP1[39];
        *state41PtrNP1 = stateNP1[40];
        *state42PtrNP1 = stateNP1[41];
        *state43PtrNP1 = stateNP1[42];
        *state44PtrNP1 = stateNP1[43];
        *state45PtrNP1 = stateNP1[44];
        *state46PtrNP1 = stateNP1[45];
        *state47PtrNP1 = stateNP1[46];
        *state48PtrNP1 = stateNP1[47];
        *state49PtrNP1 = stateNP1[48];
  }
}

template<typename ScalarT>
void updateBondLevelCDPM2Stress
(
    const ScalarT* bondLevelStrainRateXX, 
    const ScalarT* bondLevelStrainRateXY, 
    const ScalarT* bondLevelStrainRateXZ, 
    const ScalarT* bondLevelStrainRateYX, 
    const ScalarT* bondLevelStrainRateYY, 
    const ScalarT* bondLevelStrainRateYZ, 
    const ScalarT* bondLevelStrainRateZX, 
    const ScalarT* bondLevelStrainRateZY, 
    const ScalarT* bondLevelStrainRateZZ, 
    const ScalarT* bondLevelStrainXX, 
    const ScalarT* bondLevelStrainXY, 
    const ScalarT* bondLevelStrainXZ, 
    const ScalarT* bondLevelStrainYX, 
    const ScalarT* bondLevelStrainYY, 
    const ScalarT* bondLevelStrainYZ, 
    const ScalarT* bondLevelStrainZX, 
    const ScalarT* bondLevelStrainZY, 
    const ScalarT* bondLevelStrainZZ, 
    const ScalarT* bondLevelStressXXN, 
    const ScalarT* bondLevelStressXYN, 
    const ScalarT* bondLevelStressXZN, 
    const ScalarT* bondLevelStressYXN, 
    const ScalarT* bondLevelStressYYN, 
    const ScalarT* bondLevelStressYZN, 
    const ScalarT* bondLevelStressZXN, 
    const ScalarT* bondLevelStressZYN, 
    const ScalarT* bondLevelStressZZN, 
    const ScalarT* bondLevelState1N,
    const ScalarT* bondLevelState2N,
    const ScalarT* bondLevelState3N,
    const ScalarT* bondLevelState4N,
    const ScalarT* bondLevelState5N,
    const ScalarT* bondLevelState6N,
    const ScalarT* bondLevelState7N,
    const ScalarT* bondLevelState8N,
    const ScalarT* bondLevelState9N,
    const ScalarT* bondLevelState10N,
    const ScalarT* bondLevelState11N,
    const ScalarT* bondLevelState12N,
    const ScalarT* bondLevelState13N,
    const ScalarT* bondLevelState14N,
    const ScalarT* bondLevelState15N,
    const ScalarT* bondLevelState16N,
    const ScalarT* bondLevelState17N,
    const ScalarT* bondLevelState18N,
    const ScalarT* bondLevelState19N,
    const ScalarT* bondLevelState20N,
    const ScalarT* bondLevelState21N,
    const ScalarT* bondLevelState22N,
    const ScalarT* bondLevelState23N,
    const ScalarT* bondLevelState24N,
    const ScalarT* bondLevelState25N,
    const ScalarT* bondLevelState26N,
    const ScalarT* bondLevelState27N,
    const ScalarT* bondLevelState28N,
    const ScalarT* bondLevelState29N,
    const ScalarT* bondLevelState30N,
    const ScalarT* bondLevelState31N,
    const ScalarT* bondLevelState32N,
    const ScalarT* bondLevelState33N,
    const ScalarT* bondLevelState34N,
    const ScalarT* bondLevelState35N,
    const ScalarT* bondLevelState36N,
    const ScalarT* bondLevelState37N,
    const ScalarT* bondLevelState38N,
    const ScalarT* bondLevelState39N,
    const ScalarT* bondLevelState40N,
    const ScalarT* bondLevelState41N,
    const ScalarT* bondLevelState42N,
    const ScalarT* bondLevelState43N,
    const ScalarT* bondLevelState44N,
    const ScalarT* bondLevelState45N,
    const ScalarT* bondLevelState46N,
    const ScalarT* bondLevelState47N,
    const ScalarT* bondLevelState48N,
    const ScalarT* bondLevelState49N,
    const ScalarT* bondLevelInternalEnergyN,
    const ScalarT* bondLevelInelasticEnergyN,
    ScalarT* bondLevelStressXXNP1, 
    ScalarT* bondLevelStressXYNP1, 
    ScalarT* bondLevelStressXZNP1, 
    ScalarT* bondLevelStressYXNP1, 
    ScalarT* bondLevelStressYYNP1, 
    ScalarT* bondLevelStressYZNP1, 
    ScalarT* bondLevelStressZXNP1, 
    ScalarT* bondLevelStressZYNP1, 
    ScalarT* bondLevelStressZZNP1, 
    ScalarT* bondLevelState1NP1,
    ScalarT* bondLevelState2NP1,
    ScalarT* bondLevelState3NP1,
    ScalarT* bondLevelState4NP1,
    ScalarT* bondLevelState5NP1,
    ScalarT* bondLevelState6NP1,
    ScalarT* bondLevelState7NP1,
    ScalarT* bondLevelState8NP1,
    ScalarT* bondLevelState9NP1,
    ScalarT* bondLevelState10NP1,
    ScalarT* bondLevelState11NP1,
    ScalarT* bondLevelState12NP1,
    ScalarT* bondLevelState13NP1,
    ScalarT* bondLevelState14NP1,
    ScalarT* bondLevelState15NP1,
    ScalarT* bondLevelState16NP1,
    ScalarT* bondLevelState17NP1,
    ScalarT* bondLevelState18NP1,
    ScalarT* bondLevelState19NP1,
    ScalarT* bondLevelState20NP1,
    ScalarT* bondLevelState21NP1,
    ScalarT* bondLevelState22NP1,
    ScalarT* bondLevelState23NP1,
    ScalarT* bondLevelState24NP1,
    ScalarT* bondLevelState25NP1,
    ScalarT* bondLevelState26NP1,
    ScalarT* bondLevelState27NP1,
    ScalarT* bondLevelState28NP1,
    ScalarT* bondLevelState29NP1,
    ScalarT* bondLevelState30NP1,
    ScalarT* bondLevelState31NP1,
    ScalarT* bondLevelState32NP1,
    ScalarT* bondLevelState33NP1,
    ScalarT* bondLevelState34NP1,
    ScalarT* bondLevelState35NP1,
    ScalarT* bondLevelState36NP1,
    ScalarT* bondLevelState37NP1,
    ScalarT* bondLevelState38NP1,
    ScalarT* bondLevelState39NP1,
    ScalarT* bondLevelState40NP1,
    ScalarT* bondLevelState41NP1,
    ScalarT* bondLevelState42NP1,
    ScalarT* bondLevelState43NP1,
    ScalarT* bondLevelState44NP1,
    ScalarT* bondLevelState45NP1,
    ScalarT* bondLevelState46NP1,
    ScalarT* bondLevelState47NP1,
    ScalarT* bondLevelState48NP1,
    ScalarT* bondLevelState49NP1,
    ScalarT* bondLevelInternalEnergyNP1,
    ScalarT* bondLevelInelasticEnergyNP1,
    const int* neighborhoodList,
    const int numPoints, 
    const double dt
)
{
  const ScalarT* EdotXX = bondLevelStrainRateXX;
  const ScalarT* EdotXY = bondLevelStrainRateXY;
  const ScalarT* EdotXZ = bondLevelStrainRateXZ;
  const ScalarT* EdotYX = bondLevelStrainRateYX;
  const ScalarT* EdotYY = bondLevelStrainRateYY;
  const ScalarT* EdotYZ = bondLevelStrainRateYZ;
  const ScalarT* EdotZX = bondLevelStrainRateZX;
  const ScalarT* EdotZY = bondLevelStrainRateZY;
  const ScalarT* EdotZZ = bondLevelStrainRateZZ;

  const ScalarT* strainXX = bondLevelStrainXX;
  const ScalarT* strainXY = bondLevelStrainXY;
  const ScalarT* strainXZ = bondLevelStrainXZ;
  const ScalarT* strainYX = bondLevelStrainYX;
  const ScalarT* strainYY = bondLevelStrainYY;
  const ScalarT* strainYZ = bondLevelStrainYZ;
  const ScalarT* strainZX = bondLevelStrainZX;
  const ScalarT* strainZY = bondLevelStrainZY;
  const ScalarT* strainZZ = bondLevelStrainZZ;

  const ScalarT* stressXXN = bondLevelStressXXN;
  const ScalarT* stressXYN = bondLevelStressXYN;
  const ScalarT* stressXZN = bondLevelStressXZN;
  const ScalarT* stressYXN = bondLevelStressYXN;
  const ScalarT* stressYYN = bondLevelStressYYN;
  const ScalarT* stressYZN = bondLevelStressYZN;
  const ScalarT* stressZXN = bondLevelStressZXN;
  const ScalarT* stressZYN = bondLevelStressZYN;
  const ScalarT* stressZZN = bondLevelStressZZN;

  const ScalarT* state1PtrN = bondLevelState1N;
  const ScalarT* state2PtrN = bondLevelState2N;
  const ScalarT* state3PtrN = bondLevelState3N;
  const ScalarT* state4PtrN = bondLevelState4N;
  const ScalarT* state5PtrN = bondLevelState5N;
  const ScalarT* state6PtrN = bondLevelState6N;
  const ScalarT* state7PtrN = bondLevelState7N;
  const ScalarT* state8PtrN = bondLevelState8N;
  const ScalarT* state9PtrN = bondLevelState9N;
  const ScalarT* state10PtrN = bondLevelState10N;
  const ScalarT* state11PtrN = bondLevelState11N;
  const ScalarT* state12PtrN = bondLevelState12N;
  const ScalarT* state13PtrN = bondLevelState13N;
  const ScalarT* state14PtrN = bondLevelState14N;
  const ScalarT* state15PtrN = bondLevelState15N;
  const ScalarT* state16PtrN = bondLevelState16N;
  const ScalarT* state17PtrN = bondLevelState17N;
  const ScalarT* state18PtrN = bondLevelState18N;
  const ScalarT* state19PtrN = bondLevelState19N;
  const ScalarT* state20PtrN = bondLevelState20N;
  const ScalarT* state21PtrN = bondLevelState21N;
  const ScalarT* state22PtrN = bondLevelState22N;
  const ScalarT* state23PtrN = bondLevelState23N;
  const ScalarT* state24PtrN = bondLevelState24N;
  const ScalarT* state25PtrN = bondLevelState25N;
  const ScalarT* state26PtrN = bondLevelState26N;
  const ScalarT* state27PtrN = bondLevelState27N;
  const ScalarT* state28PtrN = bondLevelState28N;
  const ScalarT* state29PtrN = bondLevelState29N;
  const ScalarT* state30PtrN = bondLevelState30N;
  const ScalarT* state31PtrN = bondLevelState31N;
  const ScalarT* state32PtrN = bondLevelState32N;
  const ScalarT* state33PtrN = bondLevelState33N;
  const ScalarT* state34PtrN = bondLevelState34N;
  const ScalarT* state35PtrN = bondLevelState35N;
  const ScalarT* state36PtrN = bondLevelState36N;
  const ScalarT* state37PtrN = bondLevelState37N;
  const ScalarT* state38PtrN = bondLevelState38N;
  const ScalarT* state39PtrN = bondLevelState39N;
  const ScalarT* state40PtrN = bondLevelState40N;
  const ScalarT* state41PtrN = bondLevelState41N;
  const ScalarT* state42PtrN = bondLevelState42N;
  const ScalarT* state43PtrN = bondLevelState43N;
  const ScalarT* state44PtrN = bondLevelState44N;
  const ScalarT* state45PtrN = bondLevelState45N;
  const ScalarT* state46PtrN = bondLevelState46N;
  const ScalarT* state47PtrN = bondLevelState47N;
  const ScalarT* state48PtrN = bondLevelState48N;
  const ScalarT* state49PtrN = bondLevelState49N;

  const ScalarT* internEnergyN = bondLevelInternalEnergyN;
  const ScalarT* inelasEnergyN = bondLevelInelasticEnergyN;

  ScalarT* stressXXNP1 = bondLevelStressXXNP1;
  ScalarT* stressXYNP1 = bondLevelStressXYNP1;
  ScalarT* stressXZNP1 = bondLevelStressXZNP1;
  ScalarT* stressYXNP1 = bondLevelStressYXNP1;
  ScalarT* stressYYNP1 = bondLevelStressYYNP1;
  ScalarT* stressYZNP1 = bondLevelStressYZNP1;
  ScalarT* stressZXNP1 = bondLevelStressZXNP1;
  ScalarT* stressZYNP1 = bondLevelStressZYNP1;
  ScalarT* stressZZNP1 = bondLevelStressZZNP1;

  ScalarT* state1PtrNP1 = bondLevelState1NP1;
  ScalarT* state2PtrNP1 = bondLevelState2NP1;
  ScalarT* state3PtrNP1 = bondLevelState3NP1;
  ScalarT* state4PtrNP1 = bondLevelState4NP1;
  ScalarT* state5PtrNP1 = bondLevelState5NP1;
  ScalarT* state6PtrNP1 = bondLevelState6NP1;
  ScalarT* state7PtrNP1 = bondLevelState7NP1;
  ScalarT* state8PtrNP1 = bondLevelState8NP1;
  ScalarT* state9PtrNP1 = bondLevelState9NP1;
  ScalarT* state10PtrNP1 = bondLevelState10NP1;
  ScalarT* state11PtrNP1 = bondLevelState11NP1;
  ScalarT* state12PtrNP1 = bondLevelState12NP1;
  ScalarT* state13PtrNP1 = bondLevelState13NP1;
  ScalarT* state14PtrNP1 = bondLevelState14NP1;
  ScalarT* state15PtrNP1 = bondLevelState15NP1;
  ScalarT* state16PtrNP1 = bondLevelState16NP1;
  ScalarT* state17PtrNP1 = bondLevelState17NP1;
  ScalarT* state18PtrNP1 = bondLevelState18NP1;
  ScalarT* state19PtrNP1 = bondLevelState19NP1;
  ScalarT* state20PtrNP1 = bondLevelState20NP1;
  ScalarT* state21PtrNP1 = bondLevelState21NP1;
  ScalarT* state22PtrNP1 = bondLevelState22NP1;
  ScalarT* state23PtrNP1 = bondLevelState23NP1;
  ScalarT* state24PtrNP1 = bondLevelState24NP1;
  ScalarT* state25PtrNP1 = bondLevelState25NP1;
  ScalarT* state26PtrNP1 = bondLevelState26NP1;
  ScalarT* state27PtrNP1 = bondLevelState27NP1;
  ScalarT* state28PtrNP1 = bondLevelState28NP1;
  ScalarT* state29PtrNP1 = bondLevelState29NP1;
  ScalarT* state30PtrNP1 = bondLevelState30NP1;
  ScalarT* state31PtrNP1 = bondLevelState31NP1;
  ScalarT* state32PtrNP1 = bondLevelState32NP1;
  ScalarT* state33PtrNP1 = bondLevelState33NP1;
  ScalarT* state34PtrNP1 = bondLevelState34NP1;
  ScalarT* state35PtrNP1 = bondLevelState35NP1;
  ScalarT* state36PtrNP1 = bondLevelState36NP1;
  ScalarT* state37PtrNP1 = bondLevelState37NP1;
  ScalarT* state38PtrNP1 = bondLevelState38NP1;
  ScalarT* state39PtrNP1 = bondLevelState39NP1;
  ScalarT* state40PtrNP1 = bondLevelState40NP1;
  ScalarT* state41PtrNP1 = bondLevelState41NP1;
  ScalarT* state42PtrNP1 = bondLevelState42NP1;
  ScalarT* state43PtrNP1 = bondLevelState43NP1;
  ScalarT* state44PtrNP1 = bondLevelState44NP1;
  ScalarT* state45PtrNP1 = bondLevelState45NP1;
  ScalarT* state46PtrNP1 = bondLevelState46NP1;
  ScalarT* state47PtrNP1 = bondLevelState47NP1;
  ScalarT* state48PtrNP1 = bondLevelState48NP1;
  ScalarT* state49PtrNP1 = bondLevelState49NP1;

  ScalarT* internEnergyNP1 = bondLevelInternalEnergyNP1;
  ScalarT* inelasEnergyNP1 = bondLevelInelasticEnergyNP1;

  ScalarT deps[6];
  ScalarT epsN[6];

  ScalarT sigmaN[6];
  ScalarT sigmaNP1[6];

  ScalarT stateN[49];
  ScalarT stateNP1[49];


  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

    // All is bond level.
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
        EdotXX++, EdotXY++, EdotXZ++, 
        EdotYX++, EdotYY++, EdotYZ++, 
        EdotZX++, EdotZY++, EdotZZ++,
        strainXX++, strainXY++, strainXZ++, 
        strainYX++, strainYY++, strainYZ++, 
        strainZX++, strainZY++, strainZZ++,
        stressXXN++, stressXYN++, stressXZN++, 
        stressYXN++, stressYYN++, stressYZN++, 
        stressZXN++, stressZYN++, stressZZN++,
        ++state1PtrN, ++state2PtrN, ++state3PtrN, ++state4PtrN, ++state5PtrN,
        ++state6PtrN, ++state7PtrN, ++state8PtrN, ++state9PtrN, ++state10PtrN,
        ++state11PtrN, ++state12PtrN, ++state13PtrN, ++state14PtrN, ++state15PtrN,
        ++state16PtrN, ++state17PtrN, ++state18PtrN, ++state19PtrN, ++state20PtrN,
        ++state21PtrN, ++state22PtrN, ++state23PtrN, ++state24PtrN, ++state25PtrN,
        ++state26PtrN, ++state27PtrN, ++state28PtrN, ++state29PtrN, ++state30PtrN,
        ++state31PtrN, ++state32PtrN, ++state33PtrN, ++state34PtrN, ++state35PtrN,
        ++state36PtrN, ++state37PtrN, ++state38PtrN, ++state39PtrN, ++state40PtrN,
        ++state41PtrN, ++state42PtrN, ++state43PtrN, ++state44PtrN, ++state45PtrN,
        ++state46PtrN, ++state47PtrN, ++state48PtrN, ++state49PtrN,
        internEnergyN++, inelasEnergyN++, 
        stressXXNP1++, stressXYNP1++, stressXZNP1++, 
        stressYXNP1++, stressYYNP1++, stressYZNP1++, 
        stressZXNP1++, stressZYNP1++, stressZZNP1++,
        ++state1PtrNP1, ++state2PtrNP1, ++state3PtrNP1, ++state4PtrNP1, ++state5PtrNP1,
        ++state6PtrNP1, ++state7PtrNP1, ++state8PtrNP1, ++state9PtrNP1, ++state10PtrNP1,
        ++state11PtrNP1, ++state12PtrNP1, ++state13PtrNP1, ++state14PtrNP1, ++state15PtrNP1,
        ++state16PtrNP1, ++state17PtrNP1, ++state18PtrNP1, ++state19PtrNP1, ++state20PtrNP1,
        ++state21PtrNP1, ++state22PtrNP1, ++state23PtrNP1, ++state24PtrNP1, ++state25PtrNP1,
        ++state26PtrNP1, ++state27PtrNP1, ++state28PtrNP1, ++state29PtrNP1, ++state30PtrNP1,
        ++state31PtrNP1, ++state32PtrNP1, ++state33PtrNP1, ++state34PtrNP1, ++state35PtrNP1,
        ++state36PtrNP1, ++state37PtrNP1, ++state38PtrNP1, ++state39PtrNP1, ++state40PtrNP1,
        ++state41PtrNP1, ++state42PtrNP1, ++state43PtrNP1, ++state44PtrNP1, ++state45PtrNP1,
        ++state46PtrNP1, ++state47PtrNP1, ++state48PtrNP1, ++state49PtrNP1,
        internEnergyNP1++, inelasEnergyNP1++
        ){

        // voigt notation
        epsN[0] = *strainXX;
        epsN[1] = *strainYY;
        epsN[2] = *strainZZ;
        epsN[3] = *strainYZ;
        epsN[4] = *strainXZ;
        epsN[5] = *strainXY;

        // voigt notation
        deps[0] = *EdotXX * dt;
        deps[1] = *EdotYY * dt;
        deps[2] = *EdotZZ * dt;
        deps[3] = *EdotYZ * dt;
        deps[4] = *EdotXZ * dt;
        deps[5] = *EdotXY * dt;

        // voigt notation
        sigmaN[0] = *stressXXN;
        sigmaN[1] = *stressYYN;
        sigmaN[2] = *stressZZN;
        sigmaN[3] = *stressYZN;
        sigmaN[4] = *stressXZN;
        sigmaN[5] = *stressXYN;

        stateN[0] = *state1PtrN;
        stateN[1] = *state2PtrN;
        stateN[2] = *state3PtrN;
        stateN[3] = *state4PtrN;
        stateN[4] = *state5PtrN;
        stateN[5] = *state6PtrN;
        stateN[6] = *state7PtrN;
        stateN[7] = *state8PtrN;
        stateN[8] = *state9PtrN;
        stateN[9] = *state10PtrN;
        stateN[10] = *state11PtrN;
        stateN[11] = *state12PtrN;
        stateN[12] = *state13PtrN;
        stateN[13] = *state14PtrN;
        stateN[14] = *state15PtrN;
        stateN[15] = *state16PtrN;
        stateN[16] = *state17PtrN;
        stateN[17] = *state18PtrN;
        stateN[18] = *state19PtrN;
        stateN[19] = *state20PtrN;
        stateN[20] = *state21PtrN;
        stateN[21] = *state22PtrN;
        stateN[22] = *state23PtrN;
        stateN[23] = *state24PtrN;
        stateN[24] = *state25PtrN;
        stateN[25] = *state26PtrN;
        stateN[26] = *state27PtrN;
        stateN[27] = *state28PtrN;
        stateN[28] = *state29PtrN;
        stateN[29] = *state30PtrN;
        stateN[30] = *state31PtrN;
        stateN[31] = *state32PtrN;
        stateN[32] = *state33PtrN;
        stateN[33] = *state34PtrN;
        stateN[34] = *state35PtrN;
        stateN[35] = *state36PtrN;
        stateN[36] = *state37PtrN;
        stateN[37] = *state38PtrN;
        stateN[38] = *state39PtrN;
        stateN[39] = *state40PtrN;
        stateN[40] = *state41PtrN;
        stateN[41] = *state42PtrN;
        stateN[42] = *state43PtrN;
        stateN[43] = *state44PtrN;
        stateN[44] = *state45PtrN;
        stateN[45] = *state46PtrN;
        stateN[46] = *state47PtrN;
        stateN[47] = *state48PtrN;
        stateN[48] = *state49PtrN;

        // call the cdpm2 function
        cdpm2material_( &dt, 
                         deps, 
                         epsN,
                         sigmaN, 
                         stateN, 
                         internEnergyN, 
                         inelasEnergyN,
                         sigmaNP1, 
                         stateNP1, 
                         internEnergyNP1, 
                         inelasEnergyNP1);

        // recast the computed values back to the corresponding format
        // voigt notation
        *stressXXNP1 = sigmaNP1[0];
        *stressYYNP1 = sigmaNP1[1];
        *stressZZNP1 = sigmaNP1[2];
        *stressYZNP1 = sigmaNP1[3];
        *stressZYNP1 = sigmaNP1[3];
        *stressXZNP1 = sigmaNP1[4];
        *stressZXNP1 = sigmaNP1[4];
        *stressXYNP1 = sigmaNP1[5];
        *stressYXNP1 = sigmaNP1[5];

        *state1PtrNP1 = stateNP1[0];
        *state2PtrNP1 = stateNP1[1];
        *state3PtrNP1 = stateNP1[2];
        *state4PtrNP1 = stateNP1[3];
        *state5PtrNP1 = stateNP1[4];
        *state6PtrNP1 = stateNP1[5];
        *state7PtrNP1 = stateNP1[6];
        *state8PtrNP1 = stateNP1[7];
        *state9PtrNP1 = stateNP1[8];
        *state10PtrNP1 = stateNP1[9];
        *state11PtrNP1 = stateNP1[10];
        *state12PtrNP1 = stateNP1[11];
        *state13PtrNP1 = stateNP1[12];
        *state14PtrNP1 = stateNP1[13];
        *state15PtrNP1 = stateNP1[14];
        *state16PtrNP1 = stateNP1[15];
        *state17PtrNP1 = stateNP1[16];
        *state18PtrNP1 = stateNP1[17];
        *state19PtrNP1 = stateNP1[18];
        *state20PtrNP1 = stateNP1[19];
        *state21PtrNP1 = stateNP1[20];
        *state22PtrNP1 = stateNP1[21];
        *state23PtrNP1 = stateNP1[22];
        *state24PtrNP1 = stateNP1[23];
        *state25PtrNP1 = stateNP1[24];
        *state26PtrNP1 = stateNP1[25];
        *state27PtrNP1 = stateNP1[26];
        *state28PtrNP1 = stateNP1[27];
        *state29PtrNP1 = stateNP1[28];
        *state30PtrNP1 = stateNP1[29];
        *state31PtrNP1 = stateNP1[30];
        *state32PtrNP1 = stateNP1[31];
        *state33PtrNP1 = stateNP1[32];
        *state34PtrNP1 = stateNP1[33];
        *state35PtrNP1 = stateNP1[34];
        *state36PtrNP1 = stateNP1[35];
        *state37PtrNP1 = stateNP1[36];
        *state38PtrNP1 = stateNP1[37];
        *state39PtrNP1 = stateNP1[38];
        *state40PtrNP1 = stateNP1[39];
        *state41PtrNP1 = stateNP1[40];
        *state42PtrNP1 = stateNP1[41];
        *state43PtrNP1 = stateNP1[42];
        *state44PtrNP1 = stateNP1[43];
        *state45PtrNP1 = stateNP1[44];
        *state46PtrNP1 = stateNP1[45];
        *state47PtrNP1 = stateNP1[46];
        *state48PtrNP1 = stateNP1[47];
        *state49PtrNP1 = stateNP1[48];
    }
  }
}


// Explicit template instantiation for double
template void updateCDPM2Stress<double>
(
    const double* strainRate, 
    const double* strain, 
    const double* stressN, 
    const double* state1N,
    const double* state2N,
    const double* state3N,
    const double* state4N,
    const double* state5N,
    const double* state6N,
    const double* state7N,
    const double* state8N,
    const double* state9N,
    const double* state10N,
    const double* state11N,
    const double* state12N,
    const double* state13N,
    const double* state14N,
    const double* state15N,
    const double* state16N,
    const double* state17N,
    const double* state18N,
    const double* state19N,
    const double* state20N,
    const double* state21N,
    const double* state22N,
    const double* state23N,
    const double* state24N,
    const double* state25N,
    const double* state26N,
    const double* state27N,
    const double* state28N,
    const double* state29N,
    const double* state30N,
    const double* state31N,
    const double* state32N,
    const double* state33N,
    const double* state34N,
    const double* state35N,
    const double* state36N,
    const double* state37N,
    const double* state38N,
    const double* state39N,
    const double* state40N,
    const double* state41N,
    const double* state42N,
    const double* state43N,
    const double* state44N,
    const double* state45N,
    const double* state46N,
    const double* state47N,
    const double* state48N,
    const double* state49N,
    const double* internalEnergyN,
    const double* inelasticEnergyN,
    double* stressNP1, 
    double* state1NP1,
    double* state2NP1,
    double* state3NP1,
    double* state4NP1,
    double* state5NP1,
    double* state6NP1,
    double* state7NP1,
    double* state8NP1,
    double* state9NP1,
    double* state10NP1,
    double* state11NP1,
    double* state12NP1,
    double* state13NP1,
    double* state14NP1,
    double* state15NP1,
    double* state16NP1,
    double* state17NP1,
    double* state18NP1,
    double* state19NP1,
    double* state20NP1,
    double* state21NP1,
    double* state22NP1,
    double* state23NP1,
    double* state24NP1,
    double* state25NP1,
    double* state26NP1,
    double* state27NP1,
    double* state28NP1,
    double* state29NP1,
    double* state30NP1,
    double* state31NP1,
    double* state32NP1,
    double* state33NP1,
    double* state34NP1,
    double* state35NP1,
    double* state36NP1,
    double* state37NP1,
    double* state38NP1,
    double* state39NP1,
    double* state40NP1,
    double* state41NP1,
    double* state42NP1,
    double* state43NP1,
    double* state44NP1,
    double* state45NP1,
    double* state46NP1,
    double* state47NP1,
    double* state48NP1,
    double* state49NP1,
    double* internalEnergyNP1,
    double* inelasticEnergyNP1,
    const int numPoints, 
    const double dt
);

template void updateBondLevelCDPM2Stress<double>
(
    const double* bondLevelStrainRateXX, 
    const double* bondLevelStrainRateXY, 
    const double* bondLevelStrainRateXZ, 
    const double* bondLevelStrainRateYX, 
    const double* bondLevelStrainRateYY, 
    const double* bondLevelStrainRateYZ, 
    const double* bondLevelStrainRateZX, 
    const double* bondLevelStrainRateZY, 
    const double* bondLevelStrainRateZZ, 
    const double* bondLevelStrainXX, 
    const double* bondLevelStrainXY, 
    const double* bondLevelStrainXZ, 
    const double* bondLevelStrainYX, 
    const double* bondLevelStrainYY, 
    const double* bondLevelStrainYZ, 
    const double* bondLevelStrainZX, 
    const double* bondLevelStrainZY, 
    const double* bondLevelStrainZZ, 
    const double* bondLevelStressXXN, 
    const double* bondLevelStressXYN, 
    const double* bondLevelStressXZN, 
    const double* bondLevelStressYXN, 
    const double* bondLevelStressYYN, 
    const double* bondLevelStressYZN, 
    const double* bondLevelStressZXN, 
    const double* bondLevelStressZYN, 
    const double* bondLevelStressZZN, 
    const double* bondLevelState1N,
    const double* bondLevelState2N,
    const double* bondLevelState3N,
    const double* bondLevelState4N,
    const double* bondLevelState5N,
    const double* bondLevelState6N,
    const double* bondLevelState7N,
    const double* bondLevelState8N,
    const double* bondLevelState9N,
    const double* bondLevelState10N,
    const double* bondLevelState11N,
    const double* bondLevelState12N,
    const double* bondLevelState13N,
    const double* bondLevelState14N,
    const double* bondLevelState15N,
    const double* bondLevelState16N,
    const double* bondLevelState17N,
    const double* bondLevelState18N,
    const double* bondLevelState19N,
    const double* bondLevelState20N,
    const double* bondLevelState21N,
    const double* bondLevelState22N,
    const double* bondLevelState23N,
    const double* bondLevelState24N,
    const double* bondLevelState25N,
    const double* bondLevelState26N,
    const double* bondLevelState27N,
    const double* bondLevelState28N,
    const double* bondLevelState29N,
    const double* bondLevelState30N,
    const double* bondLevelState31N,
    const double* bondLevelState32N,
    const double* bondLevelState33N,
    const double* bondLevelState34N,
    const double* bondLevelState35N,
    const double* bondLevelState36N,
    const double* bondLevelState37N,
    const double* bondLevelState38N,
    const double* bondLevelState39N,
    const double* bondLevelState40N,
    const double* bondLevelState41N,
    const double* bondLevelState42N,
    const double* bondLevelState43N,
    const double* bondLevelState44N,
    const double* bondLevelState45N,
    const double* bondLevelState46N,
    const double* bondLevelState47N,
    const double* bondLevelState48N,
    const double* bondLevelState49N,
    const double* bondLevelInternalEnergyN,
    const double* bondLevelInelasticEnergyN,
    double* bondLevelStressXXNP1, 
    double* bondLevelStressXYNP1, 
    double* bondLevelStressXZNP1, 
    double* bondLevelStressYXNP1, 
    double* bondLevelStressYYNP1, 
    double* bondLevelStressYZNP1, 
    double* bondLevelStressZXNP1, 
    double* bondLevelStressZYNP1, 
    double* bondLevelStressZZNP1, 
    double* bondLevelState1NP1,
    double* bondLevelState2NP1,
    double* bondLevelState3NP1,
    double* bondLevelState4NP1,
    double* bondLevelState5NP1,
    double* bondLevelState6NP1,
    double* bondLevelState7NP1,
    double* bondLevelState8NP1,
    double* bondLevelState9NP1,
    double* bondLevelState10NP1,
    double* bondLevelState11NP1,
    double* bondLevelState12NP1,
    double* bondLevelState13NP1,
    double* bondLevelState14NP1,
    double* bondLevelState15NP1,
    double* bondLevelState16NP1,
    double* bondLevelState17NP1,
    double* bondLevelState18NP1,
    double* bondLevelState19NP1,
    double* bondLevelState20NP1,
    double* bondLevelState21NP1,
    double* bondLevelState22NP1,
    double* bondLevelState23NP1,
    double* bondLevelState24NP1,
    double* bondLevelState25NP1,
    double* bondLevelState26NP1,
    double* bondLevelState27NP1,
    double* bondLevelState28NP1,
    double* bondLevelState29NP1,
    double* bondLevelState30NP1,
    double* bondLevelState31NP1,
    double* bondLevelState32NP1,
    double* bondLevelState33NP1,
    double* bondLevelState34NP1,
    double* bondLevelState35NP1,
    double* bondLevelState36NP1,
    double* bondLevelState37NP1,
    double* bondLevelState38NP1,
    double* bondLevelState39NP1,
    double* bondLevelState40NP1,
    double* bondLevelState41NP1,
    double* bondLevelState42NP1,
    double* bondLevelState43NP1,
    double* bondLevelState44NP1,
    double* bondLevelState45NP1,
    double* bondLevelState46NP1,
    double* bondLevelState47NP1,
    double* bondLevelState48NP1,
    double* bondLevelState49NP1,
    double* bondLevelInternalEnergyNP1,
    double* bondLevelInelasticEnergyNP1,
    const int* neighborhoodList,
    const int numPoints, 
    const double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */

}
