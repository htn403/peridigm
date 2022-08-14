//! \file cdpm2_correspondence.h

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
#ifndef CDPM2_CORRESPONDENCE_H
#define CDPM2_CORRESPONDENCE_H

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
);

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
);

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
);

}

#endif // CDPM2_CORRESPONDENCE_H
