//! \file microplane_correspondence.cxx

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

#include "microplane_correspondence.h"
#include "correspondence.h"
//#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h> // for sqrt

extern"C" {
void inputparams_(const double *E, 
                  const double *nu, 
                  const double *k1, 
                  const double *k2, 
                  const double *k3, 
                  const double *k4);
}

extern"C" {
void setsystem_(void);
}

extern"C" {
void m7fmaterial_(const double *dt,
                  double *strainInc,
                  double *epsOld,
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

void initializeMicroplaneM7Model
(
    const double E, 
    const double nu, 
    const double k1, 
    const double k2, 
    const double k3, 
    const double k4 
)
{
  inputparams_(&E, 
               &nu, 
               &k1, 
               &k2, 
               &k3, 
               &k4);

  setsystem_();
}


template<typename ScalarT>
void updateMicroplaneM7Stress
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
    const ScalarT* state50N,
    const ScalarT* state51N,
    const ScalarT* state52N,
    const ScalarT* state53N,
    const ScalarT* state54N,
    const ScalarT* state55N,
    const ScalarT* state56N,
    const ScalarT* state57N,
    const ScalarT* state58N,
    const ScalarT* state59N,
    const ScalarT* state60N,
    const ScalarT* state61N,
    const ScalarT* state62N,
    const ScalarT* state63N,
    const ScalarT* state64N,
    const ScalarT* state65N,
    const ScalarT* state66N,
    const ScalarT* state67N,
    const ScalarT* state68N,
    const ScalarT* state69N,
    const ScalarT* state70N,
    const ScalarT* state71N,
    const ScalarT* state72N,
    const ScalarT* state73N,
    const ScalarT* state74N,
    const ScalarT* state75N,
    const ScalarT* state76N,
    const ScalarT* state77N,
    const ScalarT* state78N,
    const ScalarT* state79N,
    const ScalarT* state80N,
    const ScalarT* state81N,
    const ScalarT* state82N,
    const ScalarT* state83N,
    const ScalarT* state84N,
    const ScalarT* state85N,
    const ScalarT* state86N,
    const ScalarT* state87N,
    const ScalarT* state88N,
    const ScalarT* state89N,
    const ScalarT* state90N,
    const ScalarT* state91N,
    const ScalarT* state92N,
    const ScalarT* state93N,
    const ScalarT* state94N,
    const ScalarT* state95N,
    const ScalarT* state96N,
    const ScalarT* state97N,
    const ScalarT* state98N,
    const ScalarT* state99N,
    const ScalarT* state100N,
    const ScalarT* state101N,
    const ScalarT* state102N,
    const ScalarT* state103N,
    const ScalarT* state104N,
    const ScalarT* state105N,
    const ScalarT* state106N,
    const ScalarT* state107N,
    const ScalarT* state108N,
    const ScalarT* state109N,
    const ScalarT* state110N,
    const ScalarT* state111N,
    const ScalarT* state112N,
    const ScalarT* state113N,
    const ScalarT* state114N,
    const ScalarT* state115N,
    const ScalarT* state116N,
    const ScalarT* state117N,
    const ScalarT* state118N,
    const ScalarT* state119N,
    const ScalarT* state120N,
    const ScalarT* state121N,
    const ScalarT* state122N,
    const ScalarT* state123N,
    const ScalarT* state124N,
    const ScalarT* state125N,
    const ScalarT* state126N,
    const ScalarT* state127N,
    const ScalarT* state128N,
    const ScalarT* state129N,
    const ScalarT* state130N,
    const ScalarT* state131N,
    const ScalarT* state132N,
    const ScalarT* state133N,
    const ScalarT* state134N,
    const ScalarT* state135N,
    const ScalarT* state136N,
    const ScalarT* state137N,
    const ScalarT* state138N,
    const ScalarT* state139N,
    const ScalarT* state140N,
    const ScalarT* state141N,
    const ScalarT* state142N,
    const ScalarT* state143N,
    const ScalarT* state144N,
    const ScalarT* state145N,
    const ScalarT* state146N,
    const ScalarT* state147N,
    const ScalarT* state148N,
    const ScalarT* state149N,
    const ScalarT* state150N,
    const ScalarT* state151N,
    const ScalarT* state152N,
    const ScalarT* state153N,
    const ScalarT* state154N,
    const ScalarT* state155N,
    const ScalarT* state156N,
    const ScalarT* state157N,
    const ScalarT* state158N,
    const ScalarT* state159N,
    const ScalarT* state160N,
    const ScalarT* state161N,
    const ScalarT* state162N,
    const ScalarT* state163N,
    const ScalarT* state164N,
    const ScalarT* state165N,
    const ScalarT* state166N,
    const ScalarT* state167N,
    const ScalarT* state168N,
    const ScalarT* state169N,
    const ScalarT* state170N,
    const ScalarT* state171N,
    const ScalarT* state172N,
    const ScalarT* state173N,
    const ScalarT* state174N,
    const ScalarT* state175N,
    const ScalarT* state176N,
    const ScalarT* state177N,
    const ScalarT* state178N,
    const ScalarT* state179N,
    const ScalarT* state180N,
    const ScalarT* state181N,
    const ScalarT* state182N,
    const ScalarT* state183N,
    const ScalarT* state184N,
    const ScalarT* state185N,
    const ScalarT* state186N,
    const ScalarT* state187N,
    const ScalarT* state188N,
    const ScalarT* state189N,
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
    ScalarT* state50NP1,
    ScalarT* state51NP1,
    ScalarT* state52NP1,
    ScalarT* state53NP1,
    ScalarT* state54NP1,
    ScalarT* state55NP1,
    ScalarT* state56NP1,
    ScalarT* state57NP1,
    ScalarT* state58NP1,
    ScalarT* state59NP1,
    ScalarT* state60NP1,
    ScalarT* state61NP1,
    ScalarT* state62NP1,
    ScalarT* state63NP1,
    ScalarT* state64NP1,
    ScalarT* state65NP1,
    ScalarT* state66NP1,
    ScalarT* state67NP1,
    ScalarT* state68NP1,
    ScalarT* state69NP1,
    ScalarT* state70NP1,
    ScalarT* state71NP1,
    ScalarT* state72NP1,
    ScalarT* state73NP1,
    ScalarT* state74NP1,
    ScalarT* state75NP1,
    ScalarT* state76NP1,
    ScalarT* state77NP1,
    ScalarT* state78NP1,
    ScalarT* state79NP1,
    ScalarT* state80NP1,
    ScalarT* state81NP1,
    ScalarT* state82NP1,
    ScalarT* state83NP1,
    ScalarT* state84NP1,
    ScalarT* state85NP1,
    ScalarT* state86NP1,
    ScalarT* state87NP1,
    ScalarT* state88NP1,
    ScalarT* state89NP1,
    ScalarT* state90NP1,
    ScalarT* state91NP1,
    ScalarT* state92NP1,
    ScalarT* state93NP1,
    ScalarT* state94NP1,
    ScalarT* state95NP1,
    ScalarT* state96NP1,
    ScalarT* state97NP1,
    ScalarT* state98NP1,
    ScalarT* state99NP1,
    ScalarT* state100NP1,
    ScalarT* state101NP1,
    ScalarT* state102NP1,
    ScalarT* state103NP1,
    ScalarT* state104NP1,
    ScalarT* state105NP1,
    ScalarT* state106NP1,
    ScalarT* state107NP1,
    ScalarT* state108NP1,
    ScalarT* state109NP1,
    ScalarT* state110NP1,
    ScalarT* state111NP1,
    ScalarT* state112NP1,
    ScalarT* state113NP1,
    ScalarT* state114NP1,
    ScalarT* state115NP1,
    ScalarT* state116NP1,
    ScalarT* state117NP1,
    ScalarT* state118NP1,
    ScalarT* state119NP1,
    ScalarT* state120NP1,
    ScalarT* state121NP1,
    ScalarT* state122NP1,
    ScalarT* state123NP1,
    ScalarT* state124NP1,
    ScalarT* state125NP1,
    ScalarT* state126NP1,
    ScalarT* state127NP1,
    ScalarT* state128NP1,
    ScalarT* state129NP1,
    ScalarT* state130NP1,
    ScalarT* state131NP1,
    ScalarT* state132NP1,
    ScalarT* state133NP1,
    ScalarT* state134NP1,
    ScalarT* state135NP1,
    ScalarT* state136NP1,
    ScalarT* state137NP1,
    ScalarT* state138NP1,
    ScalarT* state139NP1,
    ScalarT* state140NP1,
    ScalarT* state141NP1,
    ScalarT* state142NP1,
    ScalarT* state143NP1,
    ScalarT* state144NP1,
    ScalarT* state145NP1,
    ScalarT* state146NP1,
    ScalarT* state147NP1,
    ScalarT* state148NP1,
    ScalarT* state149NP1,
    ScalarT* state150NP1,
    ScalarT* state151NP1,
    ScalarT* state152NP1,
    ScalarT* state153NP1,
    ScalarT* state154NP1,
    ScalarT* state155NP1,
    ScalarT* state156NP1,
    ScalarT* state157NP1,
    ScalarT* state158NP1,
    ScalarT* state159NP1,
    ScalarT* state160NP1,
    ScalarT* state161NP1,
    ScalarT* state162NP1,
    ScalarT* state163NP1,
    ScalarT* state164NP1,
    ScalarT* state165NP1,
    ScalarT* state166NP1,
    ScalarT* state167NP1,
    ScalarT* state168NP1,
    ScalarT* state169NP1,
    ScalarT* state170NP1,
    ScalarT* state171NP1,
    ScalarT* state172NP1,
    ScalarT* state173NP1,
    ScalarT* state174NP1,
    ScalarT* state175NP1,
    ScalarT* state176NP1,
    ScalarT* state177NP1,
    ScalarT* state178NP1,
    ScalarT* state179NP1,
    ScalarT* state180NP1,
    ScalarT* state181NP1,
    ScalarT* state182NP1,
    ScalarT* state183NP1,
    ScalarT* state184NP1,
    ScalarT* state185NP1,
    ScalarT* state186NP1,
    ScalarT* state187NP1,
    ScalarT* state188NP1,
    ScalarT* state189NP1,
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
  const ScalarT* state50PtrN = state50N;
  const ScalarT* state51PtrN = state51N;
  const ScalarT* state52PtrN = state52N;
  const ScalarT* state53PtrN = state53N;
  const ScalarT* state54PtrN = state54N;
  const ScalarT* state55PtrN = state55N;
  const ScalarT* state56PtrN = state56N;
  const ScalarT* state57PtrN = state57N;
  const ScalarT* state58PtrN = state58N;
  const ScalarT* state59PtrN = state59N;
  const ScalarT* state60PtrN = state60N;
  const ScalarT* state61PtrN = state61N;
  const ScalarT* state62PtrN = state62N;
  const ScalarT* state63PtrN = state63N;
  const ScalarT* state64PtrN = state64N;
  const ScalarT* state65PtrN = state65N;
  const ScalarT* state66PtrN = state66N;
  const ScalarT* state67PtrN = state67N;
  const ScalarT* state68PtrN = state68N;
  const ScalarT* state69PtrN = state69N;
  const ScalarT* state70PtrN = state70N;
  const ScalarT* state71PtrN = state71N;
  const ScalarT* state72PtrN = state72N;
  const ScalarT* state73PtrN = state73N;
  const ScalarT* state74PtrN = state74N;
  const ScalarT* state75PtrN = state75N;
  const ScalarT* state76PtrN = state76N;
  const ScalarT* state77PtrN = state77N;
  const ScalarT* state78PtrN = state78N;
  const ScalarT* state79PtrN = state79N;
  const ScalarT* state80PtrN = state80N;
  const ScalarT* state81PtrN = state81N;
  const ScalarT* state82PtrN = state82N;
  const ScalarT* state83PtrN = state83N;
  const ScalarT* state84PtrN = state84N;
  const ScalarT* state85PtrN = state85N;
  const ScalarT* state86PtrN = state86N;
  const ScalarT* state87PtrN = state87N;
  const ScalarT* state88PtrN = state88N;
  const ScalarT* state89PtrN = state89N;
  const ScalarT* state90PtrN = state90N;
  const ScalarT* state91PtrN = state91N;
  const ScalarT* state92PtrN = state92N;
  const ScalarT* state93PtrN = state93N;
  const ScalarT* state94PtrN = state94N;
  const ScalarT* state95PtrN = state95N;
  const ScalarT* state96PtrN = state96N;
  const ScalarT* state97PtrN = state97N;
  const ScalarT* state98PtrN = state98N;
  const ScalarT* state99PtrN = state99N;
  const ScalarT* state100PtrN = state100N;
  const ScalarT* state101PtrN = state101N;
  const ScalarT* state102PtrN = state102N;
  const ScalarT* state103PtrN = state103N;
  const ScalarT* state104PtrN = state104N;
  const ScalarT* state105PtrN = state105N;
  const ScalarT* state106PtrN = state106N;
  const ScalarT* state107PtrN = state107N;
  const ScalarT* state108PtrN = state108N;
  const ScalarT* state109PtrN = state109N;
  const ScalarT* state110PtrN = state110N;
  const ScalarT* state111PtrN = state111N;
  const ScalarT* state112PtrN = state112N;
  const ScalarT* state113PtrN = state113N;
  const ScalarT* state114PtrN = state114N;
  const ScalarT* state115PtrN = state115N;
  const ScalarT* state116PtrN = state116N;
  const ScalarT* state117PtrN = state117N;
  const ScalarT* state118PtrN = state118N;
  const ScalarT* state119PtrN = state119N;
  const ScalarT* state120PtrN = state120N;
  const ScalarT* state121PtrN = state121N;
  const ScalarT* state122PtrN = state122N;
  const ScalarT* state123PtrN = state123N;
  const ScalarT* state124PtrN = state124N;
  const ScalarT* state125PtrN = state125N;
  const ScalarT* state126PtrN = state126N;
  const ScalarT* state127PtrN = state127N;
  const ScalarT* state128PtrN = state128N;
  const ScalarT* state129PtrN = state129N;
  const ScalarT* state130PtrN = state130N;
  const ScalarT* state131PtrN = state131N;
  const ScalarT* state132PtrN = state132N;
  const ScalarT* state133PtrN = state133N;
  const ScalarT* state134PtrN = state134N;
  const ScalarT* state135PtrN = state135N;
  const ScalarT* state136PtrN = state136N;
  const ScalarT* state137PtrN = state137N;
  const ScalarT* state138PtrN = state138N;
  const ScalarT* state139PtrN = state139N;
  const ScalarT* state140PtrN = state140N;
  const ScalarT* state141PtrN = state141N;
  const ScalarT* state142PtrN = state142N;
  const ScalarT* state143PtrN = state143N;
  const ScalarT* state144PtrN = state144N;
  const ScalarT* state145PtrN = state145N;
  const ScalarT* state146PtrN = state146N;
  const ScalarT* state147PtrN = state147N;
  const ScalarT* state148PtrN = state148N;
  const ScalarT* state149PtrN = state149N;
  const ScalarT* state150PtrN = state150N;
  const ScalarT* state151PtrN = state151N;
  const ScalarT* state152PtrN = state152N;
  const ScalarT* state153PtrN = state153N;
  const ScalarT* state154PtrN = state154N;
  const ScalarT* state155PtrN = state155N;
  const ScalarT* state156PtrN = state156N;
  const ScalarT* state157PtrN = state157N;
  const ScalarT* state158PtrN = state158N;
  const ScalarT* state159PtrN = state159N;
  const ScalarT* state160PtrN = state160N;
  const ScalarT* state161PtrN = state161N;
  const ScalarT* state162PtrN = state162N;
  const ScalarT* state163PtrN = state163N;
  const ScalarT* state164PtrN = state164N;
  const ScalarT* state165PtrN = state165N;
  const ScalarT* state166PtrN = state166N;
  const ScalarT* state167PtrN = state167N;
  const ScalarT* state168PtrN = state168N;
  const ScalarT* state169PtrN = state169N;
  const ScalarT* state170PtrN = state170N;
  const ScalarT* state171PtrN = state171N;
  const ScalarT* state172PtrN = state172N;
  const ScalarT* state173PtrN = state173N;
  const ScalarT* state174PtrN = state174N;
  const ScalarT* state175PtrN = state175N;
  const ScalarT* state176PtrN = state176N;
  const ScalarT* state177PtrN = state177N;
  const ScalarT* state178PtrN = state178N;
  const ScalarT* state179PtrN = state179N;
  const ScalarT* state180PtrN = state180N;
  const ScalarT* state181PtrN = state181N;
  const ScalarT* state182PtrN = state182N;
  const ScalarT* state183PtrN = state183N;
  const ScalarT* state184PtrN = state184N;
  const ScalarT* state185PtrN = state185N;
  const ScalarT* state186PtrN = state186N;
  const ScalarT* state187PtrN = state187N;
  const ScalarT* state188PtrN = state188N;
  const ScalarT* state189PtrN = state189N;
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
  ScalarT* state50PtrNP1 = state50NP1;
  ScalarT* state51PtrNP1 = state51NP1;
  ScalarT* state52PtrNP1 = state52NP1;
  ScalarT* state53PtrNP1 = state53NP1;
  ScalarT* state54PtrNP1 = state54NP1;
  ScalarT* state55PtrNP1 = state55NP1;
  ScalarT* state56PtrNP1 = state56NP1;
  ScalarT* state57PtrNP1 = state57NP1;
  ScalarT* state58PtrNP1 = state58NP1;
  ScalarT* state59PtrNP1 = state59NP1;
  ScalarT* state60PtrNP1 = state60NP1;
  ScalarT* state61PtrNP1 = state61NP1;
  ScalarT* state62PtrNP1 = state62NP1;
  ScalarT* state63PtrNP1 = state63NP1;
  ScalarT* state64PtrNP1 = state64NP1;
  ScalarT* state65PtrNP1 = state65NP1;
  ScalarT* state66PtrNP1 = state66NP1;
  ScalarT* state67PtrNP1 = state67NP1;
  ScalarT* state68PtrNP1 = state68NP1;
  ScalarT* state69PtrNP1 = state69NP1;
  ScalarT* state70PtrNP1 = state70NP1;
  ScalarT* state71PtrNP1 = state71NP1;
  ScalarT* state72PtrNP1 = state72NP1;
  ScalarT* state73PtrNP1 = state73NP1;
  ScalarT* state74PtrNP1 = state74NP1;
  ScalarT* state75PtrNP1 = state75NP1;
  ScalarT* state76PtrNP1 = state76NP1;
  ScalarT* state77PtrNP1 = state77NP1;
  ScalarT* state78PtrNP1 = state78NP1;
  ScalarT* state79PtrNP1 = state79NP1;
  ScalarT* state80PtrNP1 = state80NP1;
  ScalarT* state81PtrNP1 = state81NP1;
  ScalarT* state82PtrNP1 = state82NP1;
  ScalarT* state83PtrNP1 = state83NP1;
  ScalarT* state84PtrNP1 = state84NP1;
  ScalarT* state85PtrNP1 = state85NP1;
  ScalarT* state86PtrNP1 = state86NP1;
  ScalarT* state87PtrNP1 = state87NP1;
  ScalarT* state88PtrNP1 = state88NP1;
  ScalarT* state89PtrNP1 = state89NP1;
  ScalarT* state90PtrNP1 = state90NP1;
  ScalarT* state91PtrNP1 = state91NP1;
  ScalarT* state92PtrNP1 = state92NP1;
  ScalarT* state93PtrNP1 = state93NP1;
  ScalarT* state94PtrNP1 = state94NP1;
  ScalarT* state95PtrNP1 = state95NP1;
  ScalarT* state96PtrNP1 = state96NP1;
  ScalarT* state97PtrNP1 = state97NP1;
  ScalarT* state98PtrNP1 = state98NP1;
  ScalarT* state99PtrNP1 = state99NP1;
  ScalarT* state100PtrNP1 = state100NP1;
  ScalarT* state101PtrNP1 = state101NP1;
  ScalarT* state102PtrNP1 = state102NP1;
  ScalarT* state103PtrNP1 = state103NP1;
  ScalarT* state104PtrNP1 = state104NP1;
  ScalarT* state105PtrNP1 = state105NP1;
  ScalarT* state106PtrNP1 = state106NP1;
  ScalarT* state107PtrNP1 = state107NP1;
  ScalarT* state108PtrNP1 = state108NP1;
  ScalarT* state109PtrNP1 = state109NP1;
  ScalarT* state110PtrNP1 = state110NP1;
  ScalarT* state111PtrNP1 = state111NP1;
  ScalarT* state112PtrNP1 = state112NP1;
  ScalarT* state113PtrNP1 = state113NP1;
  ScalarT* state114PtrNP1 = state114NP1;
  ScalarT* state115PtrNP1 = state115NP1;
  ScalarT* state116PtrNP1 = state116NP1;
  ScalarT* state117PtrNP1 = state117NP1;
  ScalarT* state118PtrNP1 = state118NP1;
  ScalarT* state119PtrNP1 = state119NP1;
  ScalarT* state120PtrNP1 = state120NP1;
  ScalarT* state121PtrNP1 = state121NP1;
  ScalarT* state122PtrNP1 = state122NP1;
  ScalarT* state123PtrNP1 = state123NP1;
  ScalarT* state124PtrNP1 = state124NP1;
  ScalarT* state125PtrNP1 = state125NP1;
  ScalarT* state126PtrNP1 = state126NP1;
  ScalarT* state127PtrNP1 = state127NP1;
  ScalarT* state128PtrNP1 = state128NP1;
  ScalarT* state129PtrNP1 = state129NP1;
  ScalarT* state130PtrNP1 = state130NP1;
  ScalarT* state131PtrNP1 = state131NP1;
  ScalarT* state132PtrNP1 = state132NP1;
  ScalarT* state133PtrNP1 = state133NP1;
  ScalarT* state134PtrNP1 = state134NP1;
  ScalarT* state135PtrNP1 = state135NP1;
  ScalarT* state136PtrNP1 = state136NP1;
  ScalarT* state137PtrNP1 = state137NP1;
  ScalarT* state138PtrNP1 = state138NP1;
  ScalarT* state139PtrNP1 = state139NP1;
  ScalarT* state140PtrNP1 = state140NP1;
  ScalarT* state141PtrNP1 = state141NP1;
  ScalarT* state142PtrNP1 = state142NP1;
  ScalarT* state143PtrNP1 = state143NP1;
  ScalarT* state144PtrNP1 = state144NP1;
  ScalarT* state145PtrNP1 = state145NP1;
  ScalarT* state146PtrNP1 = state146NP1;
  ScalarT* state147PtrNP1 = state147NP1;
  ScalarT* state148PtrNP1 = state148NP1;
  ScalarT* state149PtrNP1 = state149NP1;
  ScalarT* state150PtrNP1 = state150NP1;
  ScalarT* state151PtrNP1 = state151NP1;
  ScalarT* state152PtrNP1 = state152NP1;
  ScalarT* state153PtrNP1 = state153NP1;
  ScalarT* state154PtrNP1 = state154NP1;
  ScalarT* state155PtrNP1 = state155NP1;
  ScalarT* state156PtrNP1 = state156NP1;
  ScalarT* state157PtrNP1 = state157NP1;
  ScalarT* state158PtrNP1 = state158NP1;
  ScalarT* state159PtrNP1 = state159NP1;
  ScalarT* state160PtrNP1 = state160NP1;
  ScalarT* state161PtrNP1 = state161NP1;
  ScalarT* state162PtrNP1 = state162NP1;
  ScalarT* state163PtrNP1 = state163NP1;
  ScalarT* state164PtrNP1 = state164NP1;
  ScalarT* state165PtrNP1 = state165NP1;
  ScalarT* state166PtrNP1 = state166NP1;
  ScalarT* state167PtrNP1 = state167NP1;
  ScalarT* state168PtrNP1 = state168NP1;
  ScalarT* state169PtrNP1 = state169NP1;
  ScalarT* state170PtrNP1 = state170NP1;
  ScalarT* state171PtrNP1 = state171NP1;
  ScalarT* state172PtrNP1 = state172NP1;
  ScalarT* state173PtrNP1 = state173NP1;
  ScalarT* state174PtrNP1 = state174NP1;
  ScalarT* state175PtrNP1 = state175NP1;
  ScalarT* state176PtrNP1 = state176NP1;
  ScalarT* state177PtrNP1 = state177NP1;
  ScalarT* state178PtrNP1 = state178NP1;
  ScalarT* state179PtrNP1 = state179NP1;
  ScalarT* state180PtrNP1 = state180NP1;
  ScalarT* state181PtrNP1 = state181NP1;
  ScalarT* state182PtrNP1 = state182NP1;
  ScalarT* state183PtrNP1 = state183NP1;
  ScalarT* state184PtrNP1 = state184NP1;
  ScalarT* state185PtrNP1 = state185NP1;
  ScalarT* state186PtrNP1 = state186NP1;
  ScalarT* state187PtrNP1 = state187NP1;
  ScalarT* state188PtrNP1 = state188NP1;
  ScalarT* state189PtrNP1 = state189NP1;
  ScalarT* internEnergyNP1 = internalEnergyNP1;
  ScalarT* inelasEnergyNP1 = inelasticEnergyNP1;

  ScalarT deps[6];
  ScalarT epsN[6];

  ScalarT sigmaN[6];
  ScalarT sigmaNP1[6];

  ScalarT stateN[189];
  ScalarT stateNP1[189];

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
      ++state46PtrN, ++state47PtrN, ++state48PtrN, ++state49PtrN, ++state50PtrN,
      ++state51PtrN, ++state52PtrN, ++state53PtrN, ++state54PtrN, ++state55PtrN,
      ++state56PtrN, ++state57PtrN, ++state58PtrN, ++state59PtrN, ++state60PtrN,
      ++state61PtrN, ++state62PtrN, ++state63PtrN, ++state64PtrN, ++state65PtrN,
      ++state66PtrN, ++state67PtrN, ++state68PtrN, ++state69PtrN, ++state70PtrN,
      ++state71PtrN, ++state72PtrN, ++state73PtrN, ++state74PtrN, ++state75PtrN,
      ++state76PtrN, ++state77PtrN, ++state78PtrN, ++state79PtrN, ++state80PtrN,
      ++state81PtrN, ++state82PtrN, ++state83PtrN, ++state84PtrN, ++state85PtrN,
      ++state86PtrN, ++state87PtrN, ++state88PtrN, ++state89PtrN, ++state90PtrN,
      ++state91PtrN, ++state92PtrN, ++state93PtrN, ++state94PtrN, ++state95PtrN,
      ++state96PtrN, ++state97PtrN, ++state98PtrN, ++state99PtrN, ++state100PtrN,
      ++state101PtrN, ++state102PtrN, ++state103PtrN, ++state104PtrN, ++state105PtrN,
      ++state106PtrN, ++state107PtrN, ++state108PtrN, ++state109PtrN, ++state110PtrN,
      ++state111PtrN, ++state112PtrN, ++state113PtrN, ++state114PtrN, ++state115PtrN,
      ++state116PtrN, ++state117PtrN, ++state118PtrN, ++state119PtrN, ++state120PtrN,
      ++state121PtrN, ++state122PtrN, ++state123PtrN, ++state124PtrN, ++state125PtrN,
      ++state126PtrN, ++state127PtrN, ++state128PtrN, ++state129PtrN, ++state130PtrN,
      ++state131PtrN, ++state132PtrN, ++state133PtrN, ++state134PtrN, ++state135PtrN,
      ++state136PtrN, ++state137PtrN, ++state138PtrN, ++state139PtrN, ++state140PtrN,
      ++state141PtrN, ++state142PtrN, ++state143PtrN, ++state144PtrN, ++state145PtrN,
      ++state146PtrN, ++state147PtrN, ++state148PtrN, ++state149PtrN, ++state150PtrN,
      ++state151PtrN, ++state152PtrN, ++state153PtrN, ++state154PtrN, ++state155PtrN,
      ++state156PtrN, ++state157PtrN, ++state158PtrN, ++state159PtrN, ++state160PtrN,
      ++state161PtrN, ++state162PtrN, ++state163PtrN, ++state164PtrN, ++state165PtrN,
      ++state166PtrN, ++state167PtrN, ++state168PtrN, ++state169PtrN, ++state170PtrN,
      ++state171PtrN, ++state172PtrN, ++state173PtrN, ++state174PtrN, ++state175PtrN,
      ++state176PtrN, ++state177PtrN, ++state178PtrN, ++state179PtrN, ++state180PtrN,
      ++state181PtrN, ++state182PtrN, ++state183PtrN, ++state184PtrN, ++state185PtrN,
      ++state186PtrN, ++state187PtrN, ++state188PtrN, ++state189PtrN, 
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
      ++state46PtrNP1, ++state47PtrNP1, ++state48PtrNP1, ++state49PtrNP1, ++state50PtrNP1,
      ++state51PtrNP1, ++state52PtrNP1, ++state53PtrNP1, ++state54PtrNP1, ++state55PtrNP1,
      ++state56PtrNP1, ++state57PtrNP1, ++state58PtrNP1, ++state59PtrNP1, ++state60PtrNP1,
      ++state61PtrNP1, ++state62PtrNP1, ++state63PtrNP1, ++state64PtrNP1, ++state65PtrNP1,
      ++state66PtrNP1, ++state67PtrNP1, ++state68PtrNP1, ++state69PtrNP1, ++state70PtrNP1,
      ++state71PtrNP1, ++state72PtrNP1, ++state73PtrNP1, ++state74PtrNP1, ++state75PtrNP1,
      ++state76PtrNP1, ++state77PtrNP1, ++state78PtrNP1, ++state79PtrNP1, ++state80PtrNP1,
      ++state81PtrNP1, ++state82PtrNP1, ++state83PtrNP1, ++state84PtrNP1, ++state85PtrNP1,
      ++state86PtrNP1, ++state87PtrNP1, ++state88PtrNP1, ++state89PtrNP1, ++state90PtrNP1,
      ++state91PtrNP1, ++state92PtrNP1, ++state93PtrNP1, ++state94PtrNP1, ++state95PtrNP1,
      ++state96PtrNP1, ++state97PtrNP1, ++state98PtrNP1, ++state99PtrNP1, ++state100PtrNP1,
      ++state101PtrNP1, ++state102PtrNP1, ++state103PtrNP1, ++state104PtrNP1, ++state105PtrNP1,
      ++state106PtrNP1, ++state107PtrNP1, ++state108PtrNP1, ++state109PtrNP1, ++state110PtrNP1,
      ++state111PtrNP1, ++state112PtrNP1, ++state113PtrNP1, ++state114PtrNP1, ++state115PtrNP1,
      ++state116PtrNP1, ++state117PtrNP1, ++state118PtrNP1, ++state119PtrNP1, ++state120PtrNP1,
      ++state121PtrNP1, ++state122PtrNP1, ++state123PtrNP1, ++state124PtrNP1, ++state125PtrNP1,
      ++state126PtrNP1, ++state127PtrNP1, ++state128PtrNP1, ++state129PtrNP1, ++state130PtrNP1,
      ++state131PtrNP1, ++state132PtrNP1, ++state133PtrNP1, ++state134PtrNP1, ++state135PtrNP1,
      ++state136PtrNP1, ++state137PtrNP1, ++state138PtrNP1, ++state139PtrNP1, ++state140PtrNP1,
      ++state141PtrNP1, ++state142PtrNP1, ++state143PtrNP1, ++state144PtrNP1, ++state145PtrNP1,
      ++state146PtrNP1, ++state147PtrNP1, ++state148PtrNP1, ++state149PtrNP1, ++state150PtrNP1,
      ++state151PtrNP1, ++state152PtrNP1, ++state153PtrNP1, ++state154PtrNP1, ++state155PtrNP1,
      ++state156PtrNP1, ++state157PtrNP1, ++state158PtrNP1, ++state159PtrNP1, ++state160PtrNP1,
      ++state161PtrNP1, ++state162PtrNP1, ++state163PtrNP1, ++state164PtrNP1, ++state165PtrNP1,
      ++state166PtrNP1, ++state167PtrNP1, ++state168PtrNP1, ++state169PtrNP1, ++state170PtrNP1,
      ++state171PtrNP1, ++state172PtrNP1, ++state173PtrNP1, ++state174PtrNP1, ++state175PtrNP1,
      ++state176PtrNP1, ++state177PtrNP1, ++state178PtrNP1, ++state179PtrNP1, ++state180PtrNP1,
      ++state181PtrNP1, ++state182PtrNP1, ++state183PtrNP1, ++state184PtrNP1, ++state185PtrNP1,
      ++state186PtrNP1, ++state187PtrNP1, ++state188PtrNP1, ++state189PtrNP1
      ){

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
        deps[3] = *(Edot+5) * dt;
        deps[4] = *(Edot+2) * dt;
        deps[5] = *(Edot+1) * dt;

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
        stateN[49] = *state50PtrN;
        stateN[50] = *state51PtrN;
        stateN[51] = *state52PtrN;
        stateN[52] = *state53PtrN;
        stateN[53] = *state54PtrN;
        stateN[54] = *state55PtrN;
        stateN[55] = *state56PtrN;
        stateN[56] = *state57PtrN;
        stateN[57] = *state58PtrN;
        stateN[58] = *state59PtrN;
        stateN[59] = *state60PtrN;
        stateN[60] = *state61PtrN;
        stateN[61] = *state62PtrN;
        stateN[62] = *state63PtrN;
        stateN[63] = *state64PtrN;
        stateN[64] = *state65PtrN;
        stateN[65] = *state66PtrN;
        stateN[66] = *state67PtrN;
        stateN[67] = *state68PtrN;
        stateN[68] = *state69PtrN;
        stateN[69] = *state70PtrN;
        stateN[70] = *state71PtrN;
        stateN[71] = *state72PtrN;
        stateN[72] = *state73PtrN;
        stateN[73] = *state74PtrN;
        stateN[74] = *state75PtrN;
        stateN[75] = *state76PtrN;
        stateN[76] = *state77PtrN;
        stateN[77] = *state78PtrN;
        stateN[78] = *state79PtrN;
        stateN[79] = *state80PtrN;
        stateN[80] = *state81PtrN;
        stateN[81] = *state82PtrN;
        stateN[82] = *state83PtrN;
        stateN[83] = *state84PtrN;
        stateN[84] = *state85PtrN;
        stateN[85] = *state86PtrN;
        stateN[86] = *state87PtrN;
        stateN[87] = *state88PtrN;
        stateN[88] = *state89PtrN;
        stateN[89] = *state90PtrN;
        stateN[90] = *state91PtrN;
        stateN[91] = *state92PtrN;
        stateN[92] = *state93PtrN;
        stateN[93] = *state94PtrN;
        stateN[94] = *state95PtrN;
        stateN[95] = *state96PtrN;
        stateN[96] = *state97PtrN;
        stateN[97] = *state98PtrN;
        stateN[98] = *state99PtrN;
        stateN[99] = *state100PtrN;
        stateN[100] = *state101PtrN;
        stateN[101] = *state102PtrN;
        stateN[102] = *state103PtrN;
        stateN[103] = *state104PtrN;
        stateN[104] = *state105PtrN;
        stateN[105] = *state106PtrN;
        stateN[106] = *state107PtrN;
        stateN[107] = *state108PtrN;
        stateN[108] = *state109PtrN;
        stateN[109] = *state110PtrN;
        stateN[110] = *state111PtrN;
        stateN[111] = *state112PtrN;
        stateN[112] = *state113PtrN;
        stateN[113] = *state114PtrN;
        stateN[114] = *state115PtrN;
        stateN[115] = *state116PtrN;
        stateN[116] = *state117PtrN;
        stateN[117] = *state118PtrN;
        stateN[118] = *state119PtrN;
        stateN[119] = *state120PtrN;
        stateN[120] = *state121PtrN;
        stateN[121] = *state122PtrN;
        stateN[122] = *state123PtrN;
        stateN[123] = *state124PtrN;
        stateN[124] = *state125PtrN;
        stateN[125] = *state126PtrN;
        stateN[126] = *state127PtrN;
        stateN[127] = *state128PtrN;
        stateN[128] = *state129PtrN;
        stateN[129] = *state130PtrN;
        stateN[130] = *state131PtrN;
        stateN[131] = *state132PtrN;
        stateN[132] = *state133PtrN;
        stateN[133] = *state134PtrN;
        stateN[134] = *state135PtrN;
        stateN[135] = *state136PtrN;
        stateN[136] = *state137PtrN;
        stateN[137] = *state138PtrN;
        stateN[138] = *state139PtrN;
        stateN[139] = *state140PtrN;
        stateN[140] = *state141PtrN;
        stateN[141] = *state142PtrN;
        stateN[142] = *state143PtrN;
        stateN[143] = *state144PtrN;
        stateN[144] = *state145PtrN;
        stateN[145] = *state146PtrN;
        stateN[146] = *state147PtrN;
        stateN[147] = *state148PtrN;
        stateN[148] = *state149PtrN;
        stateN[149] = *state150PtrN;
        stateN[150] = *state151PtrN;
        stateN[151] = *state152PtrN;
        stateN[152] = *state153PtrN;
        stateN[153] = *state154PtrN;
        stateN[154] = *state155PtrN;
        stateN[155] = *state156PtrN;
        stateN[156] = *state157PtrN;
        stateN[157] = *state158PtrN;
        stateN[158] = *state159PtrN;
        stateN[159] = *state160PtrN;
        stateN[160] = *state161PtrN;
        stateN[161] = *state162PtrN;
        stateN[162] = *state163PtrN;
        stateN[163] = *state164PtrN;
        stateN[164] = *state165PtrN;
        stateN[165] = *state166PtrN;
        stateN[166] = *state167PtrN;
        stateN[167] = *state168PtrN;
        stateN[168] = *state169PtrN;
        stateN[169] = *state170PtrN;
        stateN[170] = *state171PtrN;
        stateN[171] = *state172PtrN;
        stateN[172] = *state173PtrN;
        stateN[173] = *state174PtrN;
        stateN[174] = *state175PtrN;
        stateN[175] = *state176PtrN;
        stateN[176] = *state177PtrN;
        stateN[177] = *state178PtrN;
        stateN[178] = *state179PtrN;
        stateN[179] = *state180PtrN;
        stateN[180] = *state181PtrN;
        stateN[181] = *state182PtrN;
        stateN[182] = *state183PtrN;
        stateN[183] = *state184PtrN;
        stateN[184] = *state185PtrN;
        stateN[185] = *state186PtrN;
        stateN[186] = *state187PtrN;
        stateN[187] = *state188PtrN;
        stateN[188] = *state189PtrN;

        // call the m7 function
        m7fmaterial_(&dt, 
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
        *state50PtrNP1 = stateNP1[49];
        *state51PtrNP1 = stateNP1[50];
        *state52PtrNP1 = stateNP1[51];
        *state53PtrNP1 = stateNP1[52];
        *state54PtrNP1 = stateNP1[53];
        *state55PtrNP1 = stateNP1[54];
        *state56PtrNP1 = stateNP1[55];
        *state57PtrNP1 = stateNP1[56];
        *state58PtrNP1 = stateNP1[57];
        *state59PtrNP1 = stateNP1[58];
        *state60PtrNP1 = stateNP1[59];
        *state61PtrNP1 = stateNP1[60];
        *state62PtrNP1 = stateNP1[61];
        *state63PtrNP1 = stateNP1[62];
        *state64PtrNP1 = stateNP1[63];
        *state65PtrNP1 = stateNP1[64];
        *state66PtrNP1 = stateNP1[65];
        *state67PtrNP1 = stateNP1[66];
        *state68PtrNP1 = stateNP1[67];
        *state69PtrNP1 = stateNP1[68];
        *state70PtrNP1 = stateNP1[69];
        *state71PtrNP1 = stateNP1[70];
        *state72PtrNP1 = stateNP1[71];
        *state73PtrNP1 = stateNP1[72];
        *state74PtrNP1 = stateNP1[73];
        *state75PtrNP1 = stateNP1[74];
        *state76PtrNP1 = stateNP1[75];
        *state77PtrNP1 = stateNP1[76];
        *state78PtrNP1 = stateNP1[77];
        *state79PtrNP1 = stateNP1[78];
        *state80PtrNP1 = stateNP1[79];
        *state81PtrNP1 = stateNP1[80];
        *state82PtrNP1 = stateNP1[81];
        *state83PtrNP1 = stateNP1[82];
        *state84PtrNP1 = stateNP1[83];
        *state85PtrNP1 = stateNP1[84];
        *state86PtrNP1 = stateNP1[85];
        *state87PtrNP1 = stateNP1[86];
        *state88PtrNP1 = stateNP1[87];
        *state89PtrNP1 = stateNP1[88];
        *state90PtrNP1 = stateNP1[89];
        *state91PtrNP1 = stateNP1[90];
        *state92PtrNP1 = stateNP1[91];
        *state93PtrNP1 = stateNP1[92];
        *state94PtrNP1 = stateNP1[93];
        *state95PtrNP1 = stateNP1[94];
        *state96PtrNP1 = stateNP1[95];
        *state97PtrNP1 = stateNP1[96];
        *state98PtrNP1 = stateNP1[97];
        *state99PtrNP1 = stateNP1[98];
        *state100PtrNP1 = stateNP1[99];
        *state101PtrNP1 = stateNP1[100];
        *state102PtrNP1 = stateNP1[101];
        *state103PtrNP1 = stateNP1[102];
        *state104PtrNP1 = stateNP1[103];
        *state105PtrNP1 = stateNP1[104];
        *state106PtrNP1 = stateNP1[105];
        *state107PtrNP1 = stateNP1[106];
        *state108PtrNP1 = stateNP1[107];
        *state109PtrNP1 = stateNP1[108];
        *state110PtrNP1 = stateNP1[109];
        *state111PtrNP1 = stateNP1[110];
        *state112PtrNP1 = stateNP1[111];
        *state113PtrNP1 = stateNP1[112];
        *state114PtrNP1 = stateNP1[113];
        *state115PtrNP1 = stateNP1[114];
        *state116PtrNP1 = stateNP1[115];
        *state117PtrNP1 = stateNP1[116];
        *state118PtrNP1 = stateNP1[117];
        *state119PtrNP1 = stateNP1[118];
        *state120PtrNP1 = stateNP1[119];
        *state121PtrNP1 = stateNP1[120];
        *state122PtrNP1 = stateNP1[121];
        *state123PtrNP1 = stateNP1[122];
        *state124PtrNP1 = stateNP1[123];
        *state125PtrNP1 = stateNP1[124];
        *state126PtrNP1 = stateNP1[125];
        *state127PtrNP1 = stateNP1[126];
        *state128PtrNP1 = stateNP1[127];
        *state129PtrNP1 = stateNP1[128];
        *state130PtrNP1 = stateNP1[129];
        *state131PtrNP1 = stateNP1[130];
        *state132PtrNP1 = stateNP1[131];
        *state133PtrNP1 = stateNP1[132];
        *state134PtrNP1 = stateNP1[133];
        *state135PtrNP1 = stateNP1[134];
        *state136PtrNP1 = stateNP1[135];
        *state137PtrNP1 = stateNP1[136];
        *state138PtrNP1 = stateNP1[137];
        *state139PtrNP1 = stateNP1[138];
        *state140PtrNP1 = stateNP1[139];
        *state141PtrNP1 = stateNP1[140];
        *state142PtrNP1 = stateNP1[141];
        *state143PtrNP1 = stateNP1[142];
        *state144PtrNP1 = stateNP1[143];
        *state145PtrNP1 = stateNP1[144];
        *state146PtrNP1 = stateNP1[145];
        *state147PtrNP1 = stateNP1[146];
        *state148PtrNP1 = stateNP1[147];
        *state149PtrNP1 = stateNP1[148];
        *state150PtrNP1 = stateNP1[149];
        *state151PtrNP1 = stateNP1[150];
        *state152PtrNP1 = stateNP1[151];
        *state153PtrNP1 = stateNP1[152];
        *state154PtrNP1 = stateNP1[153];
        *state155PtrNP1 = stateNP1[154];
        *state156PtrNP1 = stateNP1[155];
        *state157PtrNP1 = stateNP1[156];
        *state158PtrNP1 = stateNP1[157];
        *state159PtrNP1 = stateNP1[158];
        *state160PtrNP1 = stateNP1[159];
        *state161PtrNP1 = stateNP1[160];
        *state162PtrNP1 = stateNP1[161];
        *state163PtrNP1 = stateNP1[162];
        *state164PtrNP1 = stateNP1[163];
        *state165PtrNP1 = stateNP1[164];
        *state166PtrNP1 = stateNP1[165];
        *state167PtrNP1 = stateNP1[166];
        *state168PtrNP1 = stateNP1[167];
        *state169PtrNP1 = stateNP1[168];
        *state170PtrNP1 = stateNP1[169];
        *state171PtrNP1 = stateNP1[170];
        *state172PtrNP1 = stateNP1[171];
        *state173PtrNP1 = stateNP1[172];
        *state174PtrNP1 = stateNP1[173];
        *state175PtrNP1 = stateNP1[174];
        *state176PtrNP1 = stateNP1[175];
        *state177PtrNP1 = stateNP1[176];
        *state178PtrNP1 = stateNP1[177];
        *state179PtrNP1 = stateNP1[178];
        *state180PtrNP1 = stateNP1[179];
        *state181PtrNP1 = stateNP1[180];
        *state182PtrNP1 = stateNP1[181];
        *state183PtrNP1 = stateNP1[182];
        *state184PtrNP1 = stateNP1[183];
        *state185PtrNP1 = stateNP1[184];
        *state186PtrNP1 = stateNP1[185];
        *state187PtrNP1 = stateNP1[186];
        *state188PtrNP1 = stateNP1[187];
        *state189PtrNP1 = stateNP1[188];
  }
}

template<typename ScalarT>
void updateBondLevelMicroplaneM7Stress
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
    const ScalarT* bondLevelState50N,
    const ScalarT* bondLevelState51N,
    const ScalarT* bondLevelState52N,
    const ScalarT* bondLevelState53N,
    const ScalarT* bondLevelState54N,
    const ScalarT* bondLevelState55N,
    const ScalarT* bondLevelState56N,
    const ScalarT* bondLevelState57N,
    const ScalarT* bondLevelState58N,
    const ScalarT* bondLevelState59N,
    const ScalarT* bondLevelState60N,
    const ScalarT* bondLevelState61N,
    const ScalarT* bondLevelState62N,
    const ScalarT* bondLevelState63N,
    const ScalarT* bondLevelState64N,
    const ScalarT* bondLevelState65N,
    const ScalarT* bondLevelState66N,
    const ScalarT* bondLevelState67N,
    const ScalarT* bondLevelState68N,
    const ScalarT* bondLevelState69N,
    const ScalarT* bondLevelState70N,
    const ScalarT* bondLevelState71N,
    const ScalarT* bondLevelState72N,
    const ScalarT* bondLevelState73N,
    const ScalarT* bondLevelState74N,
    const ScalarT* bondLevelState75N,
    const ScalarT* bondLevelState76N,
    const ScalarT* bondLevelState77N,
    const ScalarT* bondLevelState78N,
    const ScalarT* bondLevelState79N,
    const ScalarT* bondLevelState80N,
    const ScalarT* bondLevelState81N,
    const ScalarT* bondLevelState82N,
    const ScalarT* bondLevelState83N,
    const ScalarT* bondLevelState84N,
    const ScalarT* bondLevelState85N,
    const ScalarT* bondLevelState86N,
    const ScalarT* bondLevelState87N,
    const ScalarT* bondLevelState88N,
    const ScalarT* bondLevelState89N,
    const ScalarT* bondLevelState90N,
    const ScalarT* bondLevelState91N,
    const ScalarT* bondLevelState92N,
    const ScalarT* bondLevelState93N,
    const ScalarT* bondLevelState94N,
    const ScalarT* bondLevelState95N,
    const ScalarT* bondLevelState96N,
    const ScalarT* bondLevelState97N,
    const ScalarT* bondLevelState98N,
    const ScalarT* bondLevelState99N,
    const ScalarT* bondLevelState100N,
    const ScalarT* bondLevelState101N,
    const ScalarT* bondLevelState102N,
    const ScalarT* bondLevelState103N,
    const ScalarT* bondLevelState104N,
    const ScalarT* bondLevelState105N,
    const ScalarT* bondLevelState106N,
    const ScalarT* bondLevelState107N,
    const ScalarT* bondLevelState108N,
    const ScalarT* bondLevelState109N,
    const ScalarT* bondLevelState110N,
    const ScalarT* bondLevelState111N,
    const ScalarT* bondLevelState112N,
    const ScalarT* bondLevelState113N,
    const ScalarT* bondLevelState114N,
    const ScalarT* bondLevelState115N,
    const ScalarT* bondLevelState116N,
    const ScalarT* bondLevelState117N,
    const ScalarT* bondLevelState118N,
    const ScalarT* bondLevelState119N,
    const ScalarT* bondLevelState120N,
    const ScalarT* bondLevelState121N,
    const ScalarT* bondLevelState122N,
    const ScalarT* bondLevelState123N,
    const ScalarT* bondLevelState124N,
    const ScalarT* bondLevelState125N,
    const ScalarT* bondLevelState126N,
    const ScalarT* bondLevelState127N,
    const ScalarT* bondLevelState128N,
    const ScalarT* bondLevelState129N,
    const ScalarT* bondLevelState130N,
    const ScalarT* bondLevelState131N,
    const ScalarT* bondLevelState132N,
    const ScalarT* bondLevelState133N,
    const ScalarT* bondLevelState134N,
    const ScalarT* bondLevelState135N,
    const ScalarT* bondLevelState136N,
    const ScalarT* bondLevelState137N,
    const ScalarT* bondLevelState138N,
    const ScalarT* bondLevelState139N,
    const ScalarT* bondLevelState140N,
    const ScalarT* bondLevelState141N,
    const ScalarT* bondLevelState142N,
    const ScalarT* bondLevelState143N,
    const ScalarT* bondLevelState144N,
    const ScalarT* bondLevelState145N,
    const ScalarT* bondLevelState146N,
    const ScalarT* bondLevelState147N,
    const ScalarT* bondLevelState148N,
    const ScalarT* bondLevelState149N,
    const ScalarT* bondLevelState150N,
    const ScalarT* bondLevelState151N,
    const ScalarT* bondLevelState152N,
    const ScalarT* bondLevelState153N,
    const ScalarT* bondLevelState154N,
    const ScalarT* bondLevelState155N,
    const ScalarT* bondLevelState156N,
    const ScalarT* bondLevelState157N,
    const ScalarT* bondLevelState158N,
    const ScalarT* bondLevelState159N,
    const ScalarT* bondLevelState160N,
    const ScalarT* bondLevelState161N,
    const ScalarT* bondLevelState162N,
    const ScalarT* bondLevelState163N,
    const ScalarT* bondLevelState164N,
    const ScalarT* bondLevelState165N,
    const ScalarT* bondLevelState166N,
    const ScalarT* bondLevelState167N,
    const ScalarT* bondLevelState168N,
    const ScalarT* bondLevelState169N,
    const ScalarT* bondLevelState170N,
    const ScalarT* bondLevelState171N,
    const ScalarT* bondLevelState172N,
    const ScalarT* bondLevelState173N,
    const ScalarT* bondLevelState174N,
    const ScalarT* bondLevelState175N,
    const ScalarT* bondLevelState176N,
    const ScalarT* bondLevelState177N,
    const ScalarT* bondLevelState178N,
    const ScalarT* bondLevelState179N,
    const ScalarT* bondLevelState180N,
    const ScalarT* bondLevelState181N,
    const ScalarT* bondLevelState182N,
    const ScalarT* bondLevelState183N,
    const ScalarT* bondLevelState184N,
    const ScalarT* bondLevelState185N,
    const ScalarT* bondLevelState186N,
    const ScalarT* bondLevelState187N,
    const ScalarT* bondLevelState188N,
    const ScalarT* bondLevelState189N,
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
    ScalarT* bondLevelState50NP1,
    ScalarT* bondLevelState51NP1,
    ScalarT* bondLevelState52NP1,
    ScalarT* bondLevelState53NP1,
    ScalarT* bondLevelState54NP1,
    ScalarT* bondLevelState55NP1,
    ScalarT* bondLevelState56NP1,
    ScalarT* bondLevelState57NP1,
    ScalarT* bondLevelState58NP1,
    ScalarT* bondLevelState59NP1,
    ScalarT* bondLevelState60NP1,
    ScalarT* bondLevelState61NP1,
    ScalarT* bondLevelState62NP1,
    ScalarT* bondLevelState63NP1,
    ScalarT* bondLevelState64NP1,
    ScalarT* bondLevelState65NP1,
    ScalarT* bondLevelState66NP1,
    ScalarT* bondLevelState67NP1,
    ScalarT* bondLevelState68NP1,
    ScalarT* bondLevelState69NP1,
    ScalarT* bondLevelState70NP1,
    ScalarT* bondLevelState71NP1,
    ScalarT* bondLevelState72NP1,
    ScalarT* bondLevelState73NP1,
    ScalarT* bondLevelState74NP1,
    ScalarT* bondLevelState75NP1,
    ScalarT* bondLevelState76NP1,
    ScalarT* bondLevelState77NP1,
    ScalarT* bondLevelState78NP1,
    ScalarT* bondLevelState79NP1,
    ScalarT* bondLevelState80NP1,
    ScalarT* bondLevelState81NP1,
    ScalarT* bondLevelState82NP1,
    ScalarT* bondLevelState83NP1,
    ScalarT* bondLevelState84NP1,
    ScalarT* bondLevelState85NP1,
    ScalarT* bondLevelState86NP1,
    ScalarT* bondLevelState87NP1,
    ScalarT* bondLevelState88NP1,
    ScalarT* bondLevelState89NP1,
    ScalarT* bondLevelState90NP1,
    ScalarT* bondLevelState91NP1,
    ScalarT* bondLevelState92NP1,
    ScalarT* bondLevelState93NP1,
    ScalarT* bondLevelState94NP1,
    ScalarT* bondLevelState95NP1,
    ScalarT* bondLevelState96NP1,
    ScalarT* bondLevelState97NP1,
    ScalarT* bondLevelState98NP1,
    ScalarT* bondLevelState99NP1,
    ScalarT* bondLevelState100NP1,
    ScalarT* bondLevelState101NP1,
    ScalarT* bondLevelState102NP1,
    ScalarT* bondLevelState103NP1,
    ScalarT* bondLevelState104NP1,
    ScalarT* bondLevelState105NP1,
    ScalarT* bondLevelState106NP1,
    ScalarT* bondLevelState107NP1,
    ScalarT* bondLevelState108NP1,
    ScalarT* bondLevelState109NP1,
    ScalarT* bondLevelState110NP1,
    ScalarT* bondLevelState111NP1,
    ScalarT* bondLevelState112NP1,
    ScalarT* bondLevelState113NP1,
    ScalarT* bondLevelState114NP1,
    ScalarT* bondLevelState115NP1,
    ScalarT* bondLevelState116NP1,
    ScalarT* bondLevelState117NP1,
    ScalarT* bondLevelState118NP1,
    ScalarT* bondLevelState119NP1,
    ScalarT* bondLevelState120NP1,
    ScalarT* bondLevelState121NP1,
    ScalarT* bondLevelState122NP1,
    ScalarT* bondLevelState123NP1,
    ScalarT* bondLevelState124NP1,
    ScalarT* bondLevelState125NP1,
    ScalarT* bondLevelState126NP1,
    ScalarT* bondLevelState127NP1,
    ScalarT* bondLevelState128NP1,
    ScalarT* bondLevelState129NP1,
    ScalarT* bondLevelState130NP1,
    ScalarT* bondLevelState131NP1,
    ScalarT* bondLevelState132NP1,
    ScalarT* bondLevelState133NP1,
    ScalarT* bondLevelState134NP1,
    ScalarT* bondLevelState135NP1,
    ScalarT* bondLevelState136NP1,
    ScalarT* bondLevelState137NP1,
    ScalarT* bondLevelState138NP1,
    ScalarT* bondLevelState139NP1,
    ScalarT* bondLevelState140NP1,
    ScalarT* bondLevelState141NP1,
    ScalarT* bondLevelState142NP1,
    ScalarT* bondLevelState143NP1,
    ScalarT* bondLevelState144NP1,
    ScalarT* bondLevelState145NP1,
    ScalarT* bondLevelState146NP1,
    ScalarT* bondLevelState147NP1,
    ScalarT* bondLevelState148NP1,
    ScalarT* bondLevelState149NP1,
    ScalarT* bondLevelState150NP1,
    ScalarT* bondLevelState151NP1,
    ScalarT* bondLevelState152NP1,
    ScalarT* bondLevelState153NP1,
    ScalarT* bondLevelState154NP1,
    ScalarT* bondLevelState155NP1,
    ScalarT* bondLevelState156NP1,
    ScalarT* bondLevelState157NP1,
    ScalarT* bondLevelState158NP1,
    ScalarT* bondLevelState159NP1,
    ScalarT* bondLevelState160NP1,
    ScalarT* bondLevelState161NP1,
    ScalarT* bondLevelState162NP1,
    ScalarT* bondLevelState163NP1,
    ScalarT* bondLevelState164NP1,
    ScalarT* bondLevelState165NP1,
    ScalarT* bondLevelState166NP1,
    ScalarT* bondLevelState167NP1,
    ScalarT* bondLevelState168NP1,
    ScalarT* bondLevelState169NP1,
    ScalarT* bondLevelState170NP1,
    ScalarT* bondLevelState171NP1,
    ScalarT* bondLevelState172NP1,
    ScalarT* bondLevelState173NP1,
    ScalarT* bondLevelState174NP1,
    ScalarT* bondLevelState175NP1,
    ScalarT* bondLevelState176NP1,
    ScalarT* bondLevelState177NP1,
    ScalarT* bondLevelState178NP1,
    ScalarT* bondLevelState179NP1,
    ScalarT* bondLevelState180NP1,
    ScalarT* bondLevelState181NP1,
    ScalarT* bondLevelState182NP1,
    ScalarT* bondLevelState183NP1,
    ScalarT* bondLevelState184NP1,
    ScalarT* bondLevelState185NP1,
    ScalarT* bondLevelState186NP1,
    ScalarT* bondLevelState187NP1,
    ScalarT* bondLevelState188NP1,
    ScalarT* bondLevelState189NP1,
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
  const ScalarT* state50PtrN = bondLevelState50N;
  const ScalarT* state51PtrN = bondLevelState51N;
  const ScalarT* state52PtrN = bondLevelState52N;
  const ScalarT* state53PtrN = bondLevelState53N;
  const ScalarT* state54PtrN = bondLevelState54N;
  const ScalarT* state55PtrN = bondLevelState55N;
  const ScalarT* state56PtrN = bondLevelState56N;
  const ScalarT* state57PtrN = bondLevelState57N;
  const ScalarT* state58PtrN = bondLevelState58N;
  const ScalarT* state59PtrN = bondLevelState59N;
  const ScalarT* state60PtrN = bondLevelState60N;
  const ScalarT* state61PtrN = bondLevelState61N;
  const ScalarT* state62PtrN = bondLevelState62N;
  const ScalarT* state63PtrN = bondLevelState63N;
  const ScalarT* state64PtrN = bondLevelState64N;
  const ScalarT* state65PtrN = bondLevelState65N;
  const ScalarT* state66PtrN = bondLevelState66N;
  const ScalarT* state67PtrN = bondLevelState67N;
  const ScalarT* state68PtrN = bondLevelState68N;
  const ScalarT* state69PtrN = bondLevelState69N;
  const ScalarT* state70PtrN = bondLevelState70N;
  const ScalarT* state71PtrN = bondLevelState71N;
  const ScalarT* state72PtrN = bondLevelState72N;
  const ScalarT* state73PtrN = bondLevelState73N;
  const ScalarT* state74PtrN = bondLevelState74N;
  const ScalarT* state75PtrN = bondLevelState75N;
  const ScalarT* state76PtrN = bondLevelState76N;
  const ScalarT* state77PtrN = bondLevelState77N;
  const ScalarT* state78PtrN = bondLevelState78N;
  const ScalarT* state79PtrN = bondLevelState79N;
  const ScalarT* state80PtrN = bondLevelState80N;
  const ScalarT* state81PtrN = bondLevelState81N;
  const ScalarT* state82PtrN = bondLevelState82N;
  const ScalarT* state83PtrN = bondLevelState83N;
  const ScalarT* state84PtrN = bondLevelState84N;
  const ScalarT* state85PtrN = bondLevelState85N;
  const ScalarT* state86PtrN = bondLevelState86N;
  const ScalarT* state87PtrN = bondLevelState87N;
  const ScalarT* state88PtrN = bondLevelState88N;
  const ScalarT* state89PtrN = bondLevelState89N;
  const ScalarT* state90PtrN = bondLevelState90N;
  const ScalarT* state91PtrN = bondLevelState91N;
  const ScalarT* state92PtrN = bondLevelState92N;
  const ScalarT* state93PtrN = bondLevelState93N;
  const ScalarT* state94PtrN = bondLevelState94N;
  const ScalarT* state95PtrN = bondLevelState95N;
  const ScalarT* state96PtrN = bondLevelState96N;
  const ScalarT* state97PtrN = bondLevelState97N;
  const ScalarT* state98PtrN = bondLevelState98N;
  const ScalarT* state99PtrN = bondLevelState99N;
  const ScalarT* state100PtrN = bondLevelState100N;
  const ScalarT* state101PtrN = bondLevelState101N;
  const ScalarT* state102PtrN = bondLevelState102N;
  const ScalarT* state103PtrN = bondLevelState103N;
  const ScalarT* state104PtrN = bondLevelState104N;
  const ScalarT* state105PtrN = bondLevelState105N;
  const ScalarT* state106PtrN = bondLevelState106N;
  const ScalarT* state107PtrN = bondLevelState107N;
  const ScalarT* state108PtrN = bondLevelState108N;
  const ScalarT* state109PtrN = bondLevelState109N;
  const ScalarT* state110PtrN = bondLevelState110N;
  const ScalarT* state111PtrN = bondLevelState111N;
  const ScalarT* state112PtrN = bondLevelState112N;
  const ScalarT* state113PtrN = bondLevelState113N;
  const ScalarT* state114PtrN = bondLevelState114N;
  const ScalarT* state115PtrN = bondLevelState115N;
  const ScalarT* state116PtrN = bondLevelState116N;
  const ScalarT* state117PtrN = bondLevelState117N;
  const ScalarT* state118PtrN = bondLevelState118N;
  const ScalarT* state119PtrN = bondLevelState119N;
  const ScalarT* state120PtrN = bondLevelState120N;
  const ScalarT* state121PtrN = bondLevelState121N;
  const ScalarT* state122PtrN = bondLevelState122N;
  const ScalarT* state123PtrN = bondLevelState123N;
  const ScalarT* state124PtrN = bondLevelState124N;
  const ScalarT* state125PtrN = bondLevelState125N;
  const ScalarT* state126PtrN = bondLevelState126N;
  const ScalarT* state127PtrN = bondLevelState127N;
  const ScalarT* state128PtrN = bondLevelState128N;
  const ScalarT* state129PtrN = bondLevelState129N;
  const ScalarT* state130PtrN = bondLevelState130N;
  const ScalarT* state131PtrN = bondLevelState131N;
  const ScalarT* state132PtrN = bondLevelState132N;
  const ScalarT* state133PtrN = bondLevelState133N;
  const ScalarT* state134PtrN = bondLevelState134N;
  const ScalarT* state135PtrN = bondLevelState135N;
  const ScalarT* state136PtrN = bondLevelState136N;
  const ScalarT* state137PtrN = bondLevelState137N;
  const ScalarT* state138PtrN = bondLevelState138N;
  const ScalarT* state139PtrN = bondLevelState139N;
  const ScalarT* state140PtrN = bondLevelState140N;
  const ScalarT* state141PtrN = bondLevelState141N;
  const ScalarT* state142PtrN = bondLevelState142N;
  const ScalarT* state143PtrN = bondLevelState143N;
  const ScalarT* state144PtrN = bondLevelState144N;
  const ScalarT* state145PtrN = bondLevelState145N;
  const ScalarT* state146PtrN = bondLevelState146N;
  const ScalarT* state147PtrN = bondLevelState147N;
  const ScalarT* state148PtrN = bondLevelState148N;
  const ScalarT* state149PtrN = bondLevelState149N;
  const ScalarT* state150PtrN = bondLevelState150N;
  const ScalarT* state151PtrN = bondLevelState151N;
  const ScalarT* state152PtrN = bondLevelState152N;
  const ScalarT* state153PtrN = bondLevelState153N;
  const ScalarT* state154PtrN = bondLevelState154N;
  const ScalarT* state155PtrN = bondLevelState155N;
  const ScalarT* state156PtrN = bondLevelState156N;
  const ScalarT* state157PtrN = bondLevelState157N;
  const ScalarT* state158PtrN = bondLevelState158N;
  const ScalarT* state159PtrN = bondLevelState159N;
  const ScalarT* state160PtrN = bondLevelState160N;
  const ScalarT* state161PtrN = bondLevelState161N;
  const ScalarT* state162PtrN = bondLevelState162N;
  const ScalarT* state163PtrN = bondLevelState163N;
  const ScalarT* state164PtrN = bondLevelState164N;
  const ScalarT* state165PtrN = bondLevelState165N;
  const ScalarT* state166PtrN = bondLevelState166N;
  const ScalarT* state167PtrN = bondLevelState167N;
  const ScalarT* state168PtrN = bondLevelState168N;
  const ScalarT* state169PtrN = bondLevelState169N;
  const ScalarT* state170PtrN = bondLevelState170N;
  const ScalarT* state171PtrN = bondLevelState171N;
  const ScalarT* state172PtrN = bondLevelState172N;
  const ScalarT* state173PtrN = bondLevelState173N;
  const ScalarT* state174PtrN = bondLevelState174N;
  const ScalarT* state175PtrN = bondLevelState175N;
  const ScalarT* state176PtrN = bondLevelState176N;
  const ScalarT* state177PtrN = bondLevelState177N;
  const ScalarT* state178PtrN = bondLevelState178N;
  const ScalarT* state179PtrN = bondLevelState179N;
  const ScalarT* state180PtrN = bondLevelState180N;
  const ScalarT* state181PtrN = bondLevelState181N;
  const ScalarT* state182PtrN = bondLevelState182N;
  const ScalarT* state183PtrN = bondLevelState183N;
  const ScalarT* state184PtrN = bondLevelState184N;
  const ScalarT* state185PtrN = bondLevelState185N;
  const ScalarT* state186PtrN = bondLevelState186N;
  const ScalarT* state187PtrN = bondLevelState187N;
  const ScalarT* state188PtrN = bondLevelState188N;
  const ScalarT* state189PtrN = bondLevelState189N;

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
  ScalarT* state50PtrNP1 = bondLevelState50NP1;
  ScalarT* state51PtrNP1 = bondLevelState51NP1;
  ScalarT* state52PtrNP1 = bondLevelState52NP1;
  ScalarT* state53PtrNP1 = bondLevelState53NP1;
  ScalarT* state54PtrNP1 = bondLevelState54NP1;
  ScalarT* state55PtrNP1 = bondLevelState55NP1;
  ScalarT* state56PtrNP1 = bondLevelState56NP1;
  ScalarT* state57PtrNP1 = bondLevelState57NP1;
  ScalarT* state58PtrNP1 = bondLevelState58NP1;
  ScalarT* state59PtrNP1 = bondLevelState59NP1;
  ScalarT* state60PtrNP1 = bondLevelState60NP1;
  ScalarT* state61PtrNP1 = bondLevelState61NP1;
  ScalarT* state62PtrNP1 = bondLevelState62NP1;
  ScalarT* state63PtrNP1 = bondLevelState63NP1;
  ScalarT* state64PtrNP1 = bondLevelState64NP1;
  ScalarT* state65PtrNP1 = bondLevelState65NP1;
  ScalarT* state66PtrNP1 = bondLevelState66NP1;
  ScalarT* state67PtrNP1 = bondLevelState67NP1;
  ScalarT* state68PtrNP1 = bondLevelState68NP1;
  ScalarT* state69PtrNP1 = bondLevelState69NP1;
  ScalarT* state70PtrNP1 = bondLevelState70NP1;
  ScalarT* state71PtrNP1 = bondLevelState71NP1;
  ScalarT* state72PtrNP1 = bondLevelState72NP1;
  ScalarT* state73PtrNP1 = bondLevelState73NP1;
  ScalarT* state74PtrNP1 = bondLevelState74NP1;
  ScalarT* state75PtrNP1 = bondLevelState75NP1;
  ScalarT* state76PtrNP1 = bondLevelState76NP1;
  ScalarT* state77PtrNP1 = bondLevelState77NP1;
  ScalarT* state78PtrNP1 = bondLevelState78NP1;
  ScalarT* state79PtrNP1 = bondLevelState79NP1;
  ScalarT* state80PtrNP1 = bondLevelState80NP1;
  ScalarT* state81PtrNP1 = bondLevelState81NP1;
  ScalarT* state82PtrNP1 = bondLevelState82NP1;
  ScalarT* state83PtrNP1 = bondLevelState83NP1;
  ScalarT* state84PtrNP1 = bondLevelState84NP1;
  ScalarT* state85PtrNP1 = bondLevelState85NP1;
  ScalarT* state86PtrNP1 = bondLevelState86NP1;
  ScalarT* state87PtrNP1 = bondLevelState87NP1;
  ScalarT* state88PtrNP1 = bondLevelState88NP1;
  ScalarT* state89PtrNP1 = bondLevelState89NP1;
  ScalarT* state90PtrNP1 = bondLevelState90NP1;
  ScalarT* state91PtrNP1 = bondLevelState91NP1;
  ScalarT* state92PtrNP1 = bondLevelState92NP1;
  ScalarT* state93PtrNP1 = bondLevelState93NP1;
  ScalarT* state94PtrNP1 = bondLevelState94NP1;
  ScalarT* state95PtrNP1 = bondLevelState95NP1;
  ScalarT* state96PtrNP1 = bondLevelState96NP1;
  ScalarT* state97PtrNP1 = bondLevelState97NP1;
  ScalarT* state98PtrNP1 = bondLevelState98NP1;
  ScalarT* state99PtrNP1 = bondLevelState99NP1;
  ScalarT* state100PtrNP1 = bondLevelState100NP1;
  ScalarT* state101PtrNP1 = bondLevelState101NP1;
  ScalarT* state102PtrNP1 = bondLevelState102NP1;
  ScalarT* state103PtrNP1 = bondLevelState103NP1;
  ScalarT* state104PtrNP1 = bondLevelState104NP1;
  ScalarT* state105PtrNP1 = bondLevelState105NP1;
  ScalarT* state106PtrNP1 = bondLevelState106NP1;
  ScalarT* state107PtrNP1 = bondLevelState107NP1;
  ScalarT* state108PtrNP1 = bondLevelState108NP1;
  ScalarT* state109PtrNP1 = bondLevelState109NP1;
  ScalarT* state110PtrNP1 = bondLevelState110NP1;
  ScalarT* state111PtrNP1 = bondLevelState111NP1;
  ScalarT* state112PtrNP1 = bondLevelState112NP1;
  ScalarT* state113PtrNP1 = bondLevelState113NP1;
  ScalarT* state114PtrNP1 = bondLevelState114NP1;
  ScalarT* state115PtrNP1 = bondLevelState115NP1;
  ScalarT* state116PtrNP1 = bondLevelState116NP1;
  ScalarT* state117PtrNP1 = bondLevelState117NP1;
  ScalarT* state118PtrNP1 = bondLevelState118NP1;
  ScalarT* state119PtrNP1 = bondLevelState119NP1;
  ScalarT* state120PtrNP1 = bondLevelState120NP1;
  ScalarT* state121PtrNP1 = bondLevelState121NP1;
  ScalarT* state122PtrNP1 = bondLevelState122NP1;
  ScalarT* state123PtrNP1 = bondLevelState123NP1;
  ScalarT* state124PtrNP1 = bondLevelState124NP1;
  ScalarT* state125PtrNP1 = bondLevelState125NP1;
  ScalarT* state126PtrNP1 = bondLevelState126NP1;
  ScalarT* state127PtrNP1 = bondLevelState127NP1;
  ScalarT* state128PtrNP1 = bondLevelState128NP1;
  ScalarT* state129PtrNP1 = bondLevelState129NP1;
  ScalarT* state130PtrNP1 = bondLevelState130NP1;
  ScalarT* state131PtrNP1 = bondLevelState131NP1;
  ScalarT* state132PtrNP1 = bondLevelState132NP1;
  ScalarT* state133PtrNP1 = bondLevelState133NP1;
  ScalarT* state134PtrNP1 = bondLevelState134NP1;
  ScalarT* state135PtrNP1 = bondLevelState135NP1;
  ScalarT* state136PtrNP1 = bondLevelState136NP1;
  ScalarT* state137PtrNP1 = bondLevelState137NP1;
  ScalarT* state138PtrNP1 = bondLevelState138NP1;
  ScalarT* state139PtrNP1 = bondLevelState139NP1;
  ScalarT* state140PtrNP1 = bondLevelState140NP1;
  ScalarT* state141PtrNP1 = bondLevelState141NP1;
  ScalarT* state142PtrNP1 = bondLevelState142NP1;
  ScalarT* state143PtrNP1 = bondLevelState143NP1;
  ScalarT* state144PtrNP1 = bondLevelState144NP1;
  ScalarT* state145PtrNP1 = bondLevelState145NP1;
  ScalarT* state146PtrNP1 = bondLevelState146NP1;
  ScalarT* state147PtrNP1 = bondLevelState147NP1;
  ScalarT* state148PtrNP1 = bondLevelState148NP1;
  ScalarT* state149PtrNP1 = bondLevelState149NP1;
  ScalarT* state150PtrNP1 = bondLevelState150NP1;
  ScalarT* state151PtrNP1 = bondLevelState151NP1;
  ScalarT* state152PtrNP1 = bondLevelState152NP1;
  ScalarT* state153PtrNP1 = bondLevelState153NP1;
  ScalarT* state154PtrNP1 = bondLevelState154NP1;
  ScalarT* state155PtrNP1 = bondLevelState155NP1;
  ScalarT* state156PtrNP1 = bondLevelState156NP1;
  ScalarT* state157PtrNP1 = bondLevelState157NP1;
  ScalarT* state158PtrNP1 = bondLevelState158NP1;
  ScalarT* state159PtrNP1 = bondLevelState159NP1;
  ScalarT* state160PtrNP1 = bondLevelState160NP1;
  ScalarT* state161PtrNP1 = bondLevelState161NP1;
  ScalarT* state162PtrNP1 = bondLevelState162NP1;
  ScalarT* state163PtrNP1 = bondLevelState163NP1;
  ScalarT* state164PtrNP1 = bondLevelState164NP1;
  ScalarT* state165PtrNP1 = bondLevelState165NP1;
  ScalarT* state166PtrNP1 = bondLevelState166NP1;
  ScalarT* state167PtrNP1 = bondLevelState167NP1;
  ScalarT* state168PtrNP1 = bondLevelState168NP1;
  ScalarT* state169PtrNP1 = bondLevelState169NP1;
  ScalarT* state170PtrNP1 = bondLevelState170NP1;
  ScalarT* state171PtrNP1 = bondLevelState171NP1;
  ScalarT* state172PtrNP1 = bondLevelState172NP1;
  ScalarT* state173PtrNP1 = bondLevelState173NP1;
  ScalarT* state174PtrNP1 = bondLevelState174NP1;
  ScalarT* state175PtrNP1 = bondLevelState175NP1;
  ScalarT* state176PtrNP1 = bondLevelState176NP1;
  ScalarT* state177PtrNP1 = bondLevelState177NP1;
  ScalarT* state178PtrNP1 = bondLevelState178NP1;
  ScalarT* state179PtrNP1 = bondLevelState179NP1;
  ScalarT* state180PtrNP1 = bondLevelState180NP1;
  ScalarT* state181PtrNP1 = bondLevelState181NP1;
  ScalarT* state182PtrNP1 = bondLevelState182NP1;
  ScalarT* state183PtrNP1 = bondLevelState183NP1;
  ScalarT* state184PtrNP1 = bondLevelState184NP1;
  ScalarT* state185PtrNP1 = bondLevelState185NP1;
  ScalarT* state186PtrNP1 = bondLevelState186NP1;
  ScalarT* state187PtrNP1 = bondLevelState187NP1;
  ScalarT* state188PtrNP1 = bondLevelState188NP1;
  ScalarT* state189PtrNP1 = bondLevelState189NP1;

  ScalarT* internEnergyNP1 = bondLevelInternalEnergyNP1;
  ScalarT* inelasEnergyNP1 = bondLevelInelasticEnergyNP1;

  ScalarT deps[6];
  ScalarT epsN[6];

  ScalarT sigmaN[6];
  ScalarT sigmaNP1[6];

  ScalarT stateN[189];
  ScalarT stateNP1[189];

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
        ++state46PtrN, ++state47PtrN, ++state48PtrN, ++state49PtrN, ++state50PtrN,
        ++state51PtrN, ++state52PtrN, ++state53PtrN, ++state54PtrN, ++state55PtrN,
        ++state56PtrN, ++state57PtrN, ++state58PtrN, ++state59PtrN, ++state60PtrN,
        ++state61PtrN, ++state62PtrN, ++state63PtrN, ++state64PtrN, ++state65PtrN,
        ++state66PtrN, ++state67PtrN, ++state68PtrN, ++state69PtrN, ++state70PtrN,
        ++state71PtrN, ++state72PtrN, ++state73PtrN, ++state74PtrN, ++state75PtrN,
        ++state76PtrN, ++state77PtrN, ++state78PtrN, ++state79PtrN, ++state80PtrN,
        ++state81PtrN, ++state82PtrN, ++state83PtrN, ++state84PtrN, ++state85PtrN,
        ++state86PtrN, ++state87PtrN, ++state88PtrN, ++state89PtrN, ++state90PtrN,
        ++state91PtrN, ++state92PtrN, ++state93PtrN, ++state94PtrN, ++state95PtrN,
        ++state96PtrN, ++state97PtrN, ++state98PtrN, ++state99PtrN, ++state100PtrN,
        ++state101PtrN, ++state102PtrN, ++state103PtrN, ++state104PtrN, ++state105PtrN,
        ++state106PtrN, ++state107PtrN, ++state108PtrN, ++state109PtrN, ++state110PtrN,
        ++state111PtrN, ++state112PtrN, ++state113PtrN, ++state114PtrN, ++state115PtrN,
        ++state116PtrN, ++state117PtrN, ++state118PtrN, ++state119PtrN, ++state120PtrN,
        ++state121PtrN, ++state122PtrN, ++state123PtrN, ++state124PtrN, ++state125PtrN,
        ++state126PtrN, ++state127PtrN, ++state128PtrN, ++state129PtrN, ++state130PtrN,
        ++state131PtrN, ++state132PtrN, ++state133PtrN, ++state134PtrN, ++state135PtrN,
        ++state136PtrN, ++state137PtrN, ++state138PtrN, ++state139PtrN, ++state140PtrN,
        ++state141PtrN, ++state142PtrN, ++state143PtrN, ++state144PtrN, ++state145PtrN,
        ++state146PtrN, ++state147PtrN, ++state148PtrN, ++state149PtrN, ++state150PtrN,
        ++state151PtrN, ++state152PtrN, ++state153PtrN, ++state154PtrN, ++state155PtrN,
        ++state156PtrN, ++state157PtrN, ++state158PtrN, ++state159PtrN, ++state160PtrN,
        ++state161PtrN, ++state162PtrN, ++state163PtrN, ++state164PtrN, ++state165PtrN,
        ++state166PtrN, ++state167PtrN, ++state168PtrN, ++state169PtrN, ++state170PtrN,
        ++state171PtrN, ++state172PtrN, ++state173PtrN, ++state174PtrN, ++state175PtrN,
        ++state176PtrN, ++state177PtrN, ++state178PtrN, ++state179PtrN, ++state180PtrN,
        ++state181PtrN, ++state182PtrN, ++state183PtrN, ++state184PtrN, ++state185PtrN,
        ++state186PtrN, ++state187PtrN, ++state188PtrN, ++state189PtrN,
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
        ++state46PtrNP1, ++state47PtrNP1, ++state48PtrNP1, ++state49PtrNP1, ++state50PtrNP1,
        ++state51PtrNP1, ++state52PtrNP1, ++state53PtrNP1, ++state54PtrNP1, ++state55PtrNP1,
        ++state56PtrNP1, ++state57PtrNP1, ++state58PtrNP1, ++state59PtrNP1, ++state60PtrNP1,
        ++state61PtrNP1, ++state62PtrNP1, ++state63PtrNP1, ++state64PtrNP1, ++state65PtrNP1,
        ++state66PtrNP1, ++state67PtrNP1, ++state68PtrNP1, ++state69PtrNP1, ++state70PtrNP1,
        ++state71PtrNP1, ++state72PtrNP1, ++state73PtrNP1, ++state74PtrNP1, ++state75PtrNP1,
        ++state76PtrNP1, ++state77PtrNP1, ++state78PtrNP1, ++state79PtrNP1, ++state80PtrNP1,
        ++state81PtrNP1, ++state82PtrNP1, ++state83PtrNP1, ++state84PtrNP1, ++state85PtrNP1,
        ++state86PtrNP1, ++state87PtrNP1, ++state88PtrNP1, ++state89PtrNP1, ++state90PtrNP1,
        ++state91PtrNP1, ++state92PtrNP1, ++state93PtrNP1, ++state94PtrNP1, ++state95PtrNP1,
        ++state96PtrNP1, ++state97PtrNP1, ++state98PtrNP1, ++state99PtrNP1, ++state100PtrNP1,
        ++state101PtrNP1, ++state102PtrNP1, ++state103PtrNP1, ++state104PtrNP1, ++state105PtrNP1,
        ++state106PtrNP1, ++state107PtrNP1, ++state108PtrNP1, ++state109PtrNP1, ++state110PtrNP1,
        ++state111PtrNP1, ++state112PtrNP1, ++state113PtrNP1, ++state114PtrNP1, ++state115PtrNP1,
        ++state116PtrNP1, ++state117PtrNP1, ++state118PtrNP1, ++state119PtrNP1, ++state120PtrNP1,
        ++state121PtrNP1, ++state122PtrNP1, ++state123PtrNP1, ++state124PtrNP1, ++state125PtrNP1,
        ++state126PtrNP1, ++state127PtrNP1, ++state128PtrNP1, ++state129PtrNP1, ++state130PtrNP1,
        ++state131PtrNP1, ++state132PtrNP1, ++state133PtrNP1, ++state134PtrNP1, ++state135PtrNP1,
        ++state136PtrNP1, ++state137PtrNP1, ++state138PtrNP1, ++state139PtrNP1, ++state140PtrNP1,
        ++state141PtrNP1, ++state142PtrNP1, ++state143PtrNP1, ++state144PtrNP1, ++state145PtrNP1,
        ++state146PtrNP1, ++state147PtrNP1, ++state148PtrNP1, ++state149PtrNP1, ++state150PtrNP1,
        ++state151PtrNP1, ++state152PtrNP1, ++state153PtrNP1, ++state154PtrNP1, ++state155PtrNP1,
        ++state156PtrNP1, ++state157PtrNP1, ++state158PtrNP1, ++state159PtrNP1, ++state160PtrNP1,
        ++state161PtrNP1, ++state162PtrNP1, ++state163PtrNP1, ++state164PtrNP1, ++state165PtrNP1,
        ++state166PtrNP1, ++state167PtrNP1, ++state168PtrNP1, ++state169PtrNP1, ++state170PtrNP1,
        ++state171PtrNP1, ++state172PtrNP1, ++state173PtrNP1, ++state174PtrNP1, ++state175PtrNP1,
        ++state176PtrNP1, ++state177PtrNP1, ++state178PtrNP1, ++state179PtrNP1, ++state180PtrNP1,
        ++state181PtrNP1, ++state182PtrNP1, ++state183PtrNP1, ++state184PtrNP1, ++state185PtrNP1,
        ++state186PtrNP1, ++state187PtrNP1, ++state188PtrNP1, ++state189PtrNP1,
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
        stateN[49] = *state50PtrN;
        stateN[50] = *state51PtrN;
        stateN[51] = *state52PtrN;
        stateN[52] = *state53PtrN;
        stateN[53] = *state54PtrN;
        stateN[54] = *state55PtrN;
        stateN[55] = *state56PtrN;
        stateN[56] = *state57PtrN;
        stateN[57] = *state58PtrN;
        stateN[58] = *state59PtrN;
        stateN[59] = *state60PtrN;
        stateN[60] = *state61PtrN;
        stateN[61] = *state62PtrN;
        stateN[62] = *state63PtrN;
        stateN[63] = *state64PtrN;
        stateN[64] = *state65PtrN;
        stateN[65] = *state66PtrN;
        stateN[66] = *state67PtrN;
        stateN[67] = *state68PtrN;
        stateN[68] = *state69PtrN;
        stateN[69] = *state70PtrN;
        stateN[70] = *state71PtrN;
        stateN[71] = *state72PtrN;
        stateN[72] = *state73PtrN;
        stateN[73] = *state74PtrN;
        stateN[74] = *state75PtrN;
        stateN[75] = *state76PtrN;
        stateN[76] = *state77PtrN;
        stateN[77] = *state78PtrN;
        stateN[78] = *state79PtrN;
        stateN[79] = *state80PtrN;
        stateN[80] = *state81PtrN;
        stateN[81] = *state82PtrN;
        stateN[82] = *state83PtrN;
        stateN[83] = *state84PtrN;
        stateN[84] = *state85PtrN;
        stateN[85] = *state86PtrN;
        stateN[86] = *state87PtrN;
        stateN[87] = *state88PtrN;
        stateN[88] = *state89PtrN;
        stateN[89] = *state90PtrN;
        stateN[90] = *state91PtrN;
        stateN[91] = *state92PtrN;
        stateN[92] = *state93PtrN;
        stateN[93] = *state94PtrN;
        stateN[94] = *state95PtrN;
        stateN[95] = *state96PtrN;
        stateN[96] = *state97PtrN;
        stateN[97] = *state98PtrN;
        stateN[98] = *state99PtrN;
        stateN[99] = *state100PtrN;
        stateN[100] = *state101PtrN;
        stateN[101] = *state102PtrN;
        stateN[102] = *state103PtrN;
        stateN[103] = *state104PtrN;
        stateN[104] = *state105PtrN;
        stateN[105] = *state106PtrN;
        stateN[106] = *state107PtrN;
        stateN[107] = *state108PtrN;
        stateN[108] = *state109PtrN;
        stateN[109] = *state110PtrN;
        stateN[110] = *state111PtrN;
        stateN[111] = *state112PtrN;
        stateN[112] = *state113PtrN;
        stateN[113] = *state114PtrN;
        stateN[114] = *state115PtrN;
        stateN[115] = *state116PtrN;
        stateN[116] = *state117PtrN;
        stateN[117] = *state118PtrN;
        stateN[118] = *state119PtrN;
        stateN[119] = *state120PtrN;
        stateN[120] = *state121PtrN;
        stateN[121] = *state122PtrN;
        stateN[122] = *state123PtrN;
        stateN[123] = *state124PtrN;
        stateN[124] = *state125PtrN;
        stateN[125] = *state126PtrN;
        stateN[126] = *state127PtrN;
        stateN[127] = *state128PtrN;
        stateN[128] = *state129PtrN;
        stateN[129] = *state130PtrN;
        stateN[130] = *state131PtrN;
        stateN[131] = *state132PtrN;
        stateN[132] = *state133PtrN;
        stateN[133] = *state134PtrN;
        stateN[134] = *state135PtrN;
        stateN[135] = *state136PtrN;
        stateN[136] = *state137PtrN;
        stateN[137] = *state138PtrN;
        stateN[138] = *state139PtrN;
        stateN[139] = *state140PtrN;
        stateN[140] = *state141PtrN;
        stateN[141] = *state142PtrN;
        stateN[142] = *state143PtrN;
        stateN[143] = *state144PtrN;
        stateN[144] = *state145PtrN;
        stateN[145] = *state146PtrN;
        stateN[146] = *state147PtrN;
        stateN[147] = *state148PtrN;
        stateN[148] = *state149PtrN;
        stateN[149] = *state150PtrN;
        stateN[150] = *state151PtrN;
        stateN[151] = *state152PtrN;
        stateN[152] = *state153PtrN;
        stateN[153] = *state154PtrN;
        stateN[154] = *state155PtrN;
        stateN[155] = *state156PtrN;
        stateN[156] = *state157PtrN;
        stateN[157] = *state158PtrN;
        stateN[158] = *state159PtrN;
        stateN[159] = *state160PtrN;
        stateN[160] = *state161PtrN;
        stateN[161] = *state162PtrN;
        stateN[162] = *state163PtrN;
        stateN[163] = *state164PtrN;
        stateN[164] = *state165PtrN;
        stateN[165] = *state166PtrN;
        stateN[166] = *state167PtrN;
        stateN[167] = *state168PtrN;
        stateN[168] = *state169PtrN;
        stateN[169] = *state170PtrN;
        stateN[170] = *state171PtrN;
        stateN[171] = *state172PtrN;
        stateN[172] = *state173PtrN;
        stateN[173] = *state174PtrN;
        stateN[174] = *state175PtrN;
        stateN[175] = *state176PtrN;
        stateN[176] = *state177PtrN;
        stateN[177] = *state178PtrN;
        stateN[178] = *state179PtrN;
        stateN[179] = *state180PtrN;
        stateN[180] = *state181PtrN;
        stateN[181] = *state182PtrN;
        stateN[182] = *state183PtrN;
        stateN[183] = *state184PtrN;
        stateN[184] = *state185PtrN;
        stateN[185] = *state186PtrN;
        stateN[186] = *state187PtrN;
        stateN[187] = *state188PtrN;
        stateN[188] = *state189PtrN;

        // call the m7 function
        m7fmaterial_(&dt, 
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
        *state50PtrNP1 = stateNP1[49];
        *state51PtrNP1 = stateNP1[50];
        *state52PtrNP1 = stateNP1[51];
        *state53PtrNP1 = stateNP1[52];
        *state54PtrNP1 = stateNP1[53];
        *state55PtrNP1 = stateNP1[54];
        *state56PtrNP1 = stateNP1[55];
        *state57PtrNP1 = stateNP1[56];
        *state58PtrNP1 = stateNP1[57];
        *state59PtrNP1 = stateNP1[58];
        *state60PtrNP1 = stateNP1[59];
        *state61PtrNP1 = stateNP1[60];
        *state62PtrNP1 = stateNP1[61];
        *state63PtrNP1 = stateNP1[62];
        *state64PtrNP1 = stateNP1[63];
        *state65PtrNP1 = stateNP1[64];
        *state66PtrNP1 = stateNP1[65];
        *state67PtrNP1 = stateNP1[66];
        *state68PtrNP1 = stateNP1[67];
        *state69PtrNP1 = stateNP1[68];
        *state70PtrNP1 = stateNP1[69];
        *state71PtrNP1 = stateNP1[70];
        *state72PtrNP1 = stateNP1[71];
        *state73PtrNP1 = stateNP1[72];
        *state74PtrNP1 = stateNP1[73];
        *state75PtrNP1 = stateNP1[74];
        *state76PtrNP1 = stateNP1[75];
        *state77PtrNP1 = stateNP1[76];
        *state78PtrNP1 = stateNP1[77];
        *state79PtrNP1 = stateNP1[78];
        *state80PtrNP1 = stateNP1[79];
        *state81PtrNP1 = stateNP1[80];
        *state82PtrNP1 = stateNP1[81];
        *state83PtrNP1 = stateNP1[82];
        *state84PtrNP1 = stateNP1[83];
        *state85PtrNP1 = stateNP1[84];
        *state86PtrNP1 = stateNP1[85];
        *state87PtrNP1 = stateNP1[86];
        *state88PtrNP1 = stateNP1[87];
        *state89PtrNP1 = stateNP1[88];
        *state90PtrNP1 = stateNP1[89];
        *state91PtrNP1 = stateNP1[90];
        *state92PtrNP1 = stateNP1[91];
        *state93PtrNP1 = stateNP1[92];
        *state94PtrNP1 = stateNP1[93];
        *state95PtrNP1 = stateNP1[94];
        *state96PtrNP1 = stateNP1[95];
        *state97PtrNP1 = stateNP1[96];
        *state98PtrNP1 = stateNP1[97];
        *state99PtrNP1 = stateNP1[98];
        *state100PtrNP1 = stateNP1[99];
        *state101PtrNP1 = stateNP1[100];
        *state102PtrNP1 = stateNP1[101];
        *state103PtrNP1 = stateNP1[102];
        *state104PtrNP1 = stateNP1[103];
        *state105PtrNP1 = stateNP1[104];
        *state106PtrNP1 = stateNP1[105];
        *state107PtrNP1 = stateNP1[106];
        *state108PtrNP1 = stateNP1[107];
        *state109PtrNP1 = stateNP1[108];
        *state110PtrNP1 = stateNP1[109];
        *state111PtrNP1 = stateNP1[110];
        *state112PtrNP1 = stateNP1[111];
        *state113PtrNP1 = stateNP1[112];
        *state114PtrNP1 = stateNP1[113];
        *state115PtrNP1 = stateNP1[114];
        *state116PtrNP1 = stateNP1[115];
        *state117PtrNP1 = stateNP1[116];
        *state118PtrNP1 = stateNP1[117];
        *state119PtrNP1 = stateNP1[118];
        *state120PtrNP1 = stateNP1[119];
        *state121PtrNP1 = stateNP1[120];
        *state122PtrNP1 = stateNP1[121];
        *state123PtrNP1 = stateNP1[122];
        *state124PtrNP1 = stateNP1[123];
        *state125PtrNP1 = stateNP1[124];
        *state126PtrNP1 = stateNP1[125];
        *state127PtrNP1 = stateNP1[126];
        *state128PtrNP1 = stateNP1[127];
        *state129PtrNP1 = stateNP1[128];
        *state130PtrNP1 = stateNP1[129];
        *state131PtrNP1 = stateNP1[130];
        *state132PtrNP1 = stateNP1[131];
        *state133PtrNP1 = stateNP1[132];
        *state134PtrNP1 = stateNP1[133];
        *state135PtrNP1 = stateNP1[134];
        *state136PtrNP1 = stateNP1[135];
        *state137PtrNP1 = stateNP1[136];
        *state138PtrNP1 = stateNP1[137];
        *state139PtrNP1 = stateNP1[138];
        *state140PtrNP1 = stateNP1[139];
        *state141PtrNP1 = stateNP1[140];
        *state142PtrNP1 = stateNP1[141];
        *state143PtrNP1 = stateNP1[142];
        *state144PtrNP1 = stateNP1[143];
        *state145PtrNP1 = stateNP1[144];
        *state146PtrNP1 = stateNP1[145];
        *state147PtrNP1 = stateNP1[146];
        *state148PtrNP1 = stateNP1[147];
        *state149PtrNP1 = stateNP1[148];
        *state150PtrNP1 = stateNP1[149];
        *state151PtrNP1 = stateNP1[150];
        *state152PtrNP1 = stateNP1[151];
        *state153PtrNP1 = stateNP1[152];
        *state154PtrNP1 = stateNP1[153];
        *state155PtrNP1 = stateNP1[154];
        *state156PtrNP1 = stateNP1[155];
        *state157PtrNP1 = stateNP1[156];
        *state158PtrNP1 = stateNP1[157];
        *state159PtrNP1 = stateNP1[158];
        *state160PtrNP1 = stateNP1[159];
        *state161PtrNP1 = stateNP1[160];
        *state162PtrNP1 = stateNP1[161];
        *state163PtrNP1 = stateNP1[162];
        *state164PtrNP1 = stateNP1[163];
        *state165PtrNP1 = stateNP1[164];
        *state166PtrNP1 = stateNP1[165];
        *state167PtrNP1 = stateNP1[166];
        *state168PtrNP1 = stateNP1[167];
        *state169PtrNP1 = stateNP1[168];
        *state170PtrNP1 = stateNP1[169];
        *state171PtrNP1 = stateNP1[170];
        *state172PtrNP1 = stateNP1[171];
        *state173PtrNP1 = stateNP1[172];
        *state174PtrNP1 = stateNP1[173];
        *state175PtrNP1 = stateNP1[174];
        *state176PtrNP1 = stateNP1[175];
        *state177PtrNP1 = stateNP1[176];
        *state178PtrNP1 = stateNP1[177];
        *state179PtrNP1 = stateNP1[178];
        *state180PtrNP1 = stateNP1[179];
        *state181PtrNP1 = stateNP1[180];
        *state182PtrNP1 = stateNP1[181];
        *state183PtrNP1 = stateNP1[182];
        *state184PtrNP1 = stateNP1[183];
        *state185PtrNP1 = stateNP1[184];
        *state186PtrNP1 = stateNP1[185];
        *state187PtrNP1 = stateNP1[186];
        *state188PtrNP1 = stateNP1[187];
        *state189PtrNP1 = stateNP1[188];
    }
  }
}


// Explicit template instantiation for double
template void updateMicroplaneM7Stress<double>
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
    const double* state50N,
    const double* state51N,
    const double* state52N,
    const double* state53N,
    const double* state54N,
    const double* state55N,
    const double* state56N,
    const double* state57N,
    const double* state58N,
    const double* state59N,
    const double* state60N,
    const double* state61N,
    const double* state62N,
    const double* state63N,
    const double* state64N,
    const double* state65N,
    const double* state66N,
    const double* state67N,
    const double* state68N,
    const double* state69N,
    const double* state70N,
    const double* state71N,
    const double* state72N,
    const double* state73N,
    const double* state74N,
    const double* state75N,
    const double* state76N,
    const double* state77N,
    const double* state78N,
    const double* state79N,
    const double* state80N,
    const double* state81N,
    const double* state82N,
    const double* state83N,
    const double* state84N,
    const double* state85N,
    const double* state86N,
    const double* state87N,
    const double* state88N,
    const double* state89N,
    const double* state90N,
    const double* state91N,
    const double* state92N,
    const double* state93N,
    const double* state94N,
    const double* state95N,
    const double* state96N,
    const double* state97N,
    const double* state98N,
    const double* state99N,
    const double* state100N,
    const double* state101N,
    const double* state102N,
    const double* state103N,
    const double* state104N,
    const double* state105N,
    const double* state106N,
    const double* state107N,
    const double* state108N,
    const double* state109N,
    const double* state110N,
    const double* state111N,
    const double* state112N,
    const double* state113N,
    const double* state114N,
    const double* state115N,
    const double* state116N,
    const double* state117N,
    const double* state118N,
    const double* state119N,
    const double* state120N,
    const double* state121N,
    const double* state122N,
    const double* state123N,
    const double* state124N,
    const double* state125N,
    const double* state126N,
    const double* state127N,
    const double* state128N,
    const double* state129N,
    const double* state130N,
    const double* state131N,
    const double* state132N,
    const double* state133N,
    const double* state134N,
    const double* state135N,
    const double* state136N,
    const double* state137N,
    const double* state138N,
    const double* state139N,
    const double* state140N,
    const double* state141N,
    const double* state142N,
    const double* state143N,
    const double* state144N,
    const double* state145N,
    const double* state146N,
    const double* state147N,
    const double* state148N,
    const double* state149N,
    const double* state150N,
    const double* state151N,
    const double* state152N,
    const double* state153N,
    const double* state154N,
    const double* state155N,
    const double* state156N,
    const double* state157N,
    const double* state158N,
    const double* state159N,
    const double* state160N,
    const double* state161N,
    const double* state162N,
    const double* state163N,
    const double* state164N,
    const double* state165N,
    const double* state166N,
    const double* state167N,
    const double* state168N,
    const double* state169N,
    const double* state170N,
    const double* state171N,
    const double* state172N,
    const double* state173N,
    const double* state174N,
    const double* state175N,
    const double* state176N,
    const double* state177N,
    const double* state178N,
    const double* state179N,
    const double* state180N,
    const double* state181N,
    const double* state182N,
    const double* state183N,
    const double* state184N,
    const double* state185N,
    const double* state186N,
    const double* state187N,
    const double* state188N,
    const double* state189N,
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
    double* state50NP1,
    double* state51NP1,
    double* state52NP1,
    double* state53NP1,
    double* state54NP1,
    double* state55NP1,
    double* state56NP1,
    double* state57NP1,
    double* state58NP1,
    double* state59NP1,
    double* state60NP1,
    double* state61NP1,
    double* state62NP1,
    double* state63NP1,
    double* state64NP1,
    double* state65NP1,
    double* state66NP1,
    double* state67NP1,
    double* state68NP1,
    double* state69NP1,
    double* state70NP1,
    double* state71NP1,
    double* state72NP1,
    double* state73NP1,
    double* state74NP1,
    double* state75NP1,
    double* state76NP1,
    double* state77NP1,
    double* state78NP1,
    double* state79NP1,
    double* state80NP1,
    double* state81NP1,
    double* state82NP1,
    double* state83NP1,
    double* state84NP1,
    double* state85NP1,
    double* state86NP1,
    double* state87NP1,
    double* state88NP1,
    double* state89NP1,
    double* state90NP1,
    double* state91NP1,
    double* state92NP1,
    double* state93NP1,
    double* state94NP1,
    double* state95NP1,
    double* state96NP1,
    double* state97NP1,
    double* state98NP1,
    double* state99NP1,
    double* state100NP1,
    double* state101NP1,
    double* state102NP1,
    double* state103NP1,
    double* state104NP1,
    double* state105NP1,
    double* state106NP1,
    double* state107NP1,
    double* state108NP1,
    double* state109NP1,
    double* state110NP1,
    double* state111NP1,
    double* state112NP1,
    double* state113NP1,
    double* state114NP1,
    double* state115NP1,
    double* state116NP1,
    double* state117NP1,
    double* state118NP1,
    double* state119NP1,
    double* state120NP1,
    double* state121NP1,
    double* state122NP1,
    double* state123NP1,
    double* state124NP1,
    double* state125NP1,
    double* state126NP1,
    double* state127NP1,
    double* state128NP1,
    double* state129NP1,
    double* state130NP1,
    double* state131NP1,
    double* state132NP1,
    double* state133NP1,
    double* state134NP1,
    double* state135NP1,
    double* state136NP1,
    double* state137NP1,
    double* state138NP1,
    double* state139NP1,
    double* state140NP1,
    double* state141NP1,
    double* state142NP1,
    double* state143NP1,
    double* state144NP1,
    double* state145NP1,
    double* state146NP1,
    double* state147NP1,
    double* state148NP1,
    double* state149NP1,
    double* state150NP1,
    double* state151NP1,
    double* state152NP1,
    double* state153NP1,
    double* state154NP1,
    double* state155NP1,
    double* state156NP1,
    double* state157NP1,
    double* state158NP1,
    double* state159NP1,
    double* state160NP1,
    double* state161NP1,
    double* state162NP1,
    double* state163NP1,
    double* state164NP1,
    double* state165NP1,
    double* state166NP1,
    double* state167NP1,
    double* state168NP1,
    double* state169NP1,
    double* state170NP1,
    double* state171NP1,
    double* state172NP1,
    double* state173NP1,
    double* state174NP1,
    double* state175NP1,
    double* state176NP1,
    double* state177NP1,
    double* state178NP1,
    double* state179NP1,
    double* state180NP1,
    double* state181NP1,
    double* state182NP1,
    double* state183NP1,
    double* state184NP1,
    double* state185NP1,
    double* state186NP1,
    double* state187NP1,
    double* state188NP1,
    double* state189NP1,
    double* internalEnergyNP1,
    double* inelasticEnergyNP1,
    const int numPoints, 
    const double dt
);

template void updateBondLevelMicroplaneM7Stress<double>
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
    const double* bondLevelState50N,
    const double* bondLevelState51N,
    const double* bondLevelState52N,
    const double* bondLevelState53N,
    const double* bondLevelState54N,
    const double* bondLevelState55N,
    const double* bondLevelState56N,
    const double* bondLevelState57N,
    const double* bondLevelState58N,
    const double* bondLevelState59N,
    const double* bondLevelState60N,
    const double* bondLevelState61N,
    const double* bondLevelState62N,
    const double* bondLevelState63N,
    const double* bondLevelState64N,
    const double* bondLevelState65N,
    const double* bondLevelState66N,
    const double* bondLevelState67N,
    const double* bondLevelState68N,
    const double* bondLevelState69N,
    const double* bondLevelState70N,
    const double* bondLevelState71N,
    const double* bondLevelState72N,
    const double* bondLevelState73N,
    const double* bondLevelState74N,
    const double* bondLevelState75N,
    const double* bondLevelState76N,
    const double* bondLevelState77N,
    const double* bondLevelState78N,
    const double* bondLevelState79N,
    const double* bondLevelState80N,
    const double* bondLevelState81N,
    const double* bondLevelState82N,
    const double* bondLevelState83N,
    const double* bondLevelState84N,
    const double* bondLevelState85N,
    const double* bondLevelState86N,
    const double* bondLevelState87N,
    const double* bondLevelState88N,
    const double* bondLevelState89N,
    const double* bondLevelState90N,
    const double* bondLevelState91N,
    const double* bondLevelState92N,
    const double* bondLevelState93N,
    const double* bondLevelState94N,
    const double* bondLevelState95N,
    const double* bondLevelState96N,
    const double* bondLevelState97N,
    const double* bondLevelState98N,
    const double* bondLevelState99N,
    const double* bondLevelState100N,
    const double* bondLevelState101N,
    const double* bondLevelState102N,
    const double* bondLevelState103N,
    const double* bondLevelState104N,
    const double* bondLevelState105N,
    const double* bondLevelState106N,
    const double* bondLevelState107N,
    const double* bondLevelState108N,
    const double* bondLevelState109N,
    const double* bondLevelState110N,
    const double* bondLevelState111N,
    const double* bondLevelState112N,
    const double* bondLevelState113N,
    const double* bondLevelState114N,
    const double* bondLevelState115N,
    const double* bondLevelState116N,
    const double* bondLevelState117N,
    const double* bondLevelState118N,
    const double* bondLevelState119N,
    const double* bondLevelState120N,
    const double* bondLevelState121N,
    const double* bondLevelState122N,
    const double* bondLevelState123N,
    const double* bondLevelState124N,
    const double* bondLevelState125N,
    const double* bondLevelState126N,
    const double* bondLevelState127N,
    const double* bondLevelState128N,
    const double* bondLevelState129N,
    const double* bondLevelState130N,
    const double* bondLevelState131N,
    const double* bondLevelState132N,
    const double* bondLevelState133N,
    const double* bondLevelState134N,
    const double* bondLevelState135N,
    const double* bondLevelState136N,
    const double* bondLevelState137N,
    const double* bondLevelState138N,
    const double* bondLevelState139N,
    const double* bondLevelState140N,
    const double* bondLevelState141N,
    const double* bondLevelState142N,
    const double* bondLevelState143N,
    const double* bondLevelState144N,
    const double* bondLevelState145N,
    const double* bondLevelState146N,
    const double* bondLevelState147N,
    const double* bondLevelState148N,
    const double* bondLevelState149N,
    const double* bondLevelState150N,
    const double* bondLevelState151N,
    const double* bondLevelState152N,
    const double* bondLevelState153N,
    const double* bondLevelState154N,
    const double* bondLevelState155N,
    const double* bondLevelState156N,
    const double* bondLevelState157N,
    const double* bondLevelState158N,
    const double* bondLevelState159N,
    const double* bondLevelState160N,
    const double* bondLevelState161N,
    const double* bondLevelState162N,
    const double* bondLevelState163N,
    const double* bondLevelState164N,
    const double* bondLevelState165N,
    const double* bondLevelState166N,
    const double* bondLevelState167N,
    const double* bondLevelState168N,
    const double* bondLevelState169N,
    const double* bondLevelState170N,
    const double* bondLevelState171N,
    const double* bondLevelState172N,
    const double* bondLevelState173N,
    const double* bondLevelState174N,
    const double* bondLevelState175N,
    const double* bondLevelState176N,
    const double* bondLevelState177N,
    const double* bondLevelState178N,
    const double* bondLevelState179N,
    const double* bondLevelState180N,
    const double* bondLevelState181N,
    const double* bondLevelState182N,
    const double* bondLevelState183N,
    const double* bondLevelState184N,
    const double* bondLevelState185N,
    const double* bondLevelState186N,
    const double* bondLevelState187N,
    const double* bondLevelState188N,
    const double* bondLevelState189N,
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
    double* bondLevelState50NP1,
    double* bondLevelState51NP1,
    double* bondLevelState52NP1,
    double* bondLevelState53NP1,
    double* bondLevelState54NP1,
    double* bondLevelState55NP1,
    double* bondLevelState56NP1,
    double* bondLevelState57NP1,
    double* bondLevelState58NP1,
    double* bondLevelState59NP1,
    double* bondLevelState60NP1,
    double* bondLevelState61NP1,
    double* bondLevelState62NP1,
    double* bondLevelState63NP1,
    double* bondLevelState64NP1,
    double* bondLevelState65NP1,
    double* bondLevelState66NP1,
    double* bondLevelState67NP1,
    double* bondLevelState68NP1,
    double* bondLevelState69NP1,
    double* bondLevelState70NP1,
    double* bondLevelState71NP1,
    double* bondLevelState72NP1,
    double* bondLevelState73NP1,
    double* bondLevelState74NP1,
    double* bondLevelState75NP1,
    double* bondLevelState76NP1,
    double* bondLevelState77NP1,
    double* bondLevelState78NP1,
    double* bondLevelState79NP1,
    double* bondLevelState80NP1,
    double* bondLevelState81NP1,
    double* bondLevelState82NP1,
    double* bondLevelState83NP1,
    double* bondLevelState84NP1,
    double* bondLevelState85NP1,
    double* bondLevelState86NP1,
    double* bondLevelState87NP1,
    double* bondLevelState88NP1,
    double* bondLevelState89NP1,
    double* bondLevelState90NP1,
    double* bondLevelState91NP1,
    double* bondLevelState92NP1,
    double* bondLevelState93NP1,
    double* bondLevelState94NP1,
    double* bondLevelState95NP1,
    double* bondLevelState96NP1,
    double* bondLevelState97NP1,
    double* bondLevelState98NP1,
    double* bondLevelState99NP1,
    double* bondLevelState100NP1,
    double* bondLevelState101NP1,
    double* bondLevelState102NP1,
    double* bondLevelState103NP1,
    double* bondLevelState104NP1,
    double* bondLevelState105NP1,
    double* bondLevelState106NP1,
    double* bondLevelState107NP1,
    double* bondLevelState108NP1,
    double* bondLevelState109NP1,
    double* bondLevelState110NP1,
    double* bondLevelState111NP1,
    double* bondLevelState112NP1,
    double* bondLevelState113NP1,
    double* bondLevelState114NP1,
    double* bondLevelState115NP1,
    double* bondLevelState116NP1,
    double* bondLevelState117NP1,
    double* bondLevelState118NP1,
    double* bondLevelState119NP1,
    double* bondLevelState120NP1,
    double* bondLevelState121NP1,
    double* bondLevelState122NP1,
    double* bondLevelState123NP1,
    double* bondLevelState124NP1,
    double* bondLevelState125NP1,
    double* bondLevelState126NP1,
    double* bondLevelState127NP1,
    double* bondLevelState128NP1,
    double* bondLevelState129NP1,
    double* bondLevelState130NP1,
    double* bondLevelState131NP1,
    double* bondLevelState132NP1,
    double* bondLevelState133NP1,
    double* bondLevelState134NP1,
    double* bondLevelState135NP1,
    double* bondLevelState136NP1,
    double* bondLevelState137NP1,
    double* bondLevelState138NP1,
    double* bondLevelState139NP1,
    double* bondLevelState140NP1,
    double* bondLevelState141NP1,
    double* bondLevelState142NP1,
    double* bondLevelState143NP1,
    double* bondLevelState144NP1,
    double* bondLevelState145NP1,
    double* bondLevelState146NP1,
    double* bondLevelState147NP1,
    double* bondLevelState148NP1,
    double* bondLevelState149NP1,
    double* bondLevelState150NP1,
    double* bondLevelState151NP1,
    double* bondLevelState152NP1,
    double* bondLevelState153NP1,
    double* bondLevelState154NP1,
    double* bondLevelState155NP1,
    double* bondLevelState156NP1,
    double* bondLevelState157NP1,
    double* bondLevelState158NP1,
    double* bondLevelState159NP1,
    double* bondLevelState160NP1,
    double* bondLevelState161NP1,
    double* bondLevelState162NP1,
    double* bondLevelState163NP1,
    double* bondLevelState164NP1,
    double* bondLevelState165NP1,
    double* bondLevelState166NP1,
    double* bondLevelState167NP1,
    double* bondLevelState168NP1,
    double* bondLevelState169NP1,
    double* bondLevelState170NP1,
    double* bondLevelState171NP1,
    double* bondLevelState172NP1,
    double* bondLevelState173NP1,
    double* bondLevelState174NP1,
    double* bondLevelState175NP1,
    double* bondLevelState176NP1,
    double* bondLevelState177NP1,
    double* bondLevelState178NP1,
    double* bondLevelState179NP1,
    double* bondLevelState180NP1,
    double* bondLevelState181NP1,
    double* bondLevelState182NP1,
    double* bondLevelState183NP1,
    double* bondLevelState184NP1,
    double* bondLevelState185NP1,
    double* bondLevelState186NP1,
    double* bondLevelState187NP1,
    double* bondLevelState188NP1,
    double* bondLevelState189NP1,
    double* bondLevelInternalEnergyNP1,
    double* bondLevelInelasticEnergyNP1,
    const int* neighborhoodList,
    const int numPoints, 
    const double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */

}
