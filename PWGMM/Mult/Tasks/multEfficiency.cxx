// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TMath.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/MathConstants.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"

#include "bestCollisionTable.h"
#include "Selections.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;
using namespace pwgmm::mult;

using PVMults = soa::Join<aod::PVMults, aod::MultSelections, aod::MultsExtra>;

using MultMC = aod::MultMCExtras;

using MultMCPV = soa::SmallGroups<soa::Join<PVMults, aod::Mult2MCExtras>>;

struct MultiplicityEfficiency {

	HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

	void init(InitContext const&)
	{
		AxisSpec axisNtrkPV = {140, 0.5, 140.5, "N_{trk, PV}"};
		AxisSpec axisNtrkGlobal = {140, 0.5, 140.5, "N_{trk, Global}"};

		histos.add("MultPVINELGt0", 		"MultPVINELGt0", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVBeforeTrigger", 	"MultPVBeforeTrigger", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVAfterTrigger",	"MultPVAfterTrigger", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVBeforeVertex",	"MultPVBeforeVertex", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVAfterVertex", 	"MultPVAfterVertex", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVAfterSameBunch", 	"MultPVAfterSameBunch", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVAfterGoodZVtx", 	"MultPVAfterGoodZVtx", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVAfterTOFMatch", 	"MultPVAfterTOFMatch", kTH1D, {axisNtrkPV}, false);
		histos.add("MultPVAfterVertexITSTPC", "MultPVAfterVertexITSTPC", kTH1D, {axisNtrkPV}, false);
	}

	void processPV(MultMC::iterator const& mc, MultMCPV const& pv)
	{
		if(mc.isInelGt0())
		{	
			histos.fill(HIST("MultPVINELGt0"), mc.multMCNParticlesEta08());

			for(auto& event : pv)
			{
				histos.fill(HIST("MultPVBeforeTrigger"), event.multNTracksPV());

				if(event.multSel8())
				{
					histos.fill(HIST("MultPVAfterTrigger"), event.multNTracksPV());
					histos.fill(HIST("MultPVBeforeVertex"), event.multNTracksPV());
					
					if(abs(event.multPVz()) < 10.0)
					{
						histos.fill(HIST("MultPVAfterVertex"), event.multNTracksPV());
					
					}

					if(event.selection_bit(aod::evsel::kNoSameBunchPileup))  histos.fill(HIST("MultPVAfterSameBunch"), event.multNTracksPV());
					if(event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))  histos.fill(HIST("MultPVAfterGoodZVtx"), event.multNTracksPV());
					if(event.selection_bit(aod::evsel::kIsVertexTOFmatched)) histos.fill(HIST("MultPVAfterTOFMatch"), event.multNTracksPV());
					if(event.selection_bit(aod::evsel::kIsVertexITSTPC)) 	 histos.fill(HIST("MultPVAfterVertexITSTPC"), event.multNTracksPV());
				}
			}
		}
	}
	PROCESS_SWITCH(MultiplicityEfficiency, processPV, "process PV contributrs", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{adaptAnalysisTask<MultiplicityEfficiency>(cfgc)};
}
