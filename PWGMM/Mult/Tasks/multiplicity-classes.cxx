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
#include <string>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/StaticFor.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "bestCollisionTable.h"
#include "Selections.h"

#define STRINGIFICATOR(X) #X

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

using MultData = soa::Join<aod::PVMults, aod::TPCMults, aod::FT0Mults, aod::MultSelections, aod::MultsExtra>;

using MultMC = aod::MultMCExtras;
using MultMCData = soa::SmallGroups<soa::Join<MultData, aod::Mult2MCExtras>>;

#define EVENT_TOTAL 1
#define EVENT_TRIG 2
#define EVENT_NO_PILEUP 3
#define EVENT_VTX 4

#define EVENT_ALL 1
#define EVENT_NOINT7 2
#define EVENT_SEL8 3
#define EVENT_NOSAMEBUNCH 4
#define EVENT_GOODZVTX 5
#define EVENT_TOFMATCHED 6
#define EVENT_VTXITSTPC 7
#define EVENT_NOTIMEFRAMEBORDER 8

struct MultiplicityClasses {

	HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
	
	Configurable<std::string> dataset{"dataset", "22", "year"};
	Configurable<bool> ismc{"isMC", false, "isMC"};
	
	TFile *efficiency = nullptr; 
	TFile* NpvNtpcCorr = nullptr;
	
	TH1D* efficiencyPVContributors = nullptr;

	void init(InitContext const&)
	{
		TString year = dataset.value;
		bool isMC = ismc.value;

		efficiency = TFile::Open("/data/JpsiAndPsi2SAnalysis/EventEfficiencyOutput/AllEfficiencies" + year + ".root");
		efficiencyPVContributors = (TH1D*) efficiency->Get("EfficiencyPVContributors");

		AxisSpec axisNPV = {141, 0.0, 140.0, "N_{PV}"};
		AxisSpec axisNTPC = {1001, 0.0, 1000.0, "N_{TPC}"};
		AxisSpec axisNFT0 = {141, 0.0, 140.0, "N_{FT0}"};
						
		AxisSpec axisNtruth = {141, 0.0, 140.0, "N_{truth}"};

		AxisSpec axisEventLoss = {9, 0.5, 8.5, "N_{evt}"};
		AxisSpec axisEvent = {5, 0.5, 5.5, "N_{evt}"};
		AxisSpec axisEventMC = {4, 0.5, 4.5, "N_{evt}"};
		
		AxisSpec axisVtxZ = {1000, -20.0, 20.0, "z_{vtx}"};
		
		AxisSpec axisVtxZShort = {40, -10.0, 10.0, "z_{vtx}"};
		
		AxisSpec axisVertexDiff = {1000, -5.0, 5.0, "z_{vtx} - z_{vtx, MC}"};

		histos.add("EventLossHist",		"EventLossHist", kTH1D, {axisEventLoss}, true);
		histos.add("EventHist",    		"EventHist",   	 kTH1D, {axisEvent}, true);
		
		histos.add("ZvtxHist", 	   		"ZvtxHist", 	kTH1D, {axisVtxZ}, true);

		histos.add("ZvtxTrigHist",   	"ZvtxTrigHist", kTH1D, {axisVtxZ}, true);
		
		histos.add("MultPV",			"MultPV", 	   	kTH1D, {axisNPV}, true);
		histos.add("MultPVRaw",			"MultPVRaw", 	kTH1D, {axisNPV}, true);

		histos.add("MultPVRawTimeRange",	"MultPVRawTimeRange", 	   	kTH1D, {axisNPV}, true);

		histos.add("MultTPC",			"MultPV", 	   	kTH1D, {axisNTPC}, true);
		histos.add("MultFT0",			"MultFT0", 	   	kTH1D, {axisNFT0}, true);
		
		histos.add("NPVvsNTPC",		"NPVvsNTPC", 			kTH2D, {axisNPV, axisNTPC}, true);
		histos.add("NPVvsFT0", 		"NPVvsFT0", 			kTH2D, {axisNPV, axisNFT0}, true);
				
		histos.add("NPVvsNTPCTriggered", "NPVvsNTPCTriggered",	kTH2D, {axisNPV, axisNTPC}, true);
		histos.add("NPVvsFT0Triggered",	 "NPVvsFT0Triggered", 	kTH2D, {axisNPV, axisNFT0}, true);
		
		histos.add("zVtxEvents", 		"zVtxEvents", 			kTH1D, {axisVtxZShort}, true);
		histos.add("NPVperZvtx", 		"NPVperZvtx", 			kTH1D, {axisVtxZShort}, true);

		if(isMC)
		{
			histos.add("NTruthEta10", 		"NTruthEta10", 		kTH1D, {axisNtruth}, true);
			histos.add("NTruthEta08", 		"NTruthEta08", 		kTH1D, {axisNtruth}, true);
			histos.add("NTruthEta05", 		"NTruthEta05", 		kTH1D, {axisNtruth}, true);

			histos.add("NPVVsNtruth", 	 "NtrkPVVsNtruth", 	 kTH2D, {axisNPV, axisNtruth},	true);
			histos.add("NTPCVsNtruth", 	 "NtrkPVVsNtruth", 	 kTH2D, {axisNTPC, axisNtruth},	true);

			histos.add("NPVVsNtruthTriggered",	"NtrkPVVsNtruthTriggered", 	 kTH2D, {axisNPV, axisNtruth}, 	 true);
			histos.add("NTPCVsNtruthTriggered",	"NtrkPVVsNtruthTriggered", 	 kTH2D, {axisNTPC, axisNtruth},  true);
			
			histos.add("EventHistMC",    	"EventHistMC",   kTH1D, {axisEventMC}, true);
			histos.add("ZvtxMCHist",   		"ZvtxMCHist",  	kTH1D, {axisVtxZ}, true);
			histos.add("ZvtxMCTrigHist", 	"ZvtxMCTrigHist", 			 kTH1D, {axisVtxZ}, true);
			histos.add("zVtxDiff", 			"zVtxDiff", 				 kTH1D, {axisVertexDiff}, true);
			histos.add("zVtxDiffTriggered", "zVtxDiffTriggered", kTH1D, {axisVertexDiff}, true);
			
			auto eventHistMC = histos.get<TH1>(HIST("EventHistMC"));
			auto* xAxisEventHistMC = eventHistMC->GetXaxis();
			xAxisEventHistMC->SetBinLabel(EVENT_TOTAL, "Total");
			xAxisEventHistMC->SetBinLabel(EVENT_TRIG, "INEL>0");
			xAxisEventHistMC->SetBinLabel(EVENT_NO_PILEUP, "|z_{vtx}| < 10 cm");
			xAxisEventHistMC->SetBinLabel(4, "");
		}

		auto eventLossHist = histos.get<TH1>(HIST("EventLossHist"));
		auto* xAxisEventLossHist = eventLossHist->GetXaxis();
		xAxisEventLossHist->SetBinLabel(EVENT_ALL, 					"All Events");
		xAxisEventLossHist->SetBinLabel(EVENT_NOINT7, 				"|zVtx| < 10");
		xAxisEventLossHist->SetBinLabel(EVENT_SEL8, 				"sel8");
		xAxisEventLossHist->SetBinLabel(EVENT_NOSAMEBUNCH, 			"NoSameBunchPileup");
		xAxisEventLossHist->SetBinLabel(EVENT_GOODZVTX, 			"GoodZvtxvsPV");
		xAxisEventLossHist->SetBinLabel(8, 					"");

		auto eventHist = histos.get<TH1>(HIST("EventHist"));
		auto* xAxisEventHist = eventHist->GetXaxis();
		xAxisEventHist->SetBinLabel(EVENT_TOTAL, "Total");
		xAxisEventHist->SetBinLabel(EVENT_TRIG, "Triggered Sel8");
		xAxisEventHist->SetBinLabel(EVENT_NO_PILEUP, "Pileup rejected");
		xAxisEventHist->SetBinLabel(EVENT_VTX, "|z_{vtx}| < 10 cm");
		xAxisEventHist->SetBinLabel(5, "");
	}

	void processEventData(MultData::iterator const& multData)
	{
		histos.fill(HIST("EventLossHist"), 	EVENT_ALL);
		if(fabs(multData.multPVz()) < 10.0) 						histos.fill(HIST("EventLossHist"), 	EVENT_NOINT7);
		if(multData.multSel8()) 									histos.fill(HIST("EventLossHist"),	EVENT_SEL8);
		if(multData.selection_bit(aod::evsel::kNoSameBunchPileup)) 	histos.fill(HIST("EventLossHist"), 	EVENT_NOSAMEBUNCH);
		if(multData.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) 	histos.fill(HIST("EventLossHist"), 	EVENT_GOODZVTX);
	}
	PROCESS_SWITCH(MultiplicityClasses, processEventData, "process event information", true);

	void processData(MultData::iterator const& multData)
	{
		float zVtx 			= multData.multPVz();
		int pvContributors 	= multData.multNTracksPV();
		int nTPCTracks 		= multData.multTPC();
		int nFT0	 		= multData.multFT0M();

		histos.fill(HIST("EventHist"), EVENT_TOTAL);
		histos.fill(HIST("ZvtxHist"), zVtx);
		histos.fill(HIST("NPVvsNTPC"), pvContributors, nTPCTracks);
		histos.fill(HIST("NPVvsFT0"), pvContributors, nFT0);

		if(multData.multSel8())
		{
			histos.fill(HIST("EventHist"), EVENT_TRIG);

			if(multData.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)
			&& multData.selection_bit(aod::evsel::kNoSameBunchPileup))
			{
				histos.fill(HIST("EventHist"), EVENT_NO_PILEUP);
				histos.fill(HIST("ZvtxTrigHist"), zVtx);

				if(fabs(zVtx) < 10.0)
				{
					histos.fill(HIST("EventHist"), EVENT_VTX);
					histos.fill(HIST("NPVvsNTPCTriggered"), pvContributors, nTPCTracks);
					histos.fill(HIST("NPVvsFT0Triggered"), pvContributors, nFT0);

					histos.fill(HIST("zVtxEvents"), zVtx);
					histos.fill(HIST("NPVperZvtx"), zVtx, pvContributors);

					double weightPV = efficiencyPVContributors->GetBinContent(pvContributors);

					int countPV = 1;

					while(weightPV == 0.0)
					{
						weightPV = efficiencyPVContributors->GetBinContent(pvContributors-countPV); 
						++countPV;
					}
					
					histos.fill(HIST("MultPV"), 	pvContributors, 1.0/weightPV);

					histos.fill(HIST("MultPVRaw"), 	pvContributors);
					histos.fill(HIST("MultTPC"), 	nTPCTracks);
					histos.fill(HIST("MultFT0"), 	nFT0);
				}
			}
		}
	}
	PROCESS_SWITCH(MultiplicityClasses, processData, "process data information", true);

	void processMC(MultMC::iterator const& multMC)
	{
		float zVtx = multMC.multMCPVz();

		histos.fill(HIST("EventHistMC"), EVENT_TOTAL);
		histos.fill(HIST("ZvtxMCHist"), zVtx);

		if(multMC.isInelGt0())
		{
			histos.fill(HIST("EventHistMC"), EVENT_TRIG);
			histos.fill(HIST("ZvtxMCTrigHist"), zVtx);

			if(fabs(zVtx) < 10.0)
			{
				histos.fill(HIST("EventHistMC"), EVENT_VTX);

				histos.fill(HIST("NTruthEta10"), 	multMC.multMCNParticlesEta10());
				histos.fill(HIST("NTruthEta08"), 	multMC.multMCNParticlesEta08());
				histos.fill(HIST("NTruthEta05"), 	multMC.multMCNParticlesEta05());
			}
		}
	}
	PROCESS_SWITCH(MultiplicityClasses, processMC, "process MC information", false);

	void processResponse(MultMC::iterator const& mc, MultMCData const& data)
	{
		if(mc.isInelGt0())
		{
			double zVtxMC = mc.multMCPVz();

			if(fabs(zVtxMC) < 10.0)
			{
				double mcParticlesInEta08 = mc.multMCNParticlesEta08();

				for(auto& m : data)
				{	
					double zVtxReco   = m.multPVz();
					double vtxDiff 	  = zVtxReco - zVtxMC;

					double NPv 		  = static_cast<double>(m.multNTracksPV());
					double nTPCTracks = static_cast<double>(m.multTPC());
					
					histos.fill(HIST("NPVVsNtruth"), NPv, mcParticlesInEta08);
					histos.fill(HIST("NTPCVsNtruth"), nTPCTracks, mcParticlesInEta08);
					
					histos.fill(HIST("zVtxDiff"), vtxDiff);
					
					if(m.multSel8()
					&& fabs(zVtxReco) < 10.0
					&& m.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)
					&& m.selection_bit(aod::evsel::kNoSameBunchPileup))
					{
						histos.fill(HIST("zVtxDiffTriggered"), vtxDiff);
						
						histos.fill(HIST("NPVVsNtruthTriggered"), NPv, mcParticlesInEta08);
						histos.fill(HIST("NTPCVsNtruthTriggered"), nTPCTracks, mcParticlesInEta08);
					}
				}
			}			
		}
	}
	PROCESS_SWITCH(MultiplicityClasses, processResponse, "process Response Matrix", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{adaptAnalysisTask<MultiplicityClasses>(cfgc)};
}
