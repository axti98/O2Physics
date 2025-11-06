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

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "bestCollisionTable.h"
#include "Selections.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

using MultData = soa::Join<aod::Collisions, aod::PVMults, aod::MultSelections, aod::MultsExtra>;

using MultMC = aod::MultMCExtras;
using MultMCData = soa::SmallGroups<soa::Join<MultData, aod::Mult2MCExtras>>;

using AllTracks = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra, aod::TrackSelection>;

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

struct MultiplicityDCACheck {

	HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
	Preslice<aod::Tracks> perEvent = aod::track::collisionId;

	void init(InitContext const&)
	{
		AxisSpec axisVtxZ = {40, -20.0, 20.0, "z_{vtx}"};
		AxisSpec axisVertexDiff = {1000, -5.0, 5.0, "z_{vtx} - z_{vtx, MC}"};

		AxisSpec axisTPCRows = {161, 1.0, 161.0, "TPC Rows"};
		AxisSpec TPCrowsOverClusters = {100, 0.0, 4.0, "TPC Rows / Cluster"};
		AxisSpec TPCchi2PerCluster = {100, 0.0, 8.0, "TPC #chi^{2}"};
		AxisSpec ITSchi2PerCluster = {100, 0.0, 50.0, "ITS #chi^{2}"};
		
		AxisSpec axisDCAz = {200, -10.0, 10.0, "DCA_{z}"};
		AxisSpec axisDCAxy = {200, -5.0, 5.0, "DCA_{xy}"};

		histos.add("zVtxDiff", "zVtxDiff", kTH1D, {axisVertexDiff}, true);
		
		histos.add("TPCRows", "TPCRows", kTH1D, {axisTPCRows}, true);
		histos.add("TPCRowsFromPV", "TPCRowsFromPV", kTH1D, {axisTPCRows}, true);
		histos.add("TPCRowsNotPV", "TPCRowsNotPV", kTH1D, {axisTPCRows}, true);

		histos.add("TPCrowsOverClusters", "TPCrowsOverClusters", kTH1D, {TPCrowsOverClusters}, true);
		histos.add("TPCrowsOverClustersFromPV", "TPCrowsOverClustersFromPV", kTH1D, {TPCrowsOverClusters}, true);
		histos.add("TPCrowsOverClustersNotPV", "TPCrowsOverClustersNotPV", kTH1D, {TPCrowsOverClusters}, true);

		histos.add("TPCchi2PerCluster", "TPCchi2PerCluster", kTH1D, {TPCchi2PerCluster}, true);
		histos.add("TPCchi2PerClusterFromPV", "TPCchi2PerClusterFromPV", kTH1D, {TPCchi2PerCluster}, true);
		histos.add("TPCchi2PerClusterNotPV", "TPCchi2PerClusterNotPV", kTH1D, {TPCchi2PerCluster}, true);

		histos.add("ITSchi2PerCluster", "ITSchi2PerCluster", kTH1D, {ITSchi2PerCluster}, true);
		histos.add("ITSchi2PerClusterFromPV", "ITSchi2PerClusterFromPV", kTH1D, {ITSchi2PerCluster}, true);
		histos.add("ITSchi2PerClusterNotPV", "ITSchi2PerClusterNotPV", kTH1D, {ITSchi2PerCluster}, true);

		histos.add("DCAz", "DCAz", kTH1D, {axisDCAz}, true);
		histos.add("DCAzFromPV", "DCAzFromPV", kTH1D, {axisDCAz}, true);
		histos.add("DCAzNotPV", "DCAzNotPV", kTH1D, {axisDCAz}, true);

		histos.add("DCAxy", "DCAxy", kTH1D, {axisDCAxy}, true);
		histos.add("DCAxyFromPV", "DCAxyFromPV", kTH1D, {axisDCAxy}, true);
		histos.add("DCAxyNotPV", "DCAxyNotPV", kTH1D, {axisDCAxy}, true);
	}

	void processDummy(aod::Collision const& c)
	{
		
	}
	PROCESS_SWITCH(MultiplicityDCACheck, processDummy, "process Dummy", true);

	void processPVTrackCheck(aod::Collision const& c, AllTracks const& tracks)
	{
		for(auto& track : tracks)
		{
			histos.fill(HIST("TPCRows"), track.tpcNClsCrossedRows());
			histos.fill(HIST("TPCrowsOverClusters"), track.tpcCrossedRowsOverFindableCls());
			histos.fill(HIST("TPCchi2PerCluster"), track.tpcChi2NCl());
			histos.fill(HIST("ITSchi2PerCluster"), track.itsChi2NCl());
			histos.fill(HIST("DCAz"), track.dcaZ());
			histos.fill(HIST("DCAxy"), track.dcaXY());
			
			if(track.isPVContributor())
			{
				histos.fill(HIST("TPCRowsFromPV"), track.tpcNClsCrossedRows());
				histos.fill(HIST("TPCrowsOverClustersFromPV"), track.tpcCrossedRowsOverFindableCls());
				histos.fill(HIST("TPCchi2PerClusterFromPV"), track.tpcChi2NCl());
				histos.fill(HIST("ITSchi2PerClusterFromPV"), track.itsChi2NCl());
				histos.fill(HIST("DCAzFromPV"), track.dcaZ());
				histos.fill(HIST("DCAxyFromPV"), track.dcaXY());
			} else if(!track.isPVContributor()) {
				histos.fill(HIST("TPCRowsNotPV"), track.tpcNClsCrossedRows());
				histos.fill(HIST("TPCrowsOverClustersNotPV"), track.tpcCrossedRowsOverFindableCls());
				histos.fill(HIST("TPCchi2PerClusterNotPV"), track.tpcChi2NCl());
				histos.fill(HIST("ITSchi2PerClusterNotPV"), track.itsChi2NCl());
				histos.fill(HIST("DCAzNotPV"), track.dcaZ());
				histos.fill(HIST("DCAxyNotPV"), track.dcaXY());
			}
		}
	}
	PROCESS_SWITCH(MultiplicityDCACheck, processPVTrackCheck, "process Track Cuts", true);

	void processResponse(MultMC::iterator const& mc, MultMCData const& data, AllTracks const& tracks)
	{
		double zVtxMC = mc.multMCPVz();

		if(fabs(zVtxMC) < 10.0)
		{
			for(auto& m : data)
			{
				double zVtxReco   = m.multPVz();
				double zVtxDiff = zVtxMC - zVtxReco;
				
				histos.fill(HIST("zVtxDiff"), zVtxDiff);

				if(m.multSel8()
				&& fabs(zVtxReco) < 10.0
				&& m.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)
				&& m.selection_bit(aod::evsel::kNoSameBunchPileup)
				&& m.selection_bit(aod::evsel::kIsVertexITSTPC)
				&& m.selection_bit(aod::evsel::kIsVertexTOFmatched)
				&& m.selection_bit(aod::evsel::kNoTimeFrameBorder))
				{
					auto groupedTracks = tracks.sliceBy(perEvent, m.globalIndex());

					for(auto const& track : groupedTracks)
					{
						if (std::fabs(track.eta()) < 0.8 
						&& track.tpcNClsFound() >= 80 
						&& track.tpcNClsCrossedRows() >= 100
						&& track.isGlobalTrack())
						{
							double dcaXY = track.dcaXY();
							double dcaZ = track.dcaZ();
							
							histos.fill(HIST("globalTrackDCAxyNoCut"), dcaXY);
							if(fabs(zVtxDiff) < 1.0) histos.fill(HIST("globalTrackDCAxy1cm"), dcaXY);
							if(fabs(zVtxDiff) < 0.5)histos.fill(HIST("globalTrackDCAxy05cm"), dcaXY);
							if(fabs(zVtxDiff) < 0.1)histos.fill(HIST("globalTrackDCAxy01cm"), dcaXY);

							histos.fill(HIST("globalTrackDCAzNoCut"), dcaZ);
							if(fabs(zVtxDiff) < 1.0)histos.fill(HIST("globalTrackDCAz1cm"), dcaZ);
							if(fabs(zVtxDiff) < 0.5)histos.fill(HIST("globalTrackDCAz05cm"), dcaZ);
							if(fabs(zVtxDiff) < .1)histos.fill(HIST("globalTrackDCAz01cm"), dcaZ);
						}
					}
				}
			}
		}
	}
	PROCESS_SWITCH(MultiplicityDCACheck, processResponse, "process DCA information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{adaptAnalysisTask<MultiplicityDCACheck>(cfgc)};
}
