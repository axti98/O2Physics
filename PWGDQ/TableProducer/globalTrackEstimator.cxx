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
/// \author Alexander Tiekoetter <alexander.tiekoetter@cern.ch>, Muenster
/// \brief Global Track estimator with track-to-collision association
/// \file globalTrackEstimator.cxx

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/VarManager.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <string>

namespace o2::aod 
{
namespace dqanalysisflags
{
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 8);                      //! Event decision
}

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTSA", dqanalysisflags::IsEventSelected);
}

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEventsBasic = o2::soa::Join<o2::aod::ReducedEvents, o2::aod::ReducedEventsExtended>;
using MyEventsBasicSelected = o2::soa::Join<o2::aod::ReducedEvents, o2::aod::ReducedEventsExtended, o2::aod::EventCuts>;
using MyBarrelTracks = o2::soa::Join<o2::aod::ReducedTracks, o2::aod::ReducedTracksBarrel, o2::aod::ReducedTracksBarrelPID>;

using MyDielectronCandidates = o2::soa::Join<o2::aod::Dielectrons, o2::aod::DielectronsExtra, o2::aod::DielectronsAll>;

struct GlobalTrackEstimator {
    Produces<o2::aod::ReducedEventsGlobalTrk> globalTracks;
    Produces<o2::aod::ReducedEventsMultAzimuth> azimuthalTrackMult;

    Configurable<float> cfgMaxDCAxy{"cfgMaxDCAxy", 2.0, "Maximum DCAxy"};
    Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 2.0, "Maximum DCAz"};
    Configurable<float> cfgMinPt{"cfgMinPt", 0.1, "Minimum required pT"};
    Configurable<float> cfgEtaRange{"cfgEtaRange", 0.8, "Pseudorapidity range"};
    Configurable<unsigned int> cfgMinTPCCluster{"cfgMinTPCCluster", 70, "Minimum number of crossed TPC rows"};
    Configurable<float> cfggRatioTPCClusterOverFindableCluster{"cfgRatioTPCClusterOverFindableCluster", 0.8, "Minimum ratio of crossed TPC rows over findable clusters"};
    Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "Maximum chi2 per cluster in TPC"};
    Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "Maximum chi2 per cluster in ITS"};
    Configurable<bool> cfgRequireITSibAny{"cfgRequireITSibAny", true, "Pseudorapidity range"};
    
    o2::framework::Preslice<o2::aod::ReducedTracksAssoc> preslice = o2::aod::reducedtrack_association::reducedeventId;

    void init(InitContext& context)
    {
        bool isAnyProcessEnabled = context.mOptions.get<bool>("processSkimmed") || context.mOptions.get<bool>("processSkimmedBasic") || context.mOptions.get<bool>("processSkimmedWithZdc") || context.mOptions.get<bool>("processSkimmedWithMultExtra") || context.mOptions.get<bool>("processSkimmedWithMultExtraNoQvector") || context.mOptions.get<bool>("processSkimmedWithMultExtraZdc") || context.mOptions.get<bool>("processSkimmedWithMultExtraZdcFit") || context.mOptions.get<bool>("processSkimmedWithQvectorCentr") || context.mOptions.get<bool>("processSkimmedWithQvectorCentrMerge");
        bool isDummyEnabled = context.mOptions.get<bool>("processDummy");

        if (isDummyEnabled) {
            if (isAnyProcessEnabled) {
                LOG(fatal) << "You have enabled both the dummy process and at least one normal process function! Check your config!";
            }
            return;
        }
        
        if (!isAnyProcessEnabled) {
            return;
        }
    }

    void processGlobalTracksSkimmed(MyEventsBasicSelected const& events, o2::aod::ReducedTracksAssoc const& assocs, MyBarrelTracks const& tracks)
    {
        for(const auto& event : events) {
            if (event.isEventSelected_bit(0)) {
                auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
            }
        }
    } 

    PROCESS_SWITCH(GlobalTrackEstimator, processGlobalTracksSkimmed, "Count global tracks from skimmed data model", false);


    void processAzimuthalRegions(MyDielectronCandidates const& dielectronCandidates)
    {

    }

    PROCESS_SWITCH(GlobalTrackEstimator, processAzimuthalRegions, "Count global tracks in azimuthal regions with respect to a dielectron candidate", false);

    void processDummy(MyEventsBasic const&) 
    {
        
    }

    PROCESS_SWITCH(GlobalTrackEstimator, processDummy, "Dummy process function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GlobalTrackEstimator>(cfgc)};
}