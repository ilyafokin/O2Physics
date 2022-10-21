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

/// \author ilya.fokin@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGCF/DataModel/IdentityTables.h"

#include <array>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum eIdentities {
  kEl,
  kPi,
  kKa,
  kPr
};

namespace o2::aod
{
using myCollisions = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::EvSels>;
using myTracks = soa::Join<aod::Tracks, aod::TracksExtra, o2::aod::TracksDCA, aod::TrackSelection,
                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
using myCollision = myCollisions::iterator;
using myTrack = myTracks::iterator;
} // namespace o2::aod

struct skimmerIdentity {

  Produces<aod::IdentityCollisions> fIdentityCollisions;
  Produces<aod::IdentityTracks> fIdentityTracks;

  // loose filters to reduce tree size
  Configurable<float> fMomMinLoose{"fMomMinLoose", 0.25, "minimum momentum"};
  Configurable<float> fMomMaxLoose{"fMomMaxLoose", 5., "maximum momentum"};
  Configurable<float> fEtaMaxLoose{"fEtaMaxLoose", 0.8, "maximum pseudorapidity"};
  Configurable<float> fChi2NdfMaxLoose{"fChi2NdfMaxLoose", 2.5, "maximum chi2/ndf"};
  Configurable<float> fVertexZMaxLoose{"fVertexZMaxLoose", 10, "maximum vertex z"};

  Filter etaFilter = (nabs(aod::track::eta) < fEtaMaxLoose);
  Filter chi2NdfFilter = (aod::track::tpcChi2NCl < fChi2NdfMaxLoose);

  Filter momentumRangeFilter = (aod::track::p > fMomMinLoose) && (aod::track::p < fMomMaxLoose);
  Filter vertexZFilter = (aod::collision::posZ < fVertexZMaxLoose);

  using myFilteredTracks = soa::Filtered<aod::myTracks>;
  using myFilteredTrack = myFilteredTracks::iterator;
  using myFilteredCollisions = soa::Filtered<aod::myCollisions>;
  using myFilteredCollision = myFilteredCollisions::iterator;

  void process(myFilteredCollision const& collision, myFilteredTracks const& tracks) {

    if (collision.sel7() || true) {
      for (auto const& track : tracks) {
        if (!track.hasTPC()) continue;
        float tpcSignalExpected[4] = {//0, 0, 0, 0
          track.tpcSignal() + track.tpcExpSignalDiffEl(),
          track.tpcSignal() + track.tpcExpSignalDiffPi(),
          track.tpcSignal() + track.tpcExpSignalDiffKa(),
          track.tpcSignal() + track.tpcExpSignalDiffPr()
        };
        float tpcSigmaExpected[4] = {0, 0, 0, 0
          // track.tpcExpSigmaEl(),
          // track.tpcExpSigmaPi(),
          // track.tpcExpSigmaKa(),
          // track.tpcExpSigmaPr()
        };
        fIdentityTracks(track.eta(),
                        track.phi(),
                        track.sign(),
                        track.p(),
                        track.tpcInnerParam(),
                        track.pt(),
                        track.tpcChi2NCl(),
                        track.tpcNClsCrossedRows(),
                        track.dcaXY(),
                        track.dcaZ(),
                        track.passedITSRefit(),
                        track.itsNClsInnerBarrel() > 0, // todo: this is spd+sdd instead of spd only
                        track.tpcSignal(),
                        tpcSignalExpected,
                        tpcSigmaExpected
        );
      }

      fIdentityCollisions(collision.centRun2V0M(),
                          collision.posZ());
      
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerIdentity>(cfgc)};
}