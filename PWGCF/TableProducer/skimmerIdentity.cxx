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
using myCollisions = soa::Join<aod::Collisions, aod::CentRun2V0Ms>;
using myTracks = soa::Join<aod::Tracks, aod::TracksExtra, o2::aod::TracksDCA,
                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
using myCollision = myCollisions::iterator;
using myTrack = myTracks::iterator;
} // namespace o2::aod

struct skimmerIdentity {

  Produces<aod::IdentityCollision> fIdentityCollisions;
  Produces<aod::IdentityTrack> fIdentityTracks;

  void process(aod::myCollision const& collision, aod::myTracks const& tracks) {

    for (auto const& track : tracks) {
      fIdentityTracks(track.eta(),
                      track.phi(),
                      track.sign(),
                      track.p(),
                      track.tpcInnerParam(),
                      track.pt(),
                      track.tpcChi2NCl());

    }

    fIdentityCollisions(collision.centRun2V0M(),
                        collision.posZ());
    
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerIdentity>(cfgc)};
}