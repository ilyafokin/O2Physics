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

#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <TMath.h>

#include <array>
#include <vector>

namespace o2::aod {
    namespace identitycollision {
        DECLARE_SOA_COLUMN(CentV0, cent, float);                                   //! Centrality
        DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);                               //! Vertex Z position
        DECLARE_SOA_COLUMN(CollisionW, collisionW, float[4]);
    }; // namespace identitycollision

    namespace identitytrack {
        DECLARE_SOA_COLUMN(Eta, eta, float);                                       //! Eta
        DECLARE_SOA_COLUMN(Phi, phi, float);                                       //! Phi
        DECLARE_SOA_COLUMN(Sign, sign, char);                                     //! Phi
        DECLARE_SOA_COLUMN(P, p, float);                                           //! Momentum at the vertex
        DECLARE_SOA_COLUMN(Ptot, ptot, float);                                     //! Momentum at the TPC inner wall
        DECLARE_SOA_COLUMN(Pt, pt, float);                                         //! Transverse momentum
        DECLARE_SOA_COLUMN(Chi2Ndf, chi2ndf, float);                               //! Chi2 / number of tpc clusters
        DECLARE_SOA_COLUMN(Omegas, omegas, float[4]);                              //! probabilities ω for each particle species
    }; // namespace identitytrack

    namespace identityws {
        DECLARE_SOA_COLUMN(CollisionW, collisionW, float[4]);                      //! Event W values for electrons
    }; // namespace identityws

    namespace identityprobs {
        DECLARE_SOA_COLUMN(Rhos, rhos, float[4]);                                  
        DECLARE_SOA_COLUMN(Omegas, omegas, float[4]);                              //! probabilities ω for each particle species
    } // namespace identityprobs
    
    DECLARE_SOA_TABLE(IdentityCollision, "AOD", "IDENTITYCOLL",
                    o2::soa::Index<>,
                    identitycollision::CentV0,
                    identitycollision::VertexZ);
    
    DECLARE_SOA_TABLE(IdentityTrack, "AOD", "IDENTITYTRACK",
                    o2::soa::Index<>,
                    identitytrack::Eta,
                    identitytrack::Phi,
                    identitytrack::Sign,
                    identitytrack::P,
                    identitytrack::Ptot,
                    identitytrack::Pt,
                    identitytrack::Chi2Ndf
                    );

    DECLARE_SOA_TABLE(IdentityW, "AOD", "IDENTITYW",
                    o2::soa::Index<>,
                    identityws::CollisionW);

    DECLARE_SOA_TABLE(IdentityProb, "AOD", "IDENTITYPROB",
                    o2::soa::Index<>,
                    identityprobs::Rhos,
                    identityprobs::Omegas);
    
}; // namespace o2::aod

