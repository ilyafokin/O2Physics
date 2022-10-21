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
#include <TMath.h>

#include <array>
#include <vector>

namespace o2::aod {
    namespace identitycollision {
        DECLARE_SOA_COLUMN(CentV0, cent, float);                                   //! Centrality
        DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);                               //! Vertex Z position
    }; // namespace identitycollision

    namespace identitytrack {
        DECLARE_SOA_COLUMN(Eta, eta, float);                                       //! Eta
        DECLARE_SOA_COLUMN(Phi, phi, float);                                       //! Phi
        DECLARE_SOA_COLUMN(Sign, sign, int);                                       //! Sign
        DECLARE_SOA_COLUMN(P, p, float);                                           //! Momentum at the vertex
        DECLARE_SOA_COLUMN(Ptot, ptot, float);                                     //! Momentum at the TPC inner wall
        DECLARE_SOA_COLUMN(Pt, pt, float);                                         //! Transverse momentum
        DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);                         //! Chi2 / number of tpc clusters
        DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNclsCrossedRows, float);         //! Number of TPC crossed rows
        DECLARE_SOA_COLUMN(Dcaxy, dcaxy, float);                                   //! DCA xy
        DECLARE_SOA_COLUMN(Dcaz, dcaz, float);                                     //! DCA z
        DECLARE_SOA_COLUMN(ItsRefit, itsRefit, bool);                              //! ITS refit
        DECLARE_SOA_COLUMN(ItsPixel, itsPixel, bool);                              //! Hits in ITS pixel detectors

        DECLARE_SOA_COLUMN(TpcSignal, tpcSignal, float);                                     //! tpc signal
        DECLARE_SOA_COLUMN(TpcSignalExpected, tpcSignalExpected, float[4]);                  //! expected tpc signal means for each particle species
        DECLARE_SOA_COLUMN(TpcSigmaExpected, tpcSigmaExpected, float[4]);        //! expected tpc signal widths
    }; // namespace identitytrack

    namespace identityws {
        DECLARE_SOA_COLUMN(CollisionW, collisionW, float[4]);                      //! Event W values for electrons
    }; // namespace identityws

    namespace identityprobs {
        DECLARE_SOA_COLUMN(Rhos, rhos, float[4]);                                  
        DECLARE_SOA_COLUMN(Omegas, omegas, float[4]);                              //! probabilities ω for each particle species
    } // namespace identityprobs
    
    DECLARE_SOA_TABLE(IdentityCollisions, "AOD", "IDENTITYCOLLS",
                    o2::soa::Index<>,
                    identitycollision::CentV0,
                    identitycollision::VertexZ
                    );
    
    DECLARE_SOA_TABLE(IdentityTracks, "AOD", "IDENTITYTRACKS",
                    o2::soa::Index<>,
                    identitytrack::Eta,
                    identitytrack::Phi,
                    identitytrack::Sign,
                    identitytrack::P,
                    identitytrack::Ptot,
                    identitytrack::Pt,
                    identitytrack::TpcChi2NCl,
                    identitytrack::TpcNClsCrossedRows,
                    identitytrack::Dcaxy,
                    identitytrack::Dcaz,
                    identitytrack::ItsRefit,
                    identitytrack::ItsPixel,
                    identitytrack::TpcSignal,
                    identitytrack::TpcSignalExpected,
                    identitytrack::TpcSigmaExpected
                    );

    DECLARE_SOA_TABLE(IdentityWs, "AOD", "IDENTITYWS",
                    o2::soa::Index<>,
                    identityws::CollisionW
                    );

    DECLARE_SOA_TABLE(IdentityProbs, "AOD", "IDENTITYPROBS",
                    o2::soa::Index<>,
                    identityprobs::Rhos,
                    identityprobs::Omegas
                    );
    
}; // namespace o2::aod

