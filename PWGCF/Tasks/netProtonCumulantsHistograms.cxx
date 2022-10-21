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

#include "PWGCF/Tasks/netProtonCumulants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::vector;
using std::map;

struct netProtonCumulantsHistograms {

  // event cuts
  Configurable<float> fVertexZMax{"fVertexZMax", 7, "maximum vertex z position"};
  Configurable<float> fVertexZMin{"fVertexZMin", 0.15, "maximum vertex z position"};

  // track cuts
  Configurable<float> fEtaMax{"fEtaMax", 0.8, "maximum pseudorapidity"};
  Configurable<float> fNTPCClusterMin{"fNTPCClusterMin", 80, "minimum number of TPC clusters"};
  Configurable<float> fChi2NdfMax{"fChi2NdfMax", 2.5, "maximum chi2/ndf"};
  Configurable<float> fDcaxyA{"fDcaxyA", 0.0208, "constant term in max dcaxy"}; // A + B * pT^C
  Configurable<float> fDcaxyB{"fDcaxyB", 0.04, "factor in max dcaxy"};
  Configurable<float> fDcaxyC{"fDcaxyC", -1.01, "exponent in max dcaxy"};
  Configurable<float> fDcazA{"fDcazA", 0.0208, "constant term in max dcaz"}; // A + B * pT^C
  Configurable<float> fDcazB{"fDcazB", 0.04, "factor in max dcaz"};
  Configurable<float> fDcazC{"fDcazC", -1.01, "exponent in max dcaz"};

  // tree size reduction cuts
  Configurable<float> fMomMin{"fMomMin", 0.25, "minimum momentum"};
  Configurable<float> fMomMax{"fMomMax", 5., "maximum momentum"};

  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};
  HistogramSpec hEtaPhi{"hEtaPhi", "hEtaPhi;#phi;#eta;counts", {HistType::kTH2I, {gAxis_phi, gAxis_eta}}};
  HistogramSpec hdEdxPtot{"hDedxPtot", "hDedxPtot;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dEdx}}};
  HistogramSpec hdEdxPtotBefore{"hDedxPtotBefore", "hDedxPtotBefore;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dEdx}}};

  HistogramSpec hdEdxPtotEl{"hDedxPtotEl", "hDedxPtotEl;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dEdx}}};

  HistogramSpec hCentrality{"hCentrality", "hCentrality;counts", {HistType::kTH1I, {gAxis_cent}}};
  HistogramSpec hP{"hP", "hP;counts", {HistType::kTH1I, {gAxis_ptot}}};
  HistogramSpec hPtot{"hPtot", "hPtot;counts", {HistType::kTH1I, {gAxis_ptot}}};
  HistogramSpec hChi2ndf{"hChi2ndf", "hChi2ndf;#chi^{2}/ndf;counts", {HistType::kTH1I, {gAxis_chi2ndf}}};
  HistogramSpec hNTPCCrossedRows{"hNTPCCrossedRows", "hNTPCCrossedRows;Num of crossed rows;counts", {HistType::kTH1I, {gAxis_ncrossedrows}}};
  HistogramSpec hDcaxyPt{"hDcaxyPt", "hDcaxyPt;p_{T};dca_{xy};counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dcaxy}}};
  HistogramSpec hDcazPt{"hDcazPt", "hDcazPt;p_{T};dca_{z};counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dcaz}}};

  HistogramSpec hDcaxyPtBefore{"hDcaxyPtBefore", "hDcaxyPtBefore;p_{T};dca_{xy};counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dcaxy}}};
  HistogramSpec hDcazPtBefore{"hDcazPtBefore", "hDcazPtBefore;p_{T};dca_{z};counts", {HistType::kTH2I, {gAxis_ptot, gAxis_dcaz}}};

  HistogramSpec hdEdxPtotEtaCent{"hDedxPtotEtaCent", "hDedxPtotEtaCent;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTHnI, {gAxis_ptot, gAxis_dEdx, gAxis_eta, gAxis_cent}}};
  // map<size_t, map<size_t, HistPtr>> mapDedxPtotEtaCent{}; // hDedxPtotEtaCent[iCent][iEta]

  template <typename T>
  bool collisionPassesCuts(const T& theCollision)
  {
    // vertex z cut
    if (fabs(theCollision.vertexZ()) < fVertexZMin || fabs(theCollision.vertexZ()) > fVertexZMax) {
      return false;
    }

  return true;
  }

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // return true;
    // eta cut
    if (!(fabs(theTrack.eta()) < fEtaMax)) {
      return false;
    }

    // number of tpc clusters
    if (!(theTrack.tpcNclsCrossedRows() > fNTPCClusterMin)) {
      return false;
    }

    // chi2/ndf
    if (!(theTrack.tpcChi2NCl() < fChi2NdfMax)) {
      return false;
    }

    // dcaxy
    if (!(fabs(theTrack.dcaxy()) < fDcaxyA + fDcaxyB * pow(theTrack.pt(), fDcaxyC))) {
      return false;
    }

    // dcaz
    if (!(fabs(theTrack.dcaz()) < fDcazA + fDcazB * pow(theTrack.pt(), fDcazC))) {
      return false;
    }

  return true;
  }

  void init(InitContext&){
    LOGF(info, "in init()");
    fHistogramRegistry.add(hEtaPhi);
    fHistogramRegistry.add(hdEdxPtot);
    fHistogramRegistry.add(hdEdxPtotBefore);

    fHistogramRegistry.add(hdEdxPtotEl);

    fHistogramRegistry.add(hCentrality);
    fHistogramRegistry.add(hP);
    fHistogramRegistry.add(hPtot);
    fHistogramRegistry.add(hChi2ndf);
    fHistogramRegistry.add(hNTPCCrossedRows);
    fHistogramRegistry.add(hDcaxyPt);
    fHistogramRegistry.add(hDcazPt);

    fHistogramRegistry.add(hDcaxyPtBefore);
    fHistogramRegistry.add(hDcazPtBefore);

    fHistogramRegistry.add(hdEdxPtotEtaCent);
  }

  void process(aod::IdentityCollisions::iterator const& collision, aod::IdentityTracks const& tracks) {

    if (collisionPassesCuts(collision)) {
      float cent = collision.cent();
      fHistogramRegistry.fill(HIST("hCentrality"), cent);

      for (auto const& track : tracks) {

        Float_t ptot = track.ptot();
        Float_t dEdx = track.tpcSignal();

        fHistogramRegistry.fill(HIST("hDcaxyPtBefore"), track.pt(), track.dcaxy());
        fHistogramRegistry.fill(HIST("hDcazPtBefore"), track.pt(), track.dcaz());

        fHistogramRegistry.fill(HIST("hDedxPtotBefore"), ptot, dEdx);
        if (trackPassesCuts(track)) {
          float eta = track.eta();

          fHistogramRegistry.fill(HIST("hEtaPhi"), track.phi(), track.eta());
          fHistogramRegistry.fill(HIST("hDedxPtot"), ptot, dEdx);

          fHistogramRegistry.fill(HIST("hDedxPtotEl"), ptot, track.tpcSigmaExpected()[1]);

          // LOGF(info, "expected: %f", track.tpcSigmaExpected()[1]);

          fHistogramRegistry.fill(HIST("hP"), track.p());
          fHistogramRegistry.fill(HIST("hPtot"), ptot);
          fHistogramRegistry.fill(HIST("hChi2ndf"), track.tpcChi2NCl());
          fHistogramRegistry.fill(HIST("hNTPCCrossedRows"), track.tpcNclsCrossedRows());
          fHistogramRegistry.fill(HIST("hDcaxyPt"), track.pt(), track.dcaxy());
          fHistogramRegistry.fill(HIST("hDcazPt"), track.pt(), track.dcaz());

          fHistogramRegistry.fill(HIST("hDedxPtotEtaCent"), ptot, dEdx, eta, 2.5);

        }
      }
    }

  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<netProtonCumulantsHistograms>(cfgc)};
}