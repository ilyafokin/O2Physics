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

namespace o2::aod
{
using myCollisions = soa::Join<aod::Collisions, aod::CentRun2V0Ms>;
using myTracks = soa::Join<aod::Tracks, aod::TracksExtra, o2::aod::TracksDCA,
                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
using myCollision = myCollisions::iterator;
using myTrack = myTracks::iterator;
} // namespace o2::aod

struct netProtonCumulants {

  // event cuts
  // Configurable<float> fVertexZMax{"fVertexZMax", 7, "maximum vertex z position"};
  // Configurable<float> fVertexZMin{"fVertexZMin", 0.15, "maximum vertex z position"};

  // track cuts
  Configurable<float> fEtaMax{"fEtaMax", 0.8, "maximum pseudorapidity"};
  Configurable<float> fNTPCClusterMin{"fNTPCClusterMin", 80, "minimum number of TPC clusters"};
  Configurable<float> fChi2NdfMax{"fChi2NdfMax", 2.5, "maximum chi2/ndf"};
  Configurable<float> fEtaMaxLoose{"fEtaMaxLoose", 0.8, "maximum pseudorapidity"};
  Configurable<float> fChi2NdfMaxLoose{"fChi2NdfMaxLoose", 2.5, "maximum chi2/ndf"};
  // Configurable<float> fDcaxyA{"fDcaxyA", 0.0208, "constant term in max dcaxy"}; // A + B * pT^C
  // Configurable<float> fDcaxyB{"fDcaxyB", 0.04, "factor in max dcaxy"};
  // Configurable<float> fDcaxyC{"fDcaxyC", -1.01, "exponent in max dcaxy"};
  // Configurable<float> fDcaxyA{"fDcaxyA", 0.0208, "constant term in max dcaxy"}; // A + B * pT^C
  // Configurable<float> fDcaxyB{"fDcaxyB", 0.04, "factor in max dcaxy"};
  // Configurable<float> fDcaxyC{"fDcaxyC", -1.01, "exponent in max dcaxy"};

  // tree size reduction cuts
  Configurable<float> fMomMinLoose{"fMomMinLoose", 0.25, "minimum momentum"};
  Configurable<float> fMomMaxLoose{"fMomMaxLoose", 5., "maximum momentum"};
  Configurable<float> fMomMin{"fMomMin", 0.25, "minimum momentum"};
  Configurable<float> fMomMax{"fMomMax", 5., "maximum momentum"};

  // apply cuts using filters
  Filter etaFilter = (nabs(aod::track::eta) < fEtaMaxLoose);
  Filter chi2NdfFilter = (aod::track::tpcChi2NCl < fChi2NdfMaxLoose);

  Filter momentumRangeFilter = (aod::track::p > fMomMinLoose) && (aod::track::p < fMomMaxLoose);

  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};
  HistogramSpec hEtaPhi{"hEtaPhi", "hEtaPhi;#phi;#eta;counts", {HistType::kTH2F, {gAxis_phi, gAxis_eta}}};
  HistogramSpec hdEdxPtot{"hDedxPtot", "hDedxPtot;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTH2F, {gAxis_ptot, gAxis_dEdx}}};
  HistogramSpec hdEdxPtotBefore{"hDedxPtotBefore", "hDedxPtotBefore;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTH2F, {gAxis_ptot, gAxis_dEdx}}};

  HistogramSpec hCentrality{"hCentrality", "hCentrality;counts", {HistType::kTH1I, {gAxis_cent}}};
  HistogramSpec hP{"hP", "hP;counts", {HistType::kTH1I, {gAxis_ptot}}};
  HistogramSpec hPtot{"hPtot", "hPtot;counts", {HistType::kTH1I, {gAxis_ptot}}};
  HistogramSpec hChi2ndf{"hChi2ndf", "hChi2ndf;#chi^{2}/ndf;counts", {HistType::kTH1I, {gAxis_chi2ndf}}};
  HistogramSpec hNTPCCrossedRows{"hNTPCCrossedRows", "hNTPCCrossedRows;Num of crossed rows;counts", {HistType::kTH1I, {gAxis_ncrossedrows}}};

  HistogramSpec hdEdxPtotEtaCent{"hDedxPtotEtaCent", "hDedxPtotEtaCent;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTHnI, {gAxis_ptot, gAxis_dEdx, gAxis_eta, gAxis_cent}}};
  std::map<std::map<std::shared_ptr<TH1>, size_t>, size_t> hDedxPtotEtaCent{}; // hDedxPtotEtaCent[iCent][iEta]

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // eta cut
    if (!(fabs(theTrack.eta()) < fEtaMax)) {
      return false;
    }

    // number of tpc clusters
    if (theTrack.tpcNClsCrossedRows() >= fNTPCClusterMin) {
      return false;
    }

    // chi2/ndf
    if (!(theTrack.tpcChi2NCl() < fChi2NdfMax)) {
      return false;
    }

  return true;
  }

  void init(InitContext&){
    LOGF(info, "in init()");
    fHistogramRegistry.add(hEtaPhi);
    fHistogramRegistry.add(hdEdxPtot);
    fHistogramRegistry.add(hdEdxPtotBefore);

    fHistogramRegistry.add(hCentrality);
    fHistogramRegistry.add(hP);
    fHistogramRegistry.add(hPtot);
    fHistogramRegistry.add(hChi2ndf);
    fHistogramRegistry.add(hNTPCCrossedRows);

    fHistogramRegistry.add(hdEdxPtotEtaCent);

    // for (size_t iCent = 0; iCent < centBinning.size(); iCent++) {
    //   hDedxPtotEtaCent[iCent] = std::map<std::shared_ptr<TH1>, size_t>;
    //   for (size_t iEta = 0; iEta < etaBinning.size(); iEta++) {
    //     hDedxPtotEtaCent[iCent][iEta] = nullptr;
    //     // hCharmProtonKstarDistr[iCharmPart] = registry.add<TH1>(Form("f%sProtonKstarDistr", charmParticleNames[iCharmPart].data()), Form("#it{k}* distribution of triggered p#minus%s pairs;#it{k}* (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
    //   }
    // }

  // HistogramSpec hdEdxPtotEtaCent{"hDedxPtotEtaCent", "hDedxPtotEtaCent;p_{tot};TPC dE/dx (a.u.);counts", {HistType::kTHnI, {gAxis_ptot, gAxis_dEdx, gAxis_eta, gAxis_cent}}};
  }

  using myFilteredTracks = soa::Filtered<aod::myTracks>;
  using myFilteredTrack = myFilteredTracks::iterator;

  void process(aod::myCollision const& collision, myFilteredTracks const& tracks) {

    fHistogramRegistry.fill(HIST("hCentrality"), collision.centRun2V0M());

    for (auto const& track : tracks) {

      Float_t ptot, dEdx;
      ptot = track.tpcInnerParam();
      dEdx = track.tpcSignal();

      // LOGF(info, "crossed rows: %d", track.tpcNClsCrossedRows());

      fHistogramRegistry.fill(HIST("hDedxPtotBefore"), ptot, dEdx);
      if (trackPassesCuts(track)) {
        float eta = track.eta();

        // fHistogramRegistry.fill(HIST("hEtaPhi"), {eta, phi});
        // LOGF(info, "eta: %.2f, phi: %.2f", eta, phi);
        fHistogramRegistry.fill(HIST("hEtaPhi"), track.phi(), track.eta());
        fHistogramRegistry.fill(HIST("hDedxPtot"), ptot, dEdx);

        fHistogramRegistry.fill(HIST("hP"), track.p());
        fHistogramRegistry.fill(HIST("hPtot"), ptot);
        fHistogramRegistry.fill(HIST("hChi2ndf"), track.tpcChi2NCl());
        fHistogramRegistry.fill(HIST("hNTPCCrossedRows"), track.tpcNClsCrossedRows());

        fHistogramRegistry.fill(HIST("hDedxPtotEtaCent"), ptot, dEdx, eta, 2.5);

        // std::shared_ptr<TH2D> hDedx = fHistogramRegistry.get("hDedxPtot");

      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<netProtonCumulants>(cfgc)};
}