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

struct netProtonCumulants {
  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};
  HistogramSpec hEtaPhi{"hEtaPhi", "hEtaPhi;#eta;#phi;counts", {HistType::kTH2F, {gAxis_eta, gAxis_phi}}};
  HistogramSpec hdEdxPtot{"hdEdxPtot", "hdEdxPtot;TPC dE/dx (a.u.);p_{tot};counts", {HistType::kTH2F, {gAxis_dEdx, gAxis_ptot}}};

  void init(InitContext&){
    LOGF(info, "in init()");
    fHistogramRegistry.add(hEtaPhi);
    fHistogramRegistry.add(hdEdxPtot);
  }

  void process(o2::aod::Tracks const& tracks) {
    for (auto &track : tracks) {
      Float_t eta = track.eta();
      Float_t phi = track.phi();

      // fHistogramRegistry.fill(HIST("hEtaPhi"), {eta, phi});
      // LOGF(info, "eta: %.2f, phi: %.2f", eta, phi);
      fHistogramRegistry.get<TH2F>(HIST("hEtaPhi"))->Fill(eta, phi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<netProtonCumulants>(cfgc)};
}