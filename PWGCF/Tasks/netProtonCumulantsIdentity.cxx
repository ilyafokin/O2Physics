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
using std::vector, std::map, std::string;

struct netProtonCumulantsIdentity {
  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};

  void createLineShapesTable(TString fitResultsPath);

  void init(InitContext&){
    LOGF(info, "in init()");
    createLineShapesTable("idMethodInput_fullrange.root");
  }

  void process(aod::IdentityCollisions::iterator const& collision, aod::IdentityTracks const& tracks) {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<netProtonCumulantsIdentity>(cfgc)};
}

void netProtonCumulantsIdentity::createLineShapesTable(TString fitResultsPath) {
  // Open the lookup table and initialize array
  TFile* fitResultsFile = new TFile(fitResultsPath);
  LOGF(info, "createLineShapesTable: Read line shapes from ttree %s", fitResultsPath.Data());
  TTree* treeLookUp = (TTree*) fitResultsFile->Get("treeId");
  Int_t nent = treeLookUp -> GetEntries();
  LOGF(info,"createLineShapesTable: tree is fine (%d entries), go ahead ", nent);
}