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

struct netProtonCumulantsIdentity {
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

  //lineShapesTable[etaBin][centBin][ptotBin][part](dEdx)
  // boost::multi_array<float, 5> fitParamsTable{boost::extents[etaBinning.size()][centBinning.size()][gAxis_ptot.getNbins()][eParticleType::kTot][eFitPar::kLastFitPar]};
  // boost::multi_array<TF1*, 4> lineShapesTable{boost::extents[etaBinning.size()][centBinning.size()][gAxis_ptot.getNbins()][eParticleType::kTot]};
  boost::multi_array<float, 5> fitParamsTable{boost::extents[etaBinning.size()][centBinning.size()][nBinsPtot][eParticleType::kTot][eFitPar::kLastFitPar]};
  boost::multi_array<TF1*, 4> lineShapesTable{boost::extents[etaBinning.size()][centBinning.size()][nBinsPtot][eParticleType::kTot]};
  // boost::multi_array<TF1*, 4> pdfTable{boost::extents[etaBinning.size()][centBinning.size()][gAxis_ptot.getNbins()][eParticleType::kTot]};

  Configurable<string> fFitResultsPath{"fFitResultsPath", "idMethodInput.root", "path of file containing fit results tree"};

  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};
  HistogramSpec hProbabilities{"hProbabilities", "hProbabilities;counts", {HistType::kTH1I, {gAxis_prob}}};
  HistogramSpec hWPr{"hWPr", "hWPr;counts", {HistType::kTH1I, {gAxis_WPr}}};

  HistogramSpec hdEdx{"hdEdx", "hdEdx;ptot;tpcSignal", {HistType::kTH2I, {gAxis_ptot, gAxis_dEdx}}};
  HistogramSpec hdEdxPr{"hdEdxPr", "hdEdxPr;tpcSignal;counts", {HistType::kTH2I, {gAxis_dEdx}}};

  void createLineShapesTable(string fitResultsPath);
  array<double, eParticleType::kTot> getProbabilities(double const& eta, double const& cent, double const& ptot, double const& dEdx);

  void init(InitContext&){
    LOGF(info, "in init()");
    initColors();
    createLineShapesTable(fFitResultsPath);
    fHistogramRegistry.add(hProbabilities);
    fHistogramRegistry.add(hWPr);
    fHistogramRegistry.add(hdEdx);
    fHistogramRegistry.add(hdEdxPr);
  }

  void process(aod::IdentityCollisions::iterator const& collision, aod::IdentityTracks const& tracks) {
    if (collisionPassesCuts(collision) && collision.cent() > 40. && collision.cent() < 50.) {
      array<float, kTot> Ws;
      for (int iPart = 0; iPart < kTot; iPart++) Ws[iPart] = 0.;
      for (auto const& track: tracks) {
        if (!trackPassesCuts(track)) continue;
        if (abs(track.eta()) < 0.8 && track.ptot() > 0.6 && track.ptot() < 1.5) {
          if (track.ptot() > 0.6 && track.ptot() <= 0.62) {
            fHistogramRegistry.fill(HIST("hdEdxPr"), track.tpcSignal());
          }
          const auto probs = getProbabilities(track.eta(), collision.cent(), track.ptot(), track.tpcSignal());
          // for (int iPart = 0; iPart < kTot; iPart++) {
          //   LOGF(info, "prob %d: %f", iPart, probs[iPart]);
          // }
          fHistogramRegistry.fill(HIST("hProbabilities"), probs[kPr]);
          for (int iPart = 0; iPart < kTot; iPart++) {
            Ws[iPart] += probs[iPart];
          }
          if (isnan(Ws[kPr]))
            LOGF(info, "eta %f, cent %f, ptot %f, sign %d, dEdx %f", track.eta(), collision.cent(), track.ptot(), track.sign(), track.tpcSignal());
          fHistogramRegistry.fill(HIST("hdEdx"), track.ptot(), track.tpcSignal());
        }
      }
      fHistogramRegistry.fill(HIST("hWPr"), Ws[kPr]);
      LOGF(info, "W(Pr): %f", Ws[kPr]);
    }
  }

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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<netProtonCumulantsIdentity>(cfgc)};
}

void netProtonCumulantsIdentity::createLineShapesTable(string fitResultsPath) {
  // Open the lookup table and initialize array
  TFile* fitResultsFile = new TFile(fitResultsPath.c_str());
  LOGF(info, "createLineShapesTable: Read line shapes from ttree %s", fitResultsPath.c_str());
  TTree* treeLookUp = (TTree*) fitResultsFile->Get("treeId");
  Int_t nent = treeLookUp -> GetEntries();
  LOGF(info, "createLineShapesTable: tree is fine (%d entries), go ahead ", nent);
  LOGF(info, "building lookuptable: %d x %d x %d x %d = %d functions", lineShapesTable.shape()[0],
                                                                       lineShapesTable.shape()[1],
                                                                       lineShapesTable.shape()[2],
                                                                       lineShapesTable.shape()[3],
                                                                       lineShapesTable.shape()[4]);

  Int_t sign      = 0;
  Int_t it        = 0;
  Int_t sl        = 0;
  Int_t syst      = 0;
  Int_t orig      = 0;
  Int_t corr      = -1;
  Float_t eta, cent, ptot;

  Double_t elA  = 0., elM  = 0., elSi = 0., elK = 0., elSk = 0.;
  Double_t piA  = 0., piM  = 0., piSi = 0., piK = 0., piSk = 0.;
  Double_t kaA  = 0., kaM  = 0., kaSi = 0., kaK = 0., kaSk = 0.;
  Double_t prA  = 0., prM  = 0., prSi = 0., prK = 0., prSk = 0.;
  Double_t deA  = 0., deM  = 0., deSi = 0., deK = 0., deSk = 0.;

  treeLookUp->SetBranchAddress("p"    ,&ptot);
  treeLookUp->SetBranchAddress("eta"  ,&eta);
  treeLookUp->SetBranchAddress("cent" ,&cent);
  treeLookUp->SetBranchAddress("it"   ,&it);

  treeLookUp->SetBranchAddress("syst"   ,&syst);
  treeLookUp->SetBranchAddress("orig"   ,&orig);

  treeLookUp->SetBranchAddress("sign"   ,&sign);
  treeLookUp->SetBranchAddress("sl"     ,&sl);

  treeLookUp->SetBranchAddress("elM" ,&elM);
  treeLookUp->SetBranchAddress("piM" ,&piM);
  treeLookUp->SetBranchAddress("kaM" ,&kaM);
  treeLookUp->SetBranchAddress("prM" ,&prM);
  treeLookUp->SetBranchAddress("deM" ,&deM);

  treeLookUp->SetBranchAddress("elSi" ,&elSi);
  treeLookUp->SetBranchAddress("piSi" ,&piSi);
  treeLookUp->SetBranchAddress("kaSi" ,&kaSi);
  treeLookUp->SetBranchAddress("prSi" ,&prSi);
  treeLookUp->SetBranchAddress("deSi" ,&deSi);

  treeLookUp->SetBranchAddress("elA" ,&elA);
  treeLookUp->SetBranchAddress("piA" ,&piA);
  treeLookUp->SetBranchAddress("kaA" ,&kaA);
  treeLookUp->SetBranchAddress("prA" ,&prA);
  treeLookUp->SetBranchAddress("deA" ,&deA);

  treeLookUp->SetBranchAddress("elSk" ,&elSk);
  treeLookUp->SetBranchAddress("piSk" ,&piSk);
  treeLookUp->SetBranchAddress("kaSk" ,&kaSk);
  treeLookUp->SetBranchAddress("prSk" ,&prSk);
  treeLookUp->SetBranchAddress("deSk" ,&deSk);

  treeLookUp->SetBranchAddress("elK" ,&elK);
  treeLookUp->SetBranchAddress("piK" ,&piK);
  treeLookUp->SetBranchAddress("kaK" ,&kaK);
  treeLookUp->SetBranchAddress("prK" ,&prK);
  treeLookUp->SetBranchAddress("deK" ,&deK);


  // tempFitPar[part][ipar]
  boost::multi_array<float, 2> tempFitPar{boost::extents[eParticleType::kTot][eFitPar::kLastFitPar]};
  for (int iPart = 0; iPart < tempFitPar.shape()[0]; iPart++) {
    for (int iPar = 0; iPar < tempFitPar.shape()[1]; iPar++) {
      tempFitPar[iPart][iPar] = 0.;
    }
  }

  // fill table of fit parameters using the tree
  for(int i = 0; i < nent; i++) {
    treeLookUp->GetEntry(i);

    if (it != 1) continue;
    if (orig != 0) continue;

    for (int iPart = 0; iPart < tempFitPar.shape()[0]; iPart++) {
      for (int iPar = 0; iPar < tempFitPar.shape()[1]; iPar++) {
        tempFitPar[iPart][iPar] = 0.;
      }
    }
    
    if (sign > 0) {
      tempFitPar[kEl][kAmpl] = elA; tempFitPar[kEl][kMean] = elM; tempFitPar[kEl][kSigm] = elSi; tempFitPar[kEl][kSkew] = elSk; tempFitPar[kEl][kKurt] = elK;
      tempFitPar[kPi][kAmpl] = piA; tempFitPar[kPi][kMean] = piM; tempFitPar[kPi][kSigm] = piSi; tempFitPar[kPi][kSkew] = piSk; tempFitPar[kPi][kKurt] = piK;
      tempFitPar[kKa][kAmpl] = kaA; tempFitPar[kKa][kMean] = kaM; tempFitPar[kKa][kSigm] = kaSi; tempFitPar[kKa][kSkew] = kaSk; tempFitPar[kKa][kKurt] = kaK;
      tempFitPar[kPr][kAmpl] = prA; tempFitPar[kPr][kMean] = prM; tempFitPar[kPr][kSigm] = prSi; tempFitPar[kPr][kSkew] = prSk; tempFitPar[kPr][kKurt] = prK;
    } else {
      tempFitPar[kBEl][kAmpl] = elA; tempFitPar[kBEl][kMean] = -elM; tempFitPar[kBEl][kSigm] = elSi; tempFitPar[kBEl][kSkew] = -elSk; tempFitPar[kBEl][kKurt] = elK;
      tempFitPar[kBPi][kAmpl] = piA; tempFitPar[kBPi][kMean] = -piM; tempFitPar[kBPi][kSigm] = piSi; tempFitPar[kBPi][kSkew] = -piSk; tempFitPar[kBPi][kKurt] = piK;
      tempFitPar[kBKa][kAmpl] = kaA; tempFitPar[kBKa][kMean] = -kaM; tempFitPar[kBKa][kSigm] = kaSi; tempFitPar[kBKa][kSkew] = -kaSk; tempFitPar[kBKa][kKurt] = kaK;
      tempFitPar[kBPr][kAmpl] = prA; tempFitPar[kBPr][kMean] = -prM; tempFitPar[kBPr][kSigm] = prSi; tempFitPar[kBPr][kSkew] = -prSk; tempFitPar[kBPr][kKurt] = prK;
    }

    LOGF(info, "eta %f, cent %f, ptot %f, sign %d", eta, cent, ptot, sign);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kEl][kAmpl], tempFitPar[kEl][kMean], tempFitPar[kEl][kSigm], tempFitPar[kEl][kSigm], tempFitPar[kEl][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kPi][kAmpl], tempFitPar[kPi][kMean], tempFitPar[kPi][kSigm], tempFitPar[kPi][kSigm], tempFitPar[kPi][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kKa][kAmpl], tempFitPar[kKa][kMean], tempFitPar[kKa][kSigm], tempFitPar[kKa][kSigm], tempFitPar[kKa][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kPr][kAmpl], tempFitPar[kPr][kMean], tempFitPar[kPr][kSigm], tempFitPar[kPr][kSigm], tempFitPar[kPr][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kBEl][kAmpl], tempFitPar[kBEl][kMean], tempFitPar[kBEl][kSigm], tempFitPar[kBEl][kSigm], tempFitPar[kBEl][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kBPi][kAmpl], tempFitPar[kBPi][kMean], tempFitPar[kBPi][kSigm], tempFitPar[kBPi][kSigm], tempFitPar[kBPi][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kBKa][kAmpl], tempFitPar[kBKa][kMean], tempFitPar[kBKa][kSigm], tempFitPar[kBKa][kSigm], tempFitPar[kBKa][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", tempFitPar[kBPr][kAmpl], tempFitPar[kBPr][kMean], tempFitPar[kBPr][kSigm], tempFitPar[kBPr][kSigm], tempFitPar[kBPr][kKurt]);

    for (int iPar = 0; iPar < fitParamsTable.shape()[4]; iPar++) {
      if (sign > 0) {
        for (int iPart = kEl; iPart <= kPr; iPart++) {
          fitParamsTable[getBinIndex(eta + 1e-5, etaBinning)]
                        [getBinIndex(cent + 1e-5, centBinning)]
                        [getBinIndex(ptot + 1e-5, gAxis_ptot.binEdges)]
                        [iPart][iPar] = tempFitPar[iPart][iPar];
        }
      } else {
        for (int iPart = kBEl; iPart <= kBPr; iPart++) {
          fitParamsTable[getBinIndex(eta + 1e-5, etaBinning)]
                        [getBinIndex(cent + 1e-5, centBinning)]
                        [getBinIndex(ptot + 1e-5, gAxis_ptot.binEdges)]
                        [iPart][iPar] = tempFitPar[iPart][iPar];
        }
      }
    }
  }

  LOGF(info, "done reading the tree");

  // build functions from fit parameters
  // for (int iEta = 0; iEta < lineShapesTable.shape()[0]; iEta++) {
  //   for (int iCent = 0; iCent < lineShapesTable.shape()[1]; iCent++) {
  //     for (int iPtot = 0; iPtot < lineShapesTable.shape()[2]; iPtot++) {
  //       for (int iPart = 0; iPart < lineShapesTable.shape()[3]; iPart++) {
  //         TString histName = Form("particle_%d_bin_%d_bin_%d_bin_%d", iPart, iEta, iCent, iPtot);
  //         lineShapesTable[iEta][iCent][iPtot][iPart] = new TF1("f" + histName, fitFunctionGausPar(0), gAxis_dEdx.binEdges[0], gAxis_dEdx.binEdges[gAxis_dEdx.getNbins()]);
  //         // lineShapesTable[iEta][iCent][iPtot][iPart]->SetNpx(4 * gAxis_dEdx.getNbins() / 2);
  //         lineShapesTable[iEta][iCent][iPtot][iPart]->SetLineColor(myColors[iPart]);
  //         for (int iPar = 0; iPar < fitParamsTable.shape()[4]; iPar++) {
  //           lineShapesTable[iEta][iCent][iPtot][iPart]->SetParameter(iPar, fitParamsTable[iEta][iCent][iPtot][iPart][iPar]);
  //         }
  //       }
  //     }
  //   }
  // }
}

array<double, eParticleType::kTot> netProtonCumulantsIdentity::getProbabilities(double const& eta, double const& cent, double const& ptot, double const& dEdx) {
  array<double, eParticleType::kTot> ret;
  double divisor = 0.;
  for (int i = 0; i < ret.size(); i++) ret.at(i) = 0.;

  const int etaBin = getBinIndex(eta, etaBinning);
  const int centBin = getBinIndex(cent, centBinning);
  const int ptotBin = getBinIndex(ptot, gAxis_ptot.binEdges);

  for (int iPart = 0; iPart < eParticleType::kTot; iPart++) {
    ret.at(iPart) = fitFunctionGaus(dEdx, fitParamsTable[etaBin][centBin][ptotBin][iPart]);
    if (isnan(ret.at(iPart))) {
      LOGF(info, "ret: %f, %f, %f, %f, %f", ret[0],
                                      ret[1],
                                      ret[2],
                                      ret[3],
                                      ret[4]
      );
      LOGF(info, "eta %f, cent %f, ptot %f", eta, cent, ptot);
    }
  }
  for (int iPart = 0; iPart < eParticleType::kTot; iPart++) {
    divisor += fitFunctionGaus(dEdx, fitParamsTable[etaBin][centBin][ptotBin][iPart]);
  }

  if (divisor < 1e-20) {
    LOGF(info, "divisor: %f", divisor);
    printf("ret: %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f\n", ret[0],
                                    ret[1],
                                    ret[2],
                                    ret[3],
                                    ret[4],
                                    ret[5],
                                    ret[6],
                                    ret[7],
                                    ret[8]
    );
    LOGF(info, "eta %f, cent %f, ptot %f, dEdx %f", eta, cent, ptot, dEdx);

    LOGF(info, "fit parameters:");

    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kEl][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kEl][kMean], fitParamsTable[etaBin][centBin][ptotBin][kEl][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kEl][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kEl][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kPi][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kPi][kMean], fitParamsTable[etaBin][centBin][ptotBin][kPi][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kPi][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kPi][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kKa][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kKa][kMean], fitParamsTable[etaBin][centBin][ptotBin][kKa][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kKa][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kKa][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kPr][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kPr][kMean], fitParamsTable[etaBin][centBin][ptotBin][kPr][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kPr][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kPr][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kBEl][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kBEl][kMean], fitParamsTable[etaBin][centBin][ptotBin][kBEl][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBEl][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBEl][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kBPi][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kBPi][kMean], fitParamsTable[etaBin][centBin][ptotBin][kBPi][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBPi][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBPi][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kBKa][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kBKa][kMean], fitParamsTable[etaBin][centBin][ptotBin][kBKa][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBKa][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBKa][kKurt]);
    LOGF(info, "%f, %f, %f, %f, %f", fitParamsTable[etaBin][centBin][ptotBin][kBPr][kAmpl], fitParamsTable[etaBin][centBin][ptotBin][kBPr][kMean], fitParamsTable[etaBin][centBin][ptotBin][kBPr][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBPr][kSigm], fitParamsTable[etaBin][centBin][ptotBin][kBPr][kKurt]);

  }


  for (int iPart = 0; iPart < eParticleType::kTot; iPart++) {
    ret.at(iPart) = ret.at(iPart) / divisor;
  }

  return ret;
}