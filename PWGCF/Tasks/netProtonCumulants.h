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

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "PWGCF/DataModel/IdentityTables.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>

using namespace o2::framework;
using namespace o2::constants;

std::vector<double> etaBinning = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                                   0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
std::vector<double> centBinning = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

const AxisSpec gAxis_eta{etaBinning, "#eta"};
const AxisSpec gAxis_cent{centBinning, "Centrality (%)"};
const AxisSpec gAxis_phi{400, 0., TwoPI};
const AxisSpec gAxis_dEdx{1024, 0., 1024.};
const AxisSpec gAxis_ptot{400, 0., 4.};
const AxisSpec gAxis_chi2ndf{100, 0., 10.};
const AxisSpec gAxis_ncrossedrows{160, 0., 160.};
const AxisSpec gAxis_dcaxy{100, -1., 1.};
const AxisSpec gAxis_dcaz{100, -4., 4.};

template<typename TFloat, typename TArray>
size_t getBinIndex(TFloat const& x, TArray const& arr)
{
    for (size_t idx = 0; idx < arr.size() - 1; idx++)
        if (x >= arr[idx] && x < arr[idx + 1])
            return idx;
    throw runtime_error(Form("Error in getBinIndex: %f out of range", x));
    return 0;
};
/*
void ReadFitParamsFromTree(TString paramTreeName, Int_t fitIter)
{
  //
  // Read the fit parameters from the paramtree file "ParamTree.root" which comes from PID FITS
  //
  // Open the lookup table and initialize array
  fLineShapesLookUpTable = new TFile(paramTreeName);
  cout << " ReadFitParamsFromTree.Info: Read line shapes from ttree " << endl;
  treeLookUp = (TTree*)fLineShapesLookUpTable->Get(treeLineShapes);
  Int_t nent = treeLookUp -> GetEntries();
  cout << nent <<  " ReadFitParamsFromTree.Info: tree is fine go ahead " << endl;

  Int_t sign      = 0;
  Int_t it        = 0;
  Int_t sl        = 0;
  Int_t syst      = 0;
  Int_t orig      = 0;
  Int_t corr      = -1;
  Int_t myBin[3]  = {0};  // 0; eta, 1;cent, 2;ptot
  // Double_t ptot;
  Float_t eta, cent, ptot;

  Double_t elA  = 0., elM  = 0., elSi = 0., elK = 0., elSk = 0.;
  Double_t piA  = 0., piM  = 0., piSi = 0., piK = 0., piSk = 0.;
  Double_t kaA  = 0., kaM  = 0., kaSi = 0., kaK = 0., kaSk = 0.;
  Double_t prA  = 0., prM  = 0., prSi = 0., prK = 0., prSk = 0.;
  Double_t deA  = 0., deM  = 0., deSi = 0., deK = 0., deSk = 0.;

  TBranch *fMyBinBrach = (TBranch *)treeLookUp->FindBranch("myBin");
  if (fMyBinBrach)
  { // old version of tree format
    treeLookUp->SetBranchAddress("myBin"  ,myBin);
    treeLookUp->SetBranchAddress("it"     ,&it);
  } else
  { // new version of tree format
    treeLookUp->SetBranchAddress("p"    ,&ptot);
    treeLookUp->SetBranchAddress("eta"  ,&eta);
    treeLookUp->SetBranchAddress("cent" ,&cent);
    treeLookUp->SetBranchAddress("it"   ,&it);
  }

  treeLookUp->SetBranchAddress("syst"   ,&syst); // ???
  treeLookUp->SetBranchAddress("corr"   ,&corr); // ???
  treeLookUp->SetBranchAddress("orig"   ,&orig); // ???

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

  for(Int_t i = 0; i < nent; ++i)
  {
    // myBin[0] --> Eta, myBin[1]--> Centrality, myBin[2]-->Momentum
    treeLookUp -> GetEntry(i);
    if (it != fitIter) continue;
    if (corr != fCorr) continue;
    if (!fMyBinBrach){
      myBin[0]=fhEta  -> FindBin(eta + 0.0000001) - 1;
      myBin[1]=fhCent -> FindBin(cent+ 0.0000001) - 1;
      myBin[2]=fhPtot -> FindBin(ptot+ 0.0000001) - 1;
    }
    //
    if (sign==1){
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kDe] = deA;

      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kDe] = deM;

      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kDe] = deSi;

      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kDe] = deSk;

      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kDe] = deK;
    }

    if (sign==-1){
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBDe] = deA;

      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = -elM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = -piM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = -kaM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = -prM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBDe] = -deM;

      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBDe] = deSi;

      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = -elSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = -piSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = -kaSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = -prSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBDe] = -deSk;

      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBDe] = deK;
    }
  }

  Int_t centBinRange[2] = {TMath::Max(fCentInputBin,0), TMath::Min(fCentInputBin+1,fnCentBins)};
  Int_t etaBinRange[2]  = {TMath::Max(fEtaDownBin-1,0), TMath::Min(fEtaUpBin+1,fnEtaBins)};
  // Int_t momBinRange[2]  = {TMath::Max(fpDownBin-1,0), TMath::Min(fpUpBin+1,fnMomBins)}; // TODO
  Int_t momBinRange[2]  = {0, fnMomBins};
  cout << "====================================================" << endl;
  cout << "   Reading the file is started " << endl;
  cout << "   Bin ranges where fits are read from lookup table " << endl;
  cout << "====================================================" << endl;
  cout << "cent Bin window = " << centBinRange[0] << "  " << centBinRange[1] << endl;
  cout << "eta  Bin window = " << etaBinRange[0]  << "  " << etaBinRange[1] << endl;
  cout << "mom  Bin window = " << momBinRange[0]  << "  " << momBinRange[1] << endl;
  cout << "====================================================" << endl;
  //
  // Fill TF1s into TClonesArray
  outFile->cd();
  // for (Int_t i = 0; i<fnEtaBins; i++){
  //   for (Int_t j = 0; j< fnCentBins; j++){
  //     for (Int_t k = 0; k<fnMomBins; k++){
  for (Int_t i = etaBinRange[0]; i<etaBinRange[1]; i++){
    for (Int_t j = centBinRange[0]; j<centBinRange[1]; j++){
      for (Int_t k = momBinRange[0];  k<momBinRange[1]; k++){

        TString objPtotonName     = Form("particle_0_bin_%d_bin_%d_bin_%d",i,j,k);
        TString objAntiProtonName = Form("particle_1_bin_%d_bin_%d_bin_%d",i,j,k);
        TString objOthersName     = Form("particle_2_bin_%d_bin_%d_bin_%d",i,j,k);
        //
        // proton
        fParticles[i][j][k][0] = new TF1("f" + objPtotonName,fitFunctionGausPar(0),fMindEdx,fMaxdEdx); // protons
        fParticles[i][j][k][0]->FixParameter(0,fAmpArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(1,fMeanArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(2,fSigmaArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(3,fKurtosisArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(4,fSkewArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->SetLineColor(kGreen+2);
        fParticles[i][j][k][0]->SetNpx(nBinsLineShape);
        //
        hParticles[i][j][k][0] = (TH1D*)fParticles[i][j][k][0]->GetHistogram()->Clone(objPtotonName);
        // if (fSubsample>0) hParticles[i][j][k][0]->Scale(1./nSubSample);
        hParticles[i][j][k][0]->Scale(1./nSubSample);
        hParticles[i][j][k][0]->SetName(objPtotonName);
        hParticles[i][j][k][0]->SetLineColor(kGreen+2);
        //
        // antiprotons
        fParticles[i][j][k][1] = new TF1("f" + objAntiProtonName,fitFunctionGausPar(0),fMindEdx,fMaxdEdx); // antiprotons
        fParticles[i][j][k][1]->FixParameter(0,fAmpArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(1,fMeanArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(2,fSigmaArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(3,fKurtosisArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(4,fSkewArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->SetLineColor(kYellow+2);
        fParticles[i][j][k][1]->SetNpx(nBinsLineShape);
        //
        hParticles[i][j][k][1] = (TH1D*)fParticles[i][j][k][1]->GetHistogram()->Clone(objAntiProtonName);
        // if (fSubsample>0) hParticles[i][j][k][1]->Scale(1./nSubSample);
        hParticles[i][j][k][1]->Scale(1./nSubSample);
        hParticles[i][j][k][1]->SetName(objAntiProtonName);
        hParticles[i][j][k][1]->SetLineColor(kYellow+2);

        //
        // Other particles as backgtound
        fParticles[i][j][k][2] = MergeOtherTF1s(i,j,k);
        fParticles[i][j][k][2]->SetName("f" + objOthersName);
        fParticles[i][j][k][2]->SetLineColor(kRed+2);
        fParticles[i][j][k][2]->SetNpx(nBinsLineShape);
        //
        hParticles[i][j][k][2] = (TH1D*)fParticles[i][j][k][2]->GetHistogram()->Clone(objOthersName);
        // if (fSubsample>0) hParticles[i][j][k][2]->Scale(1./nSubSample);
        hParticles[i][j][k][2]->Scale(1./nSubSample);
        hParticles[i][j][k][2]->SetName(objOthersName);
        hParticles[i][j][k][2]->SetLineColor(kRed+2);


      }
    }
  }
*/