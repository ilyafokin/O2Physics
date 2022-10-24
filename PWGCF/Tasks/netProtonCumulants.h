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
#include "TF1.h"
#include "TString.h"
#include "TColor.h"

#include "boost/multi_array.hpp"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants;
using std::vector, std::map, std::string, std::array;
using std::exp, std::abs, std::erf;

std::vector<double> etaBinning = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                                   0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
std::vector<double> centBinning = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

constexpr int nBinsPtot = 200;

const AxisSpec gAxis_eta{etaBinning, "#eta"};
const AxisSpec gAxis_cent{centBinning, "Centrality (%)"};
const AxisSpec gAxis_phi{400, 0., TwoPI};
const AxisSpec gAxis_dEdx{1024, 0., 1024.};
const AxisSpec gAxis_ptot{nBinsPtot, 0., 4.};
const AxisSpec gAxis_chi2ndf{100, 0., 10.};
const AxisSpec gAxis_ncrossedrows{160, 0., 160.};
const AxisSpec gAxis_dcaxy{100, -1., 1.};
const AxisSpec gAxis_dcaz{100, -4., 4.};
const AxisSpec gAxis_prob{200, -1, 2};
const AxisSpec gAxis_WPr{200, 0, 100};

template<typename TFloat, typename TArray>
size_t getBinIndex(TFloat const& x, TArray const& arr)
{
    for (size_t idx = 0; idx < arr.size() - 1; idx++)
        if (x >= arr[idx] && x < arr[idx + 1])
            return idx;
    throw runtime_error(Form("Error in getBinIndex: %f out of range. Maximum is %f", x, arr[arr.size() - 1]));
    return 0;
};

enum eParticleType : Int_t {
  kEl = 0,
  kPi = 1,
  kKa = 2,
  kPr = 3,
  kBEl = 4,
  kBPi = 5,
  kBKa = 6,
  kBPr = 7,
  kTot = 8
} EPart;

enum eFitPar : Int_t {
    kAmpl = 0,
    kMean = 1,
    kSigm = 2,
    kKurt = 3,
    kSkew = 4,
    kLastFitPar = 5
} EFitPar;

TString fitFunctionGausPar(Int_t start) {
  return Form("[%d]*exp(-0.5*(abs(x-[%d])/[%d])**[%d])*(1+TMath::Erf([%d]*(x-[%d])/[%d]/1.414213))",
              start, start+1, start+2, start+3, start+4, start+1, start+2);
}

array<Int_t, 10> myColors;

void initColors() {
  TColor* cTableau0 = new TColor(TColor::GetFreeColorIndex(), 76, 120, 168);
  TColor* cTableau1 = new TColor(TColor::GetFreeColorIndex(), 245, 133, 24);
  TColor* cTableau2 = new TColor(TColor::GetFreeColorIndex(), 228, 87, 86);
  TColor* cTableau3 = new TColor(TColor::GetFreeColorIndex(), 114, 183, 178);
  TColor* cTableau4 = new TColor(TColor::GetFreeColorIndex(), 84, 162, 75);
  TColor* cTableau5 = new TColor(TColor::GetFreeColorIndex(), 238, 202, 59);
  TColor* cTableau6 = new TColor(TColor::GetFreeColorIndex(), 178, 121, 162);
  TColor* cTableau7 = new TColor(TColor::GetFreeColorIndex(), 255, 157, 166);
  TColor* cTableau8 = new TColor(TColor::GetFreeColorIndex(), 157, 177, 93);
  TColor* cTableau9 = new TColor(TColor::GetFreeColorIndex(), 186, 176, 172);

  myColors[0] = cTableau0->GetNumber();
  myColors[1] = cTableau1->GetNumber();
  myColors[2] = cTableau2->GetNumber();
  myColors[3] = cTableau3->GetNumber();
  myColors[4] = cTableau4->GetNumber();
  myColors[5] = cTableau5->GetNumber();
  myColors[6] = cTableau6->GetNumber();
  myColors[7] = cTableau7->GetNumber();
  myColors[8] = cTableau8->GetNumber();
  myColors[9] = cTableau9->GetNumber();
}

template<typename TArray>
double fitFunctionGaus(double const& x, TArray par) {
  double ret = par[0] * exp(-0.5 * pow(abs(x - par[1]) / par[2], par[3])) * (1 + erf(par[4] * (x - par[1]) / (M_SQRT2 * par[2])));
  return ret;
}