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

using namespace o2::framework;
using namespace o2::constants;

std::vector<double> etaBinning = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                                   0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
std::vector<double> centBinning = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

const AxisSpec gAxis_eta{etaBinning, "#eta"};
const AxisSpec gAxis_cent{centBinning, "Centrality (%)"};
const AxisSpec gAxis_phi{400, 0., TwoPI};
const AxisSpec gAxis_dEdx{1024, 0, 1024};
const AxisSpec gAxis_ptot{400, 0., 4.};
const AxisSpec gAxis_chi2ndf{100, 0., 10.};
const AxisSpec gAxis_ncrossedrows{160, 0., 160.};

template<typename TFloat, typename TArray>
size_t getBinIndex(TFloat const& x, TArray const& arr)
{
    for (size_t idx = 0; idx < arr.size() - 1; idx++)
        if (x > arr[idx] && x < arr[idx + 1])
            return idx;
    throw runtime_error(Form("Error in getBinIndex: %f out of range", x));
    return 0;
};