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

using namespace o2::framework;
using namespace o2::constants;

const AxisSpec gAxis_eta{400, -0.9, 0.9};
const AxisSpec gAxis_phi{400, 0., TwoPI};
const AxisSpec gAxis_dEdx{1024, 0, 1024};
const AxisSpec gAxis_ptot{400, 0., 4.};