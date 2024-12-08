// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_UTILS_H
#define FFL_UTILS_H

#include <string>
#include <vector>
#include <map>

#ifdef FT_KERNEL
namespace FTK { class FFlLinkHandler; }
#else
class FFlLinkHandler;
#endif
class FaMat34;


namespace FFl
{
#ifdef FT_KERNEL
  using PartHandler = FTK::FFlLinkHandler;
#else
  using PartHandler = FFlLinkHandler;
#endif

  struct DepDOF
  {
    int node;
    int lDof;
    double coeff;
    DepDOF(int n = 1, int d = 1, double c = 1.0) : node(n), lDof(d), coeff(c) {}
  };

  using DepDOFs = std::vector<DepDOF>;
  using MPC     = std::map<short int,DepDOFs>;
  using MPCMap  = std::map<int,MPC>;

  bool convertMPCsToWAVGM(PartHandler* part, const MPCMap& mpcs);

  bool activateElmGroups(PartHandler* part, const std::string& elmGroups = "");

  bool extractBodyFromShell(PartHandler* part, const FaMat34& partCS,
                            const std::string& fname);
}

#endif
