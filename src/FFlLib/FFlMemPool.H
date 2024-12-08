// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_MEM_POOL_H
#define FFL_MEM_POOL_H

#ifdef FT_KERNEL
namespace FTK { class FFlLinkHandler; }
#else
class FFlLinkHandler;
#endif


namespace FFlMemPool
{
#ifdef FT_KERNEL
  using PartHandler = FTK::FFlLinkHandler;
#else
  using PartHandler = FFlLinkHandler;
#endif
  void deleteAllLinkMemPools();
  void deleteVisualsMemPools();
  void setAsMemPoolPart(PartHandler*);
  void freeMemPoolPart(PartHandler*);
  void resetMemPoolPart();
}

#endif
