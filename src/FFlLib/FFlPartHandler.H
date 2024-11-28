// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PART_HANDLER_H
#define FFL_PART_HANDLER_H

#ifdef FT_KERNEL
namespace FTK {
#endif
class FFlLinkHandler;
#ifdef FT_KERNEL
}
using FFlPartHandler = FTK::FFlLinkHandler;
#else
using FFlPartHandler = FFlLinkHandler;
#endif

#endif
