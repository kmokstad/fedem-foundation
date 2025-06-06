// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PNSM_H
#define FFL_PNSM_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPNSM : public FFlAttributeBase
{
public:
  FFlPNSM(int ID);
  FFlPNSM(const FFlPNSM& obj);
  virtual ~FFlPNSM() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPNSM);
  FFL_TYPE_INFO(FFlPNSM);

  FFlField<double> NSM;
  FFlField<int> isShell; // used to determine the proper unit conversion
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
