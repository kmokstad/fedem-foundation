// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPSTRC.H"

#ifdef FT_KERNEL
namespace FTK {
#endif


FFlPSTRC::FFlPSTRC(int id) : FFlAttributeBase(id)
{
  this->addField(name);
}


FFlPSTRC::FFlPSTRC(const FFlPSTRC& obj) : FFlAttributeBase(obj)
{
  this->addField(name);
  name = obj.name;
}


void FFlPSTRC::init()
{
  using TypeInfoSpec  = FFaSingelton<FFlTypeInfoSpec,FFlPSTRC>;
  using AttributeSpec = FFaSingelton<FFlFEAttributeSpec,FFlPSTRC>;

  TypeInfoSpec::instance()->setTypeName("PSTRC");
  TypeInfoSpec::instance()->setDescription("Strain coat properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_PROP);

  AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  AttributeSpec::instance()->addLegalAttribute("PTHICKREF",false);
  AttributeSpec::instance()->addLegalAttribute("PHEIGHT",false);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPSTRC::create,int,FFlAttributeBase*&));
}

#ifdef FT_KERNEL
} // namespace FTK
#endif
