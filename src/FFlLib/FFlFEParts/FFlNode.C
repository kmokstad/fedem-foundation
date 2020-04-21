// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlVertex.H"
#include "FFlLib/FFlFEParts/FFlPCOORDSYS.H"
#include "FFlLib/FFlFEResultBase.H"

#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


FFlNode::FFlNode(int id, double x, double y, double z, int s) : FFlPartBase(id)
{
  myStatus = s;
  myDOFCount = 0;

  myVertex = new FFlVertex(x,y,z);
  myVertex->ref();
  myVertex->setNode(this);

  myResults = NULL;
}


FFlNode::FFlNode(int id, const FaVec3& pos, int s) : FFlPartBase(id)
{
  myStatus = s;
  myDOFCount = 0;

  myVertex = new FFlVertex(pos);
  myVertex->ref();
  myVertex->setNode(this);

  myResults = NULL;
}


FFlNode::FFlNode(const FFlNode& otherNode) : FFlPartBase(otherNode)
{
  myStatus = otherNode.myStatus;
  myDOFCount = otherNode.myDOFCount;

  if (otherNode.myVertex)
  {
    myVertex = new FFlVertex(*otherNode.myVertex);
    myVertex->ref();
    myVertex->setNode(this);
  }
  else
    myVertex = NULL;

  myResults = NULL;
}


FFlNode::~FFlNode()
{
  if (myVertex)
  {
    myVertex->setNode(NULL);
    myVertex->unRef();
  }
  this->deleteResults();
}


void FFlNode::init()
{
  FFlNodeTypeInfoSpec::instance()->setTypeName("Node");
  FFlNodeTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::NODE);
}


void FFlNode::calculateChecksum(FFaCheckSum* cs, int precision,
                                bool includeExtNodeInfo)
{
  FFlPartBase::checksum(cs);

  if (myVertex)
    cs->add(*myVertex,precision);

  cs->add(includeExtNodeInfo && myStatus == 1);
  if (myStatus < 0 || (includeExtNodeInfo && myStatus == 13))
    cs->add((int)myStatus);

  int localCS = myLocalSystem.getID();
  if (localCS > 0) cs->add(localCS);
}


void FFlNode::convertUnits(const FFaUnitCalculator* convCal)
{
  if (myVertex) convCal->convert(*myVertex, "LENGTH");
}


int FFlNode::getVertexID() const
{
  if (myVertex) return myVertex->getRunningID();

  std::cerr <<"FFlNode::getVertexID(): No vertex set up for node "
            << this->getID() << std::endl;

  return -1;
}


const FaVec3& FFlNode::getPos() const
{
  if (myVertex) return *myVertex;

  std::cerr <<"FFlNode::getPos(): No vertex set up for node "
            << this->getID() << std::endl;

  static FaVec3 dummy;
  return dummy;
}


void FFlNode::setVertex(FFlVertex* aVertex)
{
  if (myVertex)
  {
    myVertex->setNode(NULL);
    myVertex->unRef();
  }
  myVertex = aVertex;
  myVertex->ref();
  myVertex->setNode(this);
}


bool FFlNode::setStatus(int newStat)
{
  if (myStatus == newStat)
    return false;

  myStatus = newStat;
  return true;
}


int FFlNode::getStatus(int ignore) const
{
  return myStatus > ignore ? myStatus : 0;
}


bool FFlNode::setExternal(bool ext)
{
  if (myStatus == 2 || myStatus == 3 || myStatus == (int)ext)
    return false;
  else if (myStatus == 13 && ext)
    return false;

  myStatus = ext;
  return true;
}


bool FFlNode::isFixed(int dof) const
{
  if (myStatus < 0)
    switch (dof) {
    case 1: return -myStatus & 1;
    case 2: return -myStatus & 2;
    case 3: return -myStatus & 4;
    case 4: return -myStatus & 8;
    case 5: return -myStatus & 16;
    case 6: return -myStatus & 32;
    default: return true;
    }

  return false;
}


void FFlNode::setLocalSystem(const FFlPCOORDSYS* coordSys)
{
  myLocalSystem = coordSys;
}


bool FFlNode::resolveLocalSystem(const std::map<int,FFlAttributeBase*>& possibleCSs,
                                 bool suppressErrmsg)
{
  if (myLocalSystem.resolve(possibleCSs))
  {
    FFlAttributeBase* localCS = myLocalSystem.getReference();
    if (!localCS)
      return true; // no local coordinate system for this node
    else if (localCS->getTypeName() == "PCOORDSYS")
      return true; // we have local system

    // Logic error, indicates programming error
    myLocalSystem = (FFlAttributeBase*)NULL;
    std::cerr <<"FFlNode::resolveSystem(): Invalid attribute type provided "
              << localCS->getTypeName() << std::endl;
  }
  else if (!suppressErrmsg)
    ListUI <<"\n *** Error: Failed to resolve PCOORDSYS "
           << myLocalSystem.getID() <<"\n";

  return false;
}


int FFlNode::getLocalSystemID() const
{
  return myLocalSystem.getID();
}


FFlPCOORDSYS* FFlNode::getLocalSystem() const
{
  if (myLocalSystem.isResolved())
    return static_cast<FFlPCOORDSYS*>(myLocalSystem.getReference());
  else
    return NULL;
}


void FFlNode::deleteResults()
{
  delete myResults;
  myResults = NULL;
}
