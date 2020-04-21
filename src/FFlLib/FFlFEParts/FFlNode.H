// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////
/*!
  \file FFlNode.H
  \brief Finite element node class.
*/

#ifndef FFL_NODE_H
#define FFL_NODE_H

#include "FFlLib/FFlPartBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlReference.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"

class FFlVertex;
class FFlPCOORDSYS;
class FFlAttributeBase;
class FFlFEResultBase;
class FFlFENodeResult;
class FFaUnitCalculator;
class FFaCheckSum;
class FaVec3;


/*!
  \brief Class representing a finite element node.
  \details Stores the position and some additional data for the FE node.
*/

class FFlNode : public FFlPartBase
{
public:
  FFlNode(int ID, double x, double y, double z, int s = 0);
  FFlNode(int ID, const FaVec3& pos, int s = 0);
  FFlNode(const FFlNode& node);
  virtual ~FFlNode();

  FFlNode& operator=(const FFlNode&) = delete;
  FFlNode* clone() const { return new FFlNode(*this); }

  static void init();

  FFL_TYPE_INFO(FFlNode);

  void calculateChecksum(FFaCheckSum* cs, int precision,
                         bool includeExtNodeInfo);
  void convertUnits(const FFaUnitCalculator* convCal);

  // Vertex management:

  int getVertexID() const;
  void setVertex(FFlVertex* aVertex);
  FFlVertex* getVertex() const { return myVertex; }
  const FaVec3& getPos() const;

  // Node status management:

  enum { INTERNAL = 0, EXTERNAL = 1, SLAVENODE = 2, REFNODE = 3 };

  bool setStatus(int newStat);
  bool setExternal(bool ext = true);

  int getStatus(int ignore = 0) const;
  bool isFixed(int dof = 0) const;
  bool isExternal() const   { return myStatus == 1 || myStatus == 13; }
  bool isSlaveNode() const  { return myStatus == 2 || myStatus == 3; }
  bool isRefNode() const    { return myStatus == 3; }
  bool isAttachable() const { return myStatus >= 0 && myStatus != 2; }

  // Local dofs management:

  enum { FFL_ALL_DOFS = -1, FFL_THREE_DOFS = 3, FFL_SIX_DOFS = 6 };

  void pushDOFs(int dofs)          { if (dofs > myDOFCount) myDOFCount = dofs; }
  bool hasDOFs(int dofs = 1) const { return myDOFCount >= dofs; }
  int  getMaxDOFs() const          { return myDOFCount; }

  // Local coordinate system management:

  void setLocalSystem(const FFlPCOORDSYS* coorSys);
  void setLocalSystem(int ID) { myLocalSystem = ID; }

  bool resolveLocalSystem(const std::map<int,FFlAttributeBase*>& possibleCSs,
                          bool suppressErrmsg = false);

  bool hasLocalSystem() const { return myLocalSystem.isResolved(); }
  int  getLocalSystemID() const;
  FFlPCOORDSYS* getLocalSystem() const;

  // Nodal result management:

  bool hasResults() const { return myResults != NULL; }
  FFlFENodeResult* getResults();
  void deleteResults();

protected:
  char myDOFCount; //!< Number of degrees of freedom (DOF) at this node

  //! \brief Status flag identifying what type of node this is.
  //! \details A node can have any of the following statuses:
  //! - 0 : Internal, i.e., condensed out in the superelement reduction
  //! - 1 : External, i.e., retained in the superelement reduction (supernode)
  //! - 2 : Internal slave node, can't be changed to external node
  //! - 3 : Internal WAVGM reference node, can't be changed to external node
  //! - 13 : External translational DOFs, but internal rotational DOFs
  //! - &lt;0 : Internal node with fixtures, binary digits indicate which DOFs
  char myStatus;

  FFlVertex*       myVertex;  //!< Pointer to vertex with nodal position
  FFlFEResultBase* myResults; //!< Pointer to nodal results container

  FFlReference<FFlAttributeBase> myLocalSystem; //!< Solution coordinate system

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlNode> FFlNodeTypeInfoSpec;

#endif
