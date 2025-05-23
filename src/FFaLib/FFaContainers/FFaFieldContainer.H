// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_FIELD_CONTAINER_H
#define FFA_FIELD_CONTAINER_H

#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include "FFaLib/FFaContainers/FFaQueryBase.H"
#include "FFaLib/FFaPatterns/FFaSingelton.H"
#include "FFaLib/FFaTypeCheck/FFaTypeCheck.H"
#include "FFaLib/FFaDynCalls/FFaDynCB.H"

class FFaFieldBase;
class FFaReferenceBase;
class FFaReferenceListBase;


class FFaFieldContainer : public FFaQueryBase
{
  Fmd_BASE_HEADER_INIT();

  using IntVec    = std::vector<int>;
  using FieldMap  = std::map<std::string,FFaFieldBase*>;
  using ObjectMap = std::multimap<std::string,FFaFieldContainer*>;
  using ObjectVec = std::vector<FFaFieldContainer*>;

public:
  FFaFieldContainer();
  virtual ~FFaFieldContainer();

  bool erase();

  // Field browsing

  void getFields(FieldMap& mapToFill) const;
  FFaFieldBase* getField(const std::string& fieldName) const;

  // Topology browsing

  void getReferredObjs (ObjectMap& mapToFill) const;
  void getReferringObjs(ObjectMap& mapToFill) const;
  void getReferringObjs(ObjectVec& vecToFill,
                        const std::string& contextName,
                        bool sortOnId = false) const;

  template<class T>
  void getReferringObjs(std::vector<T*>& vecToFill,
                        const std::string& context = "",
                        bool sortOnId = false) const
  {
    T* ptr;
    if (sortOnId)
    {
      // Return the objects in sorted order w.r.t. their user ID
      ObjectVec allRefs;
      this->getReferringObjs(allRefs,context,sortOnId);
      for (FFaFieldContainer* ref : allRefs)
        if ((ptr = dynamic_cast<T*>(ref)))
          vecToFill.push_back(ptr);
    }
    else
    {
      // Return the objects in arbitrary order
      // (might be different in two runs of the same problem).
      FFaFieldContainer* ref = this->getNext(context,true);
      for (ptr = NULL; ref; ref = this->getNext(context))
        if ((ptr = dynamic_cast<T*>(ref)))
          vecToFill.push_back(ptr);
    }
  }

  // Method for checking if there are any references to this object, of the type
  // of the argument \a ptr. The pointer to the first object is returned if any.
  template<class T>
  bool hasReferringObjs(T*& ptr, const std::string& context = "") const
  {
    FFaFieldContainer* ref = this->getNext(context,true);
    for (ptr = NULL; ref; ref = this->getNext(context))
      if ((ptr = dynamic_cast<T*>(ref)))
        return true;
    return false;
  }

  // Method for directly removing references to this object, and to make the
  // objects referring to this object refer to another object instead
  void releaseReferencesToMe(const std::string& contextName = "",
                             FFaFieldContainer* replacement = NULL);

  using FindCB = FFaDynCB4<FFaFieldContainer*&,int,int,const IntVec&>;

  // Methods for resolving of FFaReference and FFaReferenceList after reading

  void resolve(FindCB& findCB);
  void unresolve();

  // Method for initializing fields after resolve

  virtual void initAfterResolve() {}

  // Method returning the ID number and assemblyIDs used in IO and resolving

  virtual int getResolvedID() const { return -1; }
  virtual void getResolvedAssemblyID(IntVec&) const {}

  // Methods updating the assemblyID in references in this object.

  void updateReferences(int oldAssId, int newAssId);
  void updateReferences(const IntVec& oldAssId, const IntVec& newAssID);

  // Method for reading the value of a named field from an input stream

  bool readField(const std::string& key, std::istream& is,
		 bool datafieldsOnly = false);

  // Methods for copying all field values from another object into this

  bool copy(const FFaFieldContainer* other,
	    bool fieldsOnly = false, bool unresolve = false);
  bool copyFields(const FFaFieldContainer* obj) { return this->copy(obj,true); }

  // Method for resetting fields to their default values

  bool resetFields(const FFaFieldContainer* obj);

protected:
  virtual bool eraseOptions() { return true; }

  // Methods to add guarded pointers to other FFaFieldContainers
  // into the internal book keeping

  void addRef(FFaReferenceBase* ref);
  void addRefList(FFaReferenceListBase* list);

  // Methods to add and remove read/writable fields of "any" type
  // to/from the internal book keeping

  void addField(const std::string& identifier, FFaFieldBase* field);
  void removeField(const std::string& identifier);

  // The read/writable fields

  using FDict   = std::set<std::string>;
  using FDictIt = FDict::const_iterator;

  struct FDictLess
  {
    bool operator()(FDictIt s1, FDictIt s2) const { return *s1 < *s2; }
  };

  using FieldContainerMap  = std::map<FDictIt,FFaFieldBase*,FDictLess>;
  using FieldContainerDict = FFaSingelton<FDict,FFaFieldContainer>;

  FieldContainerMap myFields;

public:
  static void removeDictInstance() { FieldContainerDict::removeInstance(); }

private:
  friend class FFaReferenceBase;

  // Methods to mark/unmark this object as referred to by an object
  // used by FFaReferenceBase only

  void insertInRefBy  (FFaReferenceBase* ref);
  void deleteFromRefBy(FFaReferenceBase* ref);

  // Helper method to loop over all objects referring to this

  FFaFieldContainer* getNext(const std::string& context,
                             bool getFirst = false) const;

  // The internal book keeping of who is pointing to whom

  using ReferenceSet = std::vector<FFaReferenceBase*>;

  ReferenceSet                     myRefBy;
  std::list<FFaReferenceBase*>     myRefTo;
  std::list<FFaReferenceListBase*> myRefLists;
};


#define FFA_FIELD_INIT(field, defaultValue, identifier)			\
  this->addField(identifier, &field);					\
  field.setValue(defaultValue,false);					\
  field.setDefaultValue(field.getValue(),false)

#define FFA_OBSOLETE_FIELD_INIT(field, defaultValue, identifier, obj)	\
  obj->addField(identifier, &field);					\
  field.setValue(defaultValue,false);					\
  field.setDefaultValue(field.getValue(),false)

#define FFA_FIELD_DEFAULT_INIT(field, identifier)	\
  this->addField(identifier, &field)

#define FFA_OBSOLETE_FIELD_DEFAULT_INIT(field, identifier, obj)	\
  obj->addField(identifier, &field)

#define FFA_REFERENCE_INIT(reference)		\
  this->addRef(&reference);			\
  reference.setContextName(#reference)

#define FFA_REFERENCELIST_INIT(referenceList)	\
  this->addRefList(&referenceList);		\
  referenceList.setContextName(#referenceList)

#define FFA_REFERENCE_FIELD_INIT(field, reference, identifier) \
  FFA_FIELD_INIT(field, &reference, identifier);	       \
  FFA_REFERENCE_INIT(reference)

#define FFA_REFERENCELIST_FIELD_INIT(field, referenceList, identifier) \
  FFA_FIELD_INIT(field, &referenceList, identifier);		       \
  FFA_REFERENCELIST_INIT(referenceList)

#define FFA_OBSOLETE_FIELD_REMOVE(identifier, obj)	\
  obj->removeField(identifier)

#endif
