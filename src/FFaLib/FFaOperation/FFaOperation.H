// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaOperation.H
  \brief Classes that make up the core of a operation/transformation system.

  \details An operation is a kind of function with other functions as arguments.
  Trees are built by setting operations as arguments to other operations.

  \author Jacob Storen
*/

#ifndef FFA_OPERATION_H
#define FFA_OPERATION_H

#include <map>
#include <vector>
#include <string>

#include "FFaLib/FFaDynCalls/FFaDynCB.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#include "FFaLib/FFaPatterns/FFaSingelton.H"
#endif


/*!
  \brief Base class for all operations.
  \details Operations are reference counting each other,
  and deleted when their reference count reaches zero.
*/

class FFaOperationBase
{
protected:
  //! \brief Default constructor.
  FFaOperationBase() { myRefCount = 0; }
  //! \brief Empty destructor.
  virtual ~FFaOperationBase() {}

public:
  //! \brief Increments the reference counter.
  void ref() { myRefCount++; }
  //! \brief Decrements the reference counter, and deletes *this if zero.
  void unref() { if (--myRefCount <= 0) delete this; }

private:
  short int myRefCount; //!< Reference counter
};


/*!
  \brief The base class used as argument in other operations.
*/

template <class RetType>
class FFaOperation : public FFaOperationBase
{
protected:
  //! \brief Default constructor.
  FFaOperation() : IAmEvaluated(false) {}
  //! \brief Empty destructor.
  virtual ~FFaOperation() {}

public:
  //! \brief Returns current value of this operation.
  bool invoke(RetType& value)
  {
    bool ok = true;
    if (IAmEvaluated)
      value = myCachedValue;
    else
    {
      ok = this->evaluate(value);
      myCachedValue = value;
      IAmEvaluated = true;
    }

    return ok;
  }

  //! \brief Invalidates the cached value of this operation.
  virtual void invalidate() { IAmEvaluated = false; }

  //! \brief Checks if this operation will return any data.
  virtual bool hasData() const = 0;
  //! \brief Invokes the actual operation.
  virtual bool evaluate(RetType& value) = 0;

private:
  RetType myCachedValue; //!< Cached operation value, to avoid re-evaluation

protected:
  bool IAmEvaluated; //!< If \e true, the cached value will be used
};


/*!
  \brief A class used to manage and store callbacks of a certain type.
  \details This must be a base class of all operations that are supposed
  to have changeable calculations. The type of the callback is given through
  the template argument \a CB.
*/

template <class CB>
class FFaOperationContainer
{
  //! \brief Convenience type definition, name-to-operation mapping.
  using OperCBMap = std::map<std::string,CB>;

public:
  //! \brief Adds a named operation to the operation container.
  static void addOperation(const std::string& name, CB operation)
  {
    static OperCBMap ops;
    ourOperations = &ops;
    if (!name.empty())
      ops[name] = operation;
  }

  //! \brief Returns a list of all operations in this container.
  static void getOperations(std::vector<std::string>& names)
  {
    names.clear();
    if (ourOperations)
    {
      names.reserve(ourOperations->size());
      for (const typename OperCBMap::value_type& oper : *ourOperations)
        names.push_back(oper.first);
    }
  }

  //! \brief Defines the default operation.
  static void setDefaultOperation(const std::string& o) { ourDefaultOper = o; }
  //! \brief Returns the default operation.
  static const std::string& getDefaultOperation() { return ourDefaultOper; }

protected:
  //! \brief The default constructor defines the default operation.
  FFaOperationContainer() { this->selectOperation(ourDefaultOper); }

  //! \brief Selects the named operation to use.
  //! \details The given \a name has to be in the \ref ourOperations container.
  //! If not, this operation will be undefined and
  //! when invoked it will not change the output argument.
  void selectOperation(const std::string& name)
  {
    if (ourOperations && !name.empty())
    {
      typename OperCBMap::const_iterator it = ourOperations->find(name);
      if (it != ourOperations->end())
      {
        myOperationCB = it->second;
        return;
      }
    }

    // The selected operation is not defined
    myOperationCB.erase();
  }

  CB myOperationCB; //!< The selected operation

private:
  //! Container of all operations of the type \a CB.
  static OperCBMap*  ourOperations;
  //! Name of the default operation for this type.
  static std::string ourDefaultOper;
};


/*!
  \brief A class for scalar operations.
*/

template <class RetType, class PrmType = RetType>
class FFaUnaryOp : public FFaOperation<RetType>,
                   public FFaOperationContainer< FFaDynCB2<RetType&,const PrmType&> >
{
public:
  //! \brief The constructor selects the actual operation to use from \a name.
  FFaUnaryOp(FFaOperation<PrmType>* param, const std::string& name = "")
  {
    if ((myParam = param))
      myParam->ref();

    if (!name.empty())
      this->selectOperation(name);
  }

  //! \brief Invalidates the cached value of this operation.
  virtual void invalidate()
  {
    if (this->IAmEvaluated && myParam)
      myParam->invalidate();

    this->IAmEvaluated = false;
  }

  //! \brief Checks if this operation will return any data.
  virtual bool hasData() const
  {
    return myParam ? myParam->hasData() : false;
  }

  //! \brief Invokes the actual operation.
  virtual bool evaluate(RetType& value)
  {
    bool ok = false;
    if (myParam && !this->myOperationCB.empty())
    {
      PrmType param;
      ok = myParam->invoke(param);
      this->myOperationCB.invoke(value,param);
    }

    return ok;
  }

protected:
  //! \brief The destructor deletes the child operation, if any.
  virtual ~FFaUnaryOp()
  {
    if (myParam) myParam->unref();
  }

private:
  FFaOperation<PrmType>* myParam; //!< Parameter child operation

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


/*!
  \brief A class for vector operations.
*/

template <class RetType>
class FFaNToOneOp : public FFaOperation<RetType>,
                    public FFaOperationContainer< FFaDynCB2<RetType&,const std::vector<RetType>&> >
{
public:
  //! \brief The constructor selects the actual operation to use from \a name.
  FFaNToOneOp(const std::vector<FFaOperation<RetType>*>& params,
              const std::string& name = "")
  {
    myParams = params;
    for (FFaOperation<RetType>* p : myParams)
      if (p) p->ref();

    if (!name.empty())
      this->selectOperation(name);
  }

  //! \brief Invalidates the cached value of this operation.
  virtual void invalidate()
  {
    if (this->IAmEvaluated)
      for (FFaOperation<RetType>* p : myParams)
        if (p) p->invalidate();

    this->IAmEvaluated = false;
  }

  //! \brief Checks if this operation will return any data.
  virtual bool hasData() const
  {
    for (FFaOperation<RetType>* p : myParams)
      if (p && p->hasData())
        return true;

    return false;
  }

  //! \brief Invokes the actual operation.
  virtual bool evaluate(RetType& value)
  {
    bool ok = !this->myOperationCB.empty();
    if (ok)
    {
      std::vector<RetType> pVals(myParams.size());
      for (size_t i = 0; i < myParams.size(); i++)
        if (myParams[i])
          ok &= myParams[i]->invoke(pVals[i]);

      this->myOperationCB.invoke(value,pVals);
    }

    return ok;
  }

protected:
  //! \brief The destructor deletes the child operations, if any.
  virtual ~FFaNToOneOp()
  {
    for (FFaOperation<RetType>* p : myParams)
      if (p) p->unref();
  }

private:
  std::vector<FFaOperation<RetType>*> myParams; //!< Parameter child operations

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


namespace FFaOper
{
  //! \brief Defines a special result value and its associated equivalent value.
  //! \details The special result value is used by some unary operations
  //! transforming result vectors into a double value. When \a value is
  //! encountered among the input values for an operation, it will be
  //! replaced by \a equiv before the actual operation is invoked.
  void setSpecialResultValue(double value, double equiv);

#ifdef FT_USE_MEMPOOL
  //! \brief Convenience alias for the operations memory pool manager type.
  using MemPoolMgr = FFaSingelton<FFaMemPoolMgr,FFaOperationBase>;
  //! \brief Returns a pointer to the operations memory pool manager.
  FFaMemPoolMgr* getMemPoolMgr();
  //! \brief Frees the memory pool associated with the operation manager.
  void freeMemPools(bool removeMgr = false);
#endif
}

//! \cond DO_NOT_DOCUMENT
// Constructor specializations for basic types

template <> inline
FFaOperation<int>::FFaOperation() : myCachedValue(0), IAmEvaluated(false) {}

template <> inline
FFaOperation<float>::FFaOperation() : myCachedValue(0), IAmEvaluated(false) {}

template <> inline
FFaOperation<double>::FFaOperation() : myCachedValue(0), IAmEvaluated(false) {}

// Initialization of static variables

template <class CB>
std::map<std::string,CB>* FFaOperationContainer<CB>::ourOperations = NULL;

template <class CB>
std::string FFaOperationContainer<CB>::ourDefaultOper;

#ifdef FT_USE_MEMPOOL

template <class RetType, class PrmType>
FFaMemPool FFaUnaryOp<RetType,PrmType>::ourMemPool =
  FFaMemPool(sizeof(FFaUnaryOp<RetType,PrmType>),FFaOper::getMemPoolMgr());

template <class RetType>
FFaMemPool FFaNToOneOp<RetType>::ourMemPool =
  FFaMemPool(sizeof(FFaNToOneOp<RetType>),FFaOper::getMemPoolMgr());

#endif
//! \endcond
#endif
