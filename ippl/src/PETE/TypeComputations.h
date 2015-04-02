// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/


///////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//    TypeComputations.h
//
// CREATED
//    July 11, 1997
//
// DESCRIPTION
//    PETE: Portable Expression Template Engine.
//
//    This header file defines templates used to construct the
//    return types for unary, binary, and trinary operations.
//
///////////////////////////////////////////////////////////////////////////

#ifndef TYPE_COMPUTATIONS_H
#define TYPE_COMPUTATIONS_H


///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    PETE_Type2Index<T>
//
// DESCRIPTION
//    This template describes a set of trait classes that associate an
//    index with a type. Each concrete class---a type that is designed
//    to be used like a built-in type---must have a specialization
//    of this class that provides a unique index. This index is used
//    to figure out the default return type of a binary/trinary operation.
//    Specifically, the largest index is chosen.
//
///////////////////////////////////////////////////////////////////////////

template<class Type>
struct PETE_Type2Index
{
  enum { val = 666666 };
};

template<> struct PETE_Type2Index<bool>
{
  enum { val = 1 };
};

template<> struct PETE_Type2Index<char>
{
  enum { val = 2 };
};

template<> struct PETE_Type2Index<short>
{
  enum { val = 3 };
};

template<> struct PETE_Type2Index<int>
{
  enum { val = 4 };
};

template<> struct PETE_Type2Index<long>
{
  enum { val = 5 };
};

template<> struct PETE_Type2Index<float>
{
  enum { val = 6 };
};

template<> struct PETE_Type2Index<double>
{
  enum { val = 7 };
};


///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    PETE_Index2Type<Index>
//
// DESCRIPTION
//    This template describes a set of trait classes that associate an
//    a type with an index. Each concrete class---a type that is designed
//    to be used like a built-in type---must have a specialization
//    of this class that provides a type given its unique index.
//
// NOTE
//    It is good form to provide these, but they are really only necessary
//    if your compiler (and all compilers you will ever want to port
//    your code to) support partial specialization.
//
///////////////////////////////////////////////////////////////////////////

#if !defined(IPPL_USE_PARTIAL_SPECIALIZATION)

template<int Index>
struct PETE_Index2Type
{
};

template<> struct PETE_Index2Type<1>
{
  typedef bool type;
};

template<> struct PETE_Index2Type<2>
{
  typedef char type;
};

template<> struct PETE_Index2Type<3>
{
  typedef short type;
};

template<> struct PETE_Index2Type<4>
{
  typedef int type;
};

template<> struct PETE_Index2Type<5>
{
  typedef long type;
};

template<> struct PETE_Index2Type<6>
{
  typedef float type;
};

template<> struct PETE_Index2Type<7>
{
  typedef double type;
};

#endif


///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    PETEUnaryReturn<T, Op>
//
// DESCRIPTION
//    This template describes the default mechanism for calculating the
//    return type to a unary expression given the argument type T and
//    the operation type Op.
//
//    There are two sensible things one can do automatically:
//      o (Op::tag == PETE_UnaryPassThruTag) make the return type 
//        the same as the argument to the function/operation. 
//        For example, operator-(T) should return a T.
//      o (Op::tag != PETE_UnaryPassThruTag) return a type based entirely on 
//        the operation. For example, operator! always returns a bool.
//    To figure out which approach to take, PETEUnaryReturn uses the
//    tag from the operator and another template, PETE_ComputeUnaryType.
//    Consider negation (unary minus). The operator would be formed like:
//      struct OpUnaryMinus {
//        enum { tag = PETE_UnaryPassThruTag };
//      };
//    Logical negation (not) would be formed like:
//      struct OpNot {
//        enum { tag = PETE_Type2Index<bool> };
//        typedef bool type;
//      };
//    The minor redundancy (specifying the tag and the type) is required to 
//    allow easy support for compilers that may or may not support partial
//    specialization.
//
//    Special cases can be handled by directly specializing PETEUnaryReturn.
//    For example, the abs function typically returns a double when the
//    argument is a complex<double>. The appropriate specialization here
//    would be:
//      template<> struct PETEUnaryReturn< complex<double>, FnAbs > {
//        typedef double type;
//      };
//
///////////////////////////////////////////////////////////////////////////

const int PETE_UnaryPassThruTag = 0;

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

template<class T, class Op, int OpTag>
struct PETE_ComputeUnaryType
{
  typedef typename Op::type type;
};

template<class T, class Op>
struct PETE_ComputeUnaryType<T, Op, PETE_UnaryPassThruTag>
{
  typedef T type;
};

template<class T, class Op>
struct PETEUnaryReturn
{
  typedef typename PETE_ComputeUnaryType<T, Op, Op::tag>::type type;
};

#else

template<int t, int op>
struct PETE_ComputeUnaryType
{
  typedef typename 
    PETE_Index2Type<(op == PETE_UnaryPassThruTag ? t : op)>::type type;
};

template<class T, class Op>
struct PETEUnaryReturn
{
  typedef typename
    PETE_ComputeUnaryType<PETE_Type2Index<T>::val, Op::tag>::type type;
};

#endif


///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    PETEBinaryReturn<T1, T2, Op>
//
// DESCRIPTION
//    This template describes the default mechanism for calculating the
//    return type to a binary expression given the argument types T1 and
//    T2 and the operation type Op.
//
//    There are four sensible things one can do automatically:
//      o (Op::tag == PETE_BinaryPromoteTag) make the return type by 
//        promoting/converting the "simpler" type into the more "complex." 
//        For example, we typically want to do this with addition. 
//      o (Op::tag == PETE_BinaryUseLeftTag) return the type of the 
//        left-hand operand. For example, this is what happens with operator<<.
//      o (Op::tag == PETE_BinaryUseRightTag) return the type of the 
//        right-hand operand. 
//      o (otherwise) return a type based entirely on the operation. 
//        For example, operator!= always returns a bool.
//    To figure out which approach to take, PETEBinaryReturn uses the
//    tag from the operator and another template, PETE_ComputeBinaryType.
//    Consider addition. The operator would be formed like:
//      struct OpAdd {
//        enum { tag = PETE_BinaryPromoteTag };
//      };
//
//    Special cases can be handled by directly specializing PETEBinaryReturn.
//    For example, the multipplication between a matrix and a vector might do a
//    matrix/vector product, thereby returning a vector. The appropriate 
//    specialization here would be:
//      struct PETEBinaryReturn< Mat<double,3>, Vec<float,3>, OpMultipply > {
//        typedef Vector<double,3> type;
//      };
//    Notice how the element type is promoted.
//
///////////////////////////////////////////////////////////////////////////

const int PETE_BinaryPromoteTag = -2;
const int PETE_BinaryUseLeftTag = -1;
const int PETE_BinaryUseRightTag = 0;

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

// This is still harder than it has to be. There are bugs in
// the EDG front end.

template<class T1, class T2, class Op, int op>
struct PETE_ComputeBinaryType
{
  typedef typename Op::type type;
};

template<class T1, class T2, class Op>
struct PETE_ComputeBinaryType<T1, T2, Op, PETE_BinaryUseLeftTag>
{
  typedef T1 type;
};

template<class T1, class T2, class Op>
struct PETE_ComputeBinaryType<T1, T2, Op, PETE_BinaryUseRightTag>
{
  typedef T2 type;
};

template<class T1, class T2, bool lr>
struct PETE_ComputePromoteType
{
};

template<class T1, class T2>
struct PETE_ComputePromoteType<T1, T2, true>
{
  typedef T1 type;
};

template<class T1, class T2>
struct PETE_ComputePromoteType<T1, T2, false>
{
  typedef T2 type;
};

template<class T1, class T2, int t1, int t2>
struct PETE_ComputePromoteType2
{
  typedef typename
    PETE_ComputePromoteType<T1, T2, (t1 >= t2)>::type type;
};

template<class T1, class T2, class Op>
struct PETE_ComputeBinaryType<T1, T2, Op, PETE_BinaryPromoteTag>
{
  typedef typename PETE_ComputePromoteType2<T1, T2, 
    PETE_Type2Index<T1>::val, PETE_Type2Index<T2>::val>::type type;
};

template<class T1, class T2, class Op>
struct PETEBinaryReturn
{
  typedef typename PETE_ComputeBinaryType<T1, T2, Op, Op::tag>::type type;
};

#else

template<int t1, int t2, int op>
struct PETE_ComputeBinaryType
{
  typedef typename PETE_Index2Type
    <(op == PETE_BinaryPromoteTag ? 
     (t1 >= t2 ? t1 : t2) : 
     (op == PETE_BinaryUseLeftTag ? t1 : 
     (op == PETE_BinaryUseRightTag ? t2 : op)))>::type type;
};

template<class T1, class T2, class Op>
struct PETEBinaryReturn
{
  typedef typename PETE_ComputeBinaryType
    <PETE_Type2Index<T1>::val, PETE_Type2Index<T2>::val, Op::tag>::type type;
};

#endif


///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    PETETrinaryReturn<T1, T2, T3, Op>
//
// DESCRIPTION
//    This template describes the default mechanism for calculating the
//    return type to a trinary expression given the argument types T1, T2, and
//    T3 and the operation type Op. The only trinary expression supported
//    in C++ is the ?: operation. In this case, T1 should end up being bool
//    and the result of the calculation is of type Binary_Promotion(T1,T2)
//    with the value being that associated with T2 if T1's associated value 
//    turns out to be true and T3 if T1's associated value turns out to be 
//    false.
//
///////////////////////////////////////////////////////////////////////////

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

template<class T1, class T2, class T3, class Op>
struct PETETrinaryReturn
{
  typedef typename PETE_ComputeBinaryType<T2, T3, Op, Op::tag>::type type;
};

#else

template<class T1, class T2, class T3, class Op>
struct PETETrinaryReturn
{
  typedef typename PETE_ComputeBinaryType
  <PETE_Type2Index<T2>::val, PETE_Type2Index<T3>::val, Op::tag>::type type;
};

#endif


///////////////////////////////////////////////////////////////////////////
//
// UNARY OPERATORS: -, +, ~, !, Identity
//
///////////////////////////////////////////////////////////////////////////

struct OpIdentity
{
#ifdef IPPL_PURIFY
  OpIdentity() {}
  OpIdentity(const OpIdentity &) {}
  OpIdentity& operator=(const OpIdentity &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct OpUnaryMinus
{
#ifdef IPPL_PURIFY
  OpUnaryMinus() {}
  OpUnaryMinus(const OpUnaryMinus &) {}
  OpUnaryMinus& operator=(const OpUnaryMinus &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct OpUnaryPlus
{
#ifdef IPPL_PURIFY
  OpUnaryPlus() {}
  OpUnaryPlus(const OpUnaryPlus &) {}
  OpUnaryPlus& operator=(const OpUnaryPlus &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct OpBitwiseNot
{
#ifdef IPPL_PURIFY
  OpBitwiseNot() {}
  OpBitwiseNot(const OpBitwiseNot &) {}
  OpBitwiseNot& operator=(const OpBitwiseNot &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct OpNot
{
#ifdef IPPL_PURIFY
  OpNot() {}
  OpNot(const OpNot &) {}
  OpNot& operator=(const OpNot &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

template <class T>
struct OpCast
{
#ifdef IPPL_PURIFY
  OpCast() {}
  OpCast(const OpCast<T> &) {}
  OpCast& operator=(const OpCast<T> &) { return *this; }
#endif
  typedef T type;
  enum { tag = PETE_Type2Index<T>::val };
};

///////////////////////////////////////////////////////////////////////////
//
// UNARY FUNCTIONS: acos, asin, atan, ceil, cos, cosh, exp, fabs, floor,
//                  log, log10, sin, sinh, sqrt, tan, tanh
//
///////////////////////////////////////////////////////////////////////////

struct FnArcCos
{
#ifdef IPPL_PURIFY
  FnArcCos() {}
  FnArcCos(const FnArcCos &) {}
  FnArcCos& operator=(const FnArcCos &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnArcSin
{
#ifdef IPPL_PURIFY
  FnArcSin() {}
  FnArcSin(const FnArcSin &) {}
  FnArcSin& operator=(const FnArcSin &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnArcTan
{
#ifdef IPPL_PURIFY
  FnArcTan() {}
  FnArcTan(const FnArcTan &) {}
  FnArcTan& operator=(const FnArcTan &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnCeil
{
#ifdef IPPL_PURIFY
  FnCeil() {}
  FnCeil(const FnCeil &) {}
  FnCeil& operator=(const FnCeil &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnCos
{
#ifdef IPPL_PURIFY
  FnCos() {}
  FnCos(const FnCos &) {}
  FnCos& operator=(const FnCos &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnHypCos
{
#ifdef IPPL_PURIFY
  FnHypCos() {}
  FnHypCos(const FnHypCos &) {}
  FnHypCos& operator=(const FnHypCos &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnExp
{
#ifdef IPPL_PURIFY
  FnExp() {}
  FnExp(const FnExp &) {}
  FnExp& operator=(const FnExp &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnFabs
{
#ifdef IPPL_PURIFY
  FnFabs() {}
  FnFabs(const FnFabs &) {}
  FnFabs& operator=(const FnFabs &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnFloor
{
#ifdef IPPL_PURIFY
  FnFloor() {}
  FnFloor(const FnFloor &) {}
  FnFloor& operator=(const FnFloor &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnLog
{
#ifdef IPPL_PURIFY
  FnLog() {}
  FnLog(const FnLog &) {}
  FnLog& operator=(const FnLog &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnLog10
{
#ifdef IPPL_PURIFY
  FnLog10() {}
  FnLog10(const FnLog10 &) {}
  FnLog10& operator=(const FnLog10 &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnSin
{
#ifdef IPPL_PURIFY
  FnSin() {}
  FnSin(const FnSin &) {}
  FnSin& operator=(const FnSin &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnHypSin
{
#ifdef IPPL_PURIFY
  FnHypSin() {}
  FnHypSin(const FnHypSin &) {}
  FnHypSin& operator=(const FnHypSin &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnSqrt
{
#ifdef IPPL_PURIFY
  FnSqrt() {}
  FnSqrt(const FnSqrt &) {}
  FnSqrt& operator=(const FnSqrt &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnTan
{
#ifdef IPPL_PURIFY
  FnTan() {}
  FnTan(const FnTan &) {}
  FnTan& operator=(const FnTan &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnHypTan
{
#ifdef IPPL_PURIFY
  FnHypTan() {}
  FnHypTan(const FnHypTan &) {}
  FnHypTan& operator=(const FnHypTan &) { return *this; }
#endif
  enum { tag = PETE_UnaryPassThruTag };
};


///////////////////////////////////////////////////////////////////////////
//
// BINARY OPERATORS: +, -, *, /, %, <, >, <=, >=, ==, !=, &&, ||, ^, &, 
//                   |, <<, >>
//
///////////////////////////////////////////////////////////////////////////

struct OpAdd
{
#ifdef IPPL_PURIFY
  OpAdd() {}
  OpAdd(const OpAdd &) {}
  OpAdd& operator=(const OpAdd &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpSubtract
{
#ifdef IPPL_PURIFY
  OpSubtract() {}
  OpSubtract(const OpSubtract &) {}
  OpSubtract& operator=(const OpSubtract &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpMultipply
{
#ifdef IPPL_PURIFY
  OpMultipply() {}
  OpMultipply(const OpMultipply &) {}
  OpMultipply& operator=(const OpMultipply &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpDivide
{
#ifdef IPPL_PURIFY
  OpDivide() {}
  OpDivide(const OpDivide &) {}
  OpDivide& operator=(const OpDivide &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpMod
{
#ifdef IPPL_PURIFY
  OpMod() {}
  OpMod(const OpMod &) {}
  OpMod& operator=(const OpMod &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpLT
{
#ifdef IPPL_PURIFY
  OpLT() {}
  OpLT(const OpLT &) {}
  OpLT& operator=(const OpLT &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpGT
{
#ifdef IPPL_PURIFY
  OpGT() {}
  OpGT(const OpGT &) {}
  OpGT& operator=(const OpGT &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpLE
{
#ifdef IPPL_PURIFY
  OpLE() {}
  OpLE(const OpLE &) {}
  OpLE& operator=(const OpLE &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpGE
{
#ifdef IPPL_PURIFY
  OpGE() {}
  OpGE(const OpGE &) {}
  OpGE& operator=(const OpGE &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpEQ
{
#ifdef IPPL_PURIFY
  OpEQ() {}
  OpEQ(const OpEQ &) {}
  OpEQ& operator=(const OpEQ &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpNE
{
#ifdef IPPL_PURIFY
  OpNE() {}
  OpNE(const OpNE &) {}
  OpNE& operator=(const OpNE &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpAnd
{
#ifdef IPPL_PURIFY
  OpAnd() {}
  OpAnd(const OpAnd &) {}
  OpAnd& operator=(const OpAnd &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpOr
{
#ifdef IPPL_PURIFY
  OpOr() {}
  OpOr(const OpOr &) {}
  OpOr& operator=(const OpOr &) { return *this; }
#endif
  typedef bool type;
  enum { tag = PETE_Type2Index<bool>::val };
};

struct OpBitwiseXor
{
#ifdef IPPL_PURIFY
  OpBitwiseXor() {}
  OpBitwiseXor(const OpBitwiseXor &) {}
  OpBitwiseXor& operator=(const OpBitwiseXor &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpBitwiseAnd
{
#ifdef IPPL_PURIFY
  OpBitwiseAnd() {}
  OpBitwiseAnd(const OpBitwiseAnd &) {}
  OpBitwiseAnd& operator=(const OpBitwiseAnd &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpBitwiseOr
{
#ifdef IPPL_PURIFY
  OpBitwiseOr() {}
  OpBitwiseOr(const OpBitwiseOr &) {}
  OpBitwiseOr& operator=(const OpBitwiseOr &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct OpLeftShift
{
#ifdef IPPL_PURIFY
  OpLeftShift() {}
  OpLeftShift(const OpLeftShift &) {}
  OpLeftShift& operator=(const OpLeftShift &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpRightShift
{
#ifdef IPPL_PURIFY
  OpRightShift() {}
  OpRightShift(const OpRightShift &) {}
  OpRightShift& operator=(const OpRightShift &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};


///////////////////////////////////////////////////////////////////////////
//
// BINARY FUNCTIONS: copysign, ldexp, pow, fmod, atan2
//
///////////////////////////////////////////////////////////////////////////

struct FnCopysign
{
#ifdef IPPL_PURIFY
  FnCopysign() {}
  FnCopysign(const FnCopysign &) {}
  FnCopysign& operator=(const FnCopysign &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnLdexp
{
#ifdef IPPL_PURIFY
  FnLdexp() {}
  FnLdexp(const FnLdexp &) {}
  FnLdexp& operator=(const FnLdexp &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnPow
{
#ifdef IPPL_PURIFY
  FnPow() {}
  FnPow(const FnPow &) {}
  FnPow& operator=(const FnPow &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnFmod
{
#ifdef IPPL_PURIFY
  FnFmod() {}
  FnFmod(const FnFmod &) {}
  FnFmod& operator=(const FnFmod &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnArcTan2
{
#ifdef IPPL_PURIFY
  FnArcTan2() {}
  FnArcTan2(const FnArcTan2 &) {}
  FnArcTan2& operator=(const FnArcTan2 &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};


///////////////////////////////////////////////////////////////////////////
//
// ASSIGNMENT OPERATORS: =, +=, -=, *=, /=, %=, &=, ^=, |=, <<=, >>=
//
///////////////////////////////////////////////////////////////////////////

struct OpAssign
{
#ifdef IPPL_PURIFY
  OpAssign() {}
  OpAssign(const OpAssign &) {}
  OpAssign& operator=(const OpAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpAddAssign
{
#ifdef IPPL_PURIFY
  OpAddAssign() {}
  OpAddAssign(const OpAddAssign &) {}
  OpAddAssign& operator=(const OpAddAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpSubtractAssign
{
#ifdef IPPL_PURIFY
  OpSubtractAssign() {}
  OpSubtractAssign(const OpSubtractAssign &) {}
  OpSubtractAssign& operator=(const OpSubtractAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpMultipplyAssign
{
#ifdef IPPL_PURIFY
  OpMultipplyAssign() {}
  OpMultipplyAssign(const OpMultipplyAssign &) {}
  OpMultipplyAssign& operator=(const OpMultipplyAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpDivideAssign
{
#ifdef IPPL_PURIFY
  OpDivideAssign() {}
  OpDivideAssign(const OpDivideAssign &) {}
  OpDivideAssign& operator=(const OpDivideAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpModAssign
{
#ifdef IPPL_PURIFY
  OpModAssign() {}
  OpModAssign(const OpModAssign &) {}
  OpModAssign& operator=(const OpModAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpBitwiseXorAssign
{
#ifdef IPPL_PURIFY
  OpBitwiseXorAssign() {}
  OpBitwiseXorAssign(const OpBitwiseXorAssign &) {}
  OpBitwiseXorAssign& operator=(const OpBitwiseXorAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpBitwiseAndAssign
{
#ifdef IPPL_PURIFY
  OpBitwiseAndAssign() {}
  OpBitwiseAndAssign(const OpBitwiseAndAssign &) {}
  OpBitwiseAndAssign& operator=(const OpBitwiseAndAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpBitwiseOrAssign
{
#ifdef IPPL_PURIFY
  OpBitwiseOrAssign() {}
  OpBitwiseOrAssign(const OpBitwiseOrAssign &) {}
  OpBitwiseOrAssign& operator=(const OpBitwiseOrAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpLeftShiftAssign
{
#ifdef IPPL_PURIFY
  OpLeftShiftAssign() {}
  OpLeftShiftAssign(const OpLeftShiftAssign &) {}
  OpLeftShiftAssign& operator=(const OpLeftShiftAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpRightShiftAssign
{
#ifdef IPPL_PURIFY
  OpRightShiftAssign() {}
  OpRightShiftAssign(const OpRightShiftAssign &) {}
  OpRightShiftAssign& operator=(const OpRightShiftAssign &) { return *this; }
#endif
  enum { tag = PETE_BinaryUseLeftTag };
};


///////////////////////////////////////////////////////////////////////////
//
// TRINARY OPERATORS: ?: (where)
//
///////////////////////////////////////////////////////////////////////////

struct OpWhere
{
#ifdef IPPL_PURIFY
  OpWhere() {}
  OpWhere(const OpWhere &) {}
  OpWhere& operator=(const OpWhere &) { return *this; }
#endif
  enum { tag = PETE_BinaryPromoteTag };
};


///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif // TYPE_COMPUTATIONS_H

/***************************************************************************
 * $RCSfile: TypeComputations.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:28 $
 * IPPL_VERSION_ID: $Id: TypeComputations.h,v 1.1.1.1 2003/01/23 07:40:28 adelmann Exp $ 
 ***************************************************************************/
