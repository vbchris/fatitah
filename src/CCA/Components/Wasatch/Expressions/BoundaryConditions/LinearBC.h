/*
 * The MIT License
 *
 * Copyright (c) 2012 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef LinearBC_Expr_h
#define LinearBC_Expr_h

#include <expression/Expression.h>

template< typename FieldT >
class LinearBC
: public BoundaryConditionBase<FieldT>
{
  LinearBC( const Expr::Tag& indepVarTag,
              const double a,
              const double b ) :
  indepVarTag_ (indepVarTag),
  a_(a), b_(b)
  {}
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    Builder( const Expr::Tag& resultTag,
            const Expr::Tag& indepVarTag,
            const double a,
            const double b) :
    ExpressionBuilder(resultTag),
    indepVarTag_ (indepVarTag),
    a_(a), b_(b)
    {}
    Expr::ExpressionBase* build() const{ return new LinearBC(indepVarTag_, a_, b_); }
  private:
    const Expr::Tag indepVarTag_;
    const double a_, b_;
  };
  
  ~LinearBC(){}
  void advertise_dependents( Expr::ExprDeps& exprDeps ){  exprDeps.requires_expression( indepVarTag_ );}
  void bind_fields( const Expr::FieldManagerList& fml ){
    const typename Expr::FieldMgrSelector<FieldT>::type& phifm = fml.template field_manager<FieldT>();
    x_    = &phifm.field_ref( indepVarTag_    );
  }
  void evaluate();
private:
  const FieldT* x_;
  const Expr::Tag indepVarTag_;
  const double a_, b_;
};

// ###################################################################
//
//                          Implementation
//
// ###################################################################


//--------------------------------------------------------------------

template< typename FieldT >
void
LinearBC<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& f = this->value();
  const double ci = this->ci_;
  const double cg = this->cg_;
  std::vector<int>::const_iterator ia = this->flatGhostPoints_.begin(); // ia is the ghost flat index
  std::vector<int>::const_iterator ib = this->flatInteriorPoints_.begin(); // ib is the interior flat index
  for( ; ia != this->flatGhostPoints_.end(); ++ia, ++ib )
    f[*ia] = ( ( a_ * (*x_)[*ia] + b_) - ci*f[*ib] ) / cg;
}

#endif // LinearBC_Expr_h
