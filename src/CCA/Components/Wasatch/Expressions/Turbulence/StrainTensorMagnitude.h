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

#ifndef StrainTensorMagnitude_Expr_h
#define StrainTensorMagnitude_Expr_h

#include <CCA/Components/Wasatch/Operators/OperatorTypes.h>

#include <expression/Expression.h>

/**
 *  \brief obtain the tag for the strain tensor magnitude
 */
Expr::Tag straintensormagnitude_tag();
Expr::Tag wale_tensormagnitude_tag();
Expr::Tag vreman_tensormagnitude_tag();

/**
 *  \class StrainTensorMagnitude
 *  \authors Amir Biglari, Tony Saad
 *  \date   Jan, 2012. (Originally created: June, 2012).
 *  \ingroup Expressions
 *
 *  \brief given all components of the velocity, \f$u_i\f$, this calculates strain tensor magnitude, \f$|\tilde{S}|=(2\tilde{S_{kl}}\tilde{S_{kl}})^{1/2}\f$ where \f$S_{kl}=\frac{1}{2}(\frac{\partial\tilde{u}_k}{\partial x_l}+\frac{\partial\tilde{u}_l}{\partial x_k})\f$.
 *
 */
class StrainTensorMagnitude
 : public Expr::Expression<SVolField>
{
protected:
  const Expr::Tag vel1t_, vel2t_, vel3t_;
  
  //A_SURF_B_Field = A vol, B surface
  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, XVolField, SpatialOps::structured::XSurfYField >::type dudyT;
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, YVolField, SpatialOps::structured::YSurfXField >::type dvdxT; 
  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, XVolField, SpatialOps::structured::XSurfZField >::type dudzT; 
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, ZVolField, SpatialOps::structured::ZSurfXField >::type dwdxT; 
  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, YVolField, SpatialOps::structured::YSurfZField >::type dvdzT; 
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, ZVolField, SpatialOps::structured::ZSurfYField >::type dwdyT; 
  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, XVolField, SVolField >::type dudxT;
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, YVolField, SVolField >::type dvdyT;
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Gradient, ZVolField, SVolField >::type dwdzT;

  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, SpatialOps::structured::XSurfYField, SVolField >::type XYInterpT;  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, SpatialOps::structured::YSurfXField, SVolField >::type YXInterpT;

  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, SpatialOps::structured::XSurfZField, SVolField >::type XZInterpT;  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, SpatialOps::structured::ZSurfXField, SVolField >::type ZXInterpT;
  
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, SpatialOps::structured::YSurfZField, SVolField >::type YZInterpT;
  typedef SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, SpatialOps::structured::ZSurfYField, SVolField >::type ZYInterpT;

  const XVolField* vel1_;
  const YVolField* vel2_;
  const ZVolField* vel3_;
  
  const bool doX_, doY_, doZ_;
  
  const dudyT* dudyOp_;
  const dvdxT* dvdxOp_;
  const dudzT* dudzOp_;
  const dwdxT* dwdxOp_;
  const dvdzT* dvdzOp_;
  const dwdyT* dwdyOp_;
  const dudxT* dudxOp_;
  const dvdyT* dvdyOp_;
  const dwdzT* dwdzOp_;

  const XYInterpT* xyInterpOp_;
  const YXInterpT* yxInterpOp_;  

  const XZInterpT* xzInterpOp_;
  const ZXInterpT* zxInterpOp_;
  
  const YZInterpT* yzInterpOp_;
  const ZYInterpT* zyInterpOp_;
  
  StrainTensorMagnitude( const Expr::Tag& vel1tag,
                         const Expr::Tag& vel2tag,
                         const Expr::Tag& vel3tag );
 
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:

    /**
     *  \param vel1tag the first component of the velocity 
     *  \param vel2tag the second component of the velocity 
     *  \param vel3tag the third component of the velocity 
     */
    Builder( const Expr::Tag& result,
             const Expr::Tag& vel1tag,
             const Expr::Tag& vel2tag,
             const Expr::Tag& vel3tag );
    ~Builder(){}
    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag v1t_, v2t_, v3t_;
  };

  ~StrainTensorMagnitude();

  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );
  void evaluate();
};

/**
 *  \class  WaleTensorMagnitude
 *  \author Tony Saad
 *  \date   June, 2012
 *  \ingroup Expressions
 *
 *  \brief This calculates the square velocity gradient tensor. 
           This is used in the W.A.L.E. turbulent model. 
           See:
           Nicoud and Ducros, 1999, Subgrid-Scale Stress Modelling Based on the
           Square of the Velocity Gradient Tensor
 *
 */
class WaleTensorMagnitude : public StrainTensorMagnitude {

  WaleTensorMagnitude( const Expr::Tag& vel1tag,
                               const Expr::Tag& vel2tag,
                               const Expr::Tag& vel3tag);  
  public:
    class Builder : public Expr::ExpressionBuilder
    {
    public:
      
      /**
       *  \param vel1tag the first component of the velocity 
       *  \param vel2tag the second component of the velocity 
       *  \param vel3tag the third component of the velocity 
       */
      Builder( const Expr::Tag& result,
               const Expr::Tag& vel1tag,
               const Expr::Tag& vel2tag,
               const Expr::Tag& vel3tag );
      ~Builder(){}
      Expr::ExpressionBase* build() const;
      
    private:
      const Expr::Tag v1t_, v2t_, v3t_;
    };
  
  ~WaleTensorMagnitude();
  void evaluate();
};

/**
 *  \class  VremanTensorMagnitude
 *  \author Tony Saad
 *  \date   June, 2012
 *  \ingroup Expressions
 *
 *  \brief This calculates the square velocity gradient tensor. 
 This is used in the W.A.L.E. turbulent model. 
 See:
 Nicoud and Ducros, 1999, Subgrid-Scale Stress Modelling Based on the
 Square of the Velocity Gradient Tensor
 *
 */
class VremanTensorMagnitude : public StrainTensorMagnitude {
  
  VremanTensorMagnitude( const Expr::Tag& vel1tag,
                              const Expr::Tag& vel2tag,
                              const Expr::Tag& vel3tag);  
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    
    /**
     *  \param vel1tag the first component of the velocity 
     *  \param vel2tag the second component of the velocity 
     *  \param vel3tag the third component of the velocity 
     */
    Builder( const Expr::Tag& result,
            const Expr::Tag& vel1tag,
            const Expr::Tag& vel2tag,
            const Expr::Tag& vel3tag );
    ~Builder(){}
    Expr::ExpressionBase* build() const;
    
  private:
    const Expr::Tag v1t_, v2t_, v3t_;
  };
  
  ~VremanTensorMagnitude();
  void evaluate();
};
#endif // StrainTensorMagnitude_Expr_h
