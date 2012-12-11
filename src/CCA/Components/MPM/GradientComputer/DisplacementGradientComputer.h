#ifndef __DISPLACEMENT_GRADIENT_COMPUTER_H__
#define __DISPLACEMENT_GRADIENT_COMPUTER_H__

#include <CCA/Components/MPM/GradientComputer/GradientComputer.h>


namespace Uintah {

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class DisplacementGradientComputer
    \brief Class for computing displacement gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class DisplacementGradientComputer : public GradientComputer {

  public:
         
    DisplacementGradientComputer(MPMFlags* MFlag);
    DisplacementGradientComputer(const DisplacementGradientComputer* gc);
    virtual ~DisplacementGradientComputer();

    // Make a clone of the gradient computer
    DisplacementGradientComputer* clone();

    // Actually compute displacement gradient
    void computeDispGrad(const ParticleInterpolator* interp,
                         const double* oodx,
                         constParticleVariable<Point> px,
                         constParticleVariable<Matrix3> psize,
                         constParticleVariable<Matrix3> pDefGrad_old,
                         Matrix3& dispGrad_new);

  };
} // End namespace Uintah
      


#endif  // __DISPLACEMENT_GRADIENT_COMPUTER_H__

