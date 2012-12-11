#ifndef __DEFORMATION_GRADIENT_COMPUTER_H__
#define __DEFORMATION_GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/ComputeSet.h>
#include <vector>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Containers/StaticArray.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/LinearInterpolator.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Math/FastMatrix.h>
#include <CCA/Components/MPM/MPMFlags.h>

namespace Uintah {

  class Task;
  class Patch;
  class VarLabel;
  class MPMLabel;
  class MPMFlags;
  class MPMMaterial;
  class DataWarehouse;
  class ParticleSubset;
  class ParticleVariableBase;

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class DeformationGradientComputer
    \brief Class for computing deformation gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class DeformationGradientComputer {

  public:
         
    DeformationGradientComputer(MPMFlags* MFlag);
    DeformationGradientComputer(const DeformationGradientComputer* gc);
    virtual ~DeformationGradientComputer();

    // Make a clone of the gradient computer
    DeformationGradientComputer* clone();

    void computeDeformationGradientFromDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleVariable<Matrix3> &Fnew,
                                           constParticleVariable<Matrix3> &Fold,
                                           Vector dx,
                                           ParticleInterpolator* interp);

    void computeDeformationGradientFromVelocity(
                                           constNCVariable<Vector> gVel,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> psize,
                                           constParticleVariable<Matrix3> Fold,
                                           ParticleVariable<Matrix3> &Fnew,
                                           Vector dx,
                                           ParticleInterpolator* interp,
                                           const double& delT);

    void computeDeformationGradientFromTotalDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           ParticleVariable<Matrix3> &Fnew,
                                           constParticleVariable<Matrix3> &Fold,
                                           Vector dx,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleInterpolator* interp);
                                                                                
    void computeDeformationGradientFromIncrementalDisplacement(
                                           constNCVariable<Vector> IncDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> Fold,
                                           ParticleVariable<Matrix3> &Fnew,
                                           Vector dx,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleInterpolator* interp);
  protected:

    MPMLabel* lb;
    MPMFlags* flag;
    int NGP;
    int NGN;
  };

} // End namespace Uintah
      


#endif  // __DEFORMATION_GRADIENT_COMPUTER_H__

