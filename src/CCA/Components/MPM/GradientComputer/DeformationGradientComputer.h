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

  private:

    const Matrix3 Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

  public:
         
    DeformationGradientComputer(MPMFlags* MFlag);
    DeformationGradientComputer(const DeformationGradientComputer* gc);
    virtual ~DeformationGradientComputer();

    // Make a clone of the gradient computer
    DeformationGradientComputer* clone();

    // Computes and requires     
    void addComputesAndRequires(Task* task,
                                const MPMMaterial* mpm_matl,
                                const PatchSet*);

    void addComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches,
                                const bool /*recurse*/,
                                const bool SchedParent) const;

    void computeDeformationGradientExplicit(const Patch* patch,
                                            const MPMMaterial* mpm_matl,
                                            DataWarehouse& old_dw,
                                            DataWarehouse& new_dw);

    void computeDeformationGradientImplicit(const Patch* patch,
                                            const MPMMaterial* mpm_matl,
                                            DataWarehouse& old_dw,
                                            DataWarehouse& parent_old_dw,
                                            DataWarehouse& new_dw);

  protected:

    void addComputesAndRequiresExplicit(Task* task,
                                        const MPMMaterial* mpm_matl);
   
    void addComputesAndRequiresImplicit(Task* task,
                                        const MPMMaterial* mpm_matl);

    void computeDeformationGradientFromVelocity(const Matrix3& velGrad_old,
                                                const Matrix3& velGrad_new,
                                                const Matrix3& defGrad_old,
                                                const double& delT,
                                                Matrix3& defGrad_new,
                                                Matrix3& defGrad_inc);

    void computeDeformationGradientFromTotalDisplacement(const Matrix3& dispGrad_new,
                                                         const Matrix3& defGrad_old,
                                                         Matrix3& defGrad_new,
                                                         Matrix3& defGrad_inc);

    void seriesUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                     const Matrix3& defGrad_old,
                                     const double& delT,
                                     Matrix3& defGrad_new,
                                     Matrix3& defGrad_inc);

    void seriesUpdateLinearVelGrad(const Matrix3& velGrad_old,
                                   const Matrix3& velGrad_new,
                                   const Matrix3& defGrad_old,
                                   const double& delT,
                                   Matrix3& defGrad_new,
                                   Matrix3& defGrad_inc);

    void subcycleUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                       const Matrix3& defGrad_old,
                                       const double& delT,
                                       Matrix3& defGrad_new,
                                       Matrix3& defGrad_inc);

    void computeDeformationGradientFromIncrementalDisplacement(const Matrix3& dispGrad_new,
                                                               const Matrix3& defGrad_old,
                                                               Matrix3& defGrad_new,
                                                               Matrix3& defGrad_inc);

  protected:

    MPMLabel* lb;
    MPMFlags* flag;
    int NGP;
    int NGN;
  };

} // End namespace Uintah
      


#endif  // __DEFORMATION_GRADIENT_COMPUTER_H__

