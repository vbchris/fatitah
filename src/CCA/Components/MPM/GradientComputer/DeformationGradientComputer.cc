#include <CCA/Components/MPM/GradientComputer/DeformationGradientComputer.h>
#include <Core/Exceptions/InvalidValue.h>
#include <iostream>

using namespace Uintah;

DeformationGradientComputer::DeformationGradientComputer(MPMFlags* Mflag)
{
  lb = scinew MPMLabel();
  flag = Mflag;
  if(flag->d_8or27==8){
    NGN=1;
  } else{ 
    NGN=2;
  }
}

DeformationGradientComputer::DeformationGradientComputer(const DeformationGradientComputer* dg)
{
  lb = scinew MPMLabel();
  flag = cm->flag;
  NGN = cm->NGN;
  NGP = cm->NGP;
}

DeformationGradientComputer* DeformationGradientComputer::clone()
{
  return scinew DeformationGradientComputer(*this);
}

DeformationGradientComputer::~DeformationGradientComputer()
{
  delete lb;
}

void
DeformationGradientComputer::computeDeformationGradient(const Patch* patch,
                                                        const MPMMaterial* mpm_matl,
                                                        DataWarehouse& old_dw)
{
  // Get particle info and patch info
  int dwi = mpm_matl->getDWIndex();
  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
  Vector dx = patch->dCell();
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  ParticleInterpolator* interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data 
  // for vel grad and def grad calculation
  constParticleVariable<Short27> pgCode;
  constParticleVariable<double>  pMass;
  constParticleVariable<long64>  pParticleID;
  constParticleVariable<Point>   px;
  constParticleVariable<Matrix3> pDefGrad_old, pVelGrad_old, pDispGrad_old;
  constParticleVariable<Matrix3> pSize;

  // Set up variables to store new particle and grid data 
  // for vel grad and def grad calculation
  ParticleVariable<double>     pVolume_new;
  ParticleVariable<Matrix3>    pDefGrad_new, pVelGrad_new, pDispGrad_new;
  constNCVariable<Vector>      gDisp;
  constNCVariable<Vector>      gVelocity;
  constNCVariable<Vector>      GVelocity;

  // Get the old data
  new_dw->get(gVelocity, lb->gVelocityStarLabel, dwi, patch, gac, NGN);
  if (flag->d_fracture) {
    new_dw->get(pgCode,    lb->pgCodeLabel, pset);
    new_dw->get(GVelocity, lb->GVelocityStarLabel, dwi, patch, gac, NGN);
  }
  old_dw->get(px,                lb->pXLabel,                  pset);
  old_dw->get(pMass,             lb->pMassLabel,               pset);
  old_dw->get(pSize,             lb->pSizeLabel,               pset);
  old_dw->get(pDefGrad_old,      lb->pDefGradLabel,            pset);
  old_dw->get(pVelGrad_old,      lb->pVelGradLabel,            pset);
  old_dw->get(pDispGrad_old,     lb->pDispGradLabel,           pset);

  // Allocate new data
  new_dw->allocateAndPut(pVolume_new,   lb->pVolumeLabel_preReloc, pset);
  new_dw->allocateAndPut(pDefGrad_new,  lb->pDefGradLabel_preReloc,  pset);
  new_dw->allocateAndPut(pVelGrad_new,  lb->pVelGradLabel_preReloc,  pset);
  new_dw->allocateAndPut(pDispGrad_new, lb->pDispGradLabel_preReloc,  pset);
      
  // Create VelocityGradientComputer
  VelocityGradientComputer gradComp(flag);

  // Loop through particles
  double J = 1.0;
  for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
    particleIndex idx = *iter;

    // Initialize velocity gradient
    Matrix3 velGrad_new(0.0);
    Matrix3 defGradInc(0.0); // **WARNING** should be one and not zero
  
    // Compute velocity gradient
    gradComp.computeVelGrad(interpolator, oodx, pParticleID[idx], pgCode[idx], px[idx], pSize[idx], 
                            pDefGrad_old[idx], gVelocity, GVelocity, velGrad_new);

    // Update velocity gradient
    pVelGrad_new[idx] = velGrad_new;
    //if (isnan(velGrad_new.Norm())) {
    //  cerr << "particle = " << idx << " velGrad = " << velGrad_new << endl;
    //  throw InvalidValue("**ERROR**: Nan in velocity gradient value", __FILE__, __LINE__);
    //}

    // Improve upon first order estimate of deformation gradient
        bool d_taylorSeriesForDefGrad = true;
        int d_numTaylorTerms = 10;
        int num_scs = 1;
        if (d_taylorSeriesForDefGrad) {
          // Use Taylor series expansion
          // Compute mid point velocity gradient
          Matrix3 Amat = (pVelGrad_old[idx] + pVelGrad_new[idx])*(0.5*delT);
          defGradInc = Amat.Exponential(d_numTaylorTerms);
          pDefGrad_new[idx] = pDefGradInc*pDefGrad_old[idx];
        } else {
          Matrix3 F = pDefGrad_old[idx];
          double Lnorm_dt = velGrad_new.Norm()*delT;
          num_scs = max(1,2*((int) Lnorm_dt));
          if(num_scs > 1000){
            cout << "NUM_SCS = " << num_scs << endl;
          }
          double dtsc = delT/(double (num_scs));
          Matrix3 OP_tensorL_DT = Identity + velGrad_new*dtsc;
          for(int n=0;n<num_scs;n++){
            F = OP_tensorL_DT*F;
            // if(num_scs >1000){
            //   cerr << "n = " << n << endl;
            //   cerr << "F = " << F << endl;
            //   cerr << "J = " << F.Determinant() << endl << endl;
            // }
          }
          pDefGrad_new[idx] = F;
          defGradInc = pDefGrad_new[idx]*pDefGrad_old[idx].Inverse();
        }

        // Check 1: Look at Jacobian
        J = pDefGrad_new[idx].Determinant();
        if (!(J > 0.0)) {
          constParticleVariable<long64> pParticleID;
          old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
          cerr << "matl = "  << mpm_matl << " dwi = " << dwi << " particle = " << idx
               << " particleID = " << pParticleID[idx] << endl;
          cerr << "velGrad = " << velGrad_new << endl;
          cerr << "F_old = " << pDefGrad[idx]     << endl;
          cerr << "F_inc = " << pDefGradInc       << endl;
          cerr << "F_new = " << pDefGrad_new[idx] << endl;
          cerr << "J = "     << J                 << endl;
          cerr << "NUM_SCS = " << num_scs << endl;
          cerr << "**ERROR** Negative Jacobian of deformation gradient in material # ="
               << m << << " and particle " << pParticleID[idx]  << " which has mass "
               << pMass[idx] << endl;
          throw InvalidValue("**ERROR**:UCNH", __FILE__, __LINE__);
        }
      } // End of loop over particles

      // The following is used only for pressure stabilization
      CCVariable<double> J_CC;
      new_dw->allocateTemporary(J_CC,       patch);
      J_CC.initialize(0.);
      if(flag->d_doPressureStabilization) {
        CCVariable<double> vol_0_CC;
        CCVariable<double> vol_CC;
        new_dw->allocateTemporary(vol_0_CC, patch);
        new_dw->allocateTemporary(vol_CC,   patch);

        vol_0_CC.initialize(0.);
        vol_CC.initialize(0.);
  
        // First loop thru particles
        for(ParticleSubset::iterator iter = pset->begin();
            iter != pset->end(); iter++){
          particleIndex idx = *iter;

          // get the volumetric part of the deformation
          J = pDefGrad_new[idx].Determinant();

          // Get the deformed volume
          pVolume_new[idx]=(pMass[idx]/rho_orig)*J;

          IntVector cell_index;
          patch->findCell(px[idx],cell_index);

          vol_CC[cell_index]  +=pVolume_new[idx];
          vol_0_CC[cell_index]+=pMass[idx]/rho_orig;
        }

        // Compute cell centered J
        for(CellIterator iter=patch->getCellIterator(); !iter.done();iter++){
          IntVector c = *iter;
          J_CC[c]=vol_CC[c]/vol_0_CC[c];
        }

        // Second loop thru particles
        for(ParticleSubset::iterator iter = pset->begin();
            iter != pset->end(); iter++){
          particleIndex idx = *iter;
          IntVector cell_index;
          patch->findCell(px[idx],cell_index);

          // get the original volumetric part of the deformation
          J = pDefGrad_new[idx].Determinant();

          // Change F such that the determinant is equal to the average for
          // the cell
          pDefGrad_new[idx]*=cbrt(J_CC[cell_index]/J);
          defGradInc = pDefGrad_new[idx]*pDefGrad_old[idx].Inverse();

          // Update the deformed volume
          J = pdefGrad_new[idx].Determinant();
          pVolume_new[idx]= (pMass[idx]/rho_orig)*J;

          // Check 1: Look at Jacobian
          if (!(J > 0.0)) {
            cerr << "after pressure stab "          << endl;
            cerr << "matl = "  << mpm_matl          << endl;
            cerr << "F_old = " << pDefGrad[idx]     << endl;
            cerr << "F_inc = " << pDefGradInc       << endl;
            cerr << "F_new = " << pDefGrad_new[idx] << endl;
            cerr << "J = "     << J                 << endl;
            constParticleVariable<long64> pParticleID;
            old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
            cerr << "ParticleID = " << pParticleID[idx] << endl;
            cerr << "**ERROR** Negative Jacobian of deformation gradient"
                 << " in particle " << pParticleID[idx]  << " which has mass "
                 << pMass[idx] << endl;
            throw InvalidValue("**ERROR**:Negative Jacobian in UCNH",
                                __FILE__, __LINE__);
          }
        }

      } //end of pressureStabilization loop  at the patch level

}

// Use Taylor series expansion of exact solution
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::seriesUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                                         const Matrix3& defGrad_old,
                                                         const double& delT,
                                                         Matrix3& defGrad_new,
                                                         Matrix3& defGrad_inc)
{
  Matrix3 Amat = velGrad_new*delT;
  defGrad_inc = Amat.Exponential(d_numTaylorTerms);
  defGrad_new = defGrad_inc*defGrad_old;
  return;
}

// Use Taylor series expansion of exact solution
// Assume linear velocity gradient over timestep
void
DeformationGradientComputer::seriesUpdateLinearVelGrad(const Matrix3& velGrad_old,
                                                       const Matrix3& velGrad_new,
                                                       const Matrix3& defGrad_old,
                                                       const double& delT,
                                                       Matrix3& defGrad_new,
                                                       Matrix3& defGrad_inc)
{
  Matrix3 Amat = (velGrad_old + velGrad_new)*(0.5*delT);
  defGrad_inc = Amat.Exponential(d_numTaylorTerms);
  defGrad_new = defGrad_inc*defGrad_old;
  return;
}

// Use first term of series expansion of exact solution
// and subcycling.
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::subcycleUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                                           const Matrix3& defGrad_old,
                                                           const double& delT,
                                                           Matrix3& defGrad_new,
                                                           Matrix3& defGrad_inc)
{
  Matrix3 Identity; Identity.Identity();
  defGrad_new = defGrad_old;
  double Lnorm_dt = velGrad_new.Norm()*delT;
  int num_scs = max(1,2*((int) Lnorm_dt));
  if (num_scs > 1000) {
    cout << "NUM_SCS = " << num_scs << endl;
  }
  double dtsc = delT/(double (num_scs));
  Matrix3 OP_tensorL_DT = Identity + velGrad_new*dtsc;
  for(int n=0;n<num_scs;n++){
    defGrad_new = OP_tensorL_DT*defGrad_new;
    // if(num_scs >1000){
    //   cerr << "n = " << n << endl;
    //   cerr << "F = " << defGrad_new << endl;
    //   cerr << "J = " << defGrad_new.Determinant() << endl << endl;
    // }
  }
  defGradInc = defGrad_new*defGrad_old.Inverse();
  return;
}

void
DeformationGradientComputer::computeDeformationGradientFromDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleVariable<Matrix3> &Fnew,
                                           constParticleVariable<Matrix3> &Fold,
                                           Vector dx,
                                           ParticleInterpolator* interp) 
{
  Matrix3 dispGrad,Identity;
  Identity.Identity();
  vector<IntVector> ni(interp->size());
  vector<Vector> d_S(interp->size());
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
                                                                            
  for(ParticleSubset::iterator iter = pset->begin();
       iter != pset->end(); iter++){
    particleIndex idx = *iter;
                                                                            
    // Get the node indices that surround the cell
    interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);
                                                                            
    computeGrad(dispGrad, ni, d_S, oodx, gDisp);

    // Update the deformation gradient tensor to its time n+1 value.
    // Compute the deformation gradient from the displacement gradient
    Fnew[idx] = Identity + dispGrad;

    double J = Fnew[idx].Determinant();
    if (!(J > 0)) {
      ostringstream warn;
      warn << "**ERROR** : DeformationGradientComputer::computeDeformationGradientFromDisplacement" << endl << "Negative or zero determinant of Jacobian." << endl;
      warn << "     Particle = " << idx << " J = " << J << " position = " << px[idx] << endl;
      warn << "     Disp Grad = " << dispGrad << endl; 
      warn << "     F_new = " << Fnew[idx] << endl; 
      throw InvalidValue(warn.str(), __FILE__, __LINE__);
    }
  }
}

void 
DeformationGradientComputer::computeDeformationGradientFromVelocity(
                                           constNCVariable<Vector> gVel,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> psize,
                                           constParticleVariable<Matrix3> Fold,
                                           ParticleVariable<Matrix3> &Fnew,
                                           Vector dx,
                                           ParticleInterpolator* interp,
                                           const double& delT)
{
    Matrix3 velGrad,deformationGradientInc, Identity;
    Identity.Identity();
    vector<IntVector> ni(interp->size());
    vector<Vector> d_S(interp->size());
    double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    for(ParticleSubset::iterator iter = pset->begin();
      iter != pset->end(); iter++){
      particleIndex idx = *iter;

      // Get the node indices that surround the cell
      interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);

      computeGrad(velGrad, ni, d_S, oodx, gVel);

      // Compute the deformation gradient increment using the time_step
      // velocity gradient
      // F_n^np1 = dudx * dt + Identity
      deformationGradientInc = velGrad * delT + Identity;
                                                                              
      // Update the deformation gradient tensor to its time n+1 value.
      Fnew[idx] = deformationGradientInc * Fold[idx];

      double J = Fnew[idx].Determinant();
      if (!(J > 0)) {
        ostringstream warn;
        warn << "**ERROR** Negative or zero determinant of Jacobian."
             << " Particle has inverted." << endl;
        warn << "     Particle = " << idx << ", J = " << J << ", position = " << px[idx]<<endl;
        warn << "          Vel Grad = \n" << velGrad << endl; 
        warn << "          F_inc = \n" << deformationGradientInc << endl; 
        warn << "          F_old = \n" << Fold[idx] << endl; 
        warn << "          F_new = \n" << Fnew[idx] << endl; 
        warn << "          gVelocity:" << endl;
        for(int k = 0; k < flag->d_8or27; k++) {
          warn<< "             node: " << ni[k] << " vel: " << gVel[ni[k]] << endl;
        }
        
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }

    }
}

void
DeformationGradientComputer::computeDeformationGradientFromTotalDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           ParticleVariable<Matrix3> &Fnew,
                                           constParticleVariable<Matrix3>& Fold,
                                           Vector dx,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleInterpolator* interp)
{
  Matrix3 dispGrad,Identity;
  Identity.Identity();
  vector<IntVector> ni(interp->size());
  vector<double> S(interp->size());
  vector<Vector> d_S(interp->size());
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
                                                                                
  for(ParticleSubset::iterator iter = pset->begin();
       iter != pset->end(); iter++){
    particleIndex idx = *iter;
                                                                                
    // Get the node indices that surround the cell
    interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);
                                                                                
    computeGrad(dispGrad, ni, d_S, oodx, gDisp);
                                                                                
    // Update the deformation gradient tensor to its time n+1 value.
    // Compute the deformation gradient from the displacement gradient
    Fnew[idx] = Identity + dispGrad;
  }
}
                                                                                
void
DeformationGradientComputer::computeDeformationGradientFromIncrementalDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> Fold,
                                           ParticleVariable<Matrix3> &Fnew,
                                           Vector dx,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleInterpolator* interp)
{
    Matrix3 IncDispGrad,deformationGradientInc, Identity;
    Identity.Identity();
    vector<IntVector> ni(interp->size());
    vector<double> S(interp->size());
    vector<Vector> d_S(interp->size());

    double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
                                                                                
    for(ParticleSubset::iterator iter = pset->begin();
      iter != pset->end(); iter++){
      particleIndex idx = *iter;
                                                                                
      // Get the node indices that surround the cell
      interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);
                                                                                
      computeGrad(IncDispGrad, ni, d_S, oodx, gDisp);
                                                                                
      // Compute the deformation gradient increment
      deformationGradientInc = IncDispGrad + Identity;
                                                                                
      // Update the deformation gradient tensor to its time n+1 value.
      Fnew[idx] = deformationGradientInc * Fold[idx];
    }
}
