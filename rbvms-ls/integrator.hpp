// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_LS_NAVSTO_HPP
#define RBVMS_LS_NAVSTO_HPP

#include "mfem.hpp"

using namespace mfem;

namespace RBVMS
{

/** This class defines the time-dependent integrator for the
    Residual-based Variational multiscale formulation
    for incompressible Navier-Stokes flow.
*/
class IncNavStoIntegrator
{
private:

   /// Physical coefficients
   Coefficient &c_rho;
   Coefficient &c_mu;
   VectorCoefficient &c_force;
   VectorCoefficient &c_sol_u;
   Coefficient &c_sol_phi;
   Coefficient &c_suction;
   Coefficient &c_blowing;

   real_t rho0 = 1.0;
   real_t rho1 = 1000.0;

   /// Numerical parameters
   real_t Cd = 6.0;
   real_t Ct = 1.0;
   real_t Cb = 12.0;
   real_t Cn = 100.0;

   real_t kdc0 = 0.0;
   real_t kdc1 = 0.1;

   /// Discretization parameters
   real_t dt = -1.0;
   DenseMatrix Gij;
   Vector hn;

   /// Dimension data
   int dim = -1;
   Array2D<int> hmap;

   /// Physical values
   Vector u, dudt, f, grad_p, grad_phi, res_m, up, nor, traction;
   DenseMatrix flux;

   /// Solution & Residual vector
   DenseMatrix elf_u, elf_du, elv_u;

   /// Shape function data
   Vector sh_u, ushg_u, sh_p, sh_phi, ushg_phi, dupdu;
   DenseMatrix shg_u, shh_u, shg_p, shg_phi, grad_u, hess_u;

   /// Compute RBVMS stabilisation parameters
   void GetTau(real_t &tau_m, real_t &tau_c, real_t &tau_ls, real_t &cfl2,
               real_t &rho, real_t &mu, Vector &u,
               DenseMatrix &Gij);

   /// Compute Weak Dirichlet stabilisation parameters
   void GetTauB(real_t &tau_b, real_t &tau_n,
                real_t &mu, Vector &u,
                Vector &nor,
                DenseMatrix &Gij);

   /// Compute discontinuity capturing parameters
   real_t GetKdc(Vector &res,
                 DenseMatrix &grad_u,
                 DenseMatrix &Gij);

   /// Compute density
   real_t GetRho(real_t &phi, Vector &grad_phi, DenseMatrix &Gij);

   /// Compute density gradient
   real_t GetRhoGrad(real_t &phi, Vector &grad_phi, DenseMatrix &Gij);

public:
   /// Constructor
   IncNavStoIntegrator(Coefficient &rho_,
                       Coefficient &mu_,
                       VectorCoefficient &force_,
                       VectorCoefficient &sol_u,
                       Coefficient &sol_phi,
                       Coefficient &suction_,
                       Coefficient &blowing_);

   /// Set densities of the two fluids
   void SetDensities(real_t r0, real_t r1) { rho0 = r0; rho1 = r1; };

   /// Set tau parameters
   void SetTauParams(real_t Cd_, real_t Ct_) { Cd = Cd_; Ct = Ct_; };

   /// Set wbc parameters
   void SetWBCParams(real_t Cb_, real_t Cn_) { Cb = Cb_; Cn = Cn_; };

   /// Set kbc parameters
   void SetKDCParams(real_t k0, real_t k1) { kdc0 = k0; kdc1 = k1; };

   void SetInconsistentDC(real_t k0) { kdc0 = k0; };
   void SetConsistentDC(real_t k1) { kdc1 = k1; };

   /// Set the timestep size @a dt_
   void SetTimeAndStep(const real_t &t, const real_t &dt_)
   {
      dt = dt_;
      c_mu.SetTime(t);
      c_force.SetTime(t);
      c_sol_u.SetTime(t);
      c_sol_phi.SetTime(t);
      c_suction.SetTime(t);
      c_blowing.SetTime(t);
   };

   /// Assemble the local energy
   void AssembleElementEnergy(const Array<const FiniteElement *>&el,
                              ElementTransformation &Tr,
                              const Array<const Vector *> &elfun,
                              Vector &energy);

   /// Assemble the element interior residual vectors
   void AssembleElementVector(const Array<const FiniteElement *> &el,
                              ElementTransformation &Tr,
                              const Array<const Vector *> &elsol,
                              const Array<const Vector *> &elrate,
                              const Array<Vector *> &elvec,
                              real_t &cfl);

   /// Assemble the element interior gradient matrices
   void AssembleElementGrad(const Array<const FiniteElement*> &el,
                            ElementTransformation &Tr,
                            const Array<const Vector *> &elsol,
                            const Array<const Vector *> &elrate,
                            const Array2D<DenseMatrix *> &elmats);

   /// Assemble the outflow boundary residual vectors
   void AssembleOutflowVector(const Array<const FiniteElement *> &el1,
                              const Array<const FiniteElement *> &el2,
                              FaceElementTransformations &Tr,
                              const Array<const Vector *> &elfun,
                              const Array<const Vector *> &elrate,
                              const Array<Vector *> &elvect,
                              real_t &outflow,
                              bool suction = false);

   /// Assemble the outflow boundary gradient matrices
   void AssembleOutflowGrad(const Array<const FiniteElement *>&el1,
                            const Array<const FiniteElement *>&el2,
                            FaceElementTransformations &Tr,
                            const Array<const Vector *> &elfun,
                            const Array<const Vector *> &elrate,
                            const Array2D<DenseMatrix *> &elmats,
                            bool suction = false);

   /// Assemble the weak Dirichlet BC boundary residual vectors
   void AssembleWeakDirBCVector(const Array<const FiniteElement *> &el1,
                                const Array<const FiniteElement *> &el2,
                                FaceElementTransformations &Tr,
                                const Array<const Vector *> &elfun,
                                const Array<const Vector *> &elrate,
                                const Array<Vector *> &elvect,
                                bool blowing = false);

   /// Assemble the weak Dirichlet BC boundary gradient matrices
   void AssembleWeakDirBCGrad(const Array<const FiniteElement *>&el1,
                              const Array<const FiniteElement *>&el2,
                              FaceElementTransformations &Tr,
                              const Array<const Vector *> &elfun,
                              const Array<const Vector *> &elrate,
                              const Array2D<DenseMatrix *> &elmats,
                              bool blowing = false);

   /// Assemble the normal Dirichlet BC boundary residual vectors
   void AssembleNormalBCVector(const Array<const FiniteElement *> &el1,
                               const Array<const FiniteElement *> &el2,
                               FaceElementTransformations &Tr,
                               const Array<const Vector *> &elfun,
                               const Array<const Vector *> &elrate,
                               const Array<Vector *> &elvect);

   /// Assemble the normal Dirichlet BC boundary gradient matrices
   void AssembleNormalBCGrad(const Array<const FiniteElement *>&el1,
                             const Array<const FiniteElement *>&el2,
                             FaceElementTransformations &Tr,
                             const Array<const Vector *> &elfun,
                             const Array<const Vector *> &elrate,
                             const Array2D<DenseMatrix *> &elmats);
};

} // namespace RBVMS

#endif
