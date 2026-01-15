// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_NAVSTO_HPP
#define RBVMS_NAVSTO_HPP

#include "mfem.hpp"

using namespace mfem;

/** This class defines the time-dependent integrator for the
    Residual-based Variational multiscale formulation
    for incompressible Navier-Stokes flow.
*/
class IncNavStoIntegrator: public BlockTimeDepNonlinearFormIntegrator
{
private:

   // Physical parameters
   Coefficient &c_rho;
   Coefficient &c_mu;
   VectorCoefficient &c_force;
   VectorCoefficient &c_sol;
   Coefficient &c_suction;
   Coefficient &c_blowing;

   /// Numerical parameters
   real_t dt = -1.0;
   DenseMatrix Gij;
   Vector hn;

   /// Dimension data
   int dim = -1;
   Array2D<int> hmap;

   /// Physical values
   Vector u, dudt, f, grad_p, res_m, up, nor, traction;
   DenseMatrix flux;

   /// Solution & Residual vector
   DenseMatrix elf_u, elf_du, elv_u;

   /// Shape function data
   Vector sh_u, ushg_u, sh_p, dupdu;
   DenseMatrix shg_u, shh_u, shg_p, grad_u, hess_u;

   /// Compute RBVMS stabilisation parameters
   void GetTau(real_t &tau_m, real_t &tau_c, real_t &cfl2,
               real_t &rho, real_t &mu, Vector &u,
               ElementTransformation &Tr);

   /// Compute Weak Dirichlet stabilisation parameters
   void GetTauB(real_t &tau_b, real_t &tau_n,
                real_t &mu, Vector &u,
                Vector &nor,
                FaceElementTransformations &Tr);

public:
   /// Constructor
   IncNavStoIntegrator(Coefficient &rho_,
                       Coefficient &mu_,
                       VectorCoefficient &force_,
                       VectorCoefficient &sol_,
                       Coefficient &suction_,
                       Coefficient &blowing_);

   /// Set the timestep size @a dt_
   void SetTimeAndStep(const real_t &t, const real_t &dt_)
   {
      dt = dt_;
      c_mu.SetTime(t);
      c_force.SetTime(t);
      c_sol.SetTime(t);
      c_suction.SetTime(t);
      c_blowing.SetTime(t);
   };

   /// Assemble the local energy
   real_t GetElementEnergy(const Array<const FiniteElement *>&el,
                           ElementTransformation &Tr,
                           const Array<const Vector *> &elfun,
                           const Array<const Vector *> &elrate);

   /// Assemble the element constant artifical diffusion
   real_t GetElemArtDiff(const Array<const FiniteElement *> &el,
                         ElementTransformation &Tr,
                         const Array<const Vector *> &elsol,
                         const Array<const Vector *> &elrate);

   /// Assemble the element interior residual vectors
   void AssembleElementVector(const Array<const FiniteElement *> &el,
                              ElementTransformation &Tr,
                              const Array<const Vector *> &elsol,
                              const Array<const Vector *> &elrate,
                              const Array<Vector *> &elvec);
   //     real_t &cfl);

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
};

#endif
