// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef RBVMS_DIST_SOLVER_HPP
#define RBVMS_DIST_SOLVER_HPP

#include "mfem.hpp"

using namespace mfem;

/** This class manages the smooth interface functions,
    such as the smooth heaviside, smooth sign and smooth dirac functions.*/
class Heaviside
{
public:
   /// Small number to prevent division by zero
   static real_t epsilon;

   /// Half width of the interface in terms of elements
   static real_t eps;

   /// Compute element size
   static real_t h(Vector &grad_phi, DenseMatrix &Gij);
   static real_t h(Vector &grad_phi, ElementTransformation &Tr);

   /// Compute relative distance wrt zero level
   static real_t rphi(real_t &phi, real_t &h);

   /// Smooth step function, going from 0 to 1 in the interval [-1, 1]
   static real_t step(real_t &rphi);
   static real_t step(real_t &phi, Vector &grad_phi, DenseMatrix &Gij);
   static real_t step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);

   /// Smooth sign function, going from -1 to 1 in the interval [-1, 1]
   static real_t sign(real_t &rphi);
   static real_t sign(real_t &rphi, Vector &grad_phi, DenseMatrix &Gij);
   static real_t sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);

   /// Smooth dirac function, going from 0 to 0 in the interval [-1, 1]
   static real_t dirac(real_t &rphi, real_t &h);

   /// Integral of this function over phi is 1
   static real_t dirac(real_t &rphi, Vector &grad_phi, DenseMatrix &Gij);
   static real_t dirac(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
};

/** This Class defines an integrator for stabilized multi-dimensional
    convection-reaction equation.

     $(a \cdot \nabla u, v) + (k u,v)
    + \sum (a \cdot \nabla u + k u, \tau (a \cdot \nabla v + k v))_e - (f,  v)$
*/
class StabConvReactIntegrator : public NonlinearFormIntegrator
{
private:
   /// Penalty parameter for pinning the zero level-set
   real_t lambda;

   /// Discontinuity capturing parameter for the inconsistent part
   real_t kdc0;

   /// Discontinuity capturing parameter for the inconsistent part
   real_t kdc1;

   /// The function that defines the interface location
   const GridFunction *ls_gf;

   /// The stabilization parameter
   void SetDim(int dim);

   /// The stabilization parameter
   real_t GetTau(real_t &k, Vector &a, DenseMatrix &Gij);

   /// The disconituity capturing parameter
   real_t GetKdc(real_t &res, Vector &dphidx, DenseMatrix &Gij);

   /// Temporary variables
   int dim;
   Vector a, dphidx, dphidx0;
   Vector shape, trail, test;
   DenseMatrix dshape, Gij;

public:
   /// Constructor
   StabConvReactIntegrator(real_t l = 1.0,
                           real_t k0 = 0.01,
                           real_t k1 = 0.25)
   {
      lambda = l;
      kdc0 = k0;
      kdc1 = k1;
      dim = -1;
   };

   /// Destructor
   ~StabConvReactIntegrator() {};

   /// Set the penalty parameter for pinning the zero level-set
   void SetPenalty(real_t l) { lambda = l; };

   /// Set the penalty parameter for pinning the zero level-set
   void SetInconsistentDC(real_t k0) { kdc0 = k0; };

   /// Set the penalty parameter for pinning the zero level-set
   void SetConsistentDC(real_t k1) { kdc1 = k1; };

   /// Provide the Gridfunction that specifies the zero level-set
   void SetZeroLevelSet(const GridFunction *zero_level_set)
   {
      ls_gf = zero_level_set;
   }

   /// Compute the element energy --> volume in this case
   virtual real_t GetElementEnergy(const FiniteElement &el,
                                   ElementTransformation &Tr,
                                   const Vector &elfun) override;

   /// Compute the element nonlinear residual
   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun,
                                      Vector &elvect) override;

   /// Compute the element jacobian
   virtual void AssembleElementGrad(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun,
                                    DenseMatrix &elmat) override;

   /// Give element integration rule
   static const IntegrationRule &GetRule(const FiniteElement &trial_fe,
                                         const FiniteElement &test_fe,
                                         ElementTransformation &Trans);
};


/** This Class defines a redistancing algorithm based on convection.
    The convective equation
     $S_{\eps} (\phi_0)\frac{\nabla \phi}{\|\nabla \phi\|} \cdot \nabla \phi
      + \lambda \delta_{eps}(\phi_0) (\phi-phi_0)
      = S_{\eps} (\phi_0)$
    is solved using SUPG and Discontinuity capturing.
*/
class ConvectionDistanceSolver
// : public common::DistanceSolver
// This class could derive from DistanceSolver
// To avoid compile issues this is a stand alone class here
{
private:
   /// Abstract nonlinear formulation
   ParNonlinearForm form;

   /// Linear solver
   FGMRESSolver gmres;

   /// Nonlinear solver
   NewtonSolver newton_solver;

   /// Defines the element level integrals for the convection problem
   /// using SUPG and DC
   StabConvReactIntegrator integrator;

   /// Solver
   Solver *prec = nullptr;
   Vector zero, sol;

   // Volume conservation solve params
   int maxIter_vc = 10;
   real_t relTol_vc = 1e-12;
   real_t jacEps_vc = 1e-4;

   // Compute volume
   real_t ComputeVolume(ParGridFunction &distance);

   // Compute volume jacobian
   real_t ComputeVolumeJac(ParGridFunction &distance, real_t eps = 1e-8);

public:
   /// Constructor
   ConvectionDistanceSolver(ParFiniteElementSpace &space,
                            real_t lambda = 1.0);

   /// Destructor
   ~ConvectionDistanceSolver()
   {
      if (prec) { delete prec; }
   };

   // Set formulation parameters
   void SetPenalty(real_t lambda) { integrator.SetPenalty(lambda); };
   void SetInconsistentDC(real_t k0) { integrator.SetInconsistentDC(k0); };
   void SetConsistentDC(real_t k1) { integrator.SetConsistentDC(k1); };

   // Set linear solver parameters
   void SetLinearRelTol(real_t rtol) { gmres.SetRelTol(rtol); }
   void SetLinearAbsTol(real_t atol) { gmres.SetAbsTol(atol); }
   void SetLinearMaxIter(int maxiter) { gmres.SetMaxIter(maxiter); }

   void SetLinearPreconditioner(Solver &pc)
   {
      if (prec) { delete prec; }
      prec = nullptr;
      gmres.SetPreconditioner(pc);
   }

   // Set nonlinear solver parameters
   void SetNonlinearRelTol(real_t rtol) { newton_solver.SetRelTol(rtol); }
   void SetNonlinearAbsTol(real_t atol) { newton_solver.SetAbsTol(atol); }
   void SetNonlinearMaxIter(int maxiter) { newton_solver.SetMaxIter(maxiter); }
   void SetNonlinearMonitor(IterativeSolverMonitor& monitor)
   {
      newton_solver.SetPrintLevel(-1);
      newton_solver.SetMonitor(monitor);
   }

   // Set volume conservation solver parameters
   void SetVolumeConservationMaxIter(int im) { maxIter_vc = im; };
   void SetVolumeConservationRelTol(real_t tol) { relTol_vc = tol; };
   void SetVolumeConservationJacEps(real_t eps) { jacEps_vc = eps; };

   // Compute distance field for given level-set
   void ComputeScalarDistance(ParGridFunction &zero_level_set,
                              ParGridFunction &distance);

   void ComputeScalarDistance(Coefficient &zero_level_set,
                              ParGridFunction &distance);

   // Shift distance1 to have same volume as distance0
   void CorrectVolume(ParGridFunction &distance0,
                      ParGridFunction &distance1);
};


#endif
