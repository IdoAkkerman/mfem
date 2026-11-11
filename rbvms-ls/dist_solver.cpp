// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "dist_solver.hpp"

using namespace mfem;

// Default values
real_t Heaviside::epsilon = 1e-10;
real_t Heaviside::eps = 2.0;

// Compute element size
real_t Heaviside::h(Vector &grad_phi, DenseMatrix &Gij)
{
   int dim = grad_phi.Size();
   real_t gGg = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t gj = grad_phi[j];
      for (int i = 0; i < dim; i++)
      {
         gGg += Gij(i,j)*grad_phi[i]*gj;
      }
   }
   return grad_phi.Norml2()/sqrt(fmax(gGg, epsilon));
}

real_t Heaviside::h(Vector &grad_phi, ElementTransformation &Tr)
{
   // Metric tensor
   DenseMatrix Gij(grad_phi.Size());
   MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);

   return h(grad_phi, Gij);
}

// Compute relative distance wrt zero level
real_t Heaviside::rphi(real_t &phi, real_t &h)
{
   return phi/(eps*h);
}

// Smooth step function
real_t Heaviside::step(real_t &rphi)
{
   if (rphi < -1.0)
   {
      return 0.0;
   }
   else if (rphi > 1.0)
   {
      return 1.0;
   }
   else
   {
      return (1.0 + sin(M_PI*rphi/2))/2;
   }
}

real_t Heaviside::step(real_t &phi, Vector &grad_phi, DenseMatrix &Gij)
{
   real_t h = Heaviside::h(grad_phi, Gij);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::step(rphi);
}

real_t Heaviside::step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr)
{
   real_t h = Heaviside::h(grad_phi, Tr);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::step(rphi);
}

// Smooth sign function
real_t Heaviside::sign(real_t &rphi)
{
   return 2*Heaviside::step(rphi) - 1.0;
}

real_t Heaviside::sign(real_t &phi, Vector &grad_phi, DenseMatrix &Gij)
{
   return 2*Heaviside::step(phi, grad_phi, Gij) - 1.0;
}

real_t Heaviside::sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr)
{
   return 2*Heaviside::step(phi, grad_phi, Tr) - 1.0;
}

// Smooth dirac distribution
real_t Heaviside::dirac(real_t &rphi, real_t &h)
{
   if (rphi < -1.0)
   {
      return 0.0;
   }
   else if (rphi > 1.0)
   {
      return 0.0;
   }
   else
   {
      return cos(M_PI*rphi/2)*M_PI/(4*eps*h);
   }
}

real_t Heaviside::dirac(real_t &phi, Vector &grad_phi,
                        DenseMatrix &Gij)
{
   real_t h = Heaviside::h(grad_phi, Gij);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::dirac(rphi, h);
}

real_t Heaviside::dirac(real_t &phi, Vector &grad_phi,
                        ElementTransformation &Tr)
{
   real_t h = Heaviside::h(grad_phi, Tr);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::dirac(rphi, h);
}

// Define the integration rule based on FE order
const IntegrationRule &StabConvReactIntegrator::GetRule(
   const FiniteElement &trial_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans)
{
   int order = trial_fe.GetOrder() + test_fe.GetOrder();
   return IntRules.Get(trial_fe.GetGeomType(), order);
}

// Define convective tau
real_t StabConvReactIntegrator::GetTau(real_t &k, Vector &a, DenseMatrix &Gij)
{
   // Reaction part
   real_t  tau_i2 = k*k;

   // Convective part
   int dim = Gij.Width();
   for (int j = 0; j < dim; j++)
   {
      real_t aj = a[j];
      for (int i = 0; i < dim; i++)
      {
         tau_i2 += Gij(i,j)*a[i]*aj;
      }
   }

   // Momentum stabilisation parameter
   return 1.0/fmax(sqrt(tau_i2), Heaviside::epsilon);
}

// Define convective kdc
real_t StabConvReactIntegrator::GetKdc(real_t &res,
                                       Vector &dphidx,
                                       DenseMatrix &Gij)
{
   real_t h = 1.0/sqrt(Gij.Trace()/Gij.Width());
   return kdc0*h + kdc1*h*fabs(res)/fmax(dphidx.Norml2(), Heaviside::epsilon);
}

// Set the dimension of temporary variables
void StabConvReactIntegrator::SetDim(int d)
{
   if (dim != d )
   {
      dim = d;
      Gij.SetSize(dim);
      a.SetSize(dim);
      dphidx.SetSize(dim);
      dphidx0.SetSize(dim);
   }
}

// Compute the element energy --> volume in this case
real_t StabConvReactIntegrator::GetElementEnergy(const FiniteElement &el,
                                                 ElementTransformation &Trans,
                                                 const Vector &elfun)
{
   real_t w,phi;
   real_t h,rphi;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   shape.SetSize(nd);
   dshape.SetSize(nd,dim);

   const IntegrationRule *ir = NonlinearFormIntegrator::IntRule ?
                               NonlinearFormIntegrator::IntRule : &GetRule(el, el, Trans);

   real_t elvol = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;
      MultAtB(Trans.InverseJacobian(),Trans.InverseJacobian(),Gij);

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);
      phi = shape*elfun;

      // Calculate shapes
      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      // Interface params
      h = Heaviside::h(dphidx, Gij);
      rphi = Heaviside::rphi(phi,h);

      // Integral
      elvol += w*Heaviside::step(rphi);
   }
   return elvol;
}

// Compute the element nonlinear residual
void StabConvReactIntegrator::AssembleElementVector(const FiniteElement &el,
                                                    ElementTransformation &Trans,
                                                    const Vector &elfun,
                                                    Vector &elvect)
{
   real_t w,k,tau,f,phi,phi0,res,kdc;
   real_t h,rphi,Se,de;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elvect.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   test.SetSize(nd);

   const IntegrationRule *ir = NonlinearFormIntegrator::IntRule ?
                               NonlinearFormIntegrator::IntRule : &GetRule(el, el, Trans);

   elvect = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;
      MultAtB(Trans.InverseJacobian(),Trans.InverseJacobian(),Gij);

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);
      phi = shape*elfun;

      // Calculate shapes
      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      // Interface params
      phi0 = ls_gf->GetValue(Trans, ip);
      ls_gf-> GetGradient(Trans, dphidx0);

      h = Heaviside::h(dphidx0, Gij);
      rphi = Heaviside::rphi(phi0,h);
      Se = Heaviside::sign(rphi);
      de = Heaviside::dirac(rphi,h);

      // Force, reaction and convection parameters
      f = Se + lambda*de*phi0;
      k = lambda*de;
      a.Set(Se/fmax(dphidx.Norml2(),Heaviside::epsilon), dphidx);

      // Strong residual
      res = a*dphidx + k*phi - f;

      // Numerical parameters
      tau = GetTau(k, a, Gij);
      kdc = GetKdc(res, dphidx, Gij);

      // Compute test function
      dshape.Mult(a, test); // Add Convection stabilization term
      test.Add(k, shape);   // Add Reaction stabilization term
      test *= tau;          // Scale stabilization term with parameter
      test += shape;        // Galerkin term

      // Weak residual
      elvect.Add(w*res, test);

      // Artificial diffusion
      dshape.Mult(dphidx, test);
      elvect.Add(w*kdc, test);
   }
}

// Compute the element jacobian
void StabConvReactIntegrator::AssembleElementGrad(const FiniteElement &el,
                                                  ElementTransformation &Trans,
                                                  const Vector &elfun,
                                                  DenseMatrix &elmat)
{
   real_t w,k,f,phi,phi0,res,tau,kdc;
   real_t h,rphi,Se,de;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elmat.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   trail.SetSize(nd);
   test.SetSize(nd);

   const IntegrationRule *ir = NonlinearFormIntegrator::IntRule ?
                               NonlinearFormIntegrator::IntRule : &GetRule(el, el, Trans);

   elmat = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;
      MultAtB(Trans.InverseJacobian(),Trans.InverseJacobian(),Gij);

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);
      phi = shape*elfun;

      // Calculate shapes
      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      // Interface params
      phi0 = ls_gf->GetValue(Trans, ip);
      ls_gf-> GetGradient(Trans, dphidx0);

      h = Heaviside::h(dphidx0, Gij);
      rphi = Heaviside::rphi(phi0,h);
      Se = Heaviside::sign(rphi);
      de = Heaviside::dirac(rphi,h);

      // Force, reaction and convection parameters
      f = Se + lambda*de*phi0;
      k = lambda*de;
      a.Set(Se/fmax(dphidx.Norml2(), Heaviside::epsilon), dphidx);

      // Strong residual
      res = a*dphidx + k*phi - f;

      // Numerical parameters
      tau = GetTau(k, a, Gij);
      kdc = GetKdc(res, dphidx, Gij);

      // Compute trail function
      dshape.Mult(a, trail); // Add Convection term
      trail.Add(k, shape);   // Add Reaction term

      // Compute test function
      test.Set(tau, trail);   // Add stabilization term
      test += shape;          // Add Galerkin term

      // Stablization term
      AddMult_a_VWt(w, test, trail, elmat);

      // Artificial diffusion
      AddMult_a_AAt(w*kdc, dshape, elmat);
   }
}

ConvectionDistanceSolver::ConvectionDistanceSolver(ParFiniteElementSpace &space,
                                                   real_t lambda)
   : form(&space), gmres(space.GetComm()),newton_solver(space.GetComm())
{
   form.AddDomainIntegrator(&integrator);
   form.UseExternalIntegrators();

   // Set up the Jacobian solver
   gmres.iterative_mode = false;
   gmres.SetPrintLevel(-1);

   // Default values
   gmres.SetRelTol(1e-3);
   gmres.SetAbsTol(1e-12);
   gmres.SetMaxIter(100);
   prec = new HypreSmoother();
   gmres.SetPreconditioner(*prec);

   // Set up the Newton solver
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(1);
   newton_solver.SetSolver(gmres);
   newton_solver.SetOperator(form);

   // Default values
   newton_solver.SetRelTol(1e-3);
   newton_solver.SetAbsTol(1e-12);
   newton_solver.SetMaxIter(10);
}

// Solver for the distance field given a zero level-set coefficient
void ConvectionDistanceSolver::ComputeScalarDistance(ParGridFunction &phi0,
                                                     ParGridFunction &distance)
{
   integrator.SetZeroLevelSet(&phi0);

   distance.GetTrueDofs(sol);
   newton_solver.Mult(zero, sol);
   distance.Distribute(sol);
}

// Solver for the distance field given a zero level-set coefficient
void ConvectionDistanceSolver::ComputeScalarDistance(Coefficient &phi0,
                                                     ParGridFunction &distance)
{
   GridFunctionCoefficient* cf = dynamic_cast<GridFunctionCoefficient*> (&phi0);
   GridFunction* gf = const_cast<GridFunction*>(cf->GetGridFunction());
   ParGridFunction* phi0_pgf = dynamic_cast<ParGridFunction*>(gf);
   ComputeScalarDistance(*phi0_pgf, distance);
}

// Shift distance1 to have same volume as distance0
void ConvectionDistanceSolver::CorrectVolume(ParGridFunction &distance0,
                                             ParGridFunction &distance1)
{
   // Get reference vol
   real_t vol0 = ComputeVolume(distance0);

   // Get current volume
   real_t vol1 = ComputeVolume(distance1);

   // Get jacobian
   real_t jac = ComputeVolumeJac(distance1, jacEps_vc);

   // Print header
   if (Mpi::Root())
   {
      mfem::out<<"\n Volume conservation \n"
               <<std::setw(6) <<" Iter"
               <<std::setw(10)<<" Volume"
               <<std::setw(10)<<" |V-V0|"
               <<std::setw(12)<<" |V-V0|/|V0|"
               <<std::setw(10)<<" |dphi|"<<std::endl;
   }
   // Correct
   for (int it = 0; it < maxIter_vc; it++)
   {
      // Print convergence info
      if (Mpi::Root())
      {
         mfem::out<<std::setw(6)<<it
                  <<std::setw(10)<<std::defaultfloat<<std::setprecision(4)
                  <<vol1
                  <<std::setw(10)<<std::defaultfloat<<std::setprecision(2)
                  <<fabs(vol1 - vol0)
                  <<std::setw(12)<<std::defaultfloat<<std::setprecision(2)
                  <<fabs((vol1 - vol0)/vol0)
                  <<std::setw(10)<<std::defaultfloat<<std::setprecision(2)
                  <<fabs((vol1 - vol0)/jac)<<std::endl;
      }

      // Correct, check and recompute
      distance1 -= (vol1 - vol0)/jac;
      if (fabs(vol1 - vol0) < relTol_vc*vol0) { break; }
      vol1 = ComputeVolume(distance1);
   }
}

// Compute volume
real_t ConvectionDistanceSolver::ComputeVolume(ParGridFunction &distance)
{
   integrator.SetZeroLevelSet(&distance);
   return form.GetEnergy(distance);
}

// Compute volume jacobian using finite difference
real_t ConvectionDistanceSolver::ComputeVolumeJac(ParGridFunction &distance,
                                                  real_t eps)
{
   real_t vol0 = ComputeVolume(distance);
   distance += eps;
   real_t vol1 = ComputeVolume(distance);
   distance -= eps;

   return (vol1 - vol0)/eps;
}
