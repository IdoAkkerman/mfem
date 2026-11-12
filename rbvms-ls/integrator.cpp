// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "integrator.hpp"
#include "dist_solver.hpp"

using namespace mfem;
using namespace RBVMS;

// Constructor
IncNavStoIntegrator::IncNavStoIntegrator(Coefficient &rho,
                                         Coefficient &mu,
                                         VectorCoefficient &force,
                                         VectorCoefficient &sol_u,
                                         Coefficient &sol_phi,
                                         Coefficient &suction,
                                         Coefficient &blowing)
   : c_rho(rho), c_mu(mu), c_force(force), c_sol_u(sol_u), c_sol_phi(sol_phi),
     c_suction(suction), c_blowing(blowing)
{
   dim = force.GetVDim();
   u.SetSize(dim);
   dudt.SetSize(dim);
   f.SetSize(dim);
   res_m.SetSize(dim);
   up.SetSize(dim);
   traction.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
   grad_phi.SetSize(dim);
   Gij.SetSize(dim);
   nor.SetSize(dim);
   hn.SetSize(dim);
   hmap.SetSize(dim,dim);

   if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) =  hmap(1,0) =  1;
      hmap(1,1) = 2;
   }
   else if (dim == 3)
   {
      hmap(0,0) = 0;
      hmap(0,1) = hmap(1,0) = 1;
      hmap(0,2) = hmap(2,0) = 2;
      hmap(1,1) = 3;
      hmap(1,2) = hmap(2,1) = 4;
      hmap(2,2) = 5;
   }
   else
   {
      mfem_error("Only implemented for 2D and 3D");
   }
}

// Compute RBVMS stabilisation parameters
void IncNavStoIntegrator::GetTau(real_t &tau_m, real_t &tau_c,
                                 real_t &tau_ls, real_t &cfl2,
                                 real_t& rho, real_t &mu, Vector &u_,
                                 DenseMatrix &Gij_)
{
   /*   std::cout<<"rho = "<<rho<<std::endl;
      std::cout<<"mu = "<<mu<<std::endl;
      std::cout<<"u = ";  u_.Print();
      std::cout<<"u = ";  u.Print();
      std::cout<<"Gij = ";  Gij_.Print();
      std::cout<<"Gij = ";  Gij.Print();

      std::cout<<"Ct = "<<Ct<<std::endl;
      std::cout<<"Cd = "<<Cd<<std::endl;
      std::cout<<"dim = "<<dim<<std::endl;

      std::cout<<"dt= "<<dt<<std::endl;*/



   // Temporal part
   tau_m = Ct/(dt*dt);

   // Convective part
   tau_ls = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t uj = u_[j];
      for (int i = 0; i < dim; i++)
      {
         tau_ls += Gij_(i,j)*u_[i]*uj;
      }
   }
   cfl2 = tau_ls/tau_m;

   //  std::cout<<"cfl2= "<< cfl2<<std::endl;

   tau_ls += tau_m;
   tau_m = rho*rho*tau_ls;

   // Diffusive part
   real_t tmp = Cd*Cd*mu*mu;
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_m  += tmp*Gij_(i,j)*Gij_(i,j);
      }
   }

   // Momentum stabilisation parameter
   tau_m = 1.0/sqrt(tau_m);

   // Continuity stabilisation parameter
   tau_c = 1.0/(tau_m*Gij.Trace());

   // Levelset stabilisation parameter
   tau_ls = 1.0/(tau_ls);
}

// Compute Weak Dirichlet stabilisation parameters
void IncNavStoIntegrator::GetTauB(real_t &tau_b, real_t &tau_n,
                                  real_t &mu, Vector &u_,
                                  Vector &nor_,
                                  DenseMatrix &Gij_)
{
   real_t ho1 = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t nj = nor_[j];
      for (int i = 0; i < dim; i++)
      {
         ho1 += Gij_(i,j)*nor_[i]*nj;
      }
   }

   tau_b = Cb*mu*sqrt(ho1);
   tau_n = Cn;
}

// Compute discontinuity capturing parameters
real_t IncNavStoIntegrator::GetKdc(Vector &res,
                                   DenseMatrix &grad_u_,
                                   DenseMatrix &Gij_)
{
   real_t h = 1.0/sqrt(Gij_.Trace()/Gij_.Width());
   return kdc0*h*h*grad_u_.FNorm() +
          kdc1*h*res.Norml2()/fmax(grad_u_.FNorm(), Heaviside::epsilon);
}

// Compute density
real_t IncNavStoIntegrator::GetRho(real_t &phi, Vector &grad_phi_,
                                   DenseMatrix &Gij_)
{
   real_t He = Heaviside::step(phi, grad_phi_, Gij_);
   return rho0 + (rho1-rho0)*He;
}

// Compute density gradient
real_t IncNavStoIntegrator::GetRhoGrad(real_t &phi, Vector &grad_phi_,
                                       DenseMatrix &Gij_)
{
   real_t de = Heaviside::dirac(phi, grad_phi_, Gij_);
   return  (rho1-rho0)*de;
}

// Assemble the element energy
void IncNavStoIntegrator::AssembleElementEnergy(
   const Array<const FiniteElement *>&el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   Vector &energy)
{
   if (el.Size() != 3)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();
   int dof_phi = el[2]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   //bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }

   energy.SetSize(3);
   energy = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   //real_t tau_m, tau_c, tau_ls, cfl2;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);
      real_t mu = c_mu.Eval(Tr, ip);

      c_force.Eval(f, Tr, ip);

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      //real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      el[2]->CalcPhysShape(Tr, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el[2]->CalcPhysDShape(Tr, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      Vector x;
      Tr.Transform(ip, x);

      // Compute strong residual
      energy[0] += (rho*(u*u)/2)*w;
      energy[1] -= rho*(x*f)*w;
      grad_u.Symmetrize();
      energy[2] += 2*mu*grad_u.FNorm2()*w;
   }
}

// Assemble the element interior residual vectors
void IncNavStoIntegrator::AssembleElementVector(
   const Array<const FiniteElement *> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array<Vector *> &elvec,
   real_t &elem_cfl)
{
   if (el.Size() != 3)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();
   int dof_phi = el[2]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }
   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);
   elvec[2]->SetSize(dof_phi);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, tau_ls, cfl2;
   elem_cfl = 0.0;


   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);
      real_t mu = c_mu.Eval(Tr, ip);
      c_force.Eval(f, Tr, ip);

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      el[2]->CalcPhysShape(Tr, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);
      real_t dphidt = sh_phi*(*elrate[2]);

      el[2]->CalcPhysDShape(Tr, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m -= f;             // Subtract force
      res_m *= rho;
      res_m += grad_p;        // Add pressure

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int ii = 0; ii < dim; ++ii)
         {
            for (int jj = 0; jj < dim; ++jj)
            {
               res_m[jj] -= mu*(hess_u(jj,hmap(ii,ii)) +
                                hess_u(ii,hmap(jj,ii)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }
      real_t res_c = grad_u.Trace();

      // Compute stability params
      // std::cout<<"rho = "<<rho<<std::endl;
      // std::cout<<"mu = "<<mu<<std::endl;
      // std::cout<<"u = ";  u.Print();
      // std::cout<<"Gij = ";  Gij.Print();

      GetTau(tau_m, tau_c, tau_ls, cfl2, rho, mu, u, Gij);
      elem_cfl = fmax(elem_cfl, cfl2);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      // std::cout<<"tau_m = "<<tau_m<<std::endl;
      // std::cout<<"res_m = "; res_m.Print();
      //std::cout<<"up = "; up.Print();
      //  u += up;
      p -= tau_c*res_c;

      mu += GetKdc(res_m, grad_u, Gij);

      // Compute momentum weak residual
      grad_u.Mult(u,res_m);                       // Add convection (incl cross)
      res_m += dudt;                              // Add acceleration
      res_m -= f;                                 // Add force
      AddMult_a_VWt(w*rho, sh_u, res_m, elv_u);   // Add force + acc term to rhs

      flux.Diag(-p, dim);                         // Add pressure
      grad_u.Symmetrize();                        // Grad to strain
      flux.Add(2*mu,grad_u);                      // Add stress to flux
      AddMult_a_VWt(-rho, up, u, flux);           // Add SUPG to flux (incl rey)
      AddMult_a_ABt(w, shg_u, flux, elv_u);       // Add flux term to rhs

      // Compute continuity weak residual
      elvec[1]->Add(-w*res_c, sh_p);              // Add Galerkin term
      shg_p.Mult(up, sh_p);                       // PSPG help term
      elvec[1]->Add(w, sh_p);                     // Add PSPG term

      // Compute convection weak residual
      real_t res_ls = dphidt + u*grad_phi;
      elvec[2]->Add(w*res_ls, sh_phi);            // Add Galerkin term
      elvec[2]->Add(w*res_ls*tau_ls, ushg_phi);   // Add SUPG term
   }
   //elvec[0]->Print(std::cout,999);
   //elvec[1]->Print(std::cout,999);
   //elvec[2]->Print(std::cout,999);
   elem_cfl = sqrt(elem_cfl);
}

// Assemble the element interior gradient matrices
void IncNavStoIntegrator::AssembleElementGrad(
   const Array<const FiniteElement*> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();
   int dof_phi = el[2]->GetDof();

   bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_wphi = *elmats(0,2);

   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);
   DenseMatrix &mat_qphi = *elmats(1,2);

   DenseMatrix &mat_vu = *elmats(2,0);
   DenseMatrix &mat_vp = *elmats(2,1);
   DenseMatrix &mat_vphi = *elmats(2,2);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_wphi.SetSize(dof_u*dim, dof_phi);

   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);
   mat_qphi.SetSize(dof_p, dof_phi);

   mat_vu.SetSize(dof_p, dof_u*dim);
   mat_vp.SetSize(dof_p, dof_p);
   mat_vphi.SetSize(dof_p, dof_phi);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_wphi = 0.0;

   mat_qu = 0.0;
   mat_qp = 0.0;
   mat_qphi = 0.0;

   mat_vu = 0.0;
   mat_vp = 0.0;
   mat_vphi = 0.0;

   DenseMatrix mat_wu1(dof_u, dof_u);
   mat_wu1 = 0.0;
   DenseMatrix mat_wp1[dim], mat_qu1[dim], mat_up[dim], mat_wphi1[dim];
   DenseMatrix mat_vu1[dim];

   for (int i_dim = 0; i_dim < dim; ++i_dim)
   {
      mat_wp1[i_dim].SetSize(dof_u,dof_p);
      mat_qu1[i_dim].SetSize(dof_p,dof_u);
      mat_up[i_dim].SetSize(dof_u,dof_p);
      mat_wphi1[i_dim].SetSize(dof_u,dof_phi);
      mat_vu1[i_dim].SetSize(dof_phi,dof_u);
      mat_wp1[i_dim] = 0.0;
      mat_qu1[i_dim] = 0.0;
      mat_up[i_dim] = 0.0;
      mat_wphi1[i_dim] = 0.0;
      mat_vu1[i_dim] = 0.0;
   }

   int dim2 = (dim*(dim+1))/2;
   DenseMatrix  mat_uu[dim*dim];
   for (int i_dim = 0; i_dim < dim2; ++i_dim)
   {
      mat_uu[i_dim].SetSize(dof_u,dof_u);
      mat_uu[i_dim] = 0.0;
   }

   Vector vec_u1(NULL,dof_u), vec_u2(NULL,dof_u), vec_p(NULL,dof_p);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_u);

   DenseMatrix shg_uT(dim, dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, tau_ls, cfl2;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);

      real_t mu = c_mu.Eval(Tr, ip);

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      //real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      el[2]->CalcPhysShape(Tr, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);
      real_t dphidt = sh_phi*(*elrate[2]);

      el[2]->CalcPhysDShape(Tr, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);
      real_t drho = GetRhoGrad(phi, grad_phi, Gij);

      // Compute strong residual
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m -= f;             // Subtract force
      res_m *= rho;
      res_m += grad_p;        // Add pressure

      /////
      Vector bla(dim);
      grad_u.Mult(u,bla);   // Add convection
      bla += dudt;          // Add acceleration
      bla -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int ii = 0; ii < dim; ++ii)
         {
            for (int jj = 0; jj < dim; ++jj)
            {
               res_m[jj] -= mu*(hess_u(jj,hmap(ii,ii)) +
                                hess_u(ii,hmap(jj,ii)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }

      // Compute stability params
      GetTau(tau_m, tau_c, tau_ls, cfl2, rho, mu, u, Gij);

      mu += GetKdc(res_m, grad_u, Gij);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      // u += up;

      // Compute small scale jacobian
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         dupdu(j_u) = -tau_m*rho*(sh_u(j_u) + ushg_u(j_u)*dt);
      }
      //shg_u.Mult(u, ushg_u);

      // Recompute convective gradient
      MultAtB(elf_u, shg_u, grad_u);
      shg_uT.Transpose(shg_u);

      // Level - set - Block diagonal level-set block (v,phi)

      //   real_t res_ls = dphidt + u*grad_phi;
      //   elvec[2]->Add(w*res_ls, sh_phi);            // Add Galerkin term
      //   elvec[2]->Add(w*res_ls*tau_ls, ushg_phi);   // Add SUPG term

      // Galerkin terms
      AddMult_a_VVt(w, sh_phi, mat_vphi);
      AddMult_a_VWt(dt*w, sh_phi, ushg_phi, mat_vphi);

      // SUPG terms
      AddMult_a_VWt(w*tau_ls, ushg_phi, sh_phi, mat_vphi);
      AddMult_a_VVt(dt*w*tau_ls, ushg_phi, mat_vphi);

      real_t res_ls = dphidt + u*grad_phi;

      Vector col(dof_phi);
      for (int i_dim = 0; i_dim < dim; ++i_dim)
      {
         real_t tmp = dt*w*grad_phi(i_dim);
         // AddMult_a_VWt(tmp,        sh_phi,   sh_u, mat_vu1[i_dim]);
         //AddMult_a_VWt(tmp*tau_ls, ushg_phi, sh_u, mat_vu1[i_dim]);
         AddMult_a_VWt(tmp,        sh_u,sh_phi,    mat_vu1[i_dim]);
         AddMult_a_VWt(tmp*tau_ls, sh_u,ushg_phi,  mat_vu1[i_dim]);

         shg_phi.GetColumn(i_dim,col);
         //   AddMult_a_VWt(dt*w*tau_ls*res_ls, col, sh_u, mat_vu1[i_dim]);
         AddMult_a_VWt(dt*w*tau_ls*res_ls, sh_u, col,  mat_vu1[i_dim]);
      }

      // Cross term Momentem Galerkin wrt phi
      for (int i_dim = 0; i_dim < dim; ++i_dim)
      {
         AddMult_a_VVt(drho*dt*w*bla[i_dim], sh_u, mat_wphi1[i_dim]);
      }

      // Momentum - Block diagonal Velocity block (w,u)
      // Acceleration term
      AddMult_a_VVt(rho*w, sh_u, mat_wu1);

      // Convection terms
      AddMult_a_VWt(dt*rho*w, sh_u, ushg_u, mat_wu1);
      AddMult_a_VWt(-rho*w, ushg_u, dupdu, mat_wu1);

      // Diffusion term
      AddMult_a_AAt(w*mu*dt, shg_u, mat_wu1);

      // Continuity - Pressure block (q,p)
      AddMult_a_AAt(-w*tau_m, shg_p, mat_qp);

      // Off diagional block terms
      int ii = 0;
      for (int i_dim = 0; i_dim < dim; ++i_dim)
      {
         // Getting columns for outer product
         shg_p.GetColumnReference(i_dim, vec_p);
         shg_u.GetColumnReference(i_dim, vec_u1);

         // Momentum + Continuity Galerkin terms
         AddMult_a_VWt(-w, vec_u1, sh_p, mat_up[i_dim]);

         // Momentum - Pressure block (w,p)
         AddMult_a_VWt(tau_m * w, ushg_u, vec_p, mat_wp1[i_dim] );

         // Continuity - Velocity block (q,u)
         AddMult_a_VWt(w, vec_p, dupdu, mat_qu1[i_dim]);

         // Momentum - Velocity block (w,u)
         for (int j_dim = 0; j_dim < i_dim; ++j_dim)
         {
            shg_u.GetColumnReference(j_dim, vec_u2);
            AddMult_a_VWt(w*dt*(mu+tau_c), vec_u1, vec_u2, mat_uu[ii++]);
         }
         AddMult_a_VVt(w*dt*(mu+tau_c), vec_u1, mat_uu[ii++]);
      }
   }

   // Adding sub elements
   int ii = 0;
   for (int i_dim = 0; i_dim < dim; ++i_dim)
   {
      mat_wp1[i_dim] += mat_up[i_dim];
      mat_up[i_dim].Transpose();
      mat_up[i_dim] *= dt;
      mat_qu1[i_dim] += mat_up[i_dim];

      mat_wp.AddSubMatrix(i_dim * dof_u, 0, mat_wp1[i_dim]);
      mat_qu.AddSubMatrix(0, i_dim * dof_u, mat_qu1[i_dim]);
      mat_wu.AddSubMatrix(i_dim*dof_u, mat_wu1);
      for (int j_dim = 0; j_dim < i_dim; ++j_dim)
      {
         mat_wu.AddSubMatrix(i_dim * dof_u, j_dim * dof_u, mat_uu[ii]);
         mat_uu[ii].Transpose();
         mat_wu.AddSubMatrix(j_dim * dof_u, i_dim * dof_u, mat_uu[ii++]);
      }
      mat_wu.AddSubMatrix(i_dim * dof_u, i_dim * dof_u, mat_uu[ii++]);

      mat_wphi.AddSubMatrix(i_dim * dof_u, 0, mat_wphi1[i_dim]);

      mat_vu.AddSubMatrix(0, i_dim * dof_u, mat_vu1[i_dim]);
   }

}

// Assemble the outflow boundary residual vectors
void IncNavStoIntegrator
::AssembleOutflowVector(const Array<const FiniteElement *> &el1,
                        const Array<const FiniteElement *> &el2,
                        FaceElementTransformations &Tr,
                        const Array<const Vector *> &elsol,
                        const Array<const Vector *> &elrate,
                        const Array<Vector *> &elvec,
                        real_t &outflow,
                        bool suction)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);
   elvec[2]->SetSize(dof_phi);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;
   outflow = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   //Vector
   traction.SetSize(dim);
   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor); // nor = n.da
      real_t w = ip.weight;          // No weight --> taken care of by nor
      // real_t hbc = c_hbc.Eval(Tr, ip);
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      real_t un = u*nor;
      if (suction)
      {
         real_t suction_val = c_suction.Eval(*Tr.Elem1, eip);
         traction.Set(suction_val,nor);
         AddMult_a_VWt(w, sh_u, traction, elv_u);
         outflow += rho * un * w; // No weight --> taken care of by nor
      }

      if (un > 0.0) { continue; }

      c_sol_u.Eval(up, *Tr.Elem1, eip);
      up -= u;
      // up.Neg();


      AddMult_a_VWt(rho*w*un, sh_u, up, elv_u);
      //elvec[2]->Add(w*un*(phi0-phi), sh_p);
   }
}

// Assemble the outflow boundary gradient matrices
void IncNavStoIntegrator
::AssembleOutflowGrad(const Array<const FiniteElement *>&el1,
                      const Array<const FiniteElement *>&el2,
                      FaceElementTransformations &Tr,
                      const Array<const Vector *> &elsol,
                      const Array<const Vector *> &elrate,
                      const Array2D<DenseMatrix *> &elmats,
                      bool suction)
{
   int dof_u = el1[0]->GetDof();
   //int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   elmats(0,1)->SetSize(0,0);
   elmats(0,2)->SetSize(0,0);

   elmats(1,0)->SetSize(0,0);
   elmats(1,1)->SetSize(0,0);
   elmats(1,2)->SetSize(0,0);

   elmats(2,0)->SetSize(0,0);
   elmats(2,1)->SetSize(0,0);
   elmats(2,2)->SetSize(0,0);

   mat_wu = 0.0;

   sh_u.SetSize(dof_u);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder()+2;
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      // const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor); // nor = n.da
      real_t w = ip.weight;          // No weight --> taken care of by nor
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      real_t un = u*nor;

      if (un > 0.0) { continue; }

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      // Momentum - Velocity block (w,u)
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp = -rho*sh_u(j_u)*un*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += sh_u(i_u)*tmp;
            }
         }
      }

      // Jacobian due to the normal velocity
      /*   for (int dim_v = 0; dim_v < dim; ++dim_v)
         {
            real_t tmp0 = rho*nor(dim_v)*w*dt;
            for (int j_u = 0; j_u < dof_u; ++j_u)
            {
               real_t tmp1 = tmp0*sh_u(j_u);
               for (int dim_u = 0; dim_u < dim; ++dim_u)
               {
                  real_t tmp2 = tmp1*u(dim_u);
                  for (int i_u = 0; i_u < dof_u; ++i_u)
                  {
                     mat_wu(i_u + dim_u*dof_u, j_u + dim_v*dof_u) += sh_u(i_u)*tmp2;
                  }
               }
            }
         }*/
   }
}

// Assemble the weak Dirichlet BC boundary residual vectors
void IncNavStoIntegrator
::AssembleWeakDirBCVector(const Array<const FiniteElement *> &el1,
                          const Array<const FiniteElement *> &el2,
                          FaceElementTransformations &Tr,
                          const Array<const Vector *> &elsol,
                          const Array<const Vector *> &elrate,
                          const Array<Vector *> &elvec,
                          bool blowing)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);
   elvec[2]->SetSize(dof_phi);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);

   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n, mu, w, phi,rho, un;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      mu = c_mu.Eval(*Tr.Elem1, eip);
      c_sol_u.Eval(up, *Tr.Elem1, eip);
      real_t phi0 = c_sol_phi.Eval(*Tr.Elem1, eip);

      w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      if (blowing)
      {
         real_t blowing_val = c_blowing.Eval(*Tr.Elem1, eip);
         up.Add(-blowing_val, nor);
      }

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      up -= u;
      up.Neg();
      un = up*nor;

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Symmetrize();           // Grad to strain

      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      phi = sh_phi*(*elsol[2]);
      //dphidt = sh_phi*(*elrate[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      rho = GetRho(phi, grad_phi, Gij);

      GetTauB(tau_b, tau_n, mu, u, nor, Gij);

      // Traction
      grad_u.Mult(nor, traction);
      traction *= -2*mu;             // Consistency
      traction.Add(p, nor);          // Pressure
      traction.Add(tau_b,up);        // Penalty
      traction.Add(tau_n*un,nor);    // Penalty -- normal
      AddMult_a_VWt(w, sh_u, traction, elv_u);

      // Dual consistency
      MultVWt(nor, up, flux);
      flux.Symmetrize();
      AddMult_a_ABt(-w*2*mu, shg_u, flux, elv_u);

      // Continuity
      elvec[1]->Add(w*un, sh_p);

      // Convection
      un = u*nor;
      if (un > 0.0) { continue; }

      AddMult_a_VWt(-rho*w*un, sh_u, up, elv_u);

      elvec[2]->Add(w*un*(phi0-phi), sh_p);
   }
}

// Assemble the weak Dirichlet BC boundary gradient matrices
void IncNavStoIntegrator
::AssembleWeakDirBCGrad(const
                        Array<const FiniteElement *>&el1,
                        const Array<const FiniteElement *>&el2,
                        FaceElementTransformations &Tr,
                        const Array<const Vector *> &elsol,
                        const Array<const Vector *> &elrate,
                        const Array2D<DenseMatrix *> &elmats,
                        bool blowing)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_wphi = *elmats(0,2);

   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);
   DenseMatrix &mat_qphi = *elmats(1,2);

   DenseMatrix &mat_vu = *elmats(2,0);
   DenseMatrix &mat_vp = *elmats(2,1);
   DenseMatrix &mat_vphi = *elmats(2,2);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_wphi.SetSize(dof_u*dim, dof_phi);

   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);
   mat_qphi.SetSize(dof_p, dof_phi);

   mat_vu.SetSize(dof_p, dof_u*dim);
   mat_vp.SetSize(dof_p, dof_p);
   mat_vphi.SetSize(dof_p, dof_phi);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_wphi = 0.0;

   mat_qu = 0.0;
   mat_qp = 0.0;
   mat_qphi = 0.0;

   mat_vu = 0.0;
   mat_vp = 0.0;
   mat_vphi = 0.0;

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n, mu, w, phi,rho;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();
      mu = c_mu.Eval(*Tr.Elem1, eip);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      phi = sh_phi*(*elsol[2]);
      //dphidt = sh_phi*(*elrate[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      rho = GetRho(phi, grad_phi, Gij);

      GetTauB(tau_b, tau_n, mu, u, nor, Gij);

      // Momentum - Velocity block (w,u)
      // Consistency & Dual: Part 1
      real_t tmp0 = mu*w*dt;
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            real_t tmp1 = nor(j_dim)*shg_u(j_u,j_dim)*tmp0;
            real_t tmp2 = nor(j_dim)*sh_u(j_u)*tmp0;
            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               int j_dof = j_u + i_dim*dof_u;
               for (int i_u = 0; i_u < dof_u; ++i_u)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  -= (sh_u(i_u)*tmp1 + shg_u(i_u,j_dim)*tmp2);
               }
            }
         }
      }

      // Consistency & Dual: Part 2
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            int j_dof = j_u + j_dim*dof_u;
            real_t tmp1 = nor(j_dim)*sh_u(j_u)*tmp0;
            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               real_t tmp2 = nor(j_dim)*shg_u(j_u,i_dim)*tmp0;
               for (int i_u = 0; i_u < dof_u; ++i_u)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  -= (sh_u(i_u)*tmp2 +shg_u(i_u,i_dim)*tmp1);
               }
            }
         }
      }

      // Penalty
      tmp0 = tau_b*w*dt;
      for (int dim_u = 0; dim_u < dim; ++dim_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            real_t tmp1 = sh_u(j_u)*tmp0;
            int j_dof = j_u + dim_u*dof_u;
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_dof) += sh_u(i_u)*tmp1;
            }
         }
      }

      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         tmp0 = sh_u(j_u)*w*dt*tau_n;
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            real_t tmp1 = tmp0*nor(j_dim);
            int j_dof = j_u + j_dim*dof_u;
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               for (int i_dim = 0; i_dim < dim; ++i_dim)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  += sh_u(i_u)*nor(i_dim)*tmp1;
               }
            }
         }
      }

      // Momentum - Pressure block (w,p)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         tmp0 = sh_p(i_p)*w;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            real_t tmp1 = nor(dim_u)*tmp0;
            for (int j_u = 0; j_u < dof_u; ++j_u)
            {
               mat_wp(j_u + dof_u * dim_u, i_p) += sh_u(j_u)*tmp1;
            }
         }
      }

      // Continuity - Velocity block (q,u)
      for (int dim_u = 0; dim_u < dim; ++dim_u)
      {
         tmp0 = nor(dim_u)*w*dt;
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            real_t tmp1 = sh_u(j_u)*tmp0;
            for (int i_p = 0; i_p < dof_p; ++i_p)
            {
               mat_qu(i_p, j_u + dof_u * dim_u) += sh_p(i_p)*tmp1;
            }
         }
      }

      // Convection: Momentum - Velocity block (w,u)
      real_t un = u*nor;
      if (un > 0.0) { continue; }
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp = -rho*sh_u(j_u)*un*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += sh_u(i_u)*tmp;
            }
         }
      }

      //
      for (int jj = 0; jj < dof_phi; ++jj)
      {
         // real_t tmp = -un*w*dt;
         for (int ii = 0; ii < dof_phi; ++ii)
         {
            mat_vphi(ii, jj) -= sh_phi(ii)*sh_phi(jj)*un*w*dt;
         }
      }

   }
}



// Assemble the normal Dirichlet BC boundary residual vectors
void IncNavStoIntegrator
::AssembleNormalBCVector(const Array<const FiniteElement *> &el1,
                         const Array<const FiniteElement *> &el2,
                         FaceElementTransformations &Tr,
                         const Array<const Vector *> &elsol,
                         const Array<const Vector *> &elrate,
                         const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   sh_u.SetSize(dof_u);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_n = 10000.0; // TO BE MOVED
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      c_sol_u.Eval(up, *Tr.Elem1, eip);

      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      up -= u;
      up.Neg();
      real_t un = up*nor;
      AddMult_a_VWt(w*rho*tau_n*un, sh_u, nor, elv_u);
   }
}

// Assemble the normal Dirichlet BC boundary gradient matrices
void IncNavStoIntegrator
::AssembleNormalBCGrad(const Array<const FiniteElement *>&el1,
                       const Array<const FiniteElement *>&el2,
                       FaceElementTransformations &Tr,
                       const Array<const Vector *> &elsol,
                       const Array<const Vector *> &elrate,
                       const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_wphi = *elmats(0,2);

   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);
   DenseMatrix &mat_qphi = *elmats(1,2);

   DenseMatrix &mat_vu = *elmats(2,0);
   DenseMatrix &mat_vp = *elmats(2,1);
   DenseMatrix &mat_vphi = *elmats(2,2);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_wphi.SetSize(dof_u*dim, dof_phi);

   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);
   mat_qphi.SetSize(dof_p, dof_phi);

   mat_vu.SetSize(dof_p, dof_u*dim);
   mat_vp.SetSize(dof_p, dof_p);
   mat_vphi.SetSize(dof_p, dof_phi);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_wphi = 0.0;

   mat_qu = 0.0;
   mat_qp = 0.0;
   mat_qphi = 0.0;

   mat_vu = 0.0;
   mat_vp = 0.0;
   mat_vphi = 0.0;

   sh_u.SetSize(dof_u);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_n = 10000.0; // TO BE MOVED
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      //const IntegrationPoint &eip = Tr.GetElement1IntPoint();
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      // Penalty
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp0 = sh_u(j_u)*w*dt*rho*tau_n;
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            real_t tmp1 = tmp0*nor(j_dim);
            int j_dof = j_u + j_dim*dof_u;
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               for (int i_dim = 0; i_dim < dim; ++i_dim)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  += sh_u(i_u)*nor(i_dim)*tmp1;
               }
            }
         }
      }
   }
}
























// Constructor
RBVMSIntegrator::RBVMSIntegrator(Coefficient &rho,
                                 Coefficient &mu,
                                 VectorCoefficient &force)
   : c_rho(rho), c_mu(mu), c_force(force)
{
   dim = force.GetVDim();
   u.SetSize(dim);
   dudt.SetSize(dim);
   f.SetSize(dim);
   res_m.SetSize(dim);
   up.SetSize(dim);
   traction.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
   grad_phi.SetSize(dim);
   Gij.SetSize(dim);
   nor.SetSize(dim);
   hn.SetSize(dim);
   hmap.SetSize(dim,dim);

   if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) =  hmap(1,0) =  1;
      hmap(1,1) = 2;
   }
   else if (dim == 3)
   {
      hmap(0,0) = 0;
      hmap(0,1) = hmap(1,0) = 1;
      hmap(0,2) = hmap(2,0) = 2;
      hmap(1,1) = 3;
      hmap(1,2) = hmap(2,1) = 4;
      hmap(2,2) = 5;
   }
   else
   {
      mfem_error("Only implemented for 2D and 3D");
   }
}

// Compute RBVMS stabilisation parameters
void RBVMSIntegrator::GetTau(real_t &tau_m, real_t &tau_c,
                             real_t &tau_ls, real_t &cfl2,
                             real_t& rho, real_t &mu, Vector &u_,
                             DenseMatrix &Gij_)
{
   /*  std::cout<<"rho = "<<rho<<std::endl;
     std::cout<<"mu = "<<mu<<std::endl;
     std::cout<<"u = ";  u_.Print();
     std::cout<<"u = ";  u.Print();
     std::cout<<"Gij = ";  Gij_.Print();
     std::cout<<"Gij = ";  Gij.Print();

     std::cout<<"Ct = "<<Ct<<std::endl;
     std::cout<<"Cd = "<<Cd<<std::endl;
     std::cout<<"dim = "<<dim<<std::endl;
     std::cout<<"dt= "<<dt<<std::endl;*/

   // Temporal part
   tau_m = Ct/(dt*dt);

   // Convective part
   tau_ls = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t uj = u_[j];
      for (int i = 0; i < dim; i++)
      {
         tau_ls += Gij_(i,j)*u_[i]*uj;
      }
   }
   cfl2 = tau_ls/tau_m;
   //      std::cout<<"cfl2= "<<cfl2<<std::endl;
   tau_ls += tau_m;
   tau_m = rho*rho*tau_ls;

   // Diffusive part
   real_t tmp = Cd*Cd*mu*mu;
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_m  += tmp*Gij_(i,j)*Gij_(i,j);
      }
   }

   // Momentum stabilisation parameter
   tau_m = 1.0/sqrt(tau_m);

   // Continuity stabilisation parameter
   tau_c = 1.0/(tau_m*Gij.Trace());

   // Levelset stabilisation parameter
   tau_ls = 1.0/(tau_ls);
}


// Compute discontinuity capturing parameters
real_t RBVMSIntegrator::GetKdc(Vector &res,
                               DenseMatrix &grad_u_,
                               DenseMatrix &Gij_)
{
   real_t h = 1.0/sqrt(Gij_.Trace()/Gij_.Width());
   return kdc0*h*h*grad_u_.FNorm() +
          kdc1*h*res.Norml2()/fmax(grad_u_.FNorm(), Heaviside::epsilon);
}

// Compute density
real_t RBVMSIntegrator::GetRho(real_t &phi, Vector &grad_phi_,
                               DenseMatrix &Gij_)
{
   real_t He = Heaviside::step(phi, grad_phi_, Gij_);
   return rho0 + (rho1-rho0)*He;
}

// Compute density gradient
real_t RBVMSIntegrator::GetRhoGrad(real_t &phi, Vector &grad_phi_,
                                   DenseMatrix &Gij_)
{
   real_t de = Heaviside::dirac(phi, grad_phi_, Gij_);
   return  (rho1-rho0)*de;
}

// Assemble the element energy
/*void RBVMSIntegrator::AssembleElementEnergy(
   const Array<const FiniteElement *>&el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol, energy)
{
   if (el.Size() != 3)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();
   int dof_phi = el[2]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   //bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }

   energy.SetSize(3);
   energy = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   //real_t tau_m, tau_c, tau_ls, cfl2;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);
      real_t mu = c_mu.Eval(Tr, ip);

      c_force.Eval(f, Tr, ip);

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      el[2]->CalcPhysShape(Tr, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el[2]->CalcPhysDShape(Tr, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      Vector x;
      Tr.Transform(ip, x);

      // Compute strong residual
      energy[0] += (rho*(u*u)/2)*w;
      energy[1] -= rho*(x*f)*w;
      grad_u.Symmetrize();
      energy[2] += 2*mu*grad_u.FNorm2()*w;
   }
}*/

// Assemble the element interior residual vectors
void RBVMSIntegrator::AssembleElementVector(
   const Array<const FiniteElement *> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array<Vector *> &elvec)
{
   real_t elem_cfl;   // TBD

   if (el.Size() != 3)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();
   int dof_phi = el[2]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }
   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);
   elvec[2]->SetSize(dof_phi);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, tau_ls, cfl2;
   elem_cfl = 0.0;


   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);
      real_t mu = c_mu.Eval(Tr, ip);
      c_force.Eval(f, Tr, ip);

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      el[2]->CalcPhysShape(Tr, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);
      real_t dphidt = sh_phi*(*elrate[2]);

      el[2]->CalcPhysDShape(Tr, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m -= f;             // Subtract force
      res_m *= rho;
      res_m += grad_p;        // Add pressure

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int ii = 0; ii < dim; ++ii)
         {
            for (int jj = 0; jj < dim; ++jj)
            {
               res_m[jj] -= mu*(hess_u(jj,hmap(ii,ii)) +
                                hess_u(ii,hmap(jj,ii)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }
      real_t res_c = grad_u.Trace();

      // Compute stability params

      // std::cout<<"rho = "<<rho<<std::endl;
      // std::cout<<"mu = "<<mu<<std::endl;
      // std::cout<<"u = ";  u.Print();
      // std::cout<<"Gij = ";  Gij.Print();
      GetTau(tau_m, tau_c, tau_ls, cfl2, rho, mu, u, Gij);
      elem_cfl = fmax(elem_cfl, cfl2);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      // std::cout<<"tau_m = "<<tau_m<<std::endl;
      // std::cout<<"res_m = "; res_m.Print();
      // std::cout<<"up = "; up.Print();
      //  u += up;
      p -= tau_c*res_c;

      mu += GetKdc(res_m, grad_u, Gij);

      // Compute momentum weak residual
      grad_u.Mult(u,res_m);                       // Add convection (incl cross)
      res_m += dudt;                              // Add acceleration
      res_m -= f;                                 // Add force
      AddMult_a_VWt(w*rho, sh_u, res_m, elv_u);   // Add force + acc term to rhs

      flux.Diag(-p, dim);                         // Add pressure
      grad_u.Symmetrize();                        // Grad to strain
      flux.Add(2*mu,grad_u);                      // Add stress to flux
      AddMult_a_VWt(-rho, up, u, flux);           // Add SUPG to flux (incl rey)
      AddMult_a_ABt(w, shg_u, flux, elv_u);       // Add flux term to rhs

      // Compute continuity weak residual
      elvec[1]->Add(-w*res_c, sh_p);              // Add Galerkin term
      shg_p.Mult(up, sh_p);                       // PSPG help term
      elvec[1]->Add(w, sh_p);                     // Add PSPG term

      // Compute convection weak residual
      real_t res_ls = dphidt + u*grad_phi;
      elvec[2]->Add(w*res_ls, sh_phi);            // Add Galerkin term
      elvec[2]->Add(w*res_ls*tau_ls, ushg_phi);   // Add SUPG term
   }
   //elvec[0]->Print(std::cout,999);
   //elvec[1]->Print(std::cout,999);
   //elvec[2]->Print(std::cout,999);
   elem_cfl = sqrt(elem_cfl);
}

// Assemble the element interior gradient matrices
void RBVMSIntegrator::AssembleElementGrad(
   const Array<const FiniteElement*> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();
   int dof_phi = el[2]->GetDof();

   bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_wphi = *elmats(0,2);

   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);
   DenseMatrix &mat_qphi = *elmats(1,2);

   DenseMatrix &mat_vu = *elmats(2,0);
   DenseMatrix &mat_vp = *elmats(2,1);
   DenseMatrix &mat_vphi = *elmats(2,2);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_wphi.SetSize(dof_u*dim, dof_phi);

   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);
   mat_qphi.SetSize(dof_p, dof_phi);

   mat_vu.SetSize(dof_p, dof_u*dim);
   mat_vp.SetSize(dof_p, dof_p);
   mat_vphi.SetSize(dof_p, dof_phi);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_wphi = 0.0;

   mat_qu = 0.0;
   mat_qp = 0.0;
   mat_qphi = 0.0;

   mat_vu = 0.0;
   mat_vp = 0.0;
   mat_vphi = 0.0;

   DenseMatrix mat_wu1(dof_u, dof_u);
   mat_wu1 = 0.0;
   DenseMatrix mat_wp1[dim], mat_qu1[dim], mat_up[dim], mat_wphi1[dim];
   DenseMatrix mat_vu1[dim];

   for (int i_dim = 0; i_dim < dim; ++i_dim)
   {
      mat_wp1[i_dim].SetSize(dof_u,dof_p);
      mat_qu1[i_dim].SetSize(dof_p,dof_u);
      mat_up[i_dim].SetSize(dof_u,dof_p);
      mat_wphi1[i_dim].SetSize(dof_u,dof_phi);
      mat_vu1[i_dim].SetSize(dof_phi,dof_u);
      mat_wp1[i_dim] = 0.0;
      mat_qu1[i_dim] = 0.0;
      mat_up[i_dim] = 0.0;
      mat_wphi1[i_dim] = 0.0;
      mat_vu1[i_dim] = 0.0;
   }

   int dim2 = (dim*(dim+1))/2;
   DenseMatrix  mat_uu[dim*dim];
   for (int i_dim = 0; i_dim < dim2; ++i_dim)
   {
      mat_uu[i_dim].SetSize(dof_u,dof_u);
      mat_uu[i_dim] = 0.0;
   }

   Vector vec_u1(NULL,dof_u), vec_u2(NULL,dof_u), vec_p(NULL,dof_p);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_u);

   DenseMatrix shg_uT(dim, dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, tau_ls, cfl2;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);

      real_t mu = c_mu.Eval(Tr, ip);

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      // real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      el[2]->CalcPhysShape(Tr, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);
      real_t dphidt = sh_phi*(*elrate[2]);

      el[2]->CalcPhysDShape(Tr, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);
      real_t drho = GetRhoGrad(phi, grad_phi, Gij);

      // Compute strong residual
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m -= f;             // Subtract force
      res_m *= rho;
      res_m += grad_p;        // Add pressure

      /////
      Vector bla(dim);
      grad_u.Mult(u,bla);   // Add convection
      bla += dudt;          // Add acceleration
      bla -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int ii = 0; ii < dim; ++ii)
         {
            for (int jj = 0; jj < dim; ++jj)
            {
               res_m[jj] -= mu*(hess_u(jj,hmap(ii,ii)) +
                                hess_u(ii,hmap(jj,ii)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }

      // Compute stability params
      GetTau(tau_m, tau_c, tau_ls, cfl2, rho, mu, u, Gij);

      mu += GetKdc(res_m, grad_u, Gij);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      // u += up;

      // Compute small scale jacobian
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         dupdu(j_u) = -tau_m*rho*(sh_u(j_u) + ushg_u(j_u)*dt);
      }
      //shg_u.Mult(u, ushg_u);

      // Recompute convective gradient
      MultAtB(elf_u, shg_u, grad_u);
      shg_uT.Transpose(shg_u);

      // Level - set - Block diagonal level-set block (v,phi)

      //   real_t res_ls = dphidt + u*grad_phi;
      //   elvec[2]->Add(w*res_ls, sh_phi);            // Add Galerkin term
      //   elvec[2]->Add(w*res_ls*tau_ls, ushg_phi);   // Add SUPG term

      // Galerkin terms
      AddMult_a_VVt(w, sh_phi, mat_vphi);
      AddMult_a_VWt(dt*w, sh_phi, ushg_phi, mat_vphi);

      // SUPG terms
      AddMult_a_VWt(w*tau_ls, ushg_phi, sh_phi, mat_vphi);
      AddMult_a_VVt(dt*w*tau_ls, ushg_phi, mat_vphi);

      real_t res_ls = dphidt + u*grad_phi;

      Vector col(dof_phi);
      for (int i_dim = 0; i_dim < dim; ++i_dim)
      {
         real_t tmp = dt*w*grad_phi(i_dim);
         // AddMult_a_VWt(tmp,        sh_phi,   sh_u, mat_vu1[i_dim]);
         //AddMult_a_VWt(tmp*tau_ls, ushg_phi, sh_u, mat_vu1[i_dim]);
         AddMult_a_VWt(tmp,        sh_u,sh_phi,    mat_vu1[i_dim]);
         AddMult_a_VWt(tmp*tau_ls, sh_u,ushg_phi,  mat_vu1[i_dim]);

         shg_phi.GetColumn(i_dim,col);
         //   AddMult_a_VWt(dt*w*tau_ls*res_ls, col, sh_u, mat_vu1[i_dim]);
         AddMult_a_VWt(dt*w*tau_ls*res_ls, sh_u, col,  mat_vu1[i_dim]);
      }

      // Cross term Momentem Galerkin wrt phi
      for (int i_dim = 0; i_dim < dim; ++i_dim)
      {
         AddMult_a_VVt(drho*dt*w*bla[i_dim], sh_u, mat_wphi1[i_dim]);
      }

      // Momentum - Block diagonal Velocity block (w,u)
      // Acceleration term
      AddMult_a_VVt(rho*w, sh_u, mat_wu1);

      // Convection terms
      AddMult_a_VWt(dt*rho*w, sh_u, ushg_u, mat_wu1);
      AddMult_a_VWt(-rho*w, ushg_u, dupdu, mat_wu1);

      // Diffusion term
      AddMult_a_AAt(w*mu*dt, shg_u, mat_wu1);

      // Continuity - Pressure block (q,p)
      AddMult_a_AAt(-w*tau_m, shg_p, mat_qp);

      // Off diagional block terms
      int ii = 0;
      for (int i_dim = 0; i_dim < dim; ++i_dim)
      {
         // Getting columns for outer product
         shg_p.GetColumnReference(i_dim, vec_p);
         shg_u.GetColumnReference(i_dim, vec_u1);

         // Momentum + Continuity Galerkin terms
         AddMult_a_VWt(-w, vec_u1, sh_p, mat_up[i_dim]);

         // Momentum - Pressure block (w,p)
         AddMult_a_VWt(tau_m * w, ushg_u, vec_p, mat_wp1[i_dim] );

         // Continuity - Velocity block (q,u)
         AddMult_a_VWt(w, vec_p, dupdu, mat_qu1[i_dim]);

         // Momentum - Velocity block (w,u)
         for (int j_dim = 0; j_dim < i_dim; ++j_dim)
         {
            shg_u.GetColumnReference(j_dim, vec_u2);
            AddMult_a_VWt(w*dt*(mu+tau_c), vec_u1, vec_u2, mat_uu[ii++]);
         }
         AddMult_a_VVt(w*dt*(mu+tau_c), vec_u1, mat_uu[ii++]);
      }
   }

   // Adding sub elements
   int ii = 0;
   for (int i_dim = 0; i_dim < dim; ++i_dim)
   {
      mat_wp1[i_dim] += mat_up[i_dim];
      mat_up[i_dim].Transpose();
      mat_up[i_dim] *= dt;
      mat_qu1[i_dim] += mat_up[i_dim];

      mat_wp.AddSubMatrix(i_dim * dof_u, 0, mat_wp1[i_dim]);
      mat_qu.AddSubMatrix(0, i_dim * dof_u, mat_qu1[i_dim]);
      mat_wu.AddSubMatrix(i_dim*dof_u, mat_wu1);
      for (int j_dim = 0; j_dim < i_dim; ++j_dim)
      {
         mat_wu.AddSubMatrix(i_dim * dof_u, j_dim * dof_u, mat_uu[ii]);
         mat_uu[ii].Transpose();
         mat_wu.AddSubMatrix(j_dim * dof_u, i_dim * dof_u, mat_uu[ii++]);
      }
      mat_wu.AddSubMatrix(i_dim * dof_u, i_dim * dof_u, mat_uu[ii++]);

      mat_wphi.AddSubMatrix(i_dim * dof_u, 0, mat_wphi1[i_dim]);

      mat_vu.AddSubMatrix(0, i_dim * dof_u, mat_vu1[i_dim]);
   }

   mat_wp *= dt;
   mat_qp *= dt;
   mat_vp *= dt;
}


OutflowIntegrator::OutflowIntegrator(Coefficient &rho,
                                     Coefficient &mu,
                                     VectorCoefficient &sol_u)
   : c_rho(rho), c_mu(mu), c_sol_u(sol_u)
{
   dim = c_sol_u.GetVDim();
   u.SetSize(dim);
   dudt.SetSize(dim);
   f.SetSize(dim);
   res_m.SetSize(dim);
   up.SetSize(dim);
   traction.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
   grad_phi.SetSize(dim);
   Gij.SetSize(dim);
   nor.SetSize(dim);
   hn.SetSize(dim);
   hmap.SetSize(dim,dim);

   if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) =  hmap(1,0) =  1;
      hmap(1,1) = 2;
   }
   else if (dim == 3)
   {
      hmap(0,0) = 0;
      hmap(0,1) = hmap(1,0) = 1;
      hmap(0,2) = hmap(2,0) = 2;
      hmap(1,1) = 3;
      hmap(1,2) = hmap(2,1) = 4;
      hmap(2,2) = 5;
   }
   else
   {
      mfem_error("Only implemented for 2D and 3D");
   }
}

// Compute density
real_t OutflowIntegrator::GetRho(real_t &phi, Vector &grad_phi_,
                                 DenseMatrix &Gij_)
{
   real_t He = Heaviside::step(phi, grad_phi_, Gij_);
   return rho0 + (rho1-rho0)*He;
}

// Compute density gradient
real_t OutflowIntegrator::GetRhoGrad(real_t &phi, Vector &grad_phi_,
                                     DenseMatrix &Gij_)
{
   real_t de = Heaviside::dirac(phi, grad_phi_, Gij_);
   return  (rho1-rho0)*de;
}

// Assemble the outflow boundary residual vectors
void OutflowIntegrator
::AssembleFaceVector(const Array<const FiniteElement *> &el1,
                     const Array<const FiniteElement *> &el2,
                     FaceElementTransformations &Tr,
                     const Array<const Vector *> &elsol,
                     const Array<const Vector *> &elrate,
                     const Array<Vector *> &elvec)
{
   //  real_t outflow; // TBD
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);
   elvec[2]->SetSize(dof_phi);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;
   // outflow = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   //Vector
   traction.SetSize(dim);
   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor); // nor = n.da
      real_t w = ip.weight;          // No weight --> taken care of by nor
      // real_t hbc = c_hbc.Eval(Tr, ip);
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      real_t un = u*nor;

      if (un > 0.0) { continue; }

      c_sol_u.Eval(up, *Tr.Elem1, eip);
      up -= u;
      // up.Neg();


      AddMult_a_VWt(rho*w*un, sh_u, up, elv_u);
      //elvec[2]->Add(w*un*(phi0-phi), sh_p);
   }
}

// Assemble the outflow boundary gradient matrices
void OutflowIntegrator
::AssembleFaceGrad(const Array<const FiniteElement *>&el1,
                   const Array<const FiniteElement *>&el2,
                   FaceElementTransformations &Tr,
                   const Array<const Vector *> &elsol,
                   const Array<const Vector *> &elrate,
                   const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();
   //  int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   elmats(0,1)->SetSize(0,0);
   elmats(0,2)->SetSize(0,0);

   elmats(1,0)->SetSize(0,0);
   elmats(1,1)->SetSize(0,0);
   elmats(1,2)->SetSize(0,0);

   elmats(2,0)->SetSize(0,0);
   elmats(2,1)->SetSize(0,0);
   elmats(2,2)->SetSize(0,0);

   mat_wu = 0.0;

   sh_u.SetSize(dof_u);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder()+2;
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      //  const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor); // nor = n.da
      real_t w = ip.weight;          // No weight --> taken care of by nor
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      real_t un = u*nor;

      if (un > 0.0) { continue; }

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      // Momentum - Velocity block (w,u)
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp = -rho*sh_u(j_u)*un*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += sh_u(i_u)*tmp;
            }
         }
      }

      // Jacobian due to the normal velocity
      /*   for (int dim_v = 0; dim_v < dim; ++dim_v)
         {
            real_t tmp0 = rho*nor(dim_v)*w*dt;
            for (int j_u = 0; j_u < dof_u; ++j_u)
            {
               real_t tmp1 = tmp0*sh_u(j_u);
               for (int dim_u = 0; dim_u < dim; ++dim_u)
               {
                  real_t tmp2 = tmp1*u(dim_u);
                  for (int i_u = 0; i_u < dof_u; ++i_u)
                  {
                     mat_wu(i_u + dim_u*dof_u, j_u + dim_v*dof_u) += sh_u(i_u)*tmp2;
                  }
               }
            }
         }*/
   }
}


WeakBCIntegrator::WeakBCIntegrator(Coefficient &rho,
                                   Coefficient &mu,
                                   VectorCoefficient &sol_u,
                                   Coefficient &sol_phi)
   : c_rho(rho), c_mu(mu), c_sol_u(sol_u), c_sol_phi(sol_phi)
{
   dim = c_sol_u.GetVDim();
   u.SetSize(dim);
   dudt.SetSize(dim);
   f.SetSize(dim);
   res_m.SetSize(dim);
   up.SetSize(dim);
   traction.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
   grad_phi.SetSize(dim);
   Gij.SetSize(dim);
   nor.SetSize(dim);
   hn.SetSize(dim);
   hmap.SetSize(dim,dim);

   if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) =  hmap(1,0) =  1;
      hmap(1,1) = 2;
   }
   else if (dim == 3)
   {
      hmap(0,0) = 0;
      hmap(0,1) = hmap(1,0) = 1;
      hmap(0,2) = hmap(2,0) = 2;
      hmap(1,1) = 3;
      hmap(1,2) = hmap(2,1) = 4;
      hmap(2,2) = 5;
   }
   else
   {
      mfem_error("Only implemented for 2D and 3D");
   }
}

// Compute density
real_t WeakBCIntegrator::GetRho(real_t &phi, Vector &grad_phi_,
                                DenseMatrix &Gij_)
{
   real_t He = Heaviside::step(phi, grad_phi_, Gij_);
   return rho0 + (rho1-rho0)*He;
}

// Compute density gradient
real_t WeakBCIntegrator::GetRhoGrad(real_t &phi, Vector &grad_phi_,
                                    DenseMatrix &Gij_)
{
   real_t de = Heaviside::dirac(phi, grad_phi_, Gij_);
   return  (rho1-rho0)*de;
}

// Compute Weak Dirichlet stabilisation parameters
void WeakBCIntegrator::GetTauB(real_t &tau_b, real_t &tau_n,
                               real_t &mu, Vector &u_,
                               Vector &nor_,
                               DenseMatrix &Gij_)
{
   real_t ho1 = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t nj = nor_[j];
      for (int i = 0; i < dim; i++)
      {
         ho1 += Gij_(i,j)*nor_[i]*nj;
      }
   }

   tau_b = Cb*mu*sqrt(ho1);
   tau_n = Cn;
}

// Assemble the weak Dirichlet BC boundary residual vectors
void WeakBCIntegrator
::AssembleFaceVector(const Array<const FiniteElement *> &el1,
                     const Array<const FiniteElement *> &el2,
                     FaceElementTransformations &Tr,
                     const Array<const Vector *> &elsol,
                     const Array<const Vector *> &elrate,
                     const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);
   elvec[2]->SetSize(dof_phi);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);

   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);
   ushg_phi.SetSize(dof_phi);

   flux.SetSize(dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n, mu, w, phi,rho, un;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      mu = c_mu.Eval(*Tr.Elem1, eip);
      c_sol_u.Eval(up, *Tr.Elem1, eip);
      real_t phi0 = c_sol_phi.Eval(*Tr.Elem1, eip);

      w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      up -= u;
      up.Neg();
      un = up*nor;

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Symmetrize();           // Grad to strain

      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      phi = sh_phi*(*elsol[2]);
      //dphidt = sh_phi*(*elrate[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      rho = GetRho(phi, grad_phi, Gij);

      GetTauB(tau_b, tau_n, mu, u, nor, Gij);

      // Traction
      grad_u.Mult(nor, traction);
      traction *= -2*mu;             // Consistency
      traction.Add(p, nor);         // Pressure
      traction.Add(tau_b,up);        // Penalty
      traction.Add(tau_n*un,nor);   // Penalty -- normal
      AddMult_a_VWt(w, sh_u, traction, elv_u);

      // Dual consistency
      MultVWt(nor, up, flux);
      flux.Symmetrize();
      AddMult_a_ABt(-w*2*mu, shg_u, flux, elv_u);

      // Continuity
      elvec[1]->Add(w*un, sh_p);

      // Convection
      un = u*nor;
      if (un > 0.0) { continue; }

      AddMult_a_VWt(-rho*w*un, sh_u, up, elv_u);

      elvec[2]->Add(w*un*(phi0-phi), sh_p);
   }
   //elvec[0]->Print();
   // elvec[1]->Print();
   // elvec[2]->Print();
}

// Assemble the weak Dirichlet BC boundary gradient matrices
void WeakBCIntegrator
::AssembleFaceGrad(const
                   Array<const FiniteElement *>&el1,
                   const Array<const FiniteElement *>&el2,
                   FaceElementTransformations &Tr,
                   const Array<const Vector *> &elsol,
                   const Array<const Vector *> &elrate,
                   const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_wphi = *elmats(0,2);

   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);
   DenseMatrix &mat_qphi = *elmats(1,2);

   DenseMatrix &mat_vu = *elmats(2,0);
   DenseMatrix &mat_vp = *elmats(2,1);
   DenseMatrix &mat_vphi = *elmats(2,2);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_wphi.SetSize(dof_u*dim, dof_phi);

   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);
   mat_qphi.SetSize(dof_p, dof_phi);

   mat_vu.SetSize(dof_p, dof_u*dim);
   mat_vp.SetSize(dof_p, dof_p);
   mat_vphi.SetSize(dof_p, dof_phi);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_wphi = 0.0;

   mat_qu = 0.0;
   mat_qp = 0.0;
   mat_qphi = 0.0;

   mat_vu = 0.0;
   mat_vp = 0.0;
   mat_vphi = 0.0;

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n, mu, w, phi,rho;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();
      mu = c_mu.Eval(*Tr.Elem1, eip);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      phi = sh_phi*(*elsol[2]);
      //dphidt = sh_phi*(*elrate[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);
      shg_phi.Mult(u, ushg_phi);

      rho = GetRho(phi, grad_phi, Gij);

      GetTauB(tau_b, tau_n, mu, u, nor, Gij);

      // Momentum - Velocity block (w,u)
      // Consistency & Dual: Part 1
      real_t tmp0 = mu*w*dt;
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            real_t tmp1 = nor(j_dim)*shg_u(j_u,j_dim)*tmp0;
            real_t tmp2 = nor(j_dim)*sh_u(j_u)*tmp0;
            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               int j_dof = j_u + i_dim*dof_u;
               for (int i_u = 0; i_u < dof_u; ++i_u)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  -= (sh_u(i_u)*tmp1 + shg_u(i_u,j_dim)*tmp2);
               }
            }
         }
      }

      // Consistency & Dual: Part 2
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            int j_dof = j_u + j_dim*dof_u;
            real_t tmp1 = nor(j_dim)*sh_u(j_u)*tmp0;
            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               real_t tmp2 = nor(j_dim)*shg_u(j_u,i_dim)*tmp0;
               for (int i_u = 0; i_u < dof_u; ++i_u)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  -= (sh_u(i_u)*tmp2 +shg_u(i_u,i_dim)*tmp1);
               }
            }
         }
      }

      // Penalty
      tmp0 = tau_b*w*dt;
      for (int dim_u = 0; dim_u < dim; ++dim_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            real_t tmp1 = sh_u(j_u)*tmp0;
            int j_dof = j_u + dim_u*dof_u;
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_dof) += sh_u(i_u)*tmp1;
            }
         }
      }

      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         tmp0 = sh_u(j_u)*w*dt*tau_n;
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            real_t tmp1 = tmp0*nor(j_dim);
            int j_dof = j_u + j_dim*dof_u;
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               for (int i_dim = 0; i_dim < dim; ++i_dim)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  += sh_u(i_u)*nor(i_dim)*tmp1;
               }
            }
         }
      }

      // Momentum - Pressure block (w,p)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         tmp0 = sh_p(i_p)*w;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            real_t tmp1 = nor(dim_u)*tmp0;
            for (int j_u = 0; j_u < dof_u; ++j_u)
            {
               mat_wp(j_u + dof_u * dim_u, i_p) += sh_u(j_u)*tmp1;
            }
         }
      }

      // Continuity - Velocity block (q,u)
      for (int dim_u = 0; dim_u < dim; ++dim_u)
      {
         tmp0 = nor(dim_u)*w*dt;
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            real_t tmp1 = sh_u(j_u)*tmp0;
            for (int i_p = 0; i_p < dof_p; ++i_p)
            {
               mat_qu(i_p, j_u + dof_u * dim_u) += sh_p(i_p)*tmp1;
            }
         }
      }

      // Convection: Momentum - Velocity block (w,u)
      real_t un = u*nor;
      if (un > 0.0) { continue; }
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp = -rho*sh_u(j_u)*un*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += sh_u(i_u)*tmp;
            }
         }
      }

      //
      for (int jj = 0; jj < dof_phi; ++jj)
      {
         // real_t tmp = -un*w*dt;
         for (int ii = 0; ii < dof_phi; ++ii)
         {
            mat_vphi(ii, jj) -= sh_phi(ii)*sh_phi(jj)*un*w*dt;
         }
      }

   }

   mat_wp *= dt;
   mat_qp *= dt;
   mat_vp *= dt;
}

/*

// Assemble the normal Dirichlet BC boundary residual vectors
void IncNavStoIntegrator
::AssembleNormalBCVector(const Array<const FiniteElement *> &el1,
                         const Array<const FiniteElement *> &el2,
                         FaceElementTransformations &Tr,
                         const Array<const Vector *> &elsol,
                         const Array<const Vector *> &elrate,
                         const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   *elvec[2] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   sh_u.SetSize(dof_u);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_n = 10000.0; // TO BE MOVED
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      c_sol_u.Eval(up, *Tr.Elem1, eip);

      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      up -= u;
      up.Neg();
      real_t un = up*nor;
      AddMult_a_VWt(w*rho*tau_n*un, sh_u, nor, elv_u);
   }
}

// Assemble the normal Dirichlet BC boundary gradient matrices
void IncNavStoIntegrator
::AssembleNormalBCGrad(const Array<const FiniteElement *>&el1,
                       const Array<const FiniteElement *>&el2,
                       FaceElementTransformations &Tr,
                       const Array<const Vector *> &elsol,
                       const Array<const Vector *> &elrate,
                       const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();
   int dof_phi = el1[2]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_wphi = *elmats(0,2);

   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);
   DenseMatrix &mat_qphi = *elmats(1,2);

   DenseMatrix &mat_vu = *elmats(2,0);
   DenseMatrix &mat_vp = *elmats(2,1);
   DenseMatrix &mat_vphi = *elmats(2,2);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_wphi.SetSize(dof_u*dim, dof_phi);

   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);
   mat_qphi.SetSize(dof_p, dof_phi);

   mat_vu.SetSize(dof_p, dof_u*dim);
   mat_vp.SetSize(dof_p, dof_p);
   mat_vphi.SetSize(dof_p, dof_phi);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_wphi = 0.0;

   mat_qu = 0.0;
   mat_qp = 0.0;
   mat_qphi = 0.0;

   mat_vu = 0.0;
   mat_vp = 0.0;
   mat_vphi = 0.0;

   sh_u.SetSize(dof_u);
   sh_phi.SetSize(dof_phi);
   shg_phi.SetSize(dof_phi, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_n = 10000.0; // TO BE MOVED
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      real_t w = ip.weight * Tr.Weight();
      MultAtB(Tr.Elem1->InverseJacobian(),Tr.Elem1->InverseJacobian(),Gij);
      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);

      el1[2]->CalcPhysShape(*Tr.Elem1, sh_phi);
      real_t phi = sh_phi*(*elsol[2]);

      el1[2]->CalcPhysDShape(*Tr.Elem1, shg_phi);
      shg_phi.MultTranspose(*elsol[2], grad_phi);

      real_t rho = GetRho(phi, grad_phi, Gij);

      // Penalty
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp0 = sh_u(j_u)*w*dt*rho*tau_n;
         for (int j_dim = 0; j_dim < dim; ++j_dim)
         {
            real_t tmp1 = tmp0*nor(j_dim);
            int j_dof = j_u + j_dim*dof_u;
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               for (int i_dim = 0; i_dim < dim; ++i_dim)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  += sh_u(i_u)*nor(i_dim)*tmp1;
               }
            }
         }
      }
   }
}*/


