// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "integrator.hpp"

using namespace mfem;

// Constructor
IncNavStoIntegrator::IncNavStoIntegrator(Coefficient &rho,
                                         Coefficient &mu,
                                         VectorCoefficient &force,
                                         VectorCoefficient &sol,
                                         Coefficient &suction,
                                         Coefficient &blowing)
   : c_rho(rho), c_mu(mu), c_force(force), c_sol(sol),
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
void IncNavStoIntegrator::GetTau(real_t &tau_m, real_t &tau_c, real_t &cfl2,
                                 real_t& rho, real_t &mu, Vector &u,
                                 ElementTransformation &T)
{
   real_t Cd = 32.0;
   real_t Ct = 1.0;

   // Metric tensor
   MultAtB(T.InverseJacobian(),T.InverseJacobian(),Gij);
   //   real_t tr =  Gij.Trace();
   //   Gij.Diag(tr/real_t(dim), dim);

   // Temporal part
   tau_m = Ct/(dt*dt);

   // Convective part
   tau_c = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t uj = u[j];
      for (int i = 0; i < dim; i++)
      {
         tau_c += Gij(i,j)*u[i]*uj;
      }
   }
   cfl2 = tau_c/tau_m;
   tau_m += tau_c;
   tau_m *= rho*rho;

   // Diffusive part
   real_t tmp = Cd*Cd*mu*mu;
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_m  += tmp*Gij(i,j)*Gij(i,j);
      }
   }

   // Momentum stabilisation parameter
   tau_m = 1.0/sqrt(tau_m);

   // Continuity stabilisation parameter
   tau_c = 1.0/(tau_m*Gij.Trace());
}

// Compute element average artifical diffusion
real_t IncNavStoIntegrator::GetElemArtDiff(const Array<const FiniteElement *>
                                           &el,
                                           ElementTransformation &Tr,
                                           const Array<const Vector *> &elsol,
                                           const Array<const Vector *> &elrate)
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   //bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   // shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   shg_p.SetSize(dof_p, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   double num = 0.0;
   double den = 0.0;
   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      c_force.Eval(f, Tr, ip);

      // Velocity gradient in reference element
      el[0]->CalcDShape(ip, shg_u);
      shg_u.Mult(u, ushg_u);
      double ndux = grad_u.FNorm(); // Norm of gradient in reference space

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      // Metric tensor
      MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);
      //      real_t tr =  Gij.Trace();
      //      Gij.Diag(tr/real_t(dim), dim);

      // Compute numerator and denominator
      double nr = res_m.Norml2();
      num += nr*w;
      den += ndux*w;
   }
   return 0.0*num/fmax(den,1e-6);
}

// Compute Weak Dirichlet stabilisation parameters
void IncNavStoIntegrator::GetTauB(real_t &tau_b, real_t &tau_n,
                                  real_t &mu, Vector &u,
                                  Vector &nor,
                                  FaceElementTransformations &Tr)
{
   real_t Cb = 12.0;

   Tr.Elem1->InverseJacobian().Mult(nor,hn);
   tau_b = Cb*mu*hn.Norml2();
   tau_n = 100.0;
}

// Get energy
real_t IncNavStoIntegrator::GetElementEnergy(
   const Array<const FiniteElement *>&el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate)
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::GetElementEnergy"
                 " has incorrect block finite element space size!");
   }
   int dof_u = el[0]->GetDof();
   Vector tau(3);

   sh_u.SetSize(dof_u);
   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   real_t energy = 0.0;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);

      real_t w = ip.weight * Tr.Weight();

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);

      energy += w*(u*u)/2;
   }

   return energy;
}

// Assemble the element interior residual vectors
void IncNavStoIntegrator::AssembleElementVector(
   const Array<const FiniteElement *> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array<Vector *> &elvec)
//  real_t &elem_cfl)
{
   real_t elem_cfl;
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }
   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, cfl2;
   elem_cfl = 0.0;

   double mu_ad = 0.0;// GetElemArtDiff(el, Tr, elsol, elrate);

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      real_t rho = c_rho.Eval(Tr, ip);
      real_t mu = c_mu.Eval(Tr, ip);
      real_t mu_eff = mu + mu_ad;
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

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m *= rho;
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int i = 0; i < dim; ++i)
         {
            for (int j = 0; j < dim; ++j)
            {
               res_m[j] -= mu*(hess_u(j,hmap(i,i)) +
                               hess_u(i,hmap(j,i)));
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
      GetTau(tau_m, tau_c, cfl2, rho, mu, u, Tr);
      elem_cfl = fmax(elem_cfl, cfl2);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      u += up;
      p -= tau_c*res_c;

      // Compute momentum weak residual
      flux.Diag(-p, dim);                         // Add pressure
      grad_u.Symmetrize();                        // Grad to strain
      flux.Add(2*mu_eff,grad_u);                  // Add stress to flux
      AddMult_a_VVt(-rho, u, flux);               // Add convection to flux
      AddMult_a_ABt(w, shg_u, flux, elv_u);       // Add flux term to rhs
      f.Add(-rho, dudt);                          // Add Acceleration to force
      AddMult_a_VWt(-w, sh_u, f, elv_u);          // Add force + acc term to rhs

      // Compute continuity weak residual
      elvec[1]->Add(-w*res_c, sh_p);              // Add Galerkin term
      shg_p.Mult(up, sh_p);                       // PSPG help term
      elvec[1]->Add(w, sh_p);                     // Add PSPG term
   }
   //elvec[0]->Print(std::cout,888);
   //elvec[1]->Print(std::cout,888);
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

   bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_qu = 0.0;
   mat_qp = 0.0;

   DenseMatrix mat_wu1(dof_u, dof_u);
   mat_wu1 = 0.0;
   DenseMatrix mat_wp1[dim], mat_qu1[dim], mat_up[dim];
   for (int i_dim = 0; i_dim < dim; ++i_dim)
   {
      mat_wp1[i_dim].SetSize(dof_u,dof_p);
      mat_qu1[i_dim].SetSize(dof_p,dof_u);
      mat_up[i_dim].SetSize(dof_u,dof_p);
      mat_wp1[i_dim] = 0.0;
      mat_qu1[i_dim] = 0.0;
      mat_up[i_dim] = 0.0;
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

   DenseMatrix shg_uT(dim, dof_u);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, cfl2;

   //double mu_ad = 0.0;//GetElemArtDiff(el, Tr, elsol, elrate);

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();

      real_t rho = c_rho.Eval(Tr, ip);
      real_t mu = c_mu.Eval(Tr, ip);
      //real_t mu_eff = mu + mu_ad;

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);

      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m *= rho;
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int i = 0; i < dim; ++i)
         {
            for (int j = 0; j < dim; ++j)
            {
               res_m[j] -= mu*(hess_u(j,hmap(i,i)) +
                               hess_u(i,hmap(j,i)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }

      // Compute stability params
      GetTau(tau_m, tau_c, cfl2, rho, mu, u, Tr);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      u += up;

      // Compute small scale jacobian
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         dupdu(j_u) = -tau_m*rho*(sh_u(j_u) + ushg_u(j_u)*dt);
      }
      shg_u.Mult(u, ushg_u);

      // Recompute convective gradient
      MultAtB(elf_u, shg_u, grad_u);

      shg_uT.Transpose(shg_u);

      // Momentum - Block diagonal Velocity block (w,u)
      // Acceleration term
      AddMult_a_VVt(rho*w, sh_u, mat_wu1);

      // Convection terms
      AddMult_a_VWt(-dt*rho*w, ushg_u, sh_u, mat_wu1);
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
   }

   /*std::cout<<"-------------------------\n";
      std::cout<<mat_wu.FNorm()<<std::endl;
      std::cout<<mat_wp.FNorm()<<std::endl;
      std::cout<<mat_qu.FNorm()<<std::endl;
      std::cout<<mat_qp.FNorm()<<std::endl;
   std::cout<<"-------------------------\n";*/
   mat_wp *= dt;
   mat_qp *= dt;
}

// Assemble the outflow boundary residual vectors
void IncNavStoIntegrator
::AssembleFaceVector(const Array<const FiniteElement *> &el1,
                     const Array<const FiniteElement *> &el2,
                     FaceElementTransformations &Tr,
                     const Array<const Vector *> &elsol,
                     const Array<const Vector *> &elrate,
                     const Array<Vector *> &elvec)
{

   //std::cout<<"IncNavStoIntegrator::AssembleFaceVector"<<std::endl;
   real_t outflow = 0.0;
   bool suction = false;

   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;
   outflow = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   Vector traction(dim);
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

      real_t rho = c_rho.Eval(*Tr.Elem1, eip);
      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      real_t un = u*nor;
      if (suction)
      {
         real_t suction_val = c_suction.Eval(*Tr.Elem1, eip);
         traction.Set(suction_val,nor);
         AddMult_a_VWt(w, sh_u, traction, elv_u);
         outflow += rho * un * w; // No weight --> taken care of by nor
      }

      if (un < 0.0) { continue; }

      AddMult_a_VWt(rho*w*un, sh_u, u, elv_u);
   }
}

// Assemble the outflow boundary gradient matrices
void IncNavStoIntegrator
::AssembleFaceGrad(const Array<const FiniteElement *>&el1,
                   const Array<const FiniteElement *>&el2,
                   FaceElementTransformations &Tr,
                   const Array<const Vector *> &elsol,
                   const Array<const Vector *> &elrate,
                   const Array2D<DenseMatrix *> &elmats)
{
   //std::cout<<"IncNavStoIntegrator::AssembleFaceGrad"<<std::endl;
   real_t outflow = 0.0;
   bool suction = false;

   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(0,0);
   mat_qu.SetSize(0,0);
   mat_qp.SetSize(0,0);

   mat_wu = 0.0;

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder()+2;
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

      real_t rho = c_rho.Eval(*Tr.Elem1, eip);
      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);
      real_t un = u*nor;

      if (un < 0.0) { continue; }

      // Momentum - Velocity block (w,u)
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp = rho*sh_u(j_u)*un*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += sh_u(i_u)*tmp;
            }
         }
      }

      // Jacobian due to the normal velocity
      for (int dim_v = 0; dim_v < dim; ++dim_v)
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
      }
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

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      real_t rho = c_rho.Eval(*Tr.Elem1, eip);
      real_t mu = c_mu.Eval(*Tr.Elem1, eip);
      c_sol.Eval(up, *Tr.Elem1, eip);

      real_t w = ip.weight * Tr.Weight();
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
      real_t un = up*nor;

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Symmetrize();           // Grad to strain

      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);
      real_t p = sh_p*(*elsol[1]);

      GetTauB(tau_b, tau_n, mu, u, nor,Tr);

      // Traction
      grad_u.Mult(nor, traction);
      traction *= -2*mu;             // Consistency
      traction.Add(p, nor);          // Pressure
      traction.Add(tau_b,up);        // Penalty
      traction.Add(tau_n*un,nor);    // Penalty -- normal
      AddMult_a_VWt(w, sh_u, traction, elv_u);

      // Dual consistency
      MultVWt(nor,up, flux);
      flux.Symmetrize();
      AddMult_a_ABt(-w*2*mu, shg_u, flux, elv_u);

      // Continuity
      elvec[1]->Add(w*un, sh_p);

      // Convection
      un = u*nor;
      if (un < 0.0) { continue; }

      AddMult_a_VWt(rho*w*un, sh_u, up, elv_u);
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

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);

   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_qu = 0.0;
   mat_qp = 0.0;

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();
      real_t rho = c_rho.Eval(*Tr.Elem1, eip);
      real_t mu = c_mu.Eval(*Tr.Elem1, eip);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      real_t w = ip.weight * Tr.Weight();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);

      GetTauB(tau_b, tau_n, mu, u, nor,Tr);

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
         real_t tmp0 = sh_u(j_u)*w*dt*tau_n;
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
         real_t tmp0 = sh_p(i_p)*w;
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
         real_t tmp0 = nor(dim_u)*w*dt;
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
      if (un < 0.0) { continue; }
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp = rho*sh_u(j_u)*un*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            for (int i_u = 0; i_u < dof_u; ++i_u)
            {
               mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += sh_u(i_u)*tmp;
            }
         }
      }

   }
}
