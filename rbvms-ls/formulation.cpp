// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "formulation.hpp"

using namespace mfem;
using namespace RBVMS;


// Set the boundaries were strong Dirichlet BCs are imposed
void NavStoLSForm::SetStrongBC (Array<int> strong_bdr)
{
   strong_bdr.Copy(strongBCBdr);

   // Translate indices
   Array<Array<int> *> ess_bdr(3);
   Array<int> ess_bdr_u(fes[0]->GetMesh()->bdr_attributes.Max());
   Array<int> ess_bdr_p(fes[1]->GetMesh()->bdr_attributes.Max());
   Array<int> ess_bdr_phi(fes[2]->GetMesh()->bdr_attributes.Max());

   ess_bdr_u = 0;
   ess_bdr_p = 0;
   ess_bdr_phi = 0;
   for (int b = 0; b < strongBCBdr.Size(); ++b)
   {
      ess_bdr_u[strongBCBdr[b]-1] = 1;
   }

   ess_bdr[0] = &ess_bdr_u;
   ess_bdr[1] = &ess_bdr_p;
   ess_bdr[2] = &ess_bdr_phi;

   // Dummy rhs
   Array<Vector *> rhs(3);
   rhs = nullptr;

   // Enforce BCs using function form parent class
   SetEssentialBC(ess_bdr, rhs);
}

// Set the solution of the previous time step
// and the timestep size of the current solve.
void NavStoLSForm::SetTimeAndSolution(const real_t t,
                                      const real_t dt_,
                                      const Vector &x0)
{
   xs0.SetSize(block_offsets[1]-block_offsets[0]);
   xs2.SetSize(block_offsets[3]-block_offsets[2]);
   xs_true.Update(const_cast<Vector &>(x0), block_trueOffsets);

   fes[0]->GetProlongationMatrix()->Mult(
      xs_true.GetBlock(0), xs0);

   fes[2]->GetProlongationMatrix()->Mult(
      xs_true.GetBlock(2), xs2);

   dt = dt_;
   integrator.SetTimeAndStep(t,dt);
}

// Clear the gradient matrix
void NavStoLSForm::ResetGradient()
{
   hasGrad = false;
   gradCalls = 0;
}

// Block T-Vector to Vector
Vector NavStoLSForm::GetEnergies(const Vector &x) const
{
   xs_true.Update(const_cast<Vector &>(x), block_trueOffsets);
   fes[0]->GetProlongationMatrix()->Mult(xs_true.GetBlock(0), xs.GetBlock(0));
   fes[1]->GetProlongationMatrix()->Mult(xs_true.GetBlock(1), xs.GetBlock(1));
   fes[2]->GetProlongationMatrix()->Mult(xs_true.GetBlock(2), xs.GetBlock(2));

   // Actual assembly
   Array<Array<int> *>vdofs(fes.Size());
   Array<Vector *> el_x(fes.Size());
   Array<const Vector *> el_x_const(fes.Size());

   Array<const FiniteElement *> fe(fes.Size());
   ElementTransformation *T;
   Array<DofTransformation *> doftrans(fes.Size()); doftrans = nullptr;
   Mesh *mesh = fes[0]->GetMesh();

   for (int s=0; s<fes.Size(); ++s)
   {
      el_x_const[s] = el_x[s] = new Vector();
      vdofs[s] = new Array<int>;
   }

   // Domain interior
   Vector energy(3);
   Vector el_energy(3);
   energy.SetSize(3);
   energy = 0.0;
   for (int i = 0; i < fes[0]->GetNE(); ++i)
   {
      T = fes[0]->GetElementTransformation(i);
      for (int s = 0; s < fes.Size(); ++s)
      {
         doftrans[s] = fes[s]->GetElementVDofs(i, *(vdofs[s]));
         fe[s] = fes[s]->GetFE(i);
         xs.GetBlock(s).GetSubVector(*(vdofs[s]), *el_x[s]);
         if (doftrans[s])
         {
            MFEM_WARNING("NavStoLSForm::Doftrans");
            doftrans[s]->InvTransformPrimal(*el_x[s]);
         }
      }

      integrator.AssembleElementEnergy(fe, *T,
                                       el_x_const,
                                       el_energy);
      energy += el_energy;
   }

   // Communicate boundary force
   Vector tmp(energy);
   MPI_Allreduce(tmp.GetData(), energy.GetData(), energy.Size(),
                 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return energy;
}

// Block T-Vector to Block T-Vector
void NavStoLSForm::Mult(const Vector &dx, Vector &y) const
{
   // dxs_true is not modified, so const_cast is okay
   dxs_true.Update(const_cast<Vector &>(dx), block_trueOffsets);
   ys_true.Update(y, block_trueOffsets);
   xs.Update(block_offsets);
   dxs.Update(block_offsets);
   ys.Update(block_offsets);

   fes[0]->GetProlongationMatrix()->Mult(
      dxs_true.GetBlock(0), dxs.GetBlock(0));

   add(xs0,dt,dxs.GetBlock(0),xs.GetBlock(0));   // x = x0 + dt*dx

   fes[1]->GetProlongationMatrix()->Mult(
      dxs_true.GetBlock(1), xs.GetBlock(1));

   fes[2]->GetProlongationMatrix()->Mult(
      dxs_true.GetBlock(2), dxs.GetBlock(2));

   add(xs2,dt,dxs.GetBlock(2),xs.GetBlock(2));   // x = x0 + dt*dx

   // Actual assembly
   MultBlocked(xs, dxs, ys);

   // Finalize assembly
   for (int s=0; s<fes.Size(); ++s)
   {
      fes[s]->GetProlongationMatrix()->MultTranspose(
         ys.GetBlock(s), ys_true.GetBlock(s));
      ys_true.GetBlock(s).SetSubVector(*ess_tdofs[s], 0.0);
   }

   ys_true.SyncFromBlocks();
   y.SyncMemory(ys_true);

   // Communicate CFL
   real_t tmp = cfl;
   MPI_Allreduce(&tmp, &cfl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   // Communicate outflow
   tmp = outflow;
   MPI_Allreduce(&tmp, &outflow, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   // Communicate boundary force
   DenseMatrix tmpDM(bdrForce);
   MPI_Allreduce(tmpDM.GetData(), bdrForce.GetData(),
                 bdrForce.NumRows()*bdrForce.NumCols(),
                 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

// Specialized version of Mult() for BlockVectors
void NavStoLSForm::MultBlocked(const BlockVector &bx,
                               const BlockVector &bdx,
                               BlockVector &by) const
{
   Array<Array<int> *>vdofs(fes.Size());
   Array<Array<int> *>vdofs2(fes.Size());
   Array<Vector *> el_x(fes.Size());
   Array<const Vector *> el_x_const(fes.Size());
   Array<Vector *> el_dx(fes.Size());
   Array<const Vector *> el_dx_const(fes.Size());
   Array<Vector *> el_y(fes.Size());
   Array<const FiniteElement *> fe(fes.Size());
   Array<const FiniteElement *> fe2(fes.Size());
   ElementTransformation *T;
   FaceElementTransformations *Tr;
   Array<DofTransformation *> doftrans(fes.Size()); doftrans = nullptr;
   Mesh *mesh = fes[0]->GetMesh();

   by.UseDevice(true);
   by = 0.0;
   by.SyncToBlocks();
   real_t el_cfl;
   cfl = 0.0;
   for (int s=0; s<fes.Size(); ++s)
   {
      el_x_const[s] = el_x[s] = new Vector();
      el_dx_const[s] = el_dx[s] = new Vector();
      el_y[s] = new Vector();
      vdofs[s] = new Array<int>;
      vdofs2[s] = new Array<int>;
   }

   // Domain interior
   for (int i = 0; i < fes[0]->GetNE(); ++i)
   {
      T = fes[0]->GetElementTransformation(i);
      for (int s = 0; s < fes.Size(); ++s)
      {
         doftrans[s] = fes[s]->GetElementVDofs(i, *(vdofs[s]));
         fe[s] = fes[s]->GetFE(i);
         bx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_x[s]);
         bdx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_dx[s]);
         if (doftrans[s])
         {
            MFEM_WARNING("NavStoLSForm::Doftrans");
            doftrans[s]->InvTransformPrimal(*el_x[s]);
            doftrans[s]->InvTransformPrimal(*el_dx[s]);
         }
      }

      integrator.AssembleElementVector(fe, *T,
                                       el_x_const,
                                       el_dx_const,
                                       el_y,
                                       el_cfl);
      cfl = fmax(cfl, el_cfl);
      for (int s=0; s<fes.Size(); ++s)
      {
         if (el_y[s]->Size() == 0) { continue; }
         if (doftrans[s]) {doftrans[s]->TransformDual(*el_y[s]); }
         by.GetBlock(s).AddElementVector(*(vdofs[s]), *el_y[s]);
      }
   }

   // Domain boundary Outflow
   real_t el_outflow;
   outflow = 0.0;
   for (int i = 0; i < mesh->GetNBE(); ++i)
   {
      // Determine if boundary is outflow
      const int bdr_attr = mesh->GetBdrAttribute(i);
      bool outflowBC = false;
      bool suctionBC = false;
      for (int b=0; b<outflowBdr.Size(); ++b)
      {
         if ( bdr_attr == outflowBdr[b]) { outflowBC = true; }
      }
      for (int b=0; b<suctionBdr.Size(); ++b)
      {
         if ( bdr_attr == suctionBdr[b]) { suctionBC= true; }
      }
      if ( !outflowBC && !suctionBC ) { continue; }

      // Perform assembly over outflow
      Tr = mesh->GetBdrFaceTransformations(i);
      if (Tr != NULL)
      {
         for (int s=0; s<fes.Size(); ++s)
         {
            fe[s] = fes[s]->GetFE(Tr->Elem1No);
            fe2[s] = fes[s]->GetFE(Tr->Elem1No);

            fes[s]->GetElementVDofs(Tr->Elem1No, *(vdofs[s]));
            bx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_x[s]);
            bdx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_dx[s]);
            if (doftrans[s])
            {
               MFEM_WARNING("NavStoLSForm::Doftrans");
               doftrans[s]->InvTransformPrimal(*el_x[s]);
               doftrans[s]->InvTransformPrimal(*el_dx[s]);
            }
         }

         integrator.AssembleOutflowVector(fe, fe2, *Tr,
                                          el_x_const, el_dx_const, el_y,
                                          el_outflow, suctionBC);
         outflow += el_outflow;
         for (int s=0; s<fes.Size(); ++s)
         {
            if (el_y[s]->Size() == 0) { continue; }
            if (doftrans[s]) {doftrans[s]->TransformDual(*el_y[s]); }
            by.GetBlock(s).AddElementVector(*(vdofs[s]), *el_y[s]);
         }
      }
   }

   // Conservative force extraction
   fes[0]->GetProlongationMatrix()->MultTranspose(by.GetBlock(0),
                                                  ys_true.GetBlock(0));

   int nbdr = fes[0]->GetMesh()->bdr_attributes.Max();
   bdrForce.SetSize(nbdr,fes[0]->GetVDim());
   Array<int> dofs, bdr(nbdr);
   Vector vrhs;
   for (int b=0; b<nbdr; ++b)
   {
      bdr = 0; bdr[b] = 1;
      for (int v=0; v<fes[0]->GetVDim(); ++v)
      {
         fes[0]->GetEssentialTrueDofs(bdr, dofs, v);
         ys_true.GetBlock(0).GetSubVector(dofs, vrhs);
         bdrForce(b,v) = vrhs.Sum();
      }
   }

   // Domain boundary weak Dirichelet BC
   for (int i = 0; i < mesh->GetNBE(); ++i)
   {
      // Determine if boundary is outflow
      const int bdr_attr = mesh->GetBdrAttribute(i);
      bool weakBC = false;
      bool blowingBC = false;
      bool normalBC = false;
      for (int b=0; b<weakBCBdr.Size(); ++b)
      {
         if ( bdr_attr == weakBCBdr[b]) { weakBC = true; }
      }
      for (int b=0; b<blowingBdr.Size(); ++b)
      {
         if ( bdr_attr == blowingBdr[b]) { blowingBC = true; }
      }
      for (int b=0; b<normalBCBdr.Size(); ++b)
      {
         if ( bdr_attr == normalBCBdr[b]) { normalBC = true; }
      }
      if ( !weakBC && !blowingBC && !normalBC) { continue; }

      // Perform assembly over Dirichlet boundary
      Tr = mesh->GetBdrFaceTransformations(i);
      if (Tr != NULL)
      {
         for (int s=0; s<fes.Size(); ++s)
         {
            fe[s] = fes[s]->GetFE(Tr->Elem1No);
            fe2[s] = fes[s]->GetFE(Tr->Elem1No);

            fes[s]->GetElementVDofs(Tr->Elem1No, *(vdofs[s]));
            bx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_x[s]);
            bdx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_dx[s]);
            if (doftrans[s])
            {
               MFEM_WARNING("NavStoLSForm::Doftrans");
               doftrans[s]->InvTransformPrimal(*el_x[s]);
               doftrans[s]->InvTransformPrimal(*el_dx[s]);
            }
         }

         if (normalBC)
         {
            integrator.AssembleNormalBCVector(fe, fe2, *Tr,
                                              el_x_const, el_dx_const, el_y);
         }
         else
         {
            integrator.AssembleWeakDirBCVector(fe, fe2, *Tr,
                                               el_x_const, el_dx_const, el_y,
                                               blowingBC);
         }

         for (int s=0; s<fes.Size(); ++s)
         {
            if (el_y[s]->Size() == 0) { continue; }
            if (doftrans[s]) {doftrans[s]->TransformDual(*el_y[s]); }
            by.GetBlock(s).AddElementVector(*(vdofs[s]), *el_y[s]);
         }
      }
   }

   for (int s=0; s<fes.Size(); ++s)
   {
      delete vdofs2[s];
      delete vdofs[s];
      delete el_y[s];
      delete el_x[s];
      delete el_dx[s];
   }

   if (convectionOnly)
   {
      ConvectionOnlyMult(by);
   }

   by.SyncFromBlocks();
}

// Get Gradient
BlockOperator & NavStoLSForm::GetGradient(const Vector &x) const
{
   if (hasGrad && gradCalls <  maxGradCalls)
   {
      gradCalls++;
      return *pBlockGrad;
   }
   gradCalls = 0;

   if (pBlockGrad == NULL)
   {
      pBlockGrad = new BlockOperator(block_trueOffsets);
   }

   Array<const ParFiniteElementSpace *> pfes(fes.Size());

   for (int s1=0; s1<fes.Size(); ++s1)
   {
      pfes[s1] = ParFESpace(s1);

      for (int s2=0; s2<fes.Size(); ++s2)
      {
         phBlockGrad(s1,s2)->Clear();
      }
   }

   GetLocalGradient(x); // gradients are stored in 'Grads'

   if (fnfi.Size() > 0)
   {
      MFEM_ABORT("TODO: assemble contributions from shared face terms");
   }

   for (int s1=0; s1<fes.Size(); ++s1)
   {
      for (int s2=0; s2<fes.Size(); ++s2)
      {
         OperatorHandle dA(phBlockGrad(s1,s2)->Type()),
                        Ph(phBlockGrad(s1,s2)->Type()),
                        Rh(phBlockGrad(s1,s2)->Type());

         if (s1 == s2)
         {
            dA.MakeSquareBlockDiag(pfes[s1]->GetComm(), pfes[s1]->GlobalVSize(),
                                   pfes[s1]->GetDofOffsets(), Grads(s1,s1));
            Ph.ConvertFrom(pfes[s1]->Dof_TrueDof_Matrix());
            phBlockGrad(s1,s1)->MakePtAP(dA, Ph);

            OperatorHandle Ae;
            Ae.EliminateRowsCols(*phBlockGrad(s1,s1), *ess_tdofs[s1]);
         }
         else
         {
            dA.MakeRectangularBlockDiag(pfes[s1]->GetComm(),
                                        pfes[s1]->GlobalVSize(),
                                        pfes[s2]->GlobalVSize(),
                                        pfes[s1]->GetDofOffsets(),
                                        pfes[s2]->GetDofOffsets(),
                                        Grads(s1,s2));
            Rh.ConvertFrom(pfes[s1]->Dof_TrueDof_Matrix());
            Ph.ConvertFrom(pfes[s2]->Dof_TrueDof_Matrix());

            phBlockGrad(s1,s2)->MakeRAP(Rh, dA, Ph);

            phBlockGrad(s1,s2)->EliminateRows(*ess_tdofs[s1]);
            phBlockGrad(s1,s2)->EliminateCols(*ess_tdofs[s2]);
         }

         pBlockGrad->SetBlock(s1, s2, phBlockGrad(s1,s2)->Ptr());
      }
   }
   hasGrad = true;
   return *pBlockGrad;
}

// Return the local gradient matrix for the given true-dof vector x
const BlockOperator& NavStoLSForm
::GetLocalGradient(const Vector &dx) const
{
   // dxs_true is not modified, so const_cast is okay
   dxs_true.Update(const_cast<Vector &>(dx), block_trueOffsets);
   xs.Update(block_offsets);
   dxs.Update(block_offsets);

   fes[0]->GetProlongationMatrix()->Mult(
      dxs_true.GetBlock(0), dxs.GetBlock(0));

   add(xs0,dt,dxs.GetBlock(0),xs.GetBlock(0));   // x = x0 + dt*dx

   fes[1]->GetProlongationMatrix()->Mult(
      dxs_true.GetBlock(1), xs.GetBlock(1));

   fes[2]->GetProlongationMatrix()->Mult(
      dxs_true.GetBlock(2), dxs.GetBlock(2));

   add(xs2,dt,dxs.GetBlock(2),xs.GetBlock(2));   // x = x0 + dt*dx

   // (re)assemble Grad without b.c. into 'Grads'
   ComputeGradientBlocked(xs, dxs);

   delete BlockGrad;
   BlockGrad = new BlockOperator(block_offsets);

   for (int i = 0; i < fes.Size(); ++i)
   {
      for (int j = 0; j < fes.Size(); ++j)
      {
         BlockGrad->SetBlock(i, j, Grads(i, j));
      }
   }

   return *BlockGrad;
}

// Specialized version of GetGradient() for BlockVector
void NavStoLSForm
::ComputeGradientBlocked(const BlockVector &bx,
                         const BlockVector &bdx) const
{
   const int skip_zeros = 0;
   Array<Array<int> *> vdofs(fes.Size());
   Array<Array<int> *> vdofs2(fes.Size());
   Array<Vector *> el_x(fes.Size());
   Array<const Vector *> el_x_const(fes.Size());
   Array<Vector *> el_dx(fes.Size());
   Array<const Vector *> el_dx_const(fes.Size());
   Array2D<DenseMatrix *> elmats(fes.Size(), fes.Size());
   Array<const FiniteElement *>fe(fes.Size());
   Array<const FiniteElement *>fe2(fes.Size());
   ElementTransformation * T;
   FaceElementTransformations *Tr;
   Array<DofTransformation *> doftrans(fes.Size()); doftrans = nullptr;
   Mesh *mesh = fes[0]->GetMesh();

   for (int i=0; i<fes.Size(); ++i)
   {
      el_x_const[i] = el_x[i] = new Vector();
      el_dx_const[i] = el_dx[i] = new Vector();
      vdofs[i] = new Array<int>;
      vdofs2[i] = new Array<int>;
      for (int j=0; j<fes.Size(); ++j)
      {
         elmats(i,j) = new DenseMatrix();
      }
   }

   for (int i=0; i<fes.Size(); ++i)
   {
      for (int j=0; j<fes.Size(); ++j)
      {
         if (Grads(i,j) != NULL)
         {
            *Grads(i,j) = 0.0;
         }
         else
         {
            Grads(i,j) = new SparseMatrix(fes[i]->GetVSize(),
                                          fes[j]->GetVSize());
         }
      }
   }

   // Domain interior
   for (int i = 0; i < fes[0]->GetNE(); ++i)
   {
      T = fes[0]->GetElementTransformation(i);
      for (int s = 0; s < fes.Size(); ++s)
      {
         fe[s] = fes[s]->GetFE(i);
         doftrans[s] = fes[s]->GetElementVDofs(i, *vdofs[s]);
         bx.GetBlock(s).GetSubVector(*vdofs[s], *el_x[s]);
         bdx.GetBlock(s).GetSubVector(*vdofs[s], *el_dx[s]);
         if (doftrans[s])
         {
            MFEM_WARNING("NavStoLSForm::Doftrans");
            doftrans[s]->InvTransformPrimal(*el_x[s]);
            doftrans[s]->InvTransformPrimal(*el_dx[s]);
         }
      }

      integrator.AssembleElementGrad(fe, *T, el_x_const, el_dx_const, elmats);

      for (int j=0; j<fes.Size(); ++j)
      {
         for (int l=0; l<fes.Size(); ++l)
         {
            if (elmats(j,l)->Height() == 0) { continue; }
            if (doftrans[j] || doftrans[l])
            {
               MFEM_WARNING("NavStoLSForm::Doftrans");
               TransformDual(*doftrans[j], *doftrans[l], *elmats(j,l));
            }
            Grads(j,l)->AddSubMatrix(*vdofs[j], *vdofs[l],
                                     *elmats(j,l), skip_zeros);
         }
      }
   }

   // Domain boundary Outflow
   for (int i = 0; i < mesh->GetNBE(); ++i)
   {
      // Determine if boundary is outflow
      const int bdr_attr = mesh->GetBdrAttribute(i);
      bool outflowBC = false;
      bool suctionBC = false;
      for (int b=0; b<outflowBdr.Size(); ++b)
      {
         if ( bdr_attr == outflowBdr[b]) { outflowBC = true; }
      }
      for (int b=0; b<suctionBdr.Size(); ++b)
      {
         if ( bdr_attr == suctionBdr[b]) { suctionBC = true; }
      }

      if ( !outflowBC && !suctionBC ) { continue; }

      Tr = mesh->GetBdrFaceTransformations(i);
      if (Tr != NULL)
      {
         for (int s = 0; s < fes.Size(); ++s)
         {
            fe[s] = fes[s]->GetFE(Tr->Elem1No);
            fe2[s] = fe[s];

            fes[s]->GetElementVDofs(Tr->Elem1No, *vdofs[s]);
            bx.GetBlock(s).GetSubVector(*vdofs[s], *el_x[s]);

            bx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_x[s]);
            bdx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_dx[s]);
            if (doftrans[s])
            {
               MFEM_WARNING("NavStoLSForm::Doftrans");
               doftrans[s]->InvTransformPrimal(*el_x[s]);
               doftrans[s]->InvTransformPrimal(*el_dx[s]);
            }
         }

         integrator.AssembleOutflowGrad(fe, fe2, *Tr,
                                        el_x_const, el_dx_const, elmats,
                                        suctionBC);

         for (int l=0; l<fes.Size(); ++l)
         {
            for (int j=0; j<fes.Size(); ++j)
            {
               if (elmats(j,l)->Height() == 0) { continue; }
               if (doftrans[j] || doftrans[l])
               {
                  MFEM_WARNING("NavStoLSForm::Doftrans");
                  TransformDual(*doftrans[j], *doftrans[l], *elmats(j,l));
               }
               Grads(j,l)->AddSubMatrix(*vdofs[j], *vdofs[l],
                                        *elmats(j,l), skip_zeros);
            }
         }
      }
   }

   // Domain boundary weak Dirichlet BC
   for (int i = 0; i < mesh->GetNBE(); ++i)
   {
      // Determine if boundary is outflow
      const int bdr_attr = mesh->GetBdrAttribute(i);
      bool weakBC = false;
      bool blowingBC = false;
      bool normalBC = false;
      for (int b=0; b<weakBCBdr.Size(); ++b)
      {
         if ( bdr_attr == weakBCBdr[b]) { weakBC = true; }
      }
      for (int b=0; b<blowingBdr.Size(); ++b)
      {
         if ( bdr_attr == blowingBdr[b]) { blowingBC = true; }
      }
      for (int b=0; b<normalBCBdr.Size(); ++b)
      {
         if ( bdr_attr == normalBCBdr[b]) { normalBC = true; }
      }
      if ( !weakBC && !blowingBC && !normalBC) { continue; }

      Tr = mesh->GetBdrFaceTransformations(i);
      if (Tr != NULL)
      {
         for (int s = 0; s < fes.Size(); ++s)
         {
            fe[s] = fes[s]->GetFE(Tr->Elem1No);
            fe2[s] = fe[s];

            fes[s]->GetElementVDofs(Tr->Elem1No, *vdofs[s]);
            bx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_x[s]);
            bdx.GetBlock(s).GetSubVector(*(vdofs[s]), *el_dx[s]);
            if (doftrans[s])
            {
               MFEM_WARNING("NavStoLSForm::Doftrans");
               doftrans[s]->InvTransformPrimal(*el_x[s]);
               doftrans[s]->InvTransformPrimal(*el_dx[s]);
            }
         }

         if (normalBC)
         {
            integrator.AssembleNormalBCGrad(fe, fe2, *Tr,
                                            el_x_const, el_dx_const, elmats);
         }
         else
         {
            integrator.AssembleWeakDirBCGrad(fe, fe2, *Tr,
                                             el_x_const, el_dx_const, elmats,
                                             blowingBC);
         }


         for (int l=0; l<fes.Size(); ++l)
         {
            for (int j=0; j<fes.Size(); ++j)
            {
               if (elmats(j,l)->Height() == 0) { continue; }
               if (doftrans[j] || doftrans[l])
               {
                  MFEM_WARNING("NavStoLSForm::Doftrans");
                  TransformDual(*doftrans[j], *doftrans[l], *elmats(j,l));
               }
               Grads(j,l)->AddSubMatrix(*vdofs[j], *vdofs[l],
                                        *elmats(j,l), skip_zeros);
            }
         }
      }
   }

   if (convectionOnly)
   {
      ConvectionOnlyGrad();
   }


   if (!Grads(0,0)->Finalized())
   {
      for (int i=0; i<fes.Size(); ++i)
      {
         for (int j=0; j<fes.Size(); ++j)
         {
            Grads(i,j)->Finalize(skip_zeros);
         }
      }
   }

   for (int i=0; i<fes.Size(); ++i)
   {
      for (int j=0; j<fes.Size(); ++j)
      {
         delete elmats(i,j);
      }
      delete vdofs2[i];
      delete vdofs[i];
      delete el_x[i];
      delete el_dx[i];
   }

}

void NavStoLSForm::ConvectionOnlyMult(BlockVector &by) const
{
   by.GetBlock(0) = 0.0;
   by.GetBlock(1) = 0.0;
}

void NavStoLSForm::ConvectionOnlyGrad() const
{
   *Grads(0,0) = 0.0;
   *Grads(0,1) = 0.0;
   *Grads(0,2) = 0.0;

   *Grads(1,0) = 0.0;
   *Grads(1,1) = 0.0;
   *Grads(1,2) = 0.0;

   *Grads(2,0) = 0.0;
   *Grads(2,1) = 0.0;
}


