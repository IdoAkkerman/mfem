// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_LS_EVOLUTION_HPP
#define RBVMS_LS_EVOLUTION_HPP

#include "mfem.hpp"
#include "integrator.hpp"
#include "solver.hpp"

using namespace std;
using namespace mfem;

namespace RBVMS
{


/** This class is a specialized ParBlockNonlinearForm that includes timestepping
    interpolation, vis:
       add(x0,dt,dx,x);   // x = x0 + dt*dx

    Both x and dx need to be passed to the FormIntegrator.
*/
class NavStoLSForm : public ParTimeDepBlockNonlinForm
{
private:
   IncNavStoIntegrator &integrator;

   /// Numerical parameters
   real_t dt;
   mutable real_t cfl, outflow;

   /// Boundary parameters
   Array<int> strongBCBdr;
   Array<int> weakBCBdr;
   Array<int> normalBCBdr;
   Array<int> outflowBdr;
   Array<int> suctionBdr;
   Array<int> blowingBdr;

   /// Solution & Residual vector
   mutable Vector xs0, xs2;
   mutable BlockVector dxs;
   mutable BlockVector dxs_true;

   /// Conservative boundary forces
   mutable DenseMatrix bdrForce;
   mutable bool hasGrad;
   mutable int  gradCalls;
   // Convection only things
   bool convectionOnly = false;

   void ConvectionOnlyMult(BlockVector &by) const;
   void ConvectionOnlyGrad() const;

public:
   /// Constructor
   NavStoLSForm(Array<ParFiniteElementSpace *> &pfes,
                IncNavStoIntegrator &integ)
      : ParTimeDepBlockNonlinForm(pfes), integrator(integ), hasGrad(false)
   {
   }

   void SetStrongBC (Array<int> strong_bdr);
   void SetWeakBC   (Array<int> weak_bdr) { weak_bdr.Copy(weakBCBdr); };
   void SetNormalBC   (Array<int> nor_bdr) { nor_bdr.Copy(normalBCBdr);};
   void SetOutflowBC(Array<int> outflow_bdr) { outflow_bdr.Copy(outflowBdr);};
   void SetSuctionBC(Array<int> suction_bdr) { suction_bdr.Copy(suctionBdr);};
   void SetBlowingBC(Array<int> blowing_bdr) { blowing_bdr.Copy(blowingBdr);};

   void SetInconsistentDC(real_t k0) { integrator.SetInconsistentDC(k0); };
   void SetConsistentDC(real_t k1) { integrator.SetConsistentDC(k1); };

   void SetConvectionOnly () { convectionOnly = true; };

   /// Set the solution of the previous time step @a x0
   /// and the timestep size @a dt of the current solve.
   void SetTimeAndSolution(const real_t t,
                           const real_t dt,
                           const Vector &x0);

   void ResetGradient();

   int  maxGradCalls;

   /// Get the CFL-Number
   real_t GetCFL() { return cfl;};

   /// Get the Outflow
   real_t GetOutflow() { return outflow;};

   /// Get the conservative boundary forces
   DenseMatrix& GetForce() { return bdrForce;};

   /// Block T-Vector to Vector
   Vector GetEnergies(const Vector &x) const;

   /// Block T-Vector to Block T-Vector
   void Mult(const Vector &x, Vector &y) const;

   /// Specialized version of Mult() for BlockVectors
   /// Block L-Vector to Block L-Vector
   void MultBlocked(const BlockVector &bx,
                    const BlockVector &dbx,
                    BlockVector &by) const;

   virtual BlockOperator &GetGradient(const Vector &x) const;

   /// Return the local block gradient matrix for the given true-dof vector x
   const BlockOperator& GetLocalGradient(const Vector &x) const;

   /// Specialized version of GetGradient() for BlockVector
   void ComputeGradientBlocked(const BlockVector &bx,
                               const BlockVector &dbx) const;
};

} // namespace RBVMS

#endif

