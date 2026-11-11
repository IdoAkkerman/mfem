// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_SOLVER_HPP
#define RBVMS_SOLVER_HPP

#include "mfem.hpp"


using namespace std;
using namespace mfem;

namespace RBVMS
{

// Predefine class
//class ParTimeDepBlockNonlinForm;
class ParTimeDepBlockNonlinForm : public ParBlockNonlinearForm
{
   bool hasGrad = false;

public:

   // Formulation Constructor
   ParTimeDepBlockNonlinForm(Array<ParFiniteElementSpace *> &pfes)
      : ParBlockNonlinearForm(pfes)
   {
   }

   virtual void SetTimeAndSolution(const real_t t,
                                   const real_t dt,
                                   const Vector &x0) {};

   void ResetGradient()
   {
      hasGrad = false;
   };
};


/** This class provide the correct interface between the time-dependent
    block nonlinear form (defined below) and the MFEM::ODESolver.
*/
class Evolution2 : public TimeDependentOperator
{
private:
   ParTimeDepBlockNonlinForm &form;
   IterativeSolver &solver;
   Vector dudt;

public:
   /// Constructor
   Evolution2(ParTimeDepBlockNonlinForm &form,
             IterativeSolver &solver);

   /// Stub for explicit solve of time dependent problem
   virtual void Mult(const Vector &x, Vector &k) const override { k = 0.0;};

   /// Solve time dependent problem
   virtual void ImplicitSolve(const real_t dt,
                              const Vector &x,
                              Vector &k) override;

   void SetTime(const real_t t_) override { t = t_; form.ResetGradient(); };

   /// Destructor
   ~Evolution2() {}
};


/// Newton's method for solving F(x)=b for a given operator F.
/** The method GetGradient() must be implemented for the operator F.
    The preconditioner is used (in non-iterative mode) to evaluate
    the action of the inverse gradient of the operator. */
class NewtonSystemSolver : public NewtonSolver
{
private:
   Array<int> &bOffsets;
   int nvar;

   // Compute the norms for each component
   void Norms(const Vector &r,Vector& norm) const;

public:
   // Constructor, provide MPI context and component subdivision
   NewtonSystemSolver(MPI_Comm comm_, Array<int> &offsets)
      : NewtonSolver(comm_), bOffsets(offsets)
   {
      nvar = bOffsets.Size()-1;
   }

   /// Solve the nonlinear system with right-hand side @a b.
   /** If `b.Size() != Height()`, then @a b is assumed to be zero. */
   virtual void Mult(const Vector &b, Vector &x) const;
};


/// This class help monitor the convergence of the linear Krylov solve.
class GeneralResidualMonitor : public IterativeSolverMonitor
{
private:
   const std::string prefix;
   int interval;
   mutable real_t norm0;

public:
   /// Constructor
   GeneralResidualMonitor( const std::string& prefix_,
                           int print_iv)
      : prefix(prefix_), interval(-1)
   {
      if (Mpi::Root()) { interval = print_iv; }
   }

   /// Print residual
   virtual void MonitorResidual(int it,
                                real_t norm,
                                const Vector &r,
                                bool final);

};

// Custom block preconditioner for the Jacobian
class JacobianPreconditioner : public
   BlockLowerTriangularPreconditioner
{
protected:
   Array<Solver *> prec;

public:
   /// Constructor
   JacobianPreconditioner(Array<int> &offsets)
      : BlockLowerTriangularPreconditioner (offsets), prec(offsets.Size()-1)
   { prec = nullptr;};

   /// SetPreconditioners
   void SetPreconditioner(int i, Solver *pc)
   { prec[i] = pc; };

   /// Set the diagonal and off-diagonal operators
   virtual void SetOperator(const Operator &op);

   // Destructor
   virtual ~JacobianPreconditioner()
   {
      for (int i = 0; i < prec.Size(); ++i)
      {
         if (prec[i]) { delete prec[i]; }
      }
   };
};


} // namespace RBVMS

#endif

