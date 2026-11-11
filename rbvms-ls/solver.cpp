// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "solver.hpp"

using namespace mfem;
using namespace RBVMS;

// Evolution Constructor
Evolution2::Evolution2(ParTimeDepBlockNonlinForm &form,
                     IterativeSolver &solver)
   : TimeDependentOperator(form.Width(), 0.0, IMPLICIT),
     form(form), solver(solver), dudt(form.Width())
{
   solver.SetOperator(form);
   dudt = 0.0;
}

// Solve time dependent problem
void Evolution2::ImplicitSolve(const real_t dt,
                              const Vector &u0, Vector &dudt_)
{
   //form.ResetGradient();
   form.SetTimeAndSolution(t, dt, u0);
   Vector zero;
   // Initial guess id previous solution
   dudt = 0.0;
   solver.Mult(zero, dudt);
   dudt_ = dudt;
   if (Mpi::Root())
   {
      std::cout<<"\n\tTotal # Newton iterations = "
               <<solver.GetNumIterations()<<"\n\n";
   }
}

// Compute a norm for each component
void NewtonSystemSolver::Norms(const Vector &r, Vector& lnorm) const
{
   lnorm.SetSize(nvar);
   for (int i = 0; i < nvar; ++i)
   {
      Vector r_i(r.GetData() + bOffsets[i], bOffsets[i+1] - bOffsets[i]);
      lnorm[i] = sqrt(InnerProduct(MPI_COMM_WORLD, r_i, r_i));
      MFEM_VERIFY(IsFinite(lnorm[i]), "norm[" << i << "] = " << lnorm[i]);
   }
}

// Newton solver with a convergence check on each component
void NewtonSystemSolver::Mult(const Vector &b, Vector &x) const
{
   MFEM_VERIFY(oper != NULL, "the Operator is not set (use SetOperator).");
   MFEM_VERIFY(prec != NULL, "the Solver is not set (use SetSolver).");

   int it = 0;
   Vector norm0(nvar), norm(nvar), norm_goal(nvar);
   const bool have_b = (b.Size() == Height());

   if (!iterative_mode)
   {
      x = 0.0;
   }

   ProcessNewState(x);

   oper->Mult(x, r);
   if (have_b)
   {
      r -= b;
   }

   //   initial_norm
   Norms(r, norm0);
   norm = norm0;

   if (print_options.first_and_last && !print_options.iterations)
   {
      mfem::out << "Newton iteration " << std::setw(3) << it <<"\n"
                << " ||r||\n";
      for (int i = 0; i < nvar; ++i)
      {
         mfem::out<<std::setw(8)<<std::defaultfloat<<std::setprecision(4)
                  <<norm0[i]<<" %\n";
      }
   }

   for (int i = 0; i < nvar; ++i)
   {
      norm_goal[i] = std::max(rel_tol*norm0[i], abs_tol);
   }
   prec->iterative_mode = false;

   // x_{i+1} = x_i - [DF(x_i)]^{-1} [F(x_i)-b]
   for (it = 0; true; it++)
   {
      if (print_options.iterations)
      {
         mfem::out << "Newton iteration " << std::setw(3) << it <<"\n"
                   << " ||r||  \t"<< "||r||/||r_0||\n";
         for (int i = 0; i < nvar; ++i)
         {
            mfem::out<<std::setw(8)<<std::defaultfloat<<std::setprecision(4)
                     <<norm[i]<<"\t"
                     <<std::setw(8)<<std::fixed<<std::setprecision(2)
                     <<100*norm[i]/norm0[i]<<" %\n";
         }
      }
      Monitor(it, -1.0, r, x);
      converged = true;
      for (int i = 0; i < nvar; ++i)
      {
         if (norm[i] > norm_goal[i])
         {
            converged = false;
         }
      }
      if (converged) { break; }

      if (it >= max_iter)
      {
         converged = false;
         break;
      }

      grad = &oper->GetGradient(x);
      prec->SetOperator(*grad);

      if (lin_rtol_type)
      {
         AdaptiveLinRtolPreSolve(x, it, norm.Norml2());
      }

      prec->Mult(r, c); // c = [DF(x_i)]^{-1} [F(x_i)-b]

      if (lin_rtol_type)
      {
         AdaptiveLinRtolPostSolve(c, r, it, norm.Norml2());
      }

      const real_t c_scale = ComputeScalingFactor(x, b);
      if (c_scale == 0.0)
      {
         converged = false;
         break;
      }
      add(x, -c_scale, c, x);

      ProcessNewState(x);

      oper->Mult(x, r);
      if (have_b)
      {
         r -= b;
      }
      Norms(r, norm);
   }

   final_iter = it;
   final_norm = norm.Norml2();

   if (print_options.summary || (!converged && print_options.warnings) ||
       print_options.first_and_last)
   {
      mfem::out << "Newton iteration " << std::setw(3) << it <<"\n"
                << " ||r||  \t"<< "||r||/||r_0||\n";
      for (int i = 0; i < nvar; ++i)
      {
         mfem::out<<std::setw(8)<<std::defaultfloat<<std::setprecision(4)
                  <<norm[i]<<"\t"
                  <<std::setw(8)<<std::fixed<<std::setprecision(2)
                  <<100*norm[i]/norm0[i]<<" %\n";
      }
   }
   if (!converged && (print_options.summary || print_options.warnings))
   {
      mfem::out << "Newton: No convergence!\n";
   }
}


// Print residual & relative residual
void GeneralResidualMonitor::MonitorResidual(int it,
                                             real_t norm,
                                             const Vector &r,
                                             bool final)
{
   if (interval < 0) { return; }
   if (it == 0) { norm0 = norm; }

   if ( ( it%interval == 0) || final )
   {
      mfem::out<<prefix<<" iteration "<<std::setw(3)<<it
               <<std::setw(8)<<std::defaultfloat<<std::setprecision(3)
               <<": ||r|| = "<<norm
               <<std::setw(6)<<std::fixed<<std::setprecision(2)
               <<", ||r||/||r_0|| = "<<100*norm/norm0<<" %\n";
   }
}

// Set the diagonal and off-diagonal operators
void JacobianPreconditioner::SetOperator(const Operator &op)
{
   BlockOperator *jacobian = (BlockOperator *) &op;

   for (int i = 0; i < prec.Size(); ++i)
   {
      if (prec[i])
      {
         prec[i]->SetOperator(jacobian->GetBlock(i,i));
         SetDiagonalBlock(i, prec[i]);
      }
      for (int j = i+1; j < prec.Size(); ++j)
      {
         SetBlock(j,i, &jacobian->GetBlock(j,i));
      }
   }
}
