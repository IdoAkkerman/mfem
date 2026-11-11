// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "mfem.hpp"
#include <unistd.h>
#include "buildInfo.hpp"
using namespace mfem;

void line(int len)
{
   for (int b=0; b<len; ++b) { std::cout<<"-"; }
   std::cout<<"\n";
}


// Function for printing compile time and runtime information.
void printHeader()
{
      // Print header
      line(80);
      std::cout<<R"(  _________  ___  ___  ________          ________  ___  ___  ___  ________    )"<<std::endl;
      std::cout<<R"( |\___   ___\\  \|\  \|\   ___ \        |\   ____\|\  \|\  \|\  \|\   __  \   )"<<std::endl;
      std::cout<<R"( \|___ \  \_\ \  \\\  \ \  \_|\ \       \ \  \___|\ \  \\\  \ \  \ \  \|\  \  )"<<std::endl;
      std::cout<<R"(      \ \  \ \ \  \\\  \ \  \ \\ \       \ \_____  \ \   __  \ \  \ \   ____\ )"<<std::endl;
      std::cout<<R"(       \ \  \ \ \  \\\  \ \  \_\\ \       \|____|\  \ \  \ \  \ \  \ \  \___| )"<<std::endl;
      std::cout<<R"(        \ \__\ \ \_______\ \_______\        ____\_\  \ \__\ \__\ \__\ \__\    )"<<std::endl;
      std::cout<<R"(         \|__|  \|_______|\|_______|       |\_________\|__|\|__|\|__|\|__|    )"<<std::endl;
      std::cout<<R"(                                           \|_________|                       )"<<std::endl;
      std::cout<<R"(                                                                              )"<<std::endl;
}

// Function for printing compile time and runtime information.
void printSerialInfo()
{
      // Print header
      printHeader();

      // Build info
      line(80);
      std::cout<<"Compile time info\n";
      line(80);
      std::cout<< buildInfo.str() << std::endl;

      // Run info
      line(80);
      std::cout<<"Run time info"<<std::endl;
      line(80); std::cout<<std::endl;
      time_t     now = time(0);
      struct tm  tstruct = *localtime(&now);
      char       time[80], host[80];
      strftime(time, sizeof(time), "%Y-%m-%d.%X", &tstruct);
      gethostname(host,sizeof(host));

      std::cout<<"Time: "<<time<<std::endl;
      std::cout<<"Host: "<<host<<std::endl;

      std::cout<<std::endl; line(80); std::cout<<std::endl;
}

// Function for printing compile time and runtime information.
void printParallelInfo()
{
   if (Mpi::Root())
   {
      // Print header
      printHeader();

      // Build info
      line(80);
      std::cout<<"Compile time info\n";
      line(80);
      std::cout<< buildInfo.str() << std::endl;

      // Run info
      line(80);
      std::cout<<"Run time info"<<std::endl;
      line(80); std::cout<<std::endl;
      time_t     now = time(0);
      struct tm  tstruct = *localtime(&now);
      char       time[80], host[80];
      strftime(time, sizeof(time), "%Y-%m-%d.%X", &tstruct);
      gethostname(host,sizeof(host));

      std::cout<<"Time: "<<time<<std::endl;
      std::cout<<"Numer of MPI ranks "<<Mpi::WorldSize()<<std::endl;

      std::cout<<"List  of hosts\n0: "<<host<<std::endl;

      // Receive hostnames from non-root nodes
      for (int i = 1; i < Mpi::WorldSize(); i++)
      {
         MPI_Status status;
         MPI_Recv (&host, sizeof(host), MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
         std::cout<<i<<": "<<host<<std::endl;
      }

      std::cout<<std::endl; line(80); std::cout<<std::endl;
   }
   else
   {
      // Send hostnames from non-root nodes to root
      char host[80];
      gethostname(host,sizeof(host));
      MPI_Send (&host, sizeof(host), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
   }
}


void CheckBoundaries(Array<bool> &bnd_flags,
                     Array<int> &markers,
                     Array<int> &bc_bnds,
                     int mark = 1)
{
   int amax = bnd_flags.Size();
   markers.SetSize(amax);
   markers = 0;
   // Check strong boundaries
   for (int b = 0; b < bc_bnds.Size(); b++)
   {
      int bnd = bc_bnds[b]-1;
      if ( bnd < 0 || bnd >= amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flags[bnd])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flags[bnd] = true;
      markers[bnd] = mark;
   }
}

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
   virtual void SetOperator(const Operator &op)
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
   };

   // Destructor
   virtual ~JacobianPreconditioner()
   {
      for (int i = 0; i < prec.Size(); ++i)
      {
         if (prec[i]) { delete prec[i]; }
      }
   };
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
      interval = print_iv;
   }

   /// Print residual
   virtual void MonitorResidual(int it,
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
   };

};

