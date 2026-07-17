// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
//
//       --------------------------------------------------------------
//       Poincare Miniapp: Compute the Poincare estimate for a domain
//       --------------------------------------------------------------
//
// This miniapp computes the eigenvalue problem -Delta u = lambda u with homogeneous
// Neumann boundary conditions. We compute the lowest eigenmodes by discretizing
// the Laplacian and Mass operators using a FE space of the specified order, or
// an isoparametric/isogeometric space if order < 1 (quadratic for quadratic
// curvilinear mesh, NURBS for NURBS mesh, etc.)
//
// Compile with: make poincare
//
// Sample runs:  mpirun -np 4 poincare -m ../data/square-disc.mesh -rs 4
//               mpirun -np 4 poincare -m ../data/star.mesh -rs 4


#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


class MeanZeroProjector : public Operator
{
protected:
   Operator &M;
   real_t cMc;
   mutable HypreParVector c;   // mutable to make the innerproduct work
   mutable HypreParVector Mx;
   ParGridFunction c_gf;
public:

   MeanZeroProjector(Operator &M_, ParFiniteElementSpace &fespace)
      : Operator(M_.Width()), M(M_), c(&fespace), Mx(&fespace), c_gf(&fespace)
   {
      ConstantCoefficient one(1.0);
      c_gf.ProjectCoefficient(one);
      c_gf.GetTrueDofs(c);

      M.Mult(c, Mx);
      cMc = InnerProduct(c, Mx);
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      M.Mult(x, Mx);
      real_t alpha =  InnerProduct(c, Mx)/ cMc;
      add(x, -alpha, c, y);
   }
};

class MeanZeroSolver : public Solver
{
protected:
   Operator &pc;
   MeanZeroProjector proj;
   mutable Vector xp, yp;

public:

   MeanZeroSolver(Operator &pc_, Operator &M_, ParFiniteElementSpace &fespace)
      : Solver(pc_.Height()), pc(pc_), proj(M_, fespace)
   {
      xp.SetSize(pc.Height());
      yp.SetSize(pc.Height());
   }

   void SetOperator(const Operator &op) override
   {
      mfem_error("Not defined");
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      proj.Mult(x, xp);
      pc.Mult(xp, yp);
      proj.Mult(yp, y);
   }
};

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   const char *mesh_file = "../data/star.mesh";
   int ref_levels = 2;
   int order = 1;
   int seed = 75;
   bool slu_solver  = false;
   bool sp_solver = false;
   bool cpardiso_solver = false;
   bool visualization = true;
   bool visit         = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&seed, "-s", "--seed",
                  "Random seed used to initialize LOBPCG.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit",
                  "--no-visit",
                  "Enable or disable VisIt output.");
   args.Parse();
   if (slu_solver && sp_solver)
   {
      if (myid == 0)
         cout << "WARNING: Both SuperLU and STRUMPACK have been selected,"
              << " please choose either one." << endl
              << "         Defaulting to SuperLU." << endl;
      sp_solver = false;
   }
   // The command line options are also passed to the STRUMPACK
   // solver. So do not exit if some options are not recognized.
   if (!sp_solver)
   {
      if (!args.Good())
      {
         if (myid == 0)
         {
            args.PrintUsage(cout);
         }
         return 1;
      }
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // 3. Read the (serial) mesh from the given mesh file on all processors. We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 4. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement (2 by default, or
   //    specified on the command line with -rs).
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }

   // 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution (1 time by
   //    default, or specified on the command line with -rp). Once the parallel
   //    mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   // 6. Define a parallel finite element space on the parallel mesh. Here we
   //    use continuous Lagrange finite elements of the specified order. If
   //    order < 1, we instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   NURBSExtension *ext = nullptr;
   if (order > 0)
   {
      if (mesh->NURBSext)
      {
         fec = new NURBSFECollection(order);
         ext = new NURBSExtension(pmesh->NURBSext, order);
      }
      else
      {
         fec = new H1_FECollection(order, dim);
      }
   }
   else if (pmesh->GetNodes())
   {
      fec = pmesh->GetNodes()->OwnFEC();
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
   }
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, ext, fec);
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of unknowns: " << size << endl;
   }

   // 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
   //    element space. The first corresponds to the Laplacian operator -Delta,
   //    while the second is a simple mass matrix needed on the right hand side
   //    of the generalized eigenvalue problem below. The boundary conditions
   //    are implemented by elimination with special values on the diagonal to
   //    shift the Dirichlet eigenvalues out of the computational range. After
   //    serial and parallel assembly we extract the corresponding parallel
   //    matrices A and M.
   ConstantCoefficient one(1.0);

   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));
   if (pmesh->bdr_attributes.Size() == 0)
   {
      // Add a mass term if the mesh has no boundary, e.g. periodic mesh or
      // closed surface.
      a->AddDomainIntegrator(new MassIntegrator(one));
   }
   a->Assemble();
   a->Finalize();

   ParBilinearForm *m = new ParBilinearForm(fespace);
   m->AddDomainIntegrator(new MassIntegrator(one));
   m->Assemble();
   m->Finalize();

   HypreParMatrix *A = a->ParallelAssemble();
   HypreParMatrix *M = m->ParallelAssemble();

   delete a;
   delete m;

   // 8. Define and configure the LOBPCG eigensolver and the BoomerAMG
   //    preconditioner for A to be used within the solver. Set the matrices
   //    which define the generalized eigenproblem A x = lambda M x.
   Solver * precond = NULL;
   HypreEuclid *ilu = new HypreEuclid(*A);
   precond = ilu;

   // Create nullspace
   MeanZeroSolver precond_proj(*precond, *M, *fespace);
   MeanZeroProjector  proj(*M, *fespace);

   HypreLOBPCG * lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
   lobpcg->SetNumModes(1);
   lobpcg->SetRandomSeed(seed);
   lobpcg->SetPreconditioner(precond_proj);
   lobpcg->SetMaxIter(400);
   lobpcg->SetTol(1e-8);
   lobpcg->SetPrecondUsageMode(0);
   lobpcg->SetPrintLevel(1);
   lobpcg->SetMassMatrix(*M);
   lobpcg->SetOperator(*A);
   lobpcg->SetSubSpaceProjector(proj);

   // 9. Compute the eigenmodes and extract the array of eigenvalues. Define a
   //    parallel grid function to represent each of the eigenmodes returned by
   //    the solver.
   Array<real_t> eigenvalues;
   lobpcg->Solve();
   lobpcg->GetEigenvalues(eigenvalues);
   ParGridFunction x(fespace);

   if (myid == 0)
   {
      cout<<"\n\n----------------------------------------------\n";
      cout<<" Poincare constant = "<<1.0/sqrt(eigenvalues[0])<<endl;
      cout<<"----------------------------------------------\n\n\n";
   }

   // 10. Save the refined mesh and the modes in parallel. This output can be
   //     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
   {
      ostringstream mesh_name, mode_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);


      // convert eigenvector from HypreParVector to ParGridFunction
      x = lobpcg->GetEigenvector(0);

      mode_name << "mode_" << setfill('0') << setw(2) << 0 << "."
                << setfill('0') << setw(6) << myid;

      ofstream mode_ofs(mode_name.str().c_str());
      mode_ofs.precision(8);
      x.Save(mode_ofs);
      mode_name.str("");

   }

   // 11. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream mode_sock(vishost, visport);
      mode_sock.precision(8);

      // convert eigenvector from HypreParVector to ParGridFunction
      x = lobpcg->GetEigenvector(0);

      mode_sock << "parallel " << num_procs << " " << myid << "\n"
                << "solution\n" << *pmesh << x << flush
                << " Lambda = " << eigenvalues[0] << "'" << endl;

      mode_sock.close();
   }

   // Save mesh in the VisIt format
   if (visit)
   {

      VisItDataCollection visit_dc("Poincare", pmesh);
      visit_dc.RegisterField("solution", &x);
      x = lobpcg->GetEigenvector(0);
      visit_dc.SetCycle(0);
      visit_dc.Save();
   }

   // 12. Free the used memory.
   delete lobpcg;
   delete precond;
   delete M;
   delete A;

   delete fespace;
   if (order > 0)
   {
      delete fec;
   }
   delete pmesh;

   return 0;
}
