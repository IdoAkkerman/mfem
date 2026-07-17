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
//       Inverse Estimate Miniapp: Compute the inverse estimate for a domain
//       --------------------------------------------------------------
//
// This miniapp computes the eigenvalue problem -Delta u = lambda u with homogeneous
// Neumann boundary conditions. We compute the lowest eigenmodes by discretizing
// the Laplacian and Mass operators using a FE space of the specified order, or
// an isoparametric/isogeometric space if order < 1 (quadratic for quadratic
// curvilinear mesh, NURBS for NURBS mesh, etc.)
//
// Compile with: make inverse-estimate
//
// Sample runs:  mpirun -np 4 inverse-estimate -m ../data/square-disc.mesh -rs 4
//               mpirun -np 4 inverse-estimate -m ../data/star.mesh -rs 4


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


real_t RayleighQuotient(HypreParMatrix *A,
                        HypreParMatrix *M, Vector &x)
{
   // HyprePar
   Vector tmp(x.Size());

   A->Mult(x,tmp);
   real_t xAx = x*tmp;

   M->Mult(x,tmp);
   real_t xMx = x*tmp;

   return xMx/xAx;
}

real_t NormL2(HypreParMatrix *M, Vector &x)
{
   Vector tmp(x.Size());
   M->Mult(x,tmp);
   return x*tmp;
}


/// Abstract base class LinearFormIntegrator
class DiffusionRQConstraintIntegrator : public LinearFormIntegrator
{
protected:
   GridFunction &phi_gf;
   Vector shape;
   int oa, ob;
public:

   DiffusionRQConstraintIntegrator(GridFunction &gf,
                                   const IntegrationRule *ir = NULL)
      : LinearFormIntegrator(ir), phi_gf(gf), oa(3), ob(2)
   {
   }

   /** Given a particular Finite Element and a transformation (Tr)
       computes the element vector, elvect. */
   void AssembleRHSElementVect(const FiniteElement &el,
                               ElementTransformation &Tr,
                               Vector &elvect) override
   {
      int dof = el.GetDof();
      int dim = el.GetDim();

      shape.SetSize(dof);       // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;

      Vector grad(dim);

      const IntegrationRule *ir = GetIntegrationRule(el, Tr);

      if (ir == NULL)
      {
         // ir = &IntRules.Get(el.GetGeomType(),
         //                    oa * el.GetOrder() + ob + Tr.OrderW());
         ir = &IntRules.Get(el.GetGeomType(), oa * el.GetOrder() + ob);
      }

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);

         Tr.SetIntPoint (&ip);
         phi_gf.GetGradient(Tr, grad);
         real_t val = ip.weight*Tr.Weight()*grad.Norml2();

         el.CalcPhysShape(Tr, shape);

         add(elvect, val, shape, elvect);
      }
   }
   ~DiffusionRQConstraintIntegrator() { }
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
   // delete mesh;

   // 6. Define a parallel finite element space on the parallel mesh. Here we
   //    use continuous Lagrange finite elements of the specified order. If
   //    order < 1, we instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   NURBSExtension *ext = nullptr;
   NURBSExtension *ext2 = nullptr;
   if (order > 0)
   {
      if (mesh->NURBSext)
      {
         fec = new NURBSFECollection(order);
         ext = new NURBSExtension(pmesh->NURBSext, order);
         ext2 = new NURBSExtension(mesh->NURBSext, order);
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
   FiniteElementSpace *fespace_ser = new FiniteElementSpace(mesh, ext2, fec);
   /*   HYPRE_BigInt size = fespace->GlobalTrueVSize();
      if (myid == 0)
      {
         cout << "Number of unknowns: " << size << endl;
      }*/

   // 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
   //    element space. The first corresponds to the Laplacian operator -Delta,
   //    while the second is a simple mass matrix needed on the right hand side
   //    of the generalized eigenvalue problem below. The boundary conditions
   //    are implemented by elimination with special values on the diagonal to
   //    shift the Dirichlet eigenvalues out of the computational range. After
   //    serial and parallel assembly we extract the corresponding parallel
   //    matrices A and M.
   ParGridFunction x(fespace);
   ParGridFunction ca(fespace);
   ca = 0.0;
   ConstantCoefficient one(1.0);
   GridFunctionCoefficient c_cf(&ca);

   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));
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
   HypreEuclid *ilu = new HypreEuclid(*M);
   precond = ilu;

   // Create nullspace
   MeanZeroSolver precond_proj(*precond, *M, *fespace);
   MeanZeroProjector  proj(*M, *fespace);

   HypreLOBPCG * lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
   int nev = 10;
   lobpcg->SetNumModes(nev);
   lobpcg->SetRandomSeed(seed);
   lobpcg->SetPreconditioner(precond_proj);
   //lobpcg->SetPreconditioner(*precond);
   lobpcg->SetMaxIter(1400);
   lobpcg->SetTol(1e-18);
   lobpcg->SetPrecondUsageMode(0);
   lobpcg->SetPrintLevel(1);
   lobpcg->SetMassMatrix(*A);
   lobpcg->SetOperator(*M);
   lobpcg->SetSubSpaceProjector(proj);

   // 9. Compute the eigenmodes and extract the array of eigenvalues. Define a
   //    parallel grid function to represent each of the eigenmodes returned by
   //    the solver.
   Array<real_t> eigenvalues;
   lobpcg->Solve();
   lobpcg->GetEigenvalues(eigenvalues);


   if (myid == 0)
   {
      cout<<"\n\n----------------------------------------------\n";
      cout<<" inv Ev  = "<<1.0/eigenvalues[0]<<endl;
      cout<<" inv Ev  = "<<1.0/sqrt(eigenvalues[0])<<endl;
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
   VisItDataCollection *visit_dc = nullptr;
   if (visit)
   {

      visit_dc = new VisItDataCollection ("InverseEstimate", pmesh);
      visit_dc->RegisterField("solution", &x);
      visit_dc->RegisterField("coefficient", &ca);
      for (int i = 0; i < nev ; i++)
      {
         x = lobpcg->GetEigenvector(i);
         real_t rr = RayleighQuotient(A,M,x);
         ca = rr;
         visit_dc->SetCycle(i);
         visit_dc->Save();

         ParBilinearForm *ac = new ParBilinearForm(fespace);
         ac->AddDomainIntegrator(new DiffusionIntegrator(c_cf));
         ac->Assemble();
         ac->Finalize();
         HypreParMatrix *Ac = ac->ParallelAssemble();

         real_t L2 = NormL2(M,ca);

         cout<< i <<" : "<<rr<<" --> check "<<RayleighQuotient(Ac,M,
                                                               x)<<"  L2 = "<<L2<<endl;
      }
   }

   // Define LOR mesh and space
   Array<Vector *> points;
   NURBSPointSet points_lor = NURBSPointSet::DEMKO;
   fespace->GetNURBSext()->GetPointsCompr(points, points_lor);
   for (int i = 0; i < points.Size() ; i++)
   {
      points[i]->Print(cout,88888);
   }
   Mesh mesh_lor = mesh->GetLinearNURBSMesh(points);

   // Create LOR spaces
   FiniteElementCollection *fec_lor;
   NURBSExtension *ext_lor = nullptr;

   if (mesh_lor.NURBSext)
   {
      fec_lor = new NURBSFECollection(1);
      ext_lor = new NURBSExtension(mesh_lor.NURBSext, 1);
   }
   else
   {
      fec_lor = new H1_FECollection(1, dim);
   }
   FiniteElementSpace fespace_lor(&mesh_lor, ext_lor, fec_lor);

   // Define solution on LOR
   GridFunction x_lor(&fespace_lor);
   NURBSPatchMap map(ext_lor);
   const KnotVector *kv[3];
   map.SetPatchVertexMap(0, kv);

   for (int i = 0; i <kv[0]->GetNCP(); i++)
   {
      int sx = i%2? 1.0:-1.0;
      for (int j = 0; j <kv[1]->GetNCP(); j++)
      {
         int sy = j%2? 1.0:-1.0;
         x_lor[map(i,j)] = sx*sy;
      }
   }

   // Transfer solution
   GridTransfer *gt = new InterpolationGridTransfer(*fespace_ser, fespace_lor);
   const Operator &P = gt->BackwardOperator();
   P.Mult(x_lor,x);
   x = lobpcg->GetEigenvector(0);
   real_t c = RayleighQuotient(A,M,x);
   cout<< "Define solution RR = " <<" : "<<c<<endl;

   ParLinearForm v(fespace);
   v.AddDomainIntegrator(new DiffusionRQConstraintIntegrator(x));
   v.Assemble();

   ParGridFunction Mv(fespace);
   Mv = 0.0;
   PCG(*M, *precond, v, Mv, 1, 200, 1e-8, 0.0);

   real_t cv = ca*v;
   real_t vMv = v*Mv;

   cout<<"cv = "<<cv<<endl;
   cout<<"vMv = "<<vMv<<endl;

   ca = v;
   ca /= cv/vMv;

   real_t rr = RayleighQuotient(A,M,x);
   ParBilinearForm *ac = new ParBilinearForm(fespace);
   ac->AddDomainIntegrator(new DiffusionIntegrator(c_cf));
   ac->Assemble();
   ac->Finalize();
   HypreParMatrix *Ac = ac->ParallelAssemble();
   real_t rr2 =RayleighQuotient(Ac,M,x);
   cout<< "??" <<" : "<<rr<<" --> check "<<rr2<<endl;

   ca *= rr2;

   delete ac;
   delete Ac;
   ac = new ParBilinearForm(fespace);
   ac->AddDomainIntegrator(new DiffusionIntegrator(c_cf));
   ac->Assemble();
   ac->Finalize();
   Ac = ac->ParallelAssemble();

   rr2 =RayleighQuotient(Ac,M,x);
   real_t L2 = NormL2(M,ca);
   cout<< "??" <<" : "<<rr<<" --> check "<<rr2<<"  L2 = "<<L2<<endl;



   if (visit_dc)
   {
      visit_dc->SetCycle(nev+1);
      visit_dc->Save();
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
