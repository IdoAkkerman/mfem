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


class GlobalConstriant : public Operator
{
protected:
   Operator &M;
   Vector &c;
public:

   GlobalConstriant(Operator &M_, Vector &c_)
      : Operator(M_.Width()+1), M(M_), c(c_)
   {
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      Vector xr(x.GetData(), x.Size()-1);
      Vector yr(y.GetData(), y.Size()-1);

      cout<<xr.Size()<<endl;
      cout<<yr.Size()<<endl;

      M.Mult(xr, yr);
      yr.Add(x[x.Size()-1], c);
      y[y.Size()-1] = c*xr;
   }
};

class GlobalConstriant2 : public Operator
{
protected:
   Operator &M;
   Array<Vector *>  c;
public:

   GlobalConstriant2(Operator &M_, Array<Vector *> &c_)
      : Operator(M_.Width()+c_.Size()), M(M_), c(c_)
   {
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      Vector xr(x.GetData(), x.Size()-c.Size());
      Vector lambda(x.GetData()+xr.Size(), c.Size());
      Vector yr(y.GetData(), y.Size()-c.Size());
      Vector yc(y.GetData()+yr.Size(), c.Size());

      M.Mult(xr, yr);
      for (int i = 0; i < c.Size(); i++)
      {
         yr.Add(lambda[i], *c[i]);
         yc[i] = (*c[i])*xr;
      }
   }
};




real_t MatrixNorm(HypreParMatrix *M, Vector &x)
{
   Vector tmp(x.Size());
   M->Mult(x,tmp);
   return x*tmp;
}

real_t RayleighQuotient(HypreParMatrix *A,
                        HypreParMatrix *M, Vector &x)
{
   return MatrixNorm(M,x)/MatrixNorm(A,x);
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
         // ir = &IntRules.Get(el.GetGeomType(), oa * el.GetOrder() + ob);
         ir = & DiffusionIntegrator::GetRule(el,el);
      }

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);

         Tr.SetIntPoint (&ip);
         phi_gf.GetGradient(Tr, grad);
         real_t val = ip.weight*Tr.Weight()*(grad*grad);

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
   const char *ref_file  = "";
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
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
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

   // Insert knots from file
   if (strlen(ref_file) != 0)
   {
      mesh->RefineNURBSFromFile(ref_file);
   }

   // 4. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement (2 by default, or
   //    specified on the command line with -rs).
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }

   // Report on final mesh statistics
   mesh->PrintInfo();

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
      if (pmesh->NURBSext)
      {
         fec = new NURBSFECollection(order);
         ext = new NURBSExtension(pmesh->NURBSext, order);
      }
      else
      {
         fec = new H1_FECollection(order, dim);
      }
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

   FiniteElementCollection *fec2 = new H1_FECollection(order, dim);
   // FiniteElementCollection *fec2 = new L2_FECollection(1, dim);
   ParFiniteElementSpace *fespace2 = new ParFiniteElementSpace(pmesh, fec2);

   // 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
   //    element space. The first corresponds to the Laplacian operator -Delta,
   //    while the second is a simple mass matrix needed on the right hand side
   //    of the generalized eigenvalue problem below. The boundary conditions
   //    are implemented by elimination with special values on the diagonal to
   //    shift the Dirichlet eigenvalues out of the computational range. After
   //    serial and parallel assembly we extract the corresponding parallel
   //    matrices A and M.
   ParGridFunction x(fespace);
   ParGridFunction ca(fespace2);
   ca = 1.0;
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   ConstantCoefficient eps(1e-8);
   ConstantCoefficient h1_val(0.1);
   GridFunctionCoefficient c_cf(&ca);
   //AbsoluteCoefficient c_cfa(c_cf);
   MaxCoefficient c_cfa(c_cf, eps);
   Array<real_t> eigenvalues;

   real_t L2_ref = -1;
   real_t L1_ref = -1;

   // 8. Define and configure the LOBPCG eigensolver and the BoomerAMG
   //    preconditioner for A to be used within the solver. Set the matrices
   //    which define the generalized eigenproblem A x = lambda M x.
   ParBilinearForm *a = nullptr;
   ParBilinearForm *m = new ParBilinearForm(fespace);
   m->AddDomainIntegrator(new MassIntegrator(one));
   m->Assemble();
   m->Finalize();

   ParBilinearForm *m2 = new ParBilinearForm(fespace2);
   m2->AddDomainIntegrator(new MassIntegrator(one));
   m2->AddDomainIntegrator(new DiffusionIntegrator(h1_val));
   m2->Assemble();
   m2->Finalize();

   HypreParMatrix *A = nullptr;
   HypreParMatrix *M = m->ParallelAssemble();
   HypreParMatrix *M2 = m2->ParallelAssemble();

   Solver * precond = NULL;
   HypreEuclid *ilu = new HypreEuclid(*M);
   precond = ilu;

   Solver * precond2= NULL;
   HypreEuclid *ilu2 = new HypreEuclid(*M2);
   precond2 = ilu2;

   // Create nullspace
   MeanZeroSolver precond_proj(*precond, *M, *fespace);
   MeanZeroProjector  proj(*M, *fespace);

   MeanZeroSolver precond_proj2(*precond2, *M2, *fespace2);
   MeanZeroProjector  proj2(*M2, *fespace2);

   // Loop
   Array<ParGridFunction*> ax;
   VisItDataCollection *visit_dc = nullptr;
   int nev = 16;
   Vector diff(ca.Size());
   int ncluster0;
   for (int iter = 0; iter<=50; iter++)
   {
      if (myid == 0)
      {
         cout<<"\n\n----------------------------------------------\n";
         cout<<" iter  = "<<iter<<endl;
         cout<<"----------------------------------------------\n\n\n";
      }

      delete a;
      delete A;
      a = new ParBilinearForm(fespace);
      a->AddDomainIntegrator(new DiffusionIntegrator(c_cfa));
      a->Assemble();
      a->Finalize();
      A = a->ParallelAssemble();

      // 8. Define and configure the LOBPCG eigensolver and the BoomerAMG
      //    preconditioner for A to be used within the solver. Set the matrices
      //    which define the generalized eigenproblem A x = lambda M x.
      HypreLOBPCG  lobpcg(MPI_COMM_WORLD);
      lobpcg.SetNumModes(nev);
      lobpcg.SetRandomSeed(seed);
      //      lobpcg.SetPreconditioner(precond_proj);
      lobpcg.SetPreconditioner(*precond);
      lobpcg.SetMaxIter(1400);
      lobpcg.SetTol(1e-24);
      lobpcg.SetPrecondUsageMode(0);
      lobpcg.SetPrintLevel(0);
      lobpcg.SetMassMatrix(*A);
      lobpcg.SetOperator(*M);
      // lobpcg.SetSubSpaceProjector(proj);

      // 9. Compute the eigenmodes and extract the array of eigenvalues. Define a
      //    parallel grid function to represent each of the eigenmodes returned by
      //    the solver.
      lobpcg.Solve();
      lobpcg.GetEigenvalues(eigenvalues);

      // Check for clustering of eigenvalues
      real_t cluster_tol = 1e-2;
      int ii;
      for (ii = 1; ii < nev; ii++)
      {
         if ((eigenvalues[ii]/eigenvalues[0]) > 1.0 + cluster_tol) { break; }
      }
      MFEM_VERIFY(ii < nev, "Eigen space not large enough to get cluster");
      int ncluster = ii;
      cout<<"Cluster size = "<<ncluster<<endl;

      // Initialize if this is the first iteration
      if (iter == 0)
      {
         ca = eigenvalues[0];
         diff = ca;
         ncluster0 = ncluster;

         if (visit)
         {
            ax.SetSize(ncluster);
            for (ii = 0; ii < ncluster; ii++)
            {
               ax[ii] = new ParGridFunction(fespace);
            }
            visit_dc = new VisItDataCollection ("InverseEstimate", pmesh);
            visit_dc->RegisterField("coefficient", &ca);
            for (ii = 0; ii < ax.Size(); ii++)
            {
               std::ostringstream oss;
               oss << "mode_" << std::setw(3) << std::setfill('0') << ii;
               visit_dc->RegisterField(oss.str(), ax[ii]);
            }
         }
      }

      // Print eigenmodes and coefficient - if requested
      if (visit_dc)
      {
         for (ii = 0; ii < ax.Size(); ii++)
         {
            (*ax[ii]) = lobpcg.GetEigenvector(ii);
         }
         visit_dc->SetCycle(iter);
         visit_dc->Save();
      }

      // Solve constrained minimization
      int size = ncluster0;// + ncluster0;

      Array<Vector *> vv( size );
      ParLinearForm v(fespace2);
      v.AddDomainIntegrator(new DiffusionRQConstraintIntegrator(x));
      for (int i = 0; i < size ; i++)
      {
         x = lobpcg.GetEigenvector(i);
         v = 0.0;
         v.Assemble();
         vv[i] = new Vector(v);
      }

      GlobalConstriant2 Mc(*M2, vv);
      Vector sol(ca.Size()+size); sol = 0.0;
      Vector sol_r(sol.GetData(), ca.Size());
      Vector rhs(ca.Size()+size); rhs = 0.0;
      Vector rhs_r(rhs.GetData(), ca.Size());

      M2->Mult(ca, rhs_r);
      rhs *=-1;

      int max_iter = 200;
      real_t tol = 1e-8;
      GMRES(Mc,  sol, rhs, max_iter,200, tol, 0.0, 0);

      // Check if constraints are satisfied
      for (int j = 0; j < vv.Size() ; j++)
      {
         cout<<"Constriant Check =" <<(*vv[j])*sol_r<<endl;
         delete vv[j];
      }

      // Relaxation
      real_t fac = 0.5;
      ca *= 1.0 - fac;
      ca.Add(fac,sol_r);

      // Correction -- due to cllipping of coefficient
      delete a;
      a = new ParBilinearForm(fespace);
      a->AddDomainIntegrator(new DiffusionIntegrator(c_cfa));
      a->Assemble();
      a->Finalize();
      delete A;
      A = a->ParallelAssemble();
      real_t rr =RayleighQuotient(A,M,x);
      ca *= rr;

      // Report data
      real_t L2 = MatrixNorm(M2,ca);
      real_t L1 = ca.ComputeL1Error(zero);
      if (iter == 0 )
      {
         L2_ref = L2;
         L1_ref = L1;
      }
      cout<< "ITER "<<iter <<" : "<<rr<<" "<<eigenvalues[0]<<" --> check "<<rr
          <<"  L2 = "<<L2<<" <--- "<<L2_ref
          <<"  L1 = "<<L1<<" <--- "<<L1_ref<<endl;
      // Check convergence
      diff -= ca;
      real_t ca_err = diff.Norml2();
      cout<<"Diff = "<<ca_err<<endl;
      if (ca_err < 1e-10*ca.Norml2()) { break; }
      diff = ca;
   }

   real_t h = 1.0/sqrt(real_t(pmesh->GetNE()));

   cout<<" h = "<< h<<endl;
   ca /= (h*h);
   if (visit_dc)
   {
      visit_dc->SetCycle(101);
      visit_dc->Save();
   }


   // 12. Free the used memory.
   // Report on final mesh statistics
   if (fespace == fespace2)
   {
      delete fespace;
   }
   else
   {
      delete fespace;
      delete fespace2;
   }
   delete fec2;
   delete precond;
   delete precond2;
   delete m;
   delete m2;
   delete M;
   delete M2;
   delete a;
   delete A;
   delete visit_dc;
   if (order > 0)
   {
      delete fec;
   }
   delete pmesh;

   return 0;
}
