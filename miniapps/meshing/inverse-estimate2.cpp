// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you ch2n redistribute it and/or modify it under the
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
   DenseMatrix cMc;
   mutable HypreParVector c0;   // mutable to make the innerproduct work
   mutable HypreParVector cx;   // mutable to make the innerproduct work
   mutable HypreParVector cy;   // mutable to make the innerproduct work
   mutable HypreParVector Mx;
   ParGridFunction c_gf;
public:

   MeanZeroProjector(Operator &M_, ParFiniteElementSpace &fes)
      : Operator(M_.Width()), M(M_),
        c0(&fes), cx(&fes), cy(&fes), Mx(&fes), cMc(3)
   {
      ParGridFunction c_gf(&fes);
      HypreParVector c_true;

      ConstantCoefficient one(1.0);
      c_gf.ProjectCoefficient(one);
      c_gf.GetTrueDofs(c_true);
      c0 = c_true;

      CartesianXCoefficient x_cf;
      c_gf.ProjectCoefficient(x_cf);
      c_gf.GetTrueDofs(c_true);
      cx = c_true;

      CartesianYCoefficient y_cf;
      c_gf.ProjectCoefficient(y_cf);
      c_gf.GetTrueDofs(c_true);
      cy = c_true;

      M.Mult(c0, Mx);

      Vector tmp;
      tmp = c0;
      tmp.Print();

      tmp = Mx;
      tmp.Print();
      cMc(0,0) = InnerProduct(c0, Mx);
      cMc(1,0) = InnerProduct(cx, Mx);
      cMc(2,0) = InnerProduct(cy, Mx);

      M.Mult(cx, Mx);
      cMc(0,1) = InnerProduct(c0, Mx);
      cMc(1,1) = InnerProduct(cx, Mx);
      cMc(2,1) = InnerProduct(cy, Mx);

      M.Mult(cy, Mx);
      cMc(0,2) = InnerProduct(c0, Mx);
      cMc(1,2) = InnerProduct(cx, Mx);
      cMc(2,2) = InnerProduct(cy, Mx);
      cMc.Print();
      cMc.Invert();
      cMc.Print();
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      M.Mult(x, Mx);

      Vector rhs(3);

      rhs(0) =  InnerProduct(c0, Mx);
      rhs(1) =  InnerProduct(cx, Mx);
      rhs(2) =  InnerProduct(cy, Mx);

      Vector alpha(3);
      cMc.Mult(rhs, alpha);

      add(x, -alpha(0), c0, y);
      y.Add( -alpha(1), cx);
      y.Add( -alpha(2), cy);
   }
};

/*class SymmetryProjector : public Operator
{
protected:
   Operator &op;
public:

   SymmetryProjector(Operator &op_, ParFiniteElementSpace &fes)
      : Operator(), op(op_)
   {



   }

   void Mult(const Vector &x, Vector &y) const override
   {
   }
};*/


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


      Vector laps( ir->GetNPoints());
      phi_gf.GetLaplacians(Tr.ElementNo, *ir, laps);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);

         Tr.SetIntPoint (&ip);
         phi_gf.GetGradient(Tr, grad);
         real_t val = ip.weight*Tr.Weight()*laps[i]*laps[i];

         el.CalcPhysShape(Tr, shape);

         add(elvect, val, shape, elvect);
      }
   }
   ~DiffusionRQConstraintIntegrator() { }
};

class LapLapIntegrator: public BilinearFormIntegrator
{
   Coefficient *Q;
#ifndef MFEM_THREAD_SAFE
   Vector lshape;
#endif
public:
   LapLapIntegrator(Coefficient &q, const IntegrationRule *ir = nullptr)
   {
      Q = &q;
   }

   void AssembleElementMatrix(const FiniteElement &el,
                              ElementTransformation &Trans,
                              DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      // int dim = el.GetDim();
      real_t w;

#ifdef MFEM_THREAD_SAFE
      Vector lshape;
#endif
      elmat.SetSize(nd);
      lshape.SetSize(nd);

      const IntegrationRule *ir = GetIntegrationRule(el, Trans);

      elmat = 0.0;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint (&ip);

         el.CalcPhysLaplacian(Trans, lshape);

         w = Trans.Weight() * ip.weight;
         if (Q)
         {
            w *= Q -> Eval(Trans, ip);
         }

         AddMult_a_VVt(w, lshape, elmat);
      }
   }
protected:
   const IntegrationRule* GetDefaultIntegrationRule(
      const FiniteElement& trial_fe,
      const FiniteElement& test_fe,
      const ElementTransformation& trans) const override
   {
      int order = trans.OrderGrad(&trial_fe) + trans.Order() + test_fe.GetOrder();

      return &IntRules.Get(trial_fe.GetGeomType(), order+4);
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
   //    ch2n handle triangular, quadrilateral, tetrahedral, hexahedral, surface
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
   real_t h2 = 1.0/real_t(mesh->GetNE());

   // 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution (1 time by
   //    default, or specified on the command line with -rp). Once the parallel
   //    mesh is defined, the serial mesh ch2n be deleted.
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

   FiniteElementCollection *fec2 = new H1_FECollection(2*order, dim);
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
   ParGridFunction ch2(fespace2);
   ParGridFunction c_gf(fespace2);
   ch2 = 1.0;
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   ConstantCoefficient eps(1e-8);
   ConstantCoefficient h1_val(0.1);
   GridFunctionCoefficient c_cf(&ch2);
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


   ParBilinearForm *d = new ParBilinearForm(fespace);
   d->AddDomainIntegrator(new DiffusionIntegrator(one));
   d->Assemble();
   d->Finalize();

   ParBilinearForm *m2 = new ParBilinearForm(fespace2);
   m2->AddDomainIntegrator(new MassIntegrator(one));
   m2->AddDomainIntegrator(new DiffusionIntegrator(h1_val));
   m2->Assemble();
   m2->Finalize();

   HypreParMatrix *A = nullptr;
   HypreParMatrix *M = m->ParallelAssemble();
   HypreParMatrix *D = d->ParallelAssemble();
   HypreParMatrix *M2 = m2->ParallelAssemble();

   Solver * precond = NULL;
   HypreEuclid *ilu = new HypreEuclid(*M);
   HypreBoomerAMG  *amg = new HypreBoomerAMG(*M);
   precond = ilu;

   // Create nullspace
   MeanZeroSolver precond_proj(*precond, *M, *fespace);
   MeanZeroProjector  proj(*M, *fespace);

   // Loop
   Array<ParGridFunction*> ax;
   VisItDataCollection *visit_dc = nullptr;
   int nev = 1;
   Vector diff(ch2.Size());
   int ncluster0;
   real_t ev0;
   real_t Cinv = -1;
   for (int iter = 0; iter<=0; iter++)
   {
      if (myid == 0)
      {
         cout<<"\n\n----------------------------------------------\n";
         cout<<" iter  = "<<iter<<endl;
         cout<<"----------------------------------------------\n\n\n";
      }
      cout<<"Begin ch2 = "<<ch2.Min()<<" ::  "<<ch2.Max()<<endl;
      delete a;
      delete A;
      a = new ParBilinearForm(fespace);
      a->AddDomainIntegrator(new LapLapIntegrator(c_cfa));
      a->Assemble();
      a->Finalize();
      A = a->ParallelAssemble();

      // 8. Define and configure the LOBPCG eigensolver and the BoomerAMG
      //    preconditioner for A to be used within the solver. Set the matrices
      //    which define the generalized eigenproblem A x = lambda M x.
      HypreLOBPCG  lobpcg(MPI_COMM_WORLD);
      lobpcg.SetNumModes(nev);
      lobpcg.SetRandomSeed(seed);
      //lobpcg.SetPreconditioner(precond_proj);
      lobpcg.SetPreconditioner(*precond);
      lobpcg.SetMaxIter(1400);
      lobpcg.SetTol(1e-124);
      lobpcg.SetPrecondUsageMode(0);
      lobpcg.SetPrintLevel(1);
      lobpcg.SetMassMatrix(*A);
      lobpcg.SetOperator(*D);
      lobpcg.SetSubSpaceProjector(proj);

      // 9. Compute the eigenmodes and extract the array of eigenvalues. Define a
      //    parallel grid function to represent each of the eigenmodes returned by
      //    the solver.
      lobpcg.Solve();
      lobpcg.GetEigenvalues(eigenvalues);

      // Check for clustering of eigenvalues
      real_t cluster_tol = 5e-2;
      int ii;
      for (ii = 1; ii < nev; ii++)
      {
         if ((eigenvalues[ii]/eigenvalues[0]) > 1.0 + cluster_tol) { break; }
      }
     // MFEM_VERIFY(ii < nev, "Eigen space not large enough to get cluster");
      int ncluster = ii;
      //if (ncluster > 1)
      //ncluster = 4;

      cout<<"Cluster size = "<<ncluster<<"  <-- "<<ii<<endl;

      // Initialize if this is the first iteration
      if (iter == 0)
      {
         Cinv = h2/eigenvalues[0];
         ev0 = eigenvalues[0];
         ch2 = eigenvalues[0];
         diff = ch2;
         ncluster0 = ncluster;
   cout<<"ev 0 = "<<ev0<<" -- > Cinv = "<< Cinv <<endl;
         L2_ref = MatrixNorm(M2,ch2);
         L1_ref = ch2.ComputeL1Error(zero);

         cout<< "ITER "<< -1 <<" : ev = "<< 1
             <<"  "<<ch2.Min()<<" <= ch2 <= "<<ch2.Max()
             <<"  L2 = "<<L2_ref<<" <--- "<<L2_ref
             <<"  L1 = "<<L1_ref<<" <--- "<<L1_ref<<endl;

         if (visit)
         {
            ax.SetSize(ncluster);
            for (ii = 0; ii < ncluster; ii++)
            {
               ax[ii] = new ParGridFunction(fespace);
            }
            visit_dc = new VisItDataCollection ("InverseEstimate", pmesh);
            visit_dc->RegisterField("Ch^2", &ch2);
            visit_dc->RegisterField("C", &c_gf);
            for (ii = 0; ii < ax.Size(); ii++)
            {
               std::ostringstream oss;
               oss << "mode_" << std::setw(3) << std::setfill('0') << ii;
               visit_dc->RegisterField(oss.str(), ax[ii]);
            }
         }
      }
      else
      {
         ch2 *= eigenvalues[0];
      }

      // Print eigenmodes and coefficient - if requested
      if (visit_dc)
      {
         for (ii = 0; ii < ax.Size(); ii++)
         {
            (*ax[ii]) = lobpcg.GetEigenvector(ii);
         }
         c_gf = ch2;
         c_gf /= h2;

         visit_dc->SetCycle(iter);
         visit_dc->Save();
      }


      // Solve constrained minimization
      int nconstraint = ncluster;// + ncluster0;

      Array<Vector *> vv( nconstraint );
      ParLinearForm v(fespace2);
      v.AddDomainIntegrator(new DiffusionRQConstraintIntegrator(x));
      for (int i = 0; i < vv.Size() ; i++)
      {
         x = lobpcg.GetEigenvector(i);
         v = 0.0;
         v.Assemble();
         vv[i] = new Vector(v);
      }

      GlobalConstriant2 Mc(*M2, vv);
      Vector sol(ch2.Size() + nconstraint); sol = 0.0;
      Vector sol_r(sol.GetData(), ch2.Size());
      Vector rhs(ch2.Size() + nconstraint); rhs = 0.0;
      Vector rhs_r(rhs.GetData(), ch2.Size());

      M2->Mult(ch2, rhs_r);
      rhs *=-1;

      int max_iter = 200;
      real_t tol = 1e-8;
      GMRES(Mc,  sol, rhs, max_iter,200, tol, 0.0, 0);
      ch2.Print();
      sol_r.Print();

      ch2.Add(0.6, sol_r);

      // Check if constraints are satisfied
      for (int j = 0; j < vv.Size() ; j++)
      {
         cout<<"Constriant Check =" <<(*vv[j])*sol_r<<endl;
         delete vv[j];
      }

      // Check convergence
      diff -= ch2;
      real_t ch2_err = diff.Norml2();
      diff = ch2;

      // Report data
      real_t L2 = MatrixNorm(M2,ch2);
      real_t L1 = ch2.ComputeL1Error(zero);

      cout<< "ITER "<<iter <<" : ev = "<<eigenvalues[0]
          <<"  "<<ch2.Min()<<" <= ch2 <= "<<ch2.Max()
          <<"  L2 = "<<L2<<" <--- "<<L2_ref
          <<"  L1 = "<<L1<<" <--- "<<L1_ref
          <<"  Diff = "<<ch2_err<<" "<<ch2_err/ch2.Norml2()<<endl;

      if (ch2_err < 1e-4*ch2.Norml2()) { break; }
   }
   cout<<"ev 0 = "<<ev0<<" -- > Cinv = "<< Cinv <<endl;
 
   if (visit_dc)
   {
      c_gf = ch2;
      c_gf /= h2;
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
