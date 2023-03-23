//                          MFEM Example 1 - Waves Version
//
// Compile with: make waves
//
// Sample runs:  waves -m ../../data/periodic-square-nurbs.mesh -o 2 

// Description: TBD

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


/** After spatial discretization, the conduction model can be written as:
 *
 *     d^2u/dt^2 = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class WaveOperator represents the right-hand side of the above ODE.
 */
class WaveOperator : public SecondOrderTimeDependentOperator
{
protected:
   FiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   BilinearForm *T, *K, *M;

   CGSolver T_solver; // Implicit solver for T = M + fac0*K
   //HypreBoomerAMG
   GSSmoother T_prec;  // Preconditioner for the implicit solver

   ConstantCoefficient *g_inv, *cfac0;
   mutable Vector z; // auxiliary vector

public:
   WaveOperator(FiniteElementSpace &f, Array<int> &ess_bdr,double speed);

   /** Solve the Backward-Euler equation:
       d2udt2 = f(u + fac0*d2udt2,dudt + fac1*d2udt2, t),
       for the unknown d2udt2. */
   using SecondOrderTimeDependentOperator::ImplicitSolve;
   virtual void ImplicitSolve(const double fac0, const double fac1,
                              const Vector &u, const Vector &dudt, Vector &d2udt2);

   void ComputeEnergy(const Vector &u, const Vector &dudt);

   virtual ~WaveOperator();
};


WaveOperator::WaveOperator(FiniteElementSpace &f,
                           Array<int> &free_surf, double gravity)
   : SecondOrderTimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f),
     K(NULL), T(NULL), z(height)
{
   g_inv = new ConstantCoefficient(1.0/gravity);
   cfac0 = new ConstantCoefficient(-1.0);
   T = new BilinearForm(&fespace);
   T->AddDomainIntegrator(new DiffusionIntegrator(*cfac0));
   T->AddBoundaryIntegrator(new BoundaryMassIntegrator(*g_inv),free_surf);

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(1e-4);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(450);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);

   K = new BilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator());
   K->Assemble();
   K->Finalize();

   M = new BilinearForm(&fespace);
   M->AddBoundaryIntegrator(new BoundaryMassIntegrator(*g_inv),free_surf);
   M->Assemble();
   M->Finalize();
}

void WaveOperator::ImplicitSolve(const double fac0, const double fac1,
                                 const Vector &u, const Vector &dudt, Vector &d2udt2)
{
   // Solve for d2udt2 the equation:
   // K*fac0*d2udt2 + M*d2udt2 = -K*u
   if (cfac0->constant != fac0)
   {
      cfac0->constant = fac0;
      T->Assemble();
      T->Finalize();
      T_solver.SetOperator(T->SpMat());
   }
   K->Mult(u, z);
   z.Neg();
   T_solver.Mult(z, d2udt2);
}


void WaveOperator::ComputeEnergy(const Vector &u, const Vector &dudt)
{
   K->Mult(u, z);
   double Ekin = z*u;

   M->Mult(dudt, z);
   double Epot= z*dudt;

   cout<<Ekin<<" + "<<Epot<<" = "<<Ekin+Epot<<endl;
}

WaveOperator::~WaveOperator()
{
   delete T, K;
   delete g_inv,cfac0;
}

class Initial
{
public:

   double omega,k,a,h,g;

   Initial()
   {
      g = 9.81;
      k = 0.5*2*M_PI;
      a = 0.1;
      h = 2;
      ComputeOmega();
   };

   void ComputeOmega()
   {
      omega = sqrt(g*k*tanh(k*h));
   }

   double Solution(const Vector &x)
   {
      int dim = x.Size();
      double z = x[dim-1]+1;
      return (omega/k)*a*(cosh(k*(z+h))/sinh(k*h))*cos(k*x[0]);
   }

   double Rate(const Vector &x)
   {
      int dim = x.Size();
      double z = x[dim-1]+1;
      return omega*(omega/k)*a*(cosh(k*(z+h))/sinh(k*h))*sin(k*x[0]);
   }

   std::function<double(const Vector &)> GetSolFun() { return std::bind( &Initial::Solution, this, std::placeholders::_1);}
   std::function<double(const Vector &)> GetRateFun() { return std::bind( &Initial::Rate, this, std::placeholders::_1);}

   static double Solution2(const Vector &x)
   {
      double  k = 0.5*2*M_PI;
      double  a = 1;
      double h = 1;
      double  g = 9.81;
      double  omega = sqrt(g*tanh(k*h));
      int dim = x.Size();
      double z = x[dim-1];
      return (omega/k)*a*(cosh(k*(z+h))/sinh(k*h))*cos(k*x[0]);
   }

   static double Rate2(const Vector &x)
   {
      double  k = 0.5*2*M_PI;
      double  a = 1;
      double h = 2;
      double  g = 9.81;
      double  omega = sqrt(g*tanh(k*h));
      int dim = x.Size();
      double z = x[dim-1]+1.0;
      return omega*(omega/k)*a*(cosh(k*(z+h))/sinh(k*h))*sin(k*x[0]);
   }

};

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../../data/square-nurbs.mesh";
   const char *per_file  = "none";
   int ref_levels = -1;
   Array<int> master(0);
   Array<int> slave(0);
   Array<int> free_surf(0);

   Array<int> order(1);
   order[0] = 1;
   int ode_solver_type = 10;
   double t_final = 0.5;
   double dt = 1.0e-2;
   int vis_steps = 5;
   Initial initial;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly, -1 for auto.");
   args.AddOption(&per_file, "-p", "--per",
                  "Periodic BCS file.");
   args.AddOption(&master, "-pm", "--master",
                  "Master boundaries for periodic BCs");
   args.AddOption(&slave, "-ps", "--slave",
                  "Slave boundaries for periodic BCs");
   args.AddOption(&free_surf, "-fs", "--free-surf",
                  "Boundaries for free-surface BCs");

   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: [0--10] - GeneralizedAlpha(0.1 * s),\n\t"
                  "\t   11 - Average Acceleration, 12 - Linear Acceleration\n"
                  "\t   13 - CentralDifference, 14 - FoxGoodwin"); 
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");

   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file and refine the mesh.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }
   mesh->PrintInfo();
   int dim = mesh->Dimension();

   // 3. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   NURBSExtension *NURBSext = NULL;
   int own_fec = 0;

   if (mesh->NURBSext)
   {
      fec = new NURBSFECollection(order[0]);
      own_fec = 1;

      int nkv = mesh->NURBSext->GetNKV();
      if (order.Size() == 1)
      {
         int tmp = order[0];
         order.SetSize(nkv);
         order = tmp;
      }

      if (order.Size() != nkv ) { mfem_error("Wrong number of orders set."); }
      NURBSext = new NURBSExtension(mesh->NURBSext, order);

      // Read periodic BCs from file
      std::ifstream in;
      in.open(per_file, std::ifstream::in);
      if (in.is_open())
      {
         int psize;
         in >> psize;
         master.SetSize(psize);
         slave.SetSize(psize);
         master.Load(in, psize);
         slave.Load(in, psize);
         in.close();
      }
      master.Print();
      slave.Print();
      NURBSext->ConnectBoundaries(master,slave);
   }
   else if (order[0] == -1) // Isoparametric
   {
      if (mesh->GetNodes())
      {
         fec = mesh->GetNodes()->OwnFEC();
         own_fec = 0;
         cout << "Using isoparametric FEs: " << fec->Name() << endl;
      }
      else
      {
         cout <<"Mesh does not have FEs --> Assume order 1.\n";
         fec = new H1_FECollection(1, dim);
         own_fec = 1;
      }
   }
   else
   {
      if (order.Size() > 1) { cout <<"Wrong number of orders set, needs one.\n"; }
      fec = new H1_FECollection(abs(order[0]), dim);
      own_fec = 1;
   }

   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, NURBSext, fec);
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;

   // 4. Define the ODE solver used for time integration.
   SecondOrderODESolver *ode_solver;
   switch (ode_solver_type)
   {
      // Implicit methods
      case 0: ode_solver = new GeneralizedAlpha2Solver(0.0); break;
      case 1: ode_solver = new GeneralizedAlpha2Solver(0.1); break;
      case 2: ode_solver = new GeneralizedAlpha2Solver(0.2); break;
      case 3: ode_solver = new GeneralizedAlpha2Solver(0.3); break;
      case 4: ode_solver = new GeneralizedAlpha2Solver(0.4); break;
      case 5: ode_solver = new GeneralizedAlpha2Solver(0.5); break;
      case 6: ode_solver = new GeneralizedAlpha2Solver(0.6); break;
      case 7: ode_solver = new GeneralizedAlpha2Solver(0.7); break;
      case 8: ode_solver = new GeneralizedAlpha2Solver(0.8); break;
      case 9: ode_solver = new GeneralizedAlpha2Solver(0.9); break;
      case 10: ode_solver = new GeneralizedAlpha2Solver(1.0); break;

      case 11: ode_solver = new AverageAccelerationSolver(); break;
      case 12: ode_solver = new LinearAccelerationSolver(); break;
      case 13: ode_solver = new CentralDifferenceSolver(); break;
      case 14: ode_solver = new FoxGoodwinSolver(); break;

      default:
         cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         delete mesh, fespace;
         return 3;
   }

   // 5. Set the initial conditions for u. All boundaries are considered
   //    natural.
   GridFunction u_gf(fespace);
   GridFunction dudt_gf(fespace);

   if (mesh->NURBSext)
   {
      Vector coord(dim);
      for (int i = 0; i < mesh->GetNV(); i++)
      {
         mesh->GetNode(i, &coord[0]);
         u_gf[i] = initial.Solution(coord);
         dudt_gf[i] = initial.Rate(coord);
      }

   }
   else
   {
      FunctionCoefficient u_0(initial.GetSolFun());
      FunctionCoefficient dudt_0(initial.GetRateFun());

      u_gf.ProjectCoefficient(u_0);
      dudt_gf.ProjectCoefficient(dudt_0);
   }

   Vector u, dudt;
   u_gf.GetTrueDofs(u);
   dudt_gf.GetTrueDofs(dudt);

   // 6. Initialize the wave operator and the visualization.
   if (mesh->bdr_attributes.Size() == 0)
   {
      cout << "No boundaries in the mesh\n";
      delete mesh, fespace;
      return 4;
   }

   Array<int> free_surf_idx(mesh->bdr_attributes.Max());
   free_surf_idx = 0;
   for (int i = 0; i < free_surf.Size(); i++)
   {
      for (int ii = 0; ii < mesh->bdr_attributes.Size(); ii++)
      {
         if (free_surf[i] == mesh->bdr_attributes[ii]) free_surf_idx[ii] = 1;
      }
   }

   free_surf_idx.Print();
   WaveOperator oper(*fespace, free_surf_idx, 9.81);

   u_gf.SetFromTrueDofs(u);
   dudt_gf.SetFromTrueDofs(dudt);

   VisItDataCollection visit_dc("solution/wave", mesh);
   visit_dc.RegisterField("solution", &u_gf);
   visit_dc.RegisterField("rate", &dudt_gf);
   visit_dc.SetCycle(0);
   visit_dc.SetTime(0.0);
   visit_dc.Save();

   // 7. Perform time-integration (looping over the time iterations, ti, with a
   //    time-step dt).
   ode_solver->Init(oper);
   double t = 0.0;

   bool last_step = false;
   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= t_final - dt/2)
      {
         last_step = true;
      }

      ode_solver->Step(u, dudt, t, dt);
      oper.ComputeEnergy(u, dudt);

      if (last_step || (ti % vis_steps) == 0)
      {
         cout << "step " << ti << ", t = " << t << endl;

         u_gf.SetFromTrueDofs(u);
         dudt_gf.SetFromTrueDofs(dudt);

         visit_dc.SetCycle(ti);
         visit_dc.SetTime(t);
         visit_dc.Save();
      }
   }

   // 8. Free the used memory.
   delete fespace;
   if (own_fec) { delete fec; }
   delete mesh;

   return 0;
}
