// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
//   _____  ______      ____  __  _____
//   |  __ \|  _ \ \    / /  \/  |/ ____|
//   | |__) | |_) \ \  / /| \  / | (___
//   |  _  /|  _ < \ \/ / | |\/| |\___ \
//   | | \ \| |_) | \  /  | |  | |____) |
//   |_|  \_\____/   \/   |_|  |_|_____/
//
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------
#include <sys/stat.h>
#include "mfem.hpp"
#include "../util/coefficients.hpp"
#include "solver.hpp"

#include "dist_solver.hpp"
#include "formulation.hpp"
#include "integrator.hpp"

using namespace std;
using namespace mfem;

extern void printInfo();
extern void line(int len);


// Routine for checking duplicity of boundary conditions
void CheckBoundaries2(Array<bool> &bnd_flags,
                      Array<int> &bc_bnds)
{
   int amax = bnd_flags.Size();
   // Check strong boundaries
   for (int b = 0; b < bc_bnds.Size(); b++)
   {
      int bnd = bc_bnds[b];
      if ( bnd < 0 || bnd > amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flags[bnd])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flags[bnd] = true;
   }
}

// Routine for pretty printing boundary forces to screen
void PrintForce2(const Array<int>& bdr_attributes,
                 const DenseMatrix& bdrForce)
{
   // Print line lambda function
   auto pline = [](int len)
   {
      cout<<" +";
      for (int b=0; b<len; ++b) { cout<<"-"; }
      cout<<"+\n";
   };

   // Print boundary header
   int nbdr = bdr_attributes.Size();
   cout<<"\n";
   pline(10+13*nbdr);
   cout<<" | Boundary | ";
   for (int b=0; b<nbdr; ++b)
   {
      cout<<std::setw(10)<<bdr_attributes[b]<<" | ";
   }
   cout<<"\n";
   pline(10+13*nbdr);

   // Print actual forces
   char dimName[] = "xyz";
   for (int v=0; v<bdrForce.Width(); ++v)
   {
      cout<<" | Force "<<dimName[v]<<"  | ";
      for (int b=0; b<nbdr; ++b)
      {
         int bnd = bdr_attributes[b];
         cout<<std::defaultfloat<<std::setprecision(4)<<std::setw(10);
         cout<<bdrForce(bnd-1,v)<<" | ";
      }
      cout<<"\n";
   }
   pline(10+13*nbdr);
   cout<<"\n"<<std::flush;
}

// Helper class for writting global data to file
class OutputData
{
private:
   std::ofstream os;
   Array<int> bdr_attr;
   bool print;

public:
   // Constructor
   OutputData(bool prt, int dim, Array<int>& bdr_attr_, int si)
      : print(prt), bdr_attr(bdr_attr_)
   {
      if (!print) { return; }

      std::ostringstream filename;
      filename << "output_"<<std::setw(6)<<setfill('0')<<si<< ".dat";
      os.open(filename.str().c_str());

      // Header
      char dimName[] = "xyz";
      int i = 9;
      os <<"# 1: step"<<"\t"<<"2: time"<<"\t"<<"3: dt"<<"\t"
         <<"4: cfl"<<"\t"<<"5: outflow"<<"\t"
         <<"6: Ekin"<<"\t"<<"7: Epot"<<"\t"<<"8: Visc disp"<<"\t";
      for (int b=0; b<bdr_attr.Size(); ++b)
      {
         int bnd = bdr_attr[b];
         for (int v=0; v<dim; ++v)
         {
            std::ostringstream forcename;
            forcename <<i++<<": F"<<dimName[v]<<"_"<<bnd;
            os<<forcename.str()<<"\t";
         }
      }
      os<<endl;
   };

   // Print to file
   void Print(int si, real_t t, real_t dt, real_t cfl,
              real_t outflow, Vector &energy, DenseMatrix& bdrForce)
   {
      if (!print) { return; }

      int nbdr = bdr_attr.Size();
      os<<std::setw(10);
      os<<si<<"\t"<<t<<"\t"<<dt<<"\t"<<cfl<<"\t"<<outflow<<"\t";
      os<<energy[0]<<"\t"<<energy[1]<<"\t"<<energy[2]<<"\t";
      for (int b=0; b<nbdr; ++b)
      {
         int bnd = bdr_attr[b];
         for (int v=0; v<bdrForce.Width(); ++v)
         {
            os<<bdrForce(bnd-1,v)<<"\t";
         }
      }
      os<<"\n"<< std::flush;
   };

   // Destructor
   ~OutputData()
   {
      os.close();
   };

};

int main(int argc, char *argv[])
{
   // Initialize MPI and HYPRE and print info
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();
   printInfo();

   // Parse command-line options.
   OptionsParser args(argc, argv);

   // Mesh and discretization parameters
   const char *mesh_file = "../../mfem/data/inline-quad.mesh";
   const char *ref_file  = "";
   int order = 1;
   int ref_levels = 0;
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order isoparametric space.");

   // Problem parameters
   Array<int> strong_bdr;
   Array<int> weak_bdr;
   Array<int> normal_bdr;

   Array<int> outflow_bdr;
   Array<int> suction_bdr;
   Array<int> blowing_bdr;
   Array<int> master_bdr;
   Array<int> slave_bdr;

   const char *lib_file = "libfun.so";
   real_t rho_param = 1.204;
   real_t mu_param = 1.825e-5;

   args.AddOption(&strong_bdr, "-sbc", "--strong-bdr",
                  "Boundaries where Dirichelet BCs are enforced strongly.");
   args.AddOption(&weak_bdr, "-wbc", "--weak-bdr",
                  "Boundaries where Dirichelet BCs are enforced weakly.");
   args.AddOption(&normal_bdr, "-nbc", "--normal-bdr",
                  "Boundaries where Normal Dirichelet BCs are enforced weakly.");
   args.AddOption(&outflow_bdr, "-out", "--outflow-bdr",
                  "Outflow boundaries.");
   args.AddOption(&suction_bdr, "-suc", "--suction-bdr",
                  "Suction boundaries.");
   args.AddOption(&blowing_bdr, "-blow", "--blowing-bdr",
                  "Blowing boundaries.");
   args.AddOption(&master_bdr, "-mbc", "--master-bdr",
                  "Periodic master boundaries.");
   args.AddOption(&slave_bdr, "-sbc", "--slave-bdr",
                  "Periodic slave boundaries.");
   args.AddOption(&lib_file, "-l", "--lib",
                  "Library file for case specific function definitions:\n\t"
                  " - Initial condition\n\t"
                  " - Boundary condition\n\t"
                  " - Forcing\n\t"
                  " - Diffusion\n\t");
   args.AddOption(&mu_param, "-mu", "--dyn-visc",
                  "Sets the dynamic diffusion parameters, should be positive.");
   args.AddOption(&rho_param, "-rho", "--density",
                  "Sets the density parameters, should be positive.");

   // Time stepping params
   int ode_solver_type = 35;
   real_t t_final = 10.0;
   real_t dt = 0.01;
   real_t dt_max = 1.0;
   real_t dt_min = 0.0001;

   real_t cfl_target = 2.0;
   real_t dt_gain = -1.0;

   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  ODESolver::ImplicitTypes.c_str());
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--dt",
                  "Time step.");
   args.AddOption(&dt_min, "-dtmn", "--dt-min",
                  "Minimum time step size.");
   args.AddOption(&dt_max, "-dtmx", "--dt-max",
                  "Maximum time step size.");
   args.AddOption(&cfl_target, "-cfl", "--cfl-target",
                  "CFL target.");
   args.AddOption(&dt_gain, "-dtg", "--dt-gain",
                  "Gain coefficient for time step adjustment.");

   // Navier-Stokes formulation parameters
   real_t NavSto_kdc0 = 0.0;
   real_t NavSto_kdc1 = 0.1;

   args.AddOption(&NavSto_kdc0, "-nsk0", "--navsto-kdc0",
                  "Inconsistent diffusion parameter of the NS formulation.");

   args.AddOption(&NavSto_kdc1, "-nsk1", "--navsto-kdc1",
                  "Consistent diffusion parameter of tha NS formulation.");


   // RBVMS Solver parameters
   double GMRES_RelTol = 1e-3;
   int    GMRES_MaxIter = 500;
   double Newton_RelTol = 1e-3;
   int    Newton_MaxIter = 10;

   args.AddOption(&GMRES_RelTol, "-lt", "--linear-tolerance",
                  "Relative tolerance for the GMRES solver.");
   args.AddOption(&GMRES_MaxIter, "-li", "--linear-itermax",
                  "Maximum iteration count for the GMRES solver.");
   args.AddOption(&Newton_RelTol, "-nt", "--newton-tolerance",
                  "Relative tolerance for the Newton solver.");
   args.AddOption(&Newton_MaxIter, "-ni", "--newton-itermax",
                  "Maximum iteration count for the Newton solver.");

   // Redistancing formulation parameters
   real_t Redist_Lambda = 1.0;
   real_t Redist_kdc0 = 0.01;
   real_t Redist_kdc1 = 0.25;

   args.AddOption(&Redist_Lambda, "-rl", "--redist-penalty",
                  "Interface pinning penalty of the Redistancing formulation.");

   args.AddOption(&Redist_kdc0, "-rk0", "--redist-kdc0",
                  "Inconsistent diffusion parameter of Redistancing formulation.");

   args.AddOption(&Redist_kdc1, "-rk1", "--redist-kdc1",
                  "Consistent diffusion parameter of Redistancing formulation.");

   // Redistancing solver parameters
   real_t Redist_GMRES_RelTol = 1e-3;
   int Redist_GMRES_MaxIter = 100;
   real_t Redist_Newton_RelTol = 1e-3;
   int Redist_Newton_MaxIter = 5;

   args.AddOption(&Redist_GMRES_RelTol, "-rlt", "--redist-linear-tolerance",
                  "Relative tolerance for the Redistancing GMRES solver.");
   args.AddOption(&Redist_GMRES_MaxIter, "-rli", "--redist-linear-itermax",
                  "Maximum iteration count for the Redistancing GMRES solver.");
   args.AddOption(&Redist_Newton_RelTol, "-rnt", "--redist-newton-tolerance",
                  "Relative tolerance for the Redistancing Newton solver.");
   args.AddOption(&Redist_Newton_MaxIter, "-rni", "--redist-newton-itermax",
                  "Maximum iteration count for the Redistancing Newton solver.");

   // Volume conservation solver parameters
   int VolCons_MaxIter = 10;
   real_t VolCons_RelTol = 1e-12;
   real_t VolCons_JacEps = 1e-6;

   args.AddOption(&VolCons_MaxIter, "-vci", "--volcons-newton-itermax",
                  "Maximum iteration count for the volume conservation solver.");
   args.AddOption(&VolCons_RelTol, "-vct", "--volcons-newton-tolerance",
                  "Relative tolerance for the volume conservation solver.");
   args.AddOption(&VolCons_JacEps, "-vce", "--volcons-jac-eps",
                  "Finite-difference pertubation for obtaining the jacobian\n\t"
                  "for the volume conservation solver.");

   // Solution input/output params
   bool restart = false;
   int restart_interval = -1;
   real_t dt_vis = 10*dt;
   const char *vis_dir = "solution";
   args.AddOption(&restart, "-rs", "--restart", "-f", "--fresh",
                  "Restart from solution.");
   args.AddOption(&restart_interval, "-ri", "--restart-interval",
                  "Interval between archieved time steps.\n\t"
                  "For negative values output is skipped.");
   args.AddOption(&dt_vis, "-dtv", "--dt_vis",
                  "Time interval between visualization points.");
   args.AddOption(&vis_dir, "-vd", "--vis-dir",
                  "Directory for visualization files.\n\t");

   // Parse parameters
   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   // Read the mesh from the given mesh file.
   Mesh mesh(mesh_file);//, 1, 1);
   int dim = mesh.Dimension();

   // Refine mesh
   {
      if (mesh.NURBSext && (strlen(ref_file) != 0))
      {
         mesh.RefineNURBSFromFile(ref_file);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
      if (Mpi::Root()) { mesh.PrintInfo(); }
   }

   // Partition mesh
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();

   // Boundary conditions
   if (Mpi::Root())
   {
      if (strong_bdr.Size()>0) {cout<<"Strong  = "; strong_bdr.Print();}
      if (weak_bdr.Size()>0) {cout<<"Weak    = "; weak_bdr.Print();}
      if (normal_bdr.Size()>0) { cout<<"Normal  = "; normal_bdr.Print();}
      if (outflow_bdr.Size()>0) { cout<<"Outflow = "; outflow_bdr.Print() ;}
      if (suction_bdr.Size()>0) { cout<<"Suction = "; suction_bdr.Print() ;}
      if (blowing_bdr.Size()>0) { cout<<"Blowing = "; blowing_bdr.Print() ;}
      if (master_bdr.Size()>0) {cout<<"Periodic (master) = "; master_bdr.Print();}
      if (slave_bdr.Size()>0) {cout<<"Periodic (slave)  = "; slave_bdr.Print();}
   }

   Array<bool> bnd_flag(pmesh.bdr_attributes.Max()+1);
   bnd_flag = true;
   for (int b = 0; b < pmesh.bdr_attributes.Size(); b++)
   {
      bnd_flag[pmesh.bdr_attributes[b]] = false;
   }
   CheckBoundaries2(bnd_flag, strong_bdr);
   CheckBoundaries2(bnd_flag, weak_bdr);
   CheckBoundaries2(bnd_flag, normal_bdr);
   CheckBoundaries2(bnd_flag, outflow_bdr);
   CheckBoundaries2(bnd_flag, suction_bdr);
   CheckBoundaries2(bnd_flag, blowing_bdr);
   CheckBoundaries2(bnd_flag, master_bdr);
   CheckBoundaries2(bnd_flag, slave_bdr);

   MFEM_VERIFY(master_bdr.Size() == master_bdr.Size(),
               "Master-slave count do not match.");
   for (int b = 0; b < bnd_flag.Size(); b++)
   {
      MFEM_VERIFY(bnd_flag[b],
                  "Not all boundaries have a boundary condition set.");
   }

   // Select the time integrator
   unique_ptr<ODESolver> ode_solver = ODESolver::Select(ode_solver_type);
   ODESolverWithStates*  ode_solver_ws = dynamic_cast<ODESolverWithStates*>
                                         (ode_solver.get());
   int nstate = ode_solver->GetStateSize();

   if (nstate > 1 && ( restart || restart_interval > 0 ))
   {
      mfem_error("RBVMS restart not available for this time integrator \n"
                 "Time integrator can have a maximum of one statevector.");
   }

   // Define a finite element space on the mesh.
   Array<FiniteElementCollection *> fecs(3);
   if (pmesh.NURBSext)
   {
      fecs[0] = new NURBSFECollection(order);
      fecs[1] = new NURBSFECollection(order);
      fecs[2] = new NURBSFECollection(order);
   }
   else
   {
      fecs[0] = new H1_FECollection(order, dim);
      fecs[1] = new H1_FECollection(order, dim);
      fecs[2] = new H1_FECollection(order, dim);
   }


   Array<ParFiniteElementSpace *> spaces(3);
   spaces[0] = new ParFiniteElementSpace(&pmesh, fecs[0], dim,
                                         Ordering::byNODES  //, Ordering::byVDIM);
                                        );// ,master_bdr, slave_bdr);
   spaces[1] = new ParFiniteElementSpace(&pmesh, fecs[1], 1, Ordering::byNODES
                                        );//  ,master_bdr, slave_bdr);

   spaces[2] = new ParFiniteElementSpace(&pmesh, fecs[2], 1, Ordering::byNODES
                                        );//  ,master_bdr, slave_bdr);

   // Report the degree of freedoms used
   {
      Array<int> tdof(num_procs),udof(num_procs),
            ldof(num_procs),pdof(num_procs);
      tdof = 0;
      tdof[myid] = spaces[0]->TrueVSize();
      MPI_Reduce(tdof.GetData(), udof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      tdof = 0;
      tdof[myid] = spaces[1]->TrueVSize();
      MPI_Reduce(tdof.GetData(), pdof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      tdof = 0;
      tdof[myid] = spaces[2]->TrueVSize();
      MPI_Reduce(tdof.GetData(), ldof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      int udof_t = 0;//spaces[0]->GlobalTrueVSize();
      int pdof_t = 0;//spaces[1]->GlobalTrueVSize();
      int ldof_t = 0;//spaces[2]->GlobalTrueVSize();
      if (Mpi::Root())
      {
         mfem::out << "Number of finite element unknowns:\n";
         mfem::out << "\tVelocity  = "<< udof_t << endl;
         mfem::out << "\tPressure  = "<< pdof_t << endl;
         mfem::out << "\tLevel-set = "<< ldof_t << endl;
         mfem::out << "Number of finite element unknowns per partition:\n";
         mfem::out <<  "\tVelocity  = "; udof.Print(mfem::out, num_procs);
         mfem::out <<  "\tPressure  = "; pdof.Print(mfem::out, num_procs);
         mfem::out <<  "\tLevel-set = "; ldof.Print(mfem::out, num_procs);
      }
   }

   // Get vector offsets
   Array<int> bOffsets(4);
   bOffsets[0] = 0;
   bOffsets[1] = spaces[0]->TrueVSize();
   bOffsets[2] = spaces[1]->TrueVSize();
   bOffsets[3] = spaces[2]->TrueVSize();
   bOffsets.PartialSum();

   // Define the solution vector, grid function and output
   BlockVector xp(bOffsets);
   BlockVector dxp(bOffsets);
   BlockVector xp0(bOffsets);
   BlockVector xpi(bOffsets);

   // Define the gridfunctions
   ParGridFunction x_u(spaces[0]);
   ParGridFunction x_p(spaces[1]);
   ParGridFunction x_phi(spaces[2]);
   ParGridFunction x_ref(spaces[2]);

   Array<ParGridFunction*> dx_u(nstate);
   Array<ParGridFunction*> dx_p(nstate);
   Array<ParGridFunction*> dx_phi(nstate);

   for (int i = 0; i < nstate; i++)
   {
      dx_u[i] = new ParGridFunction(spaces[0]);
      dx_p[i] = new ParGridFunction(spaces[1]);
      dx_phi[i] = new ParGridFunction(spaces[2]);
   }

   // Define the visualisation output
   VisItDataCollection vdc("step", &pmesh);
   vdc.SetPrefixPath(vis_dir);
   vdc.RegisterField("u", &x_u);
   vdc.RegisterField("p", &x_p);
   vdc.RegisterField("phi", &x_phi);

   // Get the start vector(s) from file -- or from function
   real_t t;
   int si, ri, vi;
   struct stat info;
   if (restart && stat("restart/step.dat", &info) == 0)
   {
      // Read
      if (Mpi::Root())
      {
         real_t dtr;
         std::ifstream in("restart/step.dat", std::ifstream::in);
         in>>t>>si>>ri>>vi;
         in>>dtr;
         in.close();
         cout<<"Restarting from step "<<ri-1<<endl;
         if (dt_gain > 0) { dt = dtr; }
      }
      // Synchronize
      MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&si, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ri, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&vi, 1, MPI_INT, 0, MPI_COMM_WORLD);

      // Open data files
      VisItDataCollection rdc("step", &pmesh);
      rdc.SetPrefixPath("restart");
      rdc.SetPrecision(18);

      rdc.Load(ri-1);

      x_u = *rdc.GetField("u");
      x_p = *rdc.GetField("p");
      x_phi = *rdc.GetField("phi");

      x_u.GetTrueDofs(xp.GetBlock(0));
      x_p.GetTrueDofs(xp.GetBlock(1));
      x_phi.GetTrueDofs(xp.GetBlock(2));

      if (nstate == 1 && rdc.GetField("du") &&
          rdc.GetField("dp") && rdc.GetField("dphi"))
      {
         *dx_u[0] = *rdc.GetField("du");
         *dx_p[0] = *rdc.GetField("dp");
         *dx_phi[0] = *rdc.GetField("dphi");

         dx_u[0]->GetTrueDofs(dxp.GetBlock(0));
         dx_p[0]->GetTrueDofs(dxp.GetBlock(1));
         dx_phi[2]->GetTrueDofs(dxp.GetBlock(2));

         ode_solver_ws->GetState().Append(dxp);
      }
   }
   else
   {
      // Define initial condition from file
      t = 0.0; si = 0; ri = 1; vi = 1;
      LibVectorCoefficient sol_u(dim, lib_file, "sol_u");
      LibCoefficient sol_p(lib_file, "sol_p");
      LibCoefficient sol_phi(lib_file, "sol_phi");
      sol_u.SetTime(-1.0);
      sol_p.SetTime(-1.0);
      sol_phi.SetTime(-1.0);
      x_u.ProjectCoefficient(sol_u);
      x_p.ProjectCoefficient(sol_p);
      x_phi.ProjectCoefficient(sol_phi);

      x_u.GetTrueDofs(xp.GetBlock(0));
      x_p.GetTrueDofs(xp.GetBlock(1));
      x_phi.GetTrueDofs(xp.GetBlock(2));

      // Visualize initial condition
      vdc.SetCycle(0);
      vdc.SetTime(0.0);
      vdc.Save();
   }


   // Define the restart output
   VisItDataCollection rdc("step", &pmesh);
   rdc.SetPrefixPath("restart");
   rdc.SetPrecision(18);

   rdc.RegisterField("u", &x_u);
   rdc.RegisterField("p", &x_p);
   rdc.RegisterField("phi", &x_phi);
   if (nstate == 1)
   {
      rdc.RegisterField("du", dx_u[0]);
      rdc.RegisterField("dp", dx_p[0]);
      rdc.RegisterField("dphi", dx_phi[0]);
   }

   // Set up the preconditioner
   RBVMS::JacobianPreconditioner jac_prec(bOffsets);

   Solver* pc_mom  = nullptr;
   Solver* pc_cont = nullptr;
   Solver* pc_ls   = nullptr;

   HypreILU* ilu_mom = new HypreILU();
   HypreILU* ilu_cont = new HypreILU();
   HypreILU* ilu_ls = new HypreILU();

   // HypreSmoother* ilu_mom = new HypreSmoother();
   // HypreSmoother* ilu_cont = new HypreSmoother();
   // HypreSmoother* ilu_ls = new HypreSmoother();

   pc_mom  = ilu_mom;
   pc_cont = ilu_cont;
   pc_ls   = ilu_ls;

   jac_prec.SetPreconditioner(0, pc_mom);
   jac_prec.SetPreconditioner(1, pc_cont);
   jac_prec.SetPreconditioner(2, pc_ls);

   // Set up the Jacobian solver
   RBVMS::GeneralResidualMonitor j_monitor("\t\tFGMRES", 10);
   FGMRESSolver gmres(MPI_COMM_WORLD);
   gmres.iterative_mode = false;
   gmres.SetRelTol(GMRES_RelTol);
   gmres.SetMaxIter(GMRES_MaxIter);
   gmres.SetPrintLevel(-1);
   gmres.SetMonitor(j_monitor);
   gmres.SetPreconditioner(jac_prec);

   // Set up the Newton solver
   RBVMS::NewtonSystemSolver newton_solver(MPI_COMM_WORLD,bOffsets);
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(1);
   newton_solver.SetRelTol(Newton_RelTol);
   newton_solver.SetMaxIter(Newton_MaxIter);
   newton_solver.SetSolver(gmres);

   // Define the physical parameters
   LibCoefficient rho(lib_file, "rho", false, rho_param);
   LibCoefficient mu(lib_file, "mu", false, mu_param);
   LibVectorCoefficient sol_u(dim, lib_file, "sol_u");
   LibCoefficient sol_phi(lib_file, "sol_phi");
   LibVectorCoefficient force(dim, lib_file, "force");
   LibCoefficient suction(lib_file, "suction", false, 0.0);
   LibCoefficient blowing(lib_file, "blowing", false, 0.0);

   // Define weak form and evolution
   RBVMS::IncNavStoIntegrator integrator(rho, mu, force,
                                         sol_u, sol_phi,
                                         suction, blowing);
   RBVMS::NavStoLSForm form(spaces, integrator);
   RBVMS::Evolution evo(form, newton_solver);
   ode_solver->Init(evo);

   // Set boundaries in the weakform
   form.SetStrongBC (strong_bdr);
   form.SetWeakBC   (weak_bdr);
   form.SetNormalBC (normal_bdr);
   form.SetOutflowBC(outflow_bdr);
   form.SetSuctionBC(suction_bdr);
   form.SetBlowingBC(blowing_bdr);
   form.SetConvectionOnly();

   form.SetInconsistentDC(NavSto_kdc0);
   form.SetConsistentDC(NavSto_kdc1);

   ConvectionDistanceSolver dist_solver(*spaces[2]);
   dist_solver.SetPenalty(Redist_Lambda);
   dist_solver.SetInconsistentDC(Redist_kdc0);
   dist_solver.SetConsistentDC(Redist_kdc1);

   //   RBVMS::GeneralResidualMonitor dist_monitor(MPI_COMM_WORLD,
   //                                              " - Redistance",
   //                                              1);
   //   dist_solver.SetNonlinearMonitor(dist_monitor);
   dist_solver.SetLinearRelTol(Redist_GMRES_RelTol);
   dist_solver.SetLinearMaxIter(Redist_GMRES_MaxIter);
   dist_solver.SetNonlinearRelTol(Redist_Newton_RelTol);
   dist_solver.SetNonlinearMaxIter(Redist_Newton_MaxIter);

   dist_solver.SetVolumeConservationMaxIter(VolCons_MaxIter);
   dist_solver.SetVolumeConservationRelTol(VolCons_RelTol);
   dist_solver.SetVolumeConservationJacEps(VolCons_JacEps);

   // Open output file
   OutputData output(Mpi::Root(), pmesh.Dimension(), pmesh.bdr_attributes, si);

   // Loop till final time reached
   while (t < t_final)
   {
      // Print header
      if (Mpi::Root())
      {
         line(80);
         cout<<std::defaultfloat<<std::setprecision(4);
         cout<<" step = " << si << endl;
         cout<<"   dt = " << dt << endl;
         cout<<std::defaultfloat<<std::setprecision(6);;
         cout<<" time = [" << t << ", " << t+dt <<"]"<< endl;
         cout<<std::defaultfloat<<std::setprecision(4);
         line(80);
      }

      // Actual time step
      xp0 = xp;
      ode_solver->Step(xp, t, dt);

      x_ref.Distribute(xp.GetBlock(2));
      x_phi = x_ref;
      dist_solver.ComputeScalarDistance(x_ref, x_phi);

      x_ref.Distribute(xp0.GetBlock(2));
      dist_solver.CorrectVolume(x_ref, x_phi);
      x_phi.GetTrueDofs(xp.GetBlock(2));

      si++;

      // Write visualization files
      while (t >= dt_vis*vi)
      {
         // Interpolate solution
         real_t fac = (t-dt_vis*vi)/dt;

         // Report to screen
         if (Mpi::Root())
         {
            line(80);
            cout << "Visit output: " <<vi << endl;
            cout << "        Time: " <<t-dt<<" "<<t-fac*dt<<" "<<t<<endl;
            line(80);
         }

         // Copy solution in grid functions
         add (fac, xp0.GetBlock(0), (1.0-fac), xp.GetBlock(0), xpi.GetBlock(0));
         x_u.Distribute(xpi.GetBlock(0));

         add (-1.0/dt, xp0.GetBlock(1), 1.0/dt, xp.GetBlock(1), xpi.GetBlock(1));
         x_p.Distribute(xpi.GetBlock(1));

         add (fac, xp0.GetBlock(2), (1.0-fac), xp.GetBlock(2), xpi.GetBlock(2));
         x_phi.Distribute(xpi.GetBlock(2));

         // Actually write to file
         vdc.SetCycle(vi);
         vdc.SetTime(dt_vis*vi);
         vdc.Save();
         vi++;
      }

      // Change time step
      real_t dt0 = dt;
      if (dt_gain > 0)
      {
         real_t cfl = form.GetCFL();
         dt *= pow(cfl_target/cfl, dt_gain);
         dt = min(dt, dt_max);
         dt = max(dt, dt_min);
      }

      // Write restart files
      if (restart_interval > 0 && si%restart_interval == 0)
      {
         // Report to screen
         if (Mpi::Root())
         {
            line(80);
            cout << "Restart output:" << ri << endl;
            line(80);
         }

         // Copy solution in grid functions
         x_u.Distribute(xp.GetBlock(0));
         x_p.Distribute(xp.GetBlock(1));
         x_phi.Distribute(xp.GetBlock(2));
         if (nstate == 1)
         {
            ode_solver_ws->GetState().Get(0,dxp);
            dx_u[0]->Distribute(dxp.GetBlock(0));
            dx_p[0]->Distribute(dxp.GetBlock(1));
            dx_phi[0]->Distribute(dxp.GetBlock(2));
         }

         // Actually write to file
         rdc.SetCycle(ri);
         rdc.SetTime(t);
         rdc.Save();
         ri++;

         // print meta file
         if (Mpi::Root())
         {
            std::ofstream step("restart/step.dat", std::ifstream::out);
            step<<t<<"\t"<<si<<"\t"<<ri<<"\t"<<vi<<endl;
            step<<dt<<endl;
            step.close();
         }
      }

      // Postprocess solution
      real_t cfl = form.GetCFL();
      real_t outflow = form.GetOutflow();
      Vector energy = form.GetEnergies(xp);
      DenseMatrix bdrForce = form.GetForce();

      // Write to file
      output.Print(si, t, dt, cfl, outflow, energy, bdrForce);

      // Write to screen
      if (Mpi::Root())
      {
         PrintForce2(pmesh.bdr_attributes,bdrForce);
         cout<<endl;
         line(80);
         cout<<" Kinetic Energy      = "<<energy[0]<<endl;
         cout<<" Potential Energy    = "<<energy[1]<<endl;
         cout<<" Viscous Dissipation = "<<energy[2]<<endl;
         cout<<" Mass Outflow Rate   = "<<outflow<<endl;
         cout<<" CFL-Number          = "<<cfl<<endl;
         if (dt_gain > 0) { cout<<" dt  = "<<dt0<<" --> "<<dt<<endl; }
         line(80);
         cout<<endl<<endl<<std::flush;
      }
   }


   // Free the used memory.
   for (int i = 0; i < fecs.Size(); ++i)
   {
      delete fecs[i];
   }
   for (int i = 0; i < spaces.Size(); ++i)
   {
      delete spaces[i];
   }

   return 0;
}
