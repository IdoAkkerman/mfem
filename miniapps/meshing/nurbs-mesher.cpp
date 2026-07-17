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
//    ---------------------------------------------------------------------
//    NURBS Mesher Miniapp: Construct a refined NURBS mesh from a coarse one
//    ---------------------------------------------------------------------
//
// This miniapp computes bounding boxes for each element in a given mesh, and
// also computes the bounds on the determinant of the Jacobian of the
// transformation for each element. The bounding approach is based on the
// method described in:
//
// (1) Section 3 of Mittal et al., "General Field Evaluation in High-Order
//     Meshes on GPUs"
// and
// (2) Dzanic et al., "A method for bounding high-order finite element
//     functions: Applications to mesh validity and bounds-preserving limiters".
//
//
// Compile with: make mesh-bounding-boxes
//
// Sample runs:
//  nurbs-mesher -m ../../data/square-nurbs.mesh -r 3 -o 4 -visit


#include "mfem.hpp"
#include <iostream>
#include <fstream>

using namespace mfem;
using namespace std;

int main (int argc, char *argv[])
{
   // Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   Hypre::Init();

   // Set the method's default parameters.
   const char *mesh_file = "../../data/square-nurbs.mesh";
   int mesh_poly_deg     = 1;
   const char *ref_file  = "";
   int ref_levels        = -1;
   const char *out_file  = "refined-nurbs.mesh";
   bool glvis            = true;
   int visport           = 19916;
   bool visit            = false;
   const char *visit_file  = "RefinedNURBS";

   // Parse command-line options.
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&mesh_poly_deg, "-o", "--order",
                  "Polynomial degree of mesh finite element space.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly, -1 for auto.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&out_file, "-of", "--out-file",
                  "Mesh output file");
   args.AddOption(&glvis, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visport, "-p", "--send-port", "Socket for GLVis.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit",
                  "--no-visit",
                  "Enable or disable VisIt output.");
   args.AddOption(&visit_file, "-vf", "--visit-file",
                  "Visit output file");

   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // Read the mesh from the given mesh file.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   MFEM_VERIFY(mesh->NURBSext, "Mehs must be a NURBS mesh.");

   // Increase polynomial order
   int mesh_poly_deg_old  = mesh->NURBSext->GetOrder();
   if (mesh_poly_deg > mesh_poly_deg_old)
   {
      mesh->DegreeElevate(mesh_poly_deg - mesh_poly_deg_old);
   }

   // Insert knots from file
   if (strlen(ref_file) != 0)
   {
      mesh->RefineNURBSFromFile(ref_file);
   }

   // Insert knots uniformly
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   // Report on final mesh statistics
   mesh->PrintInfo();

   // Save mesh in the mfem format
   ofstream mesh_ofs("refined.mesh");
   mesh->Print(mesh_ofs);
   mesh_ofs.close();

   // Visualize using glvis
   if (glvis)
   {
      char vishost[] = "localhost";
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << flush;
   }

   // Save mesh in the VisIt format
   if (visit)
   {
      VisItDataCollection visit_dc("NURBS-Mesh", mesh);
      visit_dc.Save();
   }
}
