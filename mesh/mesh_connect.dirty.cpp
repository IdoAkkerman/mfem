// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

// Implementation of data type mesh

#include "mesh_headers.hpp"
#include "../fem/fem.hpp"
#include "../general/sort_pairs.hpp"
#include "../general/text.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include <functional>

// Include the METIS header, if using version 5. If using METIS 4, the needed
// declarations are inlined below, i.e. no header is needed.
#if defined(MFEM_USE_METIS) && defined(MFEM_USE_METIS_5)
#include "metis.h"
#endif

#ifdef MFEM_USE_GECKO
#include "graph.h"
#endif

using namespace std;

namespace mfem
{


int relabel(Array<int> &map)
{
   // Relabel dofmap
   int size = map.Size();
   Array<int> tmp(size);
   tmp = 0;
   for (int j = 0; j < size; j++)
   {
      tmp[map[j]] = 1;
   }

   int cnt = 0;
   for (int j = 0; j < size; j++)
   {
      if (tmp[j] != 0)
      {
         tmp[j] = cnt++;
      }
   }

   for (int j = 0; j < size; j++)
   {
      map[j] = tmp[map[j]];
   }

   return cnt;
}


void RelabelElements(Array<int> &map, Array<Element *> &elements)
{
   for (int i = 0; i < elements.Size(); i++)
   {
      Element *el = elements[i];
      int *v = el->GetVertices();
      int nv = el->GetNVertices();

      for (int j = 0; j < nv; j++) { v[j] = map[v[j]]; }
   }
}

void RelabelTableRow(Array<int> &map, Table *table)
{
   for (int i = 0; i < table->Size(); i++)
   {
      int *row = table->GetRow(i);
      int size = table->RowSize(i);
      for (int j = 0; j < size; j++) { row[j] = map[row[j]]; }
   }
}

void RelabelTableCol(Array<int> &map, Table *oldt, Table *newt)
{
   for (int i = 0; i < oldt->Size(); i++)
   {
      int *oldr = oldt->GetRow(i);
      int *newr = newt->GetRow(map[i]);

      int size = oldt->RowSize(i);
      for (int j = 0; j < size; j++) { newr[j] = oldr[j]; }
   }
}


void Mesh::ConnectBoundaries(int bnd0, int bnd1, Array<int> &edge_to_knot)
{

   // Array<Element *> elements;  v2v @ #1
   // Array<Element *> boundary;  v2v @ #2
   // Array<Element *> faces;     v2v @ #3

   // Table *el_to_edge;          e2e @ #6
   // Table *el_to_face;          f2f @ #8

   // Array<int> be_to_edge; 2D   --> ignore routine for 3D at this stage
   // Table *bel_to_edge;    3D   e2e @ #7
   // Array<int> be_to_face;      f2f @ #9
   // mutable Table *face_edge;   e2e @ #5    f2f @ #10
   // mutable Table *edge_vertex; v2v @ #4    e2e @ #11

   if (NURBSext || ncmesh) { return; }
   cout<<"Connecting boundaries "<<bnd0<<" "<<bnd1<<endl;


   for (int i = 0; i < faces.Size(); i++)
   {
      Element *el = faces[i];
      int *v = el->GetVertices();
      int nv = el->GetNVertices();

      for (int j = 0; j < nv; j++) { cout <<" "<<v[j]; }
      cout<<endl;
   }

   cout<<"boundary stuff"<<endl;
   for (int i = 0; i < be_to_face.Size(); i++)
   {
      int f = be_to_face[i];
      int a = bdr_attributes[f];
      cout<<i<<" "<<f<<" "<<a<<endl;
   }

   // Find faces belonging to boundaries
   Array<int> bi0(GetNBE());
   Array<int> bi1(GetNBE());

   int bc0 = 0;
   int bc1 = 0;
   for (int i = 0; i < GetNBE(); i++)
   {
      if (GetBdrAttribute(i) == bnd0) { bi0[bc0++] = be_to_face[i]; }
      if (GetBdrAttribute(i) == bnd1) { bi1[bc1++] = be_to_face[i]; }
   }
   bi0.SetSize(bc0);
   bi1.SetSize(bc1);
   cout<<"bi0 "<<endl; bi0.Print();
   cout<<"bi1 "<<endl; bi1.Print();
   // Check found faces
   if (bc0 != bc1)
   {
      cout<<"Error when connecting boundaries :"<<bnd0<<" -> "<<bnd1<<endl;
      cout<<"  - Non matching face count:"<<bc0<<" -> "<<bc1<<endl;
      mfem_error();
   }

   // Generate Vertex to Vertex map
   //  - Assuming everything lines up nicely
   Array<int> v2v(GetNV());
   for (int i = 0; i < GetNV(); i++) { v2v[i] = i; }

   int fv2fv[]= {1,0,3,2};

   for (int i = 0; i < bc0; i++)
   {
      Element *el0 = GetBdrElement(bi0[i]);
      Element *el1 = GetBdrElement(bi1[i]);
      int *v0 = el0->GetVertices();
      int *v1 = el1->GetVertices();
      int nv0 = el0->GetNVertices();
      int nv1 = el1->GetNVertices();

      if (nv0 != nv1)
      {
         cout<<"Error when connecting boundaries   :"<<bnd0<<" -> "<<bnd1<<endl;
         cout<<"Error when connecting faces        :"<<bi0[i]<<" -> "<<bi1[i]<<endl;
         cout<<"  - Non matching face vertex count :"<<nv0<<" "<<nv1<<endl;
         mfem_error();
      }


      cout<<"v0[j] ";
      for (int j = 0; j < nv0; j++) { cout<<v0[fv2fv[j]]<<" "; }
      cout<<endl;

      cout<<"v1[j] ";
      for (int j = 0; j < nv0; j++) { cout<<v1[j]<<" "; }
      cout<<endl;


      for (int j = 0; j < nv0; j++) { v2v[v1[j]] = v0[fv2fv[j]]; }

   }
   cout<<"v2v -- dirty"<<endl;
   v2v.Print();
   NumOfVertices = relabel(v2v);
   cout<<"v2v -- clean"<<endl;
   v2v.Print();

   // Apply v2v map
   RelabelElements(v2v, elements);     // #1
   RelabelElements(v2v, boundary);     // #2
   RelabelElements(v2v, faces);        // #3
   RelabelTableRow(v2v, edge_vertex);  // #4

   /*cout<<"EDGES "<<endl;
      Array<int> vert(2);
      for (int i = 0; i < NumOfEdges; i++)
      {
         cout<<i<<" : ";
         edge_vertex->GetRow(i, vert);
         vert.Print(cout);
      }*/


   // Generate Edge to Edge map
   //  - Assuming everything lines up nicely

   cout<<"face_edge"<<endl;
   face_edge->Print();



   Array<int> e2e(NumOfEdges);
   for (int i = 0; i < NumOfEdges; i++) { e2e[i] = i; }
   int fe2fe[]= {0,3,2,1};
   for (int i = 0; i < bc0; i++)
   {
      int *row0 = face_edge->GetRow(bi0[i]);
      int size0 = face_edge->RowSize(bi0[i]);

      int *row1 = face_edge->GetRow(bi1[i]);
      int size1 = face_edge->RowSize(bi1[i]);

      if (size0 != size1)
      {
         cout<<"Error when connecting boundaries :"<<bnd0<<" -> "<<bnd1<<endl;
         cout<<"Error when connecting faces      :"<<bi0[i]<<" -> "<<bi1[i]<<endl;
         cout<<"  - Non matching face edge count :"<<size0<<" "<<size1<<endl;
         mfem_error();
      }

      cout<<"edge0[j] ";
      for (int j = 0; j < size0; j++) { cout<<row0[j]<<" "; }
      cout<<endl;

      //cout<<"edge0[j] ";
      //for (int j = 0; j < size0; j++) cout<<row0[fe2fe[j]]<<" ";
      //cout<<endl;


      cout<<"edge1[j] ";
      for (int j = 0; j < size0; j++) { cout<<row1[j]<<" "; }
      cout<<endl;

      //cout<<"edge1[j] ";
      //for (int j = 0; j < size0; j++) cout<<row1[fe2fe[j]]<<" ";
      //cout<<endl;
      // for (int j = 0; j < nv0; j++) v2v[v1[j]] = v0[fv2fv[j]];
      for (int j = 0; j < size0; j++) { e2e[row1[j]] = row0[fe2fe[j]]; }
   }
   cout<<"e2e -- dirty "<<endl;
   e2e.Print();

   NumOfEdges = relabel(e2e);
   cout<<"e2e -- clean "<<endl;
   e2e.Print();

   //cout<<237<<endl; e2e.Print();
   // Apply e2e map
   RelabelTableRow(e2e, face_edge);   // #5
   RelabelTableRow(e2e, el_to_edge);  // #6
   RelabelTableRow(e2e, bel_to_edge); // #7

   cout<<"face_edge #5"<<endl;
   face_edge->Print();

   // Apply e2e map -- #11
   Table *edge_vertex_old = edge_vertex;
   edge_vertex = new Table(NumOfEdges, 2);
   RelabelTableCol(e2e, edge_vertex_old, edge_vertex);
   delete edge_vertex_old;
   {
      Array<int> tmp;
      //edge_to_knot.Print();
      edge_to_knot.Copy(tmp);
      cout<<"edge_to_knot = "<<endl;
      tmp.Print();
      edge_to_knot.SetSize(NumOfEdges);
      cout<<"e2e = "<<endl;
      e2e.Print();
      cout<<"..... "<<endl;
      //edge_to_knot e2e
      for (int j = 0; j < tmp.Size(); j++)
      {
         cout<<j<<" "<<e2e[j]<<" "<<tmp[j]<<" "<<endl;
         edge_to_knot[e2e[j]] = tmp[j];
      }
      cout<<"edge_to_knot = "<<endl;
      edge_to_knot.Print();
   }


   //Array<int> vert;
   /*cout<<"EDGES "<<endl;
      for (int i = 0; i < NumOfEdges; i++)
      {
         cout<<i<<" : ";
         edge_vertex->GetRow(i, vert);
         vert.Print(cout);
      }*/

   // Generate Face to Face map
   //  - Assuming everything lines up nicely
   Array<int> f2f(NumOfFaces);
   for (int i = 0; i < NumOfFaces; i++) { f2f[i] = i; }

   for (int i = 0; i < bc0 ; i++)
   {
      f2f[bi1[i]] = f2f[bi0[i]];
   }
   NumOfFaces = relabel(f2f);

   // Apply e2e map -- #8
   RelabelTableRow(f2f, el_to_face);

   // Apply e2e map -- #9 #11
   cout<<" #9 #11"<<endl;
   for (int i = 0; i < be_to_face.Size(); i++)
   {
      int f = be_to_face[i];
      int a = bdr_attributes[f];
      cout<<i<<" "<<f<<" "<<a<<endl;
   }

   for (int i = 0; i < be_to_face.Size(); i++)
   {
      int f = be_to_face[i];
      int a = bdr_attributes[f];
      be_to_face[i] = f2f[f];
      bdr_attributes[f2f[f]] = bdr_attributes[f];
   }
   cout<<" after"<<endl;
   for (int i = 0; i < be_to_face.Size(); i++)
   {
      int f = be_to_face[i];
      int a = bdr_attributes[f];
      cout<<i<<" "<<f<<" "<<a<<endl;
   }

   // Apply e2e map -- #10
   Table *face_edge_old = face_edge;
   face_edge = new Table(NumOfFaces, 4);
   RelabelTableCol(f2f, face_edge_old,face_edge);
   delete face_edge_old;


   // Apply e2e map -- #11 bdr_attributes
   //cout<<" be_to_face #11"<<endl;
   //be_to_face.Print();
   //cout<<" bdr_attributes #11"<<endl;
   //   bdr_attributes.Print();
   // for (int i = 0; i < be_to_face.Size(); i++)
   // {
   //    bdr_attributes[i] = bdr_attributes[f2f[i]];
   // }
   // bdr_attributes.Print();

   //cout<<"face_edge #10"<<endl;
   //face_edge->Print();


   cout <<"NumOfEdges    " <<NumOfEdges<<endl;
   cout <<"NumOfFaces    " <<NumOfFaces<<endl;
   cout <<"NumOfVertices " <<NumOfVertices<<endl;


}

}
