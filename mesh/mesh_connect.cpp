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

Mesh::Mesh(Mesh *orig_mesh, Array<int>  &bnds0, Array<int>  &bnds1)
{
bdrOrient.SetSize(0);
bdrEdgeOrient.SetSize(0,0);

   if (bnds0.Size() != bnds1.Size()) { mfem_error("Mesh::Mesh() boundary lists not of equal size"); }

   if (!orig_mesh->NURBSext)
   {
      //cout<<"Standard mesh"<<endl;
      Copy(*orig_mesh, false);
      cout<<"Skip connection boundaries"<<endl;
      //Table *tmpTable = GetFaceEdgeTable();
      // GetBdrElementFaceOrientation();
      //cout<<"55"<<endl;
      //Array<int> edge_to_knot(0);
      //for (int i = 0; i < bnds0.Size(); i++)
      //{   cout<<i<<" 58"<<endl;
      //   ConnectBoundaries(bnds0[i], bnds1[i], edge_to_knot);
      //}
   }
   else
   {
      cout<<"knot_to_edge"<<endl;
      orig_mesh->NURBSext->edge_to_knot.Print();
      cout<<"--------------------------------------"<<endl;




      Copy(*orig_mesh, false);
      orig_mesh->GetBdrElementFaceOrientation();
      orig_mesh->GenerateBdrElementEdges();

      orig_mesh->bdrOrient.Copy(bdrOrient);
      orig_mesh->bdrEdgeOrient.Copy(bdrEdgeOrient);

 cout <<" chk 1 -->  "<<orig_mesh->bdrOrient.Size()<<"   "<<orig_mesh->NumOfBdrElements<<endl;
 cout <<" chk 1b -->  "<<bdrOrient.Size()<<"   "<<NumOfBdrElements<<endl;

      //   Table *tmpTable = GetFaceEdgeTable();
      //  GetBdrElementFaceOrientation();
      NURBSext = new NURBSExtension(orig_mesh->NURBSext, bnds0, bnds1);
 cout <<" chk 2 -->  "<<orig_mesh->bdrOrient.Size()<<"   "<<orig_mesh->NumOfBdrElements<<endl;
 cout <<" chk 2b -->  "<<bdrOrient.Size()<<"   "<<NumOfBdrElements<<endl;
 cout <<" chk 2c -->  "<<bdrOrient.Size()<<"   "<<NURBSext->GetNBE()<<endl;
      cout<<"knot_to_edge"<<endl;
      NURBSext->edge_to_knot.Print();
      cout<<"--------------------------------------"<<endl;

      Dim              = NURBSext->Dimension();
      spaceDim         = Dim;
      NumOfVertices    = NURBSext->GetNV();
      NumOfElements    = NURBSext->GetNE();
      NumOfBdrElements = NURBSext->GetNBE();

      bdrOrient.SetSize(NumOfBdrElements);
      bdrEdgeOrient.SetSize(NumOfBdrElements,4);

      Array<int> marker(orig_mesh->GetNBE());
      marker = 1;

      for (int i = 0; i < orig_mesh->GetNBE(); i++)
      {
         int attr = orig_mesh->GetBdrAttribute(i);
         for (int j = 0; j < bnds0.Size(); j++)
         {
            cout<<i<<" "<<j<<endl;
            if (attr == bnds0[j]) marker[i] = 0;
            if (attr == bnds1[j]) marker[i] = 0;
         }
      }
      cout<<"marker"<<endl;
      marker.Print();
      int cnt = 0;
      for (int i = 0; i < orig_mesh->GetNBE(); i++)
      {
         if (marker[i] == 1)
         {
            bdrOrient[cnt] = orig_mesh->bdrOrient[i];
            for (int j = 0; j < 4; j++) bdrEdgeOrient(cnt,j) = orig_mesh->bdrEdgeOrient(i,j);
            cnt++;
         }
      }
      cout<<"bdrOrient"<<endl;
      bdrOrient.Print();




      //orig_mesh->bdrEdgeOrient.Copy(bdrEdgeOrient);

      cout<<"--------------------------------------"<<endl;
      cout<<"Dim              = "<<NURBSext->Dimension()<<endl;
      cout<<"NumOfVertices    = "<<NURBSext->GetNV()<<endl;
      cout<<"NumOfElements    = "<<NURBSext->GetNE()<<endl;
      cout<<"NumOfBdrElements = "<<NURBSext->GetNBE()<<endl;
      cout<<"--------------------------------------"<<endl;
      cout<<88<<endl;
      NURBSext->GetElementTopo(elements);
      cout<<90<<endl;
      NURBSext->GetBdrElementTopo(boundary);


      cout<<92<<endl;
            GenerateFaces();
      cout<<94<<endl;
            Table *tmpTable = GetFaceEdgeTable();
      cout<<96<<endl;
            GetBdrElementFaceOrientation();
      cout<<98<<endl;
           // cout<<"elements"<<endl;
            for (int i = 0; i < elements.Size(); i++)
            {
               Element *el = elements[i];
               int *v = el->GetVertices();
               int nv = el->GetNVertices();

             //  for (int j = 0; j < nv; j++) { cout<<v[j]<<" "; }
             //  cout<<endl;
            }
            cout<<"boundary"<<endl;
            for (int i = 0; i < boundary.Size(); i++)
            {
               Element *el = boundary[i];
               int *v = el->GetVertices();
               int nv = el->GetNVertices();
               cout<<GetBdrAttribute(i)<<" ";
               for (int j = 0; j < nv; j++) { cout<<v[j]<<" "; }
               cout<<endl;
            }
             cout<<"hatse1234  begin"<<endl;
            {

               int i, *v;
               STable3D *faces_tbl;

               if (el_to_face != NULL)
               {
                  delete el_to_face;
               }
               el_to_face = new Table(NumOfElements, 6);  // must be 6 for hexahedra
               faces_tbl = new STable3D(NumOfVertices);
               int cnt= 0;
               for (i = 0; i < NumOfElements; i++)
               {
                  v = elements[i]->GetVertices();


      if (i < 993)
      {
         for (int j = 0; j < 8; j++) cout<<v[j]<<" ";
         cout<<endl;
      }
      //cout<<GetElementType(i)<<endl;
      //cout<<Element::TETRAHEDRON<<endl;
      //cout<< Element::HEXAHEDRON<<endl;

                  switch (GetElementType(i))
                  {

                     case Element::TETRAHEDRON:
                     {cout<<" TET??"<<endl;
                        for (int j = 0; j < 4; j++)
                        {
                           const int *fv = tet_t::FaceVert[j];
                           el_to_face->Push(
                              i, faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]));
                        }
                        break;
                     }
                     case Element::HEXAHEDRON:
                     {
      //cout<<" QUAD"<<endl;
                        // find the face by the vertices with the smallest 3 numbers
                        // z = 0, y = 0, x = 1, y = 1, x = 0, z = 1
                        //cout<<"Element::HEXAHEDRON"<<endl;
                        for (int j = 0; j < 6; j++)
                        {
                           const int *fv = hex_t::FaceVert[j];

                           int tmp = faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
      if (i < -3){
                           cout<<j<<" : "<<v[fv[0]]<<" "<<v[fv[1]]<<" "<< v[fv[2]]<<" "<< v[fv[3]];
                           cout<<" -> "<<tmp<<endl;
      }

                           el_to_face->Push(
                              i, tmp);
                        }
                        break;
                     }
                     default:
                        MFEM_ABORT("Unexpected type of Element.");
                  }
               }

               el_to_face->Finalize();
               //cout<<"el_to_face->Finalize();"<<endl;
               //el_to_face->Print();
               NumOfFaces = faces_tbl->NumberOfElements();
               be_to_face.SetSize(NumOfBdrElements);

               //cout<<"bdr_attributes"<<endl;
               //bdr_attributes.Print();
          cout<<248<<endl;
               cnt= 0;
               for (i = 0; i < NumOfBdrElements; i++)
               {
                  cout<<"bnd = "<<i<<endl;
                  v = boundary[i]->GetVertices();
                  switch (GetBdrElementType(i))
                  {
                     case Element::TRIANGLE:
                     {
                        //            cout<<"Element::TRIANGLE"<<endl;
                        be_to_face[i] = (*faces_tbl)(v[0], v[1], v[2]);
                        break;
                     }
                     case Element::QUADRILATERAL:
                     {
                        cout<<"Element::QUADRILATERAL"<<endl;

                        cout<<i<<" "<<boundary[i]->GetAttribute()<<endl;
                        cout<<cnt++<<" : "<<v[0]<<" "<<v[1]<<" "<< v[2]<<" "<< v[3]<<endl;
                        be_to_face[i] = (*faces_tbl)(v[0], v[1], v[2], v[3]);
                        //cout<<be_to_face[i]<<endl;
                        //(*faces_tbl)
                        break;
                     }
                     default:
                        MFEM_ABORT("Unexpected type of boundary Element.");
                  }
               }
               cout<<188<<endl;

               delete faces_tbl;

            }
            cout<<"hatse1234 "<<endl;
/*


            bdrOrient.SetSize(0);
            cout<<"be_to_face"<<endl;
      be_to_face.Print();

            GetBdrElementFaceOrientation();
            cout<<"GetElementToFaceTable"<<endl;
       // STable3D *faces_tbl = GetElementToFaceTable(0);


            vertices.SetSize(NumOfVertices);

            Nodes = orig_mesh->Nodes;
            own_nodes = 0;
            // NURBSext->SetCoordsFromPatches(*Nodes);

            int vd = Nodes->VectorDim();
            for (int i = 0; i < vd; i++)
            {
               Vector vert_val;
               Nodes->GetNodalValues(vert_val, i+1);

               for (int j = 0; j < NumOfVertices; j++)
               {
                  vertices[j](i) = vert_val(j);
               }
            }*/
   }

   cout<<"Mesh::Mesh() -->  end"<<endl;
   // exit(0);

}

void Mesh::ConnectBnd2( Array<int>  &bnds0, Array<int>  &bnds1)
{
   if (bnds0.Size() != bnds1.Size()) { mfem_error("Mesh::Mesh() boundary lists not of equal size"); }

   if (!NURBSext)
   {
      cout<<"Skip connection boundaries"<<endl;
      return;
   }

   Table *tmpTable = GetFaceEdgeTable();
   GetBdrElementFaceOrientation();
   NURBSext = new NURBSExtension(NURBSext, bnds0, bnds1);

   Dim              = NURBSext->Dimension();
   spaceDim         = Dim;
   NumOfVertices    = NURBSext->GetNV();
   NumOfElements    = NURBSext->GetNE();
   NumOfBdrElements = NURBSext->GetNBE();

   cout<<"Dim              = "<<NURBSext->Dimension()<<endl;
   cout<<"NumOfVertices    = "<<NURBSext->GetNV()<<endl;
   cout<<"NumOfElements    = "<<NURBSext->GetNE()<<endl;
   cout<<"NumOfBdrElements = "<<NURBSext->GetNBE()<<endl;

   NURBSext->GetElementTopo(elements);
   NURBSext->GetBdrElementTopo(boundary);


}

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

void Mesh::ConnectBoundaries(int bnd0, int bnd1, Array<int> &edge_to_knot)
{

   if (edge_to_knot.Size() == 0)
   {
      edge_to_knot.SetSize(NumOfEdges);
      edge_to_knot = 0;
   }

   if (NURBSext || ncmesh) { return; }
   cout<<"Connecting boundaries "<<bnd0<<" "<<bnd1<<endl;

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
      //cout<<nv1<<" "<<nv0<<endl;
      if (nv0 != nv1)
      {
         cout<<"Error when connecting boundaries   :"<<bnd0<<" -> "<<bnd1<<endl;
         cout<<"Error when connecting faces        :"<<bi0[i]<<" -> "<<bi1[i]<<endl;
         cout<<"  - Non matching face vertex count :"<<nv0<<" "<<nv1<<endl;
         mfem_error();
      }
      for (int j = 0; j < nv0; j++) { v2v[v1[j]] = v0[fv2fv[j]]; }

   }
   NumOfVertices = relabel(v2v);

   // Apply v2v map
   RelabelElements(v2v, elements);     // #1
   RelabelElements(v2v, boundary);     // #2
   RelabelElements(v2v, faces);        // #3
   RelabelTableRow(v2v, edge_vertex);  // #4

   // Generate Edge to Edge map
   //  - Assuming everything lines up nicely
   Array<int> e2e(NumOfEdges);
   Array<int> e2k(NumOfEdges);
   edge_to_knot.Copy(e2k);

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
      for (int j = 0; j < size0; j++)
      {
         e2e[row1[j]] = row0[fe2fe[j]];
         e2k[row1[j]] = edge_to_knot[row0[fe2fe[j]]];
      }
   }
   cout<<"-------------------"<<endl;
   edge_to_knot.Print();

   cout<<"-------------------"<<endl;
   for (int i = 0; i < e2e.Size(); i++) { cout<<i<<" "<<e2e[i]<<endl; }


   NumOfEdges = relabel(e2e);
   cout<<"-------------------"<<endl;
   //e2e.Print();
   for (int i = 0; i < e2e.Size(); i++) { cout<<i<<" "<<e2e[i]<<endl; }

   // Apply e2e map
   RelabelTableRow(e2e, face_edge);   // #5
   RelabelTableRow(e2e, el_to_edge);  // #6
   RelabelTableRow(e2e, bel_to_edge); // #7

   // Apply e2e map -- #11
   Table *edge_vertex_old = edge_vertex;
   edge_vertex = new Table(NumOfEdges, 2);
   RelabelTableCol(e2e, edge_vertex_old, edge_vertex);
   delete edge_vertex_old;

   // Apply e2e map -- #edge_to_knot
   edge_to_knot.SetSize(NumOfEdges);
   for (int j = 0; j < e2k.Size(); j++)
   {
      edge_to_knot[e2e[j]] = e2k[j];
   }
   cout<<"-------------------"<<endl;
   edge_to_knot.Print();
   //exit(0);

   // Generate Face to Face map
   //  - Assuming everything lines up nicely
   Array<int> f2f(NumOfFaces);
   for (int i = 0; i < NumOfFaces; i++) { f2f[i] = i; }

   for (int i = 0; i < bc0 ; i++)
   {
      f2f[bi1[i]] = f2f[bi0[i]];
   }
   NumOfFaces = relabel(f2f);
   //cout<<"f2f = ";f2f.Print();
   // Apply e2e map -- #8
   RelabelTableRow(f2f, el_to_face);

   // Apply e2e map -- #9 #11
   for (int i = 0; i < be_to_face.Size(); i++)
   {
      int f = be_to_face[i];
      int a = bdr_attributes[f];
      be_to_face[i] = f2f[f];
      bdr_attributes[f2f[f]] = bdr_attributes[f];
   }

   // Apply e2e map -- #10
   Table *face_edge_old = face_edge;
   face_edge = new Table(NumOfFaces, 4);
   RelabelTableCol(f2f, face_edge_old,face_edge);
   delete face_edge_old;


   // Remove used boundaries from the list
   Table *bel_to_edge_old = bel_to_edge;
   bel_to_edge = new Table(NumOfBdrElements - 2*bc0, 4);
   int cnt = 0;
   for (int i = 0; i < NumOfBdrElements; i++)
   {
      if ((GetBdrAttribute(i) != bnd0) &&
          (GetBdrAttribute(i) != bnd1) )
      {
         boundary[cnt] = boundary[i];
         be_to_face[cnt] = be_to_face[i];
         bdrOrient[cnt] = bdrOrient[i];

         int *oldr = bel_to_edge_old->GetRow(i);
         int *newr = bel_to_edge    ->GetRow(cnt);

         int size = bel_to_edge_old->RowSize(i);
         for (int j = 0; j < size; j++) { newr[j] = oldr[j]; }
         cnt++;
      }
      else
      {
         delete boundary[i];
      }
   }
   NumOfBdrElements = cnt;
   boundary.SetSize(NumOfBdrElements);
   be_to_face.SetSize(NumOfBdrElements);
   bdrOrient.SetSize(NumOfBdrElements);


   /*
   cout <<" NumOfBdrElements "<<NumOfBdrElements << endl;
   cout<<"bel_to_edge"<<endl;
   bel_to_edge->Print();
   cout<<"-----------------------"<<endl;
   be_to_face.Print();
   cout<<"-----------------------"<<endl;

      for (int i = 0; i < NumOfBdrElements; i++)
      {
         cout<<GetBdrAttribute(i)<<endl;
      }*/


   //exit(0);

   cout<<"--------------------------------------"<<endl;
   cout<<"Connected boundaries "<<bnd0<<" "<<bnd1<<endl;
   cout<<"--------------------------------------"<<endl;
   cout <<"NumOfEdges    " <<NumOfEdges<<endl;
   cout <<"NumOfFaces    " <<NumOfFaces<<endl;
   cout <<"NumOfVertices " <<NumOfVertices<<endl;
   cout<<"--------------------------------------"<<endl;
}

}
