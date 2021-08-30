// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_PMESH
#define MFEM_PMESH

#include "../config/config.hpp"

#ifdef MFEM_USE_MPI

#include "../general/communication.hpp"
#include "../general/globals.hpp"
#include "mesh.hpp"
#include "pncmesh.hpp"
#include <iostream>

namespace mfem
{
#ifdef MFEM_USE_PUMI
class ParPumiMesh;
#endif

/// Class for parallel meshes
class ParMesh : public Mesh
{
protected:
   MPI_Comm MyComm;
   int NRanks, MyRank;

   struct Vert3
   {
      int v[3];
      Vert3() = default;
      Vert3(int v0, int v1, int v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
      void Set(int v0, int v1, int v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
      void Set(const int *w) { v[0] = w[0]; v[1] = w[1]; v[2] = w[2]; }
   };

   struct Vert4
   {
      int v[4];
      Vert4() = default;
      Vert4(int v0, int v1, int v2, int v3)
      { v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3; }
      void Set(int v0, int v1, int v2, int v3)
      { v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3; }
      void Set(const int *w)
      { v[0] = w[0]; v[1] = w[1]; v[2] = w[2]; v[3] = w[3]; }
   };

   Array<Element *> shared_edges;
   // shared face id 'i' is:
   //   * triangle id 'i',                  if i < shared_trias.Size()
   //   * quad id 'i-shared_trias.Size()',  otherwise
   Array<Vert3> shared_trias;
   Array<Vert4> shared_quads;

   /// Shared objects in each group.
   Table group_svert;
   Table group_sedge;
   Table group_stria;  // contains shared triangle indices
   Table group_squad;  // contains shared quadrilateral indices

   /// Shared to local index mapping.
   Array<int> svert_lvert;
   Array<int> sedge_ledge;
   // sface ids: all triangles first, then all quads
   Array<int> sface_lface;

   Table *face_nbr_el_to_face;
   Table  face_nbr_el_ori; // orientations for each face (from nbr processor)

   IsoparametricTransformation FaceNbrTransformation;

   // glob_elem_offset + local element number defines a global element numbering
   mutable long glob_elem_offset, glob_offset_sequence;
   void ComputeGlobalElementOffset() const;

   /// Create from a nonconforming mesh.
   ParMesh(const ParNCMesh &pncmesh);

   // Convert the local 'meshgen' to a global one.
   void ReduceMeshGen();

   // Determine sedge_ledge and sface_lface.
   void FinalizeParTopo();

   // Mark all tets to ensure consistency across MPI tasks; also mark the
   // shared and boundary triangle faces using the consistently marked tets.
   virtual void MarkTetMeshForRefinement(DSTable &v_to_v);

   /// Return a number(0-1) identifying how the given edge has been split
   int GetEdgeSplittings(Element *edge, const DSTable &v_to_v, int *middle);
   /// Append codes identifying how the given face has been split to @a codes
   void GetFaceSplittings(const int *fv, const HashTable<Hashed2> &v_to_v,
                          Array<unsigned> &codes);

   bool DecodeFaceSplittings(HashTable<Hashed2> &v_to_v, const int *v,
                             const Array<unsigned> &codes, int &pos);

   STable3D *GetFaceNbrElementToFaceTable(int ret_ftbl = 0);

   void GetFaceNbrElementTransformation(
      int i, IsoparametricTransformation *ElTr);

   void GetGhostFaceTransformation(
      FaceElementTransformations* FETr, Element::Type face_type,
      Geometry::Type face_geom);

   /// Update the groups after triangle refinement
   void RefineGroups(const DSTable &v_to_v, int *middle);

   /// Update the groups after tetrahedron refinement
   void RefineGroups(int old_nv, const HashTable<Hashed2> &v_to_v);

   void UniformRefineGroups2D(int old_nv);

   // f2qf can be NULL if all faces are quads or there are no quad faces
   void UniformRefineGroups3D(int old_nv, int old_nedges,
                              const DSTable &old_v_to_v,
                              const STable3D &old_faces,
                              Array<int> *f2qf);

   void ExchangeFaceNbrData(Table *gr_sface, int *s2l_face);

   /// Refine a mixed 2D mesh uniformly.
   virtual void UniformRefinement2D();

   /// Refine a mixed 3D mesh uniformly.
   virtual void UniformRefinement3D();

   virtual void NURBSUniformRefinement();

   /// This function is not public anymore. Use GeneralRefinement instead.
   virtual void LocalRefinement(const Array<int> &marked_el, int type = 3);

   /// This function is not public anymore. Use GeneralRefinement instead.
   virtual void NonconformingRefinement(const Array<Refinement> &refinements,
                                        int nc_limit = 0);

   virtual bool NonconformingDerefinement(Array<double> &elem_error,
                                          double threshold, int nc_limit = 0,
                                          int op = 1);

   void RebalanceImpl(const Array<int> *partition);

   void DeleteFaceNbrData();

   bool WantSkipSharedMaster(const NCMesh::Master &master) const;

   /// Fills out partitioned Mesh::vertices
   int BuildLocalVertices(const Mesh& global_mesh, const int *partitioning,
                          Array<int> &vert_global_local);

   /// Fills out partitioned Mesh::elements
   int BuildLocalElements(const Mesh& global_mesh, const int *partitioning,
                          const Array<int> &vert_global_local);

   /// Fills out partitioned Mesh::boundary
   int BuildLocalBoundary(const Mesh& global_mesh, const int *partitioning,
                          const Array<int> &vert_global_local,
                          Array<bool>& activeBdrElem,
                          Table* &edge_element);

   void FindSharedFaces(const Mesh &mesh, const int* partition,
                        Array<int>& face_group,
                        ListOfIntegerSets& groups);

   int FindSharedEdges(const Mesh &mesh, const int* partition,
                       Table* &edge_element, ListOfIntegerSets& groups);

   int FindSharedVertices(const int *partition, Table* vertex_element,
                          ListOfIntegerSets& groups);

   void BuildFaceGroup(int ngroups, const Mesh &mesh,
                       const Array<int>& face_group,
                       int &nstria, int &nsquad);

   void BuildEdgeGroup(int ngroups, const Table& edge_element);

   void BuildVertexGroup(int ngroups, const Table& vert_element);

   void BuildSharedFaceElems(int ntri_faces, int nquad_faces,
                             const Mesh &mesh, int *partitioning,
                             const STable3D *faces_tbl,
                             const Array<int> &face_group,
                             const Array<int> &vert_global_local);

   void BuildSharedEdgeElems(int nedges, Mesh &mesh,
                             const Array<int> &vert_global_local,
                             const Table *edge_element);

   void BuildSharedVertMapping(int nvert, const Table* vert_element,
                               const Array<int> &vert_global_local);

   // Similar to Mesh::GetFacesTable()
   STable3D *GetSharedFacesTable();

   /// Ensure that bdr_attributes and attributes agree across processors
   void DistributeAttributes(Array<int> &attr);

   void LoadSharedEntities(std::istream &input);

   /// If the mesh is curved, make sure 'Nodes' is ParGridFunction.
   /** Note that this method is not related to the public 'Mesh::EnsureNodes`.*/
   void EnsureParNodes();

   /// Internal function used in ParMesh::MakeRefined (and related constructor)
   void MakeRefined_(ParMesh &orig_mesh, int ref_factor, int ref_type);

   // Mark Mesh::Swap as protected, should use ParMesh::Swap to swap @a ParMesh
   // objects.
   using Mesh::Swap;

   void Destroy();

public:
   /// Default constructor. Create an empty @a ParMesh.
   ParMesh() : MyComm(0), NRanks(0), MyRank(-1), face_nbr_el_to_face(NULL),
      glob_elem_offset(-1), glob_offset_sequence(-1),
      have_face_nbr_data(false), pncmesh(NULL) { }

   /// Create a parallel mesh by partitioning a serial Mesh.
   /** The mesh is partitioned automatically or using external partitioning
       data (the optional parameter 'partitioning_[i]' contains the desired MPI
       rank for element 'i'). Automatic partitioning uses METIS for conforming
       meshes and quick space-filling curve equipartitioning for nonconforming
       meshes (elements of nonconforming meshes should ideally be ordered as a
       sequence of face-neighbors). */
   ParMesh(MPI_Comm comm, Mesh &mesh, int *partitioning_ = NULL,
           int part_method = 1);

   /** Copy constructor. Performs a deep copy of (almost) all data, so that the
       source mesh can be modified (e.g. deleted, refined) without affecting the
       new mesh. If 'copy_nodes' is false, use a shallow (pointer) copy for the
       nodes, if present. */
   explicit ParMesh(const ParMesh &pmesh, bool copy_nodes = true);

   /// Read a parallel mesh, each MPI rank from its own file/stream.
   /** The @a refine parameter is passed to the method Mesh::Finalize(). */
   ParMesh(MPI_Comm comm, std::istream &input, bool refine = true);

   /// Deprecated: see @a ParMesh::MakeRefined
   MFEM_DEPRECATED
   ParMesh(ParMesh *orig_mesh, int ref_factor, int ref_type);

   /// Move constructor. Used for named constructors.
   ParMesh(ParMesh &&mesh);

   /// Move assignment operator.
   ParMesh& operator=(ParMesh &&mesh);

   /// Explicitly delete the copy assignment operator.
   ParMesh& operator=(ParMesh &mesh) = delete;

   /// Create a uniformly refined (by any factor) version of @a orig_mesh.
   /** @param[in] orig_mesh  The starting coarse mesh.
       @param[in] ref_factor The refinement factor, an integer > 1.
       @param[in] ref_type   Specify the positions of the new vertices. The
                             options are BasisType::ClosedUniform or
                             BasisType::GaussLobatto.

       The refinement data which can be accessed with GetRefinementTransforms()
       is set to reflect the performed refinements.

       @note The constructed ParMesh is linear, i.e. it does not have nodes. */
   static ParMesh MakeRefined(ParMesh &orig_mesh, int ref_factor, int ref_type);

   /** Create a mesh by splitting each element of @a orig_mesh into simplices.
       See @a Mesh::MakeSimplicial for more details. */
   static ParMesh MakeSimplicial(ParMesh &orig_mesh);

   virtual void Finalize(bool refine = false, bool fix_orientation = false);

   virtual void SetAttributes();

   MPI_Comm GetComm() const { return MyComm; }
   int GetNRanks() const { return NRanks; }
   int GetMyRank() const { return MyRank; }

   /** Map a global element number to a local element number. If the global
       element is not on this processor, return -1. */
   int GetLocalElementNum(long global_element_num) const;

   /// Map a local element number to a global element number.
   long GetGlobalElementNum(int local_element_num) const;

   /** The following functions define global indices for all local vertices,
       edges, faces, or elements. The global indices have no meaning or
       significance for ParMesh, but can be used for purposes beyond this class.
    */
   /// AMR meshes are not supported.
   void GetGlobalVertexIndices(Array<HYPRE_Int> &gi) const;
   /// AMR meshes are not supported.
   void GetGlobalEdgeIndices(Array<HYPRE_Int> &gi) const;
   /// AMR meshes are not supported.
   void GetGlobalFaceIndices(Array<HYPRE_Int> &gi) const;
   /// AMR meshes are supported.
   void GetGlobalElementIndices(Array<HYPRE_Int> &gi) const;

   GroupTopology gtopo;

   // Face-neighbor elements and vertices
   bool             have_face_nbr_data;
   Array<int>       face_nbr_group;
   Array<int>       face_nbr_elements_offset;
   Array<int>       face_nbr_vertices_offset;
   Array<Element *> face_nbr_elements;
   Array<Vertex>    face_nbr_vertices;
   // Local face-neighbor elements and vertices ordered by face-neighbor
   Table            send_face_nbr_elements;
   Table            send_face_nbr_vertices;

   ParNCMesh* pncmesh;

   int GetNGroups() const { return gtopo.NGroups(); }

   ///@{ @name These methods require group > 0
   int GroupNVertices(int group) { return group_svert.RowSize(group-1); }
   int GroupNEdges(int group)    { return group_sedge.RowSize(group-1); }
   int GroupNTriangles(int group) { return group_stria.RowSize(group-1); }
   int GroupNQuadrilaterals(int group) { return group_squad.RowSize(group-1); }

   int GroupVertex(int group, int i)
   { return svert_lvert[group_svert.GetRow(group-1)[i]]; }
   void GroupEdge(int group, int i, int &edge, int &o);
   void GroupTriangle(int group, int i, int &face, int &o);
   void GroupQuadrilateral(int group, int i, int &face, int &o);
   ///@}

   void GenerateOffsets(int N, HYPRE_BigInt loc_sizes[],
                        Array<HYPRE_BigInt> *offsets[]) const;

   void ExchangeFaceNbrData();
   void ExchangeFaceNbrNodes();

   virtual void SetCurvature(int order, bool discont = false, int space_dim = -1,
                             int ordering = 1);

   int GetNFaceNeighbors() const { return face_nbr_group.Size(); }
   int GetNFaceNeighborElements() const { return face_nbr_elements.Size(); }
   int GetFaceNbrGroup(int fn) const { return face_nbr_group[fn]; }
   int GetFaceNbrRank(int fn) const;

   /** Similar to Mesh::GetElementFaces */
   void GetFaceNbrElementFaces(int i, Array<int> &fcs, Array<int> &cor) const;

   /** Similar to Mesh::GetFaceToElementTable with added face-neighbor elements
       with indices offset by the local number of elements. */
   Table *GetFaceToAllElementTable() const;

   /** Get the FaceElementTransformations for the given shared face (edge 2D).
       In the returned object, 1 and 2 refer to the local and the neighbor
       elements, respectively. */
   FaceElementTransformations *
   GetSharedFaceTransformations(int sf, bool fill2 = true);

   ElementTransformation *
   GetFaceNbrElementTransformation(int i)
   {
      GetFaceNbrElementTransformation(i, &FaceNbrTransformation);

      return &FaceNbrTransformation;
   }

   /// Get the size of the i-th face neighbor element relative to the reference
   /// element.
   double GetFaceNbrElementSize(int i, int type=0);

   /// Return the number of shared faces (3D), edges (2D), vertices (1D)
   int GetNSharedFaces() const;

   /// Return the local face index for the given shared face.
   int GetSharedFace(int sface) const;

   /// See the remarks for the serial version in mesh.hpp
   MFEM_DEPRECATED virtual void ReorientTetMesh();

   /// Utility function: sum integers from all processors (Allreduce).
   virtual long ReduceInt(int value) const;

   /** Load balance the mesh by equipartitioning the global space-filling
       sequence of elements. Works for nonconforming meshes only. */
   void Rebalance();

   /** Load balance a nonconforming mesh using a user-defined partition.
       Each local element 'i' is migrated to processor rank 'partition[i]',
       for 0 <= i < GetNE(). */
   void Rebalance(const Array<int> &partition);

   /// Save the mesh in a parallel mesh format.
   void ParPrint(std::ostream &out) const;

   /** Print the part of the mesh in the calling processor adding the interface
       as boundary (for visualization purposes) using the mfem v1.0 format. */
   virtual void Print(std::ostream &out = mfem::out) const;

   /// Save the ParMesh to files (one for each MPI rank). The files will be
   /// given suffixes according to the MPI rank. The mesh will be written to the
   /// files using ParMesh::Print. The given @a precision will be used for ASCII
   /// output.
   virtual void Save(const char *fname, int precision=16) const;

#ifdef MFEM_USE_ADIOS2
   /** Print the part of the mesh in the calling processor using adios2 bp
       format. */
   virtual void Print(adios2stream &out) const;
#endif

   /** Print the part of the mesh in the calling processor adding the interface
       as boundary (for visualization purposes) using Netgen/Truegrid format .*/
   virtual void PrintXG(std::ostream &out = mfem::out) const;

   /** Write the mesh to the stream 'out' on Process 0 in a form suitable for
       visualization: the mesh is written as a disjoint mesh and the shared
       boundary is added to the actual boundary; both the element and boundary
       attributes are set to the processor number.  */
   void PrintAsOne(std::ostream &out = mfem::out) const;

   /// Save the mesh as a single file (using ParMesh::PrintAsOne). The given
   /// @a precision is used for ASCII output.
   void SaveAsOne(const char *fname, int precision=16) const;

   /// Old mesh format (Netgen/Truegrid) version of 'PrintAsOne'
   void PrintAsOneXG(std::ostream &out = mfem::out);

   /** Print the mesh in parallel PVTU format. The PVTU and VTU files will be
       stored in the directory specified by @a pathname. If the directory does
       not exist, it will be created. */
   virtual void PrintVTU(std::string pathname,
                         VTKFormat format=VTKFormat::ASCII,
                         bool high_order_output=false,
                         int compression_level=0,
                         bool bdr=false);

   /// Parallel version of Mesh::Load().
   virtual void Load(std::istream &input, int generate_edges = 0,
                     int refine = 1, bool fix_orientation = true);

   /// Returns the minimum and maximum corners of the mesh bounding box. For
   /// high-order meshes, the geometry is refined first "ref" times.
   void GetBoundingBox(Vector &p_min, Vector &p_max, int ref = 2);

   void GetCharacteristics(double &h_min, double &h_max,
                           double &kappa_min, double &kappa_max);

   /// Swaps internal data with another ParMesh, including non-geometry members.
   /// See @a Mesh::Swap
   void Swap(ParMesh &other);

   /// Print various parallel mesh stats
   virtual void PrintInfo(std::ostream &out = mfem::out);

   virtual int FindPoints(DenseMatrix& point_mat, Array<int>& elem_ids,
                          Array<IntegrationPoint>& ips, bool warn = true,
                          InverseElementTransformation *inv_trans = NULL);

   /// Debugging method
   void PrintSharedEntities(const char *fname_prefix) const;

   virtual ~ParMesh();

   friend class ParNCMesh;
#ifdef MFEM_USE_PUMI
   friend class ParPumiMesh;
#endif
#ifdef MFEM_USE_ADIOS2
   friend class adios2stream;
#endif
};

}

#endif // MFEM_USE_MPI

#endif
