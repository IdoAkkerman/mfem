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

#include "fem/univariate_assembly.hpp"
#include "fem/weighted_quadrature.hpp"
#include "linalg/densemat_gen_eigensystem.hpp"
#include "linalg/fast_diagonalization.hpp"
#include "linalg/kronecker_preconditioner.hpp"
#include "mfem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#if defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

using namespace mfem;

namespace {

using Clock = std::chrono::steady_clock;

double ElapsedSec(Clock::time_point start, Clock::time_point end) {
  return std::chrono::duration<double>(end - start).count();
}

long GetPeakRssKb() {
#if defined(__linux__)
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    return ru.ru_maxrss;
  }
  return 0;
#elif defined(__APPLE__)
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    return ru.ru_maxrss / 1024;
  }
  return 0;
#else
  return 0;
#endif
}

std::string DeriveMeshName(const std::string &path) {
  const std::size_t slash = path.find_last_of("/\\");
  const std::string fname =
      (slash == std::string::npos) ? path : path.substr(slash + 1);
  const std::size_t dot = fname.find_last_of('.');
  if (dot == std::string::npos || dot == 0) {
    return fname;
  }
  return fname.substr(0, dot);
}

std::string JsonString(const std::string &s) {
  std::ostringstream oss;
  oss << '"';
  for (unsigned char c : s) {
    switch (c) {
    case '"':
      oss << "\\\"";
      break;
    case '\\':
      oss << "\\\\";
      break;
    case '\b':
      oss << "\\b";
      break;
    case '\f':
      oss << "\\f";
      break;
    case '\n':
      oss << "\\n";
      break;
    case '\r':
      oss << "\\r";
      break;
    case '\t':
      oss << "\\t";
      break;
    default:
      if (c < 0x20) {
        oss << "\\u00" << std::hex << std::setw(2) << std::setfill('0')
            << static_cast<int>(c) << std::dec;
      } else {
        oss << static_cast<char>(c);
      }
      break;
    }
  }
  oss << '"';
  return oss.str();
}

std::string JsonBool(bool val) { return val ? "true" : "false"; }

std::string JsonNumber(double val) {
  if (!std::isfinite(val)) {
    return "null";
  }
  if (val == 0.0) {
    val = 0.0;
  }
  std::ostringstream oss;
  oss << std::setprecision(17) << val;
  return oss.str();
}

std::string JsonNullableNumber(bool has_val, double val) {
  if (!has_val) {
    return "null";
  }
  return JsonNumber(val);
}

std::string JsonIntArray(const std::vector<int> &vals) {
  std::ostringstream oss;
  oss << '[';
  for (std::size_t i = 0; i < vals.size(); i++) {
    if (i > 0) {
      oss << ", ";
    }
    oss << vals[i];
  }
  oss << ']';
  return oss.str();
}

template <typename T> std::string JsonNumberArray(const std::vector<T> &vals) {
  std::ostringstream oss;
  oss << '[';
  for (std::size_t i = 0; i < vals.size(); i++) {
    if (i > 0) {
      oss << ", ";
    }
    oss << JsonNumber(static_cast<double>(vals[i]));
  }
  oss << ']';
  return oss.str();
}

class ResidualRecorderMonitor : public IterativeSolverMonitor {
public:
  std::vector<real_t> residuals;

  void Reset() override { residuals.clear(); }

  void MonitorResidual(int it, real_t norm, const Vector &r,
                       bool final) override {
    MFEM_CONTRACT_VAR(it);
    MFEM_CONTRACT_VAR(r);
    MFEM_CONTRACT_VAR(final);
    residuals.push_back(norm);
  }
};

} // namespace

int main(int argc, char *argv[]) {
  const char *mesh_file = "meshes/square-nurbs-distorted.mesh";
  const char *mesh_name_arg = "";
  const char *problem_arg = "const-1";
  const char *ref_strategy_arg = "k";
  int requested_order = -1;
  int ref_levels = 0;
  const char *assembly_arg = "legacy";
  bool wq_symmetrize = true;
  int wq_interior_pts = 2;
  const char *solver_arg = "cg";
  const char *prec_arg = "kronecker-fd";
  real_t rtol = 1e-10;
  real_t atol = 0.0;
  int max_iter = 2000;
  int gmres_restart = 50;
  int repeats = 1;
  bool verify = false;
  bool record_residuals = false;
  const char *out_file = "";
  const char *run_id = "";
  bool quiet = false;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "NURBS mesh file to use.");
  args.AddOption(&mesh_name_arg, "-mn", "--mesh-name",
                 "Logical mesh name (defaults to filename stem).");
  args.AddOption(&problem_arg, "-prob", "--problem", "PDE problem: 'const-1'.");
  args.AddOption(&ref_strategy_arg, "-rs", "--ref-strategy",
                 "Refinement strategy: 'k', 'p', or 'h'.");
  args.AddOption(&requested_order, "-o", "--order",
                 "Target NURBS order (-1 keeps mesh base order).");
  args.AddOption(&ref_levels, "-r", "--ref-levels",
                 "Number of uniform h-refinement levels.");
  args.AddOption(&assembly_arg, "-asm", "--assembly",
                 "Operator assembly method: 'legacy' or 'wq'.");
  args.AddOption(&wq_symmetrize, "-wq-sym", "--wq-symmetrize", "-no-wq-sym",
                 "--no-wq-symmetrize",
                 "Symmetrize Weighted Quadrature matrix (default: true).");
  args.AddOption(&wq_interior_pts, "-wq-pts", "--wq-interior-pts",
                 "Interior quadrature points per knot span for WQ (>= 2).");
  args.AddOption(&solver_arg, "-sol", "--solver",
                 "Krylov solver: 'cg', 'gmres', or 'minres'.");
  args.AddOption(&prec_arg, "-prec", "--prec",
                 "Preconditioner: 'none', 'jacobi', 'gs', 'kronecker-fd'.");
  args.AddOption(&rtol, "-rtol", "--rtol", "Relative solver tolerance.");
  args.AddOption(&atol, "-atol", "--atol", "Absolute solver tolerance.");
  args.AddOption(&max_iter, "-mi", "--max-iter", "Maximum Krylov iterations.");
  args.AddOption(&gmres_restart, "-kdim", "--gmres-restart",
                 "GMRES restart dimension.");
  args.AddOption(&repeats, "-rep", "--repeats",
                 "Number of timed solve repeats (>= 1).");
  args.AddOption(&verify, "-ver", "--verify", "-no-ver", "--no-verify",
                 "Run preconditioner stale-y and symmetry checks.");
  args.AddOption(&record_residuals, "-rec", "--record-residuals", "-no-rec",
                 "--no-record-residuals",
                 "Record per-iteration residuals on the final repeat.");
  args.AddOption(&out_file, "-out", "--out-file", "Output JSON file path.");
  args.AddOption(&run_id, "-id", "--run-id",
                 "External benchmark run identifier.");
  args.AddOption(&quiet, "-q", "--quiet", "-no-q", "--no-quiet",
                 "Suppress verbose solver output.");

  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem::out);
    return 1;
  }
  if (!quiet) {
    args.PrintOptions(mfem::out);
  }

  const std::string mesh_path(mesh_file);
  const std::string mesh_name = (mesh_name_arg[0] != '\0')
                                    ? std::string(mesh_name_arg)
                                    : DeriveMeshName(mesh_path);
  const std::string problem(problem_arg);
  const std::string ref_strategy(ref_strategy_arg);
  const std::string assembly(assembly_arg);
  const std::string solver_name(solver_arg);
  const std::string prec_name(prec_arg);
  const std::string out_path(out_file);
  const std::string run_id_str(run_id);

  // Validate CLI arguments up front.
  MFEM_VERIFY(problem == "const-1", "Unsupported --problem '"
                                        << problem
                                        << "'. Only 'const-1' is supported.");
  MFEM_VERIFY(ref_strategy == "k" || ref_strategy == "p" || ref_strategy == "h",
              "Unsupported --ref-strategy '"
                  << ref_strategy << "'. Expected 'k', 'p', or 'h'.");
  MFEM_VERIFY(ref_levels >= 0,
              "Invalid --ref-levels " << ref_levels << ". Must be >= 0.");
  MFEM_VERIFY(assembly == "legacy" || assembly == "wq",
              "Unsupported --assembly '" << assembly
                                         << "'. Expected 'legacy' or 'wq'.");
  MFEM_VERIFY(wq_interior_pts >= 2, "Invalid --wq-interior-pts "
                                        << wq_interior_pts
                                        << ". Must be >= 2.");
  MFEM_VERIFY(solver_name == "cg" || solver_name == "gmres" ||
                  solver_name == "minres",
              "Unsupported --solver '"
                  << solver_name << "'. Expected 'cg', 'gmres', or 'minres'.");
  MFEM_VERIFY(prec_name == "none" || prec_name == "jacobi" ||
                  prec_name == "gs" || prec_name == "kronecker-fd",
              "Unsupported --prec '"
                  << prec_name
                  << "'. Expected 'none', 'jacobi', 'gs', or 'kronecker-fd'.");
  MFEM_VERIFY(max_iter > 0,
              "Invalid --max-iter " << max_iter << ". Must be > 0.");
  MFEM_VERIFY(gmres_restart > 0,
              "Invalid --gmres-restart " << gmres_restart << ". Must be > 0.");
  MFEM_VERIFY(repeats >= 1,
              "Invalid --repeats " << repeats << ". Must be >= 1.");

  if (assembly == "wq" && !wq_symmetrize &&
      (solver_name == "cg" || solver_name == "minres")) {
    MFEM_WARNING("Unsymmetrized Weighted Quadrature (--no-wq-symmetrize) "
                 "produces a slightly nonsymmetric matrix; solver '"
                 << solver_name
                 << "' may lose orthogonality or fail to converge.");
  }

  // 1. Load NURBS mesh.
  const Clock::time_point t_mesh_start = Clock::now();
  Mesh mesh(mesh_file, 1, 1);
  const Clock::time_point t_mesh_end = Clock::now();
  const double t_mesh_load = ElapsedSec(t_mesh_start, t_mesh_end);

  const FiniteElementSpace *nodal_fes = mesh.GetNodalFESpace();
  MFEM_VERIFY(nodal_fes != nullptr && dynamic_cast<const NURBSFECollection *>(
                                          nodal_fes->FEColl()) != nullptr,
              "Input mesh '" << mesh_path << "' is not a NURBS mesh.");
  MFEM_VERIFY(mesh.GetNE() > 0, "Loaded mesh contains no elements.");
  MFEM_VERIFY(mesh.bdr_attributes.Size() > 0,
              "Loaded mesh has no boundary attributes.");

  const int dim = mesh.Dimension();
  const int base_order = nodal_fes->GetFE(0)->GetOrder();

  MFEM_VERIFY(requested_order == -1 || requested_order >= base_order,
              "Requested --order (" << requested_order
                                    << ") must be -1 or >= mesh base order ("
                                    << base_order << ").");
  if (ref_strategy == "h") {
    MFEM_VERIFY(requested_order == -1 || requested_order == base_order,
                "Refinement strategy 'h' requires --order to be -1 or equal "
                "to mesh base order ("
                    << base_order << "), got " << requested_order << ".");
  }

  const int target_order =
      (requested_order == -1) ? base_order : requested_order;
  const int rel_degree = target_order - base_order;

  // 2. Mesh refinement (k-, p-, or h-refinement).
  const Clock::time_point t_ref_start = Clock::now();
  if (ref_strategy == "k") {
    if (rel_degree > 0) {
      mesh.DegreeElevate(rel_degree, target_order);
    }
    for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
    }
  } else if (ref_strategy == "p") {
    for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
    }
    if (rel_degree > 0) {
      mesh.DegreeElevate(rel_degree, target_order);
    }
  } else // ref_strategy == "h"
  {
    for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
    }
  }
  const Clock::time_point t_ref_end = Clock::now();
  const double t_refinement = ElapsedSec(t_ref_start, t_ref_end);

  // 3. Finite element space and homogeneous Dirichlet BCs on all boundaries.
  const Clock::time_point t_fes_start = Clock::now();
  NURBSFECollection fec;
  FiniteElementSpace fespace(&mesh, &fec);

  Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  const Clock::time_point t_fes_end = Clock::now();
  const double t_fespace_setup = ElapsedSec(t_fes_start, t_fes_end);

  const int final_order = fespace.GetFE(0)->GetOrder();
  const int num_elements = mesh.GetNE();
  const int num_tdofs = fespace.GetTrueVSize();
  const int num_ess_dofs = ess_tdof_list.Size();
  const int num_free_dofs = num_tdofs - num_ess_dofs;

  // 4. RHS assembly: -Delta u = 1.
  ConstantCoefficient one(1.0);
  const Clock::time_point t_rhs_start = Clock::now();
  LinearForm b(&fespace);
  b.AddDomainIntegrator(new DomainLFIntegrator(one));
  b.Assemble();
  const Clock::time_point t_rhs_end = Clock::now();
  const double t_rhs_assembly = ElapsedSec(t_rhs_start, t_rhs_end);

  // 5. Operator assembly + linear system formation (legacy or WQ).
  std::unique_ptr<BilinearForm> a;
  GridFunction x(&fespace);
  OperatorPtr A;
  Vector B, X;

  double t_op_wq_rules_1d = 0.0;
  double t_op_wq_weights_1d = 0.0;
  double t_op_wq_metric_eval = 0.0;
  double t_op_wq_patch_assemble = 0.0;
  double t_op_wq_symmetrize = 0.0;
  double t_op_wq_bc_elim = 0.0;

  const Clock::time_point t_op_start = Clock::now();
  x = 0.0;
  if (assembly == "legacy") {
    a = std::make_unique<BilinearForm>(&fespace);
    a->SetAssemblyLevel(AssemblyLevel::LEGACY);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    a->Assemble();
    a->FormLinearSystem(ess_tdof_list, x, b, A, X, B);
  } else // assembly == "wq"
  {
    WQOptions wq_opts;
    wq_opts.symmetrize = wq_symmetrize;
    wq_opts.interior_pts_per_span = wq_interior_pts;
    WQTimings wq_timings{};

    SparseMatrix *A_wq = new SparseMatrix();
    AssembleDiffusionWQ(fespace, one, ess_tdof_list, x, b, *A_wq, X, B, wq_opts,
                        &wq_timings);
    A.Reset(A_wq, true);

    t_op_wq_rules_1d = wq_timings.rules_1d;
    t_op_wq_weights_1d = wq_timings.weights_1d;
    t_op_wq_metric_eval = wq_timings.metric_eval;
    t_op_wq_patch_assemble = wq_timings.patch_assemble;
    t_op_wq_symmetrize = wq_timings.symmetrize;
    t_op_wq_bc_elim = wq_timings.bc_elim;
  }
  const Clock::time_point t_op_end = Clock::now();
  const double t_operator_assembly = ElapsedSec(t_op_start, t_op_end);

  const SparseMatrix *A_sp = A.As<SparseMatrix>();
  const int matrix_nnz = (A_sp != nullptr) ? A_sp->NumNonZeroElems() : 0;

  // 6. Preconditioner setup.
  double t_prec_univariate = 0.0;
  double t_prec_eigendecomp = 0.0;
  double t_prec_fd_construct = 0.0;
  double t_prec_setup_total = 0.0;

  std::vector<int> full_sizes_1d;
  std::vector<int> free_sizes_1d;

  UnivariateFactors uni;
  std::vector<DenseMatrix> U_1d;
  std::vector<Vector> lambda_1d;
  std::unique_ptr<FastDiagonalizationSolver> fd_ptr;
  std::unique_ptr<Solver> prec_ptr;

  if (prec_name == "kronecker-fd") {
    const Clock::time_point t0 = Clock::now();
    AssembleUnivariateLaplace(fespace, ess_bdr, uni);
    const Clock::time_point t1 = Clock::now();
    t_prec_univariate = ElapsedSec(t0, t1);

    U_1d.resize(dim);
    lambda_1d.resize(dim);
    for (int d = 0; d < dim; d++) {
      DenseMatrixGenEigensystem ge(uni.K[d], uni.M[d]);
      ge.Eval();
      U_1d[d] = ge.Eigenvectors();
      lambda_1d[d] = ge.Eigenvalues();
    }
    const Clock::time_point t2 = Clock::now();
    t_prec_eigendecomp = ElapsedSec(t1, t2);

    fd_ptr = std::make_unique<FastDiagonalizationSolver>(U_1d, lambda_1d);
    prec_ptr = std::make_unique<KroneckerPreconditioner>(
        fespace.GetTrueVSize(), ess_tdof_list, uni.full_sizes_1d,
        uni.free_dofs_1d, uni.tensor_to_tdof, *fd_ptr);
    const Clock::time_point t3 = Clock::now();
    t_prec_fd_construct = ElapsedSec(t2, t3);

    t_prec_setup_total =
        t_prec_univariate + t_prec_eigendecomp + t_prec_fd_construct;

    full_sizes_1d.reserve(dim);
    free_sizes_1d.reserve(dim);
    for (int d = 0; d < dim; d++) {
      full_sizes_1d.push_back(uni.full_sizes_1d[d]);
      free_sizes_1d.push_back(uni.K[d].Height());
    }
  } else if (prec_name == "jacobi") {
    MFEM_VERIFY(A_sp != nullptr,
                "Jacobi preconditioner requires SparseMatrix A.");
    const Clock::time_point t0 = Clock::now();
    prec_ptr = std::make_unique<DSmoother>(*A_sp);
    const Clock::time_point t1 = Clock::now();
    t_prec_setup_total = ElapsedSec(t0, t1);
  } else if (prec_name == "gs") {
    MFEM_VERIFY(A_sp != nullptr,
                "Gauss-Seidel preconditioner requires SparseMatrix A.");
    const Clock::time_point t0 = Clock::now();
    prec_ptr = std::make_unique<GSSmoother>(*A_sp);
    const Clock::time_point t1 = Clock::now();
    t_prec_setup_total = ElapsedSec(t0, t1);
  }

  // 7. Krylov solver configuration and repeated timed solves.
  ResidualRecorderMonitor res_monitor;
  std::unique_ptr<IterativeSolver> solver;
  if (solver_name == "cg") {
    solver = std::make_unique<CGSolver>();
  } else if (solver_name == "gmres") {
    auto gmres = std::make_unique<GMRESSolver>();
    gmres->SetKDim(gmres_restart);
    solver = std::move(gmres);
  } else // solver_name == "minres"
  {
    solver = std::make_unique<MINRESSolver>();
  }

  solver->SetRelTol(rtol);
  solver->SetAbsTol(atol);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(quiet ? 0 : 1);
  solver->SetOperator(*A);
  if (prec_ptr) {
    solver->SetPreconditioner(*prec_ptr);
  }

  const Vector X_init(X);
  std::vector<double> solve_repeats;
  solve_repeats.reserve(static_cast<std::size_t>(repeats));

  for (int rep = 0; rep < repeats; rep++) {
    X = X_init;
    if (record_residuals && rep == repeats - 1) {
      res_monitor.Reset();
      solver->SetMonitor(res_monitor);
    }
    const Clock::time_point ts0 = Clock::now();
    solver->Mult(B, X);
    const Clock::time_point ts1 = Clock::now();
    solve_repeats.push_back(ElapsedSec(ts0, ts1));
  }

  if (assembly == "legacy") {
    a->RecoverFEMSolution(X, b, x);
  } else // assembly == "wq"
  {
    x.SetFromTrueDofs(X);
  }

  const bool converged = static_cast<bool>(solver->GetConverged());
  const int iterations = solver->GetNumIterations();
  const double initial_norm = static_cast<double>(solver->GetInitialNorm());
  const double final_norm = static_cast<double>(solver->GetFinalNorm());
  const double avg_reduction_factor =
      (iterations > 0 && initial_norm > 0.0 && final_norm > 0.0)
          ? std::pow(final_norm / initial_norm,
                     1.0 / static_cast<double>(iterations))
          : 0.0;

  std::vector<double> sorted_repeats = solve_repeats;
  std::sort(sorted_repeats.begin(), sorted_repeats.end());
  const std::size_t n_rep = sorted_repeats.size();
  const double t_solve =
      (n_rep % 2 == 1)
          ? sorted_repeats[n_rep / 2]
          : 0.5 * (sorted_repeats[n_rep / 2 - 1] + sorted_repeats[n_rep / 2]);
  const double t_solve_per_iter =
      (iterations > 0) ? (t_solve / static_cast<double>(iterations)) : 0.0;
  const double t_total = t_mesh_load + t_refinement + t_fespace_setup +
                         t_operator_assembly + t_rhs_assembly +
                         t_prec_setup_total + t_solve;

  const long peak_rss_kb = GetPeakRssKb();

  // 8. Optional preconditioner verification (excluded from pipeline timings).
  bool has_prec_verif = false;
  double prec_sym_err = 0.0;
  double prec_stale_y_err = 0.0;
  if (verify && prec_ptr) {
    Vector u_p(prec_ptr->Width());
    Vector v_p(prec_ptr->Width());
    Vector pu(prec_ptr->Height());
    Vector pv(prec_ptr->Height());
    Vector pu_dirty(prec_ptr->Height());
    u_p.Randomize(3);
    v_p.Randomize(4);

    pu = 0.0;
    prec_ptr->Mult(u_p, pu);
    pu_dirty = 1e6;
    prec_ptr->Mult(u_p, pu_dirty);
    pu_dirty -= pu;

    pv = 0.0;
    prec_ptr->Mult(v_p, pv);

    prec_stale_y_err = static_cast<double>(pu_dirty.Norml2());
    prec_sym_err = static_cast<double>(std::abs(pu * v_p - u_p * pv));
    has_prec_verif = true;

    if (!quiet) {
      mfem::out << "Prec stale-y err: " << prec_stale_y_err << "\n"
                << "Prec sym err    : " << prec_sym_err << "\n";
    }
  }

  if (!quiet) {
    mfem::out << "Solver: " << solver_name << ", prec: " << prec_name
              << ", converged: " << (converged ? 1 : 0)
              << ", iterations: " << iterations
              << ", final norm: " << final_norm
              << ", median solve (s): " << t_solve << "\n";
  }

  // 9. Write JSON output if requested.
  if (!out_path.empty()) {
    std::ofstream ofs(out_path);
    MFEM_VERIFY(ofs.is_open(), "Failed to open output JSON file: " << out_path);

    ofs << "{\n"
        << "  \"schema_version\": \"1.0\",\n"
        << "  \"run_id\": " << JsonString(run_id_str) << ",\n"
        << "  \"status\": \"ok\",\n"
        << "  \"config\": {\n"
        << "    \"mesh_name\": " << JsonString(mesh_name) << ",\n"
        << "    \"mesh\": " << JsonString(mesh_path) << ",\n"
        << "    \"problem\": " << JsonString(problem) << ",\n"
        << "    \"ref_strategy\": " << JsonString(ref_strategy) << ",\n"
        << "    \"requested_order\": " << requested_order << ",\n"
        << "    \"order\": " << final_order << ",\n"
        << "    \"ref_levels\": " << ref_levels << ",\n"
        << "    \"assembly\": " << JsonString(assembly) << ",\n"
        << "    \"wq_symmetrize\": " << JsonBool(wq_symmetrize) << ",\n"
        << "    \"wq_interior_pts\": " << wq_interior_pts << ",\n"
        << "    \"solver\": " << JsonString(solver_name) << ",\n"
        << "    \"prec\": " << JsonString(prec_name) << ",\n"
        << "    \"rtol\": " << JsonNumber(static_cast<double>(rtol)) << ",\n"
        << "    \"atol\": " << JsonNumber(static_cast<double>(atol)) << ",\n"
        << "    \"max_iter\": " << max_iter << ",\n"
        << "    \"gmres_restart\": " << gmres_restart << ",\n"
        << "    \"repeats\": " << repeats << ",\n"
        << "    \"verify\": " << JsonBool(verify) << ",\n"
        << "    \"record_residuals\": " << JsonBool(record_residuals) << "\n"
        << "  },\n"
        << "  \"problem_stats\": {\n"
        << "    \"dim\": " << dim << ",\n"
        << "    \"base_order\": " << base_order << ",\n"
        << "    \"order\": " << final_order << ",\n"
        << "    \"num_elements\": " << num_elements << ",\n"
        << "    \"num_tdofs\": " << num_tdofs << ",\n"
        << "    \"num_free_dofs\": " << num_free_dofs << ",\n"
        << "    \"num_ess_dofs\": " << num_ess_dofs << ",\n"
        << "    \"full_sizes_1d\": " << JsonIntArray(full_sizes_1d) << ",\n"
        << "    \"free_sizes_1d\": " << JsonIntArray(free_sizes_1d) << ",\n"
        << "    \"matrix_nnz\": " << matrix_nnz << "\n"
        << "  },\n"
        << "  \"timings_sec\": {\n"
        << "    \"mesh_load\": " << JsonNumber(t_mesh_load) << ",\n"
        << "    \"refinement\": " << JsonNumber(t_refinement) << ",\n"
        << "    \"fespace_setup\": " << JsonNumber(t_fespace_setup) << ",\n"
        << "    \"operator_assembly\": " << JsonNumber(t_operator_assembly)
        << ",\n"
        << "    \"op_wq_rules_1d\": " << JsonNumber(t_op_wq_rules_1d) << ",\n"
        << "    \"op_wq_weights_1d\": " << JsonNumber(t_op_wq_weights_1d)
        << ",\n"
        << "    \"op_wq_metric_eval\": " << JsonNumber(t_op_wq_metric_eval)
        << ",\n"
        << "    \"op_wq_patch_assemble\": "
        << JsonNumber(t_op_wq_patch_assemble) << ",\n"
        << "    \"op_wq_symmetrize\": " << JsonNumber(t_op_wq_symmetrize)
        << ",\n"
        << "    \"op_wq_bc_elim\": " << JsonNumber(t_op_wq_bc_elim) << ",\n"
        << "    \"rhs_assembly\": " << JsonNumber(t_rhs_assembly) << ",\n"
        << "    \"prec_setup_total\": " << JsonNumber(t_prec_setup_total)
        << ",\n"
        << "    \"prec_univariate_assembly\": " << JsonNumber(t_prec_univariate)
        << ",\n"
        << "    \"prec_eigendecomp\": " << JsonNumber(t_prec_eigendecomp)
        << ",\n"
        << "    \"prec_fd_construct\": " << JsonNumber(t_prec_fd_construct)
        << ",\n"
        << "    \"solve\": " << JsonNumber(t_solve) << ",\n"
        << "    \"solve_repeats\": " << JsonNumberArray(solve_repeats) << ",\n"
        << "    \"solve_per_iter\": " << JsonNumber(t_solve_per_iter) << ",\n"
        << "    \"total\": " << JsonNumber(t_total) << "\n"
        << "  },\n"
        << "  \"convergence\": {\n"
        << "    \"converged\": " << JsonBool(converged) << ",\n"
        << "    \"iterations\": " << iterations << ",\n"
        << "    \"initial_norm\": " << JsonNumber(initial_norm) << ",\n"
        << "    \"final_norm\": " << JsonNumber(final_norm) << ",\n"
        << "    \"avg_reduction_factor\": " << JsonNumber(avg_reduction_factor)
        << ",\n"
        << "    \"residual_history\": "
        << JsonNumberArray(res_monitor.residuals) << "\n"
        << "  },\n"
        << "  \"verification\": {\n"
        << "    \"enabled\": " << JsonBool(verify) << ",\n"
        << "    \"prec_sym_err\": "
        << JsonNullableNumber(has_prec_verif, prec_sym_err) << ",\n"
        << "    \"prec_stale_y_err\": "
        << JsonNullableNumber(has_prec_verif, prec_stale_y_err) << "\n"
        << "  },\n"
        << "  \"system\": {\n"
        << "    \"peak_rss_kb\": " << peak_rss_kb << "\n"
        << "  }\n"
        << "}\n";

    ofs.flush();
    MFEM_VERIFY(ofs.good(),
                "Failed while writing output JSON file: " << out_path);
  }

  return 0;
}
