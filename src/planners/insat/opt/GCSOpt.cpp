//  * Copyright (c) 2023, Ramkumar Natarajan
//  * All rights reserved.
//  *
//  * Redistribution and use in source and binary forms, with or without
//  * modification, are permitted provided that the following conditions are met:
//  *
//  *     * Redistributions of source code must retain the above copyright
//  *       notice, this list of conditions and the following disclaimer.
//  *     * Redistributions in binary form must reproduce the above copyright
//  *       notice, this list of conditions and the following disclaimer in the
//  *       documentation and/or other materials provided with the distribution.
//  *     * Neither the name of the Carnegie Mellon University nor the names of its
//  *       contributors may be used to endorse or promote products derived from
//  *       this software without specific prior written permission.
//  *
//  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
//  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  * POSSIBILITY OF SUCH DAMAGE.
//

/*!
 * \file opt.cpp
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 9/3/23
*/

#include <planners/insat/opt/GCSOpt.hpp>

#include <iostream>

namespace ps {

// Given a list of matrices, return the matrices with every column where all
// of the matrices are zero in that column, along with a boolean vector
// indicating which columns were preserved (true) or removed (false).
  std::tuple<std::vector<Eigen::MatrixXd>, VectorXb> CondenseToNonzeroColumns(
          std::vector<Eigen::MatrixXd> matrices) {
    // Validate inputs.
    DRAKE_DEMAND(matrices.size() > 0);
    const int num_cols = matrices[0].cols();
    for (const Eigen::MatrixXd& matrix : matrices) {
      DRAKE_DEMAND(matrix.cols() == num_cols);
    }

    // Find non-zero columns.
    VectorXb nonzero_cols_mask = VectorXb::Constant(num_cols, false);
    for (const Eigen::MatrixXd& matrix : matrices) {
      nonzero_cols_mask += matrix.cast<bool>().colwise().any();
    }
    const int nonzero_cols_count = nonzero_cols_mask.count();

    // Create the output, copying only the non-zero columns.
    std::vector<Eigen::MatrixXd> condensed_matrices;
    for (const Eigen::MatrixXd& matrix : matrices) {
      Eigen::MatrixXd& condensed_matrix =
              condensed_matrices.emplace_back(matrix.rows(), nonzero_cols_count);
      int condensed_col = 0;
      for (int orig_col = 0; orig_col < matrix.cols(); ++orig_col) {
        if (nonzero_cols_mask(orig_col)) {
          condensed_matrix.col(condensed_col) = matrix.col(orig_col);
          condensed_col++;
        }
      }
    }
    return std::make_tuple(condensed_matrices, nonzero_cols_mask);
  }

// Filters variables given a vector of variables along with a boolean vector
// indicating which rows were preserved (true) or removed (false).
  Eigen::VectorX<drake::symbolic::Variable> FilterVariables(
          const Eigen::VectorX<drake::symbolic::Variable>& vars,
          const VectorXb& nonzero_cols_mask) {
    Eigen::VectorX<drake::symbolic::Variable> vars_dense(nonzero_cols_mask.count());
    int row = 0;
    for (int i = 0; i < vars.size(); ++i) {
      if (nonzero_cols_mask(i)) {
        vars_dense(row++) = vars(i);
      }
    }
    return vars_dense;
  }
}

ps::GCSOpt::GCSOpt(const std::vector<HPolyhedron> &regions,
                   const std::vector<std::pair<int, int>> &edges_between_regions,
                   int order, double h_min, double h_max,
                   double path_length_weight, double time_weight,
                   Eigen::VectorXd& vel_lb, Eigen::VectorXd& vel_ub,
                   bool verbose)
        : verbose_(verbose),
          hpoly_regions_(regions),
          edges_bw_regions_(edges_between_regions),
          order_(order),
          h_min_(h_min),
          h_max_(h_max),
          vel_lb_(vel_lb),
          vel_ub_(vel_ub),
          enable_time_cost_(false),
          enable_path_length_cost_(false),
          enable_path_velocity_constraint_(false),
          gcs_(std::make_shared<drake::geometry::optimization::GraphOfConvexSets>()) {

  drake::geometry::optimization::ConvexSets regions_cs;
  for (const auto& region : regions) {
    auto cs = MakeConvexSets(region);
    regions_cs.push_back(cs[0]);
  }

  num_positions_ = regions_cs[0]->ambient_dimension();
  if (time_weight != 0) {
//  if (true) { /// FIXME: Temporary hotfix by Ram. Setting enable_time_cost_=false is segfaulting.
    time_weight_ = time_weight;
    enable_time_cost_ = true;
  }
  if (path_length_weight != 0) {
    auto path_length_weight_matrix = path_length_weight * Eigen::MatrixXd::Identity(num_positions_, num_positions_);
    path_length_weight_.resize(num_positions_, 2*num_positions_);
    path_length_weight_.block(0, 0, num_positions_, num_positions_) = path_length_weight_matrix;
    path_length_weight_.block(0, num_positions_, num_positions_, num_positions_) = -path_length_weight_matrix;

    enable_path_length_cost_ = true;
  }
  enable_path_velocity_constraint_ = true;

  if (verbose_) std::cout << "Setting up vars!" << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  setupVars();
  auto end_time = std::chrono::high_resolution_clock::now();
  if (verbose_) std::cout << "Done setting up vars!" << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;

  if (verbose_) std::cout << "Preprocessing" << std::endl;
  start_time = std::chrono::high_resolution_clock::now();
  preprocess(regions_cs, edges_between_regions);
  end_time = std::chrono::high_resolution_clock::now();
  if (verbose_) std::cout << "Done setting up vars and preprocessing!!!" << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;
}

void ps::GCSOpt::FormulateAndSetCostsAndConstraints() {

  if (!enable_path_length_cost_ && !enable_time_cost_) {
    std::runtime_error("The problem should be minimum length or minimum time or both. Nothing is enabled now!!");
  }

  if (verbose_) std::cout << "Formulating costs and constraints" << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  formulateCostsAndConstraints();
  auto end_time = std::chrono::high_resolution_clock::now();
  if (verbose_) std::cout << "Done formulating costs and constraints!!!" << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;


  if (verbose_) std::cout << "Setting up costs and constraints" << std::endl;
  start_time = std::chrono::high_resolution_clock::now();
  setupCostsAndConstraints();
  end_time = std::chrono::high_resolution_clock::now();
  if (verbose_) std::cout << "Done setting up  costs and constraints!!!" << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;

}

ps::VertexId ps::GCSOpt::AddStart(Eigen::VectorXd &start) {

  auto start_set = MakeConvexSets(drake::geometry::optimization::Point(start))[0];

  // Make time scaling set once to avoid many allocations when adding the
  // vertices to GCS.
  const drake::geometry::optimization::HPolyhedron time_scaling_set =
          drake::geometry::optimization::HPolyhedron::MakeBox(
                  drake::Vector1d(h_min_), drake::Vector1d(h_max_));

  // Add Regions with time scaling set.
  drake::geometry::optimization::ConvexSets vertex_set;
  // Assign each control point to a separate set.
  const int num_points = order_ + 1;
  vertex_set.reserve(num_points + 1);
  vertex_set.insert(vertex_set.begin(), num_points,
                    drake::geometry::optimization::ConvexSets::value_type{start_set});
  // Add time scaling set.
  if (enable_time_cost_) {
    vertex_set.emplace_back(time_scaling_set);
  }

  vertices_.emplace_back(gcs_->AddVertex(
          drake::geometry::optimization::CartesianProduct(vertex_set),
          fmt::format("{}", "start")));
  vertex_id_to_vertex_[vertices_.back()->id().get_value()-1] = vertices_.back();
  start_vtx_ = vertices_.back();
  start_vit_ = vertices_.end()-1;

  auto start_vertex = vertices_.back();
  auto upstart = start_vertex->set().MaybeGetFeasiblePoint();
  GCSVertex* start_region_vertex;
  for (auto& v : vertices_) {
    if (v->set().PointInSet(upstart.value())) {
      start_region_vertex = v;
      break;
    }
  }
  // Connect start to start region
  GCSEdge* uv_edge = gcs_->AddEdge(start_vertex, start_region_vertex);
  edges_.emplace_back(uv_edge);
  edge_id_to_edge_[edges_.back()->id().get_value()-1] = edges_.back();
  start_eit_ = edges_.end()-1;

  GCSEdge* vu_edge = gcs_->AddEdge(start_region_vertex, start_vertex);
  edges_.emplace_back(vu_edge);
  edge_id_to_edge_[edges_.back()->id().get_value()-1] = edges_.back();

  return start_vertex->id();
}

ps::VertexId ps::GCSOpt::AddGoal(Eigen::VectorXd &goal) {

  auto goal_set = MakeConvexSets(drake::geometry::optimization::Point(goal))[0];

  // Make time scaling set once to avoid many allocations when adding the
  // vertices to GCS.
  const drake::geometry::optimization::HPolyhedron time_scaling_set =
          drake::geometry::optimization::HPolyhedron::MakeBox(
                  drake::Vector1d(h_min_), drake::Vector1d(h_max_));

  // Add Regions with time scaling set.
  drake::geometry::optimization::ConvexSets vertex_set;
  // Assign each control point to a separate set.
  const int num_points = order_ + 1;
  vertex_set.reserve(num_points + 1);
  vertex_set.insert(vertex_set.begin(), num_points,
                    drake::geometry::optimization::ConvexSets::value_type{goal_set});
  // Add time scaling set.
  if (enable_time_cost_) {
    vertex_set.emplace_back(time_scaling_set);
  }

  vertices_.emplace_back(gcs_->AddVertex(
          drake::geometry::optimization::CartesianProduct(vertex_set),
          fmt::format("{}", "goal")));
  vertex_id_to_vertex_[vertices_.back()->id().get_value()-1] = vertices_.back();
  goal_vtx_ = vertices_.back();
  goal_vit_ = vertices_.end()-1;

  auto goal_vertex = vertices_.back();
  auto upgoal = goal_vertex->set().MaybeGetFeasiblePoint();
  GCSVertex* goal_region_vertex;
  for (auto& v : vertices_) {
    if (v->set().PointInSet(upgoal.value())) {
      goal_region_vertex = v;
      break;
    }
  }
  // Connect goal region to goal
  GCSEdge* uv_edge = gcs_->AddEdge(goal_region_vertex, goal_vertex);
  edges_.emplace_back(uv_edge);
  edge_id_to_edge_[edges_.back()->id().get_value()-1] = edges_.back();
  goal_eit_ = edges_.end()-1;

  GCSEdge* vu_edge = gcs_->AddEdge(goal_vertex, goal_region_vertex);
  edges_.emplace_back(vu_edge);
  edge_id_to_edge_[edges_.back()->id().get_value()-1] = edges_.back();

  return goal_vertex->id();
}

double ps::GCSOpt::LowerboundSolve(const std::vector<int>& path_ids) {
  std::vector<VertexId> path_vids;
  for (const auto& id : path_ids) {
    path_vids.push_back(vertex_id_to_vertex_[id]->id());
  }
  auto lb_soln = Solve(path_vids);
  auto lb_traj = lb_soln.first;
  auto p0 = lb_traj.value(lb_traj.start_time());
  double dist = 0;
  for (double t=lb_traj.start_time(); t<=lb_traj.end_time(); t+=1e-1) {
    auto pt = lb_traj.value(t);
    dist += (pt-p0).norm();
    p0 = pt;
  }
  return dist;
}

std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult> ps::GCSOpt::Solve(std::vector<VertexId> &path_vids) {
  Eigen::VectorXd dummy_init_guess;
  assert(dummy_init_guess.size() == 0);
  return Solve(path_vids, dummy_init_guess);
}

std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult> ps::GCSOpt::Solve(std::vector<VertexId> &path_vids,
                                                                     Eigen::VectorXd& initial_guess) {
  std::vector<EdgeId> path_eids;
  for (int i=0; i<path_vids.size()-1; ++i) {
    int64_t uid = path_vids[i].get_value();
    int64_t vid = path_vids[i+1].get_value();
    for (auto& e : edges_) {
      if (e->u().id().get_value() == uid && e->v().id().get_value() == vid) {
        path_eids.push_back(e->id());
        break;
      }
    }
  }
  return Solve(path_vids, path_eids, initial_guess);
}

std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult>
ps::GCSOpt::Solve(std::vector<VertexId> &path_vids,
                  std::vector<EdgeId> &path_eids) {
  Eigen::VectorXd dummy_init_guess;
  assert(dummy_init_guess.size() == 0);
  return Solve(path_vids, path_eids, dummy_init_guess);
}

std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult>
ps::GCSOpt::Solve(std::vector<VertexId>& path_vids,
                  std::vector<EdgeId>& path_eids,
                  Eigen::VectorXd& initial_guess) {

  if (path_vids.size() <= 0) {
    std::runtime_error("Size of Path IDs has to be positive!!");
  }

  drake::solvers::MathematicalProgram prog;

  for (const auto& vid : path_vids) {
    if (vertex_id_to_vertex_.find(vid.get_value()) == vertex_id_to_vertex_.end()) {
      std::runtime_error("Vertex ID: " + std::to_string(vid.get_value()) + " not found in the map!!");
    }
    auto& dec_vars = vertex_id_to_vertex_[vid.get_value()]->x();
    prog.AddDecisionVariables(dec_vars);
    for (auto& slack : slack_vars_) {
      prog.AddDecisionVariables(slack);
    }
    for (const auto& cost : vertex_id_to_cost_binding_[vid.get_value()]) {
      prog.AddCost(cost);
    }
    for (const auto& constraint : vertex_id_to_constraint_binding_[vid.get_value()]) {
      prog.AddConstraint(constraint);
    }
    vertex_id_to_vertex_[vid.get_value()]->set().
            AddPointInSetConstraints(&prog, dec_vars); // except the last var which is time scaling
  }

  for (const auto& eid : path_eids) {
//    auto& edge = edge_id_to_edge_[eid.get_value()];
//    const Eigen::VectorX<drake::symbolic::Variable> edge_vars =
//            drake::solvers::ConcatenateVariableRefList({edge->xu(), edge->xv()});
//    prog.AddDecisionVariables(edge_vars);
    for (const auto& constraint : edge_id_to_constraint_binding_[eid.get_value()]) {
      prog.AddConstraint(constraint);
    }
  }

  if (verbose_) {
    std::cout << "math program " << std::endl;
    std::cout << prog.to_string() << std::endl;
  }

/// Default is using SNOPT
////  drake::solvers::MathematicalProgramResult result = drake::solvers::Solve(prog);
//  auto start_time = std::chrono::high_resolution_clock::now();
//  drake::solvers::MathematicalProgramResult result;
//  if (initial_guess.size() == 0) {
//    result = drake::solvers::Solve(prog);
//  } else if (initial_guess.size() < prog.num_vars()) {
//    Eigen::VectorXd full_init_guess(prog.num_vars());
//    full_init_guess.setZero();
//    full_init_guess.head(initial_guess.size()) = initial_guess;
//    result = drake::solvers::Solve(prog, full_init_guess);
//  } else {
//    result = drake::solvers::Solve(prog, initial_guess);
//  }
//  auto end_time = std::chrono::high_resolution_clock::now();

/// Use MOSEK
//  RewriteForConvexSolver(&prog);
  auto start_time = std::chrono::high_resolution_clock::now();
  auto mosek_solver = drake::solvers::MosekSolver();
  drake::solvers::MathematicalProgramResult result;
  if (initial_guess.size() == 0) {
    result = mosek_solver.Solve(prog);
  } else if (initial_guess.size() < prog.num_vars()) {
    Eigen::VectorXd full_init_guess(prog.num_vars());
    full_init_guess.setZero();
    full_init_guess.head(initial_guess.size()) = initial_guess;
    result = mosek_solver.Solve(prog, full_init_guess);
  } else {
    result = mosek_solver.Solve(prog, initial_guess);
  }
  auto end_time = std::chrono::high_resolution_clock::now();

  if (verbose_)  std::cout << "Solving alone took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;

  if (!result.is_success()) {
    return {drake::trajectories::CompositeTrajectory<double>({}), result};
  }

// Extract the path from the edges.
  std::vector<drake::copyable_unique_ptr<drake::trajectories::Trajectory<double>>> bezier_curves;
  for (const auto& id : path_vids) {
    auto& vertex = vertex_id_to_vertex_[id.get_value()];
    const int num_control_points = order_ + 1;
    const Eigen::MatrixX<double> path_points =
            Eigen::Map<Eigen::MatrixX<double>>(result.GetSolution(vertex->x()).data(),
                                               num_positions_, num_control_points);

    double h;
    if (enable_time_cost_) {
// Extract the duration from the solution.
      h = result.GetSolution(vertex->x()).tail<1>().value();
    } else {
      h = 1;
    }
    const double start_time = bezier_curves.empty() ? 0 : bezier_curves.back()->end_time();

// Skip edges with a single control point that spend near zero time in the
// region, since zero order continuity constraint is sufficient. These edges
// would result in a discontinuous trajectory for velocities and higher
// derivatives.
    if (!(num_control_points == 1 && h_min_ == 0)) {
      bezier_curves.emplace_back(std::make_unique<drake::trajectories::BezierCurve<double>>(
              start_time, start_time + h, path_points));
    }
  }

  auto final_trajectory = drake::trajectories::CompositeTrajectory<double>(bezier_curves);

  return {final_trajectory, result};
}

void ps::GCSOpt::CleanUp() {
  gcs_->RemoveVertex(start_vtx_);
  gcs_->RemoveVertex(goal_vtx_);
  vertices_.erase(start_vit_);
  vertices_.erase(goal_vit_);
  edges_.erase(start_eit_);
  edges_.erase(goal_eit_);
}

void ps::GCSOpt::setupVars() {
  const Eigen::MatrixX<drake::symbolic::Variable> u_control =
          drake::symbolic::MakeMatrixContinuousVariable(
                  num_positions_, order_ + 1, "xu");
  const Eigen::MatrixX<drake::symbolic::Variable> v_control =
          drake::symbolic::MakeMatrixContinuousVariable(
                  num_positions_, order_ + 1, "xv");
  Eigen::Map<const Eigen::VectorX<drake::symbolic::Variable>> u_control_vars(
          u_control.data(), u_control.size());
  Eigen::Map<const Eigen::VectorX<drake::symbolic::Variable>> v_control_vars(
          v_control.data(), v_control.size());

  if (enable_time_cost_) {
    u_h_ = drake::symbolic::MakeVectorContinuousVariable(1, "Tu");
    v_h_ = drake::symbolic::MakeVectorContinuousVariable(1, "Tv");

    u_vars_ = drake::solvers::ConcatenateVariableRefList({u_control_vars, u_h_});
    v_vars_ = drake::solvers::ConcatenateVariableRefList({v_control_vars, v_h_});
  } else {
    u_vars_ = u_control_vars;
    v_vars_ = v_control_vars;
  }

//  const Eigen::VectorX<drake::symbolic::Variable> edge_vars =
//      drake::solvers::ConcatenateVariableRefList(
//          {u_control_vars, u_h_, v_control_vars, v_h_});

  u_r_trajectory_ = drake::trajectories::BezierCurve<drake::symbolic::Expression>(
          0, 1, u_control.cast<drake::symbolic::Expression>());

  v_r_trajectory_ = drake::trajectories::BezierCurve<drake::symbolic::Expression>(
          0, 1, v_control.cast<drake::symbolic::Expression>());
}

void ps::GCSOpt::preprocess(const drake::geometry::optimization::ConvexSets &regions,
                            const std::vector<std::pair<int, int>> &edges_between_regions) {
  // Make time scaling set once to avoid many allocations when adding the
  // vertices to GCS.
  const drake::geometry::optimization::HPolyhedron time_scaling_set =
          drake::geometry::optimization::HPolyhedron::MakeBox(
                  drake::Vector1d(h_min_), drake::Vector1d(h_max_));

  // Add Regions with time scaling set.
  for (size_t i = 0; i < regions.size(); ++i) {
    drake::geometry::optimization::ConvexSets vertex_set;
    // Assign each control point to a separate set.
    const int num_points = order_ + 1;
    vertex_set.reserve(num_points + 1);
    vertex_set.insert(vertex_set.begin(), num_points,
                      drake::geometry::optimization::ConvexSets::value_type{regions[i]});
    // Add time scaling set.
    if (enable_time_cost_) {
      vertex_set.emplace_back(time_scaling_set);
    }

    vertices_.emplace_back(gcs_->AddVertex(
            drake::geometry::optimization::CartesianProduct(vertex_set),
            fmt::format("{}: {}", "v" + std::to_string(i), i)));
    std::cout << "vid: " << vertices_.back()->id().get_value()-1 << " i " << i << std::endl;
    vertex_id_to_vertex_[vertices_.back()->id().get_value()-1] = vertices_.back();
    vertex_id_to_regions_[vertices_.back()->id().get_value()-1] = hpoly_regions_[i];

//    std::cout << "Added vertex with id: " << vertices_.back()->id().get_value()-1 << std::endl;
  }

  // Connect vertices with edges.
  for (const auto& [u_index, v_index] : edges_between_regions) {
    // Add edge.
    GCSVertex* u = vertices_[u_index];
    GCSVertex* v = vertices_[v_index];
    GCSEdge* uv_edge = gcs_->AddEdge(u, v);

    edges_.emplace_back(uv_edge);
    edge_id_to_edge_[edges_.back()->id().get_value()-1] = edges_.back();
  }

}

void ps::GCSOpt::formulateTimeCost() {
  // The time cost is the sum of duration variables ∑ hᵢ
  time_cost_ =
          std::make_shared<drake::solvers::LinearCost>(time_weight_ * Eigen::VectorXd::Ones(1), 0.0);
}

void ps::GCSOpt::formulatePathLengthCost() {
  /*
    We will upper bound the trajectory length by the sum of the distances
    between the control points. ∑ |weight_matrix * (rᵢ₊₁ − rᵢ)|₂
  */

  if (order_ == 0) {
    throw std::runtime_error(
        "Path length cost is not defined for a set of order 0.");
  }

  if (verbose_) {
    std::cout << "path_length_weight_: \n" << path_length_weight_ << std::endl;
  }

  const auto path_length_cost =
      std::make_shared<drake::solvers::L2NormCost>(path_length_weight_, Eigen::VectorXd::Zero(num_positions_));
  path_length_cost_.emplace_back(path_length_cost);
}

void ps::GCSOpt::formulatePathContinuityConstraint() {
  const Eigen::VectorX<drake::symbolic::Variable> edge_vars =
          drake::solvers::ConcatenateVariableRefList({u_vars_, v_vars_});

  const Eigen::VectorX<drake::symbolic::Expression> path_continuity_error =
          v_r_trajectory_.control_points().col(0) -
          u_r_trajectory_.control_points().col(order_);
  Eigen::MatrixXd M(num_positions_, edge_vars.size());
  drake::symbolic::DecomposeLinearExpressions(path_continuity_error, edge_vars, &M);
  // Condense M to only keep non-zero columns.
  const auto& [condensed_matrices, nonzero_cols_mask] =
          CondenseToNonzeroColumns({M});
  Eigen::MatrixXd M_dense = condensed_matrices[0];

  path_continuity_constraint_ = std::make_pair(
          std::make_shared<drake::solvers::LinearEqualityConstraint>(
                  M_dense, Eigen::VectorXd::Zero(num_positions_)), nonzero_cols_mask);

  if (verbose_) {
    std::cout << "path_continuity_error:\n" << path_continuity_error << std::endl;
    std::cout << "M:\n" << M << std::endl;
    std::cout << "u_r_trajectory_.control_points():\n" << u_r_trajectory_.control_points() << std::endl;
    std::cout << "v_r_trajectory_.control_points():\n" << v_r_trajectory_.control_points() << std::endl;
    std::cout << "u_vars_:\n" << u_vars_ << std::endl;
    std::cout << "v_vars_:\n" << v_vars_ << std::endl;
    std::cout << "edge_vars:\n" << edge_vars << std::endl;
    std::cout << "nonzero_cols_mask:\n" << nonzero_cols_mask << std::endl;
    std::cout << "M_dense:\n" << M_dense << std::endl;
  }
}

void ps::GCSOpt::formulateVelocityConstraint() {
//  const Eigen::MatrixX<drake::symbolic::Expression> u_rdot_control =
//          drake::dynamic_pointer_cast_or_throw<drake::trajectories::BezierCurve<drake::symbolic::Expression>>(
//                  u_r_trajectory_.MakeDerivative())
//                  ->control_points();
//
//  Eigen::MatrixXd b(u_h_.rows(), u_vars_.size());
//  DecomposeLinearExpressions(u_h_.cast<drake::symbolic::Expression>(), u_vars_, &b);
//
//  for (int i = 0; i < u_rdot_control.cols(); ++i) {
//    /// TEMPORARY (to reproduce the results of GCSTrajectoryOptimization)
//    if (i!=0 && i!=u_rdot_control.cols()-1) {
//      continue;
//    }
//
//    Eigen::MatrixXd M(num_positions_, u_vars_.size());
//    DecomposeLinearExpressions(u_rdot_control.col(i), u_vars_, &M);
//    // Condense M and b to only keep non-zero columns.
//    const auto& [condensed_matrices, nonzero_cols_mask] =
//            CondenseToNonzeroColumns({M, b});
//    Eigen::MatrixXd M_dense = condensed_matrices[0];
//    Eigen::MatrixXd b_dense = condensed_matrices[1];
//
//    Eigen::MatrixXd H(2 * num_positions_, nonzero_cols_mask.count());
//    H << M_dense - vel_ub_ * b_dense, -M_dense + vel_lb_ * b_dense;
//
//    velocity_constraint_.emplace_back(std::make_shared<drake::solvers::LinearConstraint>(
//                                              H, Eigen::VectorXd::Constant(2 * num_positions_, -kInf),
//                                              Eigen::VectorXd::Zero(2 * num_positions_)),
//                                      nonzero_cols_mask);
//  }
}

void ps::GCSOpt::formulateCostsAndConstraints() {
  if (enable_time_cost_) {
    formulateTimeCost();
  }

  if (enable_path_length_cost_) {
    formulatePathLengthCost();
  }

  formulatePathContinuityConstraint();

  if (enable_path_velocity_constraint_) {
    formulateVelocityConstraint();
  }
}

void ps::GCSOpt::addCosts(const drake::geometry::optimization::GraphOfConvexSets::Vertex *v) {
  if (enable_time_cost_) {
    CostBinding time_cost_binding =
            drake::solvers::Binding<drake::solvers::Cost>(
                    time_cost_, v->x().tail(1));

    vertex_id_to_cost_binding_[v->id().get_value()-1].emplace_back(time_cost_binding);
  }

  if (enable_path_length_cost_) {
    for (const auto& c : path_length_cost_) {
      auto control_points = Eigen::Map<const Eigen::MatrixX<drake::symbolic::Variable>>(
          v->x().data(), num_positions_, order_ + 1);

      for (int i = 0; i < control_points.cols() - 1; ++i) {
        CostBinding path_length_cost_binding =
            drake::solvers::Binding<drake::solvers::Cost>(
        c, {control_points.col(i + 1), control_points.col(i)});
        vertex_id_to_cost_binding_[v->id().get_value()-1].emplace_back(path_length_cost_binding);
      }
    }
  }
}

void ps::GCSOpt::addConstraints(const drake::geometry::optimization::GraphOfConvexSets::Vertex *v) {
  if (enable_path_velocity_constraint_) {
//    for (const auto& vc : velocity_constraint_) {
//      ConstraintBinding velocity_limits_binding =
//              drake::solvers::Binding<drake::solvers::Constraint>(
//                      vc.first, FilterVariables(v->x(), vc.second));
//
//      vertex_id_to_constraint_binding_[v->id().get_value()-1].emplace_back(velocity_limits_binding);
//    }
  }

//  if (start_.size() != 0 && v==vertex_id_to_vertex_[start_vtx_.get_value()-1]) {
//    ConstraintBinding start_point_binding =
//            drake::solvers::Binding<drake::solvers::Constraint>(
//                    start_point_constraint_.first,
//                    FilterVariables(v->x(), start_point_constraint_.second));
//
//    vertex_id_to_constraint_binding_[start_vtx_.get_value()-1].emplace_back(start_point_binding);
//  }
}

void ps::GCSOpt::addConstraints(const drake::geometry::optimization::GraphOfConvexSets::Edge *e) {
  ConstraintBinding continuity_constraint_binding =
          drake::solvers::Binding<drake::solvers::Constraint>(
                  path_continuity_constraint_.first,
                  FilterVariables(drake::solvers::ConcatenateVariableRefList({e->u().x(), e->v().x()}),
                                  path_continuity_constraint_.second));
//        std::cout << "continuity_constraint_binding: \n" << continuity_constraint_binding << std::endl;

  edge_id_to_constraint_binding_[e->id().get_value()-1].emplace_back(continuity_constraint_binding);
}

void ps::GCSOpt::setupCostsAndConstraints() {
  for (const auto* v : vertices_) {
    addCosts(v);
    addConstraints(v);
  }
  for (const auto* e : edges_) {
    addConstraints(e);
  }
  RewriteForConvexSolver();
}

/* Most convex solvers require only support linear and quadratic costs when
operating with nonlinear constraints. This removes costs and adds variables and
constraints as needed by the solvers. */
void ps::GCSOpt::RewriteForConvexSolver() {
  drake::solvers::MathematicalProgram prog;

  const double kInf = std::numeric_limits<double>::infinity();
  for (auto& vtx_costs : vertex_id_to_cost_binding_) {
    std::vector<CostBinding> new_costs;
    for (auto& cst_bind : vtx_costs.second) {
      const drake::solvers::Cost* cst = cst_bind.evaluator().get();

      if (const auto* l2c = dynamic_cast<const drake::solvers::L2NormCost*>(cst)) {
        const int m = l2c->A().rows();
        const int n = l2c->A().cols();
        auto slack = drake::symbolic::MakeVectorContinuousVariable(1, "slack");
        slack_vars_.emplace_back(slack);
        auto new_cost = std::make_shared<drake::solvers::LinearCost>(drake::Vector1d::Ones());
        auto new_cost_binding =
            drake::solvers::Binding<drake::solvers::Cost>(new_cost, slack);
        new_costs.emplace_back(new_cost_binding);
        // |Ax+b|² ≤ slack, written as a Lorentz cone with z = [slack; Ax+b].
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m + 1, n + 1);
        A(0, 0) = 1;
        A.bottomRightCorner(m, n) = l2c->A();
        Eigen::VectorXd b(m + 1);
        b << 0, l2c->b();
        auto new_con = std::make_shared<drake::solvers::LorentzConeConstraint>(A, b);
        auto new_con_binding = drake::solvers::Binding<drake::solvers::Constraint>(new_con, {slack, cst_bind.variables()});
        vertex_id_to_constraint_binding_[vtx_costs.first].emplace_back(new_con_binding);
      } else if (const auto* l1c = dynamic_cast<const drake::solvers::L1NormCost*>(cst)) {
        const int m = l1c->A().rows();
        const int n = l1c->A().cols();
        auto slack = drake::symbolic::MakeVectorContinuousVariable(m, "slack");
        slack_vars_.emplace_back(slack);
        auto new_cost = std::make_shared<drake::solvers::LinearCost>(Eigen::VectorXd::Ones(m));
        auto new_cost_binding =
            drake::solvers::Binding<drake::solvers::Cost>(new_cost, slack);
        new_costs.emplace_back(new_cost_binding);
        // Ax + b ≤ slack, written as [A,-I][x;slack] ≤ -b.
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m + n);
        A << l1c->A(), -Eigen::MatrixXd::Identity(m, m);
        vertex_id_to_constraint_binding_[vtx_costs.first].emplace_back(prog.AddLinearConstraint(A, Eigen::VectorXd::Constant(m, -kInf), -l1c->b(),
                                  {cst_bind.variables(), slack}));
        // -(Ax + b) ≤ slack, written as [A,I][x;slack] ≥ -b.
        A.rightCols(m) = Eigen::MatrixXd::Identity(m, m);
        auto new_con = std::make_shared<drake::solvers::LinearConstraint>(A, -l1c->b(), Eigen::VectorXd::Constant(m, kInf));
        auto new_con_binding = drake::solvers::Binding<drake::solvers::Constraint>(new_con, {cst_bind.variables(), slack});
        vertex_id_to_constraint_binding_[vtx_costs.first].emplace_back(new_con_binding);
      } else if (const auto* linfc = dynamic_cast<const drake::solvers::LInfNormCost*>(cst)) {
        const int m = linfc->A().rows();
        const int n = linfc->A().cols();
        auto slack = drake::symbolic::MakeVectorContinuousVariable(1, "slack");
        slack_vars_.emplace_back(slack);
        auto new_cost = std::make_shared<drake::solvers::LinearCost>(drake::Vector1d::Ones());
        auto new_cost_binding =
            drake::solvers::Binding<drake::solvers::Cost>(new_cost, slack);
        new_costs.emplace_back(new_cost_binding);
        // ∀i, aᵢᵀx + bᵢ ≤ slack, written as [A,-1][x;slack] ≤ -b.
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, n + 1);
        A << linfc->A(), Eigen::VectorXd::Constant(m, -1);
        vertex_id_to_constraint_binding_[vtx_costs.first].emplace_back(prog.AddLinearConstraint(A, Eigen::VectorXd::Constant(m, -kInf), -linfc->b(),
                                  {cst_bind.variables(), slack}));
        // ∀i, -(aᵢᵀx + bᵢ) ≤ slack, written as [A,1][x;slack] ≥ -b.
        A.col(A.cols() - 1) = Eigen::VectorXd::Ones(m);
        auto new_con = std::make_shared<drake::solvers::LinearConstraint>(A, -linfc->b(), Eigen::VectorXd::Constant(m, kInf));
        auto new_con_binding = drake::solvers::Binding<drake::solvers::Constraint>(new_con, {cst_bind.variables(), slack});
        vertex_id_to_constraint_binding_[vtx_costs.first].emplace_back(new_con_binding);
      } else if (const auto* pqc = dynamic_cast<const drake::solvers::PerspectiveQuadraticCost*>(cst)) {
        const int m = pqc->A().rows();
        const int n = pqc->A().cols();
        auto slack = drake::symbolic::MakeVectorContinuousVariable(1, "slack");
        slack_vars_.emplace_back(slack);
        auto new_cost = std::make_shared<drake::solvers::LinearCost>(drake::Vector1d::Ones());
        auto new_cost_binding =
            drake::solvers::Binding<drake::solvers::Cost>(new_cost, slack);
        new_costs.emplace_back(new_cost_binding);
        // Written as rotated Lorentz cone with z = [slack; Ax+b].
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m + 1, n + 1);
        A(0, 0) = 1;
        A.bottomRightCorner(m, n) = pqc->A();
        Eigen::VectorXd b(m + 1);
        b << 0, pqc->b();
        auto new_con = std::make_shared<drake::solvers::RotatedLorentzConeConstraint>(A, b);
        auto new_con_binding = drake::solvers::Binding<drake::solvers::Constraint>(new_con, {slack, cst_bind.variables()});
        vertex_id_to_constraint_binding_[vtx_costs.first].emplace_back(new_con_binding);
      } else {
        new_costs.emplace_back(cst_bind);
      }
    }
    vtx_costs.second = new_costs;
  }
}

/* Most convex solvers require only support linear and quadratic costs when
operating with nonlinear constraints. This removes costs and adds variables and
constraints as needed by the solvers. */
void ps::GCSOpt::RewriteForConvexSolver(drake::solvers::MathematicalProgram* prog) {
  // Use Mosek's requirements to test the program attributes.
  drake::solvers::MosekSolver mosek;
  if (mosek.AreProgramAttributesSatisfied(*prog)) {
    return;
  }

  const double kInf = std::numeric_limits<double>::infinity();

  // Loop through all unsupported costs and rewrite them into constraints.
  std::unordered_set<drake::solvers::Binding<drake::solvers::Cost>> to_remove;

  for (const auto& binding : prog->l2norm_costs()) {
    const int m = binding.evaluator()->A().rows();
    const int n = binding.evaluator()->A().cols();
    auto slack = prog->NewContinuousVariables<1>("slack");
    prog->AddLinearCost(drake::Vector1d::Ones(), slack);
    // |Ax+b|² ≤ slack, written as a Lorentz cone with z = [slack; Ax+b].
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m + 1, n + 1);
    A(0, 0) = 1;
    A.bottomRightCorner(m, n) = binding.evaluator()->A();
    Eigen::VectorXd b(m + 1);
    b << 0, binding.evaluator()->b();
    prog->AddLorentzConeConstraint(A, b, {slack, binding.variables()});
    to_remove.insert(binding);
  }

  for (const auto& binding : prog->generic_costs()) {
    const drake::solvers::Cost* cost = binding.evaluator().get();
    if (const auto* l1c = dynamic_cast<const drake::solvers::L1NormCost*>(cost)) {
      const int m = l1c->A().rows();
      const int n = l1c->A().cols();
      auto slack = prog->NewContinuousVariables(m, "slack");
      prog->AddLinearCost(Eigen::VectorXd::Ones(m), slack);
      // Ax + b ≤ slack, written as [A,-I][x;slack] ≤ -b.
      Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m + n);
      A << l1c->A(), -Eigen::MatrixXd::Identity(m, m);
      prog->AddLinearConstraint(A, Eigen::VectorXd::Constant(m, -kInf), -l1c->b(),
                                {binding.variables(), slack});
      // -(Ax + b) ≤ slack, written as [A,I][x;slack] ≥ -b.
      A.rightCols(m) = Eigen::MatrixXd::Identity(m, m);
      prog->AddLinearConstraint(A, -l1c->b(), Eigen::VectorXd::Constant(m, kInf),
                                {binding.variables(), slack});
      to_remove.insert(binding);
    } else if (const auto* linfc = dynamic_cast<const drake::solvers::LInfNormCost*>(cost)) {
      const int m = linfc->A().rows();
      const int n = linfc->A().cols();
      auto slack = prog->NewContinuousVariables<1>("slack");
      prog->AddLinearCost(drake::Vector1d::Ones(), slack);
      // ∀i, aᵢᵀx + bᵢ ≤ slack, written as [A,-1][x;slack] ≤ -b.
      Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, n + 1);
      A << linfc->A(), Eigen::VectorXd::Constant(m, -1);
      prog->AddLinearConstraint(A, Eigen::VectorXd::Constant(m, -kInf), -linfc->b(),
                                {binding.variables(), slack});
      // ∀i, -(aᵢᵀx + bᵢ) ≤ slack, written as [A,1][x;slack] ≥ -b.
      A.col(A.cols() - 1) = Eigen::VectorXd::Ones(m);
      prog->AddLinearConstraint(A, -linfc->b(), Eigen::VectorXd::Constant(m, kInf),
                                {binding.variables(), slack});
      to_remove.insert(binding);
    } else if (const auto* pqc =
        dynamic_cast<const drake::solvers::PerspectiveQuadraticCost*>(cost)) {
      const int m = pqc->A().rows();
      const int n = pqc->A().cols();
      auto slack = prog->NewContinuousVariables<1>("slack");
      prog->AddLinearCost(drake::Vector1d::Ones(), slack);
      // Written as rotated Lorentz cone with z = [slack; Ax+b].
      Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m + 1, n + 1);
      A(0, 0) = 1;
      A.bottomRightCorner(m, n) = pqc->A();
      Eigen::VectorXd b(m + 1);
      b << 0, pqc->b();
      prog->AddRotatedLorentzConeConstraint(A, b, {slack, binding.variables()});
      to_remove.insert(binding);
    }
  }
  for (const auto& b : to_remove) {
    prog->RemoveCost(b);
  }

  if (!mosek.AreProgramAttributesSatisfied(*prog)) {
    throw std::runtime_error(fmt::format(
        "SolveConvexRestriction failed to generate a convex problem: {}",
        mosek.ExplainUnsatisfiedProgramAttributes(*prog)));
  }
}



