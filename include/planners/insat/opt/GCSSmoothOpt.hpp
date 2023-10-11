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
 * \file GCSFullOpt.hpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 10/5/23
*/

#pragma once
#ifndef GCSFULLOPT_HPP
#define GCSFULLOPT_HPP

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <chrono>

#include <drake/common/trajectories/bezier_curve.h>
#include <drake/common/trajectories/bspline_trajectory.h>
#include <drake/common/trajectories/composite_trajectory.h>
#include <drake/common/symbolic/decompose.h>
#include <drake/solvers/constraint.h>
#include <drake/common/pointer_cast.h>
#include <drake/geometry/optimization/hpolyhedron.h>
#include <drake/geometry/optimization/point.h>
#include <drake/geometry/optimization/convex_set.h>
#include <drake/geometry/optimization/graph_of_convex_sets.h>
//#include <drake/planning/trajectory_optimization/gcs_trajectory_optimization.h>
#include <drake/geometry/optimization/cartesian_product.h>
#include <drake/solvers/solve.h>
#include <drake/common/trajectories/trajectory.h>
#include <drake/solvers/mosek_solver.h>

#include <common/insat/InsatTypes.hpp>
#include <common/insatxgcs/utils.hpp>

namespace ps {

  typedef drake::geometry::optimization::GraphOfConvexSets GCS;
  typedef drake::geometry::optimization::GraphOfConvexSets::Vertex GCSVertex;
  typedef drake::geometry::optimization::GraphOfConvexSets::Edge GCSEdge;
  typedef drake::geometry::optimization::GraphOfConvexSets::VertexId VertexId;
  typedef drake::geometry::optimization::GraphOfConvexSets::EdgeId EdgeId;
  typedef drake::geometry::optimization::HPolyhedron HPolyhedron;

  typedef drake::solvers::Binding<drake::solvers::Cost> CostBinding;
  typedef drake::solvers::Binding<drake::solvers::Constraint> ConstraintBinding;
  typedef std::shared_ptr<drake::solvers::Cost> CostPtr;
  typedef std::shared_ptr<drake::solvers::Constraint> ConstraintPtr;

  using VectorXb = Eigen::Matrix<bool, 1, Eigen::Dynamic>;

  const double kInf = std::numeric_limits<double>::infinity();

  struct GCSSmoothOptResult {
    drake::trajectories::CompositeTrajectory<double> traj_;
    drake::solvers::MathematicalProgramResult result_;
    drake::VectorX<drake::symbolic::Variable> vars_;
    std::vector<drake::solvers::Binding<drake::solvers::Cost>> costs_;
    std::vector< drake::solvers::Binding<drake::solvers::Constraint>> constraints_;
  };

  class GCSSmoothOpt {

  public:

    /// Sets up the GCS regions and edges
    GCSSmoothOpt(const std::vector<HPolyhedron>& regions,
           const std::vector<std::pair<int, int>>& edges_between_regions,
           int order, int continuity,
           double path_length_weight, double time_weight,
           Eigen::VectorXd& vel_lb, Eigen::VectorXd& vel_ub,
           double h_min, double h_max, double hdot_min = 1e-6,
           bool verbose=false);

    GCSSmoothOpt (const GCSSmoothOpt &)=default;
    GCSSmoothOpt & 	operator= (const GCSSmoothOpt &)=default;
    GCSSmoothOpt (GCSSmoothOpt &&)=default;
    GCSSmoothOpt & 	operator= (GCSSmoothOpt &&)=default;

    void FormulateAndSetCostsAndConstraints();

    VertexId AddStart(Eigen::VectorXd& start);

    VertexId AddGoal(Eigen::VectorXd& goal);

    std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult> Solve(std::vector<VertexId>& path_vids,
                                                         std::vector<EdgeId>& path_eids);

    std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult> Solve(std::vector<VertexId>& path_vids,
                                                         std::vector<EdgeId>& path_eids,
                                                         Eigen::VectorXd& initial_guess);

    std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult> Solve(std::vector<VertexId>& path_vids);

    std::pair<drake::trajectories::CompositeTrajectory<double>,
        drake::solvers::MathematicalProgramResult> Solve(std::vector<VertexId>& path_vids,
                                                         Eigen::VectorXd& initial_guess);

//    GCSSmoothOptResult Solve(std::vector<VertexId>& new_path_vids,
//                       std::vector<EdgeId>& new_path_eids,
//                       Eigen::VectorXd& initial_guess,
//                       drake::VectorX<drake::symbolic::Variable>& old_dec_vars,
//                       std::vector< drake::solvers::Binding< drake::solvers::Cost > >& old_costs,
//                       std::vector< drake::solvers::Binding< drake::solvers::Constraint > >& old_constraints);

    const std::shared_ptr<drake::geometry::optimization::GraphOfConvexSets> GetGCS() const {
      return gcs_;
    }

    const std::unordered_map<int64_t, GCSVertex*> GetVertexIdToVertexMap() const {
      return vertex_id_to_vertex_;
    }

    void CleanUp();

    //// TEMPORARY
    std::vector<drake::geometry::optimization::GraphOfConvexSets::Vertex*> GetVertices() const {
      return vertices_;
    }

    std::vector<drake::geometry::optimization::GraphOfConvexSets::Edge*> GetEdges() const {
      return edges_;
    }

    void formulateStartPointConstraint();

  private:

    void setupVars();
    /// Preprocess regions to add the time scaling set and create vertices
    void preprocess(const drake::geometry::optimization::ConvexSets& regions,
                    const std::vector<std::pair<int, int>>& edges_between_regions);
    void addCosts(const GCSVertex* v);
    void addConstraints(const GCSVertex* v);
    void addConstraints(const GCSEdge* e);
    void setupCostsAndConstraints();
    void formulateTimeCost();
    void formulatePathLengthCost();
    void formulateContinuityConstraint();
    void formulateVelocityConstraint();
    void formulateCostsAndConstraints();
    std::vector<Eigen::MatrixX<drake::symbolic::Expression>>
          ControlPointsOf(const Eigen::MatrixX<drake::symbolic::Variable>& mat);
// dynamic_unique_cast for unique_ptr
    template <typename To, typename From, typename Deleter>
    std::unique_ptr<To, Deleter> DynamicUniqueCast(std::unique_ptr<From, Deleter>&& p);

    bool verbose_;

    /// Basics
    int order_;
    int continuity_;
    int num_positions_;
    double h_min_;
    double h_max_;
    double hdot_min_;
    Eigen::VectorXd vel_lb_;
    Eigen::VectorXd vel_ub_;
    std::shared_ptr<drake::geometry::optimization::GraphOfConvexSets> gcs_;

    /// Terminals
    GCSVertex* start_vtx_;
    GCSVertex* goal_vtx_;
    std::vector<GCSVertex*>::iterator start_vit_;
    std::vector<GCSVertex*>::iterator goal_vit_;
    std::vector<GCSEdge*>::iterator start_eit_;
    std::vector<GCSEdge*>::iterator goal_eit_;

    /// Variables
    Eigen::VectorX<drake::symbolic::Variable> u_h_;
    Eigen::VectorX<drake::symbolic::Variable> u_vars_;
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> u_r_trajectory_;
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> u_h_trajectory_;

    Eigen::VectorX<drake::symbolic::Variable> v_h_;
    Eigen::VectorX<drake::symbolic::Variable> v_vars_;
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> v_r_trajectory_;
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> v_h_trajectory_;

    /// Vertices and edges for optimization (for matrix operations over control points)
    std::vector<drake::geometry::optimization::GraphOfConvexSets::Vertex*> vertices_;
    std::vector<drake::geometry::optimization::GraphOfConvexSets::Edge*> edges_;

    /// Costs
    std::shared_ptr<drake::solvers::Cost> time_cost_;
    std::vector<std::pair<std::shared_ptr<drake::solvers::Cost>, VectorXb>>
        path_length_cost_;

    /// Constraints
    std::vector<std::pair<std::shared_ptr<drake::solvers::Constraint>, VectorXb>> continuity_constraint_;
    std::vector<std::pair<std::shared_ptr<drake::solvers::Constraint>, VectorXb>>
        velocity_constraint_;
    std::pair<std::shared_ptr<drake::solvers::Constraint>, VectorXb> start_point_constraint_;

    /// Maps
    /// Dict for vertex id to vertex
    std::unordered_map<int64_t, drake::geometry::optimization::GraphOfConvexSets::Vertex*> vertex_id_to_vertex_;
    /// Dict for vertex id to cost binding
//    std::unordered_map<VertexId, std::vector<CostBinding>> vertex_id_to_cost_binding_;
    std::unordered_map<int64_t, std::vector<CostBinding>> vertex_id_to_cost_binding_;
    /// Dict for vertex id to constraint binding
//    std::unordered_map<VertexId, std::vector<ConstraintBinding>> vertex_id_to_constraint_binding_;
    std::unordered_map<int64_t, std::vector<ConstraintBinding>> vertex_id_to_constraint_binding_;
    /// Dict for edge id to edge
    std::unordered_map<int64_t, drake::geometry::optimization::GraphOfConvexSets::Edge*> edge_id_to_edge_;
    /// Dict for edge id to cost binding
//    std::unordered_map<EdgeId, std::vector<CostBinding>> edge_id_to_cost_binding_;
    std::unordered_map<int64_t, std::vector<CostBinding>> edge_id_to_cost_binding_;
    /// Dict for edge id to constraint binding+-
//    std::unordered_map<EdgeId, std::vector<ConstraintBinding>> edge_id_to_constraint_binding_;
    std::unordered_map<int64_t, std::vector<ConstraintBinding>> edge_id_to_constraint_binding_;

    /// Flags for enabling/disabling costs and constraints
    bool enable_time_cost_;
    double time_weight_;
    bool enable_path_length_cost_;
    Eigen::MatrixXd path_length_weight_;
    bool enable_path_velocity_constraint_;


  };
}

#endif //INSATxGCS_OPT_HPP
