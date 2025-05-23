/*
 * Copyright (c) 2023, Ramkumar Natarajan
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*!
 * \file   INSATxGCSActions.cpp
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   2/19/23
 */

#include <cstdlib>
#include "INSATxGCSAction.hpp"

namespace ps
{
  INSATxGCSAction::INSATxGCSAction(const std::string& type,
                                   ParamsType params,
                                   OptVecPtrType opt,
                                   OptType& lb_opt,
                                   bool is_expensive) : InsatAction(type, params, is_expensive),
                                                        opt_(opt), lb_opt_(lb_opt)
  {
    const auto& edges_between_regions = (*opt_)[0].GetGCS()->Edges();
    for (auto& e : edges_between_regions) {
      adjacency_list_[e->u().id().get_value()-1].push_back(e->v().id().get_value()-1);
    }
    path_length_weight_ = params_["path_length_weight"];
    time_weight_ = params_["time_weight"];
  }

  bool INSATxGCSAction::CheckPreconditions(const StateVarsType& state, int thread_id)
  {
    return true;
  }

  ActionSuccessor INSATxGCSAction::GetSuccessor(const StateVarsType& state_vars, int thread_id)
  {
    if (adjacency_list_.find(static_cast<int>(state_vars[0])) == adjacency_list_.end()) {
      // Check whether the start vertex in adjacency list
      throw std::runtime_error("State " + std::to_string(static_cast<int>(state_vars[0])) + " not found in successor map!!");
    }

    int idx = std::stoi(type_); // an integer indicating the index of the action
    if (adjacency_list_[static_cast<int>(state_vars[0])].size() > idx) {
      StateVarsType succ;
      succ.push_back(adjacency_list_[static_cast<int>(state_vars[0])][idx]);  // add the vertex index of the neighbor set
      return ActionSuccessor(true, {std::make_pair(succ, 1.0)});
    }
    return ActionSuccessor(false, {make_pair(StateVarsType(), -DINF)});
  }

  ActionSuccessor INSATxGCSAction::GetSuccessorLazy(const StateVarsType& state_vars, int thread_id)
  {

  }

  ActionSuccessor INSATxGCSAction::Evaluate(const StateVarsType& parent_state_vars, const StateVarsType& child_state_vars, int thread_id)
  {
    return GetSuccessor(parent_state_vars, thread_id);
  }

  void INSATxGCSAction::UpdateStateToSuccs() {
    adjacency_list_.clear();
    const auto& edges_between_regions = (*opt_)[0].GetGCS()->Edges();
    for (auto& e : edges_between_regions) {
      adjacency_list_[static_cast<int>(e->u().id().get_value()-1)].push_back(static_cast<int>(e->v().id().get_value()-1));
    }
  }

  double getRandomNumberBetween(double min, double max, std::mt19937& gen)
  {
    std::uniform_real_distribution<double> distr(min, max);
    return distr(gen);
  }

  double INSATxGCSAction::GetCostToSuccessor(const StateVarsType &current_state, const StateVarsType &successor_state, int thread_id)
  {
    return 1.0;
  }

  double INSATxGCSAction::getCostToSuccessor(const VecDf &current_state, const VecDf &successor_state, int thread_id)
  {
    return (successor_state-current_state).norm();
  }


  /// INSAT
  void INSATxGCSAction::setOpt(OptVecPtrType& opt)
  {
    opt_ = opt;
  }

  TrajType INSATxGCSAction::optimize(const std::vector<StateVarsType> &ancestors,
                    const StateVarsType& successor,
                    int thread_id)
  {
    const auto& vivm = (*opt_)[thread_id].GetVertexIdToVertexMap();
    std::vector<VertexId> solve_vids;
    for (auto vid : ancestors) {
      auto it = vivm.find(static_cast<int>(vid[0]));
      if (it == vivm.end()) {
        throw std::runtime_error("State with VId:" + std::to_string(static_cast<int>(vid[0])) + " not found in GCS graph!!");
      }
      solve_vids.push_back(it->second->id());
    }
    auto it = vivm.find(static_cast<int>(successor[0]));
    if (it == vivm.end()) {
      throw std::runtime_error("State with VId:" + std::to_string(static_cast<int>(successor[0])) + " not found in GCS graph!!");
    }
    solve_vids.push_back(it->second->id());

    auto soln = (*opt_)[thread_id].Solve(solve_vids);
    return TrajType (soln.first, soln.second);
  }

  TrajType INSATxGCSAction::optimize(const std::vector<int> &gcs_nodes, int thread_id) {
    const auto& vivm = (*opt_)[thread_id].GetVertexIdToVertexMap();
    std::vector<VertexId> solve_vids;
    for (auto vid : gcs_nodes) {
      auto it = vivm.find(vid);
      if (it == vivm.end()) {
        throw std::runtime_error("State with VId:" + std::to_string(vid) + " not found in GCS graph!!");
      }
      solve_vids.push_back(it->second->id());
    }

    auto soln = (*opt_)[thread_id].Solve(solve_vids);
    return TrajType (soln.first, soln.second);
  }

  double INSATxGCSAction::lowerboundCost(const std::vector<int> &gcs_nodes, int thread_id) {
    return (*opt_)[thread_id].LowerboundSolve(gcs_nodes);
//    return lb_opt_.LowerboundSolve(gcs_nodes);
  }

  double INSATxGCSAction::getCost(const TrajType &traj, int thread_id) const
  {
//    return traj.result_.get_optimal_cost();


    if (traj.size() == 0) {
      auto disc_traj = sampleTrajectory(traj.traj_, 1e-1);
      return calculateCost(disc_traj)*path_length_weight_ +
             traj.traj_.end_time()*time_weight_;
    } else {
      return calculateCost(traj.disc_traj_)*path_length_weight_ +
              traj.traj_.end_time()*time_weight_;
    }
    return calculateCost(traj.disc_traj_)*path_length_weight_ +
           traj.traj_.end_time()*time_weight_;
  }

  MatDf INSATxGCSAction::sampleTrajectory(const GCSTraj::TrajInstanceType &traj, double dt) const {
    MatDf sampled_traj;
    int i=0;
    for (double t=0.0; t<=traj.end_time(); t+=dt)
    {
      sampled_traj.conservativeResize(traj.rows(), sampled_traj.cols()+1);
      sampled_traj.col(i) = traj.value(t);
      ++i;
    }
    return sampled_traj;
  }

  double INSATxGCSAction::calculateCost(const MatDf &disc_traj) const {
    double cost = 0;
    for (int i=0; i<disc_traj.cols()-1; ++i)
    {
      cost += (disc_traj.col(i+1)-disc_traj.col(i)).norm();
    }
    return cost;
  }

  std::unordered_map<int, std::vector<int>> INSATxGCSAction::getAdjacencyList() {return adjacency_list_;}

}
