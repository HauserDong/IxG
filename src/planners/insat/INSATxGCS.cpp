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
 * \file   INSATxGCS.cpp
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   9/13/23
 */

#include <planners/insat/INSATxGCS.hpp>
#include <boost/functional/hash.hpp>
#include <common/insatxgcs/gcsbfs.hpp>

namespace ps
{


  INSATxGCS::INSATxGCS(ParamsType planner_params) :
          Planner(planner_params)
  {
    if (planner_params.find("adaptive_opt") == planner_params.end())
    {
      planner_params["adaptive_opt"] = false;
    }
  }

  void INSATxGCS::SetStartState(const StateVarsType &state_vars) {
    start_state_ptr_ = constructInsatState(state_vars);
  }

  bool INSATxGCS::Plan() {
    initialize();
    startTimer();
    while (!insat_state_open_list_.empty() && !checkTimeout())
    {
      auto state_ptr = insat_state_open_list_.min();
      insat_state_open_list_.pop();

      // Return solution if goal state is expanded
      if (isGoalState(state_ptr))
      {
        auto t_end = std::chrono::steady_clock::now();
        double t_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end-t_start_).count();
        goal_state_ptr_ = state_ptr;

        // Reconstruct and return path
        constructPlan(state_ptr);
        planner_stats_.total_time_ = 1e-9*t_elapsed;
        exit();
        return true;
      }

      expandState(state_ptr);

    }

    auto t_end = std::chrono::steady_clock::now();
    double t_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end-t_start_).count();
    planner_stats_.total_time_ = 1e-9*t_elapsed;
    return false;
  }

  TrajType INSATxGCS::getSolutionTraj() {
    return soln_traj_;
  }

  void INSATxGCS::initialize() {

    plan_.clear();

    // Initialize planner stats
    planner_stats_ = PlannerStats();

    // Initialize start state
    start_state_ptr_->SetGValue(0);
    start_state_ptr_->SetHValue(computeHeuristic(start_state_ptr_));

    // Reset goal state
    goal_state_ptr_ = NULL;

    // Reset state
    planner_stats_ = PlannerStats();

    // Reset h_min
    h_val_min_ = DINF;

    planner_stats_.num_jobs_per_thread_.resize(1, 0);
    // Initialize open list
    start_state_ptr_->SetFValue(start_state_ptr_->GetGValue() + heuristic_w_*start_state_ptr_->GetHValue());
    insat_state_open_list_.push(start_state_ptr_);

    constructInsatActions();

#if OPTIMAL
    calculateBounds();
#endif
  }

  void INSATxGCS::calculateBounds() {
    /// Calculating the lower and upper bounds for pruning
    // Run BFS on underlying GCS graph from start and from goal
    auto gcs_adjacency = insat_actions_ptrs_[0]->getAdjacencyList();
    ixg::GCSBFS gcsbfs;
    gcsbfs.SetAdjecency(gcs_adjacency);

    // Run BFS from start
    int start_state_id = static_cast<int>(start_state_ptr_->GetStateVars()[0]);
    paths_from_start_ = gcsbfs.BFSWithPaths(start_state_id);
    // Run BFS from goal
    int goal_state_id;
    for (auto& adj : gcs_adjacency) {
      StateVarsType goal_state(1, adj.first);
      if (goal_checker_(goal_state)) {
        goal_state_id = adj.first;
        paths_from_goal_ = gcsbfs.BFSWithPaths(static_cast<int>(goal_state[0]));
        break;
      }
    }

    // calculate bounds
    auto init_soln_path = paths_from_start_[goal_state_id];
    init_soln_path.push_back(goal_state_id);
    ub_path_ = init_soln_path;
    auto init_soln_traj = insat_actions_ptrs_[0]->optimize(init_soln_path);
    if (!init_soln_traj.isValid()) {
      throw std::runtime_error("Couldn't find initial solution trajectory.");
    }
    double global_ub = insat_actions_ptrs_[0]->getCost(init_soln_traj);
#if VERBOSE
    std::cout << "init_soln_path" << std::endl;
    printPath(init_soln_path);
#endif

    int count_infeasible = 0;
    for (auto& adj : gcs_adjacency) {
      ub_cost_[adj.first] = global_ub;

      // Lower bound
      // Option 1
      StateVarsType state(1, adj.first);
      lb_cost_[adj.first] = unary_heuristic_generator_(state);

      // Option 2
//      auto pfg = paths_from_goal_[adj.first];
//      if (pfg.empty()) {
//        lb_cost_[adj.first] = 0;
//        continue;
//      }
//      pfg.pop_back();
//      if (pfg.empty()) {
//        lb_cost_[adj.first] = 0;
//        continue;
//      }
//      std::reverse(pfg.begin(), pfg.end());
//      lb_cost_[adj.first] = insat_actions_ptrs_[0]->lowerboundCost(pfg);

      // Option 3
//      lb_cost_[adj.first] = pfg.size();

//      if (ub_cost_[adj.first] < lb_cost_[adj.first]) {
//        std::cout << "Upper bound " << ub_cost_[adj.first] << " is less than lower bound " << lb_cost_[adj.first] << std::endl;
//        std::abort();
//      }
    }

//    std::cout << "% infeasible UBs: "
//              << (((double)count_infeasible)/((double)gcs_adjacency.size()))*100 << std::endl;
    std::cout << "LB set successfully! lb from start: " << lb_cost_[start_state_id] << std::endl;
    std::cout << "UB set successfully! global_ub: " << global_ub << std::endl;

  }

  std::vector<InsatStatePtrType>
  INSATxGCS::getStateAncestors(const InsatStatePtrType state_ptr, bool reverse) const {
    // Get ancestors (including the current one)
    std::vector<InsatStatePtrType> ancestors;
    ancestors.push_back(state_ptr);
    auto bp = state_ptr->GetIncomingInsatEdgePtr();
    while (bp)
    {
      ancestors.push_back(bp->lowD_parent_state_ptr_);
      bp = bp->lowD_parent_state_ptr_->GetIncomingInsatEdgePtr();
    }
    if (reverse)
    {
      std::reverse(ancestors.begin(), ancestors.end());
    }
    return ancestors;
  }

  void INSATxGCS::expandState(InsatStatePtrType state_ptr) {

    if (VERBOSE) state_ptr->Print("Expanding");
    planner_stats_.num_state_expansions_++;

    state_ptr->SetVisited();

    auto ancestors = getStateAncestors(state_ptr, true);  // get ancestors

    for (auto& action_ptr: insat_actions_ptrs_)
    {
      if (action_ptr->CheckPreconditions(state_ptr->GetStateVars()))
      {
        // Evaluate the edge
        auto action_successor = action_ptr->GetSuccessor(state_ptr->GetStateVars());
#if OPTIMAL
        if (action_successor.success_) {
          /// Do not allow cycles
          for (const auto& anc : ancestors) {
            if (anc->GetStateVars()[0] == action_successor.successor_state_vars_costs_.back().first[0]) {
              action_successor.success_ = false;
              break;
            }
          }
        }
#endif
        updateState(state_ptr, ancestors, action_ptr, action_successor);
      }
    }
  }

  void INSATxGCS::updateState(InsatStatePtrType &state_ptr, std::vector<InsatStatePtrType> &ancestors,
                                 InsatActionPtrType &action_ptr, ActionSuccessor &action_successor) {
    // using current state, its ancestors, and action to update the state

    if (action_successor.success_)
    {
#if OPTIMAL
      auto successor_state_ptr = constructInsatPath(ancestors, action_successor.successor_state_vars_costs_.back().first);
#else
      auto successor_state_ptr = constructInsatState(action_successor.successor_state_vars_costs_.back().first);  // get the successor vertex index
#endif

      if (!successor_state_ptr->IsVisited())
      {
        planner_stats_.num_evaluated_edges_++;

#if OPTIMAL
        /// counting number of incoming rewirings to the same state
        ///////////////////////////////////////////////////////////
        int state_key = static_cast<int>(successor_state_ptr->GetStateVars()[0]);
        auto it = planner_stats_.num_incoming_edges_map_.find(state_key);
        if (it == planner_stats_.num_incoming_edges_map_.end())
        {
          planner_stats_.num_incoming_edges_map_[state_key] = 0;
          planner_stats_.num_incoming_edges_map_[state_key]++;
        }
        else
        {
          planner_stats_.num_incoming_edges_map_[state_key]++;
        }
        ///////////////////////////////////////////////////////////
#endif

        TrajType traj;
        double cost = 0;
        double inc_cost = 0;
        InsatStatePtrType best_anc;

        std::vector<StateVarsType> anc_states;
        for (auto& anc: ancestors)
        {
          anc_states.emplace_back(anc->GetStateVars());
        }
        traj = action_ptr->optimize(anc_states, successor_state_ptr->GetStateVars()); // optimize in the historical set and the upcoming one

        if (!traj.isValid())
        {
          return;
        }

        cost = action_ptr->getCost(traj);
        double new_g_val = cost;
        inc_cost = state_ptr->GetIncomingInsatEdgePtr()?
                action_ptr->getCost(traj) - action_ptr->getCost(state_ptr->GetIncomingInsatEdgePtr()->GetTraj()):
                action_ptr->getCost(traj);

#if OPTIMAL
//        double lb = new_g_val + lb_cost_[static_cast<int>(successor_state_ptr->GetStateVars()[0])];
        double lb = state_ptr->GetGValue() + lb_cost_[static_cast<int>(successor_state_ptr->GetStateVars()[0])];
//        if (0.2*lb > ub_cost_[static_cast<int>(successor_state_ptr->GetStateVars()[0])]) {
        if (lb > ub_cost_[static_cast<int>(successor_state_ptr->GetStateVars()[0])]) {
          planner_stats_.num_pruned_edges_++;
#if VERBOSE
          if (std::find(ub_path_.begin(), ub_path_.end(),successor_state_ptr->GetStateVars()[0])!=ub_path_.end()) {
            printPath(ancestors);
            std::cout << successor_state_ptr->GetStateVars()[0] << std::endl;
          }
            std::cout << "pruning cuz lb is " << lb << " g: " << state_ptr->GetGValue() << " h: " << lb_cost_[static_cast<int>(successor_state_ptr->GetStateVars()[0])] << std::endl;
#endif
          return;
        }
#endif

        if (successor_state_ptr->GetGValue() > new_g_val)
        {
          best_anc = ancestors[0];

          double h_val = successor_state_ptr->GetHValue();
          if (h_val == -1)
          {
            h_val = computeHeuristic(successor_state_ptr);
//            h_val = lb_cost_[static_cast<int>(successor_state_ptr->GetStateVars()[0])];
            successor_state_ptr->SetHValue(h_val);
          }

          if (h_val != DINF)
          {
            h_val_min_ = h_val < h_val_min_ ? h_val : h_val_min_;
            successor_state_ptr->SetGValue(new_g_val); //
            successor_state_ptr->SetFValue(new_g_val + heuristic_w_*h_val); //

            auto edge_ptr = new Edge(state_ptr, action_ptr, successor_state_ptr);
            edge_ptr->SetCost(inc_cost);
            successor_state_ptr->SetIncomingEdgePtr(edge_ptr);

            auto insat_edge_ptr = new InsatEdge(state_ptr, action_ptr, best_anc, successor_state_ptr);
            insat_edge_ptr->SetTraj(traj);
            insat_edge_ptr->SetTrajCost(cost);
            insat_edge_ptr->SetCost(cost);
            if (isGoalState(successor_state_ptr))
            {
              insat_edge_ptr->SetTrajCost(0);
              insat_edge_ptr->SetCost(0);
              successor_state_ptr->SetFValue(0.0);
            }
            edge_map_.insert(std::make_pair(getEdgeKey(insat_edge_ptr), insat_edge_ptr));
            successor_state_ptr->SetIncomingInsatEdgePtr(insat_edge_ptr); //

            if (insat_state_open_list_.contains(successor_state_ptr))
            {
              insat_state_open_list_.decrease(successor_state_ptr);
            }
            else
            {
              insat_state_open_list_.push(successor_state_ptr);
            }
          }
        }
      }
    }
  }

  void INSATxGCS::constructInsatActions() {
    for (auto& action_ptr : actions_ptrs_)
    {
      insat_actions_ptrs_.emplace_back(std::dynamic_pointer_cast<InsatAction>(action_ptr));
    }
  }

  InsatStatePtrType INSATxGCS::constructInsatState(const StateVarsType &state) {
    size_t key = state_key_generator_(state); // according to state[0], generate a hash number. same state[0] will have same hash number
    auto it = insat_state_map_.find(key);
    InsatStatePtrType insat_state_ptr;

    // Check if state exists in the search state map
    if (it == insat_state_map_.end())
    {
      insat_state_ptr = new InsatState(state);
      insat_state_map_.insert(std::pair<size_t, InsatStatePtrType>(key, insat_state_ptr));  // keep the hash number and the state in the map
    }
    else
    {
      insat_state_ptr = it->second;
    }

    return insat_state_ptr;
  }

  InsatStatePtrType
  INSATxGCS::constructInsatPath(std::vector<InsatStatePtrType> &ancestors, const ps::StateVarsType &state) {
    size_t key = 0;
    boost::hash_combine(key, state[0]);
    for (auto& anc: ancestors)
    {
      boost::hash_combine(key, anc->GetStateVars()[0]);
    }

    auto it = insat_state_map_.find(key);
    InsatStatePtrType insat_state_ptr;

    // Check if state exists in the search state map
    if (it == insat_state_map_.end())
    {
      insat_state_ptr = new InsatState(state);
      insat_state_map_.insert(std::pair<size_t, InsatStatePtrType>(key, insat_state_ptr));
    }
    else
    {
      insat_state_ptr = it->second;
    }

    return insat_state_ptr;
  }


  void INSATxGCS::cleanUp() {
    for (auto& state_it : insat_state_map_)
    {
      if (state_it.second)
      {
        delete state_it.second;
        state_it.second = NULL;
      }
    }
    insat_state_map_.clear();

    for (auto& edge_it : edge_map_)
    {
      if (edge_it.second)
      {
        delete edge_it.second;
        edge_it.second = NULL;
      }
    }
    edge_map_.clear();

    State::ResetStateIDCounter();
    Edge::ResetStateIDCounter();
  }

  void INSATxGCS::resetStates() {
    for (auto it = insat_state_map_.begin(); it != insat_state_map_.end(); ++it)
    {
      it->second->ResetGValue();
      it->second->ResetFValue();
      // it->second->ResetVValue();
      it->second->ResetIncomingInsatEdgePtr();
      it->second->UnsetVisited();
      it->second->UnsetBeingExpanded();
      it->second->num_successors_ = 0;
      it->second->num_expanded_successors_ = 0;
    }
  }

  void INSATxGCS::constructPlan(InsatStatePtrType &insat_state_ptr) {
    StatePtrType state_ptr = insat_state_ptr;
    Planner::constructPlan(state_ptr);

    if (insat_state_ptr->GetIncomingInsatEdgePtr())
    {
//                planner_stats_.path_cost_ = insat_state_ptr->GetIncomingInsatEdgePtr()->GetTrajCost();
      planner_stats_.path_cost_ =
              insat_actions_ptrs_[0]->getCost(insat_state_ptr->GetIncomingInsatEdgePtr()->GetTraj());
      soln_traj_ = insat_state_ptr->GetIncomingInsatEdgePtr()->GetTraj();
    }
  }

  void INSATxGCS::exit() {
    // Clear open list
    while (!insat_state_open_list_.empty())
    {
      insat_state_open_list_.pop();
    }

    cleanUp();
  }

  void INSATxGCS::printPath(std::vector<int> &path) {
    for (auto &state_id : path) {
      std::cout << state_id << " ";
    }
    std::cout << std::endl;
  }

  void INSATxGCS::printPath(std::vector<InsatStatePtrType> &path) {
    for (auto& state_ptr : path) {
      std::cout << state_ptr->GetStateVars()[0] << " ";
    }
    std::cout << std::endl;
  }
}
