//  * Copyright (c) 2024, Ramkumar Natarajan
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
 * \file LBGraph.hpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 1/15/24
*/

#ifndef IXG_LBGRAPH_HPP
#define IXG_LBGRAPH_HPP

#include <planners/insat/opt/GCSOpt.hpp>
#include <planners/insat/opt/GCSSmoothOpt.hpp>

namespace ps {

  struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
      auto hash1 = std::hash<T1>{}(p.first);
      auto hash2 = std::hash<T2>{}(p.second);

      if (hash1 != hash2) {
        return hash1 ^ hash2;
      }

      // If hash1 == hash2, their XOR is zero.
      return hash1;
    }
  };

  class LBGraph {
  public:

    LBGraph(const std::vector<HPolyhedron>& regions,
            const std::vector<std::pair<int, int>>& edges_between_regions,
            int order, int continuity,
            double path_length_weight, double time_weight,
            Eigen::VectorXd& vel_lb, Eigen::VectorXd& vel_ub,
            double h_min, double h_max, double hdot_min = 1e-6,
            bool smooth = false,
            bool verbose=false) : hpoly_regions_(regions),
                                  edges_bw_regions_(edges_between_regions) {

      if (smooth){
        gcs_ = std::make_shared<GCSSmoothOpt>(regions, edges_between_regions, order, continuity,
                                              path_length_weight, time_weight, vel_lb, vel_ub,
                                              h_min, h_max, hdot_min, verbose);
      } else {
        gcs_ = std::make_shared<GCSOpt>(regions, edges_between_regions, order, h_min, h_max,
                                        path_length_weight, time_weight, vel_lb, vel_ub,
                                        verbose);
      }

      /// Get the adjacency list of GCS graph
      std::unordered_map<int, std::vector<int>> adjacency_list_;
      for (auto& edge : edges_between_regions) {
        adjacency_list_[edge.first].push_back(edge.second);
      }

      /// Build the set of non-zero cost LB edge triplets (in_id, node_id, out_id)
      for (auto& out : adjacency_list_) {
        int ctr = out.first;
        for (int v=0; v < out.second.size(); v++) {
          for (int w=v+1; w < out.second.size(); w++) {
            lbg_opt_edges_.push_back({out.second[v], ctr, out.second[w]});
            lbg_opt_edges_.push_back({out.second[w], ctr, out.second[v]});
          }
        }
      }

      /// Solve the LB problem for each LB edge using the triplets
      for (auto& edge : lbg_opt_edges_) {
        double cost = gcs_->LowerboundSolve(edge);
        std::pair<int, int> nz_lb_edge = {edge[0], edge[2]};
        /// Save the edges with costs
        if (lb_edge_to_costs_.find(nz_lb_edge) != lb_edge_to_costs_.end()) {
          if (cost < lb_edge_to_costs_[nz_lb_edge]) {
            lb_edge_to_costs_[nz_lb_edge] = cost;
          }
        } else {
          lb_edge_to_costs_[nz_lb_edge] = cost;
        }
      }
      /// Add the zero cost edges
      for (auto& edge : edges_between_regions) {
        std::pair<int, int> zero_lb_edge = {edge.first, edge.second};
        if (lb_edge_to_costs_.find(zero_lb_edge) != lb_edge_to_costs_.end()) {
          throw std::runtime_error("Zero LB found!!");
        }
        lb_edge_to_costs_[zero_lb_edge] = 0.0;
      }

      /// Build the adjacency list of LB graph
      for (auto& edge : lb_edge_to_costs_) {
        lbg_adj_list_[edge.first.first].push_back(edge.first.second);
        lbg_adj_cost_list_[edge.first.first].push_back(edge.second);
      }
    }

    std::unordered_map<int, std::vector<int>> GetLBAdjacencyList() const {
      return lbg_adj_list_;
    }

    std::unordered_map<int, std::vector<double>> GetLBAdjacencyCostList() const {
      return lbg_adj_cost_list_;
    }

    std::vector<HPolyhedron> hpoly_regions_;
    std::vector<std::pair<int, int>> edges_bw_regions_;
    std::shared_ptr<GCSOpt> gcs_;

    std::vector<std::vector<int>> lbg_opt_edges_;
    std::unordered_map<std::pair<int, int>, double, hash_pair> lb_edge_to_costs_;

    std::unordered_map<int, std::vector<int>> lbg_adj_list_;
    std::unordered_map<int, std::vector<double>> lbg_adj_cost_list_;
  };

}

//namespace std {
//  // A hash function used to hash a pair of any kind
//  template <>
//  struct hash<pair<int, int>> {
//    size_t operator()(const pair<int, int>& p) const
//    {
//      std::size_t seed = 0;
//      boost::hash_combine(seed, p.first);
//      boost::hash_combine(seed, p.second);
//      return seed;
//
////      auto hash1 = hash<int>{}(p.first);
////      auto hash2 = hash<int>{}(p.second);
////
////      if (hash1 != hash2) {
////        return hash1 ^ hash2;
////      }
////
////      // If hash1 == hash2, their XOR is zero.
////      return hash1;
//    }
//  };
//}



#endif //IXG_LBGRAPH_HPP
