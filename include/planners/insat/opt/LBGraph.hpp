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
#include <common/Types.hpp>

#include <sstream>
#include <fstream>
#include <queue>

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

    struct Data
    {
      std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> old_id_to_new_id_;
//      std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> entry_id_;
//      std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> exit_id_;
      std::unordered_map<int, std::vector<int>> lbg_adj_list_;
      std::unordered_map<int, std::vector<double>> lbg_adj_cost_list_;
      std::unordered_map<int, std::vector<StateVarsType>> lbg_adj_state_list_;
    };

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
        gcs_->FormulateAndSetCostsAndConstraints();
      }

      /// Get the adjacency list of GCS graph
//      std::unordered_map<int, std::vector<int>> adjacency_list_;
      std::map<int, std::vector<int>> adjacency_list_;
      for (auto& edge : edges_between_regions) {
        adjacency_list_[edge.first].push_back(edge.second);
      }

      /// Build the set of non-zero cost LB edge triplets (in_id, node_id, out_id)
      int temp = 0;
      for (auto& out : adjacency_list_) {
        int ctr = out.first;
        for (int v=0; v < out.second.size(); v++) {
          for (int w=v+1; w < out.second.size(); w++) {
            lbg_opt_edges_.push_back({out.second[v], ctr, out.second[w]});
            lbg_opt_edges_.push_back({out.second[w], ctr, out.second[v]});
          }
        }
        temp++;
        if (temp > 3) {
          break;
        }
      }

      if (verbose) {
        std::cout << "LB edges: " << std::endl;
        for (auto l: lbg_opt_edges_) {
          std::cout << l[0] << " " << l[1] << " " << l[2] << std::endl;
        }
      }

      int new_id = 0;
      /// Solve the LB problem for each LB edge using the triplets
      /// Save the edge if it is new or if the cost is lower
      /// Update the map from old id to new id
      for (auto& edge : lbg_opt_edges_) {
        double cost = gcs_->LowerboundSolve(edge);
        int in_id = new_id++;
        int out_id = new_id++;
        std::pair<int, int> nz_lb_edge = {in_id, out_id};
        /// Save the edges with costs
        lb_edge_to_costs_[nz_lb_edge] = cost;

        /// this works
        data_.old_id_to_new_id_[{edge[0], edge[1]}].push_back(in_id);
        data_.old_id_to_new_id_[{edge[1], edge[2]}].push_back(out_id);

//        /// keeping track of exits and entries separately
//        data_.entry_id_[{edge[0], edge[1]}].push_back(in_id);
//        data_.exit_id_[{edge[1], edge[2]}].push_back(out_id);
      }

      std::cout << "LB edge to cost: " << std::endl;
      /// Add the zero cost edges
      for (auto& edge : edges_between_regions) {
        if (data_.old_id_to_new_id_.find(edge) == data_.old_id_to_new_id_.end()) { continue; }

        if (verbose) {
          std::cout << "Edge: " << edge.first << " " << edge.second << " pegs: ";
          for (auto peg: data_.old_id_to_new_id_[edge]) {
            std::cout << " " << peg;
          }
          std::cout << std::endl;
        }

        auto in_new_id = data_.old_id_to_new_id_[edge];
        auto out_new_id = data_.old_id_to_new_id_[edge];

//        auto in_new_id = data_.entry_id_[edge];
//        auto out_new_id = data_.exit_id_[edge];

        for (auto& in : in_new_id) {
          for (auto& out : out_new_id) {
            if (in == out) { continue;}
            std::pair<int, int> zero_lb_edge = {in, out};
            if (lb_edge_to_costs_.find(zero_lb_edge) != lb_edge_to_costs_.end()) {
              throw std::runtime_error("Zero LB found!!");
            }
            lb_edge_to_costs_[zero_lb_edge] = 0.0;
          }
        }
      }

      /// Build the adjacency list of LB graph
      for (auto& edge : lb_edge_to_costs_) {
        data_.lbg_adj_list_[edge.first.first].push_back(edge.first.second);
        data_.lbg_adj_cost_list_[edge.first.first].push_back(edge.second);
      }

      if (verbose) {
        /// print adjacency list
        std::cout << "LB Graph Adjacency List: " << std::endl;
        for (auto &node: data_.lbg_adj_list_) {
          std::cout << node.first << ": ";
          for (auto &adj: node.second) {
            std::cout << adj << " ";
          }
          std::cout << std::endl;
        }

        /// print cost list
        std::cout << "LB Graph Cost List: " << std::endl;
        for (auto &node: data_.lbg_adj_cost_list_) {
          std::cout << node.first << ": ";
          for (auto &adj: node.second) {
            std::cout << adj << " ";
          }
          std::cout << std::endl;
        }
      }
    }

    std::unordered_map<int, std::vector<int>> GetLBAdjacencyList() const {
      return data_.lbg_adj_list_;
    }

    std::unordered_map<int, std::vector<double>> GetLBAdjacencyCostList() const {
      return data_.lbg_adj_cost_list_;
    }

    void PrintLBGraphStats() const {
      int degree = 0;
      int num_edges = 0;
      for (auto& kv : data_.lbg_adj_list_) {
        if (kv.second.size() > degree) {
          degree = kv.second.size();
        }
        num_edges += kv.second.size();
      }
      std::cout << "LB Graph Stats: " << std::endl;
      std::cout << "Number of nodes: " << data_.lbg_adj_list_.size() << std::endl;
      std::cout << "Number of edges: " << num_edges << std::endl;
      std::cout << "Degree: " << degree << std::endl;
    }

    std::string createFileName() {
      return "";
    }

// Serialize Data to a file in ASCII format
    void serializeToFile(const std::string& filePath) {
      std::ofstream file(filePath);
      if (file.is_open()) {
        // Serialize data_.lbg_adj_list_
        for (const auto& entry : data_.lbg_adj_list_) {
          file << "adj " << entry.first << ':';

          // Serialize the vector elements
          for (int value : entry.second) {
            file << value << ',';
          }
          file.seekp(-1, std::ios_base::end); // Remove the last comma
          file << '\n';
        }

        // Serialize data_.lbg_adj_cost_list_
        for (const auto& entry : data_.lbg_adj_cost_list_) {
          file << "cost " << entry.first << ':';

          // Serialize the vector elements
          for (double value : entry.second) {
            file << value << ',';
          }
          file.seekp(-1, std::ios_base::end); // Remove the last comma
          file << '\n';
        }

        // Serialize data_.lbg_adj_cost_list_
//        for (const auto& entry : data_.old_id_to_new_id_) {
//          file << "idmap " << entry.first << ':';
//
//          // Serialize the vector elements
//          for (int value : entry.second) {
//            file << value << ',';
//          }
//          file.seekp(-1, std::ios_base::end); // Remove the last comma
//          file << '\n';
//        }

        file.close();
      } else {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
      }
    }

// Deserialize Data from a file in ASCII format
    Data deserializeFromFile(const std::string& filePath) {
      Data data;
      std::ifstream file(filePath);
      if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
          std::istringstream iss(line);
          std::string mapType;
          iss >> mapType;

          int key;
          char colon;
          iss >> key >> colon;

          // Deserialize vector elements
          if (mapType == "adj") {
            int value;
            std::vector<int> values;
            while (iss >> value) {
              values.push_back(value);
              char comma;
              if (!(iss >> comma) || comma != ',') {
                break;
              }
            }
            data.lbg_adj_list_[key] = values;
          } else if (mapType == "cost") {
            double value;
            std::vector<double> values;
            while (iss >> value) {
              values.push_back(value);
              char comma;
              if (!(iss >> comma) || comma != ',') {
                break;
              }
            }
            data.lbg_adj_cost_list_[key] = values;
          }
//          else if (mapType == "idmap") {
//            int value;
//            std::vector<int> values;
//            while (iss >> value) {
//              values.push_back(value);
//              char comma;
//              if (!(iss >> comma) || comma != ',') {
//                break;
//              }
//            }
//            data.old_id_to_new_id_[key] = values;
//          }
        }

        file.close();
      } else {
        std::cerr << "Error: Unable to open file for reading." << std::endl;
      }

      return data;
    }



    std::vector<HPolyhedron> hpoly_regions_;
    std::vector<std::pair<int, int>> edges_bw_regions_;
    std::shared_ptr<GCSOpt> gcs_;

    std::vector<std::vector<int>> lbg_opt_edges_;
    std::unordered_map<std::pair<int, int>, double, hash_pair> lb_edge_to_costs_;

    Data data_;
  };

  class LBGSearch {
  public:
    const int INF = 1e9;

    std::map<int, double> Dijkstra(const std::unordered_map<int, std::vector<int>>& graph,
                                             const std::unordered_map<int, std::vector<double>>& costs,
                                             int start) {

      std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, double>>> pq;
      std::map<int, double> distance;
      for (const auto& node : graph) {
        distance[node.first] = INF;
      }

      distance[start] = 0.0;
      pq.push({0.0, start});

      while (!pq.empty()) {
        int u = pq.top().second;
        int dist_u = pq.top().first;
        pq.pop();

        if (dist_u > distance[u]) {
          continue; // Skip outdated entries in priority queue
        }

        for (int i=0; i<graph.at(u).size(); ++i) {
          int v = graph.at(u)[i];
          double weight = costs.at(u)[i];

          if (distance[u] + weight < distance[v]) {
            distance[v] = distance[u] + weight;
            pq.emplace(distance[v], v);
          }
        }
      }
      return distance;
    }

  };

}


#endif //IXG_LBGRAPH_HPP
