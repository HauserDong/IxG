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
#include <utility>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <stack>

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
      std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> old_edge_to_new_id_;
      std::unordered_map<int, std::vector<int>> old_id_to_new_id_;
      std::unordered_map<int, StateVarsType> new_id_to_state_;
      std::unordered_map<int, std::vector<int>> lbg_adj_list_;
      std::unordered_map<int, std::vector<double>> lbg_adj_cost_list_;

      std::string filename_;

//      std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> entry_id_;
//      std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> exit_id_;

// Serialization function
          void serialize(std::string& fname)
          {
            std::ofstream outFile(fname);
            if (!outFile.is_open())
            {
              std::cerr << "Error opening file for writing: " << fname << std::endl;
              return;
            }

            // Serialize old_edge_to_new_id_
            outFile << old_edge_to_new_id_.size() << "\n";
            for (const auto &entry : old_edge_to_new_id_)
            {
              outFile << entry.first.first << " " << entry.first.second << " " << entry.second.size();
              for (int value : entry.second)
              {
                outFile << " " << value;
              }
              outFile << "\n";
            }

            // Serialize old_id_to_new_id_
            outFile << old_id_to_new_id_.size() << "\n";
            for (const auto &entry : old_id_to_new_id_)
            {
              outFile << entry.first << " " << entry.second.size();
              for (int value : entry.second)
              {
                outFile << " " << value;
              }
              outFile << "\n";
            }

            // Serialize new_id_to_state_
            outFile << new_id_to_state_.size() << "\n";
            for (const auto &entry : new_id_to_state_)
            {
              outFile << entry.first << " " << entry.second.size();
              for (double value : entry.second)
              {
                outFile << " " << value;
              }
              outFile << "\n";
            }

            // Serialize lbg_adj_list_
            outFile << lbg_adj_list_.size() << "\n";
            for (const auto &entry : lbg_adj_list_)
            {
              outFile << entry.first << " " << entry.second.size();
              for (int value : entry.second)
              {
                outFile << " " << value;
              }
              outFile << "\n";
            }

            // Serialize lbg_adj_cost_list_
            outFile << lbg_adj_cost_list_.size() << "\n";
            for (const auto &entry : lbg_adj_cost_list_)
            {
              outFile << entry.first << " " << entry.second.size();
              for (double value : entry.second)
              {
                outFile << " " << value;
              }
              outFile << "\n";
            }

            outFile.close();
          }

// Deserialization function
          void deserialize(std::string& fname)
          {
            std::ifstream inFile(fname);
            if (!inFile.is_open())
            {
              std::cerr << "Error opening file for reading: " << fname << std::endl;
              return;
            }

            // Deserialize old_edge_to_new_id_
            size_t oldEdgeToNewIdSize;
            inFile >> oldEdgeToNewIdSize;
            for (size_t i = 0; i < oldEdgeToNewIdSize; ++i)
            {
              std::pair<int, int> key;
              size_t vectorSize;
              inFile >> key.first >> key.second >> vectorSize;

              std::vector<int> value(vectorSize);
              for (size_t j = 0; j < vectorSize; ++j)
              {
                inFile >> value[j];
              }

              old_edge_to_new_id_.emplace(std::move(key), std::move(value));
            }

            // Deserialize old_id_to_new_id_
            size_t oldIdToNewIdSize;
            inFile >> oldIdToNewIdSize;
            for (size_t i = 0; i < oldIdToNewIdSize; ++i)
            {
              int key;
              size_t vectorSize;
              inFile >> key >> vectorSize;

              std::vector<int> value(vectorSize);
              for (size_t j = 0; j < vectorSize; ++j)
              {
                inFile >> value[j];
              }

              old_id_to_new_id_.emplace(std::move(key), std::move(value));
            }


            // Deserialize new_id_to_state_
            size_t newIdToStateSize;
            inFile >> newIdToStateSize;
            for (size_t i = 0; i < newIdToStateSize; ++i)
            {
              int key;
              size_t vectorSize;
              inFile >> key >> vectorSize;

              std::vector<double> value(vectorSize);
              for (size_t j = 0; j < vectorSize; ++j)
              {
                inFile >> value[j];
              }

              new_id_to_state_.emplace(std::move(key), std::move(value));
            }

            // Deserialize lbg_adj_list_
            size_t lbgAdjListSize;
            inFile >> lbgAdjListSize;
            for (size_t i = 0; i < lbgAdjListSize; ++i)
            {
              int key;
              size_t vectorSize;
              inFile >> key >> vectorSize;

              std::vector<int> value(vectorSize);
              for (size_t j = 0; j < vectorSize; ++j)
              {
                inFile >> value[j];
              }

              lbg_adj_list_.emplace(std::move(key), std::move(value));
            }

            // Deserialize lbg_adj_cost_list_
            size_t lbgAdjCostListSize;
            inFile >> lbgAdjCostListSize;
            for (size_t i = 0; i < lbgAdjCostListSize; ++i)
            {
              int key;
              size_t vectorSize;
              inFile >> key >> vectorSize;

              std::vector<double> value(vectorSize);
              for (size_t j = 0; j < vectorSize; ++j)
              {
                inFile >> value[j];
              }

              lbg_adj_cost_list_.emplace(std::move(key), std::move(value));
            }

            inFile.close();
          }

      };

    LBGraph(std::string& env_name,
            const std::vector<HPolyhedron>& regions,
            const std::vector<std::pair<int, int>>& edges_between_regions,
            int order, int continuity,
            double path_length_weight, double time_weight,
            Eigen::VectorXd& vel_lb, Eigen::VectorXd& vel_ub,
            double h_min, double h_max, double hdot_min = 1e-6,
            bool smooth = false,
            std::string lbg_dir = "",
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
//      verbose = true;

      data_.filename_ = env_name +
                        "_o" + std::to_string(order) +
                        "_c" + std::to_string(continuity) +
                        "_pw" + std::to_string((int)path_length_weight) +
                        "_tw" + std::to_string((int)time_weight) +
                        "_hmin" + std::to_string((int)h_min) +
                        "_hmax" + std::to_string((int)h_max) +
                        "_hdmin" + std::to_string((int)hdot_min) +
                        "_sm" + std::to_string(smooth? 1: 0);

      std::string filename = lbg_dir + data_.filename_;
      if (Load(filename)) {
        return;
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
//        temp++;
//        if (temp > 3) {break;}
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
        auto soln = gcs_->Solve(edge);
        double cost = gcs_->CalculateCost(soln);

        int in_id = new_id++;
        int out_id = new_id++;
        std::pair<int, int> nz_lb_edge = {in_id, out_id};
        /// Save the edges with costs
        lb_edge_to_costs_[nz_lb_edge] = cost;
        /// Save the new id with states
        auto p0 = soln.first.value(soln.first.start_time());
        auto pF = soln.first.value(soln.first.start_time());
        std::vector<double> _p0(p0.data(), p0.data() + p0.size());
        std::vector<double> _pF(p0.data(), p0.data() + p0.size());
        data_.new_id_to_state_[in_id] = _p0;
        data_.new_id_to_state_[out_id] = _pF;

        /// this works
        data_.old_edge_to_new_id_[{edge[0], edge[1]}].push_back(in_id);
        data_.old_edge_to_new_id_[{edge[1], edge[2]}].push_back(out_id);

        data_.old_id_to_new_id_[edge[1]].push_back(in_id);
        data_.old_id_to_new_id_[edge[1]].push_back(out_id);
//        /// keeping track of exits and entries separately
//        data_.entry_id_[{edge[0], edge[1]}].push_back(in_id);
//        data_.exit_id_[{edge[1], edge[2]}].push_back(out_id);
      }

      std::cout << "LB edge to cost: " << std::endl;
      /// Add the zero cost edges
      for (auto& edge : edges_between_regions) {
        if (data_.old_edge_to_new_id_.find(edge) == data_.old_edge_to_new_id_.end()) { continue; }

        if (verbose) {
          std::cout << "Edge: " << edge.first << " " << edge.second << " pegs: ";
          for (auto peg: data_.old_edge_to_new_id_[edge]) {
            std::cout << " " << peg;
          }
          std::cout << std::endl;
        }

        auto in_new_id = data_.old_edge_to_new_id_[edge];
        auto out_new_id = data_.old_edge_to_new_id_[edge];

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

//            VecDf s1 = Eigen::Map<VecDf, Eigen::Unaligned>(data_.new_id_to_state_[in].data(), data_.new_id_to_state_[in].size());
//            VecDf s2 = Eigen::Map<VecDf, Eigen::Unaligned>(data_.new_id_to_state_[out].data(), data_.new_id_to_state_[out].size());
//            lb_edge_to_costs_[zero_lb_edge] = (s1-s2).norm();
          }
        }
      }

      /// Build the adjacency list of LB graph
      for (auto& edge : lb_edge_to_costs_) {
        data_.lbg_adj_list_[edge.first.first].push_back(edge.first.second);
        data_.lbg_adj_cost_list_[edge.first.first].push_back(edge.second);
      }

      /// Write to file
      Save(filename);

      if (verbose) {
        /// print old id to new id
        std::cout << "GCS graph ID to LBG ID: " << std::endl;
        for (auto &node: data_.old_edge_to_new_id_) {
          std::cout << node.first.first << ", " << node.first.second << ": ";
          for (auto &adj: node.second) {
            std::cout << adj << " ";
          }
          std::cout << std::endl;
        }


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

    void Save(std::string& filename) {
      data_.serialize(filename);
    }

    bool Load(std::string& filename) {
      std::cout << "Looking for file " << filename << std::endl;
      std::ifstream inFile(filename);
      if (inFile.is_open()) {
        data_.deserialize(filename);
        std::cout << "Loaded Lower Bound Graph from file! So not constructing one." << std::endl;
        return true;
      } else {
        std::cout << "No existing file found. Constructing Lower Bound Graph!" << std::endl;
        return false;
      }
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

    std::vector<HPolyhedron> hpoly_regions_;
    std::vector<std::pair<int, int>> edges_bw_regions_;
    std::shared_ptr<GCSOpt> gcs_;

    std::vector<std::vector<int>> lbg_opt_edges_;
    std::unordered_map<std::pair<int, int>, double, hash_pair> lb_edge_to_costs_;

    Data data_;
  };

  class LBGSearch {
  public:

    LBGSearch() {}
    LBGSearch(std::string& lbg_file) {
      data_.deserialize(lbg_file);
    }


    std::map<int, double> Dijkstra(StateVarsType& start_state,
                                             int gcs_start_id) {

      int start = 2*data_.new_id_to_state_.size();
      VecDf e_start_state = Eigen::Map<VecDf, Eigen::Unaligned>(start_state.data(), start_state.size());

      for (const auto& otn : data_.old_edge_to_new_id_) {
        if (otn.first.first == gcs_start_id) {
          for (const auto& c : otn.second) {
            data_.lbg_adj_list_[start].push_back(c);
            VecDf e_new_state = Eigen::Map<VecDf, Eigen::Unaligned>(data_.new_id_to_state_[c].data(), data_.new_id_to_state_[c].size());
            double cost = (e_new_state - e_start_state).norm();
            data_.lbg_adj_cost_list_[start].push_back(cost);
          }
        }
      }

      std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, double>>> pq;
      std::map<int, double> new_dist;
      for (const auto& node : data_.lbg_adj_list_) {
        new_dist[node.first] = DINF;
      }

      new_dist[start] = 0.0;
      pq.push({0.0, start});

      while (!pq.empty()) {
        int u = pq.top().second;
        int dist_u = pq.top().first;
        pq.pop();

        if (dist_u > new_dist[u]) {
          continue; // Skip outdated entries in priority queue
        }

        for (int i=0; i<data_.lbg_adj_list_.at(u).size(); ++i) {
          int v = data_.lbg_adj_list_.at(u)[i];
          double weight = data_.lbg_adj_cost_list_.at(u)[i];

          if (new_dist[u] + weight < new_dist[v]) {
            new_dist[v] = new_dist[u] + weight;
            pq.emplace(new_dist[v], v);
          }
        }
      }

//      std::map<int, double> old_dist;
//      for (const auto& otn : data_.old_edge_to_new_id_) {
//        double d = DINF/2;
//        for (const auto n : otn.second) {
//          if (new_dist.find(n) != new_dist.end()) {
//            if (new_dist[n] < d) {
//              d = new_dist[n];
//            }
//          }
//        }
//        old_dist[otn.first.first] = d; /// @FIXME: Is this correct?
//      }

      std::map<int, double> old_dist;
      for (const auto& id : data_.old_edge_to_new_id_) {
        old_dist[id.first.first] = DINF;
        old_dist[id.first.second] = DINF;
      }
      for (const auto& otn : data_.old_edge_to_new_id_) {
        for (const auto n : otn.second) {
          if (new_dist.find(n) != new_dist.end()) {
            if (new_dist[n] < old_dist[otn.first.first]) {
              old_dist[otn.first.first] = new_dist[n];
            }
            if (new_dist[n] < old_dist[otn.first.second]) {
              old_dist[otn.first.second] = new_dist[n];
            }
          }
        }
      }

//      std::map<int, double> old_dist;
//      for (const auto& id : data_.old_id_to_new_id_) {
//        old_dist[id.first] = DINF;
//      }
//      for (const auto& otn : data_.old_id_to_new_id_) {
//        for (const auto n : otn.second) {
//          if (new_dist.find(n) != new_dist.end()) {
//            if (new_dist[n] < old_dist[otn.first]) {
//              old_dist[otn.first] = new_dist[n];
//            }
//          }
//        }
//      }

      return old_dist;
    }

    int FindNumConnectedComponents() {
      const auto& graph = data_.lbg_adj_list_;
      std::unordered_set<int> visited;
      int components = 0;

      for (const auto& pair : graph) {
        int node = pair.first;

        if (visited.find(node) == visited.end()) {
          // Start a new traversal from an unvisited node
          components++;
          std::stack<int> stack;
          stack.push(node);

          while (!stack.empty()) {
            int current = stack.top();
            stack.pop();

            if (visited.find(current) == visited.end()) {
              visited.insert(current);

              // Add unvisited neighbors to the stack
              for (int neighbor : graph.at(current)) {
                if (visited.find(neighbor) == visited.end()) {
                  stack.push(neighbor);
                }
              }
            }
          }
        }
      }

      return components;
    }


    void DFSUtil(std::unordered_map<int, std::vector<int>>& graph, int vertex, std::vector<bool>& visited) {
      visited[vertex] = true;
      for (int neighbor : graph[vertex]) {
        if (!visited[neighbor]) {
          DFSUtil(graph, neighbor, visited);
        }
      }
    }

    int countConnectedComponents() {
      auto& graph = data_.lbg_adj_list_;
      int numComponents = 0;
      int numVertices = graph.size();
      std::vector<bool> visited(numVertices, false);

      for (auto& vertex : graph) {
        if (!visited[vertex.first]) {
          DFSUtil(graph, vertex.first, visited);
          numComponents++;
        }
      }

      return numComponents;
    }

    LBGraph::Data data_;

  };

}


#endif //IXG_LBGRAPH_HPP
