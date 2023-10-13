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
 * \file gcsbfs.cpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 10/12/23
*/

#include <common/insatxgcs/gcsbfs.hpp>

namespace ixg {

  GCSBFS::GCSBFS(std::vector<std::pair<int, int>> &edges_bw_regions) {
    for (auto& edge : edges_bw_regions) {
      AddEdge(edge.first, edge.second);
    }
  }

  void GCSBFS::AddEdge(int u, int v) {
    adjacency_[u].push_back(v);
//    adjacency_[v].push_back(u);  // For an undirected graph
  }

  std::unordered_map<int, std::vector<int>> GCSBFS::BFSWithPaths(int start) {
    if (adjacency_.empty()) {
      throw std::runtime_error("Adjacency list is empty. Please set it first.");
    }

    std::queue<int> q;
    std::unordered_map<int, bool> visited;
    std::unordered_map<int, std::vector<int>> paths;

    for (const auto& kv : adjacency_) {
      int node = kv.first;
      visited[node] = false;
      paths[node] = std::vector<int>();
    }

    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
      int current = q.front();
      q.pop();

      for (const int neighbor : adjacency_[current]) {
        if (!visited[neighbor]) {
          visited[neighbor] = true;
          q.push(neighbor);

          paths[neighbor] = paths[current];
          paths[neighbor].push_back(current);
        }
      }
    }

    return paths;
  }

  void GCSBFS::SetAdjecency(const std::unordered_map<int, std::vector<int>> &adjacency) {
    adjacency_ = adjacency;
  }

}

