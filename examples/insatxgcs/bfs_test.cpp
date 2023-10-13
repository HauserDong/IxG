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
 * \file bfs_test.cpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 10/11/23
*/

#include <common/insatxgcs/gcsbfs.hpp>

auto printPath = [](int goal, const std::vector<int>& v) {
  std::cout << "goal: " << goal << std::endl;
  std::cout << "path: ";
  for (int num : v) {
    std::cout << num << " -> ";
  }
  std::cout << std::endl;
};

int main() {
  auto edges_bw_regions = utils::DeserializeEdges("../examples/insatxgcs/resources/maze2d/maze_edges.csv");
  ixg::GCSBFS graph(*edges_bw_regions);

  int start_node = 0;

  auto start_time = std::chrono::high_resolution_clock::now();
  std::unordered_map<int, std::vector<int>> paths = graph.BFSWithPaths(start_node);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::cout << "BFS completed in " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;

  std::cout << "Paths from node " << start_node << " to other nodes:" << std::endl;
  for (const auto& kv : paths) {

    printPath(kv.first, kv.second);
//    printPath(0, paths[0]);

    break;
  }

  return 0;
}
