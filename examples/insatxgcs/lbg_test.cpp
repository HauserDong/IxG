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
 * \file lbg_test.cpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 1/17/24
*/

#include <iostream>
#include <memory>
#include <cstdlib>

#include <common/insatxgcs/utils.hpp>
#include <planners/insat/opt/LBGraph.hpp>

using namespace ps;
using drake::geometry::optimization::GraphOfConvexSets;
using drake::geometry::optimization::HPolyhedron;
using utils::operator<<;

int main() {
  setenv("MOSEKLM_LICENSE_FILE", "/home/gaussian/Documents/softwares/mosektoolslinux64x86/mosek.lic", true);
  auto lic = drake::solvers::MosekSolver::AcquireLicense();

  std::string env_name = "maze2d";
  std::vector<HPolyhedron> regions = utils::DeserializeRegions("/home/gaussian/cmu_ri_phd/phd_research/temp_INSATxGCS/INSATxGCS-Planner/src/data/maze.csv");
  auto edges_bw_regions = utils::DeserializeEdges("/home/gaussian/cmu_ri_phd/phd_research/temp_INSATxGCS/INSATxGCS-Planner/src/data/maze_edges.csv");
  std::string lbg_dir = "../examples/insatxgcs/resources/" + env_name + "/lbg/";
  std::string lbg_file = lbg_dir + "maze2d_o1_c0_pw1_tw0_hmin0_hmax1_hdmin0_sm0";

  int num_positions = 2;
  int order = 1;
  int continuity = 0;
  double h_min = 1e-3;
  double h_max = 1;
  double path_len_weight = 1;
  double time_weight = 0;
  Eigen::VectorXd vel_lb = -5 * Eigen::VectorXd::Ones(num_positions);
  Eigen::VectorXd vel_ub = 5 * Eigen::VectorXd::Ones(num_positions);
  double hdot_min = 1e-6;
  bool smooth = false;
  bool verbose = false;

  LBGraph lbg(env_name, regions, *edges_bw_regions,
              order, continuity,
              path_len_weight, time_weight,
              vel_lb, vel_ub,
              h_min, h_max, hdot_min,
              smooth, lbg_dir);

  lbg.PrintLBGraphStats();
  auto graph = lbg.GetLBAdjacencyList();
  auto costs = lbg.GetLBAdjacencyCostList();

  int gcs_start_id = 0;
  auto start_state = regions[gcs_start_id].ChebyshevCenter();
  std::vector start(start_state.data(), start_state.data() + start_state.size());

  LBGSearch search(lbg_file);

  auto start_time = std::chrono::high_resolution_clock::now();
  auto graph_dist = search.Dijkstra(start, gcs_start_id);
  auto end_time = std::chrono::high_resolution_clock::now();

  std::cout << "Dijkstra took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;

  for (auto d : graph_dist) {
    if (d.second < 1e3) {
      std::cout << d.first << "\t" << d.second << std::endl;
    } else {
      std::cout << d.first << "\t" << "> 1e3" << std::endl;
    }
  }

}
