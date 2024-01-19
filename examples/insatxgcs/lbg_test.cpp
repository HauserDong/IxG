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


  std::vector<HPolyhedron> regions = utils::DeserializeRegions("/home/gaussian/cmu_ri_phd/phd_research/temp_INSATxGCS/INSATxGCS-Planner/src/data/maze.csv");
  auto edges_bw_regions = utils::DeserializeEdges("/home/gaussian/cmu_ri_phd/phd_research/temp_INSATxGCS/INSATxGCS-Planner/src/data/maze_edges.csv");

  int num_positions = 2;
  int order = 2;
  int continuity = 1;
  double h_min = 1e-3;
  double h_max = 1;
  double path_len_weight = 1;
  double time_weight = 0;
  Eigen::VectorXd vel_lb = -5 * Eigen::VectorXd::Ones(num_positions);
  Eigen::VectorXd vel_ub = 5 * Eigen::VectorXd::Ones(num_positions);
  bool verbose = false;

  LBGraph lbg(regions, *edges_bw_regions,
             order, continuity,
             path_len_weight, time_weight,
             vel_lb, vel_ub,
              h_min, h_max);

  lbg.PrintLBGraphStats();

}
