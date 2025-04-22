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
#include <iostream>
#include <memory>
#include <cstdlib>

#include <common/insatxgcs/utils.hpp>
//#include "INSATxGCS/INSATxGCSplanner.hpp"
//#include "GCS/simple_bezier.hpp"
//#include "GCS/bezier.hpp"
#include <planners/insat/opt/GCSOpt.hpp>
#include <iomanip>

using namespace ps;
using drake::geometry::optimization::GraphOfConvexSets;
using drake::geometry::optimization::HPolyhedron;
using utils::operator<<;

std::vector<HPolyhedron> ConstructSimpleRegions() {
  Eigen::MatrixXd A_bl(3, 2), A_br(3, 2), A_tl(3, 2), A_tr(3, 2);
  Eigen::VectorXd b(3);
  A_bl <<
       -1,  0,
          0, -1,
          1,  1;
  A_br <<
       1,  0,
          0, -1,
          -1,  1;
  A_tl <<
       -1,  0,
          0,  1,
          1, -1;
  A_tr <<
       1,  0,
          0,  1,
          -1, -1;
  b <<
    3,
          3,
          -1;

  return {
          HPolyhedron(A_bl, b),
          HPolyhedron(A_br, b),
          HPolyhedron(A_tl, b),
          HPolyhedron(A_tr, b)
  };
}

std::vector<VertexId> pathToVertexId(std::vector<int64_t>& path, const std::vector<GraphOfConvexSets::Vertex*>& vertices) {
  std::vector<VertexId> solve_vids;
  for (auto vid : path) {
    for (auto& v : vertices) {
      if (v->id().get_value() == vid) {
        solve_vids.push_back(v->id());
        break;
      }
    }
  }
  return solve_vids;
}

std::vector<EdgeId > pathToEdgeId(std::vector<int64_t>& path, const std::vector<GraphOfConvexSets::Edge*>& edges) {
  std::vector<EdgeId> solve_eids;
  for (int i=0; i<path.size()-1; ++i) {
    int64_t uid = path[i];
    int64_t vid = path[i+1];
    for (auto& e : edges) {
      if (e->u().id().get_value() == uid && e->v().id().get_value() == vid) {
        solve_eids.push_back(e->id());
        break;
      }
    }
  }
  return solve_eids;
}

int main() {

  auto lic = drake::solvers::MosekSolver::AcquireLicense();

  std::vector<HPolyhedron> regions = utils::DeserializeRegions("/home/gaussian/cmu_ri_phd/phd_research/temp_INSATxGCS/INSATxGCS-Planner/src/data/maze.csv");
  auto edges_bw_regions = utils::DeserializeEdges("/home/gaussian/cmu_ri_phd/phd_research/temp_INSATxGCS/INSATxGCS-Planner/src/data/maze_edges.csv");

  int num_positions = 2;
  int order = 1;
  double h_min = 1e-3;
  double h_max = 1;
  double path_len_weight = 1;
  double time_weight = 0;
  Eigen::VectorXd vel_lb = -5 * Eigen::VectorXd::Ones(num_positions);
  Eigen::VectorXd vel_ub = 5 * Eigen::VectorXd::Ones(num_positions);
  bool verbose = false;

  GCSOpt opt(regions, *edges_bw_regions,
             order, h_min, h_max, path_len_weight, time_weight,
             vel_lb, vel_ub, verbose);

//  Eigen::VectorXd start(num_positions);
//  start << 0.5, 0;
//  auto start_vertex_id = opt.AddStart(start);
//  Eigen::VectorXd goal(num_positions);
//  goal << 49.5, 50;
//  auto goal_vertex_id = opt.AddGoal(goal);

  auto vertices = opt.GetVertices();
  opt.FormulateAndSetCostsAndConstraints();
//  std::vector<VertexId> solve_vids;
  std::vector<double> costs;

//  int64_t vid_on_path[] = {0, 50, 51};
  int64_t vid_on_path[] = {2498, 2499, 2449};
//  int64_t vid_on_path[] = {0, 50, 51, 1, 2, 3, 53, 103};

//  int64_t vid_on_path[] = {2, 52, 53, 3};
//  int64_t vid_on_path[] = {1, 51, 52};
//  int64_t vid_on_path[] = {1, 51, 52, 2};
//  int64_t vid_on_path[] = {1, 51, 52, 2, 3, 4, 54, 104};
//  int64_t vid_on_path[] = {1, 51, 52, 2, 3, 4, 54, 104, 154, 204, 203, 202, 252, 302, 352 };
//  int64_t vid_on_path[] = {1, 51, 52, 2, 3, 4, 54, 104, 154, 204, 203, 202, 252, 302, 352, 351, 401, 451, 501, 551, 552, 553, 554, 504, 505, 506, 556, 557, 507, 457, 407, 406, 456, 455, 405, 355, 356, 306, 305, 304, 303, 253, 254, 255, 205, 206, 256, 257, 207, 157, 156, 155, 105, 106, 56, 57, 7, 8, 9, 59, 58, 108, 109, 110, 111, 112, 62, 61, 11, 12, 13, 63, 64, 14, 15, 16, 66, 116, 115, 114, 164, 163, 162, 161, 211, 210, 260, 259, 209, 208, 258, 308, 358, 408, 409, 410, 411, 461, 462, 512};
//  int64_t vid_on_path[] = {1, 51, 52, 2, 3, 4, 54, 104, 154, 204, 203, 202, 252, 302, 352, 351, 401, 451, 501, 551, 552, 553, 554, 504, 505, 506, 556, 557, 507, 457, 407, 406, 456, 455, 405, 355, 356, 306, 305, 304, 303, 253, 254, 255, 205, 206, 256, 257, 207, 157, 156, 155, 105, 106, 56, 57, 7, 8, 9, 59, 58, 108, 109, 110, 111, 112, 62, 61, 11, 12, 13, 63, 64, 14, 15, 16, 66, 116, 115, 114, 164, 163, 162, 161, 211, 210, 260, 259, 209, 208, 258, 308, 358, 408, 409, 410, 411, 461, 462, 512, 513, 463, 464, 414, 415, 416, 466, 516, 517, 518, 468, 467, 417, 418, 419, 369, 319, 320, 270, 271, 321, 322, 372, 422, 472, 471, 470, 520, 570, 571, 572, 622, 623, 624, 574, 524, 525, 575, 576, 626, 676, 726, 776, 777, 827, 877, 878, 828, 829, 830, 831, 881, 882, 932, 982, 1032, 1033, 1083, 1133, 1134, 1135, 1185, 1184, 1183, 1233, 1283, 1333, 1332, 1382, 1383, 1433, 1483, 1533, 1583, 1633, 1683, 1684, 1734, 1733, 1732, 1682, 1632, 1631, 1681, 1680, 1730, 1780, 1781, 1782, 1783, 1833, 1834, 1835, 1885, 1935, 1936, 1937, 1987, 2037, 2087, 2086, 2136, 2186, 2185, 2235, 2285, 2286, 2287, 2237, 2238, 2239, 2289, 2290, 2340, 2341, 2342, 2343, 2344, 2294, 2244, 2245, 2195, 2196, 2197, 2247, 2246, 2296, 2346, 2396, 2446, 2496, 2497, 2447, 2448, 2449, 2499};
//  int64_t vid_on_path[] = {start_vertex_id.get_value(), 1, 51, 52, 2, 3, 4, 54, 104, 154, 204, 203, 202, 252, 302, 352, 351, 401, 451, 501, 551, 552, 553, 554, 504, 505, 506, 556, 557, 507, 457, 407, 406, 456, 455, 405, 355, 356, 306, 305, 304, 303, 253, 254, 255, 205, 206, 256, 257, 207, 157, 156, 155, 105, 106, 56, 57, 7, 8, 9, 59, 58, 108, 109, 110, 111, 112, 62, 61, 11, 12, 13, 63, 64, 14, 15, 16, 66, 116, 115, 114, 164, 163, 162, 161, 211, 210, 260, 259, 209, 208, 258, 308, 358, 408, 409, 410, 411, 461, 462, 512, 513, 463, 464, 414, 415, 416, 466, 516, 517, 518, 468, 467, 417, 418, 419, 369, 319, 320, 270, 271, 321, 322, 372, 422, 472, 471, 470, 520, 570, 571, 572, 622, 623, 624, 574, 524, 525, 575, 576, 626, 676, 726, 776, 777, 827, 877, 878, 828, 829, 830, 831, 881, 882, 932, 982, 1032, 1033, 1083, 1133, 1134, 1135, 1185, 1184, 1183, 1233, 1283, 1333, 1332, 1382, 1383, 1433, 1483, 1533, 1583, 1633, 1683, 1684, 1734, 1733, 1732, 1682, 1632, 1631, 1681, 1680, 1730, 1780, 1781, 1782, 1783, 1833, 1834, 1835, 1885, 1935, 1936, 1937, 1987, 2037, 2087, 2086, 2136, 2186, 2185, 2235, 2285, 2286, 2287, 2237, 2238, 2239, 2289, 2290, 2340, 2341, 2342, 2343, 2344, 2294, 2244, 2245, 2195, 2196, 2197, 2247, 2246, 2296, 2346, 2396, 2446, 2496, 2497, 2447, 2448, 2449, 2499, 2500, goal_vertex_id.get_value()};
//  int64_t vid_on_path[] = {1, 52, 2, 3, 54, 104, 154, 153, 203, 252, 302, 352, 401, 451, 501, 502, 552, 553, 554, 504, 505, 506, 507, 457, 456, 406, 456, 406, 405, 355, 356, 355, 305, 304, 254, 255, 206, 256, 206, 207, 206, 156, 106, 57, 8, 59, 58, 109, 110, 111, 61, 12, 13, 63, 64, 14, 15, 65, 115, 114, 163, 162, 211, 210, 259, 209, 258, 308, 358, 409, 410, 461, 462, 463, 464, 415, 466, 467, 418, 369, 320, 321, 371, 372, 422, 421, 471, 520, 571, 572, 573, 574, 575, 626, 676, 726, 777, 827, 828, 829, 830, 881, 932, 982, 1033, 1083, 1134, 1184, 1233, 1283, 1333, 1383, 1433, 1483, 1533, 1583, 1633, 1683, 1682, 1681, 1730, 1781, 1782, 1833, 1834, 1885, 1936, 1986, 1987, 2037, 2036, 2086, 2136, 2186, 2235, 2236, 2286, 2287, 2237, 2238, 2289, 2290, 2291, 2341, 2342, 2343, 2293, 2294, 2245, 2196, 2246, 2296, 2346, 2396, 2446, 2447, 2448, 2499, 2500};
//  int64_t vid_on_path[] = {start_vertex_id.get_value(), 1, 2, 3, 54, 104, 154, 153, 203, 202, 252, 302, 301, 351, 401, 451, 501, 552, 553, 503, 504, 505, 455, 506, 507, 457, 406, 405, 354, 355, 305, 304, 253, 254, 204, 205, 206, 156, 105, 104, 105, 55, 56, 6, 7, 8, 58, 109, 110, 111, 61, 10, 11, 12, 63, 13, 14, 15, 65, 115, 114, 113, 163, 162, 111, 161, 160, 210, 209, 158, 208, 207, 258, 308, 358, 409, 410, 461, 512, 462, 463, 413, 414, 364, 415, 466, 467, 416, 417, 367, 368, 318, 319, 269, 270, 271, 321, 372, 422, 421, 471, 470, 520, 521, 572, 573, 523, 524, 474, 525, 575, 626, 676, 726, 725, 776, 827, 828, 778, 829, 830, 780, 831, 881, 932, 982, 981, 1032, 1083, 1084, 1135, 1185, 1236, 1286, 1337, 1387, 1438, 1489, 1439, 1440, 1390, 1440, 1491, 1541, 1592, 1593, 1644, 1594, 1595, 1545, 1546, 1597, 1648, 1698, 1697, 1747, 1746, 1796, 1795, 1794, 1793, 1792, 1842, 1892, 1891, 1941, 1992, 1993, 1994, 1944, 1945, 1996, 2046, 2045, 2096, 2046, 2047, 1997, 1998, 1948, 1949, 1999, 2049, 2048, 2098, 2097, 2148, 2198, 2197, 2248, 2198, 2199, 2149, 2200, 2250, 2300, 2350, 2400, 2450, 2500, goal_vertex_id.get_value()};
//  int64_t vid_on_path[] = {7699, 1, 51, 52, 2, 3, 4, 54, 104, 154, 204, 203, 202, 252, 302, 301, 351, 401, 451, 501, 551, 552, 553, 554, 504, 505, 506, 556, 557, 507, 457, 407, 406, 456, 455, 405};
//  int64_t vid_on_path[] = {355, 356, 306, 305, 304, 303, 253, 254, 255, 205, 206, 256, 257, 207, 157, 156, 155, 105, 106, 56, 57, 7, 8, 9, 59, 58, 108, 109, 110, 111, 112, 62, 61, 11, 12, 13, 63, 64, 14, 15, 16, 66, 116, 115, 114, 164, 163, 162, 161, 211, 210, 260, 259, 209, 208, 258, 308, 358, 408, 409, 410, 411, 461, 462, 512, 513, 463, 464, 414, 415, 416, 466, 516, 517, 518, 468, 467, 417, 418, 419, 369, 319, 320, 270, 271, 321, 322, 372, 422, 472, 471, 470, 520, 570, 571, 572, 622, 623, 624, 574, 524, 525, 575, 576, 626, 676, 726, 776, 777, 827, 877, 878, 828, 829, 830, 831, 881, 882, 932, 982, 1032, 1033, 1083, 1133, 1134, 1135, 1185, 1184, 1183, 1233, 1283, 1333, 1332, 1382, 1383, 1433, 1483, 1533, 1583, 1584, 1534, 1535, 1485, 1435, 1385, 1386, 1436, 1486, 1487, 1537, 1587, 1586, 1585, 1635, 1685, 1735, 1785, 1835, 1885, 1935, 1936, 1937, 1987, 2037, 2087, 2086, 2136, 2186, 2185, 2235, 2285, 2286, 2287, 2237, 2238, 2239, 2289, 2290, 2340, 2341, 2342, 2343, 2344, 2294, 2244, 2245, 2195, 2196, 2197, 2247, 2246, 2296, 2346, 2396, 2446, 2496, 2497, 2447, 2448, 2449, 2499, 2500, 7702};
//  int64_t vid_on_path[] = {7699, 1, 51, 52, 2, 3, 4, 54, 104, 154, 204, 203, 202, 252, 302, 301, 351, 401, 451, 501, 551, 552, 553, 554, 504, 505, 506, 556, 557, 507, 457, 407, 406, 456, 455, 405, 355, 356, 306, 305, 304, 303, 253, 254, 255, 205, 206, 256, 257, 207, 157, 156, 155, 105, 106, 56, 57, 7, 8, 9, 59, 58, 108, 109, 110, 111, 112, 62, 61, 11, 12, 13, 63, 64, 14, 15, 16, 66, 116, 115, 114, 164, 163, 162, 161, 211, 210, 260, 259, 209, 208, 258, 308, 358, 408, 409, 410, 411, 461, 462, 512, 513, 463, 464, 414, 415, 416, 466, 516, 517, 518, 468, 467, 417, 418, 419, 369, 319, 320, 270, 271, 321, 322, 372, 422, 472, 471, 470, 520, 570, 571, 572, 622, 623, 624, 574, 524, 525, 575, 576, 626, 676, 726, 776, 777, 827, 877, 878, 828, 829, 830, 831, 881, 882, 932, 982, 1032, 1033, 1083, 1133, 1134, 1135, 1185, 1184, 1183, 1233, 1283, 1333, 1332, 1382, 1383, 1433, 1483, 1533, 1583, 1584, 1534, 1535, 1485, 1435, 1385, 1386, 1436, 1486, 1487, 1537, 1587, 1586, 1585, 1635, 1685, 1735, 1785, 1835, 1885, 1935, 1936, 1937, 1987, 2037, 2087, 2086, 2136, 2186, 2185, 2235, 2285, 2286, 2287, 2237, 2238, 2239, 2289, 2290, 2340, 2341, 2342, 2343, 2344, 2294, 2244, 2245, 2195, 2196, 2197, 2247, 2246, 2296, 2346, 2396, 2446, 2496, 2497, 2447, 2448, 2449, 2499, 2500, 7702};

// std::vector<int64_t> path (vid_on_path, vid_on_path + sizeof(vid_on_path) / sizeof(int64_t ) );
//
//  for (auto vid : path) {
//    for (auto& v : vertices) {
//      if (v->id().get_value() == vid) {
//        solve_vids.push_back(v->id());
//        break;
//      }
//    }
//  }
//
//  opt.FormulateAndSetCostsAndConstraints();
//  auto start_time = std::chrono::high_resolution_clock::now();
//  auto soln = opt.Solve(solve_vids);
//  auto end_time = std::chrono::high_resolution_clock::now();
//
//  if (soln.second.is_success()) {
//    std::cout << "Solved from scratch in " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;
//
//    auto init_guess = soln.second.get_x_val();
////    init_guess.bottomRows((order+1)*num_positions + 1) = Eigen::VectorXd::Zero((order+1)*num_positions + 1);
//    start_time = std::chrono::high_resolution_clock::now();
//    auto soln_w_init = opt.Solve(solve_vids, init_guess);
//    end_time = std::chrono::high_resolution_clock::now();
//    std::cout << "Solved with init guess in " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;
//
//    double path_length = 0;
//    auto final_trajectory = soln.first;
//    auto p0 = final_trajectory.value(0);
//    for (double t=0; t<final_trajectory.end_time(); t+=1e-2) {
////    std::cout << "t: " << t << " " << final_trajectory.value(t).transpose() << std::endl;
////      std::cout << final_trajectory.value(t).transpose() << std::endl;
//      std::cout << final_trajectory.value(t)(0) << " " << final_trajectory.value(t)(1) << std::endl;
//      path_length += (final_trajectory.value(t) - p0).norm();
//      p0 = final_trajectory.value(t);
//    }
//    std::cout << "path length: " << path_length << std::endl;
//
////    std::cout << "Solution cost: " << soln.second.get_optimal_cost() << std::endl;
////    std::cout << "path size: "<< path.size() << " init guess size: " << init_guess.rows() << " " << init_guess.cols() << std::endl;
//
//  } else {
//    std::cerr << "Solution NOT found..." << std::endl;
//    std::runtime_error("Solution NOT found...");
//  }
//  std::cout << "convergence status: " << soln.second.get_solution_result() << std::endl;

  int sz = static_cast<int>(sizeof(vid_on_path) / sizeof(int64_t ));
  for (int i=3; i<=sz; ++i) {

    std::vector<int64_t> path (vid_on_path, vid_on_path + i);
    std::vector<VertexId> solve_vids;

    for (auto vid : path) {
      for (auto& v : vertices) {
        if (v->id().get_value() == vid) {
          solve_vids.push_back(v->id());
          break;
        }
      }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    auto soln = opt.Solve(solve_vids);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (soln.second.is_success()) {
      std::cout << "Solved from scratch in " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << "s" << std::endl;

      double path_length = 0;
      auto final_trajectory = soln.first;
      auto p0 = final_trajectory.value(0);
      for (double t=0; t<final_trajectory.end_time(); t+=1e-2) {
        path_length += (final_trajectory.value(t) - p0).norm();
        p0 = final_trajectory.value(t);
      }
      costs.push_back(path_length);
    } else {
      std::cerr << "Solution NOT found..." << std::endl;
      std::runtime_error("Solution NOT found...");
    }
    std::cout << sz << " convergence status: " << soln.second.get_solution_result() << std::endl;
  }

  for (auto c : costs) {
    std::cout << c << ", ";
  }
  std::cout << std::endl;

  return 0;
}

