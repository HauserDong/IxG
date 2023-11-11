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
 * \file GCSSmoothOpt.cpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 10/5/23
*/

#include <planners/insat/opt/GCSSmoothOpt.hpp>

#include <iostream>

namespace ps {


  void ps::GCSSmoothOpt::formulateContinuityConstraint() {
    const Eigen::VectorX<drake::symbolic::Variable> edge_vars =
        drake::solvers::ConcatenateVariableRefList({u_vars_, v_vars_});

    for (int deriv = 0; deriv < continuity_ + 1; ++deriv) {
      /// path continuity
      auto u_path_deriv = DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>
          (u_r_trajectory_.MakeDerivative(deriv));
      auto v_path_deriv = DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>
          (v_r_trajectory_.MakeDerivative(deriv));

      const Eigen::VectorX<drake::symbolic::Expression> path_continuity_error =
          v_path_deriv->control_points().front() -
          u_path_deriv->control_points().back();
      Eigen::MatrixXd M(num_positions_, edge_vars.size());
      drake::symbolic::DecomposeLinearExpressions(path_continuity_error, edge_vars, &M);
      // Condense M to only keep non-zero columns.
      const auto &[condensed_matrices, nonzero_cols_mask] =
          utils::CondenseToNonzeroColumns({M});
      Eigen::MatrixXd M_dense = condensed_matrices[0];
      continuity_constraint_.emplace_back(std::make_shared<drake::solvers::LinearEqualityConstraint>(
          M_dense, Eigen::VectorXd::Zero(num_positions_)), nonzero_cols_mask);

      /// time continuity
      auto u_time_deriv = DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>
          (u_h_trajectory_.MakeDerivative(deriv));
      auto v_time_deriv = DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>
          (v_h_trajectory_.MakeDerivative(deriv));

      const Eigen::VectorX<drake::symbolic::Expression> time_continuity_error =
          v_time_deriv->control_points().front() -
          u_time_deriv->control_points().back();
      Eigen::MatrixXd M_h(num_positions_, edge_vars.size());
      drake::symbolic::DecomposeLinearExpressions(time_continuity_error, edge_vars, &M_h);
      // Condense M to only keep non-zero columns.
      const auto &[condensed_matrices_h, nonzero_cols_mask_h] =
          utils::CondenseToNonzeroColumns({M_h});
      Eigen::MatrixXd M_h_dense = condensed_matrices_h[0];
      continuity_constraint_.emplace_back(std::make_shared<drake::solvers::LinearEqualityConstraint>(
          M_h_dense, Eigen::VectorXd::Zero(num_positions_)), nonzero_cols_mask_h);
    }
  }



  std::vector<Eigen::MatrixX<drake::symbolic::Expression>>
  ps::GCSSmoothOpt::ControlPointsOf(const Eigen::MatrixX<drake::symbolic::Variable> &mat) {
    std::vector<Eigen::MatrixX<drake::symbolic::Expression>> vec;
    for (int i = 0; i < mat.cols(); ++i)
      vec.push_back(mat.col(i));
    return vec;
  }

  template<typename To, typename From, typename Deleter>
  std::unique_ptr<To, Deleter> ps::GCSSmoothOpt::DynamicUniqueCast(std::unique_ptr<From, Deleter> &&p) {
    if (To *cast = dynamic_cast<To *>(p.get())) {
      std::unique_ptr<To, Deleter> result(cast, std::move(p.get_deleter()));
      p.release();
      return result;
    }
    // return std::unique_ptr<To, Deleter>(nullptr); // or throw std::bad_cast() if you prefer
    throw std::runtime_error("dynamic_unique_cast failed");
  }

}