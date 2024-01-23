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
    const drake::VectorX<drake::symbolic::Variable> edge_vars =
        drake::solvers::ConcatenateVariableRefList({u_vars_, v_vars_});

    for (int deriv = 0; deriv < continuity_ + 1; ++deriv) {
      /// path continuity
      auto u_path_deriv = DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>
          (u_r_trajectory_.MakeDerivative(deriv));
      auto v_path_deriv = DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>
          (v_r_trajectory_.MakeDerivative(deriv));

      const drake::VectorX<drake::symbolic::Expression> path_continuity_error =
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

      const drake::VectorX<drake::symbolic::Expression> time_continuity_error =
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



  std::vector<drake::MatrixX<drake::symbolic::Expression>>
  ps::GCSSmoothOpt::ControlPointsOf(const drake::MatrixX<drake::symbolic::Variable> &mat) {
    std::vector<drake::MatrixX<drake::symbolic::Expression>> vec;
    for (int i = 0; i < mat.cols(); ++i)
      vec.push_back(mat.col(i));
    return vec;
  }

  void GCSSmoothOpt::setupVars() {
    // formulate edge costs and constraints
    const drake::MatrixX<drake::symbolic::Variable> u_control =
        drake::symbolic::MakeMatrixContinuousVariable(num_positions_, order_ + 1, "xu");
    const drake::MatrixX<drake::symbolic::Variable> v_control =
        drake::symbolic::MakeMatrixContinuousVariable(num_positions_, order_ + 1, "xv");
    const drake::VectorX<drake::symbolic::Variable> u_duration =
        drake::symbolic::MakeVectorContinuousVariable(order_ + 1, "Tu");
    const drake::VectorX<drake::symbolic::Variable> v_duration =
        drake::symbolic::MakeVectorContinuousVariable(order_ + 1, "Tv");

    u_vars_.resize(u_control.size() + u_duration.size());
    u_vars_ <<
            Eigen::Map<const drake::VectorX<drake::symbolic::Variable>>(u_control.data(), u_control.size()),
        u_duration;
    u_r_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
        drake::math::BsplineBasis<drake::symbolic::Expression>(
            order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
        ControlPointsOf(u_control));
    drake::MatrixX<drake::symbolic::Variable> u_duration_transpose = u_duration.transpose();
    u_h_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
        drake::math::BsplineBasis<drake::symbolic::Expression>(
            order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
        ControlPointsOf(u_duration_transpose));

    drake::VectorX<drake::symbolic::Variable> edge_vars(u_control.size() + u_duration.size() + v_control.size() + v_duration.size());
    edge_vars <<
              Eigen::Map<const drake::VectorX<drake::symbolic::Variable>>(u_control.data(), u_control.size()),
        u_duration,Eigen::Map<const drake::VectorX<drake::symbolic::Variable>>(v_control.data(), v_control.size()),
        v_duration;
    v_r_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
        drake::math::BsplineBasis<drake::symbolic::Expression>(order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
        ControlPointsOf(v_control));
    drake::MatrixX<drake::symbolic::Variable> v_duration_transpose = v_duration.transpose();
    v_h_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
        drake::math::BsplineBasis<drake::symbolic::Expression>(order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
        ControlPointsOf(v_duration_transpose));

  }

  void GCSSmoothOpt::preprocess(const drake::geometry::optimization::ConvexSets &regions,
                                const std::vector<std::pair<int, int>> &edges_between_regions) {

    Eigen::MatrixXd A_time(3 * order_ + 2, order_ + 1);
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(order_, order_ + 1);
    tmp.block(0, 1, order_, order_).setIdentity();
    A_time <<
           Eigen::MatrixXd::Identity(order_ + 1, order_ + 1),
        -Eigen::MatrixXd::Identity(order_ + 1, order_ + 1),
        Eigen::MatrixXd::Identity(order_,     order_ + 1) - tmp;
    Eigen::VectorXd b_time(3 * order_ + 2);
    b_time <<
           Eigen::VectorXd::Ones(order_ + 1) * 1e3,
        Eigen::VectorXd::Zero(order_ + 1),
        -Eigen::VectorXd::Ones(order_) * hdot_min_;
    time_scaling_set_ = HPolyhedron(A_time, b_time);

    for (size_t i = 0; i < hpoly_regions_.size(); ++i) {
      vertices_.emplace_back(gcs_->AddVertex(
          hpoly_regions_[i].CartesianPower(order_ + 1).CartesianProduct(time_scaling_set_),
          fmt::format("{}: {}", "v" + std::to_string(i), i)));
      vertex_id_to_vertex_[vertices_.back()->id().get_value()] = vertices_.back();
      vertex_id_to_regions_[vertices_.back()->id().get_value()] = hpoly_regions_[i];
    }

  }

  void GCSSmoothOpt::formulateTimeCost() {
    const std::vector<drake::MatrixX<drake::symbolic::Expression>> u_time_control = u_h_trajectory_.control_points();
    drake::VectorX<drake::symbolic::Expression> segment_time = u_time_control.back() - u_time_control.front();
    Eigen::MatrixXd M(segment_time.rows(), u_vars_.size());
    DecomposeLinearExpressions(segment_time, u_vars_, &M);
    time_cost_ = std::make_shared<drake::solvers::LinearCost>(time_weight_ * M.row(0), 0.0);
  }

  void GCSSmoothOpt::formulatePathLengthIntegralCost(int integration_points) {
    std::vector<drake::symbolic::Expression> s_points(integration_points + 1);
    double step = 1.0 / integration_points;
    for (size_t i = 0; i <= integration_points; ++i) s_points[i] = i * step;
    auto u_path_deriv =
        utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(u_r_trajectory_.MakeDerivative(1));

    if (u_path_deriv->basis().order() == 1) {
      for (drake::symbolic::Expression t : {0.0, 1.0}) {
        drake::VectorX<drake::symbolic::Expression> q_ds = u_path_deriv->value(t);
        drake::VectorX<drake::symbolic::Expression> costs(num_positions_);
        for (size_t i = 0; i < num_positions_; ++i)
          costs(i) = q_ds(i);
        Eigen::MatrixXd H(costs.rows(), u_vars_.cols());
        DecomposeLinearExpressions(costs, u_vars_, &H);
        integral_cost_ =
            std::make_shared<drake::solvers::L2NormCost>(path_length_integral_weight_ * H, Eigen::VectorXd::Zero(num_positions_));
      }
    } else {
      drake::MatrixX<drake::symbolic::Expression> q_ds = u_path_deriv->vector_values(s_points);
      for (size_t i = 0; i < integration_points + 1; ++i) {
        drake::VectorX<drake::symbolic::Expression> costs(num_positions_);
        for (size_t j = 0; j < num_positions_; ++j) {
          costs(j) = (i == 0 || i == integration_points)
                     ? 0.5 * 1.0 / integration_points * q_ds(j, i)
                     : 1.0 / integration_points * q_ds(j, i);
        }
        Eigen::MatrixXd H(costs.rows(), u_vars_.cols());
        DecomposeLinearExpressions(costs, u_vars_, &H);
        integral_cost_ =
            std::make_shared<drake::solvers::L2NormCost>(path_length_integral_weight_ * H, Eigen::VectorXd::Zero(num_positions_));
      }
    }
  }

  void GCSSmoothOpt::formulatePathEnergyCost() {
    Eigen::MatrixXd weight_matrix = Eigen::MatrixXd::Identity(num_positions_, num_positions_);
    weight_matrix *= path_energy_integral_weight_;

    const std::vector<drake::MatrixX<drake::symbolic::Expression>> u_path_control =
        utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
            u_r_trajectory_.MakeDerivative(1))->control_points();
    const std::vector<drake::MatrixX<drake::symbolic::Expression>> u_time_control =
        utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
            u_h_trajectory_.MakeDerivative(1))->control_points();
    for (size_t i = 0; i < u_path_control.size(); ++i) {
      Eigen::MatrixXd A_ctrl(u_path_control[i].rows(), u_vars_.size());
      DecomposeLinearExpressions(u_path_control[i], u_vars_, &A_ctrl);
      Eigen::MatrixXd b_ctrl(u_time_control[i].rows(), u_vars_.size());
      DecomposeLinearExpressions(u_time_control[i], u_vars_, &b_ctrl);
      assert (A_ctrl.cols() == b_ctrl.cols());
      Eigen::MatrixXd H(A_ctrl.rows() + b_ctrl.rows(), A_ctrl.cols());
      H << order_ * b_ctrl, weight_matrix.array().sqrt().matrix() * A_ctrl;
      energy_cost_ = std::make_shared<drake::solvers::PerspectiveQuadraticCost>(H, Eigen::VectorXd::Zero(H.rows()));
    }
  }

  void GCSSmoothOpt::AddDerivativeRegularization(double weight_r, double weight_h, int order) {
    assert (2 <= order && order <= order_);
    std::vector<double> weights = {weight_r, weight_h};
    std::vector<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>> trajectories =
        {u_r_trajectory_, u_h_trajectory_};

    for (size_t i = 0; i < 2; ++i) {
      auto traj = trajectories[i];
      auto weight = weights[i];
      const std::vector<drake::MatrixX<drake::symbolic::Expression>> derivative_control =
          utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
              traj.MakeDerivative(order))->control_points();

      for (auto & c : derivative_control) {
        Eigen::MatrixXd A_ctrl(c.rows(), u_vars_.size());
        DecomposeLinearExpressions(c, u_vars_, &A_ctrl);
        Eigen::MatrixXd H = A_ctrl.transpose() * A_ctrl * 2 * weight / (1 + order_ - order);
        reg_cost_.emplace_back(std::make_shared<drake::solvers::QuadraticCost>(H, Eigen::VectorXd::Zero(H.rows()), 0));
      }
    }
  }

  void GCSSmoothOpt::formulateVelocityConstraint() {
    assert (vel_lb_.size() == num_positions_);
    assert (vel_ub_.size() == num_positions_);

    const std::vector<drake::MatrixX<drake::symbolic::Expression>> u_path_control =
        utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
            u_r_trajectory_.MakeDerivative(1)
        )->control_points();
    const std::vector<drake::MatrixX<drake::symbolic::Expression>> u_time_control =
        utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
            u_h_trajectory_.MakeDerivative(1))->control_points();

    for (size_t i = 0; i < u_path_control.size(); ++i) {
      Eigen::MatrixXd A_ctrl(u_path_control[i].rows(), u_vars_.size());
      DecomposeLinearExpressions(u_path_control[i], u_vars_, &A_ctrl);
      Eigen::MatrixXd b_ctrl(u_time_control[i].rows(), u_vars_.size());
      DecomposeLinearExpressions(u_time_control[i], u_vars_, &b_ctrl);

      assert (A_ctrl.cols() == b_ctrl.cols());
      Eigen::MatrixXd A_constraint(A_ctrl.rows() * 2, A_ctrl.cols());
      A_constraint <<
                   A_ctrl - vel_ub_ * b_ctrl,
          -A_ctrl + vel_lb_ * b_ctrl;
      // A_constraint <<
      //     A_ctrl - (b_ctrl.array().colwise() * ub.array()).matrix(),
      //    -A_ctrl + (b_ctrl.array().colwise() * lb.array()).matrix();
      std::shared_ptr<drake::solvers::LinearConstraint> velocity_con =
          std::make_shared<drake::solvers::LinearConstraint>(
              A_constraint,
              -Eigen::VectorXd::Ones(2*num_positions_) * std::numeric_limits<double>::infinity(),
              Eigen::VectorXd::Zero(2*num_positions_)
          );
      deriv_constraints_.emplace_back(velocity_con);
    }
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