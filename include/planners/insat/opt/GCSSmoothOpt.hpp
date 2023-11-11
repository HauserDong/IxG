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
 * \file GCSFullOpt.hpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 10/5/23
*/

#pragma once
#ifndef GCSFULLOPT_HPP
#define GCSFULLOPT_HPP

#include <planners/insat/opt/GCSOpt.hpp>
#include <drake/common/trajectories/bspline_trajectory.h>
#include <common/insatxgcs/utils.hpp>

namespace ps {

  class GCSSmoothOpt : public GCSOpt {

  public:

    /// Sets up the GCS regions and edges
    GCSSmoothOpt(const std::vector<HPolyhedron>& regions,
           const std::vector<std::pair<int, int>>& edges_between_regions,
           int order, int continuity,
           double path_length_weight, double time_weight,
           Eigen::VectorXd& vel_lb, Eigen::VectorXd& vel_ub,
           double h_min, double h_max, double hdot_min = 1e-6,
           bool verbose=false) : GCSOpt(regions, edges_between_regions, order, h_min, h_max,
                                        path_length_weight, time_weight, vel_lb, vel_ub,
                                        verbose),
                                 continuity_(continuity),
                                 hdot_min_(hdot_min) {

      assert (continuity_ < order_);

    }

    GCSSmoothOpt (const GCSSmoothOpt &)=default;
    GCSSmoothOpt & 	operator= (const GCSSmoothOpt &)=default;
    GCSSmoothOpt (GCSSmoothOpt &&)=default;
    GCSSmoothOpt & 	operator= (GCSSmoothOpt &&)=default;

    void formulateStartPointConstraint();

  protected:

    void setupVars() override {
      // formulate edge costs and constraints
      const Eigen::MatrixX<drake::symbolic::Variable> u_control =
          drake::symbolic::MakeMatrixContinuousVariable(num_positions_, order_ + 1, "xu");
      const Eigen::MatrixX<drake::symbolic::Variable> v_control =
          drake::symbolic::MakeMatrixContinuousVariable(num_positions_, order_ + 1, "xv");
      const Eigen::VectorX<drake::symbolic::Variable> u_duration =
          drake::symbolic::MakeVectorContinuousVariable(order_ + 1, "Tu");
      const Eigen::VectorX<drake::symbolic::Variable> v_duration =
          drake::symbolic::MakeVectorContinuousVariable(order_ + 1, "Tv");

      u_vars_.resize(u_control.size() + u_duration.size());
      u_vars_ <<
              Eigen::Map<const Eigen::VectorX<drake::symbolic::Variable>>(u_control.data(), u_control.size()),
          u_duration;
      u_r_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
          drake::math::BsplineBasis<drake::symbolic::Expression>(
              order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
              ControlPointsOf(u_control));
      Eigen::MatrixX<drake::symbolic::Variable> u_duration_transpose = u_duration.transpose();
      u_h_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
          drake::math::BsplineBasis<drake::symbolic::Expression>(
              order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
              ControlPointsOf(u_duration_transpose));

      Eigen::VectorX<drake::symbolic::Variable> edge_vars(u_control.size() + u_duration.size() + v_control.size() + v_duration.size());
      edge_vars <<
                Eigen::Map<const Eigen::VectorX<drake::symbolic::Variable>>(u_control.data(), u_control.size()),
          u_duration,Eigen::Map<const Eigen::VectorX<drake::symbolic::Variable>>(v_control.data(), v_control.size()),
          v_duration;
      v_r_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
          drake::math::BsplineBasis<drake::symbolic::Expression>(order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
          ControlPointsOf(v_control));
      Eigen::MatrixX<drake::symbolic::Variable> v_duration_transpose = v_duration.transpose();
      v_h_trajectory_ = drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>(
          drake::math::BsplineBasis<drake::symbolic::Expression>(order_ + 1, order_ + 1, drake::math::KnotVectorType::kClampedUniform, 0.0, 1.0),
          ControlPointsOf(v_duration_transpose));

    }


    /// Preprocess regions to add the time scaling set and create vertices
    void preprocess(const drake::geometry::optimization::ConvexSets& regions,
                            const std::vector<std::pair<int, int>>& edges_between_regions) override {

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

    void formulateTimeCost() override {
      const std::vector<Eigen::MatrixX<drake::symbolic::Expression>> u_time_control = u_h_trajectory_.control_points();
      Eigen::VectorX<drake::symbolic::Expression> segment_time = u_time_control.back() - u_time_control.front();
      Eigen::MatrixXd M(segment_time.rows(), u_vars_.size());
      DecomposeLinearExpressions(segment_time, u_vars_, &M);
      time_cost_ = std::make_shared<drake::solvers::LinearCost>(time_weight_ * M.row(0), 0.0);
    }

    void formulatePathLengthIntegralCost(int integration_points) {
      std::vector<drake::symbolic::Expression> s_points(integration_points + 1);
      double step = 1.0 / integration_points;
      for (size_t i = 0; i <= integration_points; ++i) s_points[i] = i * step;
      auto u_path_deriv =
          utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(u_r_trajectory_.MakeDerivative(1));

      if (u_path_deriv->basis().order() == 1) {
        for (drake::symbolic::Expression t : {0.0, 1.0}) {
          Eigen::VectorX<drake::symbolic::Expression> q_ds = u_path_deriv->value(t);
          Eigen::VectorX<drake::symbolic::Expression> costs(num_positions_);
          for (size_t i = 0; i < num_positions_; ++i)
            costs(i) = q_ds(i);
          Eigen::MatrixXd H(costs.rows(), u_vars_.cols());
          DecomposeLinearExpressions(costs, u_vars_, &H);
          integral_cost_ =
              std::make_shared<drake::solvers::L2NormCost>(path_length_integral_weight_ * H, Eigen::VectorXd::Zero(num_positions_));
        }
      } else {
        Eigen::MatrixX<drake::symbolic::Expression> q_ds = u_path_deriv->vector_values(s_points);
        for (size_t i = 0; i < integration_points + 1; ++i) {
          Eigen::VectorX<drake::symbolic::Expression> costs(num_positions_);
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

    void formulatePathEnergyCost() {
      Eigen::MatrixXd weight_matrix = Eigen::MatrixXd::Identity(num_positions_, num_positions_);
      weight_matrix *= path_energy_integral_weight_;

      const std::vector<Eigen::MatrixX<drake::symbolic::Expression>> u_path_control =
          utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
              u_r_trajectory_.MakeDerivative(1))->control_points();
      const std::vector<Eigen::MatrixX<drake::symbolic::Expression>> u_time_control =
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

    void AddDerivativeRegularization(double weight_r, double weight_h, int order) {
      assert (2 <= order && order <= order_);
      std::vector<double> weights = {weight_r, weight_h};
      std::vector<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>> trajectories =
          {u_r_trajectory_, u_h_trajectory_};

      for (size_t i = 0; i < 2; ++i) {
        auto traj = trajectories[i];
        auto weight = weights[i];
        const std::vector<Eigen::MatrixX<drake::symbolic::Expression>> derivative_control =
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

    void formulateVelocityLimits(const Eigen::VectorXd& lb, const Eigen::VectorXd& ub) {
      assert (lb.size() == num_positions_);
      assert (ub.size() == num_positions_);

      const std::vector<Eigen::MatrixX<drake::symbolic::Expression>> u_path_control =
      utils::DynamicUniqueCast<drake::trajectories::BsplineTrajectory<drake::symbolic::Expression>>(
          u_r_trajectory_.MakeDerivative(1)
      )->control_points();
      const std::vector<Eigen::MatrixX<drake::symbolic::Expression>> u_time_control =
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
                     A_ctrl - ub * b_ctrl,
            -A_ctrl + lb * b_ctrl;
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



    void formulateContinuityConstraint();
    std::vector<Eigen::MatrixX<drake::symbolic::Expression>>
          ControlPointsOf(const Eigen::MatrixX<drake::symbolic::Variable>& mat);
// dynamic_unique_cast for unique_ptr
    template <typename To, typename From, typename Deleter>
    std::unique_ptr<To, Deleter> DynamicUniqueCast(std::unique_ptr<From, Deleter>&& p);

    /// Basics
    int continuity_;
    double hdot_min_;

    /// weights
    double path_length_integral_weight_;
    double path_energy_integral_weight_;

    HPolyhedron time_scaling_set_;

    /// Variables
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> u_r_trajectory_;
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> u_h_trajectory_;

    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> v_r_trajectory_;
    drake::trajectories::BsplineTrajectory<drake::symbolic::Expression> v_h_trajectory_;

    /// Costs
    std::shared_ptr<drake::solvers::Cost> time_cost_;
    std::shared_ptr<drake::solvers::Cost> integral_cost_;
    std::shared_ptr<drake::solvers::Cost> energy_cost_;
    std::vector<std::shared_ptr<drake::solvers::Cost>> reg_cost_;

    /// Constraints
    std::vector<std::pair<std::shared_ptr<drake::solvers::Constraint>, VectorXb>> continuity_constraint_;
    std::vector<std::shared_ptr<drake::solvers::Constraint>> deriv_constraints_;

  };
}

#endif //INSATxGCS_OPT_HPP
