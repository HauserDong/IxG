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

    void setupVars() override;


    /// Preprocess regions to add the time scaling set and create vertices
    void preprocess(const drake::geometry::optimization::ConvexSets& regions,
                            const std::vector<std::pair<int, int>>& edges_between_regions) override;

    void formulateTimeCost() override;
    void formulatePathLengthIntegralCost(int integration_points);
    void formulatePathEnergyCost();
    void formulateVelocityConstraint() override;
    void AddDerivativeRegularization(double weight_r, double weight_h, int order);



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
