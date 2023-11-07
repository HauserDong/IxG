#ifndef INSATxGCS_PLANNER_HPP
#define INSATxGCS_PLANNER_HPP

#include <future>
#include <utility>
#include "planners/Planner.hpp"
#include <common/insat/InsatState.hpp>
#include <common/insat/InsatEdge.hpp>

namespace ps
{

  class INSATxGCS : virtual public Planner
  {
  public:

    // Typedefs
    typedef std::unordered_map<size_t, InsatStatePtrType> InsatStatePtrMapType;
    typedef smpl::intrusive_heap<InsatState, IsLesserState> InsatStateQueueMinType;

    INSATxGCS(ParamsType planner_params);;

    ~INSATxGCS() {};

    void SetStartState(const StateVarsType& state_vars);

    bool Plan();

    TrajType getSolutionTraj();

  protected:
    void initialize();

    void calculateBounds();

    std::vector<InsatStatePtrType> getStateAncestors(const InsatStatePtrType state_ptr, bool reverse=false) const;

    void expandState(InsatStatePtrType state_ptr);

    void updateState(InsatStatePtrType& state_ptr,
                     std::vector<InsatStatePtrType>& ancestors,
                     InsatActionPtrType& action_ptr,
                     ActionSuccessor& action_successor);

    void constructInsatActions();

    InsatStatePtrType constructInsatState(const StateVarsType& state);

    InsatStatePtrType constructInsatPath(std::vector<InsatStatePtrType> &ancestors, const StateVarsType& state);

    void cleanUp();

    void resetStates();

    void constructPlan(InsatStatePtrType& insat_state_ptr);

    void exit();

    /// Temporary
    void printPath(std::vector<InsatStatePtrType> &path);
    void printPath(std::vector<int> &path);

    /// Paths to every node from start and goal
    std::unordered_map<int, std::vector<int>> paths_from_start_;
    std::unordered_map<int, std::vector<int>> paths_from_goal_;
    std::vector<int> ub_path_;

    /// Map that maps state IDs to lower-bound and upper-bound costs from start to goal via them.
    std::unordered_map<int, double> lb_cost_;
    std::unordered_map<int, double> ub_cost_;

    std::vector<std::shared_ptr<InsatAction>> insat_actions_ptrs_;
    InsatStatePtrType start_state_ptr_;
    InsatStatePtrType goal_state_ptr_;
    InsatStateQueueMinType insat_state_open_list_;
    InsatStatePtrMapType insat_state_map_;
    TrajType soln_traj_;

  };

}

#endif
