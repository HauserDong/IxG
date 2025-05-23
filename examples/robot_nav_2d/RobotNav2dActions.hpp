#ifndef POINT_ROBOT_ACTION_HPP
#define POINT_ROBOT_ACTION_HPP

#include <common/Action.hpp>

namespace ps
{

class RobotNav2dAction : public Action
{

public:
    RobotNav2dAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true)
    : Action(type, params, is_expensive), map_(map) {};
    virtual bool CheckPreconditions(const StateVarsType& state, int thread_id=0); 
    ActionSuccessor GetSuccessor(const StateVarsType& state_vars, int thread_id); 
    ActionSuccessor GetSuccessorLazy(const StateVarsType& state_vars, int thread_id); 
    ActionSuccessor Evaluate(const StateVarsType& parent_state_vars, const StateVarsType& child_state_vars, int thread_id=0);
     
protected:
    bool isValidCell(int x, int y);
    bool inRange(int x, int y);
    std::vector<std::pair<int, int>> getFootPrintRectangular(int x, int y, int footprint_size);
    std::vector<double> move_dir_;
    std::vector<std::vector<int>> map_;
    std::vector<std::pair<int, int>> footprint_;
    LockType lock_;
};

class MoveUpAction : public RobotNav2dAction
{

public:
    MoveUpAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {0, 1};
    };
};

class MoveUpRightAction : public RobotNav2dAction
{

public:
    MoveUpRightAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {1, 1};
    };
};

class MoveRightAction : public RobotNav2dAction
{

public:
    MoveRightAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {1, 0};
    };
};

class MoveRightDownAction : public RobotNav2dAction
{

public:
    MoveRightDownAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {1, -1};
    };
};


class MoveDownAction : public RobotNav2dAction
{

public:
    MoveDownAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {0, -1};
    };
};

class MoveDownLeftAction : public RobotNav2dAction
{

public:
    MoveDownLeftAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {-1, -1};
    };
};


class MoveLeftAction : public RobotNav2dAction
{

public:
    MoveLeftAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {-1, 0};
    };
};


class MoveLeftUpAction : public RobotNav2dAction
{

public:
    MoveLeftUpAction(const std::string& type, ParamsType params, std::vector<std::vector<int>> map, bool is_expensive = true):
    RobotNav2dAction(type, params, map, is_expensive)
    {
        move_dir_ = {-1, 1};
    };
};

}

#endif
