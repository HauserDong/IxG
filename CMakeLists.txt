cmake_minimum_required(VERSION 3.2)

project (ixg)
add_compile_options(-std=c++17)

list(APPEND CMAKE_PREFIX_PATH /opt/drake)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(drake REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem serialization)
find_package(PkgConfig REQUIRED)
find_package(OpenCV REQUIRED)

pkg_check_modules(SBPL REQUIRED sbpl)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-std=c++17)

set(CMAKE_CXX_FLAGS "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_C_FLAGS_RELEASE "-O3")

include_directories(include
                    ${EIGEN3_INCLUDE_DIR}
                    ${drake_INCLUDE_DIRS})
include_directories(/opt/drake/include)
include_directories(${SBPL_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})

link_directories(${SBPL_LIBRARY_DIRS})
#link_directories(/opt/drake/lib)

set(INSATxGCS_SOURCES
        src/common/State.cpp
        src/common/Edge.cpp
        src/common/insat/InsatEdge.cpp
        src/common/insatxgcs/utils.cpp
        src/common/insatxgcs/gcsbfs.cpp
        src/planners/Planner.cpp
        src/planners/GepasePlanner.cpp
        src/planners/insat/INSATxGCS.cpp
        src/planners/insat/pINSATxGCS.cpp
        src/planners/insat/opt/GCSOpt.cpp
        src/planners/insat/opt/GCSSmoothOpt.cpp
        src/planners/insat/opt/LBGraph.cpp
)

set(INSAT_SOURCES
        src/common/State.cpp
        src/common/Edge.cpp
        src/common/insat/InsatEdge.cpp
        src/common/insatxgcs/utils.cpp
        src/common/insatxgcs/gcsbfs.cpp
        src/planners/Planner.cpp
        src/planners/GepasePlanner.cpp
        src/planners/WastarPlanner.cpp
        src/planners/PwastarPlanner.cpp
        src/planners/EpasePlanner.cpp
        src/planners/PasePlanner.cpp
        src/planners/MplpPlanner.cpp
        src/planners/insat/INSATxGCS.cpp
        src/planners/insat/pINSATxGCS.cpp
        src/planners/insat/opt/GCSOpt.cpp
        src/planners/insat/opt/GCSSmoothOpt.cpp
        src/planners/insat/opt/LBGraph.cpp

)

# EXECUTABLES
add_executable(run_insatxgcs
        examples/insatxgcs/run_insatxgcs.cpp
        examples/insatxgcs/INSATxGCSAction.cpp
        ${INSATxGCS_SOURCES})

target_link_libraries(run_insatxgcs
        ${drake_LIBRARIES}
        pthread)


add_executable(gcsopt_test
        examples/insatxgcs/gcsopt_test.cpp
        src/common/insatxgcs/utils.cpp
        src/planners/insat/opt/GCSOpt.cpp)

target_link_libraries(gcsopt_test
        ${drake_LIBRARIES}
        pthread)

add_executable(bfs_test
        examples/insatxgcs/bfs_test.cpp
        src/common/insatxgcs/utils.cpp
        src/common/insatxgcs/gcsbfs.cpp)

target_link_libraries(bfs_test
        ${drake_LIBRARIES}
        pthread)

add_executable(trigcs_monotonicity
        examples/insatxgcs/trigcs_monotonicity.cpp
        src/common/insatxgcs/utils.cpp
        src/planners/insat/opt/GCSOpt.cpp)

target_link_libraries(trigcs_monotonicity
        ${drake_LIBRARIES}
        pthread)

add_executable(gcsopt_monotonicity
        examples/insatxgcs/gcsopt_monotonicity.cpp
        src/common/insatxgcs/utils.cpp
        src/planners/insat/opt/GCSOpt.cpp)

target_link_libraries(gcsopt_monotonicity
        ${drake_LIBRARIES}
        pthread)

add_executable(read_itamar_logs
        examples/insatxgcs/read_itamar_logs.cpp)

target_link_libraries(read_itamar_logs
        ${Boost_LIBRARIES})

add_executable(lbg_test
        examples/insatxgcs/lbg_test.cpp
        src/common/insatxgcs/utils.cpp
        src/planners/insat/opt/GCSOpt.cpp
        src/planners/insat/opt/GCSSmoothOpt.cpp
        src/planners/insat/opt/LBGraph.cpp)

target_link_libraries(lbg_test
        ${drake_LIBRARIES}
        pthread)

add_executable(run_robot_nav_2d
        examples/robot_nav_2d/run_robot_nav_2d.cpp
        examples/robot_nav_2d/RobotNav2dActions.cpp
        ${INSAT_SOURCES})

target_link_libraries(run_robot_nav_2d
        ${drake_LIBRARIES}
        ${SBPL_LIBRARIES}
        ${OpenCV_LIBS}
        pthread)