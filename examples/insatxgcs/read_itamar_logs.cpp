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
 * \file read_itamar_logs.cpp 
 * \author Ram Natarajan (rnataraj@cs.cmu.edu)
 * \date 11/8/23
*/

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

struct Region
{
  friend class boost::serialization::access;
  std::vector<double> start;
//    unsigned int radius;
  double radius;
  std::vector<double> state; // workspace state
  std::vector<std::vector<double> > path;
  std::vector<std::vector<double> > path2;   // TODO: Do I want this here?
  // Goal to goal paths
  std::vector<std::vector<std::vector<double> >> paths;

  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & start;
    ar & radius;
    ar & state;
    ar & path;
    ar & path2;
    ar & paths;
  }
};

void writeVectorOfVectorToFile(const std::vector<Region>& regions, const std::string& filename) {
  std::ofstream file(filename); // Open the file for writing

  if (file.is_open()) {
    for (auto& r : regions) {
      auto data = r.path;

      for (const std::vector<double>& innerVector : data) {
        for (double value : innerVector) {
          file << value << " "; // Write the double values to the file
        }
        file << "\n"; // Move to the next line for the next inner vector
      }
      file << "-1 -1 -1 -1 -1 -1" << "\n"; // Write the double values to the file
    }

    file.close(); // Close the file
    std::cout << "Data has been written to the file: " << filename << std::endl;
  } else {
    std::cerr << "Failed to open the file for writing: " << filename << std::endl;
  }
}

void writeStartsAndGoalsToFile(const std::vector<Region>& regions,
                               const std::string& start_file,
                               const std::string& goal_file) {
  std::ofstream sfile(start_file); // Open the file for writing
  std::ofstream gfile(goal_file); // Open the file for writing

  if (sfile.is_open()) {
    for (auto& r : regions) {
      auto data = r.path;
      auto start = data.front();
      auto goal = data.back();
      for (double value : start) {
        sfile << value << " "; // Write the double values to the file
      }
      sfile << "\n"; // Move to the next line for the next inner vector

      for (double value : goal) {
        gfile << value << " "; // Write the double values to the file
      }
      gfile << "\n"; // Move to the next line for the next inner vector
    }

    sfile.close(); // Close the file
    gfile.close(); // Close the file
    std::cout << "Data has been written to the file: " << start_file << std::endl;
    std::cout << "Data has been written to the file: " << goal_file << std::endl;
  } else {
    std::cerr << "Failed to open the file for writing: " << start_file << std::endl;
    std::cerr << "Failed to open the file for writing: " << goal_file << std::endl;
  }
}

std::vector<Region> ReadRegions(std::string path) {
  std::vector<Region> regions;

  try {
    boost::filesystem::path myFile = path; // boost::filesystem::current_path() /
    boost::filesystem::ifstream ifs(myFile);
    boost::archive::text_iarchive ta(ifs);
    ta >> regions;
  }
  catch (...) {
    std::runtime_error("Unable to read preprocessed file");
  }

  return regions;
}

int main() {

  std::string path = "/home/gaussian/cmu_ri_phd/phd_research/ixg/logs/itamar/";
  std::string filename = "manipulator_3_pick";
  std::string filepath = path + filename + ".dat";

  std::vector<Region> regions = ReadRegions(filepath);

  std::string savepath = path+filename+".txt";
  writeVectorOfVectorToFile(regions, savepath);

  std::string startpath = path+filename+"_starts.txt";
  std::string goalpath = path+filename+"_goals.txt";
  writeStartsAndGoalsToFile(regions, startpath, goalpath);

  return 0;
}
