#include <common/insatxgcs/utils.hpp>

namespace utils {

/*
This code is copied and modifed from scipy.optimize.root_scalar's source code for brentq. 
origin: https://github.com/scipy/scipy/blob/main/scipy/optimize/Zeros/brentq.c
*/
  double Brentq(const std::function<double(double)>& func, double xa, double xb, double xtol, double rtol, int iter) {
    double xpre = xa, xcur = xb;
    double xblk = 0.0, fpre, fcur, fblk = 0.0, spre = 0.0, scur = 0.0, sbis;
    double delta;
    double stry, dpre, dblk;
    int i;
    fpre = func(xpre);
    fcur = func(xcur);
    if (fpre == 0) return xpre;
    if (fcur == 0) return xcur;
    if (signbit(fpre) == signbit(fcur))
      throw std::invalid_argument("Root is not bracketed, cannot proceed");

    for (i = 0; i < iter; i++) {
      if (fpre != 0 && fcur != 0 && (signbit(fpre) != signbit(fcur))) {
        xblk = xpre;
        fblk = fpre;
        spre = scur = xcur - xpre;
      }
      if (fabs(fblk) < fabs(fcur)) {
        xpre = xcur;
        xcur = xblk;
        xblk = xpre;

        fpre = fcur;
        fcur = fblk;
        fblk = fpre;
      }

      delta = (xtol + rtol*fabs(xcur))/2;
      sbis = (xblk - xcur)/2;
      if (fcur == 0 || fabs(sbis) < delta) return xcur;

      if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
        if (xpre == xblk) {
          stry = -fcur*(xcur - xpre)/(fcur - fpre);
        }
        else {
          dpre = (fpre - fcur)/(xpre - xcur);
          dblk = (fblk - fcur)/(xblk - xcur);
          stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre));
        }
        if (2*fabs(stry) < std::min(fabs(spre), 3*fabs(sbis) - delta)) {
          spre = scur;
          scur = stry;
        } else {
          spre = sbis;
          scur = sbis;
        }
      }
      else {
        spre = sbis;
        scur = sbis;
      }

      xpre = xcur; fpre = fcur;
      if (fabs(scur) > delta) xcur += scur;
      else xcur += (sbis > 0 ? delta : -delta);

      fcur = func(xcur);
    }
    return xcur;
  }

  std::shared_ptr<std::vector<std::pair<int,int>>> DeserializeEdges(const std::string& file_path) {
    std::vector<std::pair<int,int>> edges;

    std::ifstream file(file_path);
    if (!file.is_open())
      throw std::runtime_error("Cannot open file: " + file_path);

    size_t size; int u, v; char comma;
    file >> size;
    for (size_t idx = 0; idx < size; ++idx) {
      file >> comma >> u >> comma >> v;
      edges.push_back(std::make_pair(u, v));
    }
    return std::make_shared<std::vector<std::pair<int,int>>>(edges);
  }

  std::vector<HPolyhedron> DeserializeRegions(const std::string& file_path) {
    std::vector<HPolyhedron> regions;

    std::ifstream file(file_path);
    if (!file.is_open())
      throw std::runtime_error("Cannot open file: " + file_path);

    size_t size, rows, cols; double value; char comma;
    file >> size;

    for (size_t idx = 0; idx < size; ++idx) {
      file >> comma >> rows >> comma >> cols;
      Eigen::MatrixXd A(rows, cols);
      Eigen::VectorXd b(rows);

      for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
          file >> comma >> A(i, j);

      for (size_t i = 0; i < rows; ++i)
        file >> comma >> b(i);

      regions.push_back(HPolyhedron(A, b));
    }

    return regions;
  }

  void SerializeRegions(const std::vector<HPolyhedron>& regions, const std::string& file_path) {
    std::ofstream file(file_path);
    if (!file.is_open())
      throw std::runtime_error("Cannot open file: " + file_path);

    file << regions.size();
    for (auto & region : regions) {
      file << "," << region.A().rows() << "," << region.A().cols();
      for (size_t i = 0; i < region.A().rows(); ++i) {
        for (int j=0; j<region.A().row(i).size(); ++j) {
          file << "," << region.A().row(i)(j);
        }
      }
//        for (auto & entry : region.A().row(i))
//          file << "," << entry;
      for (int i=0; i<region.b().size(); ++i) {
        file << "," << region.b()(i);
      }
//      for (auto & entry : region.b())
//        file << "," << entry;
    }
  }

  std::tuple<std::vector<Eigen::MatrixXd>, VectorXb> CondenseToNonzeroColumns(std::vector<Eigen::MatrixXd> matrices) {
    // Validate inputs.
    DRAKE_DEMAND(matrices.size() > 0);
    const int num_cols = matrices[0].cols();
    for (const Eigen::MatrixXd& matrix : matrices) {
      DRAKE_DEMAND(matrix.cols() == num_cols);
    }

    // Find non-zero columns.
    VectorXb nonzero_cols_mask = VectorXb::Constant(num_cols, false);
    for (const Eigen::MatrixXd& matrix : matrices) {
      nonzero_cols_mask += matrix.cast<bool>().colwise().any();
    }
    const int nonzero_cols_count = nonzero_cols_mask.count();

    // Create the output, copying only the non-zero columns.
    std::vector<Eigen::MatrixXd> condensed_matrices;
    for (const Eigen::MatrixXd& matrix : matrices) {
      Eigen::MatrixXd& condensed_matrix =
          condensed_matrices.emplace_back(matrix.rows(), nonzero_cols_count);
      int condensed_col = 0;
      for (int orig_col = 0; orig_col < matrix.cols(); ++orig_col) {
        if (nonzero_cols_mask(orig_col)) {
          condensed_matrix.col(condensed_col) = matrix.col(orig_col);
          condensed_col++;
        }
      }
    }
    return std::make_tuple(condensed_matrices, nonzero_cols_mask);
  }

  drake::VectorX<drake::symbolic::Variable>
  FilterVariables(const drake::VectorX<drake::symbolic::Variable> &vars, const VectorXb &nonzero_cols_mask) {
    drake::VectorX<drake::symbolic::Variable> vars_dense(nonzero_cols_mask.count());
    int row = 0;
    for (int i = 0; i < vars.size(); ++i) {
      if (nonzero_cols_mask(i)) {
        vars_dense(row++) = vars(i);
      }
    }
    return vars_dense;
  }

} // namespace utils