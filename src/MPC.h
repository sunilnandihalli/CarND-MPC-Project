#ifndef MPC_H
#define MPC_H
#define HAVE_CSTDDEF


#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  //vector<double>
  tuple<vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,bool>
    Solve(Eigen::VectorXd state, double dt, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
