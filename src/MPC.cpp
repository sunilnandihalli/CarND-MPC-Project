#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

size_t N = 30;
static double dt = 0.05;
const double Lf = 2.67;
double ref_v = 10;
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;
size_t pd_start = a_start + N - 1; // positive direction constraint

class FG_eval {
 public:
  Eigen::VectorXd coeffs;   
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    
    for (int t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2);
      fg[0] += CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }
    for (int t = 0; t < N - 1; t++) {
      // Minimize the use of actuators.
      fg[0] += CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }
    for (int t = 0; t < N - 2; t++) {  // reduce change in actuation neighbouring timesteps
      fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2)*0.1;
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
   
    for (int t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];
      
      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      
      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      auto& c(coeffs);
      //AD<double> f1 = ((c[3]*x1+c[2])*x1+c[1])*x1+c[0] ;
      AD<double> f1 = 0;
      for(int i=c.size()-1;i>-1;i--)
	f1 = f1*x1+c[i];
      //AD<double> fdash1 = (3*c[3]*x1+2*c[2])*x1+c[1];
      AD<double> fdash1 = 0;
      for(int i=c.size()-1;i>0;i--)
	fdash1 = fdash1*x1+c[i]*i;
      AD<double> norm = CppAD::sqrt(1+fdash1*fdash1);
      
      AD<double> x1p,y1p;
      AD<double> tol=1e-9;
      AD<double> dist = v0*dt + 0.5*a0*dt*dt;
      if(CppAD::abs_geq(delta0,tol)) {
	AD<double> r = Lf/delta0;
	x1p = x0 - r*(CppAD::sin(psi1)-CppAD::sin(psi0));
	y1p = y0 - r*(CppAD::cos(psi0)-CppAD::cos(psi1));
      } else {
	x1p = x0 + dist*CppAD::cos(psi0);
	y1p = y0 + dist*CppAD::sin(psi0);
      }
      AD<double> psi1p = psi0 - dist*(delta0 / Lf );
      AD<double> v1p = v0 + a0 * dt;
      AD<double> cte1p = (f1-y1)/norm;
      AD<double> epsi1p = (fdash1*CppAD::cos(psi1)-CppAD::sin(psi1))/norm;
      fg[1 + x_start + t] = x1 - x1p;
      fg[1 + y_start + t] = y1 - y1p;
      fg[1 + psi_start + t] = psi1 - psi1p;
      fg[1 + v_start + t] = v1 - v1p;
      fg[1 + cte_start + t] = cte1 - cte1p;
      fg[1 + epsi_start + t] = epsi1 - epsi1p;
    }
  }
};

MPC::MPC() {}
MPC::~MPC() {}

tuple<vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,vector<double>,bool>
MPC::Solve(Eigen::VectorXd x0, double deltaT,Eigen::VectorXd coeffs) {
  dt = deltaT;
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  double x(x0[0]),y(x0[1]),psi(x0[2]),v(x0[3]),cte(x0[4]),epsi(x0[5]);
  size_t n_vars(x0.size()*N + 2*(N-1)),n_constraints(x0.size()*N);
  
  Dvector vars(n_vars),vars_lowerbound(n_vars),vars_upperbound(n_vars),constraints_lowerbound(n_constraints),constraints_upperbound(n_constraints);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332; // 25 deg in radians
  }
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0; // throttle
  }

  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  FG_eval fg_eval(coeffs);// object that computes objective and constraints

  std::string options; // options for IPOPT solver
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  //options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          250\n";
  CppAD::ipopt::solve_result<Dvector> solution;
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  {
    vector<double> x(N),y(N),psi(N),v(N),cte(N),epsi(N),delta(N-1),a(N-1);
    for(int i=0;i<N;i++) {
      x[i]=solution.x[x_start+i];
      y[i] = solution.x[y_start+i];
      psi[i] = solution.x[psi_start+i];
      v[i] = solution.x[v_start+i];
      cte[i] = solution.x[cte_start+i];
      epsi[i] = solution.x[epsi_start+i];
      if(i<N-1) {
	delta[i] = solution.x[delta_start+i];
	a[i] = solution.x[a_start+i];
      }
    }
    return make_tuple(x,y,psi,v,cte,epsi,delta,a,ok);
  }
}
