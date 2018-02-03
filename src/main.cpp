#define HAVE_CSTDDEF
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
const double Lf = 2.67;
const double tol = 1e-9;
// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd& coeffs, double x) {
  auto& c(coeffs);
  double result = 0.0;
  for(int i=c.size()-1;i>-1;i--)
    result = result*x+c[i];
  return result;
}

vector<double> polyeval(Eigen::VectorXd& coeffs,vector<double> xs) {
  vector<double> ys;
  for(double x:xs) ys.push_back(polyeval(coeffs,x));
  return ys;
}

double polyDerivativeEval(Eigen::VectorXd coeffs, double x) {
  auto& c(coeffs);
  double result = 0;
  for (int i=c.size()-1;i>0;i--) {
    result = result*x + c[i]*i;
  }
  return result;
}

Eigen::VectorXd polyfit(vector<double> xvals,vector<double> yvals, int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  cout<<"polyfit : ";
  for(int i=0;i<xvals.size();i++) {
    cout<<"("<<xvals[i]<<","<<yvals[i]<<") ";
  }
  cout<<endl;
  Eigen::MatrixXd A(xvals.size(), order + 1);
  Eigen::VectorXd Y(xvals.size());
  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
    Y(i) = yvals[i];
  }
  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals[j];
    }
  }
  auto Q = A.householderQr();
  auto result = Q.solve(Y);
  return result;
}

void globalKinematic(double x0,double y0,double psi0,double v0,double delta0,double a0,double dt,
		     double& x1,double& y1,double& psi1,double& v1) {
  v1 = v0 + a0*dt;
  double dist = v0*dt+0.5*a0*dt*dt;
  psi1 = psi0 - (delta0/Lf)*dist;
  if(fabs(delta0)>tol) {
    double r  =  Lf/delta0;
    x1 = x0 - r*(sin(psi1)-sin(psi0));
    y1 = y0 - r*(cos(psi0)-cos(psi1));
  } else {
    double psi_avg = 0.5*(psi0+psi1);
    x1 = x0 + dist*cos(psi_avg);
    y1 = y0 + dist*sin(psi_avg);
  }
}

tuple<vector<double>,vector<double>> toCarCoords(double px,double py,double psi,vector<double> xs,vector<double> ys) {
  cout<<" toCarCoords inp-sizes : "<<xs.size()<<" "<<ys.size()<<endl;
  double ct = cos(psi);
  double st = sin(psi);
  auto f = [&px,&py,&ct,&st](double x,double y)->tuple<double,double> {
    double dx = x-px;
    double dy = y-py;
    double x1= dx*ct+dy*st;
    double y1= -dx*st+dy*ct;
    return std::make_tuple(x1,y1);
  };
  vector<double> rxs(xs.size()),rys(xs.size());
  for(int i =0;i<xs.size();i++) {
    std::tie(rxs[i],rys[i]) = f(xs[i],ys[i]);
  }
  return std::make_tuple(rxs,rys);
}

template<class T>
void printRow(vector<T> vs) {
  return;
  for(auto& v:vs) 
    cout<<setw(15)<<v;
  cout<<endl;
}

template<class T>
void dumpRow(vector<T> vs,std::ofstream& fout) {
  std::string comma(" ");
  for(int i=0;i<vs.size();i++) {
    fout<<vs[i];
    if(i<vs.size()-1)
      fout<<comma;
  }
  fout<<endl;
}

tuple<double,double> cte_and_epsi(double x,double y,double psi,Eigen::VectorXd C) {
  double fx = polyeval(C,x);
  double fdash = polyDerivativeEval(C,x);
  double norm = sqrt(1+fdash*fdash);
  double cte = (fx-y)/norm;
  double epsi = (fdash*cos(psi)-sin(psi))/norm;
  return make_tuple(cte,epsi);
}

typedef tuple<vector<double>,vector<double>,vector<double>,vector<double>,
	      vector<double>,vector<double>,vector<double>,vector<double>,bool > solret;

void printSolution(solret& sol,Eigen::VectorXd& C,double px,double py,double psi,vector<double> ptsx,vector<double> ptsy) {
  static int id=0;
  id++;
  cout<<"---------------------------------start "<<id<<"--------------------------------------------------"<<endl;
  double dt = 0.05;
  vector<double> xs,ys,psis,vs,ctes,epsis,deltas,as,cxs,cys,fxs;
  bool ok;
  std::tie(xs,ys,psis,vs,ctes,epsis,deltas,as,ok) = sol;
  std::tie(cxs,cys) = toCarCoords(px,py,psi,xs,ys);
  fxs = polyeval(C,xs);
  int N = xs.size();
  {
    std::stringstream ss;
    ss<<"solution_"<<setw(5)<<setfill('0')<<id;
    std::ofstream fout(ss.str());
    vector<double> actes(N),aepsis(N);
    for(int i=0;i<N;i++) {
      std::tie(actes[i],aepsis[i]) = cte_and_epsi(xs[i],ys[i],psis[i],C);
    }
    printRow<string>({"xs","ax","fxs","cxs","ys","ay","cys","psis","apsi","vs","av","ctes","actes","epsis","aepsis","deltas","as"});
    dumpRow<string>({"xs","ax","fxs","cxs","ys","ay","cys","psis","apsi","vs","av","ctes","actes","epsis","aepsis","deltas","as"},fout);
    for(int i=0;i<N-1;i++) {
      double psisi = psis[i]*180.0/M_PI;
      double acte,aepsi;
      double ax,ay,apsi,av;
      if(i>0) 
	globalKinematic(xs[i-1],ys[i-1],psis[i-1],vs[i-1],deltas[i-1],as[i-1],dt,ax,ay,apsi,av);
      apsi *= 180.0/M_PI;
      printRow<double>({xs[i],ax,fxs[i],cxs[i],ys[i],ay,cys[i],psisi,apsi,vs[i],av,ctes[i],actes[i],epsis[i],aepsis[i],deltas[i],as[i]});
      dumpRow<double>({xs[i],ax,fxs[i],cxs[i],ys[i],ay,cys[i],psis[i],apsi,vs[i],av,ctes[i],actes[i],epsis[i],aepsis[i],deltas[i],as[i]},fout);
    }
    double psisi = psis[N-1]*180.0/M_PI;
    double ax,ay,apsi,av;
    globalKinematic(xs[N-2],ys[N-2],psis[N-2],vs[N-2],deltas[N-2],as[N-2],dt,ax,ay,apsi,av);
    apsi *= 180.0/M_PI;
    printRow<double>({xs[N-1],ax,fxs[N-1],cxs[N-1],ys[N-1],ay,cys[N-1],psisi,apsi,vs[N-1],av,ctes[N-1],actes[N-1],epsis[N-1],aepsis[N-1]});
    dumpRow<double>({xs[N-1],ax,fxs[N-1],cxs[N-1],ys[N-1],ay,cys[N-1],psis[N-1],apsi,vs[N-1],av,ctes[N-1],actes[N-1],epsis[N-1],aepsis[N-1],0,0},fout);
  }
  {
    std::stringstream ss;
    ss<<"way_"<<setw(5)<<setfill('0')<<id;
    std::ofstream fout(ss.str());
    vector<double> cptsx,cptsy;
    std::tie(cptsx,cptsy) = toCarCoords(px,py,psi,ptsx,ptsy);
    fout<<"x y"<<endl;
    for(int i=0;i<ptsx.size();i++) {
      fout<<ptsx[i]<<' '<<ptsy[i]<<endl;
    }
  }
  cout<<"##################################### end "<<id<<" ##############################################"<<endl;
  if(!ok)
    cout<<"solver failed";
  for(int i=1;i<N;i++)
    if(!(cxs[i]>cxs[i-1])) cout<<" order inconsistent "<<endl;
}

double samedir(double x0,double y0,double x1,double y1,double psi) {
  double dx(x1-x0),dy(y1-y0);
  double r = sqrt(dx*dx+dy*dy);
  return (dx*cos(psi)+dy*sin(psi))/r;
}

void calculateActivation(double x0,double y0,double psi0,double v0,
			 double delta0,double a0,vector<double> cptsx0,vector<double> cptsy0,
			 double& delta1,double& a1,vector<double>& xs_pred,vector<double>& ys_pred,double deltaT) {
  static MPC mpc;
  
  Eigen::VectorXd coeffs(polyfit(cptsx0,cptsy0,3)),state(6);
  cout<<" coeffs : "<<coeffs.transpose()<<endl;
  double cte0,epsi0;
  cout<<"state : "<<x0<<" "<<y0<<" "<<psi0<<" "<<v0<<" "<<delta0<<" "<<a0<<endl;
  //assert(samedir(cptsx0[0],cptsy0[0],cptsx0[1],cptsy0[1],0)>0);
  std::tie(cte0,epsi0) = cte_and_epsi(0,0,0,coeffs);
  state<<0,0,0,v0,cte0,epsi0;
  auto controls = mpc.Solve(state,deltaT,coeffs);
  printSolution(controls,coeffs,0,0,0,cptsx0,cptsy0);
  delta1 = std::get<6>(controls)[0];
  a1 = std::get<7>(controls)[0];
  xs_pred = std::get<0>(controls);
  ys_pred = std::get<1>(controls);
}

string cat(std::vector<string> x) {
  stringstream ss;
  for(auto s:x) {
    ss<<s;
  }
  return ss.str();
}
string path(string s) {
  return cat({"/home/sunil/carnd/CarND-MPC-Project/",s});
}

double seconds(std::chrono::system_clock::time_point start,std::chrono::system_clock::time_point end) {
  return std::chrono::duration<double>(end-start).count();
}

int main() {
  uWS::Hub h;
  static int id = 0;
  ofstream fout(path("debug/history.txt"));
  fout<<"type,time,lag,x,y,psi,v,delta,a,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5"<<endl;
  auto time = std::chrono::system_clock::now();
  auto start = time;
  h.onMessage([&fout,&time,&start](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    string sdata = string(data).substr(0, length);
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
	  auto now = std::chrono::system_clock::now();
	  double sim_round_trip = seconds(time,now);
	  time = now;
          vector<double> ptsx=j[1]["ptsx"],ptsy=j[1]["ptsy"];
	  double px(j[1]["x"]),py(j[1]["y"]),psi(j[1]["psi"]),
	    v(j[1]["speed"]),psi_unity(j[1]["psi_unity"]),
	    delta(j[1]["steering_angle"]),throttle(j[1]["throttle"]);
	  fout<<"recv,"<<seconds(start,now)<<","<<sim_round_trip<<","<<px<<","<<py<<","<<psi<<","<<v<<","<<delta<<","<<throttle;
	  for(auto x:ptsx) fout<<","<<x; for(auto y:ptsy) fout<<","<<y; fout<<endl; fout.flush();
	  static double delay_estimate = 0.1; // delay of 100 MilliSeconds
          static double steer_value = 0.0;
          static double throttle_value = 0.0;
	  globalKinematic(px,py,psi,v,steer_value,throttle,delay_estimate,px,py,psi,v);
	  vector<double> mpc_x_vals,mpc_y_vals,next_x_vals,next_y_vals;
	  std::tie(next_x_vals,next_y_vals) = toCarCoords(px,py,psi,ptsx,ptsy);
	  calculateActivation(px,py,psi,v,steer_value,throttle_value,next_x_vals,next_y_vals,
			      steer_value,throttle_value,mpc_x_vals,mpc_y_vals,delay_estimate);
	  this_thread::sleep_for(chrono::milliseconds(100)); // simulate latency
	  now = std::chrono::system_clock::now();
	  delay_estimate = seconds(time,now);
	  time = now;
	  fout<<"pred,"<<seconds(start,now)<<","<<delay_estimate<<","<<px<<","<<py<<","<<psi<<","<<v<<","<<steer_value<<","<<throttle_value;
	  for(auto x:next_x_vals) fout<<","<<x; for(auto y:next_y_vals) fout<<","<<y; fout<<endl; fout.flush();

	  
	  json msgJson;
	  // rescale steer_value which is in radians to -1 to 1
          msgJson["steering_angle"] = steer_value/0.436332;
	  msgJson["throttle"] = throttle_value;
	  msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;
	  msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;
	  
          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
