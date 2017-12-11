#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
/*
運動方程式
x'' = f(x', x, t)
↓x' := v
v' = f(v,x,t)
x' = k(v,x,t) = v
を解く。
*/


class RungeKutta
{
private:
  double g,l,a;//重力加速度、ブランコの紐の長さ、人間の重心の移動範囲
  double w_0;//単振り子の場合の角振動数
  double w_a;//重心の角振動数
  double f(vector<double> &x, double t);
  double k(vector<double> &x, double t);

public:
  RungeKutta();
  void cal(vector<double> &x, double t,double dt);
  double get_w0(){return w_0;}

};
RungeKutta::RungeKutta()
{

  g = 9.8;//m/t^2
  l = 1.8;//m
  a = 0.3;//m
  w_0 = sqrt(g/l);
  w_a = 2.*w_0;
}

double RungeKutta::f(vector<double> &x, double t)
{
  //x[0] = v, x[1] = x
  return -g/l * x[1];
}
double RungeKutta::k(vector<double> &x, double t)
{
  return x[0];
}

void RungeKutta::cal(vector<double> &x, double t, double dt)
{
  int size = x.size();
  vector<double> x0(size,0.), x1(size,0.), x2(size,0.), x3(size,0.);

  x0 = x;

  x1[0] = x0[0] + dt/2. * f(x0,t);
  x1[1] = x0[1] + dt/2. * k(x0,t);

  x2[0] = x0[0] + dt/2. * f(x1,t+dt/2.);
  x2[1] = x0[1] + dt/2. * k(x1,t+dt/2.);

  x3[0] = x0[0] + dt * f(x2,t+dt/2.);
  x3[1] = x0[1] + dt * k(x2,t+dt/2.);

  x[0] = x0[0] + dt/6. * (f(x0,t) + 2*f(x1,t+dt/2.) + 2*f(x2,t+dt/2.)
                          + f(x3,t+dt));
  x[1] = x0[1] + dt/6. * (k(x0,t) + 2*k(x1,t+dt/2.) + 2*k(x2,t+dt/2.)
                          + k(x3,t+dt));
}






int main(int argc, char *argv[])
{
  int np=1;
  if (argc!=np+1) {
    cout << "Enter (Nt)" << endl;
    exit (0);
  }

int Nt = atoi(argv[1]);
RungeKutta run;
double w_0 = run.get_w0();
double T = 2*M_PI/w_0;
double t_0 = 0;
double t_fin = 4.*T;
double dt = (t_fin - t_0)/Nt;

//初期値
double v_0 = 0.;
double x_0 = 0.3;

vector<double> x(2);
x[0] = v_0;
x[1] = x_0;
ofstream out("RK.dat");

out << t_0 << " " << x[0] << " " << x[1] << " " << x_0*cos(w_0 * t_0) << endl;
for(int i=1; i<=Nt; i++)
{
  double t = t_0 + dt*i;
  run.cal(x,t,dt);
  out << t << " " << x[0] << " " << x[1] << " " << x_0*cos(w_0 * t) << endl;
}
out.close();

//////////////////////////////
//gnuplotの操作
//////////////////////////////
stringstream ss_Nt;
ss_Nt << Nt; // 文字列と数値を混ぜて書ける
string Nt_st = ss_Nt.str(); // .str()でstring型の値を得る


string FILE_NAME = "Nt" + Nt_st  + "_RK.eps";
string OPEN_FILE_NAME = "open ./eps/" + FILE_NAME;
string command = "gnuplot -e \"filename = \'" + FILE_NAME + "\" plotRK.gp" ;


const char *open_gnu = NULL;
open_gnu = command.c_str();
//    cout << open_gnu;
system(open_gnu);

const char *open_eps = NULL;
open_eps = OPEN_FILE_NAME.c_str();
system(open_eps);






return 0;

}
