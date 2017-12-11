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

  double del;//重心の振動の初期位相

  double f(vector<double> &x, double t);//v' = f(v,x,t)
  double k(vector<double> &x, double t);//x' = k(v,x,t) = v

  double r(double t,double w, double del_);
  double drdt(double t,double w, double del_);

public:
  double get_r(double t);
  RungeKutta(double del_,double a_,double keisu_w_a);
  void cal(vector<double> &x, double t,double dt);
  double get_w0(){return w_0;}

};
RungeKutta::RungeKutta(double del_,double a_,double keisu_w_a)
{

  del = del_;
  g = 9.8;//m/t^2
  l = 1.8;//m
  a = a_;//m
  w_0 = sqrt(g/l);
  w_a = keisu_w_a*w_0;
}

double RungeKutta::get_r(double t)
{
  return l - a*cos(w_a*t + del);
}

double RungeKutta::r(double t,double w, double del_)
{
  return l - a*cos(w*t + del_);
}
double RungeKutta::drdt(double t, double w,double del_)
{
  return a*w*sin(w*t + del_);
}

double RungeKutta::f(vector<double> &x, double t)
{
  //x[0] = v, x[1] = x
  return -2*drdt(t,w_a,del)/r(t,w_a,del)*x[0] - g*sin(x[1])/r(t,w_a,del);
}
double RungeKutta::k(vector<double> &x, double t)
{
  return x[0];
}

//t=t_nの変数をベクトルで渡すと、t=t_n+1の値に書き換える
//計算のためにその際の時刻t、時間刻みを渡す。
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
  
int Nt = atoi(argv[1]);
double del = atof(argv[2]);
double keisu_w_a = atof(argv[3]);
double a = atof(argv[4]);
double t_fin
RungeKutta run(del,a,keisu_w_a);
double w_0 = run.get_w0();
double T = 2.*M_PI/w_0;
double t_0 = 0;
double t_fin = 6.*T;
double dt = (t_fin - t_0)/Nt;

//初期値
double v_0 = 0.;
double x_0 = 0.3;//初期値は、0.3 rad * 1.8 m ~ 50 cm くらいのところ

vector<double> x(2);
x[0] = v_0;
x[1] = x_0;
ofstream out("RK.dat");

out << t_0 << " " << x[0] << " " << x[1] << " " << x_0*cos(w_0 * t_0) << " " << run.get_r(t_0) << endl;
for(int i=1; i<=Nt; i++)
{
double t = t_0 + dt*i;
run.cal(x,t,dt);
out << t << " " << x[0] << " " << x[1] << " " << x_0*cos(w_0 * t) << " " << run.get_r(t)  << endl;
}
out.close();

//////////////////////////////
//gnuplotへの命令、図の表示
//////////////////////////////
stringstream ss_Nt,ss_del,ss_a,ss_keisu_omega_a;
ss_Nt << Nt; // 文字列と数値を混ぜて書ける
string Nt_st = ss_Nt.str(); // .str()でstring型の値を得る

ss_del << del;
string del_st = ss_del.str();

ss_a << a;
string a_st = ss_a.str();

ss_keisu_omega_a << keisu_w_a;
string keisu_omega_a_st = ss_keisu_omega_a.str();

stringstream ss_x0;
ss_x0 << x_0;
string x0_st = ss_x0.str();

stringstream ss_t_fin;
ss_t_fin << t_fin;
string t_fin_st = ss_t_fin.str();


string FILE_NAME = "Nt" + Nt_st  + "del" + del_st + "a" + a_st + "keisu" + keisu_omega_a_st + "x0" + x0_st + "t_fin" + t_fin_st +  "_RK_jusin.eps";
string OPEN_FILE_NAME = "open " + FILE_NAME;
string command = "gnuplot -e \"filename = \'" + FILE_NAME +
"\';del_ = \'" + del_st + "\';keisu_omega_a_ = \'"
 + keisu_omega_a_st + "\';a_ = \'" + a_st + "\';"
+ "\" plotRK_jusin.gp" ;


const char *open_gnu = NULL;
open_gnu = command.c_str();
//    cout << open_gnu;
system(open_gnu);

const char *open_eps = NULL;
open_eps = OPEN_FILE_NAME.c_str();
system(open_eps);

return 0;

}
