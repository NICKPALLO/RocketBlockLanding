#pragma once
#include "matrix.h"
#include "aerodynamic.h"
#include <QVector>
#include "constants.h"

struct LandData
{
    //»змен€ютс€ в процессе полета
    double t, m, gamma, psi, theta, a, B;
    vec r_c;
    vec r_m;
    vec V_c;
    double mt_ok;
    double mt_g;
    float kp;

    //Ќе измен€ютс€ в процессе полета
    double P;
    double Sm;
    double nn;
    double la;
    double mt_ok_n;
    double dm_ok;
    double mt_g_n;
    double dm_g;

    std::array<double,AmountPulses> tpn;
    std::array<double,AmountPulses> tpend;
};

double getVertV(const LandData& data);

double getRealTheta(double theta_breaking, double initialTheta, const LandData& data);

vec get_rm_from_rc(const vec& rc, double nn, double la);

M Mc_m(double nn, double la);

M Mv_c(double gamma, double psi, double theta);

M Mo_v(double a, double B);

double getThetafromVc(const vec& V_c);

vec getP(double t,const LandData& data);

vec getPcfromP(double gamma, double psi, double theta, const vec& P);

vec getRcfromR(double gamma, double psi, double theta, double a, double B, const vec& R);

vec get_gc(const vec& r_m, double nn, double la);

vec get_jc(const vec& r_m, double nn, double la);

vec get_kc(const vec& V_c, const vec& r_m, double nn, double la);

double dVx(double x, double y, double z, double Vx, double Vy, double Vz, double t, const LandData& data);

double dVy(double x, double y, double z, double Vx, double Vy, double Vz, double t, const LandData& data);

double dVz(double x, double y, double z, double Vx, double Vy, double Vz, double t, const LandData& data);

void rungeKutt(LandData& data, double h);
