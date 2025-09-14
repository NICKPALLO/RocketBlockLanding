#pragma once
#include <QVector>
#include "matrix.h"
#include "aerodynamic.h"
#ifndef PI
#define PI 3.14159265
#endif

struct ConstData
{
    double a0, a1, a2, a3, t1, t2;
    double P,Sm,D,L,mk,p_ok,p_g,Xk,Xdn_ok,Xdn_g;
};

double fy1(double t, double y1, double v, double b,const QVector<double>& time,const QVector<double>& H_t,
           const QVector<double>& m_t,const QVector<double>& V_t,const ConstData& constdata);

double fv1(double t, double v, double y1, double b,const QVector<double>& time,
           const QVector<double>& H_t, const QVector<double>& m_t, const QVector<double>& V_t,
           const QVector<double>& m_ok_t, const QVector<double>& m_g_t, const ConstData& constData);

double fb1(double b1, double b, double v1, double v, double y1, double y, const ConstData& constData);
double get_Cya(double Mah, double L, double D);
double get_Cyv(double P, double H,double Sm, double m, double V, double L, double D);
double get_Cyy(double H, double Sm, double m, double V, double L, double D);
double get_Cyb(double P, double m);
double get_Wy(double H);
double get_x1(double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r, double L);
double get_Cvv(double H, double V,double Sm,double Jz,double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r, double L);
double get_Cvb(double P,double Jz,double L,double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r);
double get_Cvy(double H, double V,double Sm,double Jz,double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r, double L);
