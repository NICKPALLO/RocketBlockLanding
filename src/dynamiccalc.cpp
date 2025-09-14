#include "dynamiccalc.h"

DynamicCalc::DynamicCalc(QObject *parent)
    : QObject{parent}
{}

void DynamicCalc::resetData()
{
    t = 0;
    y1 = 0;
    v1 = 0;
    b1 = 0;
    y = 0;
    v = 0;
    b = 0;
    vec_y1.clear();
    vec_v1.clear();
    vec_b1.clear();
    vec_y.clear();
    vec_v.clear();
    vec_b.clear();
    vec_t.clear();
}

void DynamicCalc::setData(double h_, const ConstData &constData_)
{
    h = h_;
    constData = constData_;
}

QVector<double> &DynamicCalc::get_y()
{
    return vec_y;
}

QVector<double> &DynamicCalc::get_y1()
{
    return vec_y1;
}

QVector<double> &DynamicCalc::get_v()
{
    return vec_v;
}

QVector<double> &DynamicCalc::get_v1()
{
    return vec_v1;
}

QVector<double> &DynamicCalc::get_b()
{
    return vec_b;
}

QVector<double> &DynamicCalc::get_b1()
{
    return vec_b1;
}

QVector<double> &DynamicCalc::get_t()
{
    return vec_t;
}

bool DynamicCalc::startCalculation(const QVector<double>& time, const QVector<double>& H_t,
                                   const QVector<double>& m_t, const QVector<double>& V_t,
                                   const QVector<double>& m_ok_t, const QVector<double>& m_g_t)
{
    resetData();

    double ky0,ky1,ky2,ky3;
    double kb0,kb1,kb2,kb3;
    double kv0,kv1,kv2,kv3;

    double qy0,qy1,qy2,qy3;
    double qb0,qb1,qb2,qb3;
    double qv0,qv1,qv2,qv3;

    double tkon = time[time.size()-1];
    do
    {
        vec_y1.push_back(y1);
        vec_v1.push_back(v1*57.3);
        vec_b1.push_back(b1*57.3);
        vec_y.push_back(y);
        vec_v.push_back(v*57.3);
        vec_b.push_back(b*57.3);
        vec_t.push_back(t);

        qy0=fy1(t,y1,v,b,time,H_t,m_t,V_t,constData);
        ky0=y1;

        qv0=fv1(t,v,y1,b,time,H_t,m_t,V_t,m_ok_t,m_g_t,constData);
        kv0=v1;

        qb0=fb1(b1,b,v1,v,y1,y,constData);
        kb0=b1;


        qy1=fy1(t+h/2,y1+qy0*h/2,v+kv0*h/2,b+kb0*h/2,time,H_t,m_t,V_t,constData);
        ky1=y1+qy0*h/2;

        qv1=fv1(t+h/2,v+kv0*h/2,y1+qv0*h/2,b+kb0*h/2,time,H_t,m_t,V_t,m_ok_t,m_g_t,constData);
        kv1=v1+qv0*h/2;

        qb1=fb1(b1+qb0*h/2,b+kb0*h/2,v1+qv0*h/2,v+kv0*h/2,y1+qy0*h/2,y+ky0*h/2, constData);
        kb1=b1+qb0*h/2;


        qy2=fy1(t+h/2,y1+qy1*h/2,v+kv1*h/2,b+kb1*h/2,time,H_t,m_t,V_t,constData);
        ky2=y1+qy1*h/2;

        qv2=fv1(t+h/2,v+kv1*h/2,y1+qv1*h/2,b+kb1*h/2,time,H_t,m_t,V_t,m_ok_t,m_g_t,constData);
        kv2=v1+qv1*h/2;

        qb2=fb1(b1+qb1*h/2,b+kb1*h/2,v1+qv1*h/2,v+kv1*h/2,y1+qy1*h/2,y+ky1*h/2, constData);
        kb2=b1+qb1*h/2;



        qy3=fy1(t+h,y1+qy2*h,v+kv2*h,b+kb2*h,time,H_t,m_t,V_t,constData);
        ky3=y1+qy2*h;

        qv3=fv1(t+h,v+kv2*h,y1+qv2*h,b+kb2*h,time,H_t,m_t,V_t,m_ok_t,m_g_t,constData);
        kv3=v1+qv2*h;

        qb3=fb1(b1+qb2*h,b+kb2*h,v1+qv2*h,v+kv2*h,y1+qy2*h,y+ky2*h, constData);
        kb3=b1+qb2*h;




        y1=y1+h/6*(qy0+2*qy1+2*qy2+qy3);
        y=y+h/6*(ky0+2*ky1+2*ky2+ky3);

        v1=v1+h/6*(qv0+2*qv1+2*qv2+qv3);
        v=v+h/6*(kv0+2*kv1+2*kv2+kv3);

        b1=b1+h/6*(qb0+2*qb1+2*qb2+qb3);
        b=b+h/6*(kb0+2*kb1+2*kb2+kb3);

        t+=h;
    }while (t<tkon);
    return true;
}
