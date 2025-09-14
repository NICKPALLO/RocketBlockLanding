#include"dynamicRungekutt.h"
#include"rungekutt.h"

double get_Cya(double Mah, double L, double D)
{
    double lambda = L/D;
    double Val;
    Val=(Mah*Mah>1 ? sqrt(Mah*Mah-1)/lambda : -sqrt(1-Mah*Mah)/lambda);


    double Cya=-0.11077*pow(Val,6)+0.14447*pow(Val,5)+0.05078*pow(Val,4)-0.08693*pow(Val,3)-0.01899*pow(Val,2)+0.00876*Val+0.03654;
    return Cya;
    //return Cya;
}

double get_Cyv(double P, double H,double Sm, double m, double V, double L, double D)
{
    double Mah = get_Mah(H,V);
    double Cya = get_Cya(Mah, L, D);
    double q = get_q(H,V);


    double Cyv=-(P+Cya*q*Sm)/m;
    return Cyv;
}

double get_Cyy(double H, double Sm, double m, double V, double L, double D)
{
    double Mah = get_Mah(H,V);
    double Cya = get_Cya(Mah, L, D);
    double q = get_q(H,V);

    double Cyy = Cya*q*Sm/(V*m);
    return Cyy;
}

double get_Cyb(double P, double m)
{
    return (-P/m);
}

double get_Wy(double H)
{
    double Vw[17] = {2,8,15,18,24,12,10,8,6,5,4,3,3,2,2,1,0};//y
    double Hw[17] = {0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,94};//x
    double Wy;
    if(H>=94 || H<=0)
    {
        Wy=0;
    }
    else
    {
        for(int i=0;i<16;i++)
        {
            if(H>=Hw[i] && H<Hw[i+1])
            {
                Wy=getLineValue(Vw[i],Vw[i+1],Hw[i],Hw[i+1],H);
                break;
            }
        }
    }
    return Wy;
}

double get_x1(double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r, double L)
{
    double Xok = Xdn_ok-mok/(2*p_ok*PI*r*r);
    double Xg = Xdn_g-mg/(2*p_g*PI*r*r);
    double Xrb = (Xk*mk+Xok*mok+Xg*mg)/(mk+mok+mg);
    double x1 = abs(Xrb-L/4);
    return x1;
}

double get_Cvv(double H, double V,double Sm,double Jz,double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r, double L)
{
    double Mah = get_Mah(H,V);
    double Cya = get_Cya(Mah, L, 2*r);
    double q = get_q(H,V);
    double x1 = get_x1(mk, mok, mg, p_ok, p_g, Xk, Xdn_ok, Xdn_g, r, L);

    double Cvv = -Cya*q*Sm*x1/Jz;
    return Cvv;
}

double get_Cvb(double P,double Jz,double L,double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r)
{
    double Xok = Xdn_ok-mok/(2*p_ok*PI*r*r);
    double Xg = Xdn_g-mg/(2*p_g*PI*r*r);
    double Xrb = (Xk*mk+Xok*mok+Xg*mg)/(mk+mok+mg);
    double x2 = abs(L-Xrb);

    double Cvb = P*x2/Jz;
    return Cvb;
}

double get_Cvy(double H, double V,double Sm,double Jz,double mk, double mok, double mg,double p_ok, double p_g, double Xk,double Xdn_ok, double Xdn_g, double r, double L)
{
    double Mah = get_Mah(H,V);
    double Cya = get_Cya(Mah, L, 2*r);
    double q = get_q(H,V);
    double x1 = get_x1(mk, mok, mg, p_ok, p_g, Xk, Xdn_ok, Xdn_g, r, L);

    double Cvy = Cya*q*Sm*x1/(Jz*V);
    return Cvy;
}

double fy1(double t, double y1, double v, double b,const QVector<double>& time,const QVector<double>& H_t,
           const QVector<double>& m_t,const QVector<double>& V_t,const ConstData& constdata)
{
    int iteration = 0;
    for (int i = 0;i<time.size()-1;++i)
    {
        if(t>=time[i] && t<time[i]+1)
        {
            iteration=i;
            break;
        }
    }
    double solution;
    double Cyv = get_Cyv(constdata.P, H_t[iteration], constdata.Sm, m_t[iteration], V_t[iteration],constdata.L,constdata.D);
    double Cyy = get_Cyy(H_t[iteration], constdata.Sm, m_t[iteration], V_t[iteration], constdata.L, constdata.D);
    double Cyb = get_Cyb(constdata.P, m_t[iteration]);
    double Wy = get_Wy(H_t[iteration]);

    solution = -Cyv * v - Cyy * y1 - Cyb * b + Cyy * Wy;
    return solution;
}

double fv1(double t, double v, double y1, double b,const QVector<double>& time,
           const QVector<double>& H_t, const QVector<double>& m_t, const QVector<double>& V_t,
           const QVector<double>& m_ok_t, const QVector<double>& m_g_t, const ConstData& constData)
{
    int iteration = 0;
    for (int i = 0;i<time.size()-1;++i)
    {
        if(t>=time[i] && t<time[i]+1)
        {
            iteration=i;
            break;
        }
    }
    double solution;

    double Jz = m_t[iteration]*constData.D*constData.D/16+m_t[iteration]*constData.L*constData.L/12;

    double Cvv = get_Cvv(H_t[iteration], V_t[iteration], constData.Sm, Jz, constData.mk, m_ok_t[iteration],
                         m_g_t[iteration], constData.p_ok, constData.p_g, constData.Xk, constData.Xdn_ok,
                         constData.Xdn_g, constData.D/2, constData.L);

    double Cvy = get_Cvy(H_t[iteration], V_t[iteration], constData.Sm, Jz, constData.mk, m_ok_t[iteration],
                         m_g_t[iteration], constData.p_ok, constData.p_g, constData.Xk, constData.Xdn_ok,
                         constData.Xdn_g, constData.D/2, constData.L);

    double Cvb = get_Cvb(constData.P, Jz, constData.L, constData.mk, m_ok_t[iteration], m_g_t[iteration],
                         constData.p_ok, constData.p_g, constData.Xk,constData.Xdn_ok, constData.Xdn_g, constData.D/2);

    double Wy = get_Wy(H_t[iteration]);
    //Cvy*=-1;

    solution = -Cvv * v - Cvy * y1 - Cvb * b + Cvy * Wy;
    return solution;
}

double fb1(double b1, double b, double v1, double v, double y1, double y, const ConstData& constData)
{
    double solution;
    solution = (constData.a0 * v + constData.a1 * v1 + constData.a2 * y + constData.a3 * y1 - constData.t1 * b1 - b) / constData.t2;
    return solution;
}
