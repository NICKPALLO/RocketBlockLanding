#include "aerodynamic.h"

double getLineValue(double y1, double y2, double x1, double x2, double x)
{
    return (y2 - y1) * (x - x1) / (x2 - x1) + y1;
}

double get_Mah(double& H, double& V)
{
    double a_zv[9] = {340.26, 299.45, 295.07,304.25,321.78,331.82,319.11,296.76,272.66};
    double a_zv_h;
    int H1=abs(H)/10;
    if(H1>7)
    {
        a_zv_h = 272.66;
    }
    else
    {

        a_zv_h = getLineValue(a_zv[H1],a_zv[H1+1],H1,H1+1,abs(H)/10);
    }
    double Mah = abs(V/a_zv_h);
    Mah=(Mah>5 ? 5 : Mah);
    return Mah;
}

double get_q(double& H, double& V)
{
    double p[10] = {1.225,0.41357,0.0887,0.0179,0.004,0.00108,0.00033,0.000093,0.000021,0.0000034};
    double p_h;
    int H1=abs(H)/10;
    if(H1>9)
    {
        p_h = 0;
    }
    else
    {
        p_h = getLineValue(p[H1],p[H1+1],H1,H1+1,abs(H)/10);
    }
    double q = p_h*V*V/2;
    return q;
}


vec getR(const vec& r_m, const vec& V_c, double Sm)
{
    vec R;
    double H = (r_m.size()-Rz)/1000;
    double V = V_c.size();
    double q = get_q(H,V);
    double Mah=get_Mah(H,V);

    int H1 = abs(H)/10;
    H1=(H1>8 ? 8 : H1);
    int H2 = H1+1;
    double Cax[10];
    Cax[0]=0.0043*pow(Mah,6)-0.0738*pow(Mah,5)+0.4923*pow(Mah,4)-1.5849*pow(Mah,3)+2.3747*pow(Mah,2)-1.0128*Mah+1.0795;
    Cax[1]=0.0046*pow(Mah,6)-0.0779*pow(Mah,5)+0.5175*pow(Mah,4)-1.6573*pow(Mah,3)+2.4707*pow(Mah,2)-1.0589*Mah+1.0874;
    Cax[2]=0.0049*pow(Mah,6)-0.0831*pow(Mah,5)+0.5497*pow(Mah,4)-1.7507*pow(Mah,3)+2.5968*pow(Mah,2)-1.1243*Mah+1.1040;
    Cax[3]=0.0053*pow(Mah,6)-0.0892*pow(Mah,5)+0.5872*pow(Mah,4)-1.8608*pow(Mah,3)+2.7500*pow(Mah,2)-1.2123*Mah+1.1344;
    Cax[4]=0.0044*pow(Mah,6)-0.0757*pow(Mah,5)+0.5003*pow(Mah,4)-1.5870*pow(Mah,3)+2.3132*pow(Mah,2)-0.8926*Mah+1.0791;
    Cax[5]=0.0069*pow(Mah,6)-0.1155*pow(Mah,5)+0.7430*pow(Mah,4)-2.2858*pow(Mah,3)+3.2605*pow(Mah,2)-1.3976*Mah+1.1542;
    Cax[6]=0.0059*pow(Mah,6)-0.0993*pow(Mah,5)+0.6541*pow(Mah,4)-2.0970*pow(Mah,3)+3.2189*pow(Mah,2)-1.6682*Mah+1.2931;
    Cax[7]=0.0086*pow(Mah,6)-0.1438*pow(Mah,5)+0.9394*pow(Mah,4)-2.9746*pow(Mah,3)+4.5409*pow(Mah,2)-2.5973*Mah+1.6095;
    Cax[8]=0.0122*pow(Mah,6)-0.2047*pow(Mah,5)+1.3350*pow(Mah,4)-4.2391*pow(Mah,3)+6.6119*pow(Mah,2)-4.2687*Mah+2.3045;
    Cax[9]=0.0049*pow(Mah,6)-0.0857*pow(Mah,5)+0.5859*pow(Mah,4)-1.9512*pow(Mah,3)+3.0916*pow(Mah,2)-1.6164*Mah+1.2732;

    double CaxH=getLineValue(Cax[H1],Cax[H2],H1,H2,H/10);

    R.x = -CaxH * q * Sm;

    return R;
}
