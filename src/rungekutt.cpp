#include"rungekutt.h"


double getVertV(const LandData& data)
{

    double fi = acos(sqrt((data.r_m.x*data.r_m.x+data.r_m.y*data.r_m.y))/data.r_m.size());
    fi=fi*180/PI;//�������
    double th = getThetafromVc(data.V_c);//�������
    if(data.V_c.y>0)
    {
        th=th*PI/180;
        th*=-1;
        th=th*180/PI;
    }
    double z = 180-(90-fi+th);//�������
    z=z*PI/180;
    double Vxy = sqrt(data.V_c.x*data.V_c.x+data.V_c.y*data.V_c.y);
    double VertV;
    VertV=Vxy*cos(z);
    return VertV;
}

double getRealTheta(double theta_breaking, double initialTheta, const LandData& data)
{
    double theta;
    if (data.t<data.tpn[0])
    {
        theta = ((theta_breaking-initialTheta)/data.tpn[0])*data.t+initialTheta;
    }
    else
    {
        theta = getThetafromVc(data.V_c);
        theta = 180-abs(theta);
        theta = theta > 80 ? theta : 80;
    }
    return theta;
}

double getThetafromVc(const vec& V_c)
{
    double theta;
    theta = -acos(V_c.x/sqrt(V_c.x*V_c.x+V_c.y*V_c.y));
    theta = theta * 180 / PI;
    return theta;
}

vec getP(double t,const LandData& data)
{
    vec P;
    if ((t > data.tpn[0] && t<=data.tpend[0]) || (t>data.tpn[1] && t<=data.tpend[1]) || (t>data.tpn[2] && t<=data.tpend[2]))
    {
        P.x = data.P * data.kp;
    }
    return P;
}


vec get_rm_from_rc(const vec& rc, double nn, double la)
{
	double fi_arg = 90 - la;
	vec rg;
	vec rm;
	MY Mc_g(-nn); // ������� �� �.��������� � ����������
	rg = Mc_g * rc;
	rm.x = cos(fi_arg * PI / 180) * rg.z + sin(fi_arg * PI / 180) * rg.y + Rz * cos(la * PI / 180);
	rm.y = rg.x;
	rm.z = cos(fi_arg * PI / 180) * rg.y - sin(fi_arg * PI / 180) * rg.z + Rz * sin(la * PI / 180);
	return rm;
}

M Mc_m(double& nn, double la) //��� ����� �� ������������
{
	MY Mc_g(-nn); // ������� �� �.��������� � ����������

	MZ MZ(90);
	MY MY(180 - la);

	M Mg_m(3, 3); // ������� �� �.��������� � ����������
	Mg_m = MZ * MY;
	Mg_m.T();    // ������� �� ���������� � �.���������
	M Mc_m(3, 3); // ������� �� �.��������� � �.��������� 
	Mc_m = Mg_m * Mc_g;
	return Mc_m;
}


M Mv_c(double gamma, double psi, double theta)
{
    MX Mx(gamma);
    MY My(psi);
    MZ Mz(theta);

	M Mv_c(3, 3); // ������� �� �.��������� � ���������

    Mv_c = Mz * My * Mx;
    //Mv_c = My * Mz * Mx;
	Mv_c.T();  // ������� �� ��������� � �.���������

	return Mv_c;
}

M Mo_v(double a, double B)
{
    MY My(-B);
    MZ Mz(-a);

	M Mo_v(3, 3);

    Mo_v = Mz * My;

	return Mo_v;
}


vec getPcfromP(double gamma, double psi, double theta, const vec& P)
{
	vec Pc;

	Pc = Mv_c(gamma, psi, theta) * P;

	return Pc;
}

vec getRcfromR(double gamma, double psi, double theta, double a, double B, const vec& R)
{
	vec Rc;
	Rc = Mv_c(gamma, psi, theta) * Mo_v(a, B) * R;
	return Rc;
}

vec get_gc(const vec& r_m, double nn, double la)
{
	vec gm;
	gm.x = -9.81 * r_m.x / r_m.size();
	gm.y = -9.81 * r_m.y / r_m.size();
	gm.z = -9.81 * r_m.z / r_m.size();

	MY My(90 - la);
	MZ Mz(nn);

	vec gm2;
    gm2 = Mz * My * gm;

	vec gc(gm2.y, gm2.z, gm2.x);
	return gc;
}

vec get_jc(const vec& r_m, double nn, double la)
{
	vec jl;
	jl.x = -(sqrt(r_m.x * r_m.x + r_m.y * r_m.y) * OMEGA * OMEGA);
	jl.y = 0;
	jl.z = 0;


	double lambda = acos(r_m.x / sqrt(r_m.x * r_m.x + r_m.y * r_m.y));//�������
	lambda = lambda * 180 / PI;//�������

	MZ MzL(-lambda);

	vec jm;

	jm = MzL * jl;


	MY My(90 - la);
	MZ Mz(nn);

	vec jm2;
    jm2 = Mz * My * jm;
	vec jc(jm2.y, jm2.z, jm2.x);

	return jc;
}

vec get_kc(const vec& V_c, const vec& r_m, double nn, double la)
{
    MY My(la-90); //la - 90 � �� �������� ��� ��� ������� �� -1 (������� � �������� �������)
    MZ Mz(-nn);
	M M3(3, 3);
    M3 = Mz * My;

	double lambda = acos(r_m.x / sqrt(r_m.x * r_m.x + r_m.y * r_m.y));//�������
	lambda = lambda * 180 / PI;//�������

	MZ MzL(lambda);

    vec Vc1;
	Vc1 = MzL * M3 * V_c;

    //vec Vl(Vc1.y, Vc1.z, Vc1.x);
    vec Vl(Vc1.z, Vc1.x, Vc1.y);


	if (r_m.x == 0 || r_m.y == 0)
	{
		vec kc(0, 0, 0);
		return kc;
	}
	else
	{
		vec kl;
		kl.x = -2 * OMEGA * Vl.y;
		kl.y = 2 * OMEGA * Vl.x;
		kl.z = 0;

		MzL.T();
		M3.T();

		vec kl1;

		kl1 = M3 * MzL * kl;

		vec kc(kl1.y, kl1.z, kl1.x);

		return kc;
	}
}


double dVx(double x, double y, double z, double Vx, double Vy, double Vz, double t,const LandData& data)
{	double solution;

	vec r_c(x, y, z);
	vec r_m = get_rm_from_rc(r_c,data.nn,data.la);
	vec V_c(Vx, Vy, Vz);
    double theta = getThetafromVc(data.V_c);
    vec R = getR(r_m, V_c,data.Sm);
    vec Rc = getRcfromR(data.gamma,data.psi,theta,data.a,data.B, R);
    vec P = getP(t,data);
    vec Pc = getPcfromP(data.gamma, data.psi, theta, P);
	vec gc;
	vec jc;
	vec kc;

	gc = get_gc(r_m, data.nn, data.la);
	jc = get_jc(r_m, data.nn, data.la);
	kc = get_kc(V_c, r_m, data.nn, data.la);
    solution = Pc.x / data.m + Rc.x / data.m + gc.x - jc.x - kc.x;


	return solution;
}


double dVy(double x, double y, double z, double Vx, double Vy, double Vz, double t, const LandData& data)
{
	double solution;

	vec r_c(x, y, z);
	vec r_m = get_rm_from_rc(r_c, data.nn, data.la);
	vec V_c(Vx, Vy, Vz);
    double theta = getThetafromVc(data.V_c);
    vec R = getR(r_m, V_c, data.Sm);
    vec Rc = getRcfromR(data.gamma, data.psi, theta, data.a, data.B, R);
    vec P = getP(t,data);
    vec Pc = getPcfromP(data.gamma, data.psi, theta, P);
	vec gc;
	vec jc;
	vec kc;

	gc = get_gc(r_m, data.nn, data.la);
	jc = get_jc(r_m, data.nn, data.la);
	kc = get_kc(V_c, r_m, data.nn, data.la);
    solution = Pc.y / data.m + Rc.y / data.m + gc.y - jc.y - kc.y;

	return solution;
}

double dVz(double x, double y, double z, double Vx, double Vy, double Vz, double t, const LandData& data)
{
	double solution;

	vec r_c(x, y, z);
	vec r_m = get_rm_from_rc(r_c, data.nn, data.la);
	vec V_c(Vx, Vy, Vz);
    double theta = getThetafromVc(data.V_c);
    vec R = getR(r_m, V_c, data.Sm);
    vec Rc = getRcfromR(data.gamma, data.psi, theta, data.a, data.B, R);
    vec P = getP(t,data);
    vec Pc = getPcfromP(data.gamma, data.psi, theta, P);

	vec gc;
	vec jc;
	vec kc;

	gc = get_gc(r_m, data.nn, data.la);
	jc = get_jc(r_m, data.nn, data.la);
	kc = get_kc(V_c, r_m, data.nn, data.la);

    solution = Pc.z / data.m + Rc.z / data.m+ gc.z - jc.z - kc.z;

	return solution;
}



void rungeKutt(LandData& data, double h)
{
    double K0x, K1x, K2x, K3x, q0x, q1x, q2x, q3x;
    double K0y, K1y, K2y, K3y, q0y, q1y, q2y, q3y;
    double K0z, K1z, K2z, K3z, q0z, q1z, q2z, q3z;

    K0x = data.V_c.x;
    q0x = dVx(data.r_c.x, data.r_c.y, data.r_c.z, data.V_c.x, data.V_c.y, data.V_c.z, data.t, data);


    K0y = data.V_c.y;
    q0y = dVy(data.r_c.x, data.r_c.y, data.r_c.z, data.V_c.x, data.V_c.y, data.V_c.z, data.t, data);

    K0z = data.V_c.z;
    q0z = dVz(data.r_c.x, data.r_c.y, data.r_c.z, data.V_c.x, data.V_c.y, data.V_c.z, data.t, data);

    K1x = data.V_c.x + q0x * h / 2;
    q1x = dVx(data.r_c.x + K0x * h / 2, data.r_c.y + K0y * h / 2, data.r_c.z + K0z * h / 2, data.V_c.x + q0x * h / 2, data.V_c.y + q0y * h / 2, data.V_c.z + q0z * h / 2, data.t + h / 2, data);

    K1y = data.V_c.y + q0y * h / 2;
    q1y = dVy(data.r_c.x + K0x * h / 2, data.r_c.y + K0y * h / 2, data.r_c.z + K0z * h / 2, data.V_c.x + q0x * h / 2, data.V_c.y + q0y * h / 2, data.V_c.z + q0z * h / 2, data.t + h / 2, data);

    K1z = data.V_c.z + q0z * h / 2;
    q1z = dVz(data.r_c.x + K0x * h / 2, data.r_c.y + K0y * h / 2, data.r_c.z + K0z * h / 2, data.V_c.x + q0x * h / 2, data.V_c.y + q0y * h / 2, data.V_c.z + q0z * h / 2, data.t + h / 2, data);

    K2x = data.V_c.x + q1x * h / 2;
    q2x = dVx(data.r_c.x + K1x * h / 2, data.r_c.y + K1y * h / 2, data.r_c.z + K1z * h / 2, data.V_c.x + q1x * h / 2, data.V_c.y + q1y * h / 2, data.V_c.z + q1z * h / 2, data.t + h / 2, data);


    K2y = data.V_c.y + q1y * h / 2;
    q2y = dVy(data.r_c.x + K1x * h / 2, data.r_c.y + K1y * h / 2, data.r_c.z + K1z * h / 2, data.V_c.x + q1x * h / 2, data.V_c.y + q1y * h / 2, data.V_c.z + q1z * h / 2, data.t + h / 2, data);


    K2z = data.V_c.z + q1z * h / 2;
    q2z = dVz(data.r_c.x + K1x * h / 2, data.r_c.y + K1y * h / 2, data.r_c.z + K1z * h / 2, data.V_c.x + q1x * h / 2, data.V_c.y + q1y * h / 2, data.V_c.z + q1z * h / 2, data.t + h / 2, data);


    K3x = data.V_c.x + q2x * h;
    q3x = dVx(data.r_c.x + K2x * h, data.r_c.y + K2y * h, data.r_c.z + K2z * h, data.V_c.x + q2x * h, data.V_c.y + q2y * h, data.V_c.z + q2z * h, data.t + h, data);


    K3y = data.V_c.y + q2y * h;
    q3y = dVy(data.r_c.x + K2x * h, data.r_c.y + K2y * h, data.r_c.z + K2z * h, data.V_c.x + q2x * h, data.V_c.y + q2y * h, data.V_c.z + q2z * h, data.t + h, data);

    K3z = data.V_c.z + q2z * h;
    q3z = dVz(data.r_c.x + K2x * h, data.r_c.y + K2y * h, data.r_c.z + K2z * h, data.V_c.x + q2x * h, data.V_c.y + q2y * h, data.V_c.z + q2z * h, data.t + h, data);


    vec next_r_c;
    vec next_V_c;

    next_r_c.x = data.r_c.x + h / 6 * (K0x + 2 * K1x + 2 * K2x + K3x);
    next_V_c.x = data.V_c.x + h / 6 * (q0x + 2 * q1x + 2 * q2x + q3x);

    next_r_c.y = data.r_c.y + h / 6 * (K0y + 2 * K1y + 2 * K2y + K3y);
    next_V_c.y = data.V_c.y + h / 6 * (q0y + 2 * q1y + 2 * q2y + q3y);

    next_r_c.z = data.r_c.z + h / 6 * (K0z + 2 * K1z + 2 * K2z + K3z);
    next_V_c.z = data.V_c.z + h / 6 * (q0z + 2 * q1z + 2 * q2z + q3z);

    data.r_c = next_r_c;
    data.V_c = next_V_c;
    data.t +=h;

    if ((data.t>data.tpn[0]+h && data.t<data.tpend[0]+h) || (data.t>data.tpn[1]+h && data.t<data.tpend[1]+h) || (data.t>data.tpn[2]+h && data.t<data.tpend[2]+h))
    {
        data.mt_g = data.mt_g - data.dm_g * data.kp * h;
        data.mt_ok = data.mt_ok - data.dm_ok * data.kp * h;
        data.m = data.m - (data.dm_g + data.dm_ok) * data.kp * h;
    }

	data.r_m = get_rm_from_rc(data.r_c, data.nn, data.la);

}
