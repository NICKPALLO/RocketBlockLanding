#pragma once
#include"matrix.h"
#include<vector>

#ifndef Rz 
#define Rz 6371100 // метров
#endif


//std::map<int, double> sndp{
//	{0,1.225},
//	{5,0.7421},
//	{10,0.4176},
//	{15,0.1916},
//	{20,0.08801},
//	{25,0.04048},
//	{30,0.01806},
//	{35,0.00839},
//	{40,0.004},
//	{45,0.00202},
//	{50,0.00107},
//	{55,0.00057},
//	{60,0.00031},
//	{65,0.00016},
//	{70,0.000085},
//	{75,0.000042}
//};

//std::vector<int> srtPh = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 75 };
int srtPh[16] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75 };
double srtP[16] = { 1.225 ,0.7421, 0.4176,  0.1916 , 0.08801 , 0.04048 ,0.01806 ,0.00839 ,0.004 ,0.00202 ,0.00107 , 0.00057 ,0.00031 ,0.00016 ,0.000085 ,0.000042 };


double getp(double h) 
{
	h /= 1000;
	double p;
	if (h >= 75)
	{
		p = 0;
	}
	else if (h>=0&&h<75)
	{
		for (int i = 0; i < 15; ++i)
		{
			if (h >= srtPh[i] && h < srtPh[i + 1])
			{
				p = (h - srtPh[i]) * (srtP[i + 1] - srtP[i]) / (srtPh[i + 1] - srtPh[i]) + srtP[i];
			}
		}
	}
	else {
		p = srtP[0];
	}
	return p;
}

/*vec getR(vec& r_m,vec& V_c, double& Sm, vec& Ca)
{
	vec R;

	double p = getp(r_m.size()-Rz);
	double q = (p * V_c.x * V_c.x) / 2; 


	R.x = -Ca.x * q * Sm;
	R.y = Ca.y * q * Sm;

	return R;
}
*/
