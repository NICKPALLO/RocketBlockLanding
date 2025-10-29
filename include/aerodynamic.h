#pragma once
#include "matrix.h"
#include "constants.h"

double getLineValue(double y1, double y2, double x1, double x2, double x);
double get_Mah(double& H, double& V);
double get_q(double& H,double& V);
vec getR(const vec& r_m, const vec& V_c, double Sm);
