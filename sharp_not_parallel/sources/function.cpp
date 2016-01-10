#include <cmath>
#include "function.h"

// !!!Guys, I hope that you will make name of variables with logical name!!! \\

double fi_f(double fi_d, double fi_c, double fi_u)
{
	double fi_diff = fi_d - fi_u;
	double fi_c_ = (fi_c - fi_u) / (fi_d - fi_u);
	if (std::abs(fi_diff) < 1e-5)
		return 0.5 * (fi_d + fi_c) - 0.125 * (fi_d - 2 * fi_c + fi_u);
	else
	{
		if (std::abs(fi_u-2*fi_c+fi_d)<=0.3*std::abs(fi_d-fi_u))
			return 0.5*(fi_d+fi_c)-0.125*(fi_d-2*fi_c+fi_u);
		else
		{
			if(fi_c_<=-1||fi_c_>=1.5||(0.35<=fi_c_&&fi_c_<=0.65))
				return 0.5*(fi_d+fi_c)-0.125*(fi_d-2*fi_c+fi_u);
			if(-1<fi_c_&&fi_c_<=0)
			{
				double fi_f_=0.375*fi_c_;
				return fi_u+(fi_d-fi_u)*fi_f_;
			}
			if((0<fi_c_&&fi_c_<0.35)||(0.65<fi_c_&&fi_c_<=1))
			{
				double A =(fi_d*fi_u-pow(fi_c, 2))/(fi_d-2*fi_c+fi_u);
				double B = fi_c - A;
				double exp_C = B/(fi_u-A);
				return A+B*pow(exp_C, 0.5);
			}
			if(1<fi_c_&&fi_c_<=1.5)
			{
            	double fi_f_=fi_c_;
            	return fi_u+(fi_d-fi_u)*fi_f_;
			}
		}
	}
	return 0;
}
