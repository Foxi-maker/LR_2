#include"Header.h"


double lambda = 1.;
double ThermCondExZero(double uu, double xx)
{
	return lambda / (PI*PI);
}

double ThermCondX(double uu, double xx)
{
	double x1 = 1. / 3, x2 = 3. / 4;
	double k1 = 0.5, k2 = 2.5;

	if (xx <= x1)
		return k1;
	else
		if (xx < x2)
			return k1 * (xx - x2) / (x1 - x2) + k2 * (xx - x1) / (x2 - x1);
		else
			return k2;
}

double alp = .1, bet = 1., gam = 3.;
double ThermCondU(double uu, double xx)
{
	return alp + bet * pow(uu, gam);
}

double cappa = .5, delta = 2.;
double ThermCondExThree(double uu, double xx)
{
	return cappa * pow(uu, delta);
}
