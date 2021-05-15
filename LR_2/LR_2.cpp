#include "Header.h"

std::ofstream Equation::out;

int main()
{
	double cappa = .5, delta = 2., cc = 5;
	double t_u0 = pow(delta * cc * cc / cappa, 1 / delta);

	Equation first(ThermCondExZero, std::string("example0Im"));
	//Equation first(ThermCondX, std::string("example3Ex"));
	//Equation first(ThermCondU);

	//Параметр сигма(от 0 до 1) разностной схемы 
	double sigma = 1.;

	first.Scheme(sigma);
	//first.SchemeEx3(sigma, t_u0);
}

