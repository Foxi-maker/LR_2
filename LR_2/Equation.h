#pragma once
#include "Header.h"

class Equation
{
	//Длина стерженя, время рассмотрения
	double L = 1., T = 5.;

	//Удельная теплоемкость, линейная плоность 
	//double c = 0.5, ro = 4;
	//Для примера 3 
	double c = 1., ro = 1.;

	//Коэффициент теплопроводности
	std::function<double(double, double)> K;

	//Начальные условия 
	double uSub;
	//Граничные условия
	double uLeft = 10., uRight = .3;

	//Поток для вывода данных
	static std::ofstream out;

	//Построение сетки
	//void Grid();

public:
	Equation(std::function<double(double, double)>, std::string);
	//Equation();

	void Scheme(double);

	void SchemeEx3(double, double);

	~Equation();
};

