#pragma once
#include "Header.h"

class Equation
{
	//����� ��������, ����� ������������
	double L = 1., T = 5.;

	//�������� ������������, �������� �������� 
	//double c = 0.5, ro = 4;
	//��� ������� 3 
	double c = 1., ro = 1.;

	//����������� ����������������
	std::function<double(double, double)> K;

	//��������� ������� 
	double uSub;
	//��������� �������
	double uLeft = 10., uRight = .3;

	//����� ��� ������ ������
	static std::ofstream out;

	//���������� �����
	//void Grid();

public:
	Equation(std::function<double(double, double)>, std::string);
	//Equation();

	void Scheme(double);

	void SchemeEx3(double, double);

	~Equation();
};

