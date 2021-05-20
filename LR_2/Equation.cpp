#include "Equation.h"

Equation::Equation(std::function<double(double, double)> TCX, std::string na)
{
	K = TCX;

	out.open("OUTPUT\\" + na + ".txt");
	if (!out.is_open())
	{
		std::cout << "Error! File for writing was not open!";
		exit(1);
	}
}

//sigma 0 to 1
void Equation::Scheme(double sig)
{
	//Шаг по времени и по положению
	double tau = 0.0001, h = 0.1;
	//double tau = 0.0025, h = 0.1 / 2;

	int Nx = L / h + 1;
	int Nt = T / tau + 1;

	//Сетка по положению и времени
	std::vector<std::vector<double>> grid(2, std::vector<double>(Nx));

	//Заполение значениями начальных условий
	double u0 = uRight;
	double x = 0;
	for (int pIter = 0; pIter < Nx; pIter++)
	{
		//grid[0][pIter] = sin(PI*x);
		grid[0][pIter] = u0 + x * (L - x);
		x += h;
	}

	for (const auto gx : grid[0])
	{
		out << gx << " ";
	}
	out << "\n";

	double taucroh = tau / (c * ro * h);
	double crohtau = c * ro * h / tau;
	double aPres, aNext;
	double mu, kap;

	double time = 0.;
	if (sig == 0.)
	{
		//В предположении, что сигма 0 (Явная схема)
		for (int tIter = 0; tIter < Nt - 1; tIter++)
		{
			//Заполнение краевыми значениями на новом временном слое
			time += tau;
			if (time < 0.5)
			{
				//Поток слева для 9 варинта 
				aPres = 0.5 * (K(grid[0][1], h) + K(grid[0][0], 0.));
				grid[1].front() = (crohtau / 2. + sig * uLeft + (1 - sig)*(uLeft - aPres / h * (grid[0][1] - grid[0][0]))) / (crohtau / 2. + sig * aPres / h);
				//grid[1].front() = uLeft;
			}
			else
			{
				grid[1].front() = 0.;
			}
			grid[1].back() = uRight;

			x = h;
			for (int pIter = 1; pIter < Nx - 1; pIter++)
			{
				aPres = 0.5 * (K(grid[0][pIter], x) + K(grid[0][pIter - 1], x - h));
				aNext = 0.5 * (K(grid[0][pIter + 1], x + h) + K(grid[0][pIter], x));

				//Явная разностная схема
				grid[1][pIter] = grid[0][pIter] + taucroh * (aNext * (grid[0][pIter + 1] - grid[0][pIter]) / h - aPres * (grid[0][pIter] - grid[0][pIter - 1]) / h);

				x += h;
			}
			grid[0] = grid[1];

			//if (tIter % 10 == 0)
			//{
			for (const auto gx : grid[0])
			{
				out << gx << " ";
			}
			out << "\n";
			//}
		}
	}
	else
	{
		//Массивы для прогонки
		double A, B, C, F;
		std::vector<double> alpha, betta;

		double sigh = sig / h;
		double temp;

		for (size_t tIter = 0; tIter < Nt - 1; tIter++)
		{
			//Заполнение краевыми значениями на новом временном слое
			time += tau;
			if (time < 0.5)
			{
				//Поток слева для 9 варинта 
				aPres = 0.5 * (K(grid[0][1], h) + K(grid[0][0], 0.));
				grid[1].front() = (crohtau / 2. + sig * uLeft + (1 - sig)*(uLeft - aPres / h * (grid[0][1] - grid[0][0]))) / (crohtau / 2. + sig * aPres / h);
				//grid[1].front() = uLeft;
			}
			else
			{
				grid[1].front() = 0.;
			}
			//Поток справа для 9 варинта 
			grid[1].back() = uRight;

			x = h;
			//Первая итерация
			aPres = 0.5 * (K(grid[0][1], x) + K(grid[0][0], x - h));
			aNext = 0.5 * (K(grid[0][2], x + h) + K(grid[0][1], x));

			A = sigh * aPres;
			B = sigh * aNext;
			C = A + B + crohtau;
			F = crohtau * grid[0][1] + (1 - sig) * (aNext * (grid[0][2] - grid[0][1]) / h - aPres * (grid[0][1] - grid[0][0]) / h);

			alpha.push_back(B / C);
			betta.push_back((F + A * grid[1][0]) / C);

			for (size_t pIter = 2; pIter < Nx - 2; pIter++)
			{
				x += h;
				aPres = 0.5 * (K(grid[0][pIter], x) + K(grid[0][pIter - 1], x - h));
				aNext = 0.5 * (K(grid[0][pIter + 1], x + h) + K(grid[0][pIter], x));

				A = sigh * aPres;
				B = sigh * aNext;
				C = A + B + crohtau;
				F = crohtau * grid[0][pIter] + (1 - sig) * (aNext * (grid[0][pIter + 1] - grid[0][pIter]) / h - aPres * (grid[0][pIter] - grid[0][pIter - 1]) / h);

				temp = C - A * alpha.back();

				alpha.push_back(B / temp);
				betta.push_back((A * betta.back() + F) / temp);
			}
			x += h;
			//Последняя итерация
			aPres = 0.5 * (K(grid[0][Nx - 2], x) + K(grid[0][Nx - 3], x - h));
			aNext = 0.5 * (K(grid[0][Nx - 1], x + h) + K(grid[0][Nx - 2], x));

			A = sigh * aPres;
			B = sigh * aNext;
			C = A + B + crohtau;
			F = crohtau * grid[0][Nx - 2] + (1 - sig) * (aNext * (grid[0][Nx - 1] - grid[0][Nx - 2]) / h - aPres * (grid[0][Nx - 2] - grid[0][Nx - 3]) / h);

			temp = C - A * alpha.back();

			grid[1][Nx - 2] = (A * betta.back() + F + B * grid[1][Nx - 1]) / temp;

			//Обратный ход
			for (size_t pIter = Nx - 3; pIter > 0; pIter--)
			{
				grid[1][pIter] = alpha[pIter - 1] * grid[1][pIter + 1] + betta[pIter - 1];
			}

			alpha.clear();
			betta.clear();

			grid[0] = grid[1];

			/*	if (tIter % 5 == 0)
				{*/
			for (const auto gx : grid[0])
			{
				out << gx << " ";
			}
			out << "\n";
			//}
		}
	}

	//Вывод в файл
	//for (const auto gx : grid[0])
	//{
	//	out << gx << " ";
	//}
}

//Специально для примера 3
void Equation::SchemeEx3(double sig, double u0)
{
	//Шаг по времени и по положению
	double tau = 0.0001, h = 0.2;

	size_t Nx = 51;		//число точек
	size_t Nt = 20000;

	//Сетка по положению и времени
	std::vector<std::vector<double>> grid(2, std::vector<double>(Nx));

	//Заполение значениями начальных и краевых условий
	for (int pIter = 0; pIter < Nx; pIter++)
		grid[0][pIter] = 0.;

	for (const auto gx : grid[0])
		out << gx << " ";
	out << "\n";

	double taucroh = tau / (c * ro * h);
	double crohtau = c * ro * h / tau;
	double aPres, aNext;

	//Массивы для прогонки
	double A, B, C, F;
	std::vector<double> alpha, betta;

	double sigh = sig / h;
	double temp;
	double time = 0.;
	double x = 0.;
	for (size_t tIter = 0; tIter < Nt - 1; tIter++)
	{
		time += tau;
		//aPres = 0.5 * (K(grid[0][1], h) + K(grid[0][0], 0.));
		//uLeft = u0 * pow(time, 0.5);
		//grid[1].front() = (crohtau / 2 + sig * uLeft + (1 - sig)*(u0*pow(time - tau, 0.5) - aPres / h * (grid[0][1] - grid[0][0]))) / (crohtau / 2. + sig * aPres / h);
		grid[1].front() = u0 * pow(time, 0.5);
		grid[1].back() = 0.;

		for (size_t Iter = 0; Iter < 3; Iter++)
		{
			x = h;
			//Первая итерация
			aPres = 0.5 * (K(grid[0][1], x) + K(grid[0][0], x - h));
			aNext = 0.5 * (K(grid[0][2], x + h) + K(grid[0][1], x));
			A = sigh * aPres;
			B = sigh * aNext;
			C = A + B + crohtau;
			F = crohtau * grid[0][1] + (1 - sig) * (aNext * (grid[0][2] - grid[0][1]) / h - aPres * (grid[0][1] - grid[0][0]) / h);

			alpha.push_back(B / C);
			betta.push_back((F + A * grid[1][0]) / C);

			for (size_t pIter = 2; pIter < Nx - 2; pIter++)
			{
				x += h;
				aPres = 0.5 * (K(grid[0][pIter], x) + K(grid[0][pIter - 1], x - h));
				aNext = 0.5 * (K(grid[0][pIter + 1], x + h) + K(grid[0][pIter], x));

				A = sigh * aPres;
				B = sigh * aNext;
				C = A + B + crohtau;
				F = crohtau * grid[0][pIter] + (1 - sig) * (aNext * (grid[0][pIter + 1] - grid[0][pIter]) / h - aPres * (grid[0][pIter] - grid[0][pIter - 1]) / h);

				temp = C - A * alpha.back();

				alpha.push_back(B / temp);
				betta.push_back((A * betta.back() + F) / temp);

			}
			x += h;
			//Последняя итерация
			aPres = 0.5 * (K(grid[0][Nx - 2], x) + K(grid[0][Nx - 3], x - h));
			aNext = 0.5 * (K(grid[0][Nx - 1], x + h) + K(grid[0][Nx - 2], x));
			A = sigh * aPres;
			B = sigh * aNext;
			C = A + B + crohtau;
			F = crohtau * grid[0][Nx - 2] + (1 - sig) * (aNext * (grid[0][Nx - 1] - grid[0][Nx - 2]) / h - aPres * (grid[0][Nx - 2] - grid[0][Nx - 3]) / h);

			temp = C - A * alpha.back();

			grid[1][Nx - 2] = (A * betta.back() + F + B * grid[1][Nx - 1]) / temp;

			//Обратный ход
			for (size_t pIter = Nx - 3; pIter > 0; pIter--)
			{
				grid[1][pIter] = alpha[pIter - 1] * grid[1][pIter + 1] + betta[pIter - 1];
			}

			alpha.clear();
			betta.clear();

			grid[0] = grid[1];

			//if (tIter % 100 == 0)
			//{
			//	for (const auto gx : grid[0])
			//		out << gx << " ";
			//	out << "\n";
			//}
		}

	}
	for (const auto gx : grid[0])
		out << gx << " ";
	out << "\n";
}

Equation::~Equation()
{
	out.close();
}
