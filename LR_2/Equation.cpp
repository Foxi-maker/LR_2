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
	double tau = 0.1, h = 0.2;
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
		grid[0][pIter] = sin(PI*x);
		//grid[0][pIter] = u0 + x * (L - x);
		x += h;
	}

	for (const auto gx : grid[0])
	{
		out << gx << " ";
	}
	out << "\n";

	//Занчение t_0 по варианту
	//double t0 = .5;
	//double time = tau;
	//for (int tIter = 1; tIter < Nt; tIter++)
	//{
	//	//if (time < t0)
	//	//	grid[tIter][0] = uLeft;
	//	//else
	//	//	grid[tIter][0] = 0;

	//	//grid[tIter].back() = u0;
	//	grid[tIter].front() = 0.;
	//	grid[tIter].back() = 0.;
	//	time += tau;
	//}

	double taucroh = tau / (c * ro * h);
	double crohtau = c * ro * h / tau;
	double aPres, aNext;

	double time = 0.;
	if (sig == 0.)
	{
		//В предположении, что сигма 0 (Явная схема)
		for (int tIter = 0; tIter < Nt - 1; tIter++)
		{
			//Заполнение краевыми значениями на новом временном слое
			//Разве нужно это делать на каждой итерации? (для примера не зависящего от тау - нет)
			time += tau;
			grid[1].front() = 0.;
			grid[1].back() = 0.;

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

			for (const auto gx : grid[0])
			{
				out << gx << " ";
			}
			out << "\n";
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
			grid[1].front() = 0.;
			grid[1].back() = 0.;

			x = h;
			//std::cout << "tIter: " << tIter << "\n";
			//Первая итерация
			aPres = 0.5 * (K(grid[0][1], x) + K(grid[0][0], x - h));
			aNext = 0.5 * (K(grid[0][2], x + h) + K(grid[0][1], x));
			//std::cout << "aPres: " << aPres << ", aNext:" << aNext << "\n";

			A = sigh * aPres;
			B = sigh * aNext;
			C = A + B + crohtau;
			F = crohtau * grid[0][1] + (1 - sig) * (aNext * (grid[0][2] - grid[0][1]) / h - aPres * (grid[0][1] - grid[0][0]) / h);
			//std::cout << "A: " << A << ", B: " << B << ", C: " << C << ", F: " << F << "\n";

			alpha.push_back(B / C);
			betta.push_back((F + A * grid[1][0]) / C);
			//std::cout << "alpha.back(): " << alpha.back() << ", betta.back(): " << betta.back() << "\n";

			for (size_t pIter = 2; pIter < Nx - 2; pIter++)
			{
				std::cout << "\n";
				x += h;
				aPres = 0.5 * (K(grid[0][pIter], x) + K(grid[0][pIter - 1], x - h));
				aNext = 0.5 * (K(grid[0][pIter + 1], x + h) + K(grid[0][pIter], x));
				//std::cout << "aPres: " << aPres << ", aNext:" << aNext << "\n";

				A = sigh * aPres;
				B = sigh * aNext;
				C = A + B + crohtau;
				F = crohtau * grid[0][pIter] + (1 - sig) * (aNext * (grid[0][pIter + 1] - grid[0][pIter]) / h - aPres * (grid[0][pIter] - grid[0][pIter - 1]) / h);
				//std::cout << "A: " << A << ", B: " << B << ", C: " << C << ", F: " << F << "\n";

				temp = C - A * alpha.back();

				alpha.push_back(B / temp);
				betta.push_back((A * betta.back() + F) / temp);
				//std::cout << "alpha.back(): " << alpha.back() << ", betta.back(): " << betta.back() << "\n";
			}
			x += h;
			//Последняя итерация
			aPres = 0.5 * (K(grid[0][Nx - 2], x) + K(grid[0][Nx - 3], x - h));
			aNext = 0.5 * (K(grid[0][Nx - 1], x + h) + K(grid[0][Nx - 2], x));
			//std::cout << "aPres: " << aPres << ", aNext:" << aNext << "\n";
			A = sigh * aPres;
			B = sigh * aNext;
			C = A + B + crohtau;
			F = crohtau * grid[0][Nx - 2] + (1 - sig) * (aNext * (grid[0][Nx - 1] - grid[0][Nx - 2]) / h - aPres * (grid[0][Nx - 2] - grid[0][Nx - 3]) / h);
			//std::cout << "A: " << A << ", B: " << B << ", C: " << C << ", F: " << F << "\n";

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

			for (const auto gx : grid[0])
			{
				out << gx << " ";
			}
			out << "\n";
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
	//Следить за условием устойчивости на tau для Явной схемы
	double tau = 0.0002, h = 0.2;

	size_t Nx = 50;		//число точек=50
	size_t Nt = 10001;

	//Сетка по положению и времени
	std::vector<std::vector<double>> grid(Nt, std::vector<double>(Nx));

	//Заполение значениями начальных и краевых условий
	double x = 0.;
	for (int pIter = 0; pIter < Nx; pIter++)
	{
		grid[0][pIter] = 0.;
		//x += h;
	}

	double time = tau;
	for (int tIter = 1; tIter < Nt; tIter++)
	{
		grid[tIter].front() = u0 * pow(time, 1. / 2);		// 1/delta
		grid[tIter].back() = 0.;
		time += tau;
	}

	double taucroh = tau / (c * ro * h);
	double crohtau = c * ro * h / tau;
	double aPres, aNext;

	if (sig == 0.)
	{
		//В предположении, что сигма 0 (Явная схема)
		for (int tIter = 0; tIter < Nt - 1; tIter++)
		{
			x = h;
			for (int pIter = 1; pIter < Nx - 1; pIter++)
			{
				aPres = 0.5 * (K(grid[tIter][pIter], x) + K(grid[tIter][pIter - 1], x - h));
				aNext = 0.5 * (K(grid[tIter][pIter + 1], x + h) + K(grid[tIter][pIter], x));

				//Явная разностная схема
				grid[tIter + 1][pIter] = grid[tIter][pIter] + taucroh * (aNext * (grid[tIter][pIter + 1] - grid[tIter][pIter]) / h - aPres * (grid[tIter][pIter] - grid[tIter][pIter - 1]) / h);

				x += h;
			}
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
			x = h;
			//Первая итерация
			aPres = 0.5 * (K(grid[tIter][1], x) + K(grid[tIter][0], x - h));
			aNext = 0.5 * (K(grid[tIter][2], x + h) + K(grid[tIter][1], x));
			A = sigh * aPres;
			B = sigh * aNext;
			C = sigh * aPres + sigh * aNext + crohtau;
			F = crohtau * grid[tIter][1] + (1 - sig) * (aNext * (grid[tIter][2] - grid[tIter][1]) / h - aPres * (grid[tIter][1] - grid[tIter][0]) / h);

			alpha.push_back(B / C);
			betta.push_back((F + A * grid[tIter + 1][0]) / C);

			for (size_t pIter = 2; pIter < Nx - 2; pIter++)
			{
				x += h;
				aPres = 0.5 * (K(grid[tIter][pIter], x) + K(grid[tIter][pIter - 1], x - h));
				aNext = 0.5 * (K(grid[tIter][pIter + 1], x + h) + K(grid[tIter][pIter], x));

				A = sigh * aPres;
				B = sigh * aNext;
				C = sigh * aPres + sigh * aNext + crohtau;
				F = crohtau * grid[tIter][pIter] + (1 - sig) * (aNext * (grid[tIter][pIter + 1] - grid[tIter][pIter]) / h - aPres * (grid[tIter][pIter] - grid[tIter][pIter - 1]) / h);

				temp = C - A * alpha.back();

				alpha.push_back(B / temp);
				betta.push_back((A * betta.back() + F) / temp);

			}
			x += h;
			//Последняя итерация
			aPres = 0.5 * (K(grid[tIter][Nx - 2], x) + K(grid[tIter][Nx - 3], x - h));
			aNext = 0.5 * (K(grid[tIter][Nx - 1], x + h) + K(grid[tIter][Nx - 2], x));
			A = sigh * aPres;
			B = sigh * aNext;
			C = sigh * aPres + sigh * aNext + crohtau;
			F = crohtau * grid[tIter][Nx - 2] + (1 - sig) * (aNext * (grid[tIter][Nx - 1] - grid[tIter][Nx - 2]) / h - aPres * (grid[tIter][Nx - 2] - grid[tIter][Nx - 3]) / h);

			temp = C - A * alpha.back();

			grid[tIter + 1][Nx - 2] = (A * betta.back() + F + B * grid[tIter + 1][Nx - 1]) / temp;

			//Обратный ход
			for (size_t pIter = Nx - 3; pIter > 0; pIter--)
			{
				grid[tIter + 1][pIter] = alpha[pIter - 1] * grid[tIter + 1][pIter + 1] + betta[pIter - 1];
			}

			alpha.clear();
			betta.clear();
		}
	}

	//Вывод в файл
	for (const auto gt : grid)
	{
		for (const auto gx : gt)
			out << gx << " ";
		out << "\n";
	}
}

Equation::~Equation()
{
	out.close();
}
