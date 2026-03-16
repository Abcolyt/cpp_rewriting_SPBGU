
#include "../file_h/matrix.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <type_traits>
#include <limits>

#pragma once


namespace function_optimization {

	enum class OutputMode { Verbose, Silent };

	template <OutputMode mode = OutputMode::Silent, typename Func>
	double dichotomic_minimize(Func f, double a, double b, double eps) {
		// Проверка валидности входных данных
		if (!std::isfinite(a) || !std::isfinite(b)) {
			throw std::invalid_argument("dichotomic_minimize: a and b must be finite");
		}
		if (a >= b) {
			throw std::invalid_argument("dichotomic_minimize: requires a < b");
		}
		if (eps <= 0.0 || !std::isfinite(eps)) {
			throw std::invalid_argument("dichotomic_minimize: eps must be positive and finite");
		}

		const double delta = eps / 2.0;
		const double h = 1e-8; // шаг для численной производной
		int iter = 0;
		const int max_iter = 1000;

		std::ostringstream buffer;
		if constexpr (mode == OutputMode::Verbose) {
			buffer << "\nThe dichotomy method:\n";
			buffer << "Iter\tx\t\tf(x)\t\t|f'(x)|\n";
		}

		while (b - a >= eps && iter < max_iter) {
			double c = (a + b) / 2.0;
			double x1 = c - delta;
			double x2 = c + delta;

			if (x1 < a) x1 = a;
			if (x2 > b) x2 = b;

			double f1 = f(x1);
			double f2 = f(x2);

			if constexpr (mode == OutputMode::Verbose) {
				double x_mid = (a + b) / 2.0;
				double f_mid = f(x_mid);
				double df = (f(x_mid + h) - f(x_mid - h)) / (2.0 * h);
				double norm_grad = std::abs(df);

				buffer << std::scientific << std::setprecision(8);
				buffer << iter << "\t" << x_mid << "\t" << f_mid << "\t" << norm_grad << "\n";
			}
			++iter;

			if (f1 < f2) {
				b = x2;
			}
			else {
				a = x1;
			}
		}

		double result = (a + b) / 2.0;
		if constexpr (mode == OutputMode::Verbose) {
			double f_res = f(result);
			double df_res = (f(result + h) - f(result - h)) / (2.0 * h);
			buffer << "Final: x_min = " << result << ", f(x_min) = " << f_res
				<< ", |f'(x_min)| = " << std::abs(df_res) << "\n";
			std::cout << buffer.str();
		}
		return result;
	}

	template <OutputMode mode = OutputMode::Silent, typename Func>
	double golden_section_minimize(Func f, double a, double b, double eps) {
		// Проверка валидности входных данных
		if (!std::isfinite(a) || !std::isfinite(b)) {
			throw std::invalid_argument("golden_section_minimize: a and b must be finite");
		}
		if (a >= b) {
			throw std::invalid_argument("golden_section_minimize: requires a < b");
		}
		if (eps <= 0.0 || !std::isfinite(eps)) {
			throw std::invalid_argument("golden_section_minimize: eps must be positive and finite");
		}

		const double phi = (1.0 + std::sqrt(5.0)) / 2.0; // ≈1.618

		const double h = 1e-8; // шаг для численной производной
		int iter = 0;

		std::ostringstream buffer;
		if constexpr (mode == OutputMode::Verbose) {
			buffer << "\nThe Golden Ratio method:\n";
			buffer << "Iter\tx\t\tf(x)\t\t|f'(x)|\n";
		}

		double x1 = b - (b - a) / phi;
		double x2 = a + (b - a) / phi;
		double f1 = f(x1);
		double f2 = f(x2);

		while (b - a > eps) {
			if constexpr (mode == OutputMode::Verbose) {
				double x_mid = (a + b) / 2.0;
				double f_mid = f(x_mid);
				double df = (f(x_mid + h) - f(x_mid - h)) / (2.0 * h);
				double norm_grad = std::abs(df);

				buffer << std::scientific << std::setprecision(8);
				buffer << iter << "\t" << x_mid << "\t" << f_mid << "\t" << norm_grad << "\n";
			}

			if (f1 < f2) {
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = b - (b - a) / phi;
				f1 = f(x1);
			}
			else {
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + (b - a) / phi;
				f2 = f(x2);
			}
			if constexpr (mode == OutputMode::Verbose) {
				++iter;
			}
		}

		double result = (a + b) / 2.0;
		if constexpr (mode == OutputMode::Verbose) {
			double f_res = f(result);
			double df_res = (f(result + h) - f(result - h)) / (2.0 * h);
			buffer << "Final: x_min = " << result << ", f(x_min) = " << f_res
				<< ", |f'(x_min)| = " << std::abs(df_res) << "\n";
			std::cout << buffer.str();
		}
		return result;
	}

	enum class GradMethod { User, NumericalVector, NumericalArray };

	template <typename Func>
	auto make_numerical_gradient_vector(Func f, double h = 1e-9) {
		return [f, h](const std::vector<double>& x) {
			std::vector<double> grad(x.size());
			std::vector<double> xp = x, xm = x;
			for (size_t i = 0; i < x.size(); ++i) {
				xp[i] += h; xm[i] -= h;
				grad[i] = (f(xp) - f(xm)) / (2.0 * h);
				xp[i] = x[i]; xm[i] = x[i];
			}
			return grad;
			};
	}

	template <typename Func, size_t N>
	auto make_numerical_gradient_array(Func f, double h = 1e-9) {
		return [f, h](const std::array<double, N>& x) {
			std::array<double, N> grad;
			std::array<double, N> xp = x, xm = x;
			for (size_t i = 0; i < N; ++i) {
				xp[i] += h; xm[i] -= h;
				grad[i] = (f(xp) - f(xm)) / (2.0 * h);
				xp[i] = x[i]; xm[i] = x[i];
			}
			return grad;
			};
	}

	// Вспомогательный шаблон для проверки, является ли тип std::array
	template <typename T>
	struct is_std_array : std::false_type {};

	template <typename T, size_t N>
	struct is_std_array<std::array<T, N>> : std::true_type {};

	template <typename T>
	inline constexpr bool is_std_array_v = is_std_array<T>::value;

	// ============================================================================
	// ГРАДИЕНТНЫЙ СПУСК / GRADIENT DESCENT
	// ============================================================================
	// 
	// Шаблонные параметры:
	// - method (GradMethod): Способ вычисления градиента
	//   • GradMethod::User            — пользовательский градиент (быстрее, точнее)
	//   • GradMethod::NumericalVector — численный градиент для std::vector<double>
	//   • GradMethod::NumericalArray  — численный градиент для std::array<double, N>
	//
	// - Point: Тип точки (вектор параметров)
	//   • std::vector<double>  — для переменной размерности (heap allocation)
	//   • std::array<double,N> — для фиксированной размерности (stack allocation, быстрее)
	//
	// ВАЖНО: Выбор method должен соответствовать типу Point!
	// - std::vector<double>  → GradMethod::NumericalVector или GradMethod::User
	// - std::array<double,N> → GradMethod::NumericalArray или GradMethod::User
	//
	// Производительность:
	// - Все проверки method выполняются на compile-time через if constexpr
	// - Для std::array компилятор генерирует одну функцию на каждый размер N
	// - Для std::vector генерируется одна универсальная функция
	// ============================================================================
	// Требования к типу Point:
	// - Должен иметь метод size() возвращающий std::size_t
	// - Должен поддерживать operator[] для чтения и записи
	// - Должен быть копируемым
	// Примеры: std::vector<double>, std::array<double, N>
	template <GradMethod method, OutputMode mode = OutputMode::Silent,
		typename Func, typename Point, typename GradFunc = std::nullptr_t>
	Point gradient_descent(Func f, Point x0, GradFunc user_grad = nullptr,
		double eps = 1e-6, int max_iter = 1000) {

		// Проверка валидности входных данных
		if (eps <= 0.0 || !std::isfinite(eps)) {
			throw std::invalid_argument("gradient_descent: eps must be positive and finite");
		}
		if (max_iter <= 0) {
			throw std::invalid_argument("gradient_descent: max_iter must be positive");
		}
		if constexpr (method == GradMethod::User) {
			if (user_grad == nullptr) {
				throw std::invalid_argument("gradient_descent: user_grad is required for GradMethod::User");
			}
		}

		// Выбор способа получения градиента (compile-time)
		auto get_gradient = [&](const Point& x) {
			if constexpr (method == GradMethod::User) {
				return user_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalVector) {
				static_assert(std::is_same_v<Point, std::vector<double>>,
					"NumericalVector requires std::vector<double>");
				auto num_grad = make_numerical_gradient_vector(f);
				return num_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalArray) {
				static_assert(is_std_array_v<Point>,
					"NumericalArray requires std::array<double, N>");
				constexpr size_t N = std::tuple_size_v<Point>;
				auto num_grad = make_numerical_gradient_array<Func, N>(f);
				return num_grad(x);
			}
			};

		Point x = x0;
		double fx = f(x);
		double last_norm = 0.0;

		std::ostringstream buffer;
		if constexpr (mode == OutputMode::Verbose) {
			buffer << "\nGradient descent (detailed output):\n";
			buffer << "Iter\tPoint\t\t\tf(x)\t\t||grad||\n";
			buffer << std::scientific << std::setprecision(8);
		}
		int iter = 0;

		for (; iter < max_iter; ++iter) {
			auto g = get_gradient(x);
			double norm2 = 0;
			for (auto gi : g) norm2 += gi * gi;
			last_norm = std::sqrt(norm2);

			if constexpr (mode == OutputMode::Verbose) {
				buffer << iter << "\t(";
				for (size_t i = 0; i < x.size(); ++i) {
					buffer << x[i];
					if (i + 1 < x.size()) buffer << ", ";
				}
				buffer << ")\t" << fx << "\t" << last_norm << "\n";
			}

			if (norm2 < eps * eps) break;

			// Поиск шага (backtracking)
			double alpha = 1.0;
			Point x_new;
			bool step_found = false;
			while (alpha > 1e-12) {
				x_new = x;
				// -vec-> adding
				for (size_t i = 0; i < x.size(); ++i)
				{
					x_new[i] -= alpha * g[i];
				}

				double f_new = f(x_new);
				if (f_new < fx) {
					step_found = true;
					break;
				}
				alpha *= 0.5;
			}
			if (!step_found) break;

			x = x_new;
			fx = f(x);
		}

		if constexpr (mode == OutputMode::Verbose) {
			buffer << "Final, iter=" << iter << ",\n x = (";
			for (size_t i = 0; i < x.size(); ++i) {
				buffer << x[i];
				if (i + 1 < x.size()) buffer << ", ";
			}
			buffer << "), f(x) = " << fx << ", ||grad|| = " << last_norm << "\n";
			std::cout << buffer.str();
		}

		return x;
	}

	// ============================================================================
	// МЕТОД СОПРЯЖЁННЫХ ГРАДИЕНТОВ / CONJUGATE GRADIENT METHOD
	// ============================================================================
	// 
	// Шаблонные параметры:
	// - method (GradMethod): Способ вычисления градиента
	//   • GradMethod::User            — пользовательский градиент (быстрее, точнее)
	//   • GradMethod::NumericalVector — численный градиент для std::vector<double>
	//   • GradMethod::NumericalArray  — численный градиент для std::array<double, N>
	//
	// - Point: Тип точки (вектор параметров)
	//   • std::vector<double>  — для переменной размерности (heap allocation)
	//   • std::array<double,N> — для фиксированной размерности (stack allocation, быстрее)
	//
	// ВАЖНО: Выбор method должен соответствовать типу Point!
	// - std::vector<double>  → GradMethod::NumericalVector или GradMethod::User
	// - std::array<double,N> → GradMethod::NumericalArray или GradMethod::User
	//
	// Производительность:
	// - Все проверки method выполняются на compile-time через if constexpr
	// - Для std::array компилятор генерирует одну функцию на каждый размер N
	// - Для std::vector генерируется одна универсальная функция
	// - ls_eps: точность линейного поиска (влияет на сходимость!)
	// ============================================================================
	// Требования к типу Point:
	// - Должен иметь метод size() возвращающий std::size_t
	// - Должен поддерживать operator[] для чтения и записи
	// - Должен быть копируемым
	// Примеры: std::vector<double>, std::array<double, N>
	template <GradMethod method, OutputMode mode = OutputMode::Silent,
		typename Func, typename Point, typename GradFunc = std::nullptr_t>
	Point conjugate_gradient(Func f, Point x0,
		double eps = 1e-6, int max_iter = 1000, GradFunc user_grad = nullptr,
		double ls_eps = 1e-6) // точность одномерного поиска
	{
		// Проверка валидности входных данных
		if (eps <= 0.0 || !std::isfinite(eps)) {
			throw std::invalid_argument("conjugate_gradient: eps must be positive and finite");
		}
		if (max_iter <= 0) {
			throw std::invalid_argument("conjugate_gradient: max_iter must be positive");
		}
		if (ls_eps <= 0.0 || !std::isfinite(ls_eps)) {
			throw std::invalid_argument("conjugate_gradient: ls_eps must be positive and finite");
		}
		if constexpr (method == GradMethod::User) {
			if (user_grad == nullptr) {
				throw std::invalid_argument("conjugate_gradient: user_grad is required for GradMethod::User");
			}
		}

		// В основном по Оптимизации Н. Н. Моисеева, страница 73

		// =====================================
		// || Лямбды для упрощения восприятия ||
		// =====================================

		// Лямбда для получения градиента
		auto get_gradient = [&](const Point& x) -> Point {
			if constexpr (method == GradMethod::User) {
				return user_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalVector) {
				static_assert(std::is_same_v<Point, std::vector<double>>,
					"NumericalVector requires std::vector<double>");
				auto num_grad = make_numerical_gradient_vector(f);
				return num_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalArray) {
				static_assert(is_std_array_v<Point>,
					"NumericalArray requires std::array<double, N>");
				constexpr size_t N = std::tuple_size_v<Point>;
				auto num_grad = make_numerical_gradient_array<Func, N>(f);
				return num_grad(x);
			}
			};

		// return p + alpha*p2
		auto vec_add_scaled_vec = [](const Point& p, double alpha, const Point& p2) {
			Point res = p;
			for (size_t i = 0; i < p.size(); ++i) res[i] += alpha * p2[i];
			return res;
			};

		// return -1 * a
		auto anti = [](const Point& a) -> Point {
			Point res = a;
			for (size_t i = 0; i < a.size(); ++i) res[i] = -a[i];
			return res;
			};

		// Вспомогательные операции
		auto dot = [](const Point& a, const Point& b) -> double {
			double res = 0;
			for (size_t i = 0; i < a.size(); ++i) res += a[i] * b[i];
			return res;
			};

		auto norm2 = [&](const Point& a) { return dot(a, a); };

		// ===============
		// || Реализция ||
		// ===============

		// Инициализация
		Point x = x0;
		Point g = get_gradient(x);
		Point d = anti(g);
		double g_norm2 = norm2(g);
		double g_norm = std::sqrt(g_norm2);
		double fx = f(x);

		const size_t N = x.size(); // размерность задачи

		std::ostringstream buffer;
		if constexpr (mode == OutputMode::Verbose) {
			buffer << "\nConjugate gradient method (detailed output):\n";
			buffer << "Iter\tPoint\t\t\tf(x)\t\t||grad||\n";
			buffer << std::scientific << std::setprecision(8);
		}

		int iter = 0;
		for (; iter < max_iter; ++iter) {
			if constexpr (mode == OutputMode::Verbose) {
				buffer << iter << "\t(";
				for (size_t i = 0; i < x.size(); ++i) {
					buffer << x[i];
					if (i + 1 < x.size()) buffer << ", ";
				}
				buffer << ")\t" << fx << "\t" << g_norm << "\n";
			}

			// Проверка сходимости
			if (g_norm2 < eps * eps) break;

			// Одномерная обертка
			auto line_func_from_F = [&](double alpha) -> double {
				Point x_new = x;
				for (size_t i = 0; i < x.size(); ++i) x_new[i] += alpha * d[i];
				return f(x_new);
				};

			double alpha_opt;
			double f0 = fx;
			double f1 = line_func_from_F(1.0);

			double low, high;
			// Допустим у нас почти унимодальная функция, найдем промежуток на котором она такая
			if (f1 < f0) {                                   // Если при α = 1 функция уменьшилась
				// Значит минимум находится правее, за пределами 1 => Нужно расширить интервал.
				double alpha = 1.0;                          // начинаем с α = 1
				double f_prev = f0;                           // значение функции на предыдущем шаге (сначала при α = 0)
				double f_curr = f1;                           // значение функции на текущем шаге (при α = 1)

				while (f_curr < f_prev) {                     // пока функция продолжает убывать
					alpha *= 2.0;                              // удваиваем α
					f_prev = f_curr;                           // предыдущее значение становится текущим
					f_curr = line_func_from_F(alpha);                  // вычисляем функцию при новом α

					if (alpha > 1e10) break;                   // защита от бесконечного цикла 
				}

				// После выхода из цикла f_curr >= f_prev, значит минимум находится между
				// предыдущим α (alpha/2) и текущим α (alpha)

				low = alpha / 2.0;                            // левая граница интервала
				high = alpha;                                 // правая граница интервала
			}
			else {                                             // Если при α = 1 функция не уменьшилась
				// Минимум находится между 0 и 1
				low = 0.0;
				high = 1.0;
			}

			alpha_opt = golden_section_minimize<OutputMode::Silent>(line_func_from_F, low, high, ls_eps);

			//x_new= x + α_opt*d
			Point x_new = vec_add_scaled_vec(x, alpha_opt, d);
			//for (size_t i = 0; i < x.size(); ++i) x_new[i] += alpha_opt * d[i];
			double f_new = f(x_new);

			// Обновляем текущую рабочую точку
			x = x_new;
			fx = f_new;

			// Сохраняем старый градиент
			Point g_old = g;
			double g_old_norm2 = g_norm2;

			// Новый градиент
			g = get_gradient(x);
			g_norm2 = norm2(g);
			g_norm = std::sqrt(g_norm2);

			// Вычисление β по формуле Полака–Рибьера (усложнение метода флетчера Ривса) с рестартом в
			// случае не квадратичного вида функции на n - той итерации
			double beta = 0.0;
			const bool do_restart = ((iter + 1) % N == 0) || (g_norm2 <
				1e-30);

			if (!do_restart && g_old_norm2 > 1e-30) {
				beta = g_norm2 / g_old_norm2;
				if (beta > 10.0) beta = 10.0;

			}

			// Построение новое направление
			// d_(k + 1) ​= −g_(k + 1) + β_k*​d_k
			Point d_new = vec_add_scaled_vec(anti(g), beta, d);


			// Проверка спускового свойства
			if (dot(d_new, g) >= 0) {
				d_new = anti(g);// Сброс на антиградиент
			}
			d = d_new;
		}

		if constexpr (mode == OutputMode::Verbose) {
			buffer << "Final: x = (";
			for (size_t i = 0; i < x.size(); ++i) {
				buffer << x[i];
				if (i + 1 < x.size()) buffer << ", ";
			}
			buffer << "), f(x) = " << fx << ", ||grad|| = " << (g_norm) << "\n";
			std::cout << buffer.str();
		}

		return x;
	}

	// ============================================================================
	// КВАЗИНЬЮТОНОВСКИЙ МЕТОД BFGS / QUASI-NEWTON BFGS METHOD
	// ============================================================================
	// 
	// Шаблонные параметры:
	// - method (GradMethod): Способ вычисления градиента
	//   • GradMethod::User            — пользовательский градиент (быстрее, точнее)
	//   • GradMethod::NumericalVector — численный градиент для std::vector<double>
	//   • GradMethod::NumericalArray  — численный градиент для std::array<double, N>
	//
	// - Point: Тип точки (вектор параметров)
	//   • std::vector<double>  — для переменной размерности (heap allocation)
	//   • std::array<double,N> — для фиксированной размерности (stack allocation, быстрее)
	//
	// ВАЖНО: Выбор method должен соответствовать типу Point!
	// - std::vector<double>  → GradMethod::NumericalVector или GradMethod::User
	// - std::array<double,N> → GradMethod::NumericalArray или GradMethod::User
	//
	// Производительность:
	// - Все проверки method выполняются на compile-time через if constexpr
	// - Для std::array компилятор генерирует одну функцию на каждый размер N
	// - Для std::vector генерируется одна универсальная функция
	// - c_armijo: параметр условия Армихо (0, 1), обычно 1e-4
	// ============================================================================
	// Требования к типу Point:
	// - Должен иметь метод size() возвращающий std::size_t
	// - Должен поддерживать operator[] для чтения и записи
	// - Должен быть копируемым
	// Примеры: std::vector<double>, std::array<double, N>
	template <GradMethod method, OutputMode mode = OutputMode::Silent,
		typename Func, typename Point, typename GradFunc = std::nullptr_t>
	Point bfgs_minimize(Func f, Point x0,
		double eps = 1e-6, int max_iter = 1000, GradFunc user_grad = nullptr,
		double c_armijo = 1e-4) {

		// Проверка валидности входных данных
		if (eps <= 0.0 || !std::isfinite(eps)) {
			throw std::invalid_argument("bfgs_minimize: eps must be positive and finite");
		}
		if (max_iter <= 0) {
			throw std::invalid_argument("bfgs_minimize: max_iter must be positive");
		}
		if (c_armijo <= 0.0 || c_armijo >= 1.0 || !std::isfinite(c_armijo)) {
			throw std::invalid_argument("bfgs_minimize: c_armijo must be in (0, 1)");
		}
		if constexpr (method == GradMethod::User) {
			if (user_grad == nullptr) {
				throw std::invalid_argument("bfgs_minimize: user_grad is required for GradMethod::User");
			}
		}

		// =====================================
		// || Лямбды для упрощения восприятия ||
		// =====================================

		// Лямбда для получения градиента
		auto get_gradient = [&](const Point& x) -> Point {
			if constexpr (method == GradMethod::User) {
				return user_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalVector) {
				static_assert(std::is_same_v<Point, std::vector<double>>,
					"NumericalVector requires std::vector<double>");
				auto num_grad = make_numerical_gradient_vector(f);
				return num_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalArray) {
				static_assert(is_std_array_v<Point>,
					"NumericalArray requires std::array<double, N>");
				constexpr size_t N = std::tuple_size_v<Point>;
				auto num_grad = make_numerical_gradient_array<Func, N>(f);
				return num_grad(x);
			}
			};

		// return p + alpha*p2
		auto vec_add_scaled_vec = [](const Point& p, double alpha, const Point& p2) {
			Point res = p;
			for (size_t i = 0; i < p.size(); ++i) res[i] += alpha * p2[i];
			return res;
			};

		// return -1 * a
		auto anti = [](const Point& a) -> Point {
			Point res = a;
			for (size_t i = 0; i < a.size(); ++i) res[i] = -a[i];
			return res;
			};

		// Вспомогательные операции
		auto dot = [](const Point& a, const Point& b) -> double {
			double res = 0;
			for (size_t i = 0; i < a.size(); ++i) res += a[i] * b[i];
			return res;
			};

		auto norm2 = [&](const Point& a) { return dot(a, a); };

		// ===============
		// || Реализция ||
		// ===============

		// Инициализация
		Point x = x0;
		Point g = get_gradient(x);
		double g_norm2 = norm2(g);
		double g_norm = std::sqrt(g_norm2);
		double fx = f(x);

		const size_t N = x.size(); // размерность задачи

		// Обратный гессиан – единичная матрица N×N
		matrix<double> H = matrix<double>::eye(N);

		std::ostringstream buffer;
		if constexpr (mode == OutputMode::Verbose) {
			buffer << "\nBFGS method (detailed output):\n";
			buffer << "Iter\tPoint\t\t\tf(x)\t\t||grad||\n";
			buffer << std::scientific << std::setprecision(8);
		}

		int iter = 0;
		for (; iter < max_iter; ++iter) {
			if constexpr (mode == OutputMode::Verbose) {
				buffer << iter << "\t(";
				for (size_t i = 0; i < x.size(); ++i) {
					buffer << x[i];
					if (i + 1 < x.size()) buffer << ", ";
				}
				buffer << ")\t" << fx << "\t" << g_norm << "\n";
			}

			// Проверка нормы градиента
			if (g_norm2 < eps * eps) {
				break;
			}

			// Направление спуска p = -H * g
			// Преобразуем Point в matrix<double> для умножения
			matrix<double> g_mat(N, 1);
			for (size_t i = 0; i < N; ++i) {
				g_mat[i][0] = g[i];
			}

			matrix<double> p_mat = H * g_mat * (-1.0);

			// Производная по направлению
			double directional_deriv = matrixfunction::DotProduct(g_mat, p_mat);
			if (directional_deriv >= 0.0) {
				// Направление не спусковое – сбрасываем H в единичную и продолжаем
				H = matrix<double>::eye(N);
				p_mat = H * g_mat * (-1.0);
				directional_deriv = matrixfunction::DotProduct(g_mat, p_mat);
			}

			// Линейный поиск с условием Армихо
			double alpha = 1.0;
			Point x_new;
			double f_new;
			bool step_found = false;

			// Преобразуем p_mat в Point для удобства
			auto p_mat_to_point = [&p_mat, N]() -> Point {
				Point res;
				if constexpr (std::is_same_v<Point, std::vector<double>>) {
					res.resize(N);
				}
				for (size_t i = 0; i < N; ++i) res[i] = p_mat[i][0];
				return res;
				};

			while (alpha > 1e-12) {
				// x_new = x + alpha * p
				x_new = vec_add_scaled_vec(x, alpha, p_mat_to_point());
				f_new = f(x_new);

				if (f_new <= fx + c_armijo * alpha * directional_deriv) {
					step_found = true;
					break;
				}
				alpha *= 0.5;
			}

			if (!step_found) {
				break; // не удалось найти подходящий шаг
			}

			// Новый градиент
			Point g_new = get_gradient(x_new);

			// Разности s и y
			matrix<double> s_mat(N, 1);
			matrix<double> y_mat(N, 1);
			for (size_t i = 0; i < N; ++i) {
				s_mat[i][0] = x_new[i] - x[i];
				y_mat[i][0] = g_new[i] - g[i];
			}

			double ys = matrixfunction::DotProduct(y_mat, s_mat);

			// Обновление BFGS, если условие кривизны выполнено
			if (ys > 1e-12) {
				double rho = 1.0 / ys;

				matrix<double> I = matrix<double>::eye(N);
				matrix<double> s_yT = s_mat * y_mat.transpose();   // N×N
				matrix<double> y_sT = y_mat * s_mat.transpose();   // N×N
				matrix<double> s_sT = s_mat * s_mat.transpose();   // N×N

				matrix<double> A = I - s_yT * rho;
				matrix<double> B = I - y_sT * rho;

				H = A * H * B + s_sT * rho;
			}

			// Переход к следующей итерации
			x = x_new;
			g = g_new;
			fx = f_new;
			g_norm2 = norm2(g);
			g_norm = std::sqrt(g_norm2);
		}

		if constexpr (mode == OutputMode::Verbose) {
			buffer << "Final: x = (";
			for (size_t i = 0; i < x.size(); ++i) {
				buffer << x[i];
				if (i + 1 < x.size()) buffer << ", ";
			}
			buffer << "), f(x) = " << fx << ", ||grad|| = " << g_norm << "\n";
			std::cout << buffer.str();
		}

		return x;
	}

}
