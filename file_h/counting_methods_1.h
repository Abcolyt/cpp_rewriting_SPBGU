#include "file_h/complex.h"
#include "file_h/fraction.h"
#include "file_h/polynomial.h"
#include "file_h/matrix.h"
#pragma once

//class Polynomial;

namespace counting_methods {

    void executeWithFileInput(std::function<void()> func, const char* filename) {
        FILE* file;

        if (freopen_s(&file, filename, "r", stdin) != 0) {
            std::cerr << "error file open: " << filename << std::endl;
            return; 
        }
        func();
        fclose(file);
    }

    namespace polinomial {
        void polynomial_test() {
            polynomial<Complex> plnm, b;
            std::cin >> plnm;
            using namespace polynomialfunctions;
            auto ans = plnm_roots<Complex>(plnm, Complex(LDBL_EPSILON, LDBL_EPSILON));
            for (size_t i = 0; i < ans.size(); i++)
            {
                std::cout << "root:(" << ans[i].first << ")\niteration for this root:" << ans[i].second << " \n";
            }
        }
    }

    namespace nonlinear_system_with_simple_iterations {
        using Function = std::function<double(const std::vector<double>&)>;

        // A function for solving a nonlinear system of equations using simple iterations
        std::vector<double> nonlinear_system_with_simple_iterations(
            const std::vector<Function>& functions,
            const std::vector<double>& initial_guess,
            double tolerance = 1e-20,
            int max_iterations = LONG_MAX)
        {
            std::vector<double> current_guess = initial_guess;
            int num_functions = functions.size();

            for (int iteration = 0; iteration < max_iterations; ++iteration) {
                std::vector<double> next_guess(num_functions);

                // Calculating the values for the next iteration
                for (int i = 0; i < num_functions; ++i) {
                    next_guess[i] = functions[i](current_guess);
                    std::cout << next_guess[i] << "    ";
                }
                std::cout << "\n";
                // Checking for convergence
                double max_diff = 0.0;
                for (int i = 0; i < num_functions; ++i) {
                    max_diff = std::max(max_diff, std::abs(next_guess[i] - current_guess[i]));
                }

                if (max_diff < tolerance) {
                    return next_guess;// Returning the found solution
                }

                current_guess = next_guess; // Moving on to the next iteration
            }

            throw std::runtime_error("Maximum iterations reached without convergence");
        }
        
        int run_nnssi_with_setted_nonlinear_function() {
            //We define the functions of the system of equations 13
            std::vector<Function> functions1 = {
                [](const std::vector<double>& x) { return 1 - ((1 + LDBL_EPSILON) / (2)) * std::sin(x[1] + 1); },//x=1-((1)/(2))sin(y+1)
                [](const std::vector<double>& x) { return 0.7 - (1 - LDBL_EPSILON) * cos(x[0] - 1); }                   //y=0.7-cos(x-1)
            };

            std::vector<double> initial_guess = { 0, 0 };

            try {

                std::vector<double> solution = nonlinear_system_with_simple_iterations(functions1, initial_guess);

                std::cout << "Solution: ";
                for (double value : solution) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;

            }
            catch (const std::exception& e) {
                std::cerr << "!Error: " << e.what() << std::endl;
            }

            return 0;
        }

        int run_nnssi_with_setted_linear_function() {
            //We define the functions of the system of equations 13
            std::vector<Function> functions1 = {
                [](const std::vector<double>& x) { return (24.4781 - (.0496 * x[1] + .0444 * x[2] + .0393 * x[3])) / 16; },
                [](const std::vector<double>& x) { return (26.0849 - (.0688 * x[0] + .0585 * x[2] + .0534 * x[3])) / 15.1; } ,
                [](const std::vector<double>& x) { return (27.3281 - (.0829 * x[0] + .0777 * x[1] + .0674 * x[3])) / 14.2000; },
                [](const std::vector<double>& x) { return (28.2078 - (.0970 * x[0] + .0918 * x[1] + .0867 * x[2])) / 13.3000;  }
            };

            std::vector<double> initial_guess = { 0.1, 0.1,0.1,0.1 };

            try {

                std::vector<double> solution = nonlinear_system_with_simple_iterations(functions1, initial_guess);

                std::cout << "Solution: ";
                for (double value : solution) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;

            }
            catch (const std::exception& e) {
                std::cerr << "!Error: " << e.what() << std::endl;
            }

            return 0;
        }

    }

    namespace nonlinear_system_with_the_tangent_method {
        using Function = std::function<double(const std::vector<double>&)>;
        using Jacobian = std::function<std::vector<std::vector<double>>(const std::vector<double>&)>;

        // A function for solving a nonlinear system of equations by the Newton method
        std::vector<double> nonlinear_system_with_newton(
            const std::vector<Function>& functions,
            const Jacobian& jacobian,
            const std::vector<double>& initial_guess,
            double tolerance = 1e-9,
            int max_iterations = 100000)
        {
            std::vector<double> current_guess = initial_guess;
            int num_functions = functions.size();

            for (int iteration = 0; iteration < max_iterations; ++iteration) {
                std::vector<double> function_values(num_functions);
                for (int i = 0; i < num_functions; ++i) {
                    function_values[i] = functions[i](current_guess);
                }

                //jacobian value
                std::vector<std::vector<double>> J = jacobian(current_guess);

                //Solve system J * delta = -F для delta
                std::vector<double> delta(num_functions);

                std::vector<std::vector<double>> augmented_matrix(num_functions, std::vector<double>(num_functions + 1));

                for (int i = 0; i < num_functions; ++i) {
                    for (int j = 0; j < num_functions; ++j) {
                        augmented_matrix[i][j] = J[i][j];
                    }
                    augmented_matrix[i][num_functions] = -function_values[i];
                }

                for (int i = 0; i < num_functions; ++i) {
                    double pivot = augmented_matrix[i][i];
                    for (int j = i; j <= num_functions; ++j) {
                        augmented_matrix[i][j] /= pivot;
                    }

                    for (int k = i + 1; k < num_functions; ++k) {
                        double factor = augmented_matrix[k][i];
                        for (int j = i; j <= num_functions; ++j) {
                            augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
                        }
                    }
                }

                for (int i = num_functions - 1; i >= 0; --i) {
                    delta[i] = augmented_matrix[i][num_functions];
                    for (int j = i + 1; j < num_functions; ++j) {
                        delta[i] -= augmented_matrix[i][j] * delta[j];
                    }
                }

                for (int i = 0; i < num_functions; ++i) {
                    current_guess[i] += delta[i];
                }

                double max_diff = 0.0;
                for (const auto& d : delta) {
                    max_diff = std::max(max_diff, std::abs(d));
                }

                if (max_diff < tolerance) {
                    return current_guess; // Returning the found solution
                }
            }

            throw std::runtime_error("Maximum iterations reached without convergence");
        }

        // example
        int nonlinsystem_tangent_method() {
            // We define the functions of the system of equations
            std::vector<Function> functions = {
                [](const std::vector<double>& x) { return std::tan(x[0] * x[1] + 0.4) - x[0] * x[0]; }, // Пример: f1(x1, x2) = x1^2 + x2 - 2
                [](const std::vector<double>& x) { return 0.8 * x[0] * x[0] + 2 * x[1] * x[1] - 1; } // Пример: f2(x1, x2) = x1 - x2^2
            };

            Jacobian jacobian = [](const std::vector<double>& x) {
                return std::vector<std::vector<double>>{
                    { x[1] / (std::cos((5 * x[1] * x[0] + 2) / 5) * std::cos((5 * x[1] * x[0] + 2) / 5)) - 2 * x[0], x[0] / (std::cos((5 * x[1] * x[0] + 2) / 5) * std::cos((5 * x[1] * x[0] + 2) / 5)) },     // df1/dx1, df1/dx2
                    { 8 * x[0],4 * x[1] }     // df2/dx1, df2/dx2
                };
                };

            std::vector<double> initial_guess = { 0.1, 0.1 };

            try {
                
                std::vector<double> solution = nonlinear_system_with_newton(functions, jacobian, initial_guess);

                std::cout << "Solution: \n";
                for (double value : solution) {
                    std::cout << value << " presizion:   " << std::sqrt(functions[0](solution) * functions[0](solution) + functions[1](solution) * functions[1](solution)) << "\n";
                }
                std::cout << std::endl;

            }
            catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << std::endl;
            }

            return 0;
        }
    }

    namespace gaus_method {

        matrix<double> mult(const matrix<double>& left, const matrix<double>& b) {
            //if (left.getcol() != b.getrow() || b.getcol() != 1) {
            //    throw std::invalid_argument("Invalid dimensions for multiplication.");
            //}

            matrix<double> result(left.getrow(), 1);
            for (size_t i = 0; i < left.getrow(); ++i) {
                for (size_t j = 0; j < left.getrow(); ++j) {
                    result[0][i] += left[i][j] * b[0][j];
                    std::cout << "result[0][<<" << i << "]=" << result[0][i] << "\n";
                }
            }
            return result;
        }

        void gaus_solver_linear_sistem() {
            try
            {
                matrix<double> a;
                std::cin >> a;
                std::cout << "\n" << a;

                matrix<double> b;
                std::cin >> b;
                std::cout << "\n" << b;
                std::cout << "\na^-1:\n" << (a.inverse_M());
                std::cout << "\na*b:\n" << mult((a.inverse_M()), b);

                std::cout << "\ncheck:\n" << a * b;
                
            }



            catch (std::exception ex)
            {
                std::cout << "exeption!!What:" << ex.what();
            }
            catch (...)
            {
                std::cout << "unknown error";
            }
        }
    }

    namespace holechi {
        int example() {
            matrix<double> A(4);

            A[0][0] = .1954; A[0][1] = .7700; A[0][2] = 1.3446; A[0][3] = 1.9192;
            A[1][0] = .7700; A[1][1] = 15.1728; A[1][2] = 21.9666; A[1][3] = 28.7604;
            A[2][0] = 1.3446; A[2][1] = 21.9666; A[2][2] = 74.7291; A[2][3] = 93.3867;
            A[3][0] = 1.9192; A[3][1] = 28.7604; A[3][2] = 93.3867; A[3][3] = 208.6609;

            std::cout << "input matrix:\n" << A;
            try {
                matrix<double> A_inv = (A.cholesky().inverse_M()).transpose() * A.cholesky().inverse_M();
                std::cout << " A^{-1}:\n" << A_inv;

                std::cout << "check A^{-1}:\n" << A_inv * A;
            }
            catch (const std::exception& e) {
                std::cerr << "error: " << e.what() << std::endl;
            }

            return 0;
        }
    }
}

namespace function_optimization{
	enum class OutputMode { Verbose, Silent };

	template <OutputMode mode = OutputMode::Silent,typename Func>
	double dichotomic_minimize(Func f, double a, double b, double eps) {
		const double delta = eps / 2.0;
		if (a >= b) throw std::invalid_argument("should be: a < b");
		if (eps <= 0.0) throw std::invalid_argument("should be: eps>=0");

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
		if (a >= b) throw std::invalid_argument("should be: a < b");
		if (eps <= 0.0) throw std::invalid_argument("should be: eps>=0");

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

#if 0

	template <typename Func, typename GradFunc>
	std::pair<double, double> gradient_descent_simple(
		Func f,
		GradFunc grad,
		double x0 = 0.0,
		double y0 = 0.0,
		double eps = 1e-11,
		int max_iter = 1000) {

		double x = x0;
		double y = y0;
		double fx = f(x, y);

		for (int iter = 0; iter < max_iter; ++iter) {
			std::pair<double, double> g = grad(x, y);
			double gx = g.first;
			double gy = g.second;
			double norm2 = gx * gx + gy * gy;

			if (norm2 < eps * eps) {
				break;
			}

			double alpha = 1.0;
			double new_x, new_y, new_f;
			bool step_found = false;

			while (alpha > 1e-12) {
				new_x = x - alpha * gx;
				new_y = y - alpha * gy;
				new_f = f(new_x, new_y);

				if (new_f < fx) {
					step_found = true;
					break;
				}

				alpha *= 0.5;
			}

			if (!step_found) {
				break;
			}

			x = new_x;
			y = new_y;
			fx = new_f;
		}

		return { x, y };
	}

#endif

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

	//градиентный спуск
	template <GradMethod method, OutputMode mode = OutputMode::Silent,
		typename Func, typename Point, typename GradFunc = std::nullptr_t>
	Point gradient_descent(Func f, Point x0, GradFunc user_grad = nullptr,
		double eps = 1e-6, int max_iter = 1000) {

		// Выбор способа получения градиента (compile-time)
		auto get_gradient = [&](const Point& x) {
			if constexpr (method == GradMethod::User) {
				return user_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalVector) {
				static_assert(std::is_same_v<Point, std::vector<double>>,
					"NumericalVector requires std::vector<double>");
				static auto num_grad = make_numerical_gradient_vector(f);
				return num_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalArray) {
				static_assert(is_std_array_v<Point>,
					"NumericalArray requires std::array<double, N>");
				constexpr size_t N = std::tuple_size_v<Point>;
				static auto num_grad = make_numerical_gradient_array<Func, N>(f);
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
				{x_new[i] -= alpha * g[i];}

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


#if 0  
	template <typename Func, typename GradFunc>
	std::pair<double, double> conjugate_gradient(
		Func f,
		double x0 = 0.0,
		double y0 = 0.0,
		double eps = 1e-6,
		int max_iter = 1000,
		GradFunc grad) {

		double x = x0;
		double y = y0;
		double fx = f(x, y);

		auto g = grad(x, y);
		double gx = g.first;
		double gy = g.second;
		double norm2 = gx * gx + gy * gy;

		double dx = -gx;
		double dy = -gy;

		const int N = 2;          

		for (int iter = 0; iter < max_iter; ++iter) {
			if (norm2 < eps * eps) {
				break;
			}

			double alpha = 1.0;
			double new_x, new_y, new_f;
			bool step_found = false;

			while (alpha > 1e-12) {
				new_x = x + alpha * dx;
				new_y = y + alpha * dy;
				new_f = f(new_x, new_y);

				if (new_f < fx) {
					step_found = true;
					break;
				}

				alpha /=2;   
			}

			if (!step_found) {
				break;
			}

			x = new_x;
			y = new_y;
			fx = new_f;

			double gx_old = gx;
			double gy_old = gy;
			double norm2_old = norm2;

			g = grad(x, y);
			gx = g.first;
			gy = g.second;
			norm2 = gx * gx + gy * gy;

			double beta = 0.0;
			if (norm2_old > 0.0) {
				beta = norm2 / norm2_old;
			}

			double new_dx = -gx + beta * dx;
			double new_dy = -gy + beta * dy;

			if (new_dx * gx + new_dy * gy >= 0.0) {
				new_dx = -gx;
				new_dy = -gy;
			}

			if ((iter + 1) % N == 0) {
				new_dx = -gx;
				new_dy = -gy;
			}

			dx = new_dx;
			dy = new_dy;
		}

		return { x, y };
	}
#else
	// функция метода сопряжённых градиентов
	template <GradMethod method, OutputMode mode = OutputMode::Silent,
		typename Func, typename Point, typename GradFunc = std::nullptr_t>
	Point conjugate_gradient(Func f, Point x0, 
		double eps = 1e-6, int max_iter = 1000,GradFunc user_grad = nullptr,
		double ls_eps = 1e-6) // точность одномерного поиска
	{
		// В освном по Оптимизации Н. Н. Моисеева, страница 73

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
				static auto num_grad = make_numerical_gradient_vector(f);
				return num_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalArray) {
				static_assert(is_std_array_v<Point>,
					"NumericalArray requires std::array<double, N>");
				constexpr size_t N = std::tuple_size_v<Point>;
				static auto num_grad = make_numerical_gradient_array<Func, N>(f);
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
			Point x_new = vec_add_scaled_vec(x,alpha_opt,d);
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
			if ((iter + 1) % N != 0) {
				//            <g_(k+1) ; (g_(k+1)−g_k)>
				//	β_k_^PR​ =  -----------------------		
				//				     ||g_k​||^2
				if (g_old_norm2 < 1e-30) beta = 0;
				else beta = (g_norm2 - dot(g, g_old)) / g_old_norm2; /// Эквивалентная формула: dot(g, g) - dot(g, g_old) = dot(g, g - g_old)

				if (beta < 0) beta = 0; // повышает услойчивость?
			}
			else {
				//рестарт
				beta = 0.0; 
			}

			// Построение новое направление
			// d_(k + 1) ​= −g_(k + 1) + β_k*​d_k
			Point d_new = vec_add_scaled_vec(anti(g), beta,d);


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


#endif
	template <typename Func, typename GradFunc>
	std::pair<double, double> bfgs_minimize(
		Func f,
		GradFunc grad,
		double x0,
		double y0,
		double eps = 1e-6,
		int max_iter = 1000,
		double c_armijo = 1e-4) {

		// Начальная точка как вектор-столбец
		matrix<double> x(2, 1);
		x[0][0] = x0;
		x[1][0] = y0;

		double fx = f(x[0][0], x[1][0]);

		// Градиент в начальной точке
		matrix<double> g = grad(x[0][0], x[1][0]); // ожидается размер 2×1

		// Начальное приближение обратного гессиана – единичная матрица 2×2
		matrix<double> H = matrix<double>::eye(2);

		for (int iter = 0; iter < max_iter; ++iter) {
			// Проверка нормы градиента
			double norm2 = matrixfunction::DotProduct(g, g);
			if (norm2 < eps * eps) {
				break;
			}

			// Направление спуска p = - H * g
			matrix<double> p = H * g * (-1.0);

			// Производная по направлению
			double directional_deriv = matrixfunction::DotProduct(g, p);
			if (directional_deriv >= 0.0) {
				// Направление не спусковое – сбрасываем H в единичную и продолжаем
				H = matrix<double>::eye(2);
				continue;
			}

			// Линейный поиск с условием Армихо
			double alpha = 1.0;
			matrix<double> x_new(2, 1);
			double f_new;
			bool step_found = false;

			while (alpha > 1e-12) {
				x_new = x + p * alpha;  // предполагается оператор + и умножение на скаляр
				f_new = f(x_new[0][0], x_new[1][0]);

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
			matrix<double> g_new = grad(x_new[0][0], x_new[1][0]);

			// Разности
			matrix<double> s = x_new - x;
			matrix<double> y_vec = g_new - g; // переименовал, чтобы избежать конфликта с функцией y

			double ys = matrixfunction::DotProduct(y_vec, s);

			// Обновление BFGS, если условие кривизны выполнено
			if (ys > 1e-12) {
				double rho = 1.0 / ys;

				matrix<double> I = matrix<double>::eye(2);
				matrix<double> s_yT = s * y_vec.transpose();   // 2×2
				matrix<double> y_sT = y_vec * s.transpose();   // 2×2
				matrix<double> s_sT = s * s.transpose();       // 2×2

				matrix<double> A = I - s_yT * rho;
				matrix<double> B = I - y_sT * rho;

				H = A * H * B + s_sT * rho;
			}

			// Переход к следующей итерации
			x = x_new;
			g = g_new;
			fx = f_new;
		}

		return { x[0][0], x[1][0] };
	}

	
	// Обобщённая версия BFGS для произвольной размерности
	template <GradMethod method, OutputMode mode = OutputMode::Silent,
		typename Func, typename Point, typename GradFunc = std::nullptr_t>
	Point bfgs_minimize(Func f, Point x0,
		double eps = 1e-6, int max_iter = 1000, GradFunc user_grad = nullptr,
		double c_armijo = 1e-4) {

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
				static auto num_grad = make_numerical_gradient_vector(f);
				return num_grad(x);
			}
			else if constexpr (method == GradMethod::NumericalArray) {
				static_assert(is_std_array_v<Point>,
					"NumericalArray requires std::array<double, N>");
				constexpr size_t N = std::tuple_size_v<Point>;
				static auto num_grad = make_numerical_gradient_array<Func, N>(f);
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
