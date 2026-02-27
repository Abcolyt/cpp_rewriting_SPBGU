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
	template <typename Func>
	double dichotomic_minimize(Func f, double a, double b, double eps, const double delta = eps / 2.0) {
		if (a >= b) {
			throw std::invalid_argument("should be: a < b");
		}
		if (eps <= 0.0) {
			throw std::invalid_argument("should be: eps>=0");
		}


		while (b - a >= eps) {
			double c = (a + b) / 2.0;
			double x1 = c - delta;
			double x2 = c + delta;

			if (x1 < a) x1 = a;
			if (x2 > b) x2 = b;

			double f1 = f(x1);
			double f2 = f(x2);

			if (f1 < f2) {
				b = x2;         
			}
			else {
				a = x1;          
			}
		}

		return (a + b) / 2.0;
	}

	template <typename Func>
	double golden_section_minimize(Func f, double a, double b, double eps) {
		if (a >= b) {
			throw std::invalid_argument("should be: a < b");
		}
		if (eps <= 0.0) {
			throw std::invalid_argument("should be: eps>=0");
		}

		const double phi = (1.0 + std::sqrt(5.0)) / 2.0;  // ≈ 1.618

		double x1 = b - (b - a) / phi;
		double x2 = a + (b - a) / phi;
		double f1 = f(x1);
		double f2 = f(x2);

		while (b - a > eps) {
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
		}

		return (a + b) / 2.0;
	}

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

=			if (norm2 < eps * eps) {
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

	template <typename Func, typename GradFunc>
	std::pair<double, double> conjugate_gradient(
		Func f,
		GradFunc grad,
		double x0 = 0.0,
		double y0 = 0.0,
		double eps = 1e-6,
		int max_iter = 1000) {

		double x = x0;
		double y = y0;
		double fx = f(x, y);

		auto g = grad(x, y);
		double gx = g.first;
		double gy = g.second;
		double norm2 = gx * gx + gy * gy;

		double dx = -gx;
		double dy = -gy;

		const int n = 2;          

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

			if ((iter + 1) % n == 0) {
				new_dx = -gx;
				new_dy = -gy;
			}

			dx = new_dx;
			dy = new_dy;
		}

		return { x, y };
	}

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
			double norm2 = DotProduct(g, g);
			if (norm2 < eps * eps) {
				break;
			}

			// Направление спуска p = - H * g
			matrix<double> p = H * g;
			p[0][0] = -p[0][0];
			p[1][0] = -p[1][0];

			// Производная по направлению
			double directional_deriv = DotProduct(g, p);
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

			double ys = DotProduct(y_vec, s);

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

}
