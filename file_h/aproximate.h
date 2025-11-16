#pragma once


namespace aproximate {
    namespace auxiliary_functions {
        // dot of ortagoality
        //degree of polinom
        template<typename P>std::vector<polynomial<P>> GetSystemOfOrthogonalPolynomials(const std::vector<P>& x, int n) {
            std::vector<polynomial<P>> polys;
            if (n < 0) return polys;
            polys.push_back(polynomial<P>(1));
            if (n == 0) return polys;

            //  alpha1
            P alpha1 = std::accumulate(x.begin(), x.end(), P(0)) / static_cast<P>(x.size());

            //P alpha1 = 1;

            //  q1(x) = x - alpha1
            polynomial<P> q1;
            q1.newsize(2);
            q1[0] = -alpha1;
            q1[1] = 1;
            polys.push_back(q1);


            if (n == 1) return polys;

            for (int j = 1; j < n; j++) {
                std::vector<P> qj_vals, qjm1_vals;
                for (auto xi : x) {
                    qj_vals.push_back(polys[j](xi));
                    qjm1_vals.push_back(polys[j - 1](xi));
                }

                P num_alpha = 0, den_alpha = 0;
                for (int i = 0; i < x.size(); i++) {
                    num_alpha += x[i] * qj_vals[i] * qj_vals[i];
                    den_alpha += qj_vals[i] * qj_vals[i];
                }
                P alpha_j1 = num_alpha / den_alpha;

                P num_beta = 0, den_beta = 0;
                for (int i = 0; i < x.size(); i++) {
                    num_beta += x[i] * qj_vals[i] * qjm1_vals[i];
                    den_beta += qjm1_vals[i] * qjm1_vals[i];
                }
                P beta_j = num_beta / den_beta;

                polynomial<P> x_minus_alpha;
                x_minus_alpha.newsize(2);
                x_minus_alpha[0] = -alpha_j1;
                x_minus_alpha[1] = 1;

                polynomial<P> term1 = x_minus_alpha * polys[j];
                polynomial<P> term2 = polys[j - 1] * beta_j;

                polys.push_back(term1 - term2);
            }

            return polys;
        }
        template<typename T>std::vector<std::pair<T, T>> AddNoiseToPoints(
            const std::vector<std::pair<T, T>>& points,
            int measurements_per_point = 3,
            T noise_level = T(0.1)
        ) {


            std::vector<std::pair<T, T>> noisy_points;
            noisy_points.reserve(points.size() * measurements_per_point);

            std::random_device rd;
            std::mt19937 gen(rd());

            for (const auto& point : points) {
                const T x = point.first;
                const T y_clean = point.second;

                std::uniform_real_distribution<T> dist(-noise_level, noise_level);

                for (int i = 0; i < measurements_per_point; ++i) {
                    T noise = dist(gen);
                    T y_noisy = y_clean + noise;

                    noisy_points.emplace_back(x, y_noisy);
                }
            }

            return noisy_points;
        }

    }
    using namespace auxiliary_functions;

    // // // ///
    template<typename P>polynomial<P> LeastSquaresNormal(std::vector<std::pair<P, P>> Array_xy, uint64_t degree_of_the_polynomial = 5) {
        auto V = matrix<double>::vander(convert_pairs_to_vector(Array_xy), degree_of_the_polynomial);
        matrix<double> b((Array_xy.size()), 1);
        int i = 0;
        for (auto& I : convert_pairs_to_vector(Array_xy, false)) {
            b[i][0] = I;
            i++;
        }
        auto matrx_ans = matrixfunction::sanitize_zeros(V.pseudo_inverse() * b, 1e-9);
        polynomial<double> pol_ans; pol_ans.set_deg(degree_of_the_polynomial);
        for (size_t i = 0; i < degree_of_the_polynomial; i++)
        {
            pol_ans[i] = matrx_ans[i][0];
        }
        return pol_ans;
    }

    template<typename P>
    polynomial<P> LeastSquaresOrthogonal(const std::vector<std::pair<P, P>>& points, uint64_t n) {
        std::vector<P> x, y;
        for (const auto& point : points) {
            x.push_back(point.first);
            y.push_back(point.second);
        }

        std::vector<polynomial<P>> phi = ::aproximate::GetSystemOfOrthogonalPolynomials(x, n);

        std::vector<P> c(n, 0.0);

        for (uint64_t k = 0; k < n; k++) {
            P numerator = 0.0;
            P denominator = 0.0;
            //Numerator: projection of the vector y onto the base vector φₖ.
            //Denominator : the square of the norm φₖ at discrete points.

            for (size_t i = 0; i < x.size(); i++) {
                P phi_k_x = phi[k](x[i]);
                numerator += y[i] * phi_k_x;
                denominator += phi_k_x * phi_k_x;
            }
            //k-th basic vector
            c[k] = numerator / denominator;
        }

        polynomial<P> result;
        for (uint64_t k = 0; k < n; k++) {
            result = result + (phi[k] * c[k]);
        }

        return result;
    }

    // // // ///
    namespace output_of_characteristics_for_different_data_size_parameters {
        template<typename P>
        P SumSquaredErrors(const polynomial<P>& poly, const std::vector<std::pair<P, P>>& points) {
            P total_error = 0;

            for (const auto& point : points) {
                P x = point.first;
                //P y_true = ;
                P y_pred = poly(x);
                P error = y_pred - point.second;
                total_error += error * error;
            }

            return total_error;
        }

        template<typename P = double, typename TheTypeOfFunctionThatBeingInterpolated>
        std::stringstream ShowAproximateStatistic(
            TheTypeOfFunctionThatBeingInterpolated F,
            int n = 10,
            int size = 99,
            P a = LDBL_EPSILON - 2 * M_PI,
            P b = 2 * M_PI - LDBL_EPSILON
        )
        {

            std::stringstream Ans;
            const std::string aproximate_name = "least_squares";

            int s1 = std::max(static_cast<int>(std::ceil(std::log10(n))), 3);
            int s3 = 2 + aproximate_name.size() + 2;

            std::string horizontal_line =
                "+" + std::string(s1, '-') + "+" +
                std::string(s3, '-') + "+" +
                std::string(s3 + 5, '-') + "+";

            Ans << "\n"
                << std::left
                << horizontal_line << "\n"
                << std::setw(1) << "|"
                << std::setw(s1) << "(n)"
                << std::setw(1) << "|"
                << std::setw(s3) << "R_" + aproximate_name + "_n"
                << std::setw(1) << "|"
                << std::setw(s3 + 5) << "R_" + aproximate_name + "_opt__n"
                << std::setw(1) << "|\n"
                << horizontal_line << "\n";

            // lines of function
            std::vector<std::function<P(P)>> functions_normal, functions_ortog;
            functions_normal.reserve(n + 1);
            functions_normal.push_back(F);
            functions_ortog.reserve(n + 1);
            functions_ortog.push_back(F);

            // // // //data
            std::vector<std::pair<double, double>> clean_points = ::polynomial_interpolation::nuton2::GeneratePointsEquallySufficient(size, a, b, F);
            auto noisy_points = AddNoiseToPoints(clean_points, 3, 0.2);
            noisy_points.insert(noisy_points.end(), clean_points.begin(), clean_points.end());
            // // // //

            for (uint64_t i = 1; i <= n; i++)
            {

                functions_normal.push_back(LeastSquaresNormal(noisy_points, i));
                functions_ortog.push_back(LeastSquaresOrthogonal(noisy_points, i));
                //
                auto sum_of_squared_errors_norm = SumSquaredErrors(LeastSquaresNormal(noisy_points, i), noisy_points);
                auto sum_of_squared_errors_ortog = SumSquaredErrors(LeastSquaresOrthogonal(noisy_points, i), noisy_points);


                Ans << std::left
                    << std::setw(1) << "|"
                    << std::setw(s1) << i
                    << std::setw(1) << "|"
                    << std::setw(s3) << std::setprecision(12) << sum_of_squared_errors_norm
                    << std::setw(1) << "|"
                    << std::setw(s3 + 5) << std::setprecision(12) << sum_of_squared_errors_ortog
                    << std::setw(1) << "|\n";
            }
            Ans << horizontal_line << "\n";
            std::cout << Ans.str();

#if __has_include(<SFML/Graphics.hpp>)
            DrawFunctions(functions_normal, a, b, noisy_points, aproximate_name + "_normal");
            DrawFunctions(functions_ortog, a, b, noisy_points, aproximate_name + "_orthogonal");
#else
            std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
                << "not working DrawFunctions(functions_normal)\n"
                << "not working DrawFunctions(functions_ortog)\n";
#endif
            return Ans;
        }
    }
    using namespace output_of_characteristics_for_different_data_size_parameters;

}