#include <cmath>
#pragma once

// |DifferencialMethod| = 5, 5 different sheme
// 
// RK2 <=> RungeKutta2ndOrder
// RK3 <=> RungeKutta3ndOrder1
// RK4 <=> RungeKutta4ndOrder1 <=> ClassicalRungeMechod
// RungeKutta4ndOrder2 <=> GillFormula
// 
enum class DifferencialMethod {
    RungeKutta2ndOrder = 0,
    RK2 =0,
    RungeKutta3ndOrder1 = 1,
    RK3 = 1,
    RungeKutta3ndOrder2 = 2,
    RungeKutta4ndOrder1 = 3,
    RK4 = 3,
    ClassicalRungeMechod = 3,
    RungeKutta4ndOrder2 = 4,
    GillFormula = 4
};

enum class ErrorMethod {
    RungeGlobal,
    RungeLocal
};

namespace details {
    //the values of the degrees of accuracy of the corresponding method
    const std::unordered_map<DifferencialMethod, uint16_t> order_of_accuracy_differencial_method{
    {DifferencialMethod::RungeKutta2ndOrder,2},
    {DifferencialMethod::RungeKutta3ndOrder1,3},
    {DifferencialMethod::RungeKutta3ndOrder2,3},
    {DifferencialMethod::RungeKutta4ndOrder1,4},
    {DifferencialMethod::RungeKutta4ndOrder2,4}
    };
    //Gives an initial, approximate value of the magnitude of the integration step (See the Runge Kutta method)
    template<typename Targ, typename Tresult>
    Targ GetInitialStepSize(std::function<Tresult(Targ, Tresult)> f,// Right-hand side function f(x, y)
        Targ x0,                                                    // Initial x
        Tresult y0,                                                 // Initial y
        const Targ epsilon,
        const Targ x_target,                                        // Target x value
        int64_t s
    ) {
        Targ h;
        if constexpr (std::is_arithmetic_v<Tresult>) {
            auto delta = std::pow(std::max(x0, x_target), s + 1) + std::pow(std::abs(f(x0, y0)), s + 1);
            h = static_cast<Targ>(std::pow(epsilon / delta, 1 / static_cast<Targ>(s + 1)));
            //one Euler step
            auto delta2 = std::pow(std::max(x0 + h, x_target), s + 1) + std::pow(std::abs(f(x0 + h, y0 + h * f(x0, y0))), s + 1);
            auto h2 = static_cast<Targ>(std::pow(epsilon / delta2, 1 / static_cast<Targ>(s + 1)));

            if (h > h2) {
                h = h2;
            }
        }
        else if constexpr (is_matrix<Tresult>::value) {
            auto delta = std::max(x0, x_target) + std::pow(f(x0, y0).norm(), s + 1);
            h = static_cast<Targ>(std::pow(epsilon / delta, 1 / (s + 1)));
            //one Euler step
            delta = std::pow(std::max(x0 + h, x_target), s + 1) + std::pow(f(x0 + h, y0 + h * f(x0, y0)).norm(), s + 1);
            auto h2 = static_cast<Targ>(std::pow(epsilon / delta, 1 / (s + 1)));

            if (h > h2) {
                h = h2;
            }
        }
        else {
            static_assert(false, "Unsupported type for Tresult");
        }


        return h;
    }

    // gives the optimal step size for a given rho and a previously specified
    // step that meets the accuracy requirements, and the values of the error estimate
    // using the Runge method.
    template<typename Targ>
    Targ GetEpsStepSize(Targ h, Targ epsilon, Targ runge_error, uint16_t s) {
        return (h / 2) * (std::pow(epsilon / runge_error, 1.0 / s));
    }

    // gives an error estimate using the runge method
    template<typename T>
    auto RungeError(const T& y1, const T& y2, uint64_t s) {
        if constexpr (is_matrix<T>::value) {
            //std::cout << (y1 - y2) <<'\n';
            return (y1 - y2).norm() / (std::pow(2, s) - 1);
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return std::abs(y1 - y2) / (std::pow(2, s) - 1);
        }
        else {
            static_assert(sizeof(T) == 0, "Unsupported type for RungeError");
        }
    }

    //gives a point-by-point representation of the form of an integral curve, with points equidistant from each other
    template<DifferencialMethod method, typename Targ, typename Tresult>
    std::vector<std::pair<Targ, Tresult>> SolveASystemOfOrdinaryDifferentialEquationsEqualSteps(
        std::function<Tresult(Targ, Tresult)> f,// Right-hand side function f(x, y)
        Targ x0,                      // Initial x
        Tresult y0,                   // Initial y
        Targ x_target,                      // target x value 
        const uint64_t number_of_steps,
        Targ c2 = 1.0 / 10                  // Parameter ξ for RK2 scheme
    ) {
        //const uint64_t number_of_steps=std::ceil(1.0 / GetInitialStepSize(f, x0, y0, rho, x_target, order_of_accuracy_of_the_method[method]));
        Targ& x_current = x0;
        Tresult& y_current = y0;
        Targ h = (x_target - x0) / (number_of_steps);
        if constexpr (method == DifferencialMethod::RungeKutta2ndOrder) {

            // coefficients
            Targ a21 = c2;
            Targ b2 = static_cast<Targ>(1) / (2 * c2);
            Targ b1 = 1 - static_cast<Targ>(1) / (2 * c2);

            //Ans
            std::vector<std::pair<Targ, Tresult>> solution;
            solution.push_back({ x0, y0 });


            while (x_current < x_target) {
                // RK2 stages
                Tresult k1 = h * f(x_current, y_current);
                Tresult k2 = h * f(x_current + c2 * h, y_current + a21 * k1);

                // Update solution
                y_current = y_current + b1 * k1 + b2 * k2;
                x_current += h;

                solution.push_back({ x_current, y_current });
            }
            return solution;
        }
        if constexpr (method == DifferencialMethod::RungeKutta3ndOrder1) {

            //Ans
            std::vector<std::pair<Targ, Tresult>> solution;
            solution.push_back({ x0, y0 });



            Tresult k1, k2, k3;
            while (x_current < x_target) {
                // RK2 stages
                k1 = h * f(x_current, y_current);
                k2 = h * f(x_current + (1.0 / 2) * h, y_current + (1.0 / 2) * k1);
                k3 = h * f(x_current + h, y_current - k1 + 2 * k2);

                // Update solution
                y_current = y_current + (k1 + 4 * k2 + k3) / 6;
                x_current += h;

                solution.push_back({ x_current, y_current });
            }

            return solution;
        }
        if constexpr (method == DifferencialMethod::RungeKutta3ndOrder2) {

            //Ans
            std::vector<std::pair<Targ, Tresult>> solution;
            solution.push_back({ x0, y0 });

            

            Tresult k1, k2, k3;
            while (x_current < x_target) {
                // RK2 stages
                k1 = h * f(x_current, y_current);
                k2 = h * f(x_current + (1.0 / 3) * h, y_current + (1.0 / 3) * k1);
                k3 = h * f(x_current + (2.0 / 3) * h, y_current + (2.0 / 3) * k2);

                // Update solution
                y_current = y_current + (k1 + 3 * k3) / 4;
                x_current += h;

                solution.push_back({ x_current, y_current });
            }

            return solution;
        }
        if constexpr (method == DifferencialMethod::RungeKutta4ndOrder1) {

            //Ans
            std::vector<std::pair<Targ, Tresult>> solution;
            solution.push_back({ x0, y0 });

           

            Tresult k1, k2, k3, k4;
            while (x_current < x_target) {
                // RK2 stages
                k1 = h * f(x_current, y_current);
                k2 = h * f(x_current + (1.0 / 2) * h, y_current + (1.0 / 2) * k1);
                k3 = h * f(x_current + (1.0 / 2) * h, y_current + (1.0 / 2) * k2);
                k4 = h * f(x_current + h, y_current + k3);

                // Update solution
                y_current = y_current + (k1 + 2 * (k2 + k4) + k3) / 6;
                x_current += h;

                solution.push_back({ x_current, y_current });
            }

            return solution;
        }
        if constexpr (method == DifferencialMethod::RungeKutta3ndOrder2) {

            //Ans
            std::vector<std::pair<Targ, Tresult>> solution;
            solution.push_back({ x0, y0 });

            

            Tresult k1, k2, k3, k4;

            auto qrt = (std::sqrt(2) - 1);
            while (x_current < x_target) {
                // RK2 stages
                k1 = h * f(x_current, y_current);
                k2 = h * f(x_current + (1.0 / 2) * h, y_current + (1.0 / 2) * k1);
                k3 = h * f(x_current + (1.0 / 2) * h, y_current + (1.0 / 2) * qrt * k2 + (qrt / 2) * k2);
                k4 = h * f(x_current + h, y_current - (1 / std::sqrt(2)) * k2 + (1 + 1 / std::sqrt(2)) * k3);
                // Update solution
                y_current = y_current + (k1 + 3 * k3) / 4;
                x_current += h;

                solution.push_back({ x_current, y_current });
            }

            return solution;
        }
        return {};
    }
}


//Solving the Cauchy problem while maintaining full error and accuracy
template<DifferencialMethod method, typename Targ, typename Tresult>
std::vector<std::pair<Targ, Tresult>> SolveASystemOfOrdinaryDifferentialEquationsEqualStepsWithATotalError(
    std::function<Tresult(Targ, Tresult)> f,// Right-hand side function f(x, y)
    const Targ x0,                      // Initial x
    const Tresult y0,                   // Initial y
    Targ x_target,                      // target x value 
    const Targ epsilon,                 //Accuracy
    Targ c2 = 1.0 / 10                  /*Parameter ξ for RK2 scheme*/) {
    using namespace details;
    uint64_t number_of_steps = std::ceil(1.0 / GetInitialStepSize(f, x0, y0, epsilon, x_target, order_of_accuracy_differencial_method.at(method)));

    std::vector<std::pair<Targ, Tresult>> solve1, solve2;
    do
    {
        solve1 = SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method>(f, x0, y0, x_target, number_of_steps, c2);
        solve2 = SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method>(f, x0, y0, x_target, number_of_steps * 2, c2);
        number_of_steps *= 2;
    } while (RungeError(solve1.back().second, solve2.back().second, order_of_accuracy_differencial_method.at(method)) > epsilon);

    if (RungeError(solve1.back().second, solve2.back().second, order_of_accuracy_differencial_method.at(method)) < epsilon) {
        return solve2;
    }
    else {
        return { {0,0} };
    }
}


//Solving the Cauchy problem while maintaining local error and accuracy
template<DifferencialMethod method, typename Targ, typename Tresult>
std::vector<std::pair<Targ, Tresult>> SolveASystemOfOrdinaryDifferentialEquationsDifferentStepsWithALocalError(
    std::function<Tresult(Targ, Tresult)> f,
    const Targ x0,
    const Tresult y0,
    Targ x_target,
    const Targ rho,
    Targ c2 = 1.0 / 10)
{
    using namespace details;
    const uint64_t s = order_of_accuracy_differencial_method.at(method);

    const uint64_t initial_steps = std::ceil(1.0 / GetInitialStepSize(f, x0, y0, rho, x_target, s));
    const Targ h_min = (x_target - x0) * 1e-13;
    const Targ h_max = (x_target - x0) * 0.1;

    Targ x_current = x0;
    Tresult y_current = y0;

    Targ h = (x_target - x0) / initial_steps;

    std::vector<std::pair<Targ, Tresult>> result;
    result.push_back({ x_current, y_current });

    while (x_current < x_target) {
        // We limit the step so as not to exceed the target point
        Targ current_step = std::min(h, x_target - x_current);

        auto sol_h = SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method>(f, x_current, y_current, x_current + current_step, 1, c2);
        Tresult y_h = sol_h.back().second;

        auto sol_half = SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method>(f, x_current, y_current, x_current + current_step, 2, c2);
        Tresult y_half = sol_half.back().second;

        double error = RungeError(y_h, y_half, s);

        //         |eps / 2^(s+1)      |eps           |eps * 2^s
        // ------------------------------------------------------------> oX
        //   sc4      sc2                   sc2         sc1
        if (error > rho * std::pow(2, s)) {
            // Scenario 1: the error is too large => reduce the step
            h = std::max(h / 2, h_min);
            continue; // repeat the calculation with a new step
        }
        else if (rho < error && error <= rho * std::pow(2, s)) {
            // Scenario 2: we make a more accurate decision (y_half), reduce the next step
            x_current += current_step;
            y_current = y_half;
            h = std::max(h / 2, h_min);
        }
        else if (rho / std::pow(2, s + 1) <= error && error <= rho) {
            // Scenario 3: we make a decision, leave the step unchanged
            x_current += current_step;
            y_current = y_h;
            // h remains the same
        }
        else if (error < rho / std::pow(2, s + 1)) {
            // Scenario 4: the error is small - we make a decision, increase the step
            x_current += current_step;
            y_current = y_h;
            h = std::min(h * 2, h_max);
        }

        if (x_current <= x_target) {
            result.push_back({ x_current, y_current });
        }

        if (h <= h_min) {
            break;
        }
    }

    return result;
}

//Numerical solution of the Cauchy problem solved with respect to dy/dx, with an error_param scalar
template<DifferencialMethod method, ErrorMethod error_method, typename Targ, typename Tresult>
std::vector<std::pair<Targ, Tresult>> SolveODE(
    std::function<Tresult(Targ, Tresult)> f,
    const Targ x0,
    const Tresult y0,
    const Targ x_target,
    const auto& error_param,  
    Targ c2 = 1.0 / 10)
{
    if constexpr (error_method == ErrorMethod::RungeGlobal) {
        
        return SolveASystemOfOrdinaryDifferentialEquationsEqualStepsWithATotalError<method>(
            f, x0, y0, x_target, static_cast<Targ>(error_param), c2);
    }
    else if constexpr (error_method == ErrorMethod::RungeLocal) {
        return SolveASystemOfOrdinaryDifferentialEquationsDifferentStepsWithALocalError<method>(
            f, x0, y0, x_target, static_cast<Targ>(error_param), c2);
    }
}

