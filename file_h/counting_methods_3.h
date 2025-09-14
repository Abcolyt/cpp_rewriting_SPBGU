#pragma once
#include "../file_h/matrix.h"
#include <vector>
namespace counting_methods_3 {
#if 0
    template<typename T>
	matrix<T> ÑoefficNewtonCotes(std::vector<T> nodes_in_the_integration_gap, T beginning_of_the_integration_interval, T end_of_the_integration_interval) {
        
        matrix<T>A(nodes_in_the_integration_gap.size()), b(nodes_in_the_integration_gap.size(), 1);

        for (int i = 0; i < nodes_in_the_integration_gap.size(); i++) {
            A[0][i] = 1;//nodes_in_the_integration_gap[i];
            b[i][0] = 1 / (static_cast<T>(1 + i));
            for (int j = 1; j < nodes_in_the_integration_gap.size(); j++) {
                A[j][i] = nodes_in_the_integration_gap[i] * A[j - 1][i];
            }
        }
		
		matrix<T> coeff = (end_of_the_integration_interval - beginning_of_the_integration_interval) * matrixfunction::solve_system(A,b);

		return  coeff;
 }
#endif
    


    //get c[j] : Integral(a,b) of function =~= Sum(j=1,n){c[j]*F(x_j)} 
    template<typename T>
    matrix<T> ÑoefficIntegralNewtonCotes(matrix<T> nodes_in_the_integration_gap, T beginning_of_the_integration_interval, T end_of_the_integration_interval) {
        if (nodes_in_the_integration_gap.is_vector()) {
            //std::cout << "nodes_in_the_integration_gap:\n" << nodes_in_the_integration_gap << "\n";
            if (nodes_in_the_integration_gap.getcol() == 1) { nodes_in_the_integration_gap = nodes_in_the_integration_gap.transpose(); }
            std::cout << nodes_in_the_integration_gap <<'\n';
            if (beginning_of_the_integration_interval != 0 | end_of_the_integration_interval != 1) {
                std::cout << (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow()))) << "\n";
                nodes_in_the_integration_gap = (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow())));
                nodes_in_the_integration_gap = nodes_in_the_integration_gap /(end_of_the_integration_interval - beginning_of_the_integration_interval);
            }
            uint64_t size= std::max(nodes_in_the_integration_gap.getrow(), nodes_in_the_integration_gap.getcol());
            matrix<T>A(size), b(size, 1);
            std::cout << nodes_in_the_integration_gap << '\n';
            for (int i = 0; i < size; i++) {
                A[0][i] = 1;//nodes_in_the_integration_gap[i];
                b[i][0] = 1 / (static_cast<T>(1 + i));
                for (int j = 1; j < size; j++) {
                    A[j][i] = nodes_in_the_integration_gap[i][0] * A[j - 1][i];
                }
            }

            matrix<T> coeff = (end_of_the_integration_interval - beginning_of_the_integration_interval) * matrixfunction::solve_system(A, b);

            return coeff;
        }
    }

    //get c[j] : Integral(a,b) of function =~= Sum(j=1,n){(b-a)*c[j]*F(x_j)} 
    template<typename T>
    matrix<T> ÑoefficNewtonCotes(matrix<T> nodes_in_the_integration_gap, T beginning_of_the_integration_interval, T end_of_the_integration_interval) {
        if (nodes_in_the_integration_gap.is_vector()) {
            //std::cout << "nodes_in_the_integration_gap:\n" << nodes_in_the_integration_gap << "\n";
            if (nodes_in_the_integration_gap.getcol() == 1) { nodes_in_the_integration_gap = nodes_in_the_integration_gap.transpose(); }
            std::cout << nodes_in_the_integration_gap << '\n';
            if (beginning_of_the_integration_interval != 0 | end_of_the_integration_interval != 1) {
                std::cout << (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow()))) << "\n";
                nodes_in_the_integration_gap = (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow())));
                nodes_in_the_integration_gap = nodes_in_the_integration_gap / (end_of_the_integration_interval - beginning_of_the_integration_interval);
            }
            uint64_t size = std::max(nodes_in_the_integration_gap.getrow(), nodes_in_the_integration_gap.getcol());
            matrix<T>A(size), b(size, 1);
            std::cout << nodes_in_the_integration_gap << '\n';
            for (int i = 0; i < size; i++) {
                A[0][i] = 1;//nodes_in_the_integration_gap[i];
                b[i][0] = 1 / (static_cast<T>(1 + i));
                for (int j = 1; j < size; j++) {
                    A[j][i] = nodes_in_the_integration_gap[i][0] * A[j - 1][i];
                }
            }

            matrix<T> coeff =  matrixfunction::solve_system(A, b);

            return coeff;
        }
    }


}