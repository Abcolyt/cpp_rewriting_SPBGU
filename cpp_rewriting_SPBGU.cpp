#include <iostream>
#include <sstream>
#include "complex.h"
#include "matrix.h"
#include "polynomial.h"
#include <string>
#include "fraction.h"

//I dont recommend using with fraction<polynomial<double>> /fraction<polynomial<float>>,because of the representation, the coefficients will not be integers, for example (-36.125 -144.5x -289x^2 -433.5x^3 -397.375x^4 -144.5x^5) / (0 72.25x 180.625x^2 -108.375x^3 -614.125x^4 -541.875x^5 -144.5x^6) 
template<typename T>void matrix_calc() {

    static char L = 'E';
    while (1) {


        matrix<T> a, b, ans_out;
        T ans_out_T;
        std::cout << "Enter first Number \n";
        
        std::cin >> a;
        std::cout <<"\n" << a<<"\n";
        std::cout << "Enter second Number \n";
        std::cin >> b;
        std::cout << "Enter Action('-','+','/','*','r') :\n";
        std::cin >> L;
        if (L == 'E' || L == 'e')
        {
            abort();
        }
        std::cout << std::endl;
        switch (L)
        {
        case '+':ans_out = a + b; break;
        case '-':ans_out = a - b; break;
        case '*':ans_out = a * b; break;
        case '/':ans_out = a / b; break;
        case 'd':ans_out_T=(a.determinant()); break;
        case 'r':ans_out = a.inverse_M(); break;
        default:abort();
        }
        if (L == 'd') {
            std::cout << ans_out_T << "\n";
        }
        else
        {
            std::cout << ans_out << "\n";
        }

       
    }


}

int main() {

    matrix_calc<fraction<polynomial<int>>>();
    
    system("pause");
    return 0;
}
