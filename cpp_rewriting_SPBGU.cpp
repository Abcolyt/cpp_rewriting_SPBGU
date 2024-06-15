#include <iostream>
#include <assert.h>
#include <sstream>
#include "complex.h"
#include "matrix.h"
#include "polynomial.h"
#include <string>
#include "fraction.h"


void complex_calc() {
    
    static char L = 'E';
    while (1) {
   
        
        Complex a,b,c;
        std::cout << "Enter first Number \n";
        std::cin >> a;
        std::cout << "Enter second Number \n";
        std::cin >> b;
 
        std::cout << "Enter Action('-','+','/','*','E') :";
        
        std::cin >> L;
        if (L == 'E'||L=='e')
        {
            abort();
        }
        std::cout << std::endl;
        switch (L)
        {
        case '+':c=a+b; break;
        case '-':c = a - b; break;
        case '*':c = a * b; break;
        case '/':c = a / b; break;
        default:abort();
        }
        std::cout << c<<"\n";
    }
 
    
}
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
    try
    {
        matrix_calc<fraction<polynomial<int>>>();
    }
    catch (std::exception ex)
    {
        std::cout << "exeption!!What:" << ex.what();
    }
    catch (...)
    {
        std::cout << "unknown error";
    }
    
    system("pause");
    return 0;
}
