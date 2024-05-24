#include <iostream>
#include <assert.h>
#include <sstream>
#include "complex.h"
#include "matrix.h"
#include "polynomial.h"
#include <string>
#include "fraction.h"

void TESTs() {
    Complex a, b, c, d;
    // "+"
    a = Complex(1, 1);
    b = Complex(-1, -1);
    c = Complex(0, 0);
    assert((a + b) == c);

    // "-"
    a = Complex(1, 1);
    b = Complex(-1, -1);
    c = Complex(2, 2);
    assert((a - b) == c);

    // "*"
    a = Complex(1, 1);
    b = Complex(1, 1);
    c = Complex(0, 2);
    assert((a * b) == c);

    // "/"
    a = Complex(1, 1);
    b = Complex(1, 1);
    c = Complex(1, 0);
    assert((a / b) == c);
    a = Complex(6, 8);
    b = Complex(5, 15);
    c = Complex(0.6, -0.2);
    assert((a / b) == c);

    // "==" and "="
    a = Complex(1, 1);
    b = Complex(10, 10);
    c = Complex(10, 10);
    assert((a = b) == c);

    // "!="
    a = Complex(1, 1);
    b = Complex(10, 15);
    c = Complex(10, 5);
    assert((a = b) != c);

    // "<<" 
    a = Complex(3, 4);
    std::stringstream so;
    so << a;
    assert(so.str() == "Real(:3)Imaginary(:4)");

    // ">>" 
    //std::stringstream si("Real(:3)Imaginary(:4)");
    //si >> a; 
    //assert(a.real_Get() == 3 && a.imag_Get() == 4);

}
void complex_calc() {
    
    static char L = 'E';
    while (1) {
   
        
        Complex a,b,c;
        std::cout << "Enter first Number \n";
        std::cin >> a;
        std::cout << "Enter second Number \n";
        std::cin >> b;
 
        std::cout << "Enter Action('-','+','/','*','E','T') :";
        
        std::cin >> L;
        //getchar();
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
        case 'T':TESTs(); break;
        case 't':TESTs(); break;
        default:abort();
        }
        std::cout << c<<"\n";
    }
 
    
}
void backspace(int count)
{
    for (int i = 0; i < count; i++)
        printf_s("\b \b");
}

template<typename T>void matrix_calc() {

    static char L = 'E';
    while (1) {


        matrix<T> a, b, ans_out;
        T ans_out_T;
        std::cout << "Enter first Number \n";
        
        std::cin >> a;
        std::cout <<"\n" << a<<"\n";
         "Enter second Number \n";
        std::cin >> b;
        std::cout << "Enter Action('-','+','/','*','r','T') :\n";

        std::cin >> L;
        //getchar();
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

template<typename P>void polynomial_matrix_calc() {

    static char L = 'E';
    while (1) {


        matrix<polynomial<P>> a, b, ans_out;
        polynomial<P> ans_out_T;
        std::cout << "Enter first Number \n";

        std::cin >> a;
        std::cout << "\n" << a << "\n";
        "Enter second Number \n";
        std::cin >> b;
        std::cout << "Enter Action('-','+','/','*','r','T') :\n";

        std::cin >> L;
        //getchar();
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
        case 'd':ans_out_T = (a.determinant()); break;
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
    fraction<polynomial<int64_t>> a(1,1), b;
    std::cout << "a=" << a << "\nb=" << b << "\n";
    std::cin >> a;
    a.reduce();
    std::cout << "a=" << a << "\n";
    std::cin >> b;
    std::cout << "b=" <<b << "\n";
    std::cout <<"a="<<a<<"b="<<b <<"a / b=" << (a / b)<<"(a == 0)"<< (a == 1);

    //polynomial_matrix_calc<double>();
    //matrix_calc<fraction<polynomial<double>>>();


    system("pause");
    return 0;
}
