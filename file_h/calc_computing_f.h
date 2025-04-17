#pragma once
#include "../file_h/complex.h"
#include "../file_h/fraction.h"
#include "../file_h/polynomial.h"
#include "../file_h/matrix.h"

namespace calc_computing_f
{
    void complex_calc() {

        static char L = 'E';
        while (1) {


            Complex a, b, c;
            std::cout << "Enter first Number \n";
            std::cin >> a;
            std::cout << "Enter second Number \n";
            std::cin >> b;

            std::cout << "Enter Action('-','+','/','*','E') :";

            std::cin >> L;
            if (L == 'E' || L == 'e')
            {
                break;
            }
            std::cout << std::endl;
            switch (L)
            {
            case '+':c = a + b; break;
            case '-':c = a - b; break;
            case '*':c = a * b; break;
            case '/':c = a / b; break;
            default:throw std::invalid_argument("unknown symbol");
            }
            std::cout << c << "\n";
        }


    }
    template<typename T>void matrix_calc() {

        static char L = 'E';
        while (1) {


            matrix<T> a, b, ans_out;
            T ans_out_T;
            std::cout << "Enter first Number \n";

            std::cin >> a;
            std::cout << "\n" << a << "\n";
            std::cout << "Enter second Number \n";
            std::cin >> b;
            std::cout << "\n" << b << "\n";
            std::cout << "Enter Action('-','+','/','*','r','E','d') :\n";
            std::cin >> L;
            if (L == 'E' || L == 'e')
            {
                break;
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
            default:throw std::invalid_argument("unknown symbol");
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
    void main_calc_menu() {
        std::cout << "Enter calc type \n";

        std::cout << "mfpint-matrix_calc<fraction<polynomial<int>>>()\n";
        std::cout << "mfpdouble-matrix_calc<fraction<polynomial<double>>>()\n";
        std::cout << "mpint-matrix_calc<polynomial<int>>()\n";
        std::cout << "mpdouble-matrix_calc<polynomial<double>>()\n";
        std::cout << "mdouble-matrix_calc<double>()\n";
        std::cout << "mint-matrix_calc<int>()\n";
        std::cout << "cmpl-complex_calc()\n";
        std::cout << "\ncalc type:";

        std::string type;
        std::cin >> type;

        if (type == "complex_calc()" || type == "cmpl")complex_calc();
        else if (type == "matrix_calc<fraction<polynomial<int>>>()" || type == "mfpint" || type == "matfrpolint") matrix_calc<fraction<polynomial<int>>>();
        else if (type == "matrix_calc<fraction<polynomial<double>>>()" || type == "mfpdouble" || type == "matfrpoldouble") matrix_calc<fraction<polynomial<double>>>();
        else if (type == "matrix_calc<polynomial<int>>()" || type == "mpint" || type == "matpolint") matrix_calc<polynomial<int>>();
        else if (type == "matrix_calc<polynomial<double>>()" || type == "mpdouble" || type == "matpoldouble") matrix_calc<polynomial<double>>();
        else if (type == "matrix_calc<double>()" || type == "mdouble" || type == "matdouble") matrix_calc<double>();
        else if (type == "matrix_calc<int>()" || type == "mint" || type == "matint") matrix_calc<int>();
        else throw std::invalid_argument("unknown calculator type or an error in the name");
    }
    //global computing function 

    void calc_global() {
        try
        {
            std::string type;
        restart:
            main_calc_menu();
            std::cout << "Do you want to get out?(yes/no)\n:";
            std::cin >> type;
            if (type == "yes" || type == "y" || !(type == "no" || type == "n"))std::exit(0);
            if (type == "no" || type == "n")goto restart;

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
