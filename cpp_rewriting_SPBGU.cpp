#include <iostream>
#include <sstream>
#include <string>
#include <fraction.h>
#include <polynomial.h>
#include <complex.h>
#include <matrix.h>


#include <vector>

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

template<typename P>
P abs(P arg) {
    if (arg < 0) { arg = arg * (-1); }
    return arg;
}

template<typename P>std::pair<P, int> solve_tangents(polynomial<P> plnm,P x0, double e) {
    int iterations = 0;
    P x_old = 5, x_new = x0;
    while (abs(x_new - x_old) > e) {
        x_old = x_new;
        x_new = x_old +  (polynomialfunctions::f_polyn_x0_<P>(plnm, x_old) / (polynomialfunctions::f_polyn_x0_<P>(polynomialfunctions::derivate(plnm), x_old)))*(-1);
        iterations++;
    }

    return { x_new, iterations };
}



template<typename P>
std::vector<std::pair<P, int>> plnm_roots(polynomial<P> plnm,P x0) {
    std::vector<std::pair<P,int>> ans_roots;
    polynomial<P> b;
    //1
    auto ans = solve_tangents<P>(plnm,x0, LDBL_EPSILON);
    b.newsize(2);
    b[1] = 1;
    b[0] = (ans.first)*(-1);
    
    plnm = (plnm / b);
    std::cout << "plnm:"<<plnm<<"\nnew plnm"<< b << "\nnew iter" << '1' << '\n';
    ans_roots.push_back(ans);

    //
    
    //2
    ans = solve_tangents<P>(plnm, x0, LDBL_EPSILON);
    b.newsize(2);
    b[1] = 1;
    b[0] = (ans.first) * (-1);
    plnm = (plnm / b);
    std::cout << "plnm:" << plnm << "\nnew plnm" << b << "new iter" << '2' << '\n';
    ans_roots.push_back(ans);

    //

    //3
    ans = solve_tangents<P>(plnm, x0, LDBL_EPSILON);
    b.newsize(2);
    b[1] = 1;
    b[0] = (ans.first) * (-1);
    plnm = (plnm / b);
    std::cout << "plnm:" << plnm << "\nnew plnm" << b << "new iter" << '3' << '\n';
    ans_roots.push_back(ans);

    //
    return ans_roots;
}

void polynomial_test() {
    polynomial<Complex> plnm, b;
    std::cin >> plnm;
    auto ans = plnm_roots<Complex>(plnm,Complex(LDBL_EPSILON,LDBL_EPSILON));
    for (size_t i = 0; i < ans.size(); i++)
    {
        std::cout << "root:(" << ans[i].first << ")\niteration for this root" << ans[i].second << "\n";
    }

    //std::cout<<polynomialfunctions::derivate(a)<<"\n"<<a<<"\n"<< polynomialfunctions::f_polyn_x0_(a, 5.25)<<std::endl;
    /*auto ans = solve_tangents(plnm, LDBL_EPSILON);
    std::cout << "root:" << ans.first << "iter:" << ans.second<<'\n';
    b.newsize(2);
    b[1] = 1;
    b[0] = -ans.first;
    std::cout << b<<'\n' << (plnm / b);*/

}
int main() {
         
    polynomial_test();
    system("pause");
    return 0;
}
