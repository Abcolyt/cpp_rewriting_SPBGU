
#include <iostream>
#include <cassert>
#include <sstream>
#pragma warning(disable : 4996)
class Complex {
private:
    double R;
    double I;
public:
    friend std::ostream& operator<<(std::ostream& stream, const Complex& number);
    Complex(double r, double i) : R(r), I(i) {}
    Complex(double r) : Complex(r, 0) {}
    Complex() : Complex(0, 0) {}
    Complex(const Complex& other) : R(other.R), I(other.I) {}
    
    double real_Get();
    double imag_Get();
    void real_Set(double r);
    void imag_Set(double i);
    
    void set_RI(double a, double b);
 
    Complex conjugate() { return Complex(R, -I); }//gives the conjugate complex number
 
    Complex operator+(const Complex& other) {
        return Complex(R + other.R, I + other.I);
    }
 
    Complex operator-(const Complex& other) {
        return Complex(R - other.R, I - other.I);
    }
 
    Complex operator*(const Complex& other) {
        return Complex(R * other.R - I * other.I, I * other.R + R * other.I);
    }
 
    Complex operator/(const Complex& other) {
        double divisor = other.R * other.R + other.I * other.I;
        
        return Complex((this->R * other.R + this->I * other.I) / divisor,(this->I * other.R - this->R * other.I) / divisor);
    }
 
    Complex& operator = (const Complex& other)
    {
        this->R = other.R;
        this->I = other.I;
        return *this;
    }
 
    bool operator==(const Complex& other) {
        return (R == other.R) && (I == other.I);
    }
    bool operator!=(const Complex& other) {
        return (R != other.R)|(I != other.I);
    }
 
};
 
std::ostream& operator<<(std::ostream& stream, const Complex& number) {
    stream << "Real(:" << number.R << ")Imaginary(:" << number.I << ")";
    return stream;
 }
 
std::istream& operator >> (std::istream& instream, Complex& number)
{
    double R,I;
    instream >> R >> I;
    number.set_RI(R,I);
 
    return instream;
}
 
 
double Complex::real_Get() {
    return R;
}
 
double Complex::imag_Get() {
    return I;
}
 
void Complex::real_Set(double r) {
    R = r;
}
 
void Complex::imag_Set(double i) {
    I = i;
}
 
void Complex::set_RI(double a, double b) {
    R = a;
    I = b;
}
 
void TESTs() {
    Complex a, b, c,d;
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
    a = Complex(6,8);
    b = Complex(5,15);
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
    a=Complex(3, 4);
    std::stringstream so;
    so << a;
    assert(so.str() == "Real(:3)Imaginary(:4)");
 
    // ">>" 
    std::stringstream si("Real(:3)Imaginary(:4)");
    si >> a; 
    assert(a.real_Get() == 3 && a.imag_Get() == 4);
    
}
void complex_calc() {
    
    static char L = 'E';
    while (1) {
        double R, I;
        
        Complex a,b,c;
        std::cout << "Enter first Number \n";
        std::cin >> a;
        std::cout << "Enter second Number \n";
        std::cin >> b;
 
        std::cout << "Enter Action('-','+','/','*','E','T') :";
        
        std::cin >> L;
        //getchar();
        if ((L == 'E') ||(L=='e'))
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

class Matrix
{
public:
    Matrix();
    ~Matrix();

private:
    unsigned int line, column;
};

Matrix::Matrix()
{
}

Matrix::~Matrix()
{
}
void matrix_calc() {

    static char L = 'E';
    while (1) {
        double R, I;

        Matrix a, b, c;
        std::cout << "Enter first Number \n";
        std::cin >> a;
        std::cout << "Enter second Number \n";
        std::cin >> b;

        std::cout << "Enter Action('-','+','/','*','E','T') :";

        std::cin >> L;
        //getchar();
        if ((L == 'E') || (L == 'e'))
        {
            abort();
        }
        std::cout << std::endl;
        switch (L)
        {
        case '+':c = a + b; break;
        case '-':c = a - b; break;
        case '*':c = a * b; break;
        case 'x':c = a[b]; break;//scalar multiplication
        case '/':c = a / b; break;
        case 'T':TESTs(); break;
        case 't':TESTs(); break;
        default:abort();
        }
        std::cout << c << "\n";
    }


}









int main() {
    complex_calc();
    system("pause");
    return 0;
}
