#include <sstream>
#pragma once
enum class output_mode
{
    FULL,
    ABBREVIATED,
    SHORT
};
class Complex {
private:
    double R;
    double I;
public:
    // I/O OPERATIONS
    
    //default output mode
    static output_mode default_output_mode;
    //class realization output mode
    output_mode outm_E = default_output_mode;
    friend std::ostream& operator<<(std::ostream& stream, const Complex& number);
    friend std::istream& operator>>(std::istream& instream, Complex& number);


    //CONSTRUCTORS\DESTRUCTORS

    Complex(double r, double i);
    Complex(double r);
    Complex();
    Complex(const Complex& other);



    //ARITHMETIC OPERATORS

    Complex operator+(const Complex& other);
    Complex operator-(const Complex& other);
    Complex operator*(const Complex& other);
    Complex operator/(const Complex& other);


    //SPECIAL METHODS

    double real_Get();
    double imag_Get();
    void real_Set(double r);
    void imag_Set(double i);
    void set_RI(double real, double imaginary);
    //for example: 2+3i and 2-3i 
    Complex conjugate();
    Complex& operator=(const Complex& other);


    //COMPARISON OPERATORS

    bool operator==(const Complex& other);
    bool operator!=(const Complex& other);
};
output_mode Complex::default_output_mode = output_mode::FULL;
