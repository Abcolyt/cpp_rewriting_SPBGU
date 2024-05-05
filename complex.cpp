#include "complex.h"

Complex::Complex(double r, double i) : R(r), I(i) {}
Complex::Complex(double r) : Complex(r, 0) {}
Complex::Complex() : Complex(0, 0) {}
Complex::Complex(const Complex& other) : R(other.R), I(other.I) {}

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

Complex Complex::conjugate() {
    return Complex(R, -I);
}

Complex Complex::operator+(const Complex& other) {
    return Complex(R + other.R, I + other.I);
}

Complex Complex::operator-(const Complex& other) {
    return Complex(R - other.R, I - other.I);
}

Complex Complex::operator*(const Complex& other) {
    return Complex(R * other.R - I * other.I, I * other.R + R * other.I);
}

Complex Complex::operator/(const Complex& other) {
    double divisor = other.R * other.R + other.I * other.I;
    return Complex((R * other.R + I * other.I) / divisor, (I * other.R - R * other.I) / divisor);
}

Complex& Complex::operator=(const Complex& other) {
    R = other.R;
    I = other.I;
    return *this;
}

bool Complex::operator==(const Complex& other) {
    return (R == other.R) && (I == other.I);
}

bool Complex::operator!=(const Complex& other) {
    return (R != other.R) || (I != other.I);
}

std::ostream& operator<<(std::ostream& stream, const Complex& number) {
    stream << "Real(" << number.R << ") Imaginary(" << number.I << ")";
    return stream;
}

std::istream& operator>>(std::istream& instream, Complex& number) {
    double R, I;
    instream >> R >> I;
    number.set_RI(R, I);

    return instream;
}