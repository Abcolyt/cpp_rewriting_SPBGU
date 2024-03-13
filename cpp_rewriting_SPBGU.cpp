
#include "stdafx.h"
#include <iostream>
#include <cassert>

class Complex {
private:
	double R;
	double I;
public:
	Complex(double r, double i) : R(r), I(i) {}
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
	// "=="
	a = Complex(1, 1);
	b = Complex(10, 10);
	c = Complex(10, 10);
	assert((a = b) == c);
	//"="
	a = Complex(1, 1);
	b = Complex(10, 15);
	c = Complex(10, 5);
	assert((a = b) != c);
}

int main() {
	TESTs();
	return 0;
}
