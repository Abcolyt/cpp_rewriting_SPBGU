#include <sstream>
#pragma once
class Complex {
private:
    double R;
    double I;
public:
    friend std::ostream& operator<<(std::ostream& stream, const Complex& number);
    friend std::istream& operator>>(std::istream& instream, Complex& number);

    Complex(double r, double i);
    Complex(double r);
    Complex();
    Complex(const Complex& other);

    double real_Get();
    double imag_Get();
    void real_Set(double r);
    void imag_Set(double i);

    void set_RI(double a, double b);

    Complex conjugate();

    Complex operator+(const Complex& other);
    Complex operator-(const Complex& other);
    Complex operator*(const Complex& other);
    Complex operator/(const Complex& other);

    Complex& operator=(const Complex& other);

    bool operator==(const Complex& other);
    bool operator!=(const Complex& other);
//new


    // ���������� ��������� <
    bool operator<(const Complex& other) const {
        return this->magnitude() < other.magnitude();
    }

    // ���������� ��������� < ��� ��������� � �����
    bool operator<(double value) const {
        return this->magnitude() < value;
    }

    // ���������� ��������� >
    bool operator>(const Complex& other) const {
        return this->magnitude() > other.magnitude();
    }

    // ���������� ��������� > ��� ��������� � �����
    bool operator>(double value) const {
        return this->magnitude() > value;
    }

    // ����� ��� ���������� ������ ������������ �����
    double magnitude() const {
        return std::sqrt(R * R + I * I);
    }

};

