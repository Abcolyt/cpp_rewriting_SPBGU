#pragma once
#include"fraction.h"
template<typename T>fraction<T>& fraction<T>::reduce()
{
    //if (denominator == 0)return fraction < T>(1);
    //
    //std::max(denominator, numerator);
    //std::min(denominator, numerator);
    T M = std::max(denominator, numerator);
    T m = std::min(denominator, numerator);
    T gcd = findGCD(M,m);
    std::cout << "\ngcd(out7):" << gcd << "\n";
    //std::cout << "\ngcd:" << gcd[0] << "\n";
    std::cout << "\n(numerator/gcd):" << (numerator/gcd) << "\n";
    std::cout << "\n(denominator/ gcd):" << (denominator / gcd) << "\n";
    numerator = (numerator/gcd);
   // std::cout << "\nnumerator:" << numerator << "\n";
    this->denominator = (denominator/ gcd);
    std::cout << "\ndenominator:" << denominator << "\n";
    //std::cout << "\n(*this):" << (*this)<< "\n";
    return *this;
}

template<typename T>T fraction<T>::findGCD(T& a, T& b)const
{
    
    if (b == 0) {
        std::cout << "\ngcd:a:" << a << " b: " << b << " \n";
        return a;
    }
    else {
        T c = (a % b);
        return findGCD(b, c);
    }
}
template<typename T>std::ostream& operator<<(std::ostream& os, const fraction<T>& fraction)
{
    os << fraction.numerator << " / (" << fraction.denominator << ")";
    return os;
}

template<typename T>std::istream& operator>>(std::istream& is, fraction<T>& fraction)
{
    //char slash;
    is >> fraction.numerator; std::cout << "/\n";  is >> fraction.denominator;
    //(fraction).reduce();
    //std::cout << fraction<<"\n";
    std::cout << "\nall\n";
    return is;
}

template<typename T>fraction<T> fraction<T>::operator+(const fraction<T>& other)const
{
    fraction<T> result;
    result.numerator = numerator * other.denominator + other.numerator * denominator;
    result.denominator = denominator * other.denominator;
    return result.reduce();
}

template<typename T>fraction<T> fraction<T>::operator*(const fraction<T>& other)const
{
    fraction<T> result;
    result.numerator = numerator * other.numerator;
    result.denominator = denominator * other.denominator;
    return result.reduce();
}

template<typename T>fraction<T> fraction<T>::operator/(const fraction<T>& other) const 
{
    fraction<T> result;
    result.numerator = this->numerator * other.denominator;
    result.denominator = this->denominator * other.numerator;

    return result.reduce();
}

template<typename T>fraction<T> fraction<T>::operator-(const fraction<T>& other) const
{
    fraction<T> result;
    result.numerator = (this->numerator * other.denominator) - (other.numerator * this->denominator);
    result.denominator = this->denominator * other.denominator;
    return result.reduce();
}
template<typename T>fraction<T> fraction<T>::operator=(const fraction other)const {
    this->numerator = other.numerator;
    this->denominator = other.denominator;

}

template<typename T>fraction<T> fraction<T>::operator=(const T other)const {
    if (T==0)
    {
        numerator = 0;
        denominator = 1;
    }
    else if (T==1)
    {
        numerator = 1;
        denominator = 1;
    }
    else {
       throw std::invalid_argument("the ability to convert to an ordinary fraction is not available");
    }
    return *this;
}

template<typename T>bool fraction<T>::operator==(const fraction& other) const {
    return (this->numerator * other.denominator) == (other.numerator * this->denominator);
}

template<typename T>bool fraction<T>::operator!=(const fraction<T>& other) const {
    return !(*this == other);
}
template<typename T>bool fraction<T>::operator==(const T other) const {
    if (other == 0 and numerator == 0) { 
        return true; 
    }
    else {
        if (fraction<T>(1) ==(*this))
        {
            return true;
        }
    
    }
    return false;

}