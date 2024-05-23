#pragma once
#include"fraction.h"
template<typename T>fraction<T> fraction<T>::reduce()
{
    //if (denominator == 0)return fraction < T>(1);
    T& gcd = findGCD(denominator, numerator);
    std::cout << "\ngcd(out7):" << gcd << "\n";
    //std::cout << "\ngcd:" << gcd[0] << "\n";
    numerator = (numerator/gcd);
    denominator = (denominator/ gcd);
    std::cout << "\ngcd:" << (*this)<< "\n";
    return *this;
}

template<typename T>T& fraction<T>::findGCD(T& a, T& b)const
{
    if (b == 0) {
        //std::cout << "\ngcd:a:" << a << " b: " << b << " \n";
        return a;
    }
    else {
        return findGCD(b, (a % b));
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
    return is;
}

template<typename T>fraction<T> fraction<T>::operator+(const fraction<T>& other)const
{
    T newNumerator = numerator * other.denominator + other.numerator * denominator;
    T newDenominator = denominator * other.denominator;
    return fraction<T>(newNumerator, newDenominator);
}

template<typename T>fraction<T> fraction<T>::operator*(const fraction<T>& other)const
{
    T newNumerator = numerator * other.numerator;
    T newDenominator = denominator * other.denominator;
    return fraction<T>(newNumerator, newDenominator);
}

template<typename T>fraction<T> fraction<T>::operator/(const fraction<T>& other) const {
    
    T newNumerator = this->numerator * other.denominator;
    T newDenominator = this->denominator * other.numerator;
    return ((fraction<T>(newNumerator, newDenominator)).reduce());
}

template<typename T>fraction<T> fraction<T>::operator-(const fraction<T>& other) const {
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