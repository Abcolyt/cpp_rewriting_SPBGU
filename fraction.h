#pragma once
#include<iostream>

template <typename T>class fraction;
template <typename T>std::ostream& operator<<(std::ostream& os, const fraction<T>& fraction);
template <typename T>std::istream& operator>>(std::istream& is, fraction<T>& fraction);

template <typename T>class fraction {
private:
    T numerator;
    T denominator;

public:
    fraction() :fraction(0, 1) {};
    fraction(T num):fraction(num,1){}//поправить
    fraction(T num, T denom)
    {
        denominator = denom;
        numerator = num;
        reduce();
        //std::cout<<"\nreduce:"<<reduce()<<"\n";
    }

    fraction<T>& reduce();

    T findGCD(T& a, T& b)const;

    friend std::ostream& operator<<<>(std::ostream& os, const fraction<T>& fraction);

    friend std::istream& operator>><>(std::istream& is, fraction<T>& fraction);

    fraction operator+(const fraction& other)const;
    fraction operator*(const fraction& other)const;
    fraction operator/(const fraction & other)const;
    fraction operator-(const fraction& other)const;

    fraction operator=(const fraction other);
    fraction operator=(const T other)const;
    bool operator==(const fraction& other)const;
    bool operator==(const T other) const;//it only works for comparison with 1 and 0
    bool operator!=(const fraction& other) const;
    
    T getNumerator() { return numerator; }
    T getDenominator() { return denominator; }
    
};

#include "fraction.cpp"