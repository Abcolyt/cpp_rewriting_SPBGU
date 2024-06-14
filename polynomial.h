#pragma once
#include<iostream>
#include <sstream>
enum class output_mode
{
	FULL,
	ABBREVIATED,
	SHORT
};
template<typename P> class polynomial;
template<typename P> std::ostream& operator<<(std::ostream& out, const polynomial<P>& plnm);
template<typename P> std::istream& operator>> (std::istream& in, polynomial<P>& plnm);

template<typename P> class polynomial
{
private:

public:

	output_mode outm_E;



	P* ptr;
	uint64_t deg;
	polynomial() :polynomial(0) {};
	polynomial(P number);
	~polynomial();
	polynomial(const polynomial<P>& other);
	polynomial<P> operator+(const polynomial<P>& other)const;//addition of polynomials
	polynomial<P> operator-(const polynomial<P>& other)const;
	polynomial<P> operator/(const polynomial<P>& other)const;//binary division of a polynomial by a polynomial
	polynomial<P> operator%(const polynomial<P>& other)const;//the remainder of the division of a polynomial by a polynomial
	polynomial<P> operator*(const polynomial<P>& other)const;//binary polynomial multiplication 
	polynomial<P> operator*(const P& other)const; //binary polynomial multiplication by an element from the field

	friend std::ostream& operator<<<>(std::ostream& out, const polynomial<P>& plnm); // overloading the output operator
	friend std::istream& operator>><>(std::istream& in, polynomial<P>& plnm); // overloading the input operator

	polynomial<P> operator>>(const uint64_t power)const;//coefficient shift (increase). or (plnm*x^n)
	polynomial<P> operator<<(const uint64_t power)const;//coefficient shift(decrease). or the whole part of (plnm* x^(-n))

	P& operator[](uint64_t index);//the access operator to the polynomial coefficient
	polynomial<P>& operator=(const P other);
	polynomial<P>& operator=(const polynomial<P>& other);

	bool operator==(const polynomial<P>& other)const;// 
	bool operator==(int64_t zero)const;//isZero()?(all values==0 so polynomial have zero coefficent in any position)
	bool operator!=(int64_t zero)const { return not((*this) == zero); }
	bool operator!=(const polynomial<P>& other)const;

    bool operator>(const polynomial<P>& other)const;
	bool operator<(const polynomial<P>& other)const;
	

	polynomial cutbag()const;
	void newsize(uint64_t size);

};

#include "polynomial.cpp"
