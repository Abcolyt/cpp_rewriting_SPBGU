#pragma once
#include<iostream>
#include <sstream>
#include <complex.h>
//the degree of detail of the output
extern enum class output_mode
{
	FULL,
	ABBREVIATED,
	SHORT,
	NO
};
template<typename P> class polynomial;
template<typename P> std::ostream& operator<<(std::ostream& out, const polynomial<P>& plnm);
template<typename P> std::istream& operator>>(std::istream& in, polynomial<P>& plnm);
template<typename P> class polynomial
{
private:
	//degree of the polynomial
	uint64_t deg;
public:
	//an array of coefficients
	P* ptr;
	//default output mode
	static output_mode default_output_mode;
	//class realization output mode
	output_mode outm_E;


	//CONSTRUCTORS\DESTRUCTORS
	polynomial(std::string constructor_mode, int polynomial_degree) {

	}
	polynomial() :polynomial(0) {};
	polynomial(P number, output_mode mode) :polynomial(number) { outm_E = mode; }
	polynomial(P number);
	polynomial(const polynomial<P>& other, output_mode mode) :polynomial(other) { outm_E = mode; }
	polynomial(const polynomial<P>& other);
	~polynomial();


	//DATA ACCESS
	
    //get the degree of the polynomial
	uint64_t get_deg() { return deg; };
	//set the degree of the polynomial
	void set_deg(uint64_t newdeg) { deg = newdeg; };


	//ARITHMETIC OPERATORS
	 
	//addition of polynomials
	polynomial<P> operator+(const polynomial<P>& other)const;
	//subtruction of polynomials
	polynomial<P> operator-(const polynomial<P>& other)const;
	//binary division of a polynomial by a polynomial
	polynomial<P> operator/(const polynomial<P>& other)const;
	//the remainder of the division of a polynomial by a polynomial
	polynomial<P> operator%(const polynomial<P>& other)const;
	//binary polynomial multiplication 
	polynomial<P> operator*(const polynomial<P>& other)const;
	//binary polynomial multiplication by an element from the field
	polynomial<P> operator*(const P& other)const;


    // I/O OPERATIONS

	//the output operator (with the degree of detail specified in the outm_E field)
	friend std::ostream& operator<<<>(std::ostream& out, const polynomial<P>& plnm); 
	//the input operator
	friend std::istream& operator>><>(std::istream& in, polynomial<P>& plnm);


	//CHARACTER SHIFTS

	//coefficient shift (increase). or (plnm*x^n)
	polynomial<P> operator>>(const uint64_t power)const;
	//coefficient shift(decrease). or the whole part of (plnm* x^(-n))
	polynomial<P> operator<<(const uint64_t power)const;


	//COMPARISON OPERATORS
	
	//is it true if the polynomials are equal
	bool operator==(const polynomial<P>& other)const;
	//is it false if the polynomials are equal
	bool operator!=(const polynomial<P>& other)const;
	//is it true if this->deg > other.deg;
	bool operator>(const polynomial<P>& other)const;
	//is it true if this->deg < other.deg;
	bool operator<(const polynomial<P>& other)const;


	//COMPARISON OPERATORS WITH ZERO

	//isZero()?(all values==0 so polynomial have zero coefficent in any position)
	bool operator==(int64_t zero)const;
	//!isZero()?(all values==0 so polynomial have zero coefficent in any position)
	bool operator!=(int64_t zero)const { return not((*this) == zero); }
	

	//SPECIAL METHODS

    //the access operator to the polynomial coefficient
	P& operator[](uint64_t index);
	//equalization operator(low coefficient=other )
	polynomial<P>& operator=(const P other);
	//the equalization operator
	polynomial<P>& operator=(const polynomial<P>& other);
	//reduction of the polynomial (due to the higher zero coefficients)
	polynomial cutbag()const;
	//reallocate memory with zeros (old data is not saved)
	void newsize(uint64_t size);
	//set the output mode via a variable of type: uint64_t newmode
	void output_mode_set(uint64_t newmode);
	//set the output mode via a variable of type: uint64_t newmode 
	void output_mode_set(output_mode new_outm_E);
	//set the output mode via a variable of type: std::string newmode
	void output_mode_set(std::string newmode);

};
template<typename P> output_mode polynomial<P>::default_output_mode = output_mode::SHORT;

#include "../_cpp_realisation_file/polynomial.cpp"
namespace polynomialfunctions {

	template<typename P>inline polynomial<P> derivate(polynomial<P> polynom) {
		polynomial<P> ans;
		ans.newsize(polynom.get_deg());
		if (polynom.get_deg() > 0) {
			for (uint64_t i = polynom.get_deg() - 1; i >= 1; i--) {
				ans[i - 1] = polynom[i] * i;
			}
			ans = ans.cutbag();
		}
		else {
			ans = 0;
		}
		return ans;
	}
	template<typename P>inline P f_polyn_x0_(polynomial<P> polynom, P x0) {
		P ans(0);
		for (uint64_t i = 0; i < polynom.get_deg(); i++) {
			//std::cout << polynom[i] + x0 << '\n' << std::endl;
		}
		for (int i = polynom.get_deg() - 1; i >= 0; i--) {
			ans = polynom[i] + (ans)*x0;
		}
		return ans;
	}

}