#pragma once
#include<iostream>
#include <sstream>
#include <complex.h>
#include <vector>

template<typename T>
std::vector<T> convert_pairs_to_vector(
	const std::vector<std::pair<T, T>>& input,
	bool take_first_element
);
//the degree of detail of the output
enum class output_mode
{
	FULL,
	ABBREVIATED,
	SHORT
};
template<typename P> class polynomial;
namespace polynomialfunctions {
	template<typename P>inline polynomial<P> derivate(const polynomial<P>& polynom);
	template<typename P>inline P f_polyn_x0_(const polynomial<P>& polynom, const P& x0);
	template<typename P>P abs(P arg);
	template<typename P>P max(P arg1,P arg2);
	template<typename P>std::pair<P, int> solve_tangents(polynomial<P> plnm, P x0, double e, uint64_t max_iter_number = 1'000'000);
	template<typename P>std::vector<std::pair<P, int>> plnm_roots(polynomial<P> plnm, P x0 = FLT_EPSILON);
	template<typename P>polynomial<P> filter_large_epsilon(polynomial<P>, P eps);

}
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
	uint64_t get_deg()const { return deg; };
	//set the degree of the polynomial
	void set_deg(uint64_t newdeg) {
		
		polynomial<P> newpol;
		newpol.newsize(newdeg);
		for (uint64_t i = 0; i < std::min(deg,newdeg); i++)
			newpol[i] = ptr[i];
		/*for (uint64_t i = std::min(deg, newdeg); i < newdeg; i++)
			newpol[i] = 0;*/
		(*this) = newpol;
	};


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
	
	//NEW start
	//binary polynomial addition with an element from the field
	polynomial<P> operator+(const P& other)const {
		polynomial<P> ans(*this + static_cast<polynomial<P>>(other));
		return ans;
	}
	friend polynomial<P> operator+(const P& scalar, const polynomial<P>& poly) { return poly + scalar; }
	polynomial<P> operator-(const P& other)const {
		polynomial<P> ans(*this - static_cast<polynomial<P>>(other));
		return ans;
	}
	friend polynomial<P> operator-(const P& scalar, const polynomial<P>& poly) { return (poly - scalar)* ( - 1); }
	//NEW end


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
	P& operator[](uint64_t index)const;
	//equalization operator(low coefficient=other )
	polynomial<P>& operator=(const P other);
	//the equalization operator
	polynomial<P>& operator=(const polynomial<P>& other);
	//reduction of the polynomial (due to the higher zero coefficients)
	polynomial cutbag()const;
	//reallocate memory with zeros (old data is not saved)
	void newsize(uint64_t size);
	//set the output mode via a variable of type: uint64_t newmode
	polynomial<P> output_mode_set(uint64_t newmode);
	//set the output mode via a variable of type: uint64_t newmode 
	polynomial<P>& output_mode_set(output_mode new_outm_E);
	//set the output mode via a variable of type: std::string newmode
	polynomial<P> output_mode_set(std::string newmode);
	
	//NEW
	//return Pol(x0)
	P operator()(const P& x) const {
		return polynomialfunctions::f_polyn_x0_(*this, x);
	}
	void the_root_constructor(const std::vector<P>& array_xy);
	polynomial<P> get_first_derrivate()const { return polynomialfunctions::derivate(*this); }
	P maximum_abs(P a, P b)const;
	P maximum()const;

	std::vector<std::pair<P, int>> plnm_roots( P x0 = FLT_EPSILON) const {
		return polynomialfunctions::plnm_roots(*this);
	};
};
template<typename P> output_mode polynomial<P>::default_output_mode = output_mode::SHORT;

#include "../_cpp_realisation_file/polynomial.cpp"



template<typename P>class counting_polynomial :public polynomial<P>
{
public:
	P operator()(const P& x) const {
		return polynomialfunctions::f_polyn_x0_(*this, x);
	}

	//ARITHMETIC OPERATORS
	
	//binary polynomial addition with an element from the field
	counting_polynomial<P> operator+(const P& other)const {
		*this = *this[0] + other;
		return *this;
	}
private:

};

