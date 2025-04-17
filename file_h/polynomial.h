#pragma once
#include<iostream>
#include <sstream>
#include <complex.h>
#include <vector>


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
	template<typename P>std::pair<P, int> solve_tangents(polynomial<P> plnm, P x0, double e, uint64_t max_iter_number = 1'000'000);
	template<typename P>std::vector<std::pair<P, int>> plnm_roots(polynomial<P> plnm, P x0 = FLT_EPSILON);
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
	void output_mode_set(uint64_t newmode);
	//set the output mode via a variable of type: uint64_t newmode 
	polynomial<P>& output_mode_set(output_mode new_outm_E);
	//set the output mode via a variable of type: std::string newmode
	void output_mode_set(std::string newmode);
	
	//NEW
	//return Pol(x0)
	P operator()(const P& x) const {
		return polynomialfunctions::f_polyn_x0_(*this, x);
	}
	void the_root_constructor(const std::vector<P>& array_xy);
	polynomial<P> get_first_derrivate()const { return polynomialfunctions::derivate(*this); }
	P maximum(P a, P b)const;
	P maximum()const;
};
template<typename P> output_mode polynomial<P>::default_output_mode = output_mode::SHORT;

#include "../_cpp_realisation_file/polynomial.cpp"
namespace polynomialfunctions {

	template<typename P>inline polynomial<P> derivate(const polynomial<P>& polynom) {
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
	template<typename P>inline P f_polyn_x0_(const polynomial<P>& polynom,const P& x0) {
		P ans(0);
		/*for (uint64_t i = 0; i < polynom.get_deg(); i++) {
			std::cout << '\n' << polynom[i]   << std::endl;
		}*/
		for (uint64_t i = polynom.get_deg() ; i > 0; --i) {
			//std::cout << "\n ans:" << ans<<"i=" << i - 1 << std::endl;
			ans = polynom[i - 1] + (ans)*x0;
		}
		//std::cout << "\n ans:" << ans;
		return ans;
	}

	template<typename P>P abs(P arg) {
		if (arg < 0) { arg = arg * (-1); }
		return arg;
	}

	template<typename P>std::pair<P, int> solve_tangents(polynomial<P> plnm, P x0, double e, uint64_t max_iter_number ) {
		int iterations = 0;
		P x_old = 5, x_new = x0;
		auto derivate = polynomialfunctions::derivate(plnm);
		while (abs(x_new - x_old) > e and iterations < max_iter_number) {
			auto T = abs(x_new - x_old);
			x_old = x_new;
			x_new = x_old + (polynomialfunctions::f_polyn_x0_<P>(plnm, x_old) / (polynomialfunctions::f_polyn_x0_<P>(derivate, x_old))) * (-1);
			iterations++;
		}
		return { x_new, iterations };
	}

	template<typename P>std::vector<std::pair<P, int>> plnm_roots(polynomial<P> plnm, P x0)
	{
		std::vector<std::pair<P, int>> ans_roots;
		polynomial<P> b;
		size_t deg = plnm.get_deg() - 1;
		//std::cout << plnm << "  deg" << plnm.get_deg() << "\n";

		for (size_t i = 0; i < deg; i++)
		{
			//std::cout << "plnm:" << plnm << "\nnew plnm" << b << "\nnew iter" << "'" << i << "' \n";
			auto ans = solve_tangents<P>(plnm, x0, LDBL_EPSILON);
			//std::cout << "ans" << solve_tangents<P>(plnm, x0, LDBL_EPSILON).first << "\n";
			//std::cout <<"root["<<i<<"]=" << ans.first << '\n';
			b.newsize(2);
			b[1] = 1;
			b[0] = (ans.first) * (-1);

			plnm = (plnm / b);
			//std::cout << "plnm after:" << plnm<<"\n";
			ans_roots.push_back(ans);
		}
		return ans_roots;
	}

}


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

