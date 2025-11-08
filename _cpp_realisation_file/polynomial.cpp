#pragma once
#include"../file_h/polynomial.h"

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
	template<typename P>inline P f_polyn_x0_(const polynomial<P>& polynom, const P& x0) {
		P ans(0);

		for (uint64_t i = polynom.get_deg(); i > 0; --i) {
			ans = polynom[i - 1] + (ans)*x0;
		}
		return ans;
	}

	template<typename P>P abs(P arg) {
		if (arg < 0) { arg = arg * (-1); }
		return arg;
	}
	template<typename P>P max(P arg1, P arg2) {
		if (arg1 <= arg2) {
			return arg2;
		}
		else
			return arg1;
	}
	template<typename P>std::pair<P, int> solve_tangents(polynomial<P> plnm, P x0, double e, uint64_t max_iter_number) {
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

	template<typename P>
	std::vector<std::pair<P, int>> plnm_roots(polynomial<P> plnm, P x0) {
		polynomial<P> original_plnm = plnm; // Сохраняем исходный полином
		std::vector<std::pair<P, int>> ans_roots;
		polynomial<P> b;
		size_t deg = plnm.get_deg();

		for (size_t i = 0; i < deg - 1; i++) {
			auto ans = solve_tangents<P>(plnm, x0, LDBL_EPSILON);
			auto refined_ans = solve_tangents<P>(original_plnm, ans.first, LDBL_EPSILON);

			ans_roots.push_back(refined_ans);

			b.newsize(2);
			b[1] = 1;
			b[0] = (refined_ans.first)* ( - 1.0);
			plnm = original_plnm; 
			for (auto& root : ans_roots) {
				polynomial<P> divisor; divisor.newsize(2);
				divisor[1] = 1;
				divisor[0] = (root.first) * (-1.0);
				plnm = plnm / divisor; 
			}
		}
		return ans_roots;
	}
	template<typename P>polynomial<P> filter_large_epsilon(polynomial<P> pol, P eps) {
		for (uint64_t i = 0; i < pol.get_deg(); i++)
		{
				if (-eps < pol[i] && pol[i] < eps){
					pol[i] = 0;
			}
		}
		return pol.cutbag();

	}

}

template<typename P>polynomial<P> polynomial<P>::operator>>(const uint64_t power)const
{
	polynomial<P> result;result.newsize(deg + power);
	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i + power] = ptr[i];
		
	}
	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator<<(const uint64_t power)const
{
	polynomial result;
	if(deg >= power)
	 result.newsize(deg - power);
	else result.newsize(0);
	for (uint64_t i = 0; i < result.get_deg(); i++)
	{
		result.ptr[i] = ptr[i + power];
	}
	return result;
}

template<typename P>bool polynomial<P>::operator==(const polynomial<P>& other)const 
{
	if (deg != other.deg)return false;
	else
	{
		for (uint64_t i = 0; i < deg; i++)
		{
			if (ptr[i] != other.ptr[i])return false;
		}
	}
	return true;
}

template<typename P>bool polynomial<P>::operator==(int64_t zero)const
{
	if (zero != 0) {
		return 0;
	}
	for (uint64_t i = 0; i < this->deg; i++)
	{
		if (ptr[i] != 0)
		{
			return 0;
		}
	}
	return 1;


}

template<typename P>bool polynomial<P>::operator!=(const polynomial<P>& other)const 
{
	if (deg != other.deg)return true;
	else
	{
		for (uint64_t i = 0; i < deg; i++)
		{
			if (ptr[i] != other.ptr[i])return true;
		}
	}
	return false;
}

template<typename P>polynomial<P>::polynomial(P number) {
	
	
	outm_E = default_output_mode;
	deg = 1;
	ptr = new P[deg];
	ptr[0] = number;



}

template<typename P>void polynomial<P>::newsize(uint64_t size)
{
	deg = size;
	if (size) {
		delete[] ptr;
		ptr = new P[size];
	}
	else
		ptr = nullptr;
	for (uint64_t i = 0; i < size; i++)
		ptr[i] = 0;
}
template<typename P>polynomial<P>::polynomial(const polynomial<P>& other)
{
	
	outm_E = default_output_mode;
	deg = other.deg; 
	ptr = new P[deg + 1]; 

	for (uint64_t i = 0; i < deg; ++i)
	{
		ptr[i] = other.ptr[i]; 
	}
}
template<typename P>polynomial<P>::~polynomial()
{
	delete[] ptr;
}
template<typename P>polynomial<P>& polynomial<P>::operator=(const P other)
{
	
		delete[] ptr;
		deg = 1;
		ptr = new P[deg];
		ptr[0] = other;
	
	return *this;
}
template<typename P>polynomial<P>& polynomial<P>::operator=(const polynomial<P>& other)
{
	if (this != &other)
	{
		delete[] ptr;
		deg = other.deg;
		ptr = new P[deg];

		for (uint64_t i = 0; i < deg; ++i)
		{
			ptr[i] = other.ptr[i];
		}
	}

	return *this;
}

template<typename P>std::istream& operator>>(std::istream& in, polynomial<P>& plnm)
{
	uint64_t new_size;
	in >> new_size;

	polynomial<P> temp; temp.newsize(new_size);

	for (uint64_t i = 0; i < new_size; ++i)
	{
		in >> temp.ptr[i];
	}

	plnm = temp;

	return in;
}
//
//template<typename P>polynomial<P> polynomial<P>::output_mode_set(std::string newmode) {
//	switch (newmode)
//	{
//	case "FULL": outm_E = output_mode::FULL;			 break;
//	case "ABBREVIATED":outm_E = output_mode::ABBREVIATED;break;
//	case "SHORT": outm_E = output_mode::SHORT;			 break;
//	default:break;
//	}
//}
template<typename P>polynomial<P> polynomial<P>::output_mode_set(uint64_t newmode) {
	switch (newmode)
	{
	case 0: outm_E = output_mode::FULL;			break;
	case 1: outm_E = output_mode::ABBREVIATED;  break;
	case 2: outm_E = output_mode::SHORT;		break;
	default:break;
	}
	return *this;
}
template<typename P>polynomial<P>& polynomial<P>::output_mode_set(output_mode new_outm_E) {
	switch (new_outm_E)
	{
	case output_mode::FULL:			outm_E = output_mode::FULL;			break;
	case output_mode::ABBREVIATED:  outm_E = output_mode::ABBREVIATED;  break;
	case output_mode::SHORT:		outm_E = output_mode::SHORT;		break;
	default:break;
	}
	return *this;
}

template<typename P>std::ostream& operator<<(std::ostream& out, const polynomial<P>& plnm)
{
	if (output_mode::FULL == plnm.outm_E) {
		out << "Degree: " << plnm.deg << ", Coefficients: ";
	}
	

	for (uint64_t i = 0; i < plnm.deg; ++i)
	{
		if (plnm.ptr[i] < 0 && (output_mode::FULL == plnm.outm_E || output_mode::ABBREVIATED == plnm.outm_E)) {
			out << "(" << plnm.ptr[i] << ")";
		}
		else {
			out << plnm.ptr[i];
		}

		if (i == 1)
		{
			out << "x";
		}
		else if (i > 1) {
			out << "x^" << i;
		}
		if (1+i< plnm.deg) {

			if (output_mode::FULL == plnm.outm_E) {
				out << " + ";
			}
			else if (output_mode::ABBREVIATED == plnm.outm_E) {
				out << "+";
			}
			else if (output_mode::SHORT == plnm.outm_E)
			{
				out << " ";
			}
		}
	}
	if (plnm.deg == 0)
	{
		out << "0";
	}

	return out;
}

template<typename P>polynomial<P> polynomial<P>::operator+(const polynomial<P>& other) const
{
	uint64_t new_size = std::max(deg, other.deg);
	polynomial<P> result; result.newsize(new_size);

	for (uint64_t i = 0; i < new_size; i++)
	{
		if (i < deg && i < other.deg)
		{
			result.ptr[i] = ptr[i] + other.ptr[i];
		}
		else if (i < deg)
		{
			result.ptr[i] = ptr[i];
		}
		else if (i < other.deg)
		{
			result.ptr[i] = other.ptr[i];
		}

	}

	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator-(const polynomial<P>& other) const
{
	uint64_t new_size = std::max(deg, other.deg);
	polynomial<P> result; result.newsize(new_size);

	for (uint64_t i = 0; i < new_size; i++)
	{
		if (i < deg && i < other.deg)
		{
			result.ptr[i] = ptr[i] - other.ptr[i];
		}
		else if (i < deg)
		{
			result.ptr[i] = ptr[i];
		}
		else if (i < other.deg)
		{
			result.ptr[i] = (other.ptr[i])*(-1);
		}

	}

	return result;
}


template<typename P>polynomial<P> polynomial<P>::operator*(const polynomial<P>& other)const
{
	uint64_t new_size = deg + other.deg;
	polynomial<P> result; result.newsize(new_size);
	for (uint64_t i = 0; i < deg; i++)
	{
		for (uint64_t j = 0; j < other.deg; j++)
		{
			result.ptr[i + j] = result.ptr[i + j] + ptr[i] * other.ptr[j];
		}
	}
	return result.cutbag();
}

template<typename P>polynomial<P> polynomial<P>::operator*(const P& other)const
{

	polynomial<P> result; result.newsize(deg);

	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i] = ptr[i] * other;
	}

	return result;

}

template<typename P>polynomial<P> polynomial<P>::operator/(const polynomial<P>& other)const
{
	polynomial<P> result(1), ans;
	if (other == 0) {
		throw std::exception("parametr 2 is zero");
	}
	if (deg < other.deg)
	{
		ans = 0;
	}
	else if ((*this) == other) {
		ans = 1;

		return ans;
	}
	else
	{
		result = (*this);
		uint64_t difference = deg - other.deg;
		for (int64_t i = difference; i >= 0; i--) {
			ans = ans >> 1;
			ans.ptr[0] = (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
			result = result - (other >> i) * (ans.ptr[0]);
		}
	}

	return ans.cutbag();
}

template<typename P>polynomial<P> polynomial<P>::operator%(const polynomial<P>& other)const
{
	polynomial<P> result(1), ans;
	if (other == 0) {
		result = (*this);

	}
	else if (deg < other.deg)
	{
		result.ptr[0] = 0;
	}
	else
	{
		result = (*this);
		uint64_t difference = deg - other.deg;
		for (int64_t i = difference; i >= 0; i--) {

			ans = ans >> 1;
			ans.ptr[0] = (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
			result = result - (other >> i) * (ans.ptr[0]);
		}


	}

	return result.cutbag();
}
template<typename P>bool  polynomial<P>::operator>(const polynomial<P>& other)const
{
	return this->deg > other.deg;
}

template<typename P>bool  polynomial<P>::operator<(const polynomial<P>& other)const
{
	return this->deg <  other.deg;
}
template<typename P>P& polynomial<P>::operator[](uint64_t index) const {
	if (index >= deg) {
		throw std::out_of_range("Index out of range");
	}
	return ptr[index];
}
template<typename P>polynomial<P> polynomial<P>::cutbag() const {
	uint64_t new_deg = deg;
	while (new_deg > 0 && ptr[new_deg - 1] == 0) {
		new_deg--;
	}

	polynomial<P> new_poly;
	new_poly.newsize(new_deg);
	std::copy(ptr, ptr + new_deg, new_poly.ptr);

	return new_poly;
}

template<typename P>void polynomial<P>::the_root_constructor(const std::vector<P>& array_x) {
	polynomial<P> a(1);
	for (uint64_t i = 0; i < array_x.size(); i++)
	{
		polynomial<P> b(1);
		b = 1; b = b >> 1;
		b = b - array_x[i];
		a = a * b;
	}
	*this = a;
}

template<typename P>P polynomial<P>::maximum_abs(P a, P b) const {
	using namespace polynomialfunctions;
	P ans_max = 0;
	std::vector<std::pair<P, int>> V = polynomialfunctions::plnm_roots(this->get_first_derrivate(), DBL_EPSILON);
	//std::cout << this<<" " << *this << '\n';
	for (auto& i : V)
	{
		if (a < i.first && i.first < b) {
			//std::cout << "i.first:" << i.first << " (*this)(i.first): " << ((*this)(i.first)) << "abs:" << polynomialfunctions::abs((*this)(i.first))<< " \n";
			ans_max = polynomialfunctions::max(ans_max,abs((*this)(i.first) ));
		}
	}
	ans_max =max(ans_max, polynomialfunctions::abs((*this)(a)));
	ans_max =max(ans_max, polynomialfunctions::abs((*this)(b)));
	//std::cout <<"ans_max:" << ans_max << '\n';
	return ans_max;
}
