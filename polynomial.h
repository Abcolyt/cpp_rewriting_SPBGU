#pragma once

#include<iostream>
#include <sstream>
template<typename P> class polynomial;
template<typename P> std::ostream& operator<<  (std::ostream& out, const polynomial<P>& plnm);
template<typename P> std::istream& operator>> (std::istream& in, polynomial<P>& plnm);
template<typename P> class polynomial
{
public:

	P* ptr;
	uint64_t deg;
	polynomial() :polynomial(0) {};
	polynomial(uint64_t size);
	~polynomial();
	polynomial operator+(const polynomial<P>& other)const;//addition of polynomials
	polynomial operator-(const polynomial<P>& other)const;
	polynomial operator/(const polynomial<P>& other)const;//binary division of a polynomial by a polynomial
	//polynomial operator%(const polynomial<P>& other)const;//the remainder of the division of a polynomial by a polynomial
	polynomial operator*(const polynomial<P>& other)const;//binary polynomial multiplication 
	polynomial operator*(const P& other)const; //binary polynomial multiplication by an element from the field

	friend std::ostream& operator<<<>(std::ostream& out, const polynomial<P>& plnm); // overloading the output operator
	friend std::istream& operator>><>(std::istream& in, polynomial<P>& plnm); // overloading the input operator

	polynomial<P> operator>>(const uint64_t power)const;//coefficient shift (increase). or (plnm*x^n)
	polynomial<P> operator<<(const uint64_t power)const;//coefficient shift(decrease). or the whole part of (plnm* x^(-n))

	//P& operator[](uint64_t index);//the access operator to the polynomial coefficient

	polynomial<P>& operator=(const polynomial<P>& other);

	bool operator==(const polynomial<P>& other)const;
	bool operator!=(const polynomial<P>& other)const;
};

template<typename P>polynomial<P> polynomial<P>::operator>>(const uint64_t power)const
{
	polynomial result(deg + power);
	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i + power] = ptr[i];
	}
	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator<<(const uint64_t power)const
{
	polynomial result(deg - power);
	for (uint64_t i = 0; i <= (deg - power); i++)
	{
		result.ptr[i] = ptr[i + power];
	}
	return result;
}

template<typename P>bool polynomial<P>::operator==(const polynomial<P>& other)const //no working //C2666
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

template<typename P>bool polynomial<P>::operator!=(const polynomial<P>& other)const //no working //C2666
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

template<typename P>polynomial<P>::polynomial(uint64_t size) {
	deg = size;
	if (size)
		ptr = new P[size];
	else
		ptr = nullptr;
	for (uint64_t i = 0; i < size; i++)
		ptr[i] = 0;

}

template<typename P>polynomial<P>::~polynomial()
{
	delete[] ptr;
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

	polynomial<P> temp(new_size);

	for (uint64_t i = 0; i < new_size; ++i)
	{
		in >> temp.ptr[i];
	}

	plnm = temp;

	return in;
}

template<typename P>std::ostream& operator<<(std::ostream& out, const polynomial<P>& plnm)
{
	out << "Degree: " << plnm.deg << ", Coefficients: ";

	for (uint64_t i = 0; i < plnm.deg; ++i)
	{
		out << plnm.ptr[i];
		if (i == 1)
		{
			out << "x";
		}
		else if (i > 1) {
			out << "x^" << i;
		}
		out << " ";

	}

	out << std::endl;

	return out;
}

template<typename P>polynomial<P> polynomial<P>::operator+(const polynomial<P>& other) const
{
	uint64_t new_size = std::max(deg, other.deg);
	polynomial<P> result(new_size);

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
	polynomial<P> result(new_size);

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
			result.ptr[i] = -other.ptr[i];
		}

	}

	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator*(const polynomial<P>& other)const
{
	uint64_t new_size = deg + other.deg + (deg == 0 || other.deg == 0 ? 0 : -1);
	polynomial<P> result(new_size);
	for (uint64_t i = 0; i < deg; i++)
	{
		for (uint64_t j = 0; j < other.deg; j++)
		{
			result.ptr[i + j] = result.ptr[i + j] + ptr[i] * other.ptr[j];
		}
	}
	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator*(const P& other)const
{

	polynomial<P> result(deg);

	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i] = ptr[i] * other;
	}

	return result;

}
//no working

template<typename P>polynomial<P> polynomial<P>::operator/(const polynomial<P>& other)const
{
	polynomial<P> result(1);
	if (deg < other.deg)
	{
		result.ptr[0] = 0;
	}
	else if ((*this) == other) {
		result.ptr[0] = 1;
	}
	else
	{
		result = (*this);
		uint64_t difference = deg - other.deg;
		for (int64_t i = difference; i >= 0; i--) {
			std::cout << "i:" << i << "\n";
			std::cout << "koef:" << (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1]))) << "\n";
			std::cout << "result.ptr[other.deg + i]" << result.ptr[other.deg + i - 1] << "\n";
			std::cout << "other >> i" << (other >> i) << "\n";
			result = result - (other >> i) * (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
			std::cout << "res:" << result << "\n";
		}


	}

	return result;
}



#if 0
int main() {
	using T = int;
	polynomial<T> a;
	std::cin >> a;
	polynomial<T> b;
	std::cin >> b;

	//bool c = (a == b);
	std::cout << a << " " << b << " " << a * b;

	system("pause");
	return 0;
}
#endif