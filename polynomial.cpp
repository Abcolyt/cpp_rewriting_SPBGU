#pragma once
#include"polynomial.h"
template<typename P>polynomial<P> polynomial<P>::operator>>(const uint64_t power)const
{
	if (power == 0)return *this;
	polynomial result(deg + power);
	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i + power] = ptr[i];
	}
	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator<<(const uint64_t power)const
{
	if (power == 0)return *this;
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
		ptr[i] = (P)0;

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
	//del
	//std::cout << "\nother in operator=" << other<<"\n";
	/*if (other == 1) {
		delete[] ptr;
		deg = other;
		ptr = new P[deg];
		ptr[0] = other;
	}	*/
	//if (other == 0) {
	//	delete[] ptr;
	//	deg = other;
	//	ptr = nullptr;
	//	
	//}

	return *this;
}
template<typename P>polynomial<P>& polynomial<P>::operator=(const polynomial<P>& other)
{
	if (this != &other)
	{
		std::cout << "\nother:" << other << "\n";
		std::cout << "\nptr:" << ptr << "\n";
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
	//it is not necessary if there will be inside another class
	// 
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
	if (plnm.deg == 0)
	{
		out << "0";
	}
	//out << std::endl;

	return out;
}

template<typename P>polynomial<P>& polynomial<P>::operator+(const polynomial<P>& other) const
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

template<typename P>polynomial<P>& polynomial<P>::operator-(const polynomial<P>& other) const
{
	uint64_t new_size = std::max(deg, other.deg);
	polynomial<P> result(new_size),zero(1);

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
			result.ptr[i] = zero.ptr[i] - other.ptr[i];
		}

	}

	return result;
}

template<typename P>polynomial<P>& polynomial<P>::operator*(const polynomial<P>& other)const
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

template<typename P>polynomial<P>& polynomial<P>::operator*(const P& other)const
{

	polynomial<P> result(deg);

	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i] = ptr[i] * other;
	}

	return result;

}

template<typename P>polynomial<P>& polynomial<P>::operator/(const polynomial<P>& other)const
{
	polynomial<P> result(1), ans;
	if (deg < other.deg)
	{
		result.ptr[0] = (P)0;
	}
	else if ((*this) == other) {
		result.ptr[0] = (P)1;
	}
	else
	{
		result = (*this);
		uint64_t difference = deg - other.deg;
		for (int64_t i = difference; i >= 0; i--) {

			ans = ans >> 1;
			ans.ptr[0] = (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
			result = result - (other >> i) * (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
		}


	}

	return ans;
}

template<typename P>polynomial<P>& polynomial<P>::operator%(const polynomial<P>& other)const
{
	polynomial<P> result(1), ans;
	if (deg < other.deg)
	{
		result.ptr[0] = (P)0;
	}
	else if ((*this) == other) {
		result.ptr[0] = (P)1;
	}
	else
	{
		result = (*this);
		uint64_t difference = deg - other.deg;
		for (int64_t i = difference; i >= 0; i--) {

			ans = ans >> 1;
			ans.ptr[0] = (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
			result = result - (other >> i) * (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
		}


	}

	return result;
}

