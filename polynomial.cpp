#pragma once
#include"polynomial.h"
template<typename P>polynomial<P> polynomial<P>::operator>>(const uint64_t power)const
{
	polynomial result;result.newsize(deg + power);
	for (uint64_t i = 0; i < deg; i++)
	{
		result.ptr[i + power] = ptr[i];
	}
	return result;
}

template<typename P>polynomial<P> polynomial<P>::operator<<(const uint64_t power)const
{
	polynomial result; result.newsiz(deg - power);
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

template<typename P>polynomial<P>::polynomial(P number) {
	deg = 1;
	ptr = new P[deg];
	ptr[0] = number;
	/*
	* 	deg = size;
	if (size)
		ptr = new P[size];
	else
		ptr = nullptr;
	for (uint64_t i = 0; i < size; i++)
		ptr[i] = 0;
	*/


}

template<typename P>void polynomial<P>::newsize(uint64_t size)
{
	deg = size;
	if (size)
		ptr = new P[size];
	else
		ptr = nullptr;
	for (uint64_t i = 0; i < size; i++)
		ptr[i] = 0;
}
template<typename P>polynomial<P>::polynomial(const polynomial<P>& other)
{
	deg = other.deg; // Копируем степень полинома
	ptr = new P[deg + 1]; // Выделяем память под массив коэффициентов

	for (uint64_t i = 0; i <= deg; ++i)
	{
		ptr[i] = other.ptr[i]; // Копируем коэффициенты из другого полинома
	}
}
template<typename P>polynomial<P>::~polynomial()
{
	delete[] ptr;
}
template<typename P>polynomial<P>& polynomial<P>::operator=(const P other)
{
	//del
	//std::cout << "\nother in operator=" << other<<"\n";
	
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
			result.ptr[i] = -other.ptr[i];
		}

	}

	return result;
}


template<typename P>polynomial<P> polynomial<P>::operator*(const polynomial<P>& other)const
{
	uint64_t new_size = deg + other.deg + (deg == 0 || other.deg == 0 ? 0 : -1);
	polynomial<P> result; result.newsize(new_size);
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
		//std::cout << "\nans in op/" << ans << "\n";
		return ans;
	}
	else
	{
		result = (*this);
		uint64_t difference = deg - other.deg;
		for (int64_t i = difference; i >= 0; i--) {
			std::cout << "\nans1" << ans << "\n";
			ans = ans >> 1;
			std::cout << "\nans2" << ans << "\n";
			ans.ptr[0] = (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
			std::cout << "\nans3" << ans << "\n";
			result = result - (other >> i) * (ans.ptr[0]);
			std::cout << "\nres:" << result << "\n";
		}


	}
	std::cout << "ans:" << ans;
	return ans;
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
			result = result - (other >> i) * (((result.ptr[other.deg + i - 1] / other.ptr[other.deg - 1])));
		}


	}

	return result;
}

template<typename P>bool  polynomial<P>::operator>(const polynomial<P>& other)const
{
	return this->deg > other.deg;
}

template<typename P>bool  polynomial<P>::operator<(const polynomial<P>& other)const
{
	return this->deg <  other.deg;
}