#include "stdafx.h"
#include <iostream>
#include <cassert>
#include <sstream>
#pragma warning(disable : 4996)


class Matrix
{
public:
	Matrix(uint64_t h, uint64_t w):sizex(w), sizey(h){};
	~Matrix();

	Matrix operator-()const;//unary matrix inverse finding(if exists)

	Matrix operator+(const Matrix& other) const;//binary matrix addition
	Matrix operator-(const Matrix& other) const;//binary matrix subtraction
	Matrix operator*(const Matrix& other) const;//binary matrix multiplication
	Matrix operator/(const Matrix& other) const;//binary matrix division(if not a singular matrix on the left)

	Matrix operator*(const uint64_t other)const;////binary matrix multiplication by scalar 

	
	friend std::ostream& operator<<(std::ostream& out, const Matrix& p); //overloading the output operator
	friend std::istream& operator>>(std::istream& in, Matrix& p);//overloading the input operator

	Matrix& operator = (const Matrix& other)//
	//assignment operator overload
	{
		return;
	}

	Matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string acrion, Matrix& other);//overload(), extended union operator
	Matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string acrion);//overload(),advanced editing operator
private:
	uint64_t determinant() const;//finding the determinant if there is one
	uint64_t *ptr;
	uint64_t sizex;
	uint64_t sizey;

};

Matrix::Matrix(uint64_t h, uint64_t w)
{
}

Matrix::~Matrix()
{
}

Matrix Matrix::operator+(const Matrix& other)const
{
	return Matrix(1, 2);
}
Matrix Matrix::operator-(const Matrix& other)const
{
	return Matrix(2,1)
}
Matrix Matrix::operator*(const Matrix& other)const
{
	return Matrix(1, 2);
}
Matrix Matrix::operator/(const Matrix& other)const
{
	return Matrix(1, 2);
}

Matrix Matrix::operator*(const uint64_t other)const
{
	return Matrix(2, 1)
};

Matrix Matrix::operator-()const
{
	return Matrix(2, 1)
}

Matrix Matrix::operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, std::string acrion,Matrix& other)

{
	return Matrix(2, 1)
}
Matrix Matrix::operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, std::string acrion)
{
	return Matrix(2, 1);
};

uint64_t Matrix::determinant() const

{
	return 0;
};


double** _M_malloc(int str, int column)//memory allocation function in one piece
{
	double** M = (double**)malloc(column * sizeof(double*) + column * str * sizeof(double));
	double* start = (double*)((char*)M + column * sizeof(double*));
	for (int i = 0; i < column; i++)
	{
		M[i] = ((start)+i * str);
	}
	return M;
}


template <class T1> void complex_calc(T1 o) {

	static char L = 'E';
	while (1) {
		double R, I;

		T1 a, b, c;
		std::cout << "Enter first Number \n";
		std::cin >> a;
		std::cout << "Enter second Number \n";
		std::cin >> b;

		std::cout << "Enter Action('-','+','/','*','E','T') :";

		std::cin >> L;
		//getchar();
		if (L == 'E' | L == 'e')
		{
			abort();
		}
		std::cout << std::endl;
		switch (L)
		{
		case '+':c = a + b; break;//для одинаковой размерности
		case '-':c = a - b; break;
		case '*':c = a * b; break;// Для матриц при размере одной из матриц 1 на 1 воспринимать как умножение на скаляр
		case '/':c = a / b; break;//Для матриц умножить на обратное число
		case 'T':TESTs(); break;
		case 't':TESTs(); break;
		default:abort();
		}
		std::cout << c << "\n";
	}


}



int main()
{

}
