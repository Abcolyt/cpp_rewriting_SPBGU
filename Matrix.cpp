#include "matrix.h"

/*
double** _M_malloc(int str, int column)
{
	double** M = (double**)malloc(column * sizeof(double*) + column * str * sizeof(double));
	double* start = (double*)((char*)M + column * sizeof(double*));
	for (int i = 0; i < column; i++)
	{
		M[i] = ((start)+i * str);
	}
	return M;
}
*/


template <typename T> T** _M_malloc(int rows, int columns)//memory allocation function in one piece
{
	T** M = new T * [columns + rows * columns];
	T* start = (T*)((char*)M + columns * sizeof(T*));
	for (int i = 0; i < columns; i++) {
		M[i] = start + i * rows;
	}
	return M;
}

#if 0
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
#endif