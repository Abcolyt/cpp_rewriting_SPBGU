// cpp_rewriting_SPBGU.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//





#define _CRT_SECURE_NO_WARNINGS 
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<locale.h>
#include<stdbool.h>
#include"complex_number.cpp"

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include<math.h>
#include<cassert>
class Complex

{
private:
	double R;
	double I;
public:
	Complex(double r, double i) : R(r), I(i) {}
	Complex() : Complex(0, 0) {}
	Complex(const Complex& other) { this->R = other.R; this->I = other.I; }

	double real_Get();
	double imag_Get();
	void real_Set(double r);
	void imag_Set(double i);
	std::pair<double, double> get_RI();
	void set_RI(std::pair<double, double> a);
	void set_RI(double a, double b);


	Complex& conjugate() { Complex out((const_cast<Complex*>(this))->R, (const_cast<Complex*>(this))->I *= -1); return out; }

	Complex& operator+(const Complex& other)
	{
		Complex output;
		output.R = this->R+ other.R;
		output.I = this->I+ other.I;
		return output;
	}
	Complex& operator-(const Complex& other)
	{
		Complex output;
		output.R = this->R  - other.R;
		output.I = this->R  - other.I;
		return output;
	}
	Complex& operator=(const Complex& other)
	{
		this->R = other.R;
		this->I = other.I;
		return *this;
	}
	bool operator==(const Complex& other) { 
		return this->R == other.R && this->I == other.I ? 1 : 0; }

	Complex& operator*(const Complex& other) {
		Complex out;
		out.R = this->R * other.R - this->I * other.I;
		out.I = this->I * other.R + other.I * this->R;
		return out;
	}

	
	Complex& operator/(const Complex& other)  {
		//(a+bi)/(c+di)= (a+bi)*(c-di) / ((c*c+d*d))
		double divisor = (other.R) * (other.R) + (other.I) * (other.I);
		
		Complex other_(other.R, other.I * (-1));
		Complex out(((*this * other_).R) / divisor, ((*this * other_).I) / divisor);
		
		return out;

	}
};

double Complex::real_Get()
{
	return this->R;
}
double Complex::imag_Get()
{
	return this->I;
}
void Complex::real_Set(double r)
{
	this->R = r;
}
void Complex::imag_Set(double i)
{
	this->I = i;
}
std::pair<double, double> Complex::get_RI()
{
	return std::make_pair(this->R, this->I);
}
void Complex::set_RI(double a, double b)
{
	this->R = a;
	this->I = b;
}
void Complex::set_RI(std::pair<double, double> a)
{
	this->R = a.first;
	this->I = a.second;
}

void TESTs() {
	Complex a, b, c;

	a = Complex(1, 1);
	b = Complex(-1, -1);
	c = Complex(0, 0);
	assert((a + b) == (c));



	a = Complex(1, 1);
	b = Complex(-1, -1);
	c = Complex(2, 2);
	assert((a - b) == (c));


	a = Complex(1, 1);
	b = Complex(1, 1);
	c = Complex(0, 2);

	
	assert((a * b) == (c));

	a = Complex(1, 1);
	b = Complex(1, 1);
	c = Complex(1, 0);
	Complex d((a / b));
	
	assert( d == c);
	assert( (a / b) == c);
	
}



#if 0
//+
double** _M_malloc(int str, int column) {
	double** M = (double**)malloc(column * sizeof(double*) + column * str * sizeof(double));
	double* start = (double*)((char*)M + column * sizeof(double*));
	for (int i = 0; i < column; i++)
	{
		M[i] = ((start)+i * str);
	}
	return M;
}
class Matrix
{

public:
	Matrix() :Matrix(0, 0) {  };
	Matrix(int str, int column);
	void print_M();
	void scanf_M();
	Matrix(const Matrix& other);
	Matrix& operator=(const Matrix& other);
	Matrix& operator+(const Matrix& other);
	Matrix* get_cut_M(int nocolumn, int nostring);
	//double det_Mn_n_recursion(const Matrix *A) {
	//	double sum = 0;
	//	if ((this->number_of_str == 1) && (this->number_of_columns == 1)) {
	//		return (this->Mptr[0][0]);
	//	}
	//	if (this->number_of_str == 2) {
	//		return (this->Mptr[0][0] * this->Mptr[1][1] - this->Mptr[0][1] * this->Mptr[1][0]);
	//	}
	//	for (int i = 0; i < this->number_of_columns; i++) {
	//		//decomposition by 1 =[0] string
	//
	//		sum += ((i % 2 == 0) ? 1 : -1) * (this->Mptr[i][0]) * det_Mn_n_recursion(get_cut_M(i, 0));
	//		free(M_in.Mptr);
	//	}
	//	return sum;
	//}


	~Matrix();

private:
	int number_of_str, number_of_columns;//str,colm
	double** Mptr;
};
Matrix::Matrix(int str, int column)
{
	std::cout << "constructor"<<this<<" " << str << "   " << column << "\n";
	if (str == 0 || column == 0) {
		this->number_of_str = 0; this->number_of_columns = 0; this->Mptr = NULL;
	}
	else {


		this->number_of_str = str;
		this->number_of_columns = column;
		if (NULL != (this->Mptr)) { free((this->Mptr)); }
		this->Mptr = _M_malloc(str, column);
		for (int i = 0; i < (this->number_of_str); i++)
		{
			for (int j = 0; j < (this->number_of_columns); j++)
			{
				this->Mptr[j][i] = 0;
			}
		}
	}
}
void Matrix::scanf_M() {
if (NULL != this->Mptr) {


	for (int i = 0; i < (this->number_of_str); i++)
	{
		for (int j = 0; j < (this->number_of_columns); j++)
		{
			scanf("%lf", &(this->Mptr[j][i]));
		}

	}

}
}
inline Matrix::Matrix(const Matrix& other) {
	this->number_of_columns = other.number_of_columns;
	this->number_of_str = other.number_of_str;
	this->Mptr = _M_malloc(other.number_of_str, other.number_of_columns);
	for (int i = 0; i < this->number_of_columns; i++) {
		for (int j = 0; j < this->number_of_str; j++) {
			this->Mptr[i][j] = other.Mptr[i][j];
		}
	}
}
inline Matrix& Matrix::operator=(const Matrix& other) {
	std::cout << "operator=" << this << " " <<&other <<"\n";
	if (this != &other) {
		if (NULL != this->Mptr)free(this->Mptr);
	}
	this->number_of_str=other.number_of_str;
	this->number_of_columns=other.number_of_columns;
	this->Mptr=_M_malloc(other.number_of_str,other.number_of_columns);
	for (int i = 0; i < other.number_of_str; i++)
	{
		for (int j = 0; j < other.number_of_columns; j++)
		{	
			this->Mptr[i][j] = other.Mptr[i][j];
		}

	}
	return *this;
}
inline Matrix& Matrix::operator+(const Matrix& other) {
	if (NULL != other.Mptr &&this->number_of_str==other.number_of_str && this->number_of_columns==other.number_of_columns) {
		Matrix sum(this->number_of_str,this->number_of_columns);
		for (int i = 0; i < other.number_of_str; i++)
		{
			for (int j = 0; j < other.number_of_columns; j++)
			{
			sum.Mptr[i][j] = this->Mptr[i][j] + other.Mptr[i][j];
			}

		}
		return sum;
	}
	
}
inline Matrix* Matrix::get_cut_M(int nocolumn, int nostring)
{
	Matrix* M_in = new Matrix(this->number_of_str - 1, this->number_of_columns - 1);
	//initM0(M_in, M_out->number_of_str - 1, M_out->number_of_columns - 1);
	for (int i, i_ = i = 0; i < (M_in->number_of_str); i_++, i++)
	{
		if (i == nocolumn) i_++;
		for (int j, j_ = j = 0; j < (M_in->number_of_columns); j_++, j++)
		{
			if (j == nostring)  j_++;
			(M_in->Mptr)[i][j] = this->Mptr[i_][j_];
		}
	}
	return M_in;
}

void Matrix::print_M()
{

	if (this->Mptr != NULL) {
		printf("![col][str]\n");
		for (int i = 0; i < (this->number_of_str); i++)
		{
			for (int j = 0; j < (this->number_of_columns); j++)
			{
				printf("M[%d][%d]=%lf\t", j, i, this->Mptr[j][i]);
			}printf("\n");
		}
		printf("ptr=%p number_of_str=%d number_of_columns= %d\n", this->Mptr, this->number_of_str, this->number_of_columns);
	}
}

Matrix::~Matrix()
{
	std::cout << "destructor" << this << "\n";
	if (NULL != (this->Mptr))free((this->Mptr));
}






//матрица через классы





/*

typedef struct matrix
{
	int number_of_str, number_of_columns;//str,colm
	double** Mptr;
}matrix;
#define matrix(out) 
#define Dmatrix(out) matrix* out = (matrix*)malloc(sizeof(matrix)); out->number_of_str = 0; out->number_of_columns = 0; out->Mptr = NULL;
matrix* get_dM() {
	matrix* out = (matrix*)malloc(sizeof(matrix)); out->number_of_str = 0; out->number_of_columns = 0; out->Mptr = NULL;
}

//+
void print_M(matrix* M) {

	if (M->Mptr != NULL) {
		printf("![col][str]\n");
		for (int i = 0; i < (M->number_of_str); i++)
		{
			for (int j = 0; j < (M->number_of_columns); j++)
			{
				printf("M[%d][%d]=%lf\t", j, i, M->Mptr[j][i]);
			}printf("\n");
		}
		printf("ptr=%p number_of_str=%d number_of_columns= %d\n", M->Mptr, M->number_of_str, M->number_of_columns);
	}
}
//+
void scanf_M(matrix* M) {
	if (NULL != M->Mptr) {


		for (int i = 0; i < (M->number_of_str); i++)
		{
			for (int j = 0; j < (M->number_of_columns); j++)
			{
				scanf("%lf", &(M->Mptr[j][i]));
			}

		}

	}
}


extern void initM0(matrix* M, int str, int column);







//+
void get_copy_M(const matrix* M_out, matrix* M_in)//âîçâðàùàåò êîïèþ íà òîò 2 ìåðíûé ìàññèâ êîòîðûé áûë ïåðåäàí 1 ïàðàìåòðîì
{
	initM0(M_in, M_out->number_of_str, M_out->number_of_str);


	for (int i = 0; i < (M_in->number_of_str); i++)
	{
		for (int j = 0; j < (M_in->number_of_columns); j++)
		{
			M_in->Mptr[i][j] = M_out->Mptr[i][j];
		}
	}


}


//+
matrix* get_cut_M(const matrix* M_out, matrix* M_in, int nocolumn, int nostring)
{

	initM0(M_in, M_out->number_of_str - 1, M_out->number_of_columns - 1);
	for (int i, i_ = i = 0; i < (M_in->number_of_str); i_++, i++)
	{
		if (i == nocolumn) i_++;
		for (int j, j_ = j = 0; j < (M_in->number_of_columns); j_++, j++)
		{
			if (j == nostring)  j_++;
			(M_in->Mptr)[i][j] = M_out->Mptr[i_][j_];
		}
	}
	return M_in;
}

//+
void initM0(matrix* M, int str, int column)
{
	M->number_of_str = str;
	M->number_of_columns = column;
	if (NULL != (M->Mptr))free((M->Mptr));
	M->Mptr = _M_malloc(str, column);
	for (int i = 0; i < (M->number_of_str); i++)
	{
		for (int j = 0; j < (M->number_of_columns); j++)
		{
			M->Mptr[j][i] = 0;
		}
	}
}
//+

double det_Mn_n_1_iter(const double* M, const int n) //((double*)M, not working with struct matrix object
{
	double plus_ans = 0, minus_ans = 0;
	for (int i = 0; i < n; i++)
	{
		double pre_plus_ans = 1;
		for (int j = 0; j < n; j++) {
			pre_plus_ans *= M[((j + i) % n) * n + j];
		}
		plus_ans += pre_plus_ans;
	}

	for (int i = 0; i < n; i++)
	{
		double pre_minus_ans = 1;
		for (int j = 0; j < n; j++) {
			pre_minus_ans *= M[((j + i) % n) * n + (n - j - 1)];
		}
		minus_ans += pre_minus_ans;
	}

	return plus_ans - minus_ans;
}
//+
double det_Mn_n2_iter(const matrix* const M) {

	double plus_ans = 0, minus_ans = 0;
	for (int i = 0; i < M->number_of_columns; i++)
	{
		double pre_plus_ans = 1;
		for (int j = 0; j < M->number_of_str; j++)
		{
			pre_plus_ans *= M->Mptr[j][((j + i) % M->number_of_columns)];
		}
		plus_ans += pre_plus_ans;
	}

	for (int i = 0; i < M->number_of_columns; i++)
	{
		double pre_minus_ans = 1;
		for (int j = 0; j < M->number_of_str; j++)
		{
			pre_minus_ans *= (M->Mptr[(M->number_of_columns) - j - 1][(j + i) % M->number_of_columns]);
		}
		minus_ans += pre_minus_ans;
	}

	return plus_ans - minus_ans;
}


double det_Mn_n_recursion(const matrix* M) {
	double sum = 0;
	if ((M->number_of_str == 1) & (M->number_of_columns == 1)) {
		return (M->Mptr[0][0]);
	}
	if (M->number_of_str == 2) {
		return (M->Mptr[0][0] * M->Mptr[1][1] - M->Mptr[0][1] * M->Mptr[1][0]);
	}
	for (int i = 0; i < M->number_of_columns; i++) {
		//decomposition by 1 =[0] string
		matrix(M_in);
		sum += ((i % 2 == 0) ? 1 : -1) * (M->Mptr[i][0]) * det_Mn_n_recursion(get_cut_M(M, &M_in, i, 0));
		free(M_in.Mptr);
	}
	return sum;
}
void swap(double* n, double* m) {
	double t = *n; *n = *m; *m = t;
}
void transp_M(matrix* M) {

	for (int i = 0; i < M->number_of_columns; i++)
	{
		for (int j = i; j < M->number_of_str; j++)
		{
			swap(&(M->Mptr[i][j]), &(M->Mptr[j][i]));
		}
	}

}
void transp_M_all(matrix* M) {
	matrix(M_new)

		for (int i = 0; i < M->number_of_columns; i++)
		{
			for (int j = i; j < M->number_of_str; j++)
			{
				swap(&(M->Mptr[i][j]), &(M->Mptr[j][i]));
			}
		}

}
void get_UnM(const matrix* M_out, matrix* M_in) {

	initM0(M_in, M_out->number_of_str, M_out->number_of_columns);
	for (int i = 0; i < M_in->number_of_str; i++)
	{
		for (int j = 0; j < M_in->number_of_columns; j++)
		{
			matrix(M_preans)
				get_cut_M(M_out, &M_preans, j, i);
			//printf("str=%j\tcol=%i\n",j,i);
			//print_M(&M_preans);
			//printf("det[%j][%i]=%lf",j,i, det_Mn_n_recursion(&M_preans));
			(M_in->Mptr[j][i]) = (((i + j) % 2 == 1) ? -1 : 1) * det_Mn_n_recursion(&M_preans);
			free(M_preans.Mptr);
		}



	}
	transp_M(M_in);

}


//+
const matrix* scalar_M_division(const matrix* M, double divisor) {
	for (int i = 0; i < M->number_of_columns; i++)
	{
		for (int j = 0; j < M->number_of_str; j++)
		{
			(M->Mptr[i][j]) /= divisor;
		}
	}
	return M;
}
//-
matrix* get_inverse_M(const matrix* M, matrix* b) {

	initM0(b, M->number_of_columns, M->number_of_str);
	scalar_M_division(b, det_Mn_n_recursion(M));



	return b;
}

//+
matrix* multiply_M(const matrix* const M, const matrix* const M2)//the first parameter has multiplication by columns, the second by rows
{
	if (M2->number_of_columns != M->number_of_str)
		return NULL;
	matrix* out = get_dM();
	initM0(out, M->number_of_str, M2->number_of_columns);
	for (int i = 0; i < M->number_of_str; i++)
	{
		for (int j = 0; j < M2->number_of_columns; j++)
		{
			for (int k = 0; k < M2->number_of_str; k++)
			{
				out->Mptr[i][j] += M2->Mptr[i][k] * M->Mptr[k][j];
			}
		}
	}
	return out;
}
enum TypeLogik
{
	CramerMatrixSolver = 0,
	reversefinder = 1,
};

typedef struct TpLg
{
	enum TypeLogik log;
} TpLg;
*/
int main() {
	Matrix a, b, c(2,2);
	c.scanf_M();
	a = b = c;
	a.print_M();
	a = (a + a);
	a.print_M();
	return 0;
}
/*
int main()
{

	TpLg l1;
	printf("CramerMatrixSolver(0) or reversefinder(1)?\n");
	scanf("%d", &l1.log);
	printf("%d", l1.log);



	switch (l1.log)
	{
	case reversefinder: {
		int size = 2;
		printf("enter Marix size(int type ,the square matrix): ");
		scanf("%d", &size);
		matrix(a)
			matrix(b)
			initM0(&a, size, size);
		initM0(&b, size, size);
		printf("enter matrix elements");
		scanf_M(&a);
		print_M(&a);
		get_UnM(&a, &b);
		if (det_Mn_n_recursion(&a) != 0)
		{
			scalar_M_division(&b, det_Mn_n_recursion(&a));
			printf("the inverse matrix:\n");
			print_M(&b);
			printf("the original matrix multiplied by the inverse:\n");
			print_M(multiply_M(&a, &b));
		}
		else printf("error(det (M) =0)\n The matrix does not exist!!\n");
		break;
	}
	case CramerMatrixSolver:
	{
		int size = 2;

		printf("enter Marix size(int type,the square matrix): ");
		scanf("%d", &size);
		matrix(a)
			matrix(b)
			initM0(&a, size, size);
		initM0(&b, size, size);
		printf("enter matrix elements\n");
		scanf_M(&a);
		printf("+\n");
		if (det_Mn_n_recursion(&a) != 0) {
			matrix(free_members);
			initM0(&free_members, 1, size);
			get_copy_M(&a, &b);
			scanf_M(&free_members);
			for (int j = 0; j < size; j++)
			{
				get_copy_M(&a, &b);
				for (int i = 0; i < size; i++) {

					b.Mptr[j][i] = free_members.Mptr[i][0];//вставляем в массив а столбец свободных коэффициентов

				}
				printf("%lf\n", (det_Mn_n_recursion(&b) / det_Mn_n_recursion(&a)));
			}

			break;
		}
		else printf("error(det (M) =0)\n");

	}
	default:
		printf("%s", "unknown command");
		return 666666;
	}








	system("pause");
	return 0;
}
*/

#endif


int main() {
	TESTs();
	return 0;
}
