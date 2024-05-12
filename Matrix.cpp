#include "matrix.h"

template<typename T>matrix<T>::~matrix()
{
    delete[] ptr;
}

template<typename T>matrix<T> matrix<T>::sqprediag(const uint64_t S)const
{
    matrix<T> temp(S, S);
    for (int i = S - 1; i >= 0; i--)
    {
        temp[i][S - 1 - i] = 1;
    }
    return temp;
}

template<typename T>matrix<T> matrix<T>::transpose() const {

    matrix<T> temp(rowsize, colsize);

    for (size_t i = 0; i < colsize; ++i) {
        for (size_t j = 0; j < rowsize; ++j) {
            temp[j][i] = (*this)[i][j];
        }
    }
    return temp;
}

template<typename T>matrix<T> matrix<T>::inverse_M()const
{
    if (colsize != rowsize)
    {
        throw std::invalid_argument("the number of rows in the matrix is not equal to the number of columns");
    }
    uint64_t sizeM = (*this).colsize;
    matrix<T> a = (*this), reversM, E_lower_diagonal, t1;
    E_lower_diagonal = a.sqprediag(sizeM);
    reversM = E_lower_diagonal * E_lower_diagonal;

    t1 = a.to_uptrng(reversM);

    t1 = ((t1 * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();
    reversM = ((reversM * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();

    t1 = t1.to_uptrng(reversM);

    t1 = ((t1 * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();
    reversM = ((reversM * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();

    for (uint64_t i = 0; i < sizeM; i++)
    {
        T koef = t1[i][i];
        if (koef == 0)
        {
            throw std::invalid_argument("the matrix is not reversible");
        }
        for (uint64_t k = 0; k < sizeM; k++)
        {
            t1[i][k] = t1[i][k] / koef;
            reversM[i][k] = reversM[i][k] / koef;

        }

    }
    return reversM;
}

template<typename T>matrix<T> matrix<T>::to_uptrng(matrix<T>& other)const
{
    if ((this->colsize == this->rowsize) && (this->colsize == (other.colsize)) && ((other.colsize) == other.rowsize))
    {



        matrix<T> temp((*this));
        for (uint64_t i = 0; i < rowsize - 1; i++) {
            for (uint64_t j = i + 1; j < rowsize; j++) {

                while (temp[i][i] == 0 && temp[j][i] == 0) {
                    j++;
                    if (j >= rowsize) {
                        throw "писец";

                        //return хз что ;
                    }
                }

                if (temp[i][i] == 0) {
                    for (uint64_t k = 0; k < rowsize; k++) {
                        std::swap(temp[i][k], temp[j][k]);
                        std::swap(other[i][k], other[j][k]);

                    }

                }

                T koef1 = temp[j][i] / temp[i][i];
                for (uint64_t k = 0/*=i*/; k < rowsize; k++) {
                    temp[j][k] = temp[j][k] - temp[i][k] * koef1;
                    other[j][k] = other[j][k] - other[i][k] * koef1;
                }

            }

        }
        return temp;
    }
    throw std::invalid_argument("[  ((this->colsize == this->rowsize) && (this->colsize == (other.colsize)) && ((other.colsize) == other.rowsize)) ]==0");
    abort();
}

template<typename T>matrix<T> matrix<T>::to_uptrng()const
{
    matrix<T> temp((*this));

    T det = 1;
    for (uint64_t i = 0; i < rowsize - 1; i++) {
        for (uint64_t j = i + 1; j < rowsize; j++) {

            while (temp[i][i] == 0 && temp[j][i] == 0) {
                j++;
                if (j >= rowsize) {
                    throw std::invalid_argument("the matrix is irreducible to the triangular form");


                }
            }

            if (temp[i][i] == 0) {
                for (uint64_t k = 0; k < rowsize; k++) {
                    std::swap(temp[i][k], temp[j][k]);
                    det = det * (-1);
                }

            }

            T koef = temp[j][i] / temp[i][i];

            for (uint64_t k = i; k < rowsize; k++) {
                temp[j][k] = temp[j][k] - temp[i][k] * koef;

            }

        }

    }
    for (uint64_t i = 0; i < rowsize - 1; i++) {
        temp[0][i] = (temp[0][i]) * det;
    }
    return temp;
}

template <typename T>T matrix<T>::determinant() const {
    if (rowsize != colsize) {
        //throw()
        return 0;
    }

    matrix<T> temp((*this).to_uptrng());

    T det = 1;
    for (uint64_t i = 0; i < rowsize; i++) {
        det *= temp[i][i];
    }

    return det;
}

template<typename T> std::ostream& operator<<(std::ostream& out, const matrix<T>& mtrx)
{
    out << "sizex:" << mtrx.getcol() << "sizey:" << mtrx.getrow() << "\n";
    for (uint64_t i = 0; i < mtrx.getcol(); i++)
    {
        for (uint64_t j = 0; j < mtrx.getrow(); j++)
        {
            out << "[" << i << "][" << j << "] = \t" << mtrx[i][j] << "\t | ";
        }
        out << "\n";
    }
    return out;
}

template<typename T>std::istream& operator>>(std::istream& in, matrix<T>& p)
{
    for (uint64_t i = 0; i < p.getcol(); i++) {
        for (uint64_t j = 0; j < p.getrow(); j++) {
            in >> p[i][j];
        }
    }
    return in;
}

template<typename T>matrix<T>& matrix<T>::operator=(const matrix<T>& other) // assignment operator overload
{
    if (this != &other)
    {

        delete[] ptr;
        ptr = new T[other.getcol() * other.getrow()];
        this->colsize = other.getcol();
        this->rowsize = other.getrow();

        for (uint64_t i = 0; i < other.getrow(); i++)
        {
            for (uint64_t j = 0; j < other.getcol(); j++)
            {
                (*this)[i][j] = other[i][j];
            }
        }
    }

    return *this;
}

template<typename T>matrix<T> matrix<T>::operator+(const matrix<T>& other) const
{
    matrix<T> result(other.getcol(), other.getrow());
    if ((other.getcol()) == ((*this).getcol()) && (other.getrow() == (*this).getrow())) {

        for (uint64_t i = 0; i < other.getrow(); i++)
        {
            for (uint64_t j = 0; j < other.getcol(); j++)
            {
                result[i][j] = (other[i][j] + (*this)[i][j]);
            }

        }

    }
    else {
        throw std::invalid_argument("((other.getcol()) == ((*this).getcol()) && (other.getrow() == (*this).getrow())) == false");
    }
    return result;

}

template<typename T>matrix<T> matrix<T>::operator-(const matrix<T>& other) const
{
    matrix<T> result(other.getcol(), other.getrow());
    if ((other.getcol()) == ((*this).getcol()) && (other.getrow() == (*this).getrow())) {

        for (uint64_t i = 0; i < other.getrow(); i++)
        {
            for (uint64_t j = 0; j < other.getcol(); j++)
            {
                result[i][j] = (other[i][j] - (*this)[i][j]);
            }

        }

    }
    return result;

}

template<typename T>matrix<T> matrix<T>::operator*(const T& other) const
{
    matrix<T> result(this->getcol(), this->getrow());


    for (uint64_t i = 0; i < this->getrow(); i++)
    {
        for (uint64_t j = 0; j < this->getcol(); j++)
        {
            result[i][j] = other * ((*this)[i][j]);
        }

    }


    return result;
}

template<typename T>matrix<T> matrix<T>::operator*(const matrix<T>& other) const
{
    if (rowsize != other.colsize)//problematic condition???
    {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");

    }

    matrix<T> result(colsize, other.rowsize);

    for (uint64_t i = 0; i < colsize; ++i) {
        for (uint64_t j = 0; j < other.rowsize; ++j) {
            for (uint64_t k = 0; k < rowsize; ++k) {
                result[i][j] = result[i][j] + ((*this)[i][k] * other[k][j]);
            }
        }
    }

    return result;
}

template<typename T>matrix<T> matrix<T>::operator/(const matrix<T>& other) const
{
    if (other.determinant() != 0) {
        return (*this) * (other.inverse_M());
    }
    throw std::invalid_argument("the matrix on the right is not reversible");
}


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