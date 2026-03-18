#include <../file_h/matrix.h>

template <typename T>
class derived_matrix : public matrix<T> {
private:
	enum class StorageMethod {
		None = 666,									//не хранит

		One_dimensional_array_of_strings = 0,	//по строкам 
		one_dimensional_array_of_columns = 1,	//по столбцам

		row_major_order = 0,					//по строкам 
		column_major_oreder = 1		    		//по столбцам
	};
	static StorageMethod BaseOrderStorage = StorageMethod::row_major_orde;
	StorageMethod OrderStorage = BaseOrderStorage;

public:

	StorageMethod GetOrderStorage()const {
		return OrderStorage;
	}
	void SetOrderStorage(StorageMethod OrderStorage)const {
		this->OrderStorage = OrderStorage;
	}
	//CONSTRUCTORS\DESTRUCTORS
	//the copying constructor
	derived_matrix(const matrix<T>& mtrx) :
		matrix<T>()	{
		//  опируем размеры и режим вывода из исходной матрицы
		this->colsize = mtrx.getcol();
		this->rowsize = mtrx.getrow();
		this->out_mode = mtrx.get_output_mode(); // предполагаем наличие геттера
												 // ¬ыдел€ем пам€ть под данные
		this->ptr = new T[this->colsize * this->rowsize];

		const matrix<T>* p = dynamic_cast<const matrix<T>*>(&mtrx);
		if (p) {
			this->OrderStorage = BaseOrderStorage;
		}
		else {
			this->OrderStorage = mtrx.OrderStorage;
		}

		if (this->OrderStorage == StorageMethod::row_major_order) {
			for (uint64_t i = 0; i < this->colsize * this->rowsize; ++i)
				this->ptr[i] = mtrx.ptr[i];
		}
		else { // column-major
			for (uint64_t i = 0; i < this->colsize; ++i) {
				for (uint64_t j = 0; j < this->rowsize; ++j) {
					(*this)[j][i] = mtrx[i][j]; // обратите внимание на перестановку индексов
				}
			}
		}

	}




	//square matrix constructor
	derived_matrix(uint64_t size_diag) :matrix<T>(size_diag)  {}
	// # + задание OrderStorage 
	//constructor with dimensions
	derived_matrix(uint64_t colsize, uint64_t rowsize);
	// # + задание OrderStorage 
	//the default constructor
	derived_matrix();
	// # + задание OrderStorage 
	template<typename U>matrix(const derived_matrix<U>& other);
	// # + задание OrderStorage 
	derived_matrix(std::initializer_list<std::initializer_list<T>> init);
	// # проверить
	derived_matrix(matrix<T>&& other) noexcept;
	// # + задание OrderStorage (приведение к виду по умолчанию?)


	//DATA ACCESS

	//to index row access operator
	T* operator[](uint64_t index) const {
		if (OrderStorage == StorageMethod::row_major_order)
			return this->ptr + index * this->rowsize;
		else if (OrderStorage == StorageMethod::column_major_order)
			return this->ptr + index * this->colsize;
		else
			throw std::runtime_error("Unknown storage order");
	}	
	// # подумать над альтернативо выкидывани€ ошибок

	//ARITHMETIC OPERATORS
	// # выбрать способ хранени€ по умолчни€(наслудетс€ из коснтурктора)


	//r-value moving
	derived_matrix<T>& operator=(derived_matrix<T>&& other)
	// # + задание OrderStorage (приведение к виду по умолчанию?)



	//DATA ACCESS
	/////// ===========================
	/////// # исправить типы на новые
	/////// ===========================
	matrix<T>& set_output_mode(output_mode mode) { this->out_mode = mode; return *this; }

	matrix<T> get_column(uint64_t col) const;
	matrix<T> get_row(uint64_t row) const;

	//ARITHMETIC OPERATORS

	// the unary operator returns a matrix with inverse (multiplied by minus 1 ) elements
	matrix<T> operator-() const;
	// binary matrix addition
	matrix<T> operator+(const matrix<T>& other) const;
	// binary matrix subtraction
	matrix<T> operator-(const matrix<T>& other) const;
	// binary matrix multiplication
	matrix<T> operator*(const matrix<T>& other) const;
	// binary matrix multiplication by  an element from the field
	matrix<T> operator*(const T& other) const;
	// assignment operator overload
	matrix<T>& operator=(const matrix<T>& other);
	//r-value moving
	matrix<T>& operator=(matrix<T>&& other) noexcept {
		if (this != &other) {
			delete[] ptr;
			ptr = other.ptr;
			colsize = other.colsize;
			rowsize = other.rowsize;
			out_mode = other.out_mode;
			other.ptr = nullptr;
			other.colsize = 0;
			other.rowsize = 0;
		}
		return *this;
	}
	template <typename U>matrix<T>& operator=(const matrix<U>& other) {
		if (static_cast<const void*>(this) != static_cast<const void*>(&other)) {
			delete[] ptr;
			colsize = other.getcol();
			rowsize = other.getrow();
			allocateMemory();
			for (uint64_t i = 0; i < colsize; ++i) {
				for (uint64_t j = 0; j < rowsize; ++j) {
					(*this)[i][j] = static_cast<T>(other[i][j]);
				}
			}
		}
		return *this;
	}
	// binary matrix division(if not a singular matrix on the left)
	matrix<T> operator/(const matrix<T>& other) const;
	matrix<T> operator/(const T& other) const;

	// I/O OPERATIONS

	// overloading the output operator
	template<typename T>friend std::ostream& operator<<<>(std::ostream& out, const matrix<T>& p);
	// overloading the input operator
	template<typename T>friend std::istream& operator>><>(std::istream& in, matrix<T>& p);

	//SPECIAL METHODS
	// The p norm method for a matrix
	T norm(int p = 2) const;
	//return of the upper triangular matrix after transformations
	matrix to_uptrng()const;
	//bringing to the upper triangular view together with the "other" matrix
	matrix to_uptrng(matrix<T>& other)const;
	//matrix transposition
	matrix transpose()const;
	// finding the determinant if there is one
	T determinant() const;
	//return of the square matrix from 1 to the lower diagonal
	//  <-----S---->
	//  0 0  ..  0 1     |
	//  0 0  ..  1 0     |
	//  . ..   ...  .. ..     S                    
	//  0 1  ..  0 0     |
	//  1 0  ..  0 0     |
	matrix sqprediag(const uint64_t S)const;
	//return of the inverse matrix  
	matrix inverse_M()const;
	//The Cholesky decomposition
	matrix<T> cholesky() const;
	//LUP decomposition
	LUResult<T> LUP() const;
	//replacement by a matrix of zero T elements
	static matrix<T> zeros(uint64_t colsize, uint64_t rowsize);

	//replacement by a matrix of single T elements
	//  <---rows--->
	//  1 1  ..  1 1     |
	//  1 1  ..  1 1     |
	//  . ..   ...  .. ..    cols                   
	//  1 1  ..  1 1     |
	//  1 1  ..  1 1     |
	static matrix<T> ones(uint64_t colsize, uint64_t rowsize);
	//The identity matrix
	//  <-----S---->
	//  1 0  ..  0 0     |
	//  0 1  ..  0 0     |
	//  . ..   ...  .. ..     S                    
	//  0 0  ..  1 0     |
	//  0 0  ..  0 1     |
	static matrix<T> eye(uint64_t S);
	//Random matrix
	// min < R < max
	//  <---rows--->
	//  R R  ..  R R     |
	//  R R  ..  R R     |
	//  . ..   ...  .. ..    cols                   
	//  R R  ..  R R     |
	//  R R  ..  R R     |
	static matrix<T> random(size_t rows, size_t cols, T min, T max);

	//Random matrix
	// min < R < max
	//  <---rows--->
	//  R 0  ..  0 0     |
	//  0 R  ..  0 0     |
	//  . ..   ...  .. ..     cols                    
	//  0 0  ..  R 0     |
	//  0 0  ..  0 R     |
	static matrix<T> randomDiagonal(size_t n, T min, T max);

	// Constructed from vector x of length n, with m columns (if m=0 then m=n)
	// Matrix dimensions: n rows x m columns
	// Element at row i, column j: x[i]^j  [increasing powers]
	//  <--- m columns --->
	//  x0^0  x0^1  x0^2  ...  x0^(m-1)   | row 0
	//  x1^0  x1^1  x1^2  ...  x1^(m-1)   | row 1
	//  ...                                | ...
	//  xn^0  xn^1  xn^2  ...  xn^(m-1)   | row n-1
	static matrix<T> vander(const std::vector<T>& x, uint64_t m = 0);

	// (A'A)^-1 * A'
	matrix<T> left_pseudo_reverse()const;
	// A'*( A*A')^-1
	matrix<T> right_pseudo_reverse()const;

	// General pseudo-inverse matrix (automatic selection)
	matrix<T> pseudo_inverse()const;

	// Method for calculating the maximum eigenvalue with initial vector
	T max_eigenvalue(const matrix<T>& initial_vec, double epsilon = 1e-6, int max_iter = 1000) const;

	//
	//mat:
	//cols: 4, rows: 4
	// 68 -94  73 -32
	// -5  -7  32  39
	//-69 -79 -91  43
	// 81  40   6  49
	//sub:
	//cols: 3, rows: 3
	// -94  73 -32
	// -7  32  39
	// -79 -91  43
	//
	matrix<T> submatrix(uint64_t start_row, uint64_t start_col, uint64_t rows, uint64_t cols) const;

	//LOGIC OPERATIONS
	//Is a vector a row or a vector a column
	bool operator==(const matrix<T>& other) const {
		if (this->getcol() == other.getcol() && this->getrow() == other.getrow()) {

			for (size_t i = 0; i < this->getcol(); i++)
			{
				for (size_t j = 0; j < this->getrow(); j++)
				{
					if ((*this)[j][i] != other[j][i]) {
						return false;
					}
				}
			}
			return true;
		}
		else {
			return false;
		}
	}
};


