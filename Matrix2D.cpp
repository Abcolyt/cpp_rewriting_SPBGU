#pragma once
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <utility>
#include <random>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


// ============================================================================
// ==========================  CLASS Matrix2D  ================================
// ============================================================================
// Новый независимый класс матрицы с поддержкой двух форматов хранения:
// - RowMajor (построчное)
// - ColumnMajor (по столбцовое)
// ============================================================================

template<typename T>
class Matrix2D {
public:
	/// Формат хранения данных в памяти
	enum class StorageOrder {
		RowMajor,    // построчное хранение (как в базовом matrix<T>)
		ColumnMajor  // по столбцовое хранение
	};

private:
	T* ptr;           ///< Указатель на данные матрицы
	uint64_t cols;    ///< Количество столбцов
	uint64_t rows;    ///< Количество строк
	StorageOrder storage;  ///< Формат хранения данных

	/// Внутренний доступ к элементу по индексам (с учётом формата хранения)
	/// @param row индекс строки
	/// @param col индекс столбца
	/// @return ссылка на элемент матрицы
	T& at(uint64_t row, uint64_t col);
	const T& at(uint64_t row, uint64_t col) const;

	/// Выделение памяти под матрицу (с инициализацией нулями)
	void allocateMemory();

public:
	// ==================== КОНСТРУКТОРЫ / ДЕСТРУКТОР ====================

	/// Конструктор по умолчанию — создаёт пустую матрицу 0×0
	Matrix2D();

	/// Конструктор с указанием размеров и формата хранения
	/// @param cols количество столбцов
	/// @param rows количество строк
	/// @param storage формат хранения (по умолчанию RowMajor)
	Matrix2D(uint64_t cols, uint64_t rows, StorageOrder storage = StorageOrder::RowMajor);

	/// Конструктор квадратной матрицы
	/// @param size размер матрицы (количество строк и столбцов)
	/// @param storage формат хранения (по умолчанию RowMajor)
	Matrix2D(uint64_t size, StorageOrder storage = StorageOrder::RowMajor);

	/// Конструктор инициализации списком значений
	/// @param init список списков значений для инициализации
	/// @param storage формат хранения (по умолчанию RowMajor)
	Matrix2D(std::initializer_list<std::initializer_list<T>> init, StorageOrder storage = StorageOrder::RowMajor);

	/// Конструктор копирования — создаёт копию другой матрицы
	/// @param other матрица для копирования
	Matrix2D(const Matrix2D& other);

	/// Конструктор перемещения — перемещает данные из другой матрицы
	/// @param other матрица для перемещения
	Matrix2D(Matrix2D&& other) noexcept;

	/// Деструктор — освобождает выделенную память
	~Matrix2D();

	// ==================== ОПЕРАТОРЫ ПРИСВАИВАНИЯ ====================

	/// Оператор присваивания копированием
	/// @param other матрица для копирования
	/// @return ссылка на текущую матрицу
	Matrix2D& operator=(const Matrix2D& other);

	/// Оператор присваивания перемещением
	/// @param other матрица для перемещения
	/// @return ссылка на текущую матрицу
	Matrix2D& operator=(Matrix2D&& other) noexcept;

	// ==================== ДОСТУП К ДАННЫМ ====================

	/// Получить количество столбцов матрицы
	/// @return количество столбцов
	uint64_t getcol() const;

	/// Получить количество строк матрицы
	/// @return количество строк
	uint64_t getrow() const;

	/// Получить текущий формат хранения данных
	/// @return формат хранения (RowMajor или ColumnMajor)
	StorageOrder getStorageOrder() const;

	/// Установить новый формат хранения (с конвертацией данных)
	/// @param order новый формат хранения
	void setStorageOrder(StorageOrder order);

	/// Оператор доступа к строке по индексу
	/// @param row индекс строки
	/// @return указатель на начало строки
	/// @note Для ColumnMajor создаётся временная копия строки (неэффективно)
	T* operator[](uint64_t row);

	/// Оператор доступа к строке по индексу (const-версия)
	/// @param row индекс строки
	/// @return указатель на начало строки
	const T* operator[](uint64_t row) const;

	/// Оператор доступа к элементу по индексам (рекомендуется для ColumnMajor)
	/// @param row индекс строки
	/// @param col индекс столбца
	/// @return ссылка на элемент
	T& operator()(uint64_t row, uint64_t col);

	/// Оператор доступа к элементу по индексам (const-версия)
	/// @param row индекс строки
	/// @param col индекс столбца
	/// @return константная ссылка на элемент
	const T& operator()(uint64_t row, uint64_t col) const;

	/// Получить столбец матрицы как матрицу размера rows×1
	/// @param col индекс столбца
	/// @return матрица-столбец
	/// @throws std::out_of_range если индекс выходит за границы
	Matrix2D getColumn(uint64_t col) const;

	/// Получить строку матрицы как матрицу размера 1×cols
	/// @param row индекс строки
	/// @return матрица-строка
	/// @throws std::out_of_range если индекс выходит за границы
	Matrix2D getRow(uint64_t row) const;

	// ==================== УНАРНЫЕ ОПЕРАТОРЫ ====================

	/// Унарный минус — возвращает матрицу с противоположным знаком элементов
	/// @return новая матрица с элементами -ptr[i]
	Matrix2D operator-() const;

	// ==================== БИНАРНЫЕ АРИФМЕТИЧЕСКИЕ ОПЕРАТОРЫ ====================

	/// Сложение матриц
	/// @param other матрица для сложения
	/// @return новая матрица — результат сложения
	/// @throws std::invalid_argument если размеры матриц не совпадают
	Matrix2D operator+(const Matrix2D& other) const;

	/// Вычитание матриц
	/// @param other матрица для вычитания
	/// @return новая матрица — результат вычитания
	/// @throws std::invalid_argument если размеры матриц не совпадают
	Matrix2D operator-(const Matrix2D& other) const;

	/// Умножение матриц
	/// @param other матрица для умножения
	/// @return новая матрица — результат умножения
	/// @throws std::invalid_argument если число строк первой не равно числу столбцов второй
	Matrix2D operator*(const Matrix2D& other) const;

	/// Умножение матрицы на скаляр
	/// @param scalar скаляр для умножения
	/// @return новая матрица с умноженными элементами
	Matrix2D operator*(const T& scalar) const;

	/// Деление матрицы на скаляр
	/// @param scalar скаляр для деления
	/// @return новая матрица с разделёнными элементами
	/// @throws std::invalid_argument если скаляр равен нулю
	Matrix2D operator/(const T& scalar) const;

	/// Дружественный оператор для умножения скаляра на матрицу
	/// @param scalar скаляр для умножения
	/// @param mat матрица для умножения
	/// @return новая матрица — результат умножения
	friend Matrix2D operator*(const T& scalar, const Matrix2D& mat) {
		return mat * scalar;
	}

	// ==================== СПЕЦИАЛЬНЫЕ МЕТОДЫ ====================

	/// Транспонирование матрицы — строки становятся столбцами
	/// @return новая транспонированная матрица
	Matrix2D transpose() const;

	/// Вычисление нормы матрицы (Lp-норма)
	/// @param p порядок нормы (1, 2, INT_MAX для бесконечной нормы)
	/// @return значение нормы матрицы
	T norm(int p = 2) const;

	/// Вычисление определителя матрицы (только для квадратных)
	/// @return значение определителя
	/// @throws std::invalid_argument если матрица не квадратная
	T determinant() const;

	/// Вычисление обратной матрицы (только для невырожденных квадратных)
	/// @return новая обратная матрица
	/// @throws std::invalid_argument если матрица вырожденная или не квадратная
	Matrix2D inverse() const;

	/// Создать единичную матрицу заданного размера
	/// @param size размер матрицы
	/// @param storage формат хранения (по умолчанию RowMajor)
	/// @return единичная матрица с единицами на главной диагонали
	static Matrix2D eye(uint64_t size, StorageOrder storage = StorageOrder::RowMajor);

	/// Создать нулевую матрицу заданного размера
	/// @param cols количество столбцов
	/// @param rows количество строк
	/// @param storage формат хранения (по умолчанию RowMajor)
	/// @return матрица, заполненная нулями
	static Matrix2D zeros(uint64_t cols, uint64_t rows, StorageOrder storage = StorageOrder::RowMajor);

	/// Создать матрицу из единиц заданного размера
	/// @param cols количество столбцов
	/// @param rows количество строк
	/// @param storage формат хранения (по умолчанию RowMajor)
	/// @return матрица, заполненная единицами
	static Matrix2D ones(uint64_t cols, uint64_t rows, StorageOrder storage = StorageOrder::RowMajor);

	/// Создать случайную матрицу заданного размера
	/// @param cols количество столбцов
	/// @param rows количество строк
	/// @param min минимальное значение элемента
	/// @param max максимальное значение элемента
	/// @param storage формат хранения (по умолчанию RowMajor)
	/// @return матрица со случайными элементами
	static Matrix2D random(uint64_t cols, uint64_t rows, T min, T max, StorageOrder storage = StorageOrder::RowMajor);

	/// Создать случайную диагональную матрицу
	/// @param size размер матрицы
	/// @param min минимальное значение элемента
	/// @param max максимальное значение элемента
	/// @param storage формат хранения (по умолчанию RowMajor)
	/// @return диагональная матрица со случайными элементами на диагонали
	static Matrix2D randomDiagonal(uint64_t size, T min, T max, StorageOrder storage = StorageOrder::RowMajor);

	/// Создать матрицу Вандермонда из вектора значений
	/// @param x вектор значений
	/// @param m количество столбцов (степени от 0 до m-1), по умолчанию равно размеру x
	/// @param storage формат хранения (по умолчанию RowMajor)
	/// @return матрица Вандермонда размера n×m
	static Matrix2D vander(const std::vector<T>& x, uint64_t m = 0, StorageOrder storage = StorageOrder::RowMajor);

	/// Получить подматрицу — часть исходной матрицы
	/// @param startRow индекс начальной строки
	/// @param startCol индекс начального столбца
	/// @param numRows количество строк подматрицы
	/// @param numCols количество столбцов подматрицы
	/// @return новая матрица — подматрица
	/// @throws std::out_of_range если размеры выходят за границы
	Matrix2D submatrix(uint64_t startRow, uint64_t startCol, uint64_t numRows, uint64_t numCols) const;

	/// Поменять местами две строки матрицы
	/// @param i индекс первой строки
	/// @param j индекс второй строки
	void swapRows(size_t i, size_t j);

	// ==================== ОПЕРАТОРЫ СРАВНЕНИЯ ====================

	/// Проверка на равенство двух матриц
	/// @param other матрица для сравнения
	/// @return true если матрицы одинакового размера и все элементы равны
	bool operator==(const Matrix2D& other) const;

	/// Проверка на неравенство двух матриц
	/// @param other матрица для сравнения
	/// @return true если матрицы различаются размером или хотя бы одним элементом
	bool operator!=(const Matrix2D& other) const;

	// ==================== ВВОД / ВЫВОД ====================

	/// Оператор вывода матрицы в поток
	/// @param out выходной поток
	/// @param m матрица для вывода
	/// @return ссылка на поток
	template<typename U>
	friend std::ostream& operator<<(std::ostream& out, const Matrix2D<U>& m);

	/// Оператор ввода матрицы из потока
	/// @param in входной поток
	/// @param m матрица для ввода
	/// @return ссылка на поток
	template<typename U>
	friend std::istream& operator>>(std::istream& in, Matrix2D<U>& m);
};


// ============================================================================
// ==================== РЕАЛИЗАЦИИ МЕТОДОВ Matrix2D ===========================
// ============================================================================

// ------------------- Внутренние методы -------------------

template<typename T>
T& Matrix2D<T>::at(uint64_t row, uint64_t col) {
	if (storage == StorageOrder::RowMajor)
		return ptr[row * cols + col];
	else
		return ptr[col * rows + row];
}

template<typename T>
const T& Matrix2D<T>::at(uint64_t row, uint64_t col) const {
	if (storage == StorageOrder::RowMajor)
		return ptr[row * cols + col];
	else
		return ptr[col * rows + row];
}

template<typename T>
void Matrix2D<T>::allocateMemory() {
	delete[] ptr;
	ptr = new T[cols * rows]();  // инициализация нулями
}

// ------------------- Конструкторы / Деструктор -------------------

template<typename T>
Matrix2D<T>::Matrix2D()
	: ptr(nullptr), cols(0), rows(0), storage(StorageOrder::RowMajor) {
}

template<typename T>
Matrix2D<T>::Matrix2D(uint64_t cols, uint64_t rows, StorageOrder storage)
	: cols(cols), rows(rows), storage(storage) {
	ptr = new T[cols * rows]();
}

template<typename T>
Matrix2D<T>::Matrix2D(uint64_t size, StorageOrder storage)
	: cols(size), rows(size), storage(storage) {
	ptr = new T[size * size]();
}

template<typename T>
Matrix2D<T>::Matrix2D(std::initializer_list<std::initializer_list<T>> init, StorageOrder storage)
	: rows(init.size()), cols(0), storage(storage) {

	if (rows > 0) {
		cols = init.begin()->size();
		for (const auto& row : init) {
			if (row.size() != cols) {
				throw std::invalid_argument("All rows must have the same length");
			}
		}
	}

	ptr = new T[cols * rows]();

	uint64_t i = 0;
	for (const auto& row : init) {
		uint64_t j = 0;
		for (const auto& value : row) {
			at(i, j) = value;
			j++;
		}
		i++;
	}
}

template<typename T>
Matrix2D<T>::Matrix2D(const Matrix2D& other)
	: cols(other.cols), rows(other.rows), storage(other.storage) {
	ptr = new T[cols * rows];
	for (uint64_t i = 0; i < cols * rows; ++i) {
		ptr[i] = other.ptr[i];
	}
}

template<typename T>
Matrix2D<T>::Matrix2D(Matrix2D&& other) noexcept
	: ptr(other.ptr), cols(other.cols), rows(other.rows), storage(other.storage) {
	other.ptr = nullptr;
	other.cols = 0;
	other.rows = 0;
}

template<typename T>
Matrix2D<T>::~Matrix2D() {
	delete[] ptr;
}

// ------------------- Операторы присваивания -------------------

template<typename T>
Matrix2D<T>& Matrix2D<T>::operator=(const Matrix2D& other) {
	if (this != &other) {
		delete[] ptr;
		cols = other.cols;
		rows = other.rows;
		storage = other.storage;
		ptr = new T[cols * rows];
		for (uint64_t i = 0; i < cols * rows; ++i) {
			ptr[i] = other.ptr[i];
		}
	}
	return *this;
}

template<typename T>
Matrix2D<T>& Matrix2D<T>::operator=(Matrix2D&& other) noexcept {
	if (this != &other) {
		delete[] ptr;
		ptr = other.ptr;
		cols = other.cols;
		rows = other.rows;
		storage = other.storage;
		other.ptr = nullptr;
		other.cols = 0;
		other.rows = 0;
	}
	return *this;
}

// ------------------- Доступ к данным -------------------

template<typename T>
uint64_t Matrix2D<T>::getcol() const {
	return cols;
}

template<typename T>
uint64_t Matrix2D<T>::getrow() const {
	return rows;
}

template<typename T>
typename Matrix2D<T>::StorageOrder Matrix2D<T>::getStorageOrder() const {
	return storage;
}

template<typename T>
void Matrix2D<T>::setStorageOrder(StorageOrder order) {
	if (order != storage) {
		T* newPtr = new T[cols * rows];
		if (storage == StorageOrder::RowMajor) {
			// Было row-major, станет column-major
			for (uint64_t i = 0; i < cols; ++i) {
				for (uint64_t j = 0; j < rows; ++j) {
					newPtr[j * cols + i] = ptr[i * rows + j];
				}
			}
		}
		else {
			// Было column-major, станет row-major
			for (uint64_t i = 0; i < cols; ++i) {
				for (uint64_t j = 0; j < rows; ++j) {
					newPtr[i * rows + j] = ptr[j * cols + i];
				}
			}
		}
		delete[] ptr;
		ptr = newPtr;
		storage = order;
	}
}

template<typename T>
T* Matrix2D<T>::operator[](uint64_t row) {
	if (storage == StorageOrder::RowMajor) {
		return ptr + row * cols;
	}
	else {
		// Для column-major создаём временную строку
		static T* tempRow = nullptr;
		static uint64_t tempRowSize = 0;
		if (tempRowSize != cols) {
			delete[] tempRow;
			tempRow = new T[cols];
			tempRowSize = cols;
		}
		for (uint64_t c = 0; c < cols; ++c) {
			tempRow[c] = at(row, c);
		}
		return tempRow;
	}
}

template<typename T>
const T* Matrix2D<T>::operator[](uint64_t row) const {
	if (storage == StorageOrder::RowMajor) {
		return ptr + row * cols;
	}
	else {
		// Для column-major создаём временную строку
		static T* tempRow = nullptr;
		static uint64_t tempRowSize = 0;
		if (tempRowSize != cols) {
			delete[] tempRow;
			tempRow = new T[cols];
			tempRowSize = cols;
		}
		for (uint64_t c = 0; c < cols; ++c) {
			tempRow[c] = at(row, c);
		}
		return tempRow;
	}
}

template<typename T>
T& Matrix2D<T>::operator()(uint64_t row, uint64_t col) {
	return at(row, col);
}

template<typename T>
const T& Matrix2D<T>::operator()(uint64_t row, uint64_t col) const {
	return at(row, col);
}

template<typename T>
Matrix2D<T> Matrix2D<T>::getColumn(uint64_t col) const {
	if (col >= cols) {
		throw std::out_of_range("Column index out of bounds");
	}
	Matrix2D result(1, rows, storage);
	for (uint64_t row = 0; row < rows; ++row) {
		result.at(0, row) = at(row, col);
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::getRow(uint64_t row) const {
	if (row >= rows) {
		throw std::out_of_range("Row index out of bounds");
	}
	Matrix2D result(cols, 1, storage);
	for (uint64_t col = 0; col < cols; ++col) {
		result.at(col, 0) = at(row, col);
	}
	return result;
}

// ------------------- Унарные операторы -------------------

template<typename T>
Matrix2D<T> Matrix2D<T>::operator-() const {
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = -ptr[i];
	}
	return result;
}

// ------------------- Бинарные арифметические операторы -------------------

template<typename T>
Matrix2D<T> Matrix2D<T>::operator+(const Matrix2D& other) const {
	if (cols != other.cols || rows != other.rows) {
		throw std::invalid_argument("Matrix dimensions must match for addition");
	}
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = ptr[i] + other.ptr[i];
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator-(const Matrix2D& other) const {
	if (cols != other.cols || rows != other.rows) {
		throw std::invalid_argument("Matrix dimensions must match for subtraction");
	}
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = ptr[i] - other.ptr[i];
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator*(const Matrix2D& other) const {
	if (rows != other.cols) {
		throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
	}
	Matrix2D result(cols, other.rows, storage);
	for (uint64_t i = 0; i < result.cols; ++i) {
		for (uint64_t j = 0; j < result.rows; ++j) {
			T sum = 0;
			for (uint64_t k = 0; k < rows; ++k) {
				sum += at(k, i) * other.at(j, k);
			}
			result.at(j, i) = sum;
		}
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator*(const T& scalar) const {
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = ptr[i] * scalar;
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator/(const T& scalar) const {
	if (scalar == T(0)) {
		throw std::invalid_argument("Division by zero");
	}
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = ptr[i] / scalar;
	}
	return result;
}

// ------------------- Специальные методы -------------------

template<typename T>
Matrix2D<T> Matrix2D<T>::transpose() const {
	Matrix2D result(rows, cols, storage);
	for (uint64_t i = 0; i < cols; ++i) {
		for (uint64_t j = 0; j < rows; ++j) {
			result.at(j, i) = at(i, j);
		}
	}
	return result;
}

template<typename T>
T Matrix2D<T>::norm(int p) const {
	T result = 0;
	const uint64_t size = cols * rows;

	if (p == 1) {
		for (uint64_t i = 0; i < size; ++i) {
			result += std::abs(ptr[i]);
		}
	}
	else if (p == 2) {
		for (uint64_t i = 0; i < size; ++i) {
			result += ptr[i] * ptr[i];
		}
		result = std::sqrt(result);
	}
	else if (p == INT_MAX) {
		for (uint64_t i = 0; i < size; ++i) {
			T abs_val = std::abs(ptr[i]);
			if (abs_val > result) {
				result = abs_val;
			}
		}
	}
	else {
		for (uint64_t i = 0; i < size; ++i) {
			result += std::pow(std::abs(ptr[i]), p);
		}
		result = std::pow(result, 1.0 / p);
	}
	return result;
}

template<typename T>
T Matrix2D<T>::determinant() const {
	if (cols != rows) {
		throw std::invalid_argument("Determinant requires square matrix");
	}

	Matrix2D M = *this;
	T det = T(1);
	const size_t n = rows;

	if (n == 0) return T(1);

	for (size_t k = 0; k < n; ++k) {
		size_t maxRow = k;
		T maxVal = std::abs(M.at(k, k));
		for (size_t i = k + 1; i < n; ++i) {
			T val = std::abs(M.at(k, i));
			if (val > maxVal) {
				maxVal = val;
				maxRow = i;
			}
		}

		if (maxVal == T(0)) {
			return T(0);
		}

		if (maxRow != k) {
			for (size_t j = 0; j < n; ++j) {
				T temp = M.at(k, j);
				M.at(k, j) = M.at(maxRow, j);
				M.at(maxRow, j) = temp;
			}
			det = -det;
		}

		det *= M.at(k, k);

		for (size_t i = k + 1; i < n; ++i) {
			T factor = M.at(k, i) / M.at(k, k);
			for (size_t j = k; j < n; ++j) {
				M.at(j, i) -= factor * M.at(j, k);
			}
		}
	}

	return det;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::inverse() const {
	if (cols != rows) {
		throw std::invalid_argument("Inverse requires square matrix");
	}

	T det = determinant();
	if (det == T(0)) {
		throw std::invalid_argument("Matrix is singular, cannot invert");
	}

	Matrix2D augmented(cols, rows * 2, storage);

	for (uint64_t i = 0; i < cols; ++i) {
		for (uint64_t j = 0; j < rows; ++j) {
			augmented.at(j, i) = at(j, i);
		}
	}

	for (uint64_t i = 0; i < cols; ++i) {
		for (uint64_t j = 0; j < rows; ++j) {
			augmented.at(j, i + rows) = (i == j) ? T(1) : T(0);
		}
	}

	for (uint64_t k = 0; k < cols; ++k) {
		size_t maxRow = k;
		T maxVal = std::abs(augmented.at(k, k));
		for (size_t i = k + 1; i < rows; ++i) {
			T val = std::abs(augmented.at(k, i));
			if (val > maxVal) {
				maxVal = val;
				maxRow = i;
			}
		}

		if (maxVal == T(0)) {
			throw std::invalid_argument("Matrix is singular");
		}

		if (maxRow != k) {
			for (uint64_t j = 0; j < cols * 2; ++j) {
				T temp = augmented.at(k, j);
				augmented.at(k, j) = augmented.at(maxRow, j);
				augmented.at(maxRow, j) = temp;
			}
		}

		T pivot = augmented.at(k, k);
		for (uint64_t j = 0; j < cols * 2; ++j) {
			augmented.at(k, j) /= pivot;
		}

		for (uint64_t i = 0; i < rows; ++i) {
			if (i != k) {
				T factor = augmented.at(k, i);
				for (uint64_t j = 0; j < cols * 2; ++j) {
					augmented.at(i, j) -= factor * augmented.at(k, j);
				}
			}
		}
	}

	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols; ++i) {
		for (uint64_t j = 0; j < rows; ++j) {
			result.at(j, i) = augmented.at(j, i + rows);
		}
	}

	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::eye(uint64_t size, StorageOrder storage) {
	Matrix2D result(size, size, storage);
	for (uint64_t i = 0; i < size; ++i) {
		result.at(i, i) = T(1);
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::zeros(uint64_t cols, uint64_t rows, StorageOrder storage) {
	return Matrix2D(cols, rows, storage);
}

template<typename T>
Matrix2D<T> Matrix2D<T>::ones(uint64_t cols, uint64_t rows, StorageOrder storage) {
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = T(1);
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::random(uint64_t cols, uint64_t rows, T min, T max, StorageOrder storage) {
	Matrix2D result(cols, rows, storage);
	std::random_device rd;
	std::mt19937 gen(rd());

	if constexpr (std::is_integral_v<T>) {
		std::uniform_int_distribution<T> dist(min, max);
		for (uint64_t i = 0; i < cols * rows; ++i) {
			result.ptr[i] = dist(gen);
		}
	}
	else {
		std::uniform_real_distribution<T> dist(min, max);
		for (uint64_t i = 0; i < cols * rows; ++i) {
			result.ptr[i] = dist(gen);
		}
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::randomDiagonal(uint64_t size, T min, T max, StorageOrder storage) {
	Matrix2D result(size, size, storage);
	std::random_device rd;
	std::mt19937 gen(rd());

	if constexpr (std::is_integral_v<T>) {
		std::uniform_int_distribution<T> dist(min, max);
		for (uint64_t i = 0; i < size; ++i) {
			result.at(i, i) = dist(gen);
		}
	}
	else {
		std::uniform_real_distribution<T> dist(min, max);
		for (uint64_t i = 0; i < size; ++i) {
			result.at(i, i) = dist(gen);
		}
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::vander(const std::vector<T>& x, uint64_t m, StorageOrder storage) {
	uint64_t n = x.size();
	if (m == 0) m = n;

	Matrix2D result(n, m, storage);
	for (uint64_t i = 0; i < n; ++i) {
		result.at(i, 0) = T(1);
		for (uint64_t j = 1; j < m; ++j) {
			result.at(i, j) = result.at(i, j - 1) * x[i];
		}
	}
	return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::submatrix(uint64_t startRow, uint64_t startCol, uint64_t numRows, uint64_t numCols) const {
	if (startRow + numRows > rows || startCol + numCols > cols) {
		throw std::out_of_range("Submatrix dimensions exceed matrix bounds");
	}

	Matrix2D result(numCols, numRows, storage);
	for (uint64_t i = 0; i < numCols; ++i) {
		for (uint64_t j = 0; j < numRows; ++j) {
			result.at(j, i) = at(startRow + j, startCol + i);
		}
	}
	return result;
}

template<typename T>
void Matrix2D<T>::swapRows(size_t i, size_t j) {
	for (size_t col = 0; col < cols; ++col) {
		T temp = at(i, col);
		at(i, col) = at(j, col);
		at(j, col) = temp;
	}
}

// ------------------- Операторы сравнения -------------------

template<typename T>
bool Matrix2D<T>::operator==(const Matrix2D& other) const {
	if (cols != other.cols || rows != other.rows) {
		return false;
	}
	for (uint64_t i = 0; i < cols * rows; ++i) {
		if (ptr[i] != other.ptr[i]) {
			return false;
		}
	}
	return true;
}

template<typename T>
bool Matrix2D<T>::operator!=(const Matrix2D& other) const {
	return !(*this == other);
}

// ------------------- Ввод / Вывод -------------------

template<typename T>
std::ostream& operator<<(std::ostream& out, const Matrix2D<T>& m) {
	out << "Matrix2D[" << m.cols << "x" << m.rows << "] ("
		<< (m.storage == Matrix2D<T>::StorageOrder::RowMajor ? "RowMajor" : "ColumnMajor")
		<< ")\n";
	for (uint64_t i = 0; i < m.rows; ++i) {
		for (uint64_t j = 0; j < m.cols; ++j) {
			out << std::setw(10) << m.at(i, j) << " ";
		}
		out << "\n";
	}
	return out;
}

template<typename T>
std::istream& operator>>(std::istream& in, Matrix2D<T>& m) {
	uint64_t cols, rows;
	std::cout << "Enter number of columns: ";
	in >> cols;
	std::cout << "Enter number of rows: ";
	in >> rows;

	m = Matrix2D<T>(cols, rows);
	for (uint64_t i = 0; i < rows; ++i) {
		for (uint64_t j = 0; j < cols; ++j) {
			std::cout << "[" << i << "][" << j << "]=";
			in >> m.at(i, j);
		}
	}
	return in;
}

// ------------------- Свободные операторы -------------------

template<typename T>
Matrix2D<T> operator*(const T& scalar, const Matrix2D<T>& mat) {
	return mat * scalar;
}
