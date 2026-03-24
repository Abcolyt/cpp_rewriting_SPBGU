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
#include <immintrin.h>  // AVX-интринсики


// ============================================================================
// ==========================  CLASS Matrix2D  ================================
// ============================================================================

template<typename T>
class Matrix2D {
public:
	/// Формат хранения данных в памяти
	enum class StorageOrder {
		RowMajor,    // построчное хранение
		ColumnMajor  // по столбцовое хранение
	};

private:
	T* ptr;           ///< Указатель на данные матрицы
	uint64_t cols;    ///< Количество столбцов
	uint64_t rows;    ///< Количество строк
	StorageOrder storage;  ///< Формат хранения данных

	/// Выделение памяти под матрицу
	void allocateMemory();

public:
	// ==================== КОНСТРУКТОРЫ / ДЕСТРУКТОР ====================

	Matrix2D();
	Matrix2D(uint64_t cols, uint64_t rows, StorageOrder storage = StorageOrder::RowMajor);
	Matrix2D(uint64_t size, StorageOrder storage = StorageOrder::RowMajor);
	Matrix2D(std::initializer_list<std::initializer_list<T>> init, StorageOrder storage = StorageOrder::RowMajor);
	Matrix2D(const Matrix2D& other);
	Matrix2D(Matrix2D&& other) noexcept;
	~Matrix2D();

	// ==================== ДОСТУП К ЭЛЕМЕНТАМ ====================

	/// Доступ к элементу по индексам
	T& at(uint64_t row, uint64_t col);
	const T& at(uint64_t row, uint64_t col) const;

	// ==================== ОПЕРАТОРЫ ПРИСВАИВАНИЯ ====================

	Matrix2D& operator=(const Matrix2D& other);
	Matrix2D& operator=(Matrix2D&& other) noexcept;

	// ==================== ДОСТУП К ДАННЫМ ====================

	uint64_t getcol() const;
	uint64_t getrow() const;
	StorageOrder getStorageOrder() const;
	void setStorageOrder(StorageOrder order);

	/// Получить прямой доступ к данным (для оптимизированных операций)
	/// @param row индекс строки
	/// @param col индекс столбца
	/// @return указатель на данные в памяти (последовательно для данного формата хранения)
	/// @note Для RowMajor: возвращает указатель на начало строки row
	/// @note Для ColumnMajor: возвращает указатель на начало столбца col
	const double* getDataPtr(uint64_t row, uint64_t col) const;

	T* operator[](uint64_t row);
	const T* operator[](uint64_t row) const;
	T& operator()(uint64_t row, uint64_t col);
	const T& operator()(uint64_t row, uint64_t col) const;

	// ==================== АРИФМЕТИЧЕСКИЕ ОПЕРАТОРЫ ====================

	Matrix2D operator-() const;
	Matrix2D operator+(const Matrix2D& other) const;
	Matrix2D operator-(const Matrix2D& other) const;
	Matrix2D operator*(const Matrix2D& other) const;
	Matrix2D operator*(const T& scalar) const;
	Matrix2D operator/(const T& scalar) const;

	friend Matrix2D operator*(const T& scalar, const Matrix2D& mat) {
		return mat * scalar;
	}

	// ==================== СТАТИЧЕСКИЕ МЕТОДЫ ====================

	/// Создать случайную матрицу
	static Matrix2D random(uint64_t cols, uint64_t rows, T min, T max, StorageOrder storage = StorageOrder::RowMajor);

	// ==================== ОПЕРАТОРЫ СРАВНЕНИЯ ====================

	bool operator==(const Matrix2D& other) const;
	bool operator!=(const Matrix2D& other) const;
};


// ============================================================================
// ==================== РЕАЛИЗАЦИИ МЕТОДОВ ====================================
// ============================================================================

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
	ptr = new T[cols * rows]();
}

template<typename T>
Matrix2D<T>::Matrix2D()
	: ptr(nullptr), cols(0), rows(0), storage(StorageOrder::RowMajor) {}

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
			for (uint64_t i = 0; i < cols; ++i) {
				for (uint64_t j = 0; j < rows; ++j) {
					newPtr[j * cols + i] = ptr[i * rows + j];
				}
			}
		}
		else {
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

template<>
const double* Matrix2D<double>::getDataPtr(uint64_t row, uint64_t col) const {
	if (storage == StorageOrder::RowMajor) {
		// Для RowMajor: строки хранятся последовательно
		return ptr + row * cols;
	}
	else {
		// Для ColumnMajor: столбцы хранятся последовательно
		return ptr + col * rows;
	}
}

template<typename T>
T* Matrix2D<T>::operator[](uint64_t row) {
	if (storage == StorageOrder::RowMajor) {
		return ptr + row * cols;
	}
	else {
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
Matrix2D<T> Matrix2D<T>::operator-() const {
	Matrix2D result(cols, rows, storage);
	for (uint64_t i = 0; i < cols * rows; ++i) {
		result.ptr[i] = -ptr[i];
	}
	return result;
}

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
	if (cols != other.rows) {
		throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
	}
	Matrix2D result(other.cols, rows, storage);
	for (uint64_t i = 0; i < rows; ++i) {
		for (uint64_t j = 0; j < other.cols; ++j) {
			T sum = 0;
			for (uint64_t k = 0; k < cols; ++k) {
				sum += at(i, k) * other.at(k, j);
			}
			result.at(i, j) = sum;
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

// ============================================================================
// ==================== СВОБОДНЫЕ ФУНКЦИИ =====================================
// ============================================================================

/// Умножение матриц с явным указанием результата
/// @param A первая матрица (левый множитель)
/// @param B вторая матрица (правый множитель)
/// @param result матрица для сохранения результата (должна быть правильного размера: A.rows × B.cols)
template<typename T>
void multiply(const Matrix2D<T>* A, const Matrix2D<T>* B, Matrix2D<T>* result) {
	if (A == nullptr || B == nullptr || result == nullptr) {
		throw std::invalid_argument("Null pointer passed to multiply");
	}
	
	if (A->getcol() != B->getrow()) {
		throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
	}
	
	// Проверяем размер result
	if (result->getrow() != A->getrow() || result->getcol() != B->getcol()) {
		// Пересоздаём result правильного размера
		*result = Matrix2D<T>(B->getcol(), A->getrow(), A->getStorageOrder());
	}
	
	uint64_t rowsA = A->getrow();
	uint64_t colsA = A->getcol();
	uint64_t colsB = B->getcol();
	
	for (uint64_t i = 0; i < rowsA; ++i) {
		for (uint64_t j = 0; j < colsB; ++j) {
			T sum = 0;
			for (uint64_t k = 0; k < colsA; ++k) {
				sum += A->at(i, k) * B->at(k, j);
			}
			result->at(i, j) = sum;
		}
	}
}

// ============================================================================
// AVX-версия умножения матриц (только для double)
// ============================================================================

/// Умножение матриц с использованием AVX инструкций + FMA
/// @param A первая матрица (левый множитель), должна быть в RowMajor формате
/// @param B вторая матрица (правый множитель), должна быть в ColumnMajor формате
/// @param result матрица для сохранения результата
/// @note B должна храниться в ColumnMajor для эффективной загрузки столбцов AVX-регистром
void multiplyAVX(const Matrix2D<double>* A, const Matrix2D<double>* B, Matrix2D<double>* result) {
	if (A == nullptr || B == nullptr || result == nullptr) {
		throw std::invalid_argument("Null pointer passed to multiplyAVX");
	}
	
	if (A->getcol() != B->getrow()) {
		throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
	}
	
	// Проверяем форматы хранения
	if (A->getStorageOrder() != Matrix2D<double>::StorageOrder::RowMajor) {
		throw std::invalid_argument("Matrix A must be in RowMajor format for AVX multiplication");
	}
	if (B->getStorageOrder() != Matrix2D<double>::StorageOrder::ColumnMajor) {
		throw std::invalid_argument("Matrix B must be in ColumnMajor format for AVX multiplication");
	}
	
	// Проверяем размер result
	if (result->getrow() != A->getrow() || result->getcol() != B->getcol()) {
		*result = Matrix2D<double>(B->getcol(), A->getrow(), A->getStorageOrder());
	}
	
	const uint64_t M = A->getrow();      // строки A
	const uint64_t K = A->getcol();      // столбцы A = строки B
	const uint64_t N = B->getcol();      // столбцы B
	
	// Получаем прямые указатели на данные один раз!
	const double* dataA = A->getDataPtr(0, 0);  // начало данных A
	const double* dataB = B->getDataPtr(0, 0);  // начало данных B
	
	// AVX + FMA умножение: C[i][j] = sum_k(A[i][k] * B[k][j])
	for (uint64_t i = 0; i < M; ++i) {
		// Указатель на начало i-й строки A (RowMajor: строки последовательно)
		const double* rowA = dataA + i * K;
		
		for (uint64_t j = 0; j < N; ++j) {
			// Указатель на начало j-го столбца B (ColumnMajor: столбцы последовательно)
			const double* colB = dataB + j * K;
			
			// Инициализируем аккумулятор нулём
			__m256d vsum = _mm256_setzero_pd();
			
			// Обрабатываем по 4 элемента AVX-регистром с FMA
			uint64_t k = 0;
			for (; k + 4 <= K; k += 4) {
				// Загружаем 4 элемента из строки A и столбца B
				__m256d va = _mm256_loadu_pd(rowA + k);
				__m256d vb = _mm256_loadu_pd(colB + k);
				
				// FMA: vsum = va * vb + vsum (одна инструкция!)
				vsum = _mm256_fmadd_pd(va, vb, vsum);
			}
			
			// Горизонтальная сумма аккумулятора
			__m128d vlow = _mm256_castpd256_pd128(vsum);        // [a0, a1]
			__m128d vhigh = _mm256_extractf128_pd(vsum, 1);     // [a2, a3]
			__m128d vsum128 = _mm_add_pd(vlow, vhigh);          // [a0+a2, a1+a3]
			__m128d vshuf = _mm_shuffle_pd(vsum128, vsum128, 1); // [a1+a3, a0+a2]
			__m128d vfinal = _mm_add_sd(vsum128, vshuf);        // [a0+a1+a2+a3, ...]
			
			double sum = _mm_cvtsd_f64(vfinal);
			
			// Обрабатываем "хвост" — оставшиеся 1-3 элемента скалярно
			for (; k < K; ++k) {
				sum += rowA[k] * colB[k];
			}
			
			result->at(i, j) = sum;
		}
	}
}
