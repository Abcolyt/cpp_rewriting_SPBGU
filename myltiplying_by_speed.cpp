#include <../file_h/matrix.h>

template <typename T>
class derived_matrix : public matrix<T> {
private:
	enum class StorageMethod {
		None,									//не хранит

		One_dimensional_array_of_strings = 0,	//по строкам 
		one_dimensional_array_of_columns = 1,	//по столбцам

		row_major_order = 0,					//по строкам 
		column_major_oreder	=1		    		//по столбцам
	};
	StorageMethod OrderStorage = StorageMethod::None;
	StorageMethod GetOrderStorage()const {
		return OrderStorage;
	}
	void SetOrderStorage(StorageMethod OrderStorage)const{
		this->OrderStorage = OrderStorage;
	}
public:
	//CONSTRUCTORS\DESTRUCTORS
	//the copying constructor
	matrix(const matrix<T>& mtrx);
	//# + транспонирование (в конец?)
	//square matrix constructor
	matrix(uint64_t size_diag) :matrix(size_diag, size_diag) {}
	// # + задание OrderStorage 
	//constructor with dimensions
	matrix(uint64_t colsize, uint64_t rowsize);
	// # + задание OrderStorage 
	//the default constructor
	matrix();
	// # + задание OrderStorage 
	template <typename U>matrix(const matrix<U>& other);
	// # + задание OrderStorage 
	matrix(std::initializer_list<std::initializer_list<T>> init);
	// # проверить
	matrix(matrix<T>&& other) noexcept;
	// # + задание OrderStorage (приведение к виду по умолчанию?)


	//DATA ACCESS

	//to index row access operator
	T* operator[](const uint64_t row_index) const { return ptr + row_index * rowsize; }
	//№ + учет типа хранения

	//ARITHMETIC OPERATORS
	// # выбрать способ хранения по умолчния(наслудется из коснтурктора)

	//r-value moving
	matrix<T>& operator=(matrix<T>&& other)
	// # + задание OrderStorage (приведение к виду по умолчанию?)

};


