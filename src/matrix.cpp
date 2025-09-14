#include"matrix.h"

vec::vec(const vec& other)
{
	this->x = other.x;
	this->y = other.y;
	this->z = other.z;
}

vec& vec::operator = (const vec& other)
{
	this->x = other.x;
	this->y = other.y;
	this->z = other.z;
	return *this;
}

double vec::size() const
{
	return sqrt(x * x + y * y + z * z);
}

M::M(int row, int column) : row(row), column(column)
{//первое число - строка, вторая столбец
	arr = new double* [row];
	for (int i = 0; i < row; ++i)
	{
		arr[i] = new double[column]();
	}
}

M::~M()
{
	for (int i = 0; i < row; ++i)
	{
		delete[] arr[i];
	}
	delete[] arr;
}

M::M(const M& other)
{
	row = other.row;
	column = other.column;
	arr = new double* [row];
	for (int i = 0; i < row; ++i)
	{
		arr[i] = new double[column];
	}
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			arr[i][j] = other.arr[i][j];
		}
	}
}

// double* M::operator[] (int i)
// {
// 	return arr[i];
// }

double* M::operator[] (int i) const
{
	return arr[i];
}

M& M::operator = (const M& other)
{
	for (int i = 0; i < row; ++i)
	{
		delete[] arr[i];
	}
	delete[] arr;

	row = other.row;
	column = other.column;
	arr = new double* [row];
	for (int i = 0; i < row; ++i)
	{
		arr[i] = new double[column];
	}
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			arr[i][j] = other.arr[i][j];
		}
	}
	return *this;
}

M M::operator* (const M& other)
{
	if (this->column != other.row)
	{
        throw "column A != row B";
	}
	M C(this->row, other.column);
	for (int i = 0; i < C.row; ++i)
	{
		for (int j = 0; j < C.column; ++j)
		{
			for (int r = 0; r < this->column; ++r)
			{
				C[i][j] = C[i][j] + (this->arr[i][r] * other.arr[r][j]);
			}
		}
	}
	return C;
}

void M::T()
{
	M C(this->column, this->row);
	for (int i = 0; i < C.row; ++i)
	{
		for (int j = 0; j < C.column; ++j)
		{
			C[i][j] = this->arr[j][i];
		}
	}
	*this = C;
}

vec M::operator*(const vec& other)
{
	vec c(0, 0, 0);
	c.x = this->arr[0][0] * other.x + this->arr[0][1] * other.y + this->arr[0][2] * other.z;
	c.y = this->arr[1][0] * other.x + this->arr[1][1] * other.y + this->arr[1][2] * other.z;
	c.z = this->arr[2][0] * other.x + this->arr[2][1] * other.y + this->arr[2][2] * other.z;
	return c;
}

int M::Row() const
{
	return row;
}

int M::Column() const
{
	return column;
}

MX::MX(double fi) : M(3, 3)
{
	arr[0][0] = 1;
	arr[0][1] = 0;
	arr[0][2] = 0;
	arr[1][0] = 0;
	arr[1][1] = cos(fi * PI / 180);
	arr[1][2] = sin(fi * PI / 180);
	arr[2][0] = 0;
	arr[2][1] = -sin(fi * PI / 180);
	arr[2][2] = cos(fi * PI / 180);
}

MY::MY(double fi) : M(3, 3)
{
	arr[0][0] = cos(fi * PI / 180);
	arr[0][1] = 0;
	arr[0][2] = -sin(fi * PI / 180);
	arr[1][0] = 0;
	arr[1][1] = 1;
	arr[1][2] = 0;
	arr[2][0] = sin(fi * PI / 180);
	arr[2][1] = 0;
	arr[2][2] = cos(fi * PI / 180);
}

MZ::MZ(double fi) : M(3, 3)
{
	arr[0][0] = cos(fi * PI / 180);
	arr[0][1] = sin(fi * PI / 180);
	arr[0][2] = 0;
	arr[1][0] = -sin(fi * PI / 180);;
	arr[1][1] = cos(fi * PI / 180);
	arr[1][2] = 0;
	arr[2][0] = 0;
	arr[2][1] = 0;
	arr[2][2] = 1;
}


void showM(const M& a)
{
	for (int i = 0; i < a.Row(); ++i)
	{
		for (int j = 0; j < a.Column(); ++j)
		{
			std::cout << a[i][j] << "\t\t";
		}
		std::cout << std::endl;
	}
}

void showV(const vec& a)
{
	std::cout << a.x << std::endl << a.y << std::endl << a.z << std::endl;
}
