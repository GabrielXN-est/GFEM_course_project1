#include "lin_alg.h"

// Funções auxiliares
std::size_t Vector_aux_get_vec_size_f_Mat (const Matrix& m)
{
    if (m.mat[0].size() == 1)
        return m.mat.size();
    else if (m.mat.size() == 1)
        return m.mat[0].size();
    else
        throw std::invalid_argument("Matrix must have one row or one column to convert to Vector");
}

// Construtortes
Vector::Vector (Matrix m) : 
    vec(Vector_aux_get_vec_size_f_Mat (m))
{
    if (m.mat[0].size() == 1)
    {
        for (std::size_t i {0}; i < m.mat.size(); i++)
        {
            vec[i] = m.mat[i][0];
        }
    }
    else if (m.mat.size() == 1)
    {
        for (std::size_t i {0}; i < m.mat[0].size(); i++)
        {
            vec[i] = m.mat[0][i];
        }
    }
    else
        throw std::invalid_argument("Matrix must have one row or one column to convert to Vector");
}

// Getters
const std::size_t Vector::size() const
{
    return vec.size();
}

//operadores
double& Vector::operator[] (std::size_t index)
{
    if (index >= vec.size())
        throw std::out_of_range("Vector index out of range");
    return vec[index];
}

double Vector::get (std::size_t index) const
{
    if (index >= vec.size())
        throw std::out_of_range("Vector index out of range");
    return vec[index];
}

// operadores
// transpor o vetor
Vector Vector::T() const
{
    Vector transposed = *this;
    transposed.transposed = !transposed.transposed;
    return transposed;
}

void Vector::operator+= (const Vector& other)
{
    if (vec.size() != other.vec.size())
        throw std::invalid_argument("Vectors must have the same size for addition");
    for (std::size_t i {0}; i < vec.size(); i++)
    {
        vec[i] += other.vec[i];
    }
}

Vector Vector::operator+ (const Vector& other) const
{
    if (size() != other.size())
        throw std::invalid_argument("Vectors must have the same size for addition");

    Vector result(size());

    for (std::size_t i {0}; i < size(); i++)
    {
        result.vec[i] = get(i) + other.get(i);
    }
    return std::move(result);
}

Vector Vector::operator- () const
{
    Vector result(size());
    for (std::size_t i {0}; i < size(); i++)
        {result[i] = -get(i);}
    return std::move(result);
}

Matrix Vector::operator* (const double& v2)
{
    // multiplicação por escalar
    Matrix result(size(), 1);
    for (std::size_t i {0}; i < size(); i++)
    {
        result[i][0] = get(i) * v2;
    }
    return result;
}

Matrix Vector::operator* (const Matrix& v2)
{
    Matrix m1 {*this};
    return m1.operator*(v2);
}