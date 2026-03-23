#include "lin_alg.h"

// Matrix
Matrix::Matrix (Vector v) : mat {((!v.transposed)? 
    std::vector<std::vector<double>>(v.size(), std::vector<double>(1, 0.))
    : std::vector<std::vector<double>>(1, std::vector<double>(v.size(), 0.)))}
{
    if (!v.transposed)
        for (std::size_t i {0}; i < v.size(); i++)
        {
            mat[i][0] = v.get(i);
        }
    else
        for (std::size_t i {0}; i < v.size(); i++)
        {
            mat[0][i] = v.get(i);
        }
}

std::vector<double>& Matrix::operator[] (std::size_t index)
{
    if (index >= mat.size())
        throw std::out_of_range("Vector index out of range");
    return mat[index];
}

void Matrix::operator+= (const Matrix& other)
{
    if (mat.size() != other.mat.size() || mat[0].size() != other.mat[0].size())
        throw std::invalid_argument("Matrices must have the same dimensions for addition");
    for (std::size_t i {0}; i < mat.size(); i++)
    {
        for (std::size_t j {0}; j < mat[0].size(); j++)
        {
            mat[i][j] += other.mat[i][j];
        }
    }
}

Matrix Matrix::operator+ (const Matrix& other) const
{
    if (mat.size() != other.mat.size() || mat[0].size() != other.mat[0].size())
        throw std::invalid_argument("Matrices must have the same dimensions for addition");
    
    Matrix result(mat.size(), mat[0].size());
    
    for (std::size_t i {0}; i < mat.size(); i++)
    {
        for (std::size_t j {0}; j < mat[0].size(); j++)
            {result.mat[i][j] = mat[i][j] + other.mat[i][j];}
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const 
{
    size_t rows = this->mat.size();
    size_t cols = other.mat[0].size();  // Assumindo matrizes retangulares
    size_t inner = this->mat[0].size();  // Dimensão interna
    Matrix result(rows, cols);  // Construtor Matrix(rows, cols) necessário
    
    if (inner != other.mat.size())
        throw std::invalid_argument("Incompatible matrix sizes for multiplication");

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) 
        {
            for (size_t k = 0; k < inner; ++k)
            {
                result.mat[i][j] += this->mat[i][k] * other.mat[k][j];
            }
        }
    }
    return result;
}

Matrix Matrix::operator*(const Vector& other) const 
{
    size_t rows = this->mat.size();
    size_t inner = this->mat[0].size(); 
    size_t cols {};
    if (!other.transposed)
    {
        cols = 1;  // Assumindo matrizes retangulares
        if (inner != other.size())
            throw std::invalid_argument("Incompatible matrix sizes for multiplication");
        Matrix result(rows, cols);  // Construtor Matrix(rows, cols) necessário
    }
    else
    {
        cols = other.size();  // Assumindo matrizes retangulares
        if (inner != 1)
            throw std::invalid_argument("Incompatible matrix sizes for multiplication");
    }
    
    Matrix result(rows, cols);  // Construtor Matrix(rows, cols) necessário

    for (size_t i = 0; i < rows; ++i) 
    {
        for (size_t j = 0; j < other.size(); ++j)
        {
            if (!other.transposed)
                    result.mat[i][0] += this->mat[i][j] * other.vec[j];
            else
                result.mat[i][j] += this->mat[i][0] * other.vec[j];
        }
    }

    return result;
}

Matrix Matrix::operator*(double other) const 
{
    Matrix result(this->mat.size(), this->mat[0].size());
    for (std::size_t i {0}; i < mat.size(); i++)
    {
        for (std::size_t j {0}; j < mat[0].size(); j++)
        {
            result.mat[i][j] = this->mat[i][j] * other;
        }
    }
    return result;
}

double Matrix::determinant()
{
    if (mat.size() !=  mat[0].size())
        throw std::logic_error("Matriz não quadrada");
    
    if (mat.size() == 1)
        return mat[0][0];
    else
    {
        double det {0};
        Matrix submat(mat.size() - 1, mat[0].size() - 1);
        for (std::size_t j {0}; j < mat[0].size(); j++)
        {
            if (mat[0][j] != 0)
            {
                for (std::size_t i {1}; i < mat.size(); i++)
                {
                        for (std::size_t k {0}; k < mat[0].size(); k++)
                        {
                            if (k < j)
                                submat[i-1][k] = mat[i][k];
                            else if (k > j)
                                submat[i-1][k-1] = mat[i][k];
                        }
                }
                det += (j % 2 == 0 ? 1 : -1) * mat[0][j] * submat.determinant();
            }
        }
        return det;
    }
}

Matrix operator* (double n, Matrix&& m)
{
    return std::move(m*n);
}