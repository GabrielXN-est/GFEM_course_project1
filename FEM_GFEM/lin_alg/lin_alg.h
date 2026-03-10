#ifndef LIN_ALH_H
#define LIN_ALH_H

#include <vector>
#include <fstream>

// vetores e matrizes
class Matrix
{
    public:
    std::vector<std::vector<double>> mat {};

    // Contrutores
    Matrix (std::size_t lin, std::size_t col) : mat(lin, std::vector<double>(col, 0.0))
    {}
    Matrix (std::vector<std::vector<double>> M) : mat{M}
    {}

    Matrix () {}

    //operadores
    std::vector<double>& operator[] (std::size_t index)
    {
        if (index >= mat.size())
            throw std::out_of_range("Vector index out of range");
        return mat[index];
    }

    void operator+= (const Matrix& other)
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

    Matrix operator*(const Matrix& other) const 
    {
        size_t rows = this->mat.size();
        size_t cols = other.mat[0].size();  // Assumindo matrizes retangulares
        size_t inner = this->mat[0].size();  // Dimensão interna
        Matrix result(rows, cols);  // Construtor Matrix(rows, cols) necessário
        
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

    double determinant()
    {
        if (mat.size() == 1 && mat[0].size() == 1)
            return mat[0][0];
        else
            throw std::logic_error("Not implemented for this matrix size");
    }
};

class Vector
{
    public:
    // se transposed = true: vetor linha
    // se transposed = false: vetor coluna
    std::vector<double> vec {};
    bool transposed {false};

    double& operator[] (std::size_t index)
    {
        if (index >= vec.size())
            throw std::out_of_range("Vector index out of range");
        return vec[index];
    }

    const std::size_t size() const
    {
        return vec.size();
    }

    // transpor o vetor
    Vector T() const
    {
        Vector transposed = *this;
        transposed.transposed = !transposed.transposed;
        return transposed;
    }

    Vector (std::size_t size) : vec(static_cast<int>(size)){}
    Vector (std::vector<double> v) : vec{v}{}
    Vector (){}

    Vector (Matrix m) : vec(m.mat.size())
    {
        if (m.mat[0].size() != 1)
            throw std::invalid_argument("Matrix must have one column to convert to Vector");
        for (std::size_t i {0}; i < m.mat.size(); i++)
        {
            vec[i] = m.mat[i][0];
        }
    }

    void operator+= (const Vector& other)
    {
        if (vec.size() != other.vec.size())
            throw std::invalid_argument("Vectors must have the same size for addition");
        for (std::size_t i {0}; i < vec.size(); i++)
        {
            vec[i] += other.vec[i];
        }
    }

    template <typename T1, typename T2>
    friend Matrix operator* (T1& v1, T2& v2)
    {
        if constexpr (std::is_base_of_v<Vector, T2>)
        {
            // vetor linha por coluna
            if (v1.transposed && !v2.transposed)
            {
                if (v1.size() != v2.size())
                    throw std::invalid_argument("Incompatible vector sizes for multiplication");
                double result {0.0};
                for (std::size_t i {0}; i < v1.size(); i++)
                {
                    result += v1[i] * v2[i];
                }
                return {std::vector<std::vector<double>>{{{result}}}};
            }

            // vetor coluna por linha
            else if (!v1.transposed && v2.transposed)
            {
                std::size_t size {v1.size()};
                if (size != v2.size())
                    throw std::invalid_argument("Incompatible vector sizes for multiplication");
                Matrix result(size, size);
                for (std::size_t i {0}; i < size; i++)
                {
                    for (std::size_t j {0}; j < size; j++)
                    {
                        result[i][j]= v1[i] * v2[j];
                    }
                }
                return result;
            }
            else
            {
                throw std::invalid_argument("Invalid vector orientations for multiplication");
            }
        }
        else if constexpr (std::is_floating_point_v<T2>)
        {
            // multiplicação por escalar
            Matrix result(v1.size(), 1);
            for (std::size_t i {0}; i < v1.size(); i++)
            {
                result[i] = v1[i] * v2;
            }
            return result;
        }
        else
        {
            throw std::invalid_argument("Unsupported type for multiplication");
        }
    }

    template <typename T1, typename T2>
    friend Matrix operator* (T1& v1, T2& v2)
    {
        if constexpr (std::is_base_of_v<Vector, T2>)
        {
            // vetor linha por coluna
            if (v1.transposed && !v2.transposed)
            {
                if (v1.size() != v2.size())
                    throw std::invalid_argument("Incompatible vector sizes for multiplication");
                double result {0.0};
                for (std::size_t i {0}; i < v1.size(); i++)
                {
                    result += v1[i] * v2[i];
                }
                return {std::vector<std::vector<double>>{{{result}}}};
            }

            // vetor coluna por linha
            else if (!v1.transposed && v2.transposed)
            {
                std::size_t size {v1.size()};
                if (size != v2.size())
                    throw std::invalid_argument("Incompatible vector sizes for multiplication");
                Matrix result(size, size);
                for (std::size_t i {0}; i < size; i++)
                {
                    for (std::size_t j {0}; j < size; j++)
                    {
                        result[i][j]= v1[i] * v2[j];
                    }
                }
                return result;
            }
            else
            {
                throw std::invalid_argument("Invalid vector orientations for multiplication");
            }
        }
        else if constexpr (std::is_floating_point_v<T2>)
        {
            // multiplicação por escalar
            Matrix result(v1.size(), 1);
            for (std::size_t i {0}; i < v1.size(); i++)
            {
                result[i] = v1[i] * v2;
            }
            return result;
        }
        else
        {
            throw std::invalid_argument("Unsupported type for multiplication");
        }
    }

    template <typename T1, typename T2>
    friend Matrix operator* (T1& v1, T2&& v2)
    {
        if constexpr (std::is_base_of_v<Vector, T2>)
        {
            // vetor linha por coluna
            if (v1.transposed && !v2.transposed)
            {
                if (v1.size() != v2.size())
                    throw std::invalid_argument("Incompatible vector sizes for multiplication");
                double result {0.0};
                for (std::size_t i {0}; i < v1.size(); i++)
                {
                    result += v1[i] * v2[i];
                }
                return {std::vector<std::vector<double>>{{{result}}}};
            }

            // vetor coluna por linha
            else if (!v1.transposed && v2.transposed)
            {
                std::size_t size {v1.size()};
                if (size != v2.size())
                    throw std::invalid_argument("Incompatible vector sizes for multiplication");
                Matrix result(size, size);
                for (std::size_t i {0}; i < size; i++)
                {
                    for (std::size_t j {0}; j < size; j++)
                    {
                        result[i][j]= v1[i] * v2[j];
                    }
                }
                return result;
            }
            else
            {
                throw std::invalid_argument("Invalid vector orientations for multiplication");
            }
        }
        else if constexpr (std::is_floating_point_v<T2>)
        {
            // multiplicação por escalar
            Matrix result(v1.size(), 1);
            for (std::size_t i {0}; i < v1.size(); i++)
            {
                result[i][0] = v1[i] * v2;
            }
            return result;
        }
        else
        {
            throw std::invalid_argument("Unsupported type for multiplication");
        }
    }

    /*  
    Vector& operator*= (double n)
    {
        for (double& i: vec)
        {
            i *= n;
        }
        return *this;
    }*/

    Vector operator* (double n)
    {
        Vector result (size());
        for (std::size_t i {0}; i < size(); i++)
        {
            result[i] = vec[i] * n;
        }
        return result;
    }
};

Matrix operator* (double n, Matrix&& m);

#endif