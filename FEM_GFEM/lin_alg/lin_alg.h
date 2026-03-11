#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <vector>
#include <fstream>

class Vector;

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
    Matrix (Vector v); // estudar sintaxe
    Matrix () {}

    //setters
    void clear() {mat.clear();}

    //getters
    std::vector<double>& operator[] (std::size_t index);

    //operadores
    void operator+= (const Matrix& other);

    Matrix operator*(const Matrix& other) const;

    Matrix operator*(const Vector& other) const;

    Matrix operator*(double other) const;
    
    // extraí valor de matriz 1x1
    double determinant();
};

class Vector
{
    public:
    std::vector<double> vec {};

    // se transposed = true: vetor linha
    // se transposed = false: vetor coluna
    bool transposed {false};

    // Contrutores
    Vector (std::size_t size) : vec(static_cast<int>(size)){}
    Vector (std::vector<double> v) : vec{v}{}
    Vector (){}
    Vector (Matrix m);

    //Setter
    void clear() {vec.clear();}

    //Getters
    const std::size_t size() const;

    double& operator[] (std::size_t index);

    double get (std::size_t index) const;
    
    //operadores
    // transpor o vetor
    Vector T() const;

    void operator+= (const Vector& other);

    Vector operator+ (const Vector& other) const;

    Vector operator- (const Vector& other) const {return operator+(-other);}

    Vector operator- () const;

    Matrix operator* (const double& v2);
    Matrix operator* (const Matrix& v2);

    // multiplicação de vetores (classes filhas de Vector)
    template <typename T2>
    std::enable_if_t<std::is_base_of_v<Vector, T2>, Matrix> operator* (const T2& v2)
    {
        // vetor linha por coluna
        if (transposed && !v2.transposed)
        {
            if (size() != v2.size())
                throw std::invalid_argument("Incompatible vector sizes for multiplication");
            double result {0.0};
            for (std::size_t i {0}; i < size(); i++)
            {
                result += get(i) * v2.get(i);
            }
            return {std::vector<std::vector<double>>{{{result}}}};
        }
        // vetor coluna por linha
        else if (!transposed && v2.transposed)
        {
            if (size() != v2.size())
                throw std::invalid_argument("Incompatible vector sizes for multiplication");
            Matrix result(size(), size());
            for (std::size_t i {0}; i < size(); i++)
            {
                for (std::size_t j {0}; j < size(); j++)
                {
                    result[i][j]= get(i) * v2.get(j);
                }
            }
            return result;
        }
        else
        {
            throw std::invalid_argument("Invalid vector orientations for multiplication");
        }
    }
};

Matrix operator* (double n, Matrix&& m);

Vector Gauss_elimination(Matrix K, Vector F);

#endif