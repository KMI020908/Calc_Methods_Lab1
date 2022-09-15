#ifndef MATRIX_H
#define MATRIX_H

#include<vector>
#include<iostream>

template<typename Type>
class SquareMatrix {
private:
    std::vector<std::vector<Type>> m_data;

public:
    SquareMatrix() : m_data({}) {}

    SquareMatrix(const size_t& n) {
        Reset(n);
    }

    void Reset(const size_t& n) {
        m_data.resize(n);
        for (size_t i = 0; i < n; ++i) {
            m_data.at(i).resize(n);
        }
    }

    int At(const size_t& row, const size_t& col) const {
        return m_data.at(row).at(col);
    }

    double& At(const size_t& row, const size_t& col) {
        return m_data.at(row).at(col);
    }

    int Size() const {
        return m_data.size();
    }

    friend bool SymmetricMatrixQ(const SquareMatrix& matrix)
    {
        for (int i = 0; i < matrix.m_data.size(); i++)
        {
            for (int j = 0; j < matrix.m_data.size(); j++)
                if (matrix.At(i, j) != matrix.At(j, i))
                    return false;
        }
        return true;
    }


    friend std::ostream& operator<<(std::ostream& out,const SquareMatrix& matrix) {
        for (int i = 0; i < matrix.m_data.size(); i++)
        {
            for (int j = 0; j < matrix.m_data.size(); j++)
            {
                out << matrix.At(i, j) << ' ';
            }
            std::cout << '\n';
        }
        return out;
    }
};

#endif