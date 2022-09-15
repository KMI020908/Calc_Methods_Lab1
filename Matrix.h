#ifndef MATRIX_H
#define MATRIX_H

#include<vector>
#include<iostream>

template<typename Type>
class SquareMatrix {
private:
    std::vector<std::vector<Type>> m_data;

public:
    SquareMatrix() : m_data({}){}

    SquareMatrix(const size_t n){
        reset(n);
    }

    void reset(const size_t n) {
        m_data.resize(n);
        for (size_t i = 0; i < n; ++i){
            this->at(i).resize(n);
        }
    }

    Type getValue(const size_t& row, const size_t& col) const {
        return m_data.at(row).at(col);
    }

    int size() const {
        return m_data.size();
    }

    friend std::ostream& operator<<(std::ostream& out,const SquareMatrix& matrix) {
        for (int i = 0; i < matrix.m_data.size(); i++)
        {
            for (int j = 0; j < matrix.m_data.size(); j++)
            {
                out << matrix.at(i, j) << ' ';
            }
            std::cout << '\n';
        }
        return out;
    }
};

#endif