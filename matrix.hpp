// Copyright 2018 Your Name <your_email>

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <math.h>
#include <type_traits>
#include <limits>

template <class T>
class Matrix {
  static_assert(std::is_arithmetic<T>::value, "Not arithmetic type");
  int rows, columns;
  T** m;

 public:
  Matrix() {
    rows = 0;
    columns = 0;
    m = nullptr;
  }
  Matrix(int rows, int columns) {
    this->rows = rows;
    this->columns = columns;
    m = new T*[rows];
    for (int i = 0; i < rows; i++) {
      m[i] = new T[columns];
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m[i][j] = 0;  //что-то не так с обнулением элементов
      }
    }
  }

  Matrix(const Matrix& c) {
    rows = c.rows;
    columns = c.columns;
    m = new T*[rows];
    for (int i = 0; i < rows; i++) {
      m[i] = new T[columns];
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m[i][j] = c[i][j];
      }
    }
  }

  T* operator[](int index) const { return m[index]; }

  bool operator!=(const Matrix& N) const { return !(*this == N); }

  Matrix& operator=(const Matrix& c) {
    for (int i = 0; i < rows; i++) {
      delete[] m[i];
    }
    delete[] m;

    rows = c.rows;
    columns = c.columns;
    m = new T*[rows];
    for (int i = 0; i < rows; i++) {
      m[i] = new T[columns];
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m[i][j] = c[i][j];
      }
    }

    return *this;
  }

  bool operator==(const Matrix& E) const {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        if (std::is_floating_point<T>::value) {
          if (abs(m[i][j] - E.m[i][j]) > std::numeric_limits<T>::epsilon())
            return false;
        } else {
          if (m[i][j] != E[i][j]) return false;
        }
      }
    }
    return true;
  }

  Matrix operator+(const Matrix& R) const {  // R- правый операнд
    if (columns != R.columns || rows != R.rows) {
      Matrix<T> error;
      return error;
    }
    Matrix<T> sum(rows, columns);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        sum.m[i][j] = m[i][j] + R.m[i][j];
      }
    }
    return sum;
  }
  Matrix operator-(const Matrix& r) const {  // r-правый операнд
    if (columns != r.columns || rows != r.rows) {
      Matrix<T> error;
      return error;
    }
    Matrix<T> raz(rows, columns);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        raz.m[i][j] = m[i][j] - r.m[i][j];
      }
    }
    return raz;
  }
  Matrix operator*(const Matrix& F) {
    if (columns != F.rows) {
      Matrix<T> error;
      return error;
    }
    Matrix<T> proiz(rows, F.columns);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < F.columns; j++) {
        for (int k = 0; k < columns; k++) proiz[i][j] += m[i][k] * F.m[k][j];
      }
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < F.columns; j++) {
        if (abs(proiz[i][j]) < 10 * std::numeric_limits<double>::epsilon())
          proiz[i][j] = 0;
        if (abs(abs(proiz[i][j]) - round(proiz[i][j])) <
            10 * std::numeric_limits<double>::epsilon())
          proiz[i][j] = round(proiz[i][j]);
      }
    }
    return proiz;
  }
  Matrix Inverse() const {  //переделать операторы != и == для типа double
    Matrix<T> A(rows, 2 * rows);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < rows; j++) {
        A.m[i][j] = m[i][j];
      }
    }
    for (int i = 0; i < rows; i++) {
      A.m[i][rows + i] = 1;
    }

    for (int i = 0; i < rows; i++) {
      double maxEl = abs(A.m[i][i]);
      int maxRow = i;
      for (int k = i + 1; k < rows; k++) {
        if (abs(A.m[k][i]) > maxEl) {
          maxEl = A.m[k][i];
          maxRow = k;
        }
      }

      for (int k = i; k < 2 * rows; k++) {
        double tmp = A.m[maxRow][k];
        A.m[maxRow][k] = A.m[i][k];
        A.m[i][k] = tmp;
      }

      for (int k = i + 1; k < rows; k++) {
        double c = -A.m[k][i] / A.m[i][i];
        for (int j = i; j < 2 * rows; j++) {
          if (i == j) {
            A.m[k][j] = 0;
          } else {
            A.m[k][j] += c * A.m[i][j];
          }
        }
      }
    }

    for (int i = rows - 1; i >= 0; i--) {
      for (int k = rows; k < 2 * rows; k++) {
        A.m[i][k] /= A.m[i][i];
      }
      for (int rowMod = i - 1; rowMod >= 0; rowMod--) {
        for (int columMod = rows; columMod < 2 * rows; columMod++) {
          A.m[rowMod][columMod] -= A.m[i][columMod] * A.m[rowMod][i];
        }
      }
    }

    Matrix<T> Inv(rows, rows);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < rows; j++) {
        Inv.m[i][j] = A.m[i][j + rows];
      }
    }
    return Inv;
  }

  int Rows() const { return rows; }
  int Cols() const { return columns; }

  ~Matrix() {
    if (m != nullptr) {
      for (int i = 0; i < rows; i++) {
        delete[] m[i];
      }
      delete[] m;
    }
  }
};

#endif  // INCLUDE_MATRIX_HPP_
