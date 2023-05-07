// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <fstream>
#include <vector>

template <typename T> inline void writeScalar(std::ofstream &out,T &scalar) {
  if (out.is_open()) {
    out.write(reinterpret_cast<char *>(&scalar), sizeof(T));
  } else {
    throw std::runtime_error("Output file for scalar is not open!");
  }
}

template <typename T>
inline void writeScalar(std::ofstream &out, T &&scalar) {
  if (out.is_open()) {
    out.write(reinterpret_cast<char *>(&scalar), sizeof(T));
  } else {
    throw std::runtime_error("Output file for scalar is not open!");
  }
}

template <typename T> inline void readScalar(std::ifstream &out, T &scalar) {
  if (out.is_open()) {
    out.read(reinterpret_cast<char *>(&scalar), sizeof(T));
  } else {
    throw std::runtime_error("Output file for scalar is not open!");
  }
}



template <class EigenMatrix>
inline void writeEigenMatrix(std::ofstream &out, EigenMatrix &Matrix) {
  if (out.is_open()) {
    typename EigenMatrix::Index rows = Matrix.rows(), cols = Matrix.cols();
    writeScalar(out, rows);
    writeScalar(out, cols);
    out.write(reinterpret_cast<char *>(Matrix.data()),
              rows * cols * sizeof(typename EigenMatrix::Scalar));
  } else {
    throw std::runtime_error("Output file for Eigen matrix is not open!");
  }
}

template <class EigenMatrix>
inline void readEigenMatrix(std::ifstream &in, EigenMatrix &Matrix) {
  if (in.is_open()) {
    typename EigenMatrix::Index rows, cols;
    readScalar(in, rows);
    readScalar(in, cols);
    Matrix.resize(rows, cols);
    in.read(reinterpret_cast<char *>(Matrix.data()),
              rows * cols * sizeof(typename EigenMatrix::Scalar));
  } else {
    throw std::runtime_error("Output file for Eigen matrix is not open!");
  }
}

template<class EigenSparseMatrix>
inline void writeEigenSparseMatrix(std::ofstream &out, EigenSparseMatrix &Matrix)
{
  if (out.is_open()) {
    typedef typename EigenSparseMatrix::StorageIndex IND;
    typedef typename EigenSparseMatrix::Scalar SCA;
    IND rows, cols, nnzs, outerSize, innerSize;
    rows = Matrix.rows();
    cols = Matrix.cols();
    nnzs = Matrix.nonZeros();
    outerSize = Matrix.outerSize();
    innerSize = Matrix.innerSize();

    writeScalar(out, rows);
    writeScalar(out, cols);
    writeScalar(out, nnzs);
    writeScalar(out, outerSize);
    writeScalar(out, innerSize);
    
    out.write(reinterpret_cast<char *>(Matrix.outerIndexPtr()),
              outerSize * sizeof(IND));
    out.write(reinterpret_cast<char *>(Matrix.innerIndexPtr()),
              nnzs * sizeof(IND));
    out.write(reinterpret_cast<char *>(Matrix.valuePtr()), sizeof(SCA) * nnzs);
  } else {
    throw std::runtime_error("Output file for Eigen sparse matrix is not open!");
  }
}

template <class EigenSparseMatrix>
inline void readEigenSparseMatrix(std::ifstream &in,
                                   EigenSparseMatrix &Matrix) {
  if (in.is_open()) {
    typedef typename EigenSparseMatrix::StorageIndex IND;
    typedef typename EigenSparseMatrix::Scalar SCA;

    IND rows, cols, nnzs, outerSize, innerSize;

    
    readScalar(in, rows);
    readScalar(in, cols);
    readScalar(in, nnzs);
    readScalar(in, outerSize);
    readScalar(in, innerSize);

    Matrix.resize(rows, cols);
    Matrix.makeCompressed();
    Matrix.resizeNonZeros(nnzs);

    in.read(reinterpret_cast<char *>(Matrix.outerIndexPtr()),
            outerSize * sizeof(IND));
    in.read(reinterpret_cast<char *>(Matrix.innerIndexPtr()),
            nnzs * sizeof(IND));
    in.read(reinterpret_cast<char *>(Matrix.valuePtr()),
            sizeof(SCA) * nnzs);

    Matrix.finalize();
  } else {
    throw std::runtime_error(
        "Output file for Eigen sparse matrix is not open!");
  }
}

template <class EigenSparseVector>
inline void writeEigenSparseVector(std::ofstream &out,
                                   EigenSparseVector &Vector) {
  if (out.is_open()) {
    typedef typename EigenSparseVector::StorageIndex IND;
    //typedef typename EigenSparseVector::Scalar SCA;
    IND rows, nnzs, innerSize;
    rows = Vector.rows();
    nnzs = Vector.nonZeros();
    innerSize = Vector.innerSize();

    writeScalar(out, rows);
    writeScalar(out, nnzs);
    writeScalar(out, innerSize);
    
    out.write(reinterpret_cast<char *>(Vector.innerIndexPtr()),
              nnzs * sizeof(IND));
    
  } else {
    throw std::runtime_error(
        "Output file for Eigen sparse matrix is not open!");
  }
}

template <class EigenSparseVector>
inline void readEigenSparseVector(std::ifstream &in,
                                  EigenSparseVector &Vector) {
  if (in.is_open()) {
    typedef typename EigenSparseVector::StorageIndex IND;
    //typedef typename EigenSparseVector::Scalar SCA;

    IND rows, nnzs, innerSize;

    readScalar(in, rows);
    readScalar(in, nnzs);
    readScalar(in, innerSize);

    Vector.resize(rows);
    Vector.resizeNonZeros(nnzs);
    
    in.read(reinterpret_cast<char *>(Vector.innerIndexPtr()),
            nnzs * sizeof(IND));

    Vector.finalize();
  } else {
    throw std::runtime_error(
        "Output file for Eigen sparse matrix is not open!");
  }
}


template <typename T>
inline void writeStdVector(std::ofstream &out, std::vector<T> &vec) {
  if (out.is_open()) {
    std::size_t vsize = vec.size();
    out.write(reinterpret_cast<char *>(&vsize), sizeof(std::size_t));
    out.write(reinterpret_cast<char *>(&vec[0]), vsize * sizeof(T));
  } else {
    throw std::runtime_error("Output file for vector is not open!");
  }
}

template <typename T>
inline void readStdVector(std::ifstream &out, std::vector<T> &vec) {
  if (out.is_open()) {
    std::size_t vsize;
    out.read(reinterpret_cast<char *>(&vsize), sizeof(std::size_t));
    vec.resize(vsize);
    out.read(reinterpret_cast<char *>(&vec[0]), vsize * sizeof(T));
  } else {
    throw std::runtime_error("Output file for vector is not open!");
  }
}