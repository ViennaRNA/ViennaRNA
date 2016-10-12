#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <assert.h>

template <class T>
class Matrix {
private:
    long m_;
    long n_;
    T *mtrx_;

public:
    Matrix(long m, long n) : m_(m), n_(n) {
        mtrx_=new T[m*n];
    }

    ~Matrix() {
        delete mtrx_;
    }

    inline const long xDim() const {
        return m_;
    }

    inline const long yDim() const {
        return n_;
    }

    inline const T& getAt(long x,long y) const {
        assert(x<m_ && y<n_);
        return mtrx_[y*m_+x];
    }

    inline void setAt(long x,long y, const T &val) {
        assert(x<m_ && y<n_);
        mtrx_[y*m_+x]=val;
    }
};

#endif
