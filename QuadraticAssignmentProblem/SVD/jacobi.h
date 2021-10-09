#ifndef QAP_SVD_JACOBI_H
#define QAP_SVD_JACOBI_H

#include <cmath>

#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"

namespace QAPSolver {
    namespace internal {
        class JacobiRotation
        {
        public:

            JacobiRotation() : m_c_(0), m_s_(0) {}

            JacobiRotation(const double c, const double s) : m_c_(c), m_s_(s) {}

            double& c() { return m_c_; }
            double c() const { return m_c_; }
            double& s() { return m_s_; }
            double s() const { return m_s_; }

            JacobiRotation operator*(const JacobiRotation& other)
            {
                // removed conj
                return JacobiRotation(m_c_ * other.m_c_ - m_s_ * other.m_s_,
                    m_c_ * other.m_s_ + m_s_ * other.m_c_);
            }

            JacobiRotation transpose() const { return JacobiRotation(m_c_, -m_s_); }

            JacobiRotation adjoint() const { return JacobiRotation(m_c_, -m_s_); }

            inline bool makeJacobi(const SquareMatrix&, int p, int q);
            bool makeJacobi(const double& x, const double& y, const double& z);

        protected:

            double m_c_, m_s_;
        };

        inline bool JacobiRotation::makeJacobi(const SquareMatrix& m, int p, int q)
        {
            // why do u need only 3 things here?
            return makeJacobi(m(p, p), m(p, q), m(q, q));
        }

    }
    /****************************************************************************************
    *   Implementation of MatrixBase methods
    ****************************************************************************************/

    // j goes on left
    inline void SquareMatrix::applyOnTheLeft(int p, int q, const internal::JacobiRotation& j)
    {
        //for each of the two rows
        int size = this->size();
        for (int i = 0; i < size; i++) {
            double row_q_entry = (*this)(q, i);
            double row_p_entry = (*this)(p, i);
            double jc = j.c();
            double js = j.s();
            (*this)(q, i) = jc * row_q_entry + js * row_p_entry;
            (*this)(p, i) = -js * row_q_entry + jc * row_p_entry;
        }
    }

    // j goes on right
    inline void SquareMatrix::applyOnTheRight(int p, int q, const internal::JacobiRotation& j)
    {
        //for each of the two rows
        int size = this->size();
        for (int i = 0; i < size; i++) {
            double column_q_entry = (*this)(i, q);
            double column_p_entry = (*this)(i, p);
            double jc = j.c();
            double js = j.s();
            (*this)(i, q) = jc * column_q_entry - js * column_p_entry;
            (*this)(i, p) = js * column_q_entry + jc * column_p_entry;
        }
    }

}

#endif // EIGEN_JACOBI_H