#ifndef QAP_SVD_REALSVD2X2_H
#define QAP_SVD_REALSVD2X2_H

#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"
#include "QuadraticAssignmentProblem/SVD/jacobi.h"

namespace QAPSolver {
    namespace internal {
        void real_2x2_jacobi_svd(const SquareMatrix& matrix, int p, int q,
            JacobiRotation* j_left,
            JacobiRotation* j_right);
    }
}

#endif