// TO DO: change SquareMatrix m to an std::array of size 4.

#include "QuadraticAssignmentProblem/SVD/realsvd2x2.h"

void QAPSolver::internal::real_2x2_jacobi_svd(const SquareMatrix& matrix, int p, int q,
    JacobiRotation* j_left,
    JacobiRotation* j_right)
{
    SquareMatrix m(2);

    // IS THIS FLIPPED AROUND? p and q are flipped. this might mean apply on left and right are flipped?
    // also this should be changed to an array std::array<>(m) = {matrix(p, p), ... }
    m(0, 0) = matrix(q, q);
    m(1, 0) = matrix(p, q);
    m(0, 1) = matrix(q, p);
    m(1, 1) = matrix(p, p);

    JacobiRotation rot1;
    double t = m(0, 0) + m(1, 1);
    double d = m(1, 0) - m(0, 1);

    if (abs(d) < (std::numeric_limits<double>::min)())
    {
        rot1.s() = 0;
        rot1.c() = 1;
    }
    else
    {
        // If d!=0, then t/d cannot overflow because the magnitude of the
        // entries forming d are not too small compared to the ones forming t.
        double u = t / d;
        double tmp = sqrt(1 + u * u);
        rot1.s() = 1 / tmp;
        rot1.c() = u / tmp;
    }

    // below, p is assigned to 0 and q is assigned to 1
    m.applyOnTheLeft(0, 1, rot1);


    j_right->makeJacobi(m, 0, 1);
    *j_left = rot1 * j_right->transpose();
}