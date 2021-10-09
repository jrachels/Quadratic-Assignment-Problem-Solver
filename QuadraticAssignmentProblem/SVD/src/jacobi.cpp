
#include "QuadraticAssignmentProblem/SVD/jacobi.h"

bool QAPSolver::internal::JacobiRotation::makeJacobi(const double& x, const double& y, const double& z)
{
    double deno = 2 * std::abs(y);
    if (deno < (std::numeric_limits<double>::min)())
    {
        m_c_ = 1;
        m_s_ = 0;
        return false;
    }
    else
    {
        double tau = (x - z) / deno;
        double w = sqrt(tau * tau + 1);
        double t;
        if (tau > 0)
        {
            t = 1.0 / (tau + w);
        }
        else
        {
            t = 1.0 / (tau - w);
        }
        double sign_t = t > 0 ? 1 : -1;
        double n = 1 / std::sqrt(t * t + 1);
        m_s_ = -sign_t * (y / std::abs(y)) * std::abs(t) * n;
        m_c_ = n;
        return true;
    }
}