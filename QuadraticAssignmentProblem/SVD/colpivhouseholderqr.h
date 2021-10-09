#ifndef QAP_SVD_COLPIVHOUSEHOLDERQR_H
#define QAP_SVD_COLPIVHOUSEHOLDERQR_H

#include <vector>

#include "QuadraticAssignmentProblem/Matrices/matrix.h"

namespace QAPSolver {
	namespace internal {
		std::vector<int> ColPivHouseHolderQR(Matrix& m_qr);
	}
}

#endif