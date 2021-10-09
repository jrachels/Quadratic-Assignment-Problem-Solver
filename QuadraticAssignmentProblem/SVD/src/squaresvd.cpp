#include "QuadraticAssignmentProblem/SVD/squaresvd.h"
#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"

double QAPSolver::internal::SquareSVD::TraceNorm() {
	double total = 0;
	int num_singular_values = m_singularValues_.size();
	for (int i = 0; i < num_singular_values; i++) {
		total += m_singularValues_[i];
	}
	return total;
}

// SVD Algorithm based on Eigen's Jacobi SVD
QAPSolver::internal::SquareSVD::SquareSVD(QAPSolver::SquareMatrix m_workMatrix, int size) : m_singularValues_(size) {


	// using 
	// static constexpr double precision = 2 * std::numeric_limits<double>::epsilon();
	// causes infinite loop on some inputs
	static constexpr double precision = 2* std::numeric_limits<double>::epsilon();

	static constexpr double considerAsZero = (std::numeric_limits<double>::min)();

	// Scaling factor
	double scale = m_workMatrix.MaxAbsOfEntries();
	if (scale == 0) {
		scale = 1;
	}

	m_workMatrix.ScaleBy(1.0 / scale);

	double maxDiagEntry = m_workMatrix.MaxAbsOfDiagonal();


	bool finished = false;
	while (!finished)
	{
		finished = true;

		for (int p = 1; p < size; ++p)
		{
			for (int q = 0; q < p; ++q)
			{
				double threshold = std::max<double>(considerAsZero, precision * maxDiagEntry);
				if ((std::abs(m_workMatrix(p, q)) > threshold) || (std::abs(m_workMatrix(q, p)) > threshold))
				{
					finished = false;

					JacobiRotation j_left, j_right;
					real_2x2_jacobi_svd(m_workMatrix, p, q, &j_left, &j_right);

					m_workMatrix.applyOnTheLeft(p, q, j_left);

					m_workMatrix.applyOnTheRight(p, q, j_right);

					// keep track of the largest diagonal coefficient
					maxDiagEntry = std::max<double>(maxDiagEntry, std::max<double>(abs(m_workMatrix(p, p)), abs(m_workMatrix(q, q))));

				}
			}
		}
	}

	for (int i = 0; i < size; ++i)
	{
		m_singularValues_[i] = std::abs(m_workMatrix(i, i) * scale);

	}
}
