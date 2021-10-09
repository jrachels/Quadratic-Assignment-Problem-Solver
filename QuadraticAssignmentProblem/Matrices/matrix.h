#ifndef QAP_MATRICES_MATRIX_H_
#define QAP_MATRICES_MATRIX_H_

#include <vector>

namespace QAPSolver {
	class Matrix
	{
	public:
		explicit Matrix(int m, int n) : height_(m), width_(n), mat_(m* n) { //returns matrix of all 0s
			mat_.shrink_to_fit();
		}

		Matrix() = default;

		double& operator() (int i, int j) {
			// stored as [column1, column2, column3,..., columnn] for cache reasons.
			// remember to copy column first!
			return mat_[i + j * height_];
		}

		double operator() (int i, int j) const {
			// stored as [column1, column2, column3,..., columnn] for cache reasons.
			// remember to copy column first!
			return mat_[i + j * height_];
		}

		void ScaleBy(double scale);

		int size1() const {
			return height_;
		}

		int size2() const {
			return width_;
		}

		double MaxAbsOfEntries() const;

		Matrix Transpose() const;

	private:

		// forced const
		int height_;
		int width_;
		// changeable
		std::vector<double> mat_;
	};

	namespace internal {
		double NormOfColumn(const Matrix& matrix, int column_index);
		double NormOfTailOfColumn(const Matrix& matrix, int column_index, int starting_row_index);
		void makeHouseHolderInColumn(Matrix& matrix, int column_index, double& tau, double& beta);
		void applyHouseHolderOnLeft(Matrix& matrix, int column_index_of_householder, const double tau/*, double* workspace*/);

	}
}
#endif
