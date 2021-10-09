#include "QuadraticAssignmentProblem/Matrices/matrix.h"

#include <cmath>

void QAPSolver::Matrix::ScaleBy(double scale) {
	int size = mat_.size();
	for (int i = 0; i < size; i++) {
		double temp = mat_[i];
		mat_[i] = temp * scale;
	}
	return;
}

double QAPSolver::Matrix::MaxAbsOfEntries() const {
	double max = 0;
	for (double entry : mat_) {
		if (std::abs(entry) > max) {
			max = entry;
		}
	}
	return max;
}

QAPSolver::Matrix QAPSolver::Matrix::Transpose() const {
	Matrix transpose(width_, height_);
	// should test if I should width or height first for speed
	for (int i = 0; i < width_; ++i) {
		for (int j = 0; j < height_; ++j) {
			transpose(i, j) = (this)->operator()(j, i);
		}
	}
	return transpose;
}

double QAPSolver::internal::NormOfColumn(const Matrix& matrix, int column_index) {
	return NormOfTailOfColumn(matrix, column_index, 0);
}

double QAPSolver::internal::NormOfTailOfColumn(const Matrix& matrix, int column_index, int starting_row_index) {
	double squared_norm = 0;
	int height = matrix.size1();
	for (int i = starting_row_index; i < height; i++) {
		double entry = matrix(i, column_index);
		squared_norm += entry * entry;
	}
	return std::sqrt(squared_norm);
}

void QAPSolver::internal::makeHouseHolderInColumn(Matrix& matrix, int column_index, double& tau, double& beta) {
	int num_rows = matrix.size1();

	// compute squarednorm of the tail
	double tailSqNorm = 0;
	for (int i = column_index + 1; i < num_rows; i++) {
		tailSqNorm += std::pow(matrix(i, column_index), 2);
	}

	double c0 = matrix(column_index, column_index);

	constexpr double tol = (std::numeric_limits<double>::min)();

	if (tailSqNorm <= tol) {
		tau = 0;
		beta = c0;
		// set the rest to 0
		for (int i = column_index + 1; i < num_rows; i++) {
			matrix(i, column_index) = 0;
		}
	}
	else {
		beta = std::sqrt((c0 * c0) + tailSqNorm);
		if (c0 >= 0) {
			beta = -beta;
		}
		// set the tail
		{
			const double temp_const = c0 - beta;
			for (int i = column_index + 1; i < num_rows; i++) {
				double temp_entry = matrix(i, column_index);
				matrix(i, column_index) = temp_entry / (c0 - beta);
			}
		}
		tau = (beta - c0) / beta;
	}
}

void QAPSolver::internal::applyHouseHolderOnLeft(Matrix& matrix, int column_index_of_householder, const double tau/*, double* workspace*/) {

	int rows = matrix.size1();
	int cols = matrix.size2();

	int columns_remaining = cols - column_index_of_householder - 1;
	//int rows_below = rows - k - 1;


	// Draw a picture to understand. See Eigen documentation for Householder.h
	std::vector<double> tmp(columns_remaining);
	for (int column_shift = 0; column_shift < columns_remaining; column_shift++) {
		for (int row_of_matrix = column_index_of_householder + 1; row_of_matrix < rows; row_of_matrix++) {
			tmp[column_shift] += matrix(row_of_matrix, column_index_of_householder) * matrix(row_of_matrix, column_index_of_householder + 1 + column_shift);
		}
		tmp[column_shift] += matrix(column_index_of_householder, column_index_of_householder + 1 + column_shift);
	}

	for (int column_shift = 0; column_shift < columns_remaining; column_shift++) {
		matrix(column_index_of_householder, column_index_of_householder + 1 + column_shift) -= tau * tmp[column_shift];
	}

	// now figure out bottom
	for (int column_shift = 0; column_shift < columns_remaining; column_shift++) {
		for (int row_of_matrix = column_index_of_householder + 1; row_of_matrix < rows; row_of_matrix++) {
			matrix(row_of_matrix, column_index_of_householder + 1 + column_shift) -= tau * matrix(row_of_matrix, column_index_of_householder) * tmp[column_shift];
		}
	}
}

