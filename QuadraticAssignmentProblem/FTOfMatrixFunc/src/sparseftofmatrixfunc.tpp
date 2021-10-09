#ifndef QAP_FTOFMATRIXFUNC_SRC_SPARSEFTOFMATRIXFUNC_TPP_
#define QAP_FTOFMATRIXFUNC_SRC_SPARSEFTOFMATRIXFUNC_TPP_

#include "QuadraticAssignmentProblem/Matrices/matrixconcept.h"
#include "QuadraticAssignmentProblem/FourierTransform/youngorthogonalfouriertransform.h"
#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"

template<QAPSolver::MatrixConcept M>
QAPSolver::FourierTransform::SparseFTOfMatrixFunc::SparseFTOfMatrixFunc(M& matrix, YoungOrthogonalRepresentationConstants& yor_constants) : YoungOrthogonalFourierTransform(yor_constants) {

	assert(matrix.size1() == matrix.size2()); // this should have already been checked
	//assert(matrix.size1() > 3);
	// The size of the input matrix is an important parameter throughout all of
	// these algorithms, so we pass it to the functions to use locally.
	const int size = matrix.size1();
	//yor_constants_.GrowYORConstants(size);
	FFTRowFuncs fft_row_funcs{ matrix, size };
	ft_at_identity_tableau_ = FFTAtIdentityTableau{ fft_row_funcs, size };
	ft_at_one_two_cycle_tableau_ = FFTAtOneTwoCycleTableau{ fft_row_funcs, size };
	ft_at_two_two_cycle_tableau_ = FFTAtTwoTwoCycleTableau{ fft_row_funcs, size };
	ft_at_one_three_cycle_tableau_ = FFTAtOneThreeCycleTableau{ fft_row_funcs, size };
}

template<QAPSolver::MatrixConcept M>
QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTRowFuncs::FFTRowFuncs(M& matrix, int size) {
	ft_row_funcs_.resize(size);
	ft_row_funcs_.shrink_to_fit();
	// Compute the relevant columns of the representation matrices of transpositions
	std::vector<std::vector<double>> representations_of_transpositions = RepresentationsOfTranspositions(size);//, yor_constants_);
	// Compute the column component of the FT using Clausen's formula.
	for (int row_iterator = 0; row_iterator < size; row_iterator++) {
		// initialize to 0
		std::get<0>(ft_row_funcs_[row_iterator]).resize(size - 2);
		std::get<0>(ft_row_funcs_[row_iterator]).shrink_to_fit();
		std::get<1>(ft_row_funcs_[row_iterator]) = 0;
		for (int column_iterator = 0; column_iterator < row_iterator; column_iterator++) {
			double matrix_entry = matrix(row_iterator, column_iterator);
			for (int k = 0; k < size - 2; k++) {
				std::get<0>(ft_row_funcs_[row_iterator])[k] += representations_of_transpositions[column_iterator][k] * matrix_entry;
			}
			std::get<1>(ft_row_funcs_[row_iterator]) += matrix_entry;
		}


		for (int column_iterator = row_iterator + 1; column_iterator < (size - 1); column_iterator++) {
			double matrix_entry = matrix(row_iterator, column_iterator);
			for (int k = 0; k < size - 2; k++) {
				std::get<0>(ft_row_funcs_[row_iterator])[k] += representations_of_transpositions[column_iterator][k] * matrix_entry;
			}
			std::get<1>(ft_row_funcs_[row_iterator]) += matrix_entry;
		}

		// move case outside loop to avoid checking this condition every time
		if (row_iterator != (size - 1))
		{
			double matrix_entry = matrix(row_iterator, size - 1);
			for (int k = 0; k < size - 2; k++) {
				std::get<0>(ft_row_funcs_[row_iterator])[k] += representations_of_transpositions[row_iterator][k] * matrix_entry;
			}
			std::get<1>(ft_row_funcs_[row_iterator]) += matrix_entry;
		}
	}
}

#endif