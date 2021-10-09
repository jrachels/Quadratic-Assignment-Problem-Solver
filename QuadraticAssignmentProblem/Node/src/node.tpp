#ifndef QAP_NODE_SRC_NODE_TPP_
#define QAP_NODE_SRC_NODE_TPP_

#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"
#include "QuadraticAssignmentProblem/Matrices/matrixconcept.h"
#include "QuadraticAssignmentProblem/FourierTransform/youngorthogonalfouriertransform.h"

template<QAPSolver::MatrixConcept M>
QAPSolver::FourierTransform::Node::Node(const M& matrix1, const M& matrix2, YoungOrthogonalRepresentationConstants& yor_constants) : YoungOrthogonalFourierTransform(yor_constants) {
	// error checking to make sure they are the right size
	//int size = matrix1.size1();
	size = matrix1.size1();

	SparseFTOfMatrixFunc ft_of_matrix_func_1{ matrix1, yor_constants };
	SparseFTOfMatrixFunc ft_of_matrix_func_2{ matrix2, yor_constants };

	// each entry of a FFTOfMatrixFunc is off by a factor of (n-2)!, so we will need to add that back in
	double missing_factor = internal::factorial(size - 2);

	ft_at_identity_tableau_ = FTAtIdentityTableau{ ft_of_matrix_func_1.ft_at_identity_tableau_, ft_of_matrix_func_2.ft_at_identity_tableau_, missing_factor };
	ft_at_one_two_cycle_tableau_ = FTAtOneTwoCycleTableau{ ft_of_matrix_func_1.ft_at_one_two_cycle_tableau_, ft_of_matrix_func_2.ft_at_one_two_cycle_tableau_, missing_factor };
	ft_at_two_two_cycle_tableau_ = FTAtTwoTwoCycleTableau{ ft_of_matrix_func_1.ft_at_two_two_cycle_tableau_, ft_of_matrix_func_2.ft_at_two_two_cycle_tableau_, missing_factor };
	ft_at_one_three_cycle_tableau_ = FTAtOneThreeCycleTableau{ ft_of_matrix_func_1.ft_at_one_three_cycle_tableau_, ft_of_matrix_func_2.ft_at_one_three_cycle_tableau_, missing_factor };
}

#endif