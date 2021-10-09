// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

// TO DO: Handle instances of size = 1 and size = 2.
// TO DO: SVD algorithm sometimes, but rarely, loops infinitely

// An instance of QAPSolver is used to solve many QAP instances.
//
// Step 1: Initialize matrices of a matrix class that satisfies MatrixConcept.
// Note: SquareMatrix and Matrix, included in this library, satisfy the MatrixConcept.
//
//	QAPSolver::SquareMatrix mat1(3);
// 
// mat1(0, 0) = 0;
// mat1(0, 1) = 29;
// mat1(0, 2) = -13;
// mat1(1, 0) = 2;
// mat1(1, 1) = 0;
// mat1(1, 2) = 8;
// mat1(2, 0) = 17;
// mat1(2, 1) = -12;
// mat1(2, 2) = 0;
//
// Step 2: Initialize a QAPSolver instance.
//
// QAPSolver::QAPSolver qap_solver;
//
// Step 3: Pass an instance of the QAP to qap_solver.
// Note: perm is a Permutation defined in Permutation.cpp.
// total is a double.
//
// auto [perm, total] = qap_solver(mat1, mat2);
//
// Step 4: Convert perm to a vector using member function perm_to_vec.
//
// std::vector<int> vec = perm.perm_to_vec();
//
// See examples in qapsolvertests.cpp

#ifndef QAP_QAPSOLVER_QAPSOLVER_H_
#define QAP_QAPSOLVER_QAPSOLVER_H_

#include <tuple>
#include <memory>
#include <iostream>

#include "QuadraticAssignmentProblem/FourierTransform/src/youngorthogonalrepresentationconstants.h"
#include "QuadraticAssignmentProblem/Matrices/matrixconcept.h"
#include "QuadraticAssignmentProblem/Matrices/dimensionmismatch.h"
#include "QuadraticAssignmentProblem/Permutations/permutation.h"
#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"
#include "QuadraticAssignmentProblem/Node/node.h"

namespace QAPSolver {

	class QAPSolver
	{
	public:

		template<MatrixConcept M>
		[[nodiscard("QAP output not saved by class instance")]] std::tuple<Permutation, double> operator()(const M& matrix1, const M& matrix2) {
			// grow yor_constants
			int size = matrix1.size1();
			if (size != matrix1.size2()) {
				// fine accessing multiple times because this is is cold branch
				throw DimensionMismatch(matrix1.size1(), matrix1.size2(), matrix2.size1(), matrix2.size2());
			}
			if (!(size == matrix2.size1() && size == matrix2.size2())) {
				throw DimensionMismatch(matrix1.size1(), matrix1.size2(), matrix2.size1(), matrix2.size2());
			}

			YoungOrthogonalRepresentationConstants yor_constants(size + 1);

			//yor_constants_.GrowYORConstants(size);

			FourierTransform::Node FT_Objective_Func;
			try {
				//FT_Objective_Func = FourierTransform::Node{ FourierTransform::SparseFTOfMatrixFunc(matrix1, yor_constants), FourierTransform::SparseFTOfMatrixFunc(matrix2, yor_constants), yor_constants};
				FT_Objective_Func = FourierTransform::Node{ matrix1, matrix2, yor_constants };
			}
			catch (const std::out_of_range& e) {
				std::cout << e.what() << std::endl;
				std::cout << "Out of range error accessing matrix. size1() and size2() don't match operator() accessor domain." << std::endl;
				throw;
			}
			catch (...) {
				std::cout << "Unknown error constructing Root Node." << std::endl;
				throw;
			}

			// Structured Binding
			auto [optimal_permutation, optimal_objective_value] = BestLeaf(FT_Objective_Func, size, -INFINITY);

			// add in the diagonal values
			for (int i = 0; i < size; i++) {
				optimal_objective_value += matrix1(i, i) * matrix2(i, i);
			}

			return { Permutation(*optimal_permutation, size), optimal_objective_value };
		};



	private:
		class Node;

		std::tuple<std::unique_ptr<PermutationBuilder>, double> BestLeaf(const FourierTransform::Node& Restricted_FT_Objective_Func, const int size_left, double max_so_far);

	};

}

#endif