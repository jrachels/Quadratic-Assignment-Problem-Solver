//TO DO: When you compute the children of a node, you are computing all of the children at once. You should be able to use the intermediate values you compute from one child
// to make computing the other children faster.
// TO DO: to make it so that both Node and SparseFTOfMatrixFunc have the same YoungOrthogonalRepresentationConstants&, they should inherit from a Fourier Transform class that
// contains only a static YOR_Constants_ member
// TO DO: store ft_one_two_cycle_right_ as it's transpose so I don't have to compute the transpose in the svd algorithm. requires updating the svd algorithm
// TO DO: every child of the same parent gets  the same row bound for twotwocycle and onethreecycle components. I don't need to recompute it. 
// TO DO: Improve time complexity of computing bound of OneTwo cycle by (1) computing QR decomposition of sums of dyadics as you move down the tree and (2) compute product of sum of dyadics if you get too many entries
// ^ I don't think this is possible because dyadic matrices don't commute. reducing complexity requires that f^((n-1, 1)) matrices be low rank. Prove.


#ifndef QAP_NODE_H_
#define QAP_NODE_H_

#include <stdexcept>
#include <cmath>

#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"
#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"
#include "QuadraticAssignmentProblem/Matrices/matrixconcept.h"
#include "QuadraticAssignmentProblem/SVD/tracenormsumofdyadicmatrices.h"
#include "QuadraticAssignmentProblem/FourierTransform/youngorthogonalfouriertransform.h"

namespace QAPSolver {
	namespace Tests {
		class NodeTester;
	}
	namespace FourierTransform {

		namespace internal {
			// maybe these should be in an unnamed namespace rather than an internal namespace

			// Modify vector in a way equivalent to multiplying by the representation matrix of 
			// a left rotation of elements j through i>j.
			template<class FTOfQAPObjectiveFuncComponent>
			requires requires (FTOfQAPObjectiveFuncComponent& component, const int size, const int k) {
				component.ApplyAdjacentTranspositionTranspose(size, k);
			}
			void ApplyRotationTranspose(FTOfQAPObjectiveFuncComponent& component, const int size, const int j, const int i) { 
				assert(j <= i);

				for (int k = j; k < i; k++) {
					component.ApplyAdjacentTranspositionTranspose(size, k);
				}
			}
			double factorial(int n);
		}

		class Node : protected YoungOrthogonalFourierTransform {
			friend class SparseFTOfMatrixFunc;
			friend class QAPSolver::Tests::NodeTester;
		public:
			template<QAPSolver::MatrixConcept M>
			Node(const M& matrix1, const M& matrix2, YoungOrthogonalRepresentationConstants& yor_constants);

			Node(const Node& extended_ft, int index_of_rotation);

			Node() = default;

			double Bound() const;

			double Value() const;

		private:

			int size;

			//static YoungOrthogonalRepresentationConstants& yor_constants_;

			/*******************************************************/
			// BEGIN TWO TWO CYCLE *******************************/
			/*******************************************************/

			class FTAtTwoTwoCycleTableau {
				template<class FTOfQAPObjectiveFuncComponent>
				requires requires (FTOfQAPObjectiveFuncComponent& component, const int size, const int k) {
					component.ApplyAdjacentTranspositionTranspose(size, k);
				}
				friend void internal::ApplyRotationTranspose(FTOfQAPObjectiveFuncComponent& component, const int size, const int j, const int i);

			public:

				FTAtTwoTwoCycleTableau(const SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau& ft_matrix_two_two_cycle_1, const SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau& ft_matrix_two_two_cycle_2, double missing_factor);

				static FTAtTwoTwoCycleTableau FindChildComponentFromRotations(FTAtTwoTwoCycleTableau& ft_parent, int parent_size);

				FTAtTwoTwoCycleTableau() = default;

				double Bound_Contribution() const;

				FTAtTwoTwoCycleTableau Rotate(const int index_of_rotation, const int size) const;

				std::vector<double> ft_two_two_cycle_column_;
				std::vector<double> ft_two_two_cycle_row_;

			private:

				// Issue: this code is partially duplicated from the fftofmatrixfunc class
				void ApplyAdjacentTranspositionTranspose(int size, int i);
			};
			Node::FTAtTwoTwoCycleTableau ft_at_two_two_cycle_tableau_;

			/*******************************************************/
			// BEGIN ONE THREE CYCLE *******************************/
			/*******************************************************/

			class FTAtOneThreeCycleTableau {
				template<class FTOfQAPObjectiveFuncComponent>
				requires requires (FTOfQAPObjectiveFuncComponent& component, const int size, const int k) {
					component.ApplyAdjacentTranspositionTranspose(size, k);
				}
				friend void internal::ApplyRotationTranspose(FTOfQAPObjectiveFuncComponent& component, const int size, const int j, const int i);
			public:

				FTAtOneThreeCycleTableau(const SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau& ft_matrix_one_three_cycle_1, const SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau& ft_matrix_one_three_cycle_2, double missing_factor);

				static FTAtOneThreeCycleTableau FindChildComponentFromRotations(FTAtOneThreeCycleTableau& ft_parent, int parent_size);

				FTAtOneThreeCycleTableau() = default;

				double Bound_Contribution() const;

				FTAtOneThreeCycleTableau Rotate(const int index_of_rotation, const int size) const;

				std::vector<double> ft_one_three_cycle_column_;
				std::vector<double> ft_one_three_cycle_row_;

				// Issue: this code is partially duplicated from the fftofmatrixfunc class
				void ApplyAdjacentTranspositionTranspose(int size, int i);
			};
			Node::FTAtOneThreeCycleTableau ft_at_one_three_cycle_tableau_;

			/*******************************************************/
			// BEGIN ONE TWO CYCLE *******************************/
			/*******************************************************/

			// TO DO: store ft_one_two_cycle_right_ as it's transpose so I don't have to compute the transpose in the svd algorithm. requires updating the svd algorithm

			class FTAtOneTwoCycleTableau {
				template<class FTOfQAPObjectiveFuncComponent>
				requires requires (FTOfQAPObjectiveFuncComponent& component, const int size, const int k) {
					component.ApplyAdjacentTranspositionTranspose(size, k);
				}
				friend void internal::ApplyRotationTranspose(FTOfQAPObjectiveFuncComponent& component, const int size, const int j, const int i);
			public:

				FTAtOneTwoCycleTableau(const SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau& ft_matrix_one_two_cycle_1, const SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau& ft_matrix_one_two_cycle_2, double missing_factor);

				FTAtOneTwoCycleTableau() = default;


				double Bound_Contribution() const;

				static FTAtOneTwoCycleTableau FindChildComponentFromRotations(FTAtOneTwoCycleTableau& ft_one_two_cycle_parent, Node::FTAtTwoTwoCycleTableau& ft_two_two_cycle_parent,
					Node::FTAtOneThreeCycleTableau& ft_one_three_cycle_parent, int parent_size);

				// Issue: may produce too many copies. may produce 1 less copy if i reorder the basis elements in fftofmatrixfunc
				// should produce a smaller ft?
				FTAtOneTwoCycleTableau Rotate(const int index_of_rotation, const int size) const;

				SquareMatrix ft_one_two_cycle_;

			private:
				// Issue: this code is partially duplicated from the fftofmatrixfunc class
				void ApplyAdjacentTranspositionTranspose(int size, int k);
			};
			Node::FTAtOneTwoCycleTableau ft_at_one_two_cycle_tableau_;

			/*******************************************************/
			// BEGIN IDENTITY COMPONENT ****************************/
			/*******************************************************/

			class FTAtIdentityTableau {
				template<class FTOfQAPObjectiveFuncComponent>
				requires requires (FTOfQAPObjectiveFuncComponent& component, const int size, const int k) {
					component.ApplyAdjacentTranspositionTranspose(size, k);
				}
				friend void internal::ApplyRotationTranspose(FTOfQAPObjectiveFuncComponent& component, const int size, const int j, const int i);
			public:
				FTAtIdentityTableau(const SparseFTOfMatrixFunc::FFTAtIdentityTableau& ft_matrix_identity_tableau_1, const SparseFTOfMatrixFunc::FFTAtIdentityTableau& ft_matrix_identity_tableau_2, double missing_factor) {
					ft_identity_tableau_ = missing_factor * (ft_matrix_identity_tableau_1.ft_at_identity_tableau) * (ft_matrix_identity_tableau_2.ft_at_identity_tableau);
				}
				FTAtIdentityTableau() : ft_identity_tableau_(0) {};


				static FTAtIdentityTableau FindChildComponentFromRotations(FTAtIdentityTableau& ft_identity_tableau_parent, Node::FTAtOneTwoCycleTableau& ft_one_two_cycle_parent, int parent_size);

				double Bound_Contribution() const {
					return std::abs(ft_identity_tableau_);
				}

				// Issue: may produce too many copies. may produce 1 less copy if i reorder the basis elements in fftofmatrixfunc
				// should produce a smaller ft?
				FTAtIdentityTableau Rotate(const int index_of_rotation, const int size) const;

				double ft_identity_tableau_;

			private:
				// Issue: this code is partially duplicated from the fftofmatrixfunc class
				void ApplyAdjacentTranspositionTranspose(int size, int k) {
					// for each column
					return;
				}
			};
			Node::FTAtIdentityTableau ft_at_identity_tableau_;
		};


	}
}

#include "QuadraticAssignmentProblem/Node/src/node.tpp"

#endif