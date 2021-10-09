// TO DO: consider changing order  of basis vectors. (seems good 10/3)

#ifndef QAP_FTOFMATRIXFUNC_SPARSEFTOFMATRIXFUNC_H_
#define QAP_FTOFMATRIXFUNC_SPARSEFTOFMATRIXFUNC_H_

#include <vector>
#include <cassert>
#include <cmath>
#include <tuple>
#include <limits>

#include "QuadraticAssignmentProblem/Matrices/matrixconcept.h"
#include "QuadraticAssignmentProblem/FourierTransform/youngorthogonalfouriertransform.h"

namespace QAPSolver {
	namespace Tests {
		class SparseFTOfMatrixFuncTests;
	}
	namespace FourierTransform {
		namespace internal {
			// Modify vector in a way equivalent to multiplying by the representation matrix of 
			// a transposition of element i with element j<i.
			template<class FTOfMatrixComponent>
			requires requires (std::vector<double>& intermediate_vec, int size, int k) { FTOfMatrixComponent::ApplyAdjacentTransposition(intermediate_vec, size, k); }
			void ApplyTransposition(std::vector<double>& intermediate_vec, int size, int j, int i) {
				assert(j <= i);
				for (int k = i - 1; k > j - 1; k--) {
					FTOfMatrixComponent::ApplyAdjacentTransposition(intermediate_vec, size, k);
				}
				for (int k = j + 1; k < i; k++) {
					FTOfMatrixComponent::ApplyAdjacentTransposition(intermediate_vec, size, k);
				}
			}
		}

		class Node;


		// Efficiently computes the Fourier Transform of a matrix function f.
		// Given a square matrix mat of side length mat.size1(), f:S_n -> R^n
		// such that f(sigma)=mat[sigma(n), sigma(n-1)]. This Fourier Transform
		// can then be used to quickly compute FourierTransformQAPObjectiveFunc.
		// Time complexity is O(n^3), which is faster than Clausen's FFT
		// by a factor of n. Efficiency comes from exploiting sparsity. This FT is 
		// zero except at 4 partitions, (n), (n-1, 1), (n-2, 2), and (n-2, 1, 1), 
		// with the last 3 components being sparse. 

		// NOTE: Every entry is off by a factor of (n-2)!. This factor is  left out
		// for precision reasons and added back in when FourierTransformQAPObjectiveFunc
		// is constructed.

		// NOTE: Throughout this document I follow the convention that "transpositions"
		// are not necessarily "adjacent transpositions." 

		// While the FT of a matrix function is useful and interesting, for the 
		// purposes of this application, this is only a helper object. 
		// TO DO: Consider magnitude issues.
		// TO DO: Check working problem sizes (size = 1, 2, 3.)
		// TO DO: implement a pretty print member function(?)
		// TO DO: Throw error if can't access element of Matrix?
		// TO DO: Decide on a consistent naming scheme (FFT vs FT, adding sparse to names. redundancy. rule of thumb:
		// outermost class should have adjectives that describe all members. FFTOfMatrixFunc should be 
		// SparseFTOfMatrixFunc)
		// TO DO: Reorder internals of classes. put data members together
		// TO DO: Make constructor private?
		// TO DO: Put YorConstants in their own class. Make friend of FFTOfMatrixFunc and FTOfQAPOBjFun
		// TO DO: create constructor for SparseFTOfMatrixFunc that doesn't require a YoungOrthogonalRepresentationConstants

		class SparseFTOfMatrixFunc : protected YoungOrthogonalFourierTransform
		{

			friend class Node;
			//friend class SparseFTOfMatrixFuncTests;
			friend class QAPSolver::Tests::SparseFTOfMatrixFuncTests;
			//extern friend class QAPSolver::SparseFTOfMatrixFuncTests;
		public:
			template<QAPSolver::MatrixConcept M>
			SparseFTOfMatrixFunc(M& matrix, YoungOrthogonalRepresentationConstants& yor_constants);

		private:


			// Computes the only nonzero entries of the only nonzero component, (n-2, 1), of an
			// intermediate Fourier Transform "restricted to sigma_{n-1}" used in constructing
			// FFTAtOneTwoCycleTableau, FFTAtTwoTwoCycleTableau, and FFTAtOneThreeCycleTableau.
			// fft_row_funcs_[i] corresponds to fixing (n) to row i+1 (matrix starts at row 1) and allowing 
			// (n-1) to vary across the columns according to the permutation sigma_{n-1}.

			// In theory, ft_row_funcs_first_component[i] and ft_row_funcs_first_component[i] 
			// together form this i-th fourier transform.

			// Not stored as a component of the Fourier Transform of Matrix Functions.

			// TO DO: Combine the loops computing the two components into one.
			class FFTRowFuncs {
				friend class QAPSolver::Tests::SparseFTOfMatrixFuncTests;
			public:
				template<QAPSolver::MatrixConcept M>
				FFTRowFuncs(M& matrix, int size);

			private:

				// Compute the relevant column of the (n-1) representation matrices 
				// corresponding to the (n-1) transpositions.
				// TO DO: move multiplication by particular adjacency matrix to helper function.
				static std::vector<std::vector<double>> RepresentationsOfTranspositions(int size);

			public:
				std::vector<std::tuple<std::vector<double>, double>> ft_row_funcs_;
			};



			// This is this Fourier component corresponding to the (n) partition. Since representations
			// at this partition are the identity, the fourier transform is just the sum of the 
			// entries of the matrix times (n-2)! (with the (n-2)! factor excluded).
			class FFTAtIdentityTableau {
			public:
				FFTAtIdentityTableau(FFTRowFuncs& ft_row_funcs, int size) : ft_at_identity_tableau(0) {
					for (int i = 0; i < size; i++) {
						ft_at_identity_tableau += std::get<1>(ft_row_funcs.ft_row_funcs_[i]);
					}
				}
				FFTAtIdentityTableau() = default;
				double ft_at_identity_tableau;
			};
			SparseFTOfMatrixFunc::FFTAtIdentityTableau ft_at_identity_tableau_;


			class FFTAtOneTwoCycleTableau {
				template<class FTOfMatrixComponent>
				requires requires (std::vector<double>& intermediate_vec, int size, int k) { FTOfMatrixComponent::ApplyAdjacentTransposition(intermediate_vec, size, k); }
				friend void internal::ApplyTransposition(std::vector<double>& intermediate_vec, int size, int j, int i);
			public:
				FFTAtOneTwoCycleTableau(FFTRowFuncs& ft_row_funcs, int size);

				FFTAtOneTwoCycleTableau() = default;
				// First component is the (n-1) column, second component is the (n) column.
				std::tuple<std::vector<double>, std::vector<double>> sparse_ft_at_one_two_cycle_tableau_;
			private:
				// Modify vector in a way equivalent to multiplying by the representation matrix of 
				// an adjacent transposition of element i with element i+1.
				// TO DO: remove i==1 case from this function since it is only called at the end. handle outside of for loop.
				static void ApplyAdjacentTransposition(std::vector<double>& intermediate_vec, int size, int i);

			};
			SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau ft_at_one_two_cycle_tableau_;


			// Computes the component of Fourier Transform of a matrix function corresponding to the
			// Young Tableaux with Two Two cycles, or a row 2 of length 2. Basis vectors are identified
			// by row2:[i][j] corresponding to a standard tableaux for some i and j, with i and j possibly 
			// blank (indeterminant). Basis vectors are ordered by row2:[i_1][j_1] < row2:[i_2][j_2] if
			// j_1 > j_2 and row2:[i_1][j_1] < row2:[i_2][j_1] if i_1 < i_2. 
			class FFTAtTwoTwoCycleTableau {
				template<class FTOfMatrixComponent>
				requires requires (std::vector<double>& intermediate_vec, int size, int k) { FTOfMatrixComponent::ApplyAdjacentTransposition(intermediate_vec, size, k); }
				friend void internal::ApplyTransposition(std::vector<double>& intermediate_vec, int size, int j, int i);
			public:
				FFTAtTwoTwoCycleTableau(FFTRowFuncs& ft_row_funcs, int size);

				FFTAtTwoTwoCycleTableau() = default;
				std::vector<double> sparse_ft_at_two_two_cycle_tableau_;
			private:
				// Modify vector in a way equivalent to multiplying by the representation matrix of 
				// an adjacent transposition of element i with element i+1.
				// TO DO: remove i==1 case from this function since it is only called at the end. handle outside of for loop.
				static void ApplyAdjacentTransposition(std::vector<double>& intermediate_vec, int size, int i);
			};
			SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau ft_at_two_two_cycle_tableau_;


			// Computes the component of Fourier Transform of a matrix function corresponding to the
			// Young Tableaux with One Three cycle, or a row 2 of length 1 and a row 3 of length 1. 
			// Basis vectors are identified by row2:[i]row3:[j] corresponding to a standard tableaux 
			// for some i and j, with i and j possibly blank (indeterminant). Basis vectors are ordered 
			// by row2:[i_1]row3:[j_1] < row2:[i_2]row3:[j_2] if
			// j_1 > j_2 and row2:[i_1]row3:[j_1] < row2:[i_2]row3:[j_1] if i_1 < i_2. 
			class FFTAtOneThreeCycleTableau {
				template<class FTOfMatrixComponent>
				requires requires (std::vector<double>& intermediate_vec, int size, int k) { FTOfMatrixComponent::ApplyAdjacentTransposition(intermediate_vec, size, k); }
				friend void internal::ApplyTransposition(std::vector<double>& intermediate_vec, int size, int j, int i);
			public:
				FFTAtOneThreeCycleTableau(FFTRowFuncs& ft_row_funcs, int size);
				FFTAtOneThreeCycleTableau() = default;
				std::vector<double> sparse_ft_at_one_three_cycle_tableau_;
			private:

				// Modify vector in a way equivalent to multiplying by the representation matrix of 
				// an adjacent transposition of element i with element i+1.
				// TO DO: remove i==1 case from this function since it is only called at the end. handle outside of for loop.
				static void ApplyAdjacentTransposition(std::vector<double>& intermediate_vec, int size, int i);
			};
			SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau ft_at_one_three_cycle_tableau_;



		};
	}
}

#include "QuadraticAssignmentProblem/FTOfMatrixFunc/src/sparseftofmatrixfunc.tpp"

#endif