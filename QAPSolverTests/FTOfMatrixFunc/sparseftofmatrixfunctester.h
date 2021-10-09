
#ifndef QAPSOLVERTESTS_FTOFMATRIXFUNC_SPARSEFTOFMATRIXFUNCTESTS_H
#define QAPSOLVERTESTS_FTOFMATRIXFUNC_SPARSEFTOFMATRIXFUNCTESTS_H

#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"
#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"

//  REQUIRE_THAT( computed, Catch::Approx(known) );
namespace QAPSolver {
	namespace Tests {
		class SparseFTOfMatrixFuncTests {
		public:
			explicit SparseFTOfMatrixFuncTests(QAPSolver::SquareMatrix& mat) : mat_(mat), yor_constants_(mat.size()), sparse_ft_(mat, yor_constants_) {}

			std::vector<std::tuple<std::vector<double>, double>> FTOfRowFunctions() const {
				return QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTRowFuncs(mat_, mat_.size()).ft_row_funcs_;
			}

			std::vector<std::vector<double>> RepresentationsOfTranspositions(int size) {
				// GO BACK AND MAKE THIS FUNCTION PUBLIC
				return QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTRowFuncs::RepresentationsOfTranspositions(size);
			}

			double IdentityComponent() {
				return sparse_ft_.ft_at_identity_tableau_.ft_at_identity_tableau;
			}

			// for the column with row2:[n]
			std::vector<double> OneTwoCycleComponentN() {
				return std::get<1>(sparse_ft_.ft_at_one_two_cycle_tableau_.sparse_ft_at_one_two_cycle_tableau_);
			}

			// for the column with row2:[n-1]
			std::vector<double> OneTwoCycleComponentNminus() {
				return std::get<0>(sparse_ft_.ft_at_one_two_cycle_tableau_.sparse_ft_at_one_two_cycle_tableau_);
			}

			std::vector<double> TwoTwoCycleComponent() {
				return sparse_ft_.ft_at_two_two_cycle_tableau_.sparse_ft_at_two_two_cycle_tableau_;
			}

			std::vector<double> OneThreeCycleComponent() {
				return sparse_ft_.ft_at_one_three_cycle_tableau_.sparse_ft_at_one_three_cycle_tableau_;
			}

		private:
			QAPSolver::SquareMatrix& mat_;
			YoungOrthogonalRepresentationConstants yor_constants_;
			QAPSolver::FourierTransform::SparseFTOfMatrixFunc sparse_ft_;
		};
	}
}

#endif