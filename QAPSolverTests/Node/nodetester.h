#ifndef QAPSOLVERTESTS_NODE_NODETESTS_H
#define QAPSOLVERTESTS_NODE_NODETESTS_H

#include "QuadraticAssignmentProblem/Node/node.h"
#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"

namespace QAPSolver {
	namespace Tests {
		class NodeTester {
		public:
			explicit NodeTester(SquareMatrix& mat_left, SquareMatrix& mat_right) : node_(mat_left, mat_right, yor_constants_.GrowYORConstants(mat_left.size())) { /*yor_constants_.GrowYORConstants(mat_left.size()); */ }

			// This would throw an error if I haven't called the above constructor previously since Node would only have the default constructed yor_constants.
			// However, Node is an input, so I must have called Node previously with a larget node before calling this one.
			explicit NodeTester(NodeTester& entended_ft, int index_of_rotation) : node_(entended_ft.node_, index_of_rotation) {}

			std::vector<double> TwoTwoCycleComponentLeft() {
				return node_.ft_at_two_two_cycle_tableau_.ft_two_two_cycle_column_;
			}

			std::vector<double> TwoTwoCycleComponentRight() {
				return node_.ft_at_two_two_cycle_tableau_.ft_two_two_cycle_row_;
			}

			std::vector<double> OneThreeCycleComponentLeft() {
				return node_.ft_at_one_three_cycle_tableau_.ft_one_three_cycle_column_;
			}

			std::vector<double> OneThreeCycleComponentRight() {
				return node_.ft_at_one_three_cycle_tableau_.ft_one_three_cycle_row_;
			}

			SquareMatrix OneTwoCycleComponent() {
				return node_.ft_at_one_two_cycle_tableau_.ft_one_two_cycle_;
			}

			double IdentityComponent() {
				return node_.ft_at_identity_tableau_.ft_identity_tableau_;
			}

			double Bound() {
				return node_.Bound();
			}

		private:
			static YoungOrthogonalRepresentationConstants yor_constants_;
			QAPSolver::FourierTransform::Node node_;
		};
	}
}

#endif