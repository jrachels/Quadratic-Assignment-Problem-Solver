// TO DO: try a heap instead of sorting. need to be tested
// TO DO: consider that there might be an issue with storing a unique_ptr to an element of an array (node_ member of ChildNode)

#include <cassert>
#include <tuple>

#include "QuadraticAssignmentProblem/Include/qapsolver.h"
#include "QuadraticAssignmentProblem/Permutations/permutationbuilder.h"
#include "QuadraticAssignmentProblem/Node/node.h"


class ChildNode {
public:
	ChildNode(QAPSolver::FourierTransform::Node& node, double bound, int rotation_index) : node_(std::make_unique<QAPSolver::FourierTransform::Node>(node)), bound_(bound), rotation_index_(rotation_index) {}

	bool operator< (const ChildNode& node2)
	{
		return bound_ > node2.bound_;
	}

	bool operator== (const ChildNode& node2)
	{
		return bound_ == node2.bound_;
	}

	double Bound() {
		return bound_;
	}

	int RotationIndex() {
		return rotation_index_;
	}

	QAPSolver::FourierTransform::Node FT() {
		return *node_;
	}

private:

	std::unique_ptr<QAPSolver::FourierTransform::Node> node_;

	double bound_;

	int rotation_index_;
};


std::tuple<std::unique_ptr<QAPSolver::PermutationBuilder>, double> QAPSolver::QAPSolver::BestLeaf(const FourierTransform::Node& Restricted_FT_Objective_Func, const int size_left, double max_so_far) {
	// TO DO: figure out how to better reuse sigmas. 
	//PermutationBuilder sigma; // i need to make this on the heap?
	std::unique_ptr<PermutationBuilder> sigma{}; // i want this to start off as null, then become default if gucchi
	if (size_left == 1) {
		double value = Restricted_FT_Objective_Func.Value();
		if (value > max_so_far) {
			std::unique_ptr<PermutationBuilder> identity_sigma = std::make_unique<PermutationBuilder>();
			return std::make_tuple(std::move(identity_sigma), value);
		}
		else {
			return (std::make_tuple(std::move(sigma), max_so_far)); // return null
		}
	}
	else {
		std::vector<FourierTransform::Node> derivednodes;
		derivednodes.reserve(size_left);
		std::vector<ChildNode> derivednodereferences;
		derivednodereferences.reserve(size_left);
		for (int i = 0; i < size_left; i++) { //paper says 1 to k
			derivednodes.push_back(FourierTransform::Node(Restricted_FT_Objective_Func, i + 1));// already reserved
			derivednodereferences.push_back(ChildNode(derivednodes[i], derivednodes[i].Bound(), i + 1));//compute \hat{f}_i. put in std
		}
		// this might be better as a heap because I don't need to sort everything. I only need the nodes with the first few largest bounds.
		// need to test this
		std::sort(derivednodereferences.begin(), derivednodereferences.end()); // may not work since copy constructor has been deleted

		for (int i = 0; i < size_left; i++) {
			if (derivednodereferences[i].Bound() <= max_so_far) {
				break;
			}
			else {
				auto [tao, new_max] = BestLeaf(derivednodereferences[i].FT(), size_left - 1, max_so_far);
				if (tao) /* check if tao is null */ {
					sigma = std::move(tao);
					(*sigma).AddRotation(derivednodereferences[i].RotationIndex(), size_left);
					max_so_far = new_max;
				}
			}
		}
		return std::make_tuple(std::move(sigma), max_so_far);
	}
}
