#include "QAPSolverTests/catch.hpp"

#include "QuadraticAssignmentProblem/Permutations/permutation.h"

TEST_CASE("Build a permuation and check the result") {
	QAPSolver::PermutationBuilder my_perm_builder{};

	my_perm_builder.AddRotation(1, 2);
	my_perm_builder.AddRotation(3, 3);
	my_perm_builder.AddRotation(4, 4);
	my_perm_builder.AddRotation(3, 5);
	my_perm_builder.AddRotation(2, 6);


	QAPSolver::Permutation my_perm{ my_perm_builder, 6 };

	REQUIRE(my_perm.perm_to_vec() == std::vector<int>{3, 1, 5, 6, 4, 2});

}