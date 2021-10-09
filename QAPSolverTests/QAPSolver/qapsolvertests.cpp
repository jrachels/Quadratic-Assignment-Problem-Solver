#include "QAPSolverTests/catch.hpp"
#include "QAPSolverTests/Node/nodetester.h"

#include "QuadraticAssignmentProblem/Include/qapsolver.h"

TEST_CASE("Instance of the Quadratic Assignment Problem of size 3") {

	QAPSolver::SquareMatrix mat1(3);

	mat1(0, 0) = 0;
	mat1(0, 1) = 29;
	mat1(0, 2) = -13;
	mat1(1, 0) = 2;
	mat1(1, 1) = 0;
	mat1(1, 2) = 8;
	mat1(2, 0) = 17;
	mat1(2, 1) = -12;
	mat1(2, 2) = 0;


	QAPSolver::SquareMatrix mat2(3);

	mat2(0, 0) = 0;
	mat2(0, 1) = 5;
	mat2(0, 2) = -2;
	mat2(1, 0) = -28;
	mat2(1, 1) = 0;
	mat2(1, 2) = -8;
	mat2(2, 0) = 11;
	mat2(2, 1) = 31;
	mat2(2, 2) = 0;

	QAPSolver::QAPSolver qap_solver;
	auto [perm, total] = qap_solver(mat1, mat2);
	std::vector<int> vec = perm.perm_to_vec();
	std::vector<int> correct_perm = std::vector<int>({ 3, 2, 1 });
	int max_perm_length = std::max(vec.size(), correct_perm.size());
	for (int i = 0; i < max_perm_length; ++i) {
		REQUIRE(vec[i] == correct_perm[i]);
	}
	REQUIRE(total == 422);

}

TEST_CASE("Instance of the Quadratic Assignment Problem of size 4") {

	QAPSolver::SquareMatrix mat1(4);

	mat1(0, 0) = 0;
	mat1(0, 1) = 29;
	mat1(0, 2) = -13;
	mat1(0, 3) = -14;
	mat1(1, 0) = 2;
	mat1(1, 1) = 0;
	mat1(1, 2) = 8;
	mat1(1, 3) = 4;
	mat1(2, 0) = 17;
	mat1(2, 1) = -12;
	mat1(2, 2) = 0;
	mat1(2, 3) = 10;
	mat1(3, 0) = -11;
	mat1(3, 1) = -18;
	mat1(3, 2) = -19;
	mat1(3, 3) = 0;

	QAPSolver::SquareMatrix mat2(4);

	mat2(0, 0) = 0;
	mat2(0, 1) = 5;
	mat2(0, 2) = -2;
	mat2(0, 3) = 8;
	mat2(1, 0) = -28;
	mat2(1, 1) = 0;
	mat2(1, 2) = -8;
	mat2(1, 3) = -1;
	mat2(2, 0) = 11;
	mat2(2, 1) = 31;
	mat2(2, 2) = 0;
	mat2(2, 3) = 4;
	mat2(3, 0) = -3;
	mat2(3, 1) = -22;
	mat2(3, 2) = 18;
	mat2(3, 3) = 0;

	QAPSolver::QAPSolver qap_solver;
	auto [perm, total] = qap_solver(mat1, mat2);
	std::vector<int> vec = perm.perm_to_vec();
	std::vector<int> correct_perm = std::vector<int>({ 3, 4, 2, 1 });
	int max_perm_length = std::max(vec.size(), correct_perm.size());
	for (int i = 0; i < max_perm_length; ++i) {
		REQUIRE(vec[i] == correct_perm[i]);
	}
	REQUIRE(total == 1986);

}


TEST_CASE("Instance of the Quadratic Assignment Problem of size 7") {
	QAPSolver::SquareMatrix left_mat(7);

	left_mat(0, 0) = 0;
	left_mat(0, 1) = -24;
	left_mat(0, 2) = 32;
	left_mat(0, 3) = -26;
	left_mat(0, 4) = 32;
	left_mat(0, 5) = 6;
	left_mat(0, 6) = -2;
	left_mat(1, 0) = -14;
	left_mat(1, 1) = 0;
	left_mat(1, 2) = -15;
	left_mat(1, 3) = -6;
	left_mat(1, 4) = 31;
	left_mat(1, 5) = 16;
	left_mat(1, 6) = 30;
	left_mat(2, 0) = -24;
	left_mat(2, 1) = -12;
	left_mat(2, 2) = 0;
	left_mat(2, 3) = 18;
	left_mat(2, 4) = 32;
	left_mat(2, 5) = -22;
	left_mat(2, 6) = -22;
	left_mat(3, 0) = -9;
	left_mat(3, 1) = -3;
	left_mat(3, 2) = -23;
	left_mat(3, 3) = 0;
	left_mat(3, 4) = 22;
	left_mat(3, 5) = -23;
	left_mat(3, 6) = 19;
	left_mat(4, 0) = -15;
	left_mat(4, 1) = 29;
	left_mat(4, 2) = 30;
	left_mat(4, 3) = -12;
	left_mat(4, 4) = 0;
	left_mat(4, 5) = 3;
	left_mat(4, 6) = 17;
	left_mat(5, 0) = 11;
	left_mat(5, 1) = -24;
	left_mat(5, 2) = -9;
	left_mat(5, 3) = -32;
	left_mat(5, 4) = 27;
	left_mat(5, 5) = 0;
	left_mat(5, 6) = -6;
	left_mat(6, 0) = -11;
	left_mat(6, 1) = 32;
	left_mat(6, 2) = 5;
	left_mat(6, 3) = 9;
	left_mat(6, 4) = 6;
	left_mat(6, 5) = 5;
	left_mat(6, 6) = 0;

	QAPSolver::SquareMatrix right_mat(7);

	right_mat(0, 0) = 0;
	right_mat(0, 1) = -14;
	right_mat(0, 2) = 17;
	right_mat(0, 3) = 31;
	right_mat(0, 4) = 24;
	right_mat(0, 5) = -1;
	right_mat(0, 6) = 2;
	right_mat(1, 0) = 20;
	right_mat(1, 1) = 0;
	right_mat(1, 2) = -16;
	right_mat(1, 3) = 16;
	right_mat(1, 4) = 18;
	right_mat(1, 5) = 26;
	right_mat(1, 6) = -26;
	right_mat(2, 0) = -14;
	right_mat(2, 1) = 9;
	right_mat(2, 2) = 0;
	right_mat(2, 3) = -21;
	right_mat(2, 4) = 31;
	right_mat(2, 5) = 31;
	right_mat(2, 6) = -19;
	right_mat(3, 0) = -19;
	right_mat(3, 1) = -26;
	right_mat(3, 2) = -26;
	right_mat(3, 3) = 0;
	right_mat(3, 4) = 0;
	right_mat(3, 5) = 5;
	right_mat(3, 6) = 31;
	right_mat(4, 0) = -9;
	right_mat(4, 1) = 26;
	right_mat(4, 2) = 29;
	right_mat(4, 3) = -30;
	right_mat(4, 4) = 0;
	right_mat(4, 5) = 11;
	right_mat(4, 6) = 29;
	right_mat(5, 0) = -5;
	right_mat(5, 1) = -3;
	right_mat(5, 2) = -29;
	right_mat(5, 3) = -25;
	right_mat(5, 4) = -19;
	right_mat(5, 5) = 0;
	right_mat(5, 6) = 1;
	right_mat(6, 0) = 27;
	right_mat(6, 1) = -5;
	right_mat(6, 2) = -31;
	right_mat(6, 3) = 25;
	right_mat(6, 4) = -27;
	right_mat(6, 5) = 0;
	right_mat(6, 6) = 0;

	QAPSolver::QAPSolver qap_solver;
	auto [perm, total] = qap_solver(left_mat, right_mat);
	std::vector<int> vec = perm.perm_to_vec();
	std::vector<int> correct_perm = std::vector<int>({ 2, 3, 6, 4, 1, 5, 7 });
	int max_perm_length = std::max(vec.size(), correct_perm.size());
	for (int i = 0; i < max_perm_length; ++i) {
		REQUIRE(vec[i] == correct_perm[i]);
	}
	REQUIRE(total == Approx(9314));

}