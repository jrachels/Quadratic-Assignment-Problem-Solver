#include "QAPSolverTests/catch.hpp"
#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"
#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"
#include "QAPSolverTests/FTOfMatrixFunc/sparseftofmatrixfunctester.h"

TEST_CASE("Computing the Fourier Transform of a Matrix Func") {
	// initialize matrix

	// ALL ENTRIES ARE MISSING A FACTOR OF n!

	QAPSolver::SquareMatrix mat(4);
	mat(0, 0) = 0;
	mat(0, 1) = 1;
	mat(0, 2) = 2;
	mat(0, 3) = -3;
	mat(1, 0) = -6;
	mat(1, 1) = 0;
	mat(1, 2) = 2;
	mat(1, 3) = 4;
	mat(2, 0) = 1;
	mat(2, 1) = 2;
	mat(2, 2) = 0;
	mat(2, 3) = 9;
	mat(3, 0) = -2;
	mat(3, 1) = 3;
	mat(3, 2) = 3;
	mat(3, 3) = 0;

	QAPSolver::Tests::SparseFTOfMatrixFuncTests sparse_ft_of_matrix_func_tester(mat);

	//"Computing the FT of the row functions/Computing the restricted FT"
	// missing factor of (n-2)! = 2
	std::vector<std::tuple<std::vector<double>, double>> known_ft_of_row_functions;
	known_ft_of_row_functions.reserve(4);
	known_ft_of_row_functions.shrink_to_fit();
	known_ft_of_row_functions.push_back(std::make_tuple<std::vector<double>, double>({ 2 * std::sqrt(3), 3 }, 0));
	known_ft_of_row_functions.push_back(std::make_tuple<std::vector<double>, double>({ 5.0 * std::sqrt(3), 3 }, 0));
	known_ft_of_row_functions.push_back(std::make_tuple<std::vector<double>, double>({ std::sqrt(3) / 2.0, 7.5 }, 12));
	known_ft_of_row_functions.push_back(std::make_tuple<std::vector<double>, double>({ 5.0 * std::sqrt(3) / 2.0, 2.5 }, 4));
	std::vector<std::tuple<std::vector<double>, double>> computed_ft_of_row_functions = sparse_ft_of_matrix_func_tester.FTOfRowFunctions();
	for (int i = 0; i < std::max<int>(computed_ft_of_row_functions.size(), known_ft_of_row_functions.size()); ++i) {
		for (int j = 0; j < std::max<int>(std::get<0>(computed_ft_of_row_functions[i]).size(), std::get<0>(known_ft_of_row_functions[i]).size()); ++j) {
			REQUIRE(std::get<0>(computed_ft_of_row_functions[i])[j] == Approx(std::get<0>(known_ft_of_row_functions[i])[j]));
		}
		REQUIRE(std::get<1>(computed_ft_of_row_functions[i]) == Approx(std::get<1>(known_ft_of_row_functions[i])));
	}

	// Identity Component
	REQUIRE(sparse_ft_of_matrix_func_tester.IdentityComponent() == 16);

	// One Two Cycle Component

	// check first column with row2:[3]
	std::vector<double> known_one_two_cycle_component_Nminus = { 13.0 * std::sqrt(3) / 2.0, 11.5, 6 * std::sqrt(2) };
	std::vector<double> computed_one_two_cycle_component_Nminus = sparse_ft_of_matrix_func_tester.OneTwoCycleComponentNminus();
	for (int i = 0; i < std::max<int>(known_one_two_cycle_component_Nminus.size(), computed_one_two_cycle_component_Nminus.size()); ++i) {
		REQUIRE(known_one_two_cycle_component_Nminus[i] == Approx(computed_one_two_cycle_component_Nminus[i]));
	}

	// check second column with row2:[4]
	std::vector<double> known_one_two_cycle_component_N = { 0, 8 * sqrt(2), 0 };
	std::vector<double> computed_one_two_cycle_component_N = sparse_ft_of_matrix_func_tester.OneTwoCycleComponentN();
	for (int i = 0; i < std::max<int>(known_one_two_cycle_component_N.size(), computed_one_two_cycle_component_N.size()); ++i) {
		REQUIRE(known_one_two_cycle_component_N[i] == Approx(computed_one_two_cycle_component_N[i]));
	}

	// Two Two Cycle Component

	// check column with row2:[3][4]
	std::vector<double> known_two_two_cycle_component = { 11.0 * sqrt(3) / 2.0, 2.5 };
	std::vector<double> computed_two_two_cycle_component = sparse_ft_of_matrix_func_tester.TwoTwoCycleComponent();
	for (int i = 0; i < std::max<int>(known_two_two_cycle_component.size(), computed_two_two_cycle_component.size()); ++i) {
		REQUIRE(known_two_two_cycle_component[i] == Approx(computed_two_two_cycle_component[i]));
	}

	// One Three Cycle Component

	// check column with row2:[3]row3:[4]
	std::vector<double> known_one_three_cycle_component = { -7 * std::sqrt(3) / 2, -6.5, -2 * std::sqrt(6) };
	std::vector<double> computed_one_three_cycle_component = sparse_ft_of_matrix_func_tester.OneThreeCycleComponent();
	for (int i = 0; i < std::max<int>(known_one_three_cycle_component.size(), computed_one_three_cycle_component.size()); ++i) {
		REQUIRE(known_one_three_cycle_component[i] == Approx(computed_one_three_cycle_component[i]));
	}

}