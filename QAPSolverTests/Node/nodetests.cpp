// TO DO: Test rotations before testing children

#include "QAPSolverTests/catch.hpp"
#include "QAPSolverTests/Node/nodetester.h"
#include "QAPSolverTests/FTOfMatrixFunc/sparseftofmatrixfunctester.h"

QAPSolver::YoungOrthogonalRepresentationConstants QAPSolver::Tests::NodeTester::yor_constants_(0);


TEST_CASE("Computing the Fourier Transform of a QAP objective function of size 4 and its children") {

	QAPSolver::SquareMatrix left_mat(4);

	left_mat(0, 0) = 0;
	left_mat(0, 1) = 29;
	left_mat(0, 2) = -13;
	left_mat(0, 3) = -14;
	left_mat(1, 0) = 2;
	left_mat(1, 1) = 0;
	left_mat(1, 2) = 8;
	left_mat(1, 3) = 4;
	left_mat(2, 0) = 17;
	left_mat(2, 1) = -12;
	left_mat(2, 2) = 0;
	left_mat(2, 3) = 10;
	left_mat(3, 0) = -11;
	left_mat(3, 1) = -18;
	left_mat(3, 2) = -19;
	left_mat(3, 3) = 0;

	QAPSolver::SquareMatrix right_mat(4);

	right_mat(0, 0) = 0;
	right_mat(0, 1) = 5;
	right_mat(0, 2) = -2;
	right_mat(0, 3) = 8;
	right_mat(1, 0) = -28;
	right_mat(1, 1) = 0;
	right_mat(1, 2) = -8;
	right_mat(1, 3) = -1;
	right_mat(2, 0) = 11;
	right_mat(2, 1) = 31;
	right_mat(2, 2) = 0;
	right_mat(2, 3) = 4;
	right_mat(3, 0) = -3;
	right_mat(3, 1) = -22;
	right_mat(3, 2) = 18;
	right_mat(3, 3) = 0;



	//QAPSolver::Tests::SparseFTOfMatrixFuncTests left_tester(left_mat);
	//QAPSolver::Tests::SparseFTOfMatrixFuncTests right_tester(right_mat);

	QAPSolver::Tests::NodeTester node_of_objective_func(left_mat, right_mat);

	// Check Construction of node from matrices

	// Identity Component
	REQUIRE(node_of_objective_func.IdentityComponent() == -442);
	//  REQUIRE_THAT( computed, Catch::Approx(known) );
	REQUIRE_THAT(node_of_objective_func.TwoTwoCycleComponentLeft(), Catch::Approx(std::vector<double>({19 * sqrt(3), 83})));

	REQUIRE_THAT(node_of_objective_func.TwoTwoCycleComponentRight(), Catch::Approx(std::vector<double>({-21 * sqrt(3), -8})));

	REQUIRE_THAT(node_of_objective_func.OneThreeCycleComponentLeft(), Catch::Approx(std::vector<double>({-79 * sqrt(3)/3, -29, -154*sqrt(6)/3})));

	REQUIRE_THAT(node_of_objective_func.OneThreeCycleComponentRight(), Catch::Approx(std::vector<double>({-61 * sqrt(3) / 3, 56, -7 * sqrt(6) / 3})));

	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0,0) == -903);
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1,0) == Approx(-2255*sqrt(3)/3));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2,0) == Approx(5042*sqrt(6)/3));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0,1) == Approx(476*sqrt(3)));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1,1) == Approx(-2428/3.0));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2,1) == Approx(-15668*sqrt(2)/3));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0,2) == Approx(-131*sqrt(6)));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1,2) == Approx(-1037*sqrt(2)/3));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2,2) == Approx(4246/3.0));

	// Check rotated node. 1 is chosen because it requires using all of the transposition matrices (1 2), (2 3), and (3 4).

	QAPSolver::Tests::NodeTester rotate_node_to_1{ node_of_objective_func, 1 };

	REQUIRE(rotate_node_to_1.TwoTwoCycleComponentLeft().size() == 0);

	REQUIRE(rotate_node_to_1.TwoTwoCycleComponentRight().size() == 0);

	REQUIRE_THAT(rotate_node_to_1.OneThreeCycleComponentLeft(), Catch::Approx(std::vector<double>({-27*sqrt(6)/2})));

	REQUIRE_THAT(rotate_node_to_1.OneThreeCycleComponentRight(), Catch::Approx(std::vector<double>({ -7 * sqrt(6) / 3 })));

	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(0, 0) == Approx(1806));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(1, 0) == Approx(1723*sqrt(3)));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(0, 1) == Approx(-1605*sqrt(3)));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(1, 1) == Approx(-3644));

	REQUIRE(rotate_node_to_1.IdentityComponent() == Approx(-95));

}


TEST_CASE("Computing the Fourier Transform of a QAP objective function of size 6 and its children") {
	QAPSolver::SquareMatrix left_mat(6);

	left_mat(0, 0) = 0;
	left_mat(0, 1) = -12;
	left_mat(0, 2) = 4;
	left_mat(0, 3) = 28;
	left_mat(0, 4) = -4;
	left_mat(0, 5) = -22;
	left_mat(1, 0) = 14;
	left_mat(1, 1) = 0;
	left_mat(1, 2) = -1;
	left_mat(1, 3) = 14;
	left_mat(1, 4) = -30;
	left_mat(1, 5) = 4;
	left_mat(2, 0) = -2;
	left_mat(2, 1) = -13;
	left_mat(2, 2) = 0;
	left_mat(2, 3) = 9;
	left_mat(2, 4) = -14;
	left_mat(2, 5) = -6;
	left_mat(3, 0) = 32;
	left_mat(3, 1) = 28;
	left_mat(3, 2) = -2;
	left_mat(3, 3) = 0;
	left_mat(3, 4) = 6;
	left_mat(3, 5) = 26;
	left_mat(4, 0) = -4;
	left_mat(4, 1) = -24;
	left_mat(4, 2) = -11;
	left_mat(4, 3) = 3;
	left_mat(4, 4) = 0;
	left_mat(4, 5) = -6;
	left_mat(5, 0) = 0;
	left_mat(5, 1) = 4;
	left_mat(5, 2) = 24;
	left_mat(5, 3) = 13;
	left_mat(5, 4) = -27;
	left_mat(5, 5) = 0;

	QAPSolver::SquareMatrix right_mat(6);

	right_mat(0, 0) = 0;
	right_mat(0, 1) = 6;
	right_mat(0, 2) = -29;
	right_mat(0, 3) = -5;
	right_mat(0, 4) = -2;
	right_mat(0, 5) = -25;
	right_mat(1, 0) = -5;
	right_mat(1, 1) = 0;
	right_mat(1, 2) = -29;
	right_mat(1, 3) = -14;
	right_mat(1, 4) = -27;
	right_mat(1, 5) = 2;
	right_mat(2, 0) = 11;
	right_mat(2, 1) = -11;
	right_mat(2, 2) = 0;
	right_mat(2, 3) = 13;
	right_mat(2, 4) = -4;
	right_mat(2, 5) = 3;
	right_mat(3, 0) = 19;
	right_mat(3, 1) = 22;
	right_mat(3, 2) = -31;
	right_mat(3, 3) = 0;
	right_mat(3, 4) = -21;
	right_mat(3, 5) = 28;
	right_mat(4, 0) = -12;
	right_mat(4, 1) = 9;
	right_mat(4, 2) = -24;
	right_mat(4, 3) = 32;
	right_mat(4, 4) = 0;
	right_mat(4, 5) = -29;
	right_mat(5, 0) = -8;
	right_mat(5, 1) = 12;
	right_mat(5, 2) = 26;
	right_mat(5, 3) = -17;
	right_mat(5, 4) = 20;
	right_mat(5, 5) = 0;




	//QAPSolver::Tests::SparseFTOfMatrixFuncTests left_tester(left_mat);
	//QAPSolver::Tests::SparseFTOfMatrixFuncTests right_tester(right_mat);

	QAPSolver::Tests::NodeTester node_of_objective_func(left_mat, right_mat);

	// Check Construction of node from matrices

	// Identity Component
	REQUIRE(node_of_objective_func.IdentityComponent() == -66960);
	//  REQUIRE_THAT( computed, Catch::Approx(known) );
	REQUIRE_THAT(node_of_objective_func.TwoTwoCycleComponentLeft(), Catch::Approx(std::vector<double>({ 1075.17, 883.659, - 33.5659, - 186, - 518.768, 660.989, - 36.5148, - 30.9839, - 643.988 })).margin(1));

	REQUIRE_THAT(node_of_objective_func.TwoTwoCycleComponentRight(), Catch::Approx(std::vector<double>({ 45.5895, 56.7501, - 15.7071, - 10.8333, 7.45356, 15.0616, 9.73729, 10.328, 0.745356 })).margin(1));

	REQUIRE_THAT(node_of_objective_func.OneThreeCycleComponentLeft(), Catch::Approx(std::vector<double>({ -258.042, 385.597, - 633.62, - 702, 334.626, - 37.5659, - 83.4841, 528, - 602.754, 339.411 })).margin(1));

	REQUIRE_THAT(node_of_objective_func.OneThreeCycleComponentRight(), Catch::Approx(std::vector<double>({ -14.7046, 44.0917, - 65.4534, 62.5, 33.3077, 6.4846, 53.4425, 9, - 38.9711, - 23.3345 })).margin(1));

	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0, 0) == Approx(-26334).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0, 1) == Approx(98695.7).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0, 2) == Approx(-28673.7).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0, 3) == Approx(26574.5).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(0, 4) == Approx(2844.32).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1, 0) == Approx(5913.22).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1, 1) == Approx(-30186).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1, 2) == Approx(-28264.5).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1, 3) == Approx(839.111).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(1, 4) == Approx(-29269.2).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2, 0) == Approx(8513.2).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2, 1) == Approx(7375.83).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2, 2) == Approx(179157).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2, 3) == Approx(-41911.1).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(2, 4) == Approx(139240).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(3, 0) == Approx(-40718.4).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(3, 1) == Approx(126677).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(3, 2) == Approx(-156477).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(3, 3) == Approx(63084.6).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(3, 4) == Approx(-88119.9).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(4, 0) == Approx(-8602.67).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(4, 1) == Approx(35819.1).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(4, 2) == Approx(6105.73).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(4, 3) == Approx(5646.56).margin(1));
	REQUIRE(node_of_objective_func.OneTwoCycleComponent()(4, 4) == Approx(13694.4).margin(1));


	// Check rotated node. 1 is chosen because it requires using all of the transposition matrices (1 2), (2 3), and (3 4).

	QAPSolver::Tests::NodeTester rotate_node_to_1{ node_of_objective_func, 1 };


	REQUIRE_THAT(rotate_node_to_1.TwoTwoCycleComponentLeft(), Catch::Approx(std::vector<double>({ 69.7653, - 108.444, - 129.263, 297.445, 53.6656 })).margin(1));

	REQUIRE_THAT(rotate_node_to_1.TwoTwoCycleComponentRight(), Catch::Approx(std::vector<double>({ 7.45356, 15.0616, 9.73729, 10.328, 0.745356 })).margin(1));

	REQUIRE_THAT(rotate_node_to_1.OneThreeCycleComponentLeft(), Catch::Approx(std::vector<double>({ 160.083, - 149.071, - 156.006, - 13.3333, - 34.641, - 174.42})).margin(1));

	REQUIRE_THAT(rotate_node_to_1.OneThreeCycleComponentRight(), Catch::Approx(std::vector<double>({ 33.3077, 6.4846, 53.4425, 9, - 38.9711, - 23.3345 })).margin(1));

	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(0, 0) == Approx(262).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(1, 0) == Approx(8787.27).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(2, 0) == Approx(-4657.3).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(3, 0) == Approx(-18638.5).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(0, 1) == Approx(-6923.58).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(1, 1) == Approx(11676.7).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(2, 1) == Approx(34388.5).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(3, 1) == Approx(-29732.2).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(0, 2) == Approx(-15831.5).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(1, 2) == Approx(26906.6).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(2, 2) == Approx(-43898.2).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(3, 2) == Approx(20120.8).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(0, 3) == Approx(10728).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(1, 3) == Approx(-3681.61).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(2, 3) == Approx(20185.3).margin(1));
	REQUIRE(rotate_node_to_1.OneTwoCycleComponent()(3, 3) == Approx(-5598.5).margin(1));

	REQUIRE(rotate_node_to_1.IdentityComponent() == Approx(-23076).margin(1));
}