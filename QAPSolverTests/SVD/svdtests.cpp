#include "QAPSolverTests/catch.hpp"

#include "QuadraticAssignmentProblem/SVD/squaresvd.h"
#include "QuadraticAssignmentProblem/SVD/tracenormsumofdyadicmatrices.h"

TEST_CASE("Trace Norm of 4 by 4 matrix") {
    QAPSolver::SquareMatrix square_mat(4);

    square_mat(0, 0) = -16;
    square_mat(0, 1) = 23;
    square_mat(0, 2) = 25;
    square_mat(0, 3) = -29;
    square_mat(1, 0) = -13;
    square_mat(1, 1) = -23;
    square_mat(1, 2) = 24;
    square_mat(1, 3) = -20;
    square_mat(2, 0) = -25;
    square_mat(2, 1) = -11;
    square_mat(2, 2) = -14;
    square_mat(2, 3) = -6;
    square_mat(3, 0) = -12;
    square_mat(3, 1) = 2;
    square_mat(3, 2) = 4;
    square_mat(3, 3) = 27;

    REQUIRE(QAPSolver::internal::SquareSVD(square_mat).TraceNorm() == Approx(145.0189));
}

TEST_CASE("Trace Norm of Dyadic Matrices") {
    QAPSolver::Matrix left_mat(5, 3);
    QAPSolver::Matrix right_mat(3, 5);
    QAPSolver::SquareMatrix  product_mat(5);

    left_mat(0, 0) = -30;
    left_mat(0, 1) = -32;
    left_mat(0, 2) = -2;
    left_mat(1, 0) = -28;
    left_mat(1, 1) = 11;
    left_mat(1, 2) = -5;
    left_mat(2, 0) = 9;
    left_mat(2, 1) = -24;
    left_mat(2, 2) = 23;
    left_mat(3, 0) = 16;
    left_mat(3, 1) = 25;
    left_mat(3, 2) = -25;
    left_mat(4, 0) = 6;
    left_mat(4, 1) = -17;
    left_mat(4, 2) = 29;

    right_mat(0, 0) = 29;
    right_mat(0, 1) = -21;
    right_mat(0, 2) = -30;
    right_mat(0, 3) = -20;
    right_mat(0, 4) = 32;
    right_mat(1, 0) = -5;
    right_mat(1, 1) = 10;
    right_mat(1, 2) = -29;
    right_mat(1, 3) = 7;
    right_mat(1, 4) = 2;
    right_mat(2, 0) = 23;
    right_mat(2, 1) = 15;
    right_mat(2, 2) = -7;
    right_mat(2, 3) = 30;
    right_mat(2, 4) = 25;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 3; k++) {
                product_mat(i, j) += left_mat(i, k) * right_mat(k, j);
            }
        }
    }

    REQUIRE(QAPSolver::internal::TraceNormSumOfDyadicMatrices(left_mat, right_mat) == Approx(QAPSolver::internal::SquareSVD(product_mat).TraceNorm()));

}
