#include "QuadraticAssignmentProblem/FTOfMatrixFunc/sparseftofmatrixfunc.h"

// Compute the relevant column of the (n-1) representation matrices 
// corresponding to the (n-1) transpositions.
// TO DO: move multiplication by particular adjacency matrix to helper function.
std::vector<std::vector<double>> QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTRowFuncs::RepresentationsOfTranspositions(int size) {
	std::vector<std::vector<double>> representations_of_transpositions(size - 1);  // there are (n-1) columns
	representations_of_transpositions[size - 2].resize(size - 2);  // Each column is of length (n-2)
	representations_of_transpositions[size - 2][size - 3] = 1;  // Representation of identity permutation is identity
	representations_of_transpositions[size - 3] = representations_of_transpositions[size - 2];  // Copy forward to next transposition

	// Computes the column of representation matrices for all transpositions of (n-1) with (i+1)
	// except for (i+1)=1, because this case requires special care. 
	// Transpositions are written as a product of adjacent transpositions.
	// (i j) with i<j becomes (j j-1)(j-1 j-2)...(i+2 i+1) (i+1 i) (i+1 i+2)...(j-2 j-1)(j-1 j)
	for (int i = size - 3; i > 0; i--) {
		// Multiply by representation matrix of adjacency transposition of elements i and (i+1).
		representations_of_transpositions[i][i - 1] = representations_of_transpositions[i][i] * (yor_constants_.off_diagonal(i + 1));
		representations_of_transpositions[i][i] = -representations_of_transpositions[i][i] * (yor_constants_.diagonal(i + 1));
		// smart copy forward
		representations_of_transpositions[i - 1] = representations_of_transpositions[i];
		// Multiply by remaining adjacency matrices
		for (int j = i; j < size - 3; j++) {
			double a = representations_of_transpositions[i][j];
			double b = representations_of_transpositions[i][j + 1];
			representations_of_transpositions[i][j] = (yor_constants_.diagonal(j + 2)) * a + (yor_constants_.off_diagonal(j + 2)) * b;
			representations_of_transpositions[i][j + 1] = -(yor_constants_.diagonal(j + 2)) * b + (yor_constants_.off_diagonal(j + 2)) * a;
		}
	}

	// Finish the final, special case. 
	representations_of_transpositions[0][0] = -representations_of_transpositions[0][0];
	for (int j = 0; j < size - 3; j++) {
		double a = representations_of_transpositions[0][j];
		double b = representations_of_transpositions[0][j + 1];
		representations_of_transpositions[0][j] = (yor_constants_.diagonal(j + 2)) * a + (yor_constants_.off_diagonal(j + 2)) * b;
		representations_of_transpositions[0][j + 1] = -(yor_constants_.diagonal(j + 2)) * b + (yor_constants_.off_diagonal(j + 2)) * a;
	}

	return representations_of_transpositions;
}

QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau::FFTAtOneTwoCycleTableau(FFTRowFuncs& ft_row_funcs, int size) {
	// Handle the first component
	std::get<0>(sparse_ft_at_one_two_cycle_tableau_).resize(size - 1);
	std::get<0>(sparse_ft_at_one_two_cycle_tableau_).shrink_to_fit();
	for (int i = 0; i < size; i++) {
		// copy entries to temporary
		std::vector<double> temporary(size - 1);
		for (int j = 0; j < size - 2; j++) {
			temporary[j] = std::get<0>(ft_row_funcs.ft_row_funcs_[i])[j];
		}
		// apply transpositions
		internal::ApplyTransposition<FFTAtOneTwoCycleTableau>(temporary, size, i + 1, size);
		// add to total
		for (int j = 0; j < size - 1; j++) {
			std::get<0>(sparse_ft_at_one_two_cycle_tableau_)[j] += temporary[j];
		}
	}

	// Handle the second component
	std::get<1>(sparse_ft_at_one_two_cycle_tableau_).resize(size - 1);
	std::get<1>(sparse_ft_at_one_two_cycle_tableau_).shrink_to_fit();
	for (int i = 0; i < size; i++) {
		// copy entries to temporary
		std::vector<double> temporary(size - 1);
		temporary[size - 2] = std::get<1>(ft_row_funcs.ft_row_funcs_[i]);
		// apply transpositions
		internal::ApplyTransposition<FFTAtOneTwoCycleTableau>(temporary, size, i + 1, size);
		// add to total
		for (int j = 0; j < size - 1; j++) {
			std::get<1>(sparse_ft_at_one_two_cycle_tableau_)[j] += temporary[j];
		}
	}
}

void QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau::ApplyAdjacentTransposition(std::vector<double>& intermediate_vec, int size, int i) {
	assert(size > i);
	if (i == 1) {
		double temp = intermediate_vec[0];
		intermediate_vec[0] = -temp;
	}
	else {
		double a = intermediate_vec[i - 2];
		double b = intermediate_vec[i - 1];
		intermediate_vec[i - 2] = (yor_constants_.diagonal(i)) * a + (yor_constants_.off_diagonal(i)) * b;
		intermediate_vec[i - 1] = (yor_constants_.off_diagonal(i)) * a - (yor_constants_.diagonal(i)) * b;
	}
}

QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau::FFTAtTwoTwoCycleTableau(FFTRowFuncs& ft_row_funcs, int size) {
	if (size > 3){
		const int ft_component_size = ((size) * (size - 3)) / 2;
		sparse_ft_at_two_two_cycle_tableau_.resize(ft_component_size);
		sparse_ft_at_two_two_cycle_tableau_.shrink_to_fit();
		for (int i = 0; i < size; i++) {
			// copy entries to temporary
			std::vector<double> temporary(ft_component_size);
			for (int j = 0; j < size - 2; j++) {
				temporary[j] = std::get<0>(ft_row_funcs.ft_row_funcs_[i])[j];
			}
			// apply transpositions
			internal::ApplyTransposition<FFTAtTwoTwoCycleTableau>(temporary, size, i + 1, size);
			// add to total
			for (int j = 0; j < ft_component_size; j++) {
				sparse_ft_at_two_two_cycle_tableau_[j] += temporary[j];
			}
		}
	}
}

// Modify vector in a way equivalent to multiplying by the representation matrix of 
// an adjacent transposition of element i with element i+1.
// TO DO: remove i==1 case from this function since it is only called at the end. handle outside of for loop.
void QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau::ApplyAdjacentTransposition(std::vector<double>& intermediate_vec, int size, int i) {
	assert(size > i);
	// do op for row1:[1]row2:[2]
	if (i == 1) {
		int index = 0;
		for (int c = size - 2; c > 1; index += c--) {
			double temp = intermediate_vec[index];
			intermediate_vec[index] = -temp;
		}
	}
	else {
		{ // entering scope
			// do op for row2:[i][] and row2:[i+1][] pairs
			int index = i - 2; // index of [i][n]
			for (int c = size - 2; c > i - 1; index += c--) {
				//save
				double entry_i = intermediate_vec[index]; // row2:[i]row3:[]
				double entry_i_1 = intermediate_vec[index + 1]; //row2:[i+1]row3:[]
				//swap
				intermediate_vec[index] = entry_i * (yor_constants_.diagonal(i)) + entry_i_1 * (yor_constants_.off_diagonal(i));
				intermediate_vec[index + 1] = entry_i * (yor_constants_.off_diagonal(i)) - entry_i_1 * (yor_constants_.diagonal(i));
			}
		} // exiting scope
		// do op for row1:[1][3]row2:[2][4]
		if (i == 3) {
			int index = ((size) * (size - 3) / 2) - 2;
			double temp = intermediate_vec[index];
			intermediate_vec[index] = -temp;
		}
		// do op for row2:[][i] and row2:[][i+1] pairs
		else if (i > 3) {
			const int a = (size + i - 3) * (size - i) / (2); // points to entry  row2:[2][i]
			const int b = (size + i - 2) * (size - i - 1) / (2); // points to entry row2:[2][i+1]
			for (int c = 0; c < i - 2; c++) {
				// save
				double entry_i = intermediate_vec[a + c]; // row2:[][i]
				double entry_i_1 = intermediate_vec[b + c]; // row3:[][i+1]
				// do op
				intermediate_vec[a + c] = entry_i * (yor_constants_.diagonal(i-2)) + entry_i_1 * (yor_constants_.off_diagonal(i-2));
				intermediate_vec[b + c] = entry_i * (yor_constants_.off_diagonal(i-2)) - entry_i_1 * (yor_constants_.diagonal(i-2));
			}
		}
	}
}


QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau::FFTAtOneThreeCycleTableau(FFTRowFuncs& ft_row_funcs, int size) {
	const int ft_component_size = ((size - 1) * (size - 2)) / 2;
	sparse_ft_at_one_three_cycle_tableau_.resize(ft_component_size);
	sparse_ft_at_one_three_cycle_tableau_.shrink_to_fit();
	for (int i = 0; i < size; i++) {
		// copy entries to temporary
		std::vector<double> temporary(ft_component_size);
		for (int j = 0; j < size - 2; j++) {
			temporary[j] = std::get<0>(ft_row_funcs.ft_row_funcs_[i])[j];
		}
		// apply transpositions
		internal::ApplyTransposition<FFTAtOneThreeCycleTableau>(temporary, size, i + 1, size);
		// add to total

		for (int j = 0; j < ft_component_size; j++) {
			sparse_ft_at_one_three_cycle_tableau_[j] += temporary[j];
		}
	}
}

// Modify vector in a way equivalent to multiplying by the representation matrix of 
// an adjacent transposition of element i with element i+1.
// TO DO: remove i==1 case from this function since it is only called at the end. handle outside of for loop.
void QAPSolver::FourierTransform::SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau::ApplyAdjacentTransposition(std::vector<double>& intermediate_vec, int size, int i) {
	assert(size > i);
	// do op for row2:[1]row3:[2]
	if (i == 1) {
		int index = 0;
		for (int c = size - 2; c > 0; index += c--) {
			double temp = intermediate_vec[index];
			intermediate_vec[index] = -temp;
		}
	}
	else {
		{// entering scope
			// do op for row2:[]row3:[i] and row2:[]row3:[i+1] pairs
			const int a = (size + i - 3) * (size - i) / (2); // points to entry  row2:[2]row3[i]
			const int b = (size + i - 2) * (size - i - 1) / (2); // points to entry row2:[2]row3[i+1]
			for (int c = 0; c < i - 2; c++) {
				// save
				double entry_i = intermediate_vec[a + c]; // row3:[i]
				double entry_i_1 = intermediate_vec[b + c]; // row3:[i+1]
				// do op
				intermediate_vec[a + c] = entry_i * (yor_constants_.diagonal(i)) + entry_i_1 * (yor_constants_.off_diagonal(i));
				intermediate_vec[b + c] = entry_i * (yor_constants_.off_diagonal(i)) - entry_i_1 * (yor_constants_.diagonal(i));
			}

			// do op for row2:[i]row3:[i+1]
			double temp = intermediate_vec[b + i - 2];
			intermediate_vec[b + i - 2] = -temp;
		}// exiting scope
		{
			// do op for row2:[i]row3:[] and row2:[i+1]row3:[] pairs
			int index = i - 2; // index of [i][n]
			for (int c = size - 2; c > i - 1; index += c--) {
				//save
				double entry_i = intermediate_vec[index]; // row2:[i]row3:[]
				double entry_i_1 = intermediate_vec[index + 1]; //row2:[i+1]row3:[]
				//swap
				intermediate_vec[index] = entry_i * (yor_constants_.diagonal(i)) + entry_i_1 * (yor_constants_.off_diagonal(i));
				intermediate_vec[index + 1] = entry_i * (yor_constants_.off_diagonal(i)) - entry_i_1 * (yor_constants_.diagonal(i));
			}
		}
	}
}