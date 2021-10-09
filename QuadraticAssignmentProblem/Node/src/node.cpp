
#include "QuadraticAssignmentProblem/Node/node.h"

double QAPSolver::FourierTransform::internal::factorial(int n) {
	assert(n >= 0);
	if (n == 0) {
		return 1;
	}
	else {
		return n * factorial(n - 1);
	}
}

QAPSolver::FourierTransform::Node::Node(const Node& extended_ft, int index_of_rotation) : size(extended_ft.size - 1) {
	int parent_size = extended_ft.size;
	//size = extended_ft.size - 1;
	assert(index_of_rotation <= parent_size);

	FTAtIdentityTableau rotated_ft_at_identity_tableau = extended_ft.ft_at_identity_tableau_.Rotate(index_of_rotation, parent_size);
	FTAtOneTwoCycleTableau rotated_ft_at_one_two_cycle_tableau = extended_ft.ft_at_one_two_cycle_tableau_.Rotate(index_of_rotation, parent_size);
	FTAtTwoTwoCycleTableau rotated_ft_at_two_two_cycle_tableau = extended_ft.ft_at_two_two_cycle_tableau_.Rotate(index_of_rotation, parent_size);
	FTAtOneThreeCycleTableau rotated_ft_at_one_three_cycle_tableau_ = extended_ft.ft_at_one_three_cycle_tableau_.Rotate(index_of_rotation, parent_size);

	ft_at_identity_tableau_ = FTAtIdentityTableau::FindChildComponentFromRotations(rotated_ft_at_identity_tableau, rotated_ft_at_one_two_cycle_tableau, parent_size);
	ft_at_one_two_cycle_tableau_ = FTAtOneTwoCycleTableau::FindChildComponentFromRotations(rotated_ft_at_one_two_cycle_tableau, rotated_ft_at_two_two_cycle_tableau, rotated_ft_at_one_three_cycle_tableau_, parent_size);
	ft_at_two_two_cycle_tableau_ = FTAtTwoTwoCycleTableau::FindChildComponentFromRotations(rotated_ft_at_two_two_cycle_tableau, parent_size);
	ft_at_one_three_cycle_tableau_ = FTAtOneThreeCycleTableau::FindChildComponentFromRotations(rotated_ft_at_one_three_cycle_tableau_, parent_size);
}

double QAPSolver::FourierTransform::Node::Bound() const {

	double total_bound = 0;

	total_bound += ft_at_identity_tableau_.Bound_Contribution();
	total_bound += ft_at_one_two_cycle_tableau_.Bound_Contribution();
	total_bound += ft_at_two_two_cycle_tableau_.Bound_Contribution();
	total_bound += ft_at_one_three_cycle_tableau_.Bound_Contribution();

	return total_bound / internal::factorial(size);

}

double QAPSolver::FourierTransform::Node::Value() const {
	if (size == 1) {
		return ft_at_identity_tableau_.ft_identity_tableau_;
	}
	else {
		throw std::logic_error("Attempted to extract Value from Node of size > 1.\n");
	}
}

/*******************************************************/
// BEGIN TWO TWO CYCLE *******************************/
/*******************************************************/

QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau::FTAtTwoTwoCycleTableau(const SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau& ft_matrix_two_two_cycle_1, const SparseFTOfMatrixFunc::FFTAtTwoTwoCycleTableau& ft_matrix_two_two_cycle_2, double missing_factor) {
	ft_two_two_cycle_column_ = ft_matrix_two_two_cycle_1.sparse_ft_at_two_two_cycle_tableau_;
	int temp_length = ft_two_two_cycle_column_.size();
	for (int i = 0; i < temp_length; i++) {
		ft_two_two_cycle_column_[i] *= missing_factor; // choose this component since it will be the one rotated
	}
	ft_two_two_cycle_row_ = ft_matrix_two_two_cycle_2.sparse_ft_at_two_two_cycle_tableau_;
}

QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau::FindChildComponentFromRotations(FTAtTwoTwoCycleTableau& ft_parent, int parent_size) {
	// default construct a new FTAtOneThreeCycleTableau
	FTAtTwoTwoCycleTableau child{};
	// compute size of vects
	if (parent_size > 4) {
		int child_length = std::max((parent_size - 1) * (parent_size - 4) / 2, 0);

		// set size of vects in FTAtOneThreeCycleTableau
		child.ft_two_two_cycle_column_.resize(child_length);
		child.ft_two_two_cycle_row_.resize(child_length);
		//  iterate over corresponding elements of ft_parent and multiply by the constant 
		double multiplicative_constant = (parent_size - 3) / static_cast<double>(((parent_size - 1) * (parent_size - 4)));
		int shift = parent_size - 2;
		for (int i = 0; i < child_length; i++) {
			child.ft_two_two_cycle_column_[i] = multiplicative_constant * (ft_parent.ft_two_two_cycle_column_[i + shift]);
			child.ft_two_two_cycle_row_[i] = (ft_parent.ft_two_two_cycle_row_[i + shift]); // only multiple by multiplicative constant on one matrix
		}
	}
	return child;
}

double QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau::Bound_Contribution() const {
	//double bound_contribution = 0;
	double column_bound_temp = 0;
	double row_bound_temp = 0;
	for (double x : ft_two_two_cycle_column_) {
		column_bound_temp += x * x;
	}
	for (double x : ft_two_two_cycle_row_) {
		row_bound_temp += x * x; // TO DO: every child of a parent gets  the same row bound for this component. I don't need to recompute it. 
	}
	return (ft_two_two_cycle_column_.size()) * (std::sqrt(column_bound_temp) * std::sqrt(row_bound_temp));
}

QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau::Rotate(const int index_of_rotation, const int size) const {
	// create a copy
	FTAtTwoTwoCycleTableau copy{};
	copy.ft_two_two_cycle_column_ = (this->ft_two_two_cycle_column_);
	copy.ft_two_two_cycle_row_ = (this->ft_two_two_cycle_row_);
	// apply rotation
	internal::ApplyRotationTranspose(copy, size, index_of_rotation, size);
	// return copy
	return copy;
}

void QAPSolver::FourierTransform::Node::FTAtTwoTwoCycleTableau::ApplyAdjacentTranspositionTranspose(int size, int i) {
	// for each column

	assert(size > i);
	// do op for row1:[1]row2:[2]
	if (i == 1) {
		int index = 0;
		for (int c = size - 2; c > 1; index += c--) {
			double temp = ft_two_two_cycle_column_[index];
			ft_two_two_cycle_column_[index] = -temp;
		}
	}
	else {
		{ // entering scope
			// do op for row2:[i][] and row2:[i+1][] pairs
			int index = i - 2; // index of [i][n]
			for (int c = size - 2; c > i - 1; index += c--) {
				//save
				double entry_i = ft_two_two_cycle_column_[index]; // row2:[i]row3:[]
				double entry_i_1 = ft_two_two_cycle_column_[index + 1]; //row2:[i+1]row3:[]
				//swap
				ft_two_two_cycle_column_[index] = entry_i * yor_constants_.diagonal(i) + entry_i_1 * yor_constants_.off_diagonal(i);
				ft_two_two_cycle_column_[index + 1] = entry_i * yor_constants_.off_diagonal(i) - entry_i_1 * yor_constants_.diagonal(i);
			}
		} // exiting scope
		// do op for row1:[1][3]row2:[2][4]
		if (i == 3) {
			int index = ((size) * (size - 3) / 2) - 2;
			double temp = ft_two_two_cycle_column_[index];
			ft_two_two_cycle_column_[index] = -temp;
		}
		// do op for row2:[][i] and row2:[][i+1] pairs
		else if (i > 3) {
			const int a = (size + i - 3) * (size - i) / (2); // points to entry  row2:[2][i]
			const int b = (size + i - 2) * (size - i - 1) / (2); // points to entry row2:[2][i+1]
			for (int c = 0; c < i - 2; c++) {
				// save
				double entry_i = ft_two_two_cycle_column_[a + c]; // row2:[][i]
				double entry_i_1 = ft_two_two_cycle_column_[b + c]; // row3:[][i+1]
				// do op
				ft_two_two_cycle_column_[a + c] = entry_i * yor_constants_.diagonal(i-2) + entry_i_1 * yor_constants_.off_diagonal(i-2);
				ft_two_two_cycle_column_[b + c] = entry_i * yor_constants_.off_diagonal(i-2) - entry_i_1 * yor_constants_.diagonal(i-2);
			}
		}
	}
}

/*******************************************************/
// BEGIN ONE THREE CYCLE *******************************/
/*******************************************************/

QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau::FTAtOneThreeCycleTableau(const SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau& ft_matrix_one_three_cycle_1, const SparseFTOfMatrixFunc::FFTAtOneThreeCycleTableau& ft_matrix_one_three_cycle_2, double missing_factor) {
	ft_one_three_cycle_column_ = ft_matrix_one_three_cycle_1.sparse_ft_at_one_three_cycle_tableau_;
	int temp_length = ft_one_three_cycle_column_.size();
	for (int i = 0; i < temp_length; i++) {
		ft_one_three_cycle_column_[i] *= missing_factor; // choose this component since it will be the one rotated
	}
	ft_one_three_cycle_row_ = ft_matrix_one_three_cycle_2.sparse_ft_at_one_three_cycle_tableau_;
}

QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau::FindChildComponentFromRotations(FTAtOneThreeCycleTableau& ft_parent, int parent_size) {
	// default construct a new FTAtOneThreeCycleTableau
	FTAtOneThreeCycleTableau child{};
	// compute size of vects
	if (parent_size > 3) {
		int child_length = (parent_size - 2) * (parent_size - 3) / 2;

		// set size of vects in FTAtOneThreeCycleTableau
		child.ft_one_three_cycle_column_.resize(std::max(child_length, 0));
		child.ft_one_three_cycle_row_.resize(std::max(child_length, 0));
		//  iterate over corresponding elements of ft_parent and multiply by the constant 
		double multiplicative_constant = (parent_size - 1) / static_cast<double>((parent_size * (parent_size - 3)));
		int shift = parent_size - 2;
		for (int i = 0; i < child_length; i++) {
			child.ft_one_three_cycle_column_[i] = multiplicative_constant * (ft_parent.ft_one_three_cycle_column_[i + shift]);
			child.ft_one_three_cycle_row_[i] = (ft_parent.ft_one_three_cycle_row_[i + shift]); // only multiple by multiplicative constant on one matrix
		}
	}
	//else {
	//	child.ft_one_three_cycle_column_.resize(std::max(child_length, 0));
	//	child.ft_one_three_cycle_row_.resize(std::max(child_length, 0));
	//}
	return child;
}

double QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau::Bound_Contribution() const {
	//double bound_contribution = 0;
	double column_bound_temp = 0;
	double row_bound_temp = 0;
	for (double x : ft_one_three_cycle_column_) {
		column_bound_temp += x * x;
	}
	for (double x : ft_one_three_cycle_row_) {
		row_bound_temp += x * x; // TO DO: every child of a parent gets  the same row bound for this component. I don't need to recompute it. 
	}
	//bound_contribution += std::sqrt(column_bound_temp) + std::sqrt(row_bound_temp);
	return (ft_one_three_cycle_row_.size()) * (std::sqrt(column_bound_temp) * std::sqrt(row_bound_temp));
}

QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau::Rotate(const int index_of_rotation, const int size) const {
	// create a copy
	FTAtOneThreeCycleTableau copy{};
	copy.ft_one_three_cycle_column_ = (this->ft_one_three_cycle_column_);
	copy.ft_one_three_cycle_row_ = (this->ft_one_three_cycle_row_);
	// apply rotation
	internal::ApplyRotationTranspose(copy, size, index_of_rotation, size);
	// return copy
	return copy;
}


// Issue: this code is partially duplicated from the fftofmatrixfunc class
void QAPSolver::FourierTransform::Node::FTAtOneThreeCycleTableau::ApplyAdjacentTranspositionTranspose(int size, int i) {

	assert(size > i);
	if (size > 2) {
// do op for row2:[1]row3:[2]
if (i == 1) {
	int index = 0;
	for (int c = size - 2; c > 0; index += c--) {
		double temp = ft_one_three_cycle_column_[index];
		ft_one_three_cycle_column_[index] = -temp;
	}
}
else {
	{// entering scope
		// do op for row2:[]row3:[i] and row2:[]row3:[i+1] pairs
		const int a = (size + i - 3) * (size - i) / (2); // points to entry  row2:[2]row3[i]
		const int b = (size + i - 2) * (size - i - 1) / (2); // points to entry row2:[2]row3[i+1]
		for (int c = 0; c < i - 2; c++) {
			// save
			double entry_i = ft_one_three_cycle_column_[a + c]; // row3:[i]
			double entry_i_1 = ft_one_three_cycle_column_[b + c]; // row3:[i+1]
			// do op
			ft_one_three_cycle_column_[a + c] = entry_i * yor_constants_.diagonal(i) + entry_i_1 * yor_constants_.off_diagonal(i);
			ft_one_three_cycle_column_[b + c] = entry_i * yor_constants_.off_diagonal(i) - entry_i_1 * yor_constants_.diagonal(i);
		}

		// do op for row2:[i]row3:[i+1]
		double temp = ft_one_three_cycle_column_[b + i - 2];
		ft_one_three_cycle_column_[b + i - 2] = -temp;
	}// exiting scope
	{
		// do op for row2:[i]row3:[] and row2:[i+1]row3:[] pairs
		int index = i - 2; // index of [i][n]
		for (int c = size - 2; c > i - 1; index += c--) {
			//save
			double entry_i = ft_one_three_cycle_column_[index]; // row2:[i]row3:[]
			double entry_i_1 = ft_one_three_cycle_column_[index + 1]; //row2:[i+1]row3:[]
			//swap
			ft_one_three_cycle_column_[index] = entry_i * yor_constants_.diagonal(i) + entry_i_1 * yor_constants_.off_diagonal(i);
			ft_one_three_cycle_column_[index + 1] = entry_i * yor_constants_.off_diagonal(i) - entry_i_1 * yor_constants_.diagonal(i);
		}
	}
}
	}
}

/*******************************************************/
// BEGIN ONE TWO CYCLE *******************************/
/*******************************************************/


QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau::FTAtOneTwoCycleTableau(const SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau& ft_matrix_one_two_cycle_1, const SparseFTOfMatrixFunc::FFTAtOneTwoCycleTableau& ft_matrix_one_two_cycle_2, double missing_factor) {
	int size = std::get<0>(ft_matrix_one_two_cycle_1.sparse_ft_at_one_two_cycle_tableau_).size();

	ft_one_two_cycle_ = SquareMatrix(size);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			ft_one_two_cycle_(i, j) = std::get<0>(ft_matrix_one_two_cycle_1.sparse_ft_at_one_two_cycle_tableau_)[i] * std::get<0>(ft_matrix_one_two_cycle_2.sparse_ft_at_one_two_cycle_tableau_)[j];
			ft_one_two_cycle_(i, j) += std::get<1>(ft_matrix_one_two_cycle_1.sparse_ft_at_one_two_cycle_tableau_)[i] * std::get<1>(ft_matrix_one_two_cycle_2.sparse_ft_at_one_two_cycle_tableau_)[j];
			ft_one_two_cycle_(i, j) *= missing_factor;
		}
	}

}

double QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau::Bound_Contribution() const {
	return (ft_one_two_cycle_.size1()) * ft_one_two_cycle_.TraceNorm();
}

QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau::FindChildComponentFromRotations(FTAtOneTwoCycleTableau& ft_one_two_cycle_parent, Node::FTAtTwoTwoCycleTableau& ft_two_two_cycle_parent,
	Node::FTAtOneThreeCycleTableau& ft_one_three_cycle_parent, int parent_size) {

	// default construct a new FTAtOneTwoCycleTableau
	FTAtOneTwoCycleTableau child{};
	// compute size of vects
	int child_length = (parent_size - 2); // parent_size is n

	child.ft_one_two_cycle_ = SquareMatrix(std::max(child_length, 0));

	// determine multiplicative constants
	// spread this out over multiple for loops to be cache friendly.

	// parent_size = n
	double multiplicative_constant_one_two_cycle;
	if (parent_size > 2) {
		multiplicative_constant_one_two_cycle = (parent_size - 1) / static_cast<double>(((parent_size) * (parent_size - 2)));
	}
	else {
		multiplicative_constant_one_two_cycle = (parent_size - 1) / (static_cast<double>(parent_size)); // this is irrelevant because the child has length 0? TO DO: remove else?
	}

	for (int i = 0; i < child_length; i++) {
		for (int j = 0; j < child_length; j++) {
			child.ft_one_two_cycle_(i, j) = multiplicative_constant_one_two_cycle * ft_one_two_cycle_parent.ft_one_two_cycle_(i, j);
		}
	}

	if (parent_size > 3) { // should be "if (parent_size > 3)"?

		double multiplicative_constant_two_two_cycle = (parent_size - 3) / static_cast<double>(2 * (parent_size - 2));
		// spread this out over multiple for loops to be cache friendly

		for (int i = 0; i < child_length; i++) {
			for (int j = 0; j < child_length; j++){
				child.ft_one_two_cycle_(i, j) += multiplicative_constant_two_two_cycle * ft_two_two_cycle_parent.ft_two_two_cycle_column_[i] * \
					ft_two_two_cycle_parent.ft_two_two_cycle_row_[j];
			}
		}
	}
	if (parent_size > 2) {
		double  multiplicative_constant_one_three_cycle = (parent_size - 1) / static_cast<double>(2 * parent_size);

		for (int i = 0; i < child_length; ++i) {
			for (int j = 0; j < child_length; ++j) {
				child.ft_one_two_cycle_(i, j) += multiplicative_constant_one_three_cycle * ft_one_three_cycle_parent.ft_one_three_cycle_column_[i] * \
					ft_one_three_cycle_parent.ft_one_three_cycle_row_[j];
			}
		}
	}

	return child;
}

// Issue: may produce too many copies. may produce 1 less copy if i reorder the basis elements in fftofmatrixfunc
// should produce a smaller ft?
QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau::Rotate(const int index_of_rotation, const int size) const {
	// create a copy
	FTAtOneTwoCycleTableau copy{};
	copy.ft_one_two_cycle_ = (this->ft_one_two_cycle_);
	// apply rotation
	internal::ApplyRotationTranspose(copy, size, index_of_rotation, size);
	// return copy

	return copy;
}

// Issue: this code is partially duplicated from the fftofmatrixfunc class
void QAPSolver::FourierTransform::Node::FTAtOneTwoCycleTableau::ApplyAdjacentTranspositionTranspose(int size, int k) {
	// TO DO: This was written when FTAtOneTwoCycleTableau was stored as a product of matrices. It can be optimized now.
	for (int j = 0; j < ft_one_two_cycle_.size(); j++) {
		if (k == 1) {
			if (size > 1) {
				double temp = ft_one_two_cycle_(0, j);
				ft_one_two_cycle_(0, j) = -temp;
			}
		}
		else {
			if (size > 2) {
				double a = ft_one_two_cycle_(k - 2, j);
				double b = ft_one_two_cycle_(k - 1, j);
				ft_one_two_cycle_(k - 2, j) = yor_constants_.diagonal(k) * a + yor_constants_.off_diagonal(k) * b;
				ft_one_two_cycle_(k - 1, j) = yor_constants_.off_diagonal(k) * a - yor_constants_.diagonal(k) * b;
			}
		}
	}
}

/*******************************************************/
// BEGIN IDENTITY COMPONENT ****************************/
/*******************************************************/

QAPSolver::FourierTransform::Node::FTAtIdentityTableau QAPSolver::FourierTransform::Node::FTAtIdentityTableau::FindChildComponentFromRotations(FTAtIdentityTableau& ft_identity_tableau_parent, Node::FTAtOneTwoCycleTableau& ft_one_two_cycle_parent, int parent_size) {
	// default construct a new FTAtIdentityTableau
	FTAtIdentityTableau child{};
	// initialize to 0
	child.ft_identity_tableau_ = 0;
	child.ft_identity_tableau_ += ft_identity_tableau_parent.ft_identity_tableau_ / (parent_size);
	//child.ft_identity_tableau_ += (static_cast<double>(parent_size - 1) / static_cast<double>(parent_size)) * ft_one_two_cycle_parent.ft_one_two_cycle_(parent_size - 2, parent_size - 2);
	if (parent_size > 0) {
		child.ft_identity_tableau_ += (static_cast<double>(parent_size - 1) / static_cast<double>(parent_size)) * ft_one_two_cycle_parent.ft_one_two_cycle_(parent_size - 2, parent_size - 2);

	}
	return child;
}

// Issue: may produce too many copies. may produce 1 less copy if i reorder the basis elements in fftofmatrixfunc
// should produce a smaller ft?
QAPSolver::FourierTransform::Node::FTAtIdentityTableau QAPSolver::FourierTransform::Node::FTAtIdentityTableau::Rotate(const int index_of_rotation, const int size) const {
	// create a copy
	FTAtIdentityTableau copy{};
	copy.ft_identity_tableau_ = (this->ft_identity_tableau_);
	// apply rotation
	internal::ApplyRotationTranspose(copy, size, index_of_rotation, size);
	// return copy
	return copy;
}