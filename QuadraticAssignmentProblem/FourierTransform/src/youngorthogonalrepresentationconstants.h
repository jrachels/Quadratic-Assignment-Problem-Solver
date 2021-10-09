#ifndef YoungOrthogonalRepresentationConstants_H_
#define YoungOrthogonalRepresentationConstants_H_

#include <vector>
#include <tuple>
#include <array>

namespace QAPSolver {
	// Collection of constants that are entries of the
	// Young Orthogonal Representation Matrices. Entry [i]
	// holds the constants corresponding to an axial distance
	// of i. Entry [0] is (-infinity, -infinity), but unused.
	class YoungOrthogonalRepresentationConstants
	{
	private:
		std::vector< std::array< double, 2>> constants_;
	public:
		explicit YoungOrthogonalRepresentationConstants(int i) {
			GrowYORConstants(i);
		}

		YoungOrthogonalRepresentationConstants() = default;

		// Computes YOR constants required for any current input matrix
		// and stores for later use in yor_constants_.
		// TO DO: consider order of operations for double precision arithmetic 
		// TO DO: consider immediately converting i to double so repeated implicit conversion
		// isn't necessary
		YoungOrthogonalRepresentationConstants& GrowYORConstants(int i);

		double current_size();

		double operator()(int i, int j) const {
			return constants_[i][j];
		}


		[[nodiscard]] double diagonal(int i) const {
			//return std::get<0>(constants_[i]);// [0] ;
			//return (constants_[i])[0];
			return std::get<0>(constants_[i]);
		}

		[[nodiscard]] double off_diagonal(int i) const {
			return std::get<1>(constants_[i]);
			//return std::get<1>(constants_[i]);// [1] ;
			//return (constants_[i])[1];
		}
	};
}

#endif

