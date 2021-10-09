#ifndef QAP_MATRICES_DIMENSIONMISMATCH_H_
#define QAP_MATRICES_DIMENSIONMISMATCH_H_

#include <stdexcept>
#include <array>

namespace QAPSolver {
	class DimensionMismatch : public std::runtime_error {
	public:
		const std::array<int, 2> first_matrix_sizes_;
		const std::array<int, 2> second_matrix_sizes_;

		DimensionMismatch(int a, int b, int c, int d) : first_matrix_sizes_({ a, b }), second_matrix_sizes_({ c, d }),
			runtime_error("Precondition failed: Input matrices are not square or have different dimensions sizes.") {}
	};
}

#endif