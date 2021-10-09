#ifndef QAP_MATRICES_NOTSQUARE_H_
#define QAP_MATRICES_NOTSQUARE_H_

#include <stdexcept>

namespace QAPSolver {
	class NotSquare : public std::runtime_error {
	public:
		const char* what() const noexcept
		{
			return "Precondition failed: Input matrix not Square.";
		}
	};
}

#endif
