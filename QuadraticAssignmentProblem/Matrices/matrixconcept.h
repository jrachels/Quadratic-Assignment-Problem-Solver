// TO DO: Matrix concept should allow you to specific the required member functions (accessor and size) and not apply specific names of them.
// even better if they can have common default names like size1 if the size or accessor functions aren't specified.

#ifndef QAP_MATRICES_MATRIXCONCEPT_H_
#define QAP_MATRICES_MATRIXCONCEPT_H_

#include <concepts>

namespace QAPSolver {

	template<typename T>
	concept MatrixConcept = requires (T a, int b, int c) { //intellisense issue
		{a(b, c)} -> std::convertible_to<double>;
		{a.size1()}->std::convertible_to<size_t>;
		{a.size2()}->std::convertible_to<size_t>;
	};

}

#endif