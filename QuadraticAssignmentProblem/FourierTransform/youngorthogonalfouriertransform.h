
#ifndef QAP_FOURIERTRANSFORM_YOUNGORTHOGONALFOURIERTRANSFORM_H_
#define QAP_FOURIERTRANSFORM_YOUNGORTHOGONALFOURIERTRANSFORM_H_

#include "QuadraticAssignmentProblem/FourierTransform/src/youngorthogonalrepresentationconstants.h"

namespace QAPSolver {
	namespace FourierTransform {

		// Caches some constant coefficients used in the computation of fourier transforms.

		class YoungOrthogonalFourierTransform {
		public:
			YoungOrthogonalFourierTransform(YoungOrthogonalRepresentationConstants& yor_constants) {
				if (yor_constants_.current_size() < yor_constants.current_size()) {
					yor_constants_ = yor_constants;
				}
			}
			YoungOrthogonalFourierTransform() = default;
		protected:
			static YoungOrthogonalRepresentationConstants yor_constants_;
		};
	}
}

#endif