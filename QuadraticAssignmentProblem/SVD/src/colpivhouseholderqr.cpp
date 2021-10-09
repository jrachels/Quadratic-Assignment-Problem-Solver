
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "QuadraticAssignmentProblem/SVD/colpivhouseholderqr.h"

// ONLY WORKS IF MORE ROWS THAN COLUMNS

std::vector<int> QAPSolver::internal::ColPivHouseHolderQR(Matrix& m_qr) {

    int rows = m_qr.size1();
    int cols = m_qr.size2();
    assert(rows >= cols);
    int size = std::min<int>(rows, cols);

    std::vector<double> m_hCoeffs(size); // do i need to keep track of this?

    std::vector<int> what_location_is_this_col_index(cols);
    std::vector<int> what_col_is_in_this_col_location(cols);
    for (int i = 0; i < cols; ++i) {
        what_col_is_in_this_col_location[i] = i;
    }
    int number_of_transpositions = 0;

    std::vector<double> m_colNormsUpdated(cols);
    std::vector<double> m_colNormsDirect(cols);
    double max_norm = 0;
    for (int k = 0; k < cols; ++k) {
        double norm = internal::NormOfColumn(m_qr, k);
        m_colNormsDirect[k] = norm;
        m_colNormsUpdated[k] = norm;
        if (norm > max_norm) {
            max_norm = norm;
        }
    }

    //std::max_element(m_colNormsUpdated.begin(), m_colNormsUpdated.end()); Why did i write this?
    double threshold_helper = std::pow((max_norm * std::numeric_limits<double>::epsilon()) / static_cast<double>(rows), 2);
    double norm_downdate_threshold = std::sqrt(std::numeric_limits<double>::epsilon());

    int m_nonzero_pivots = size; 
    double m_maxpivot = 0;

    for (int k = 0; k < size; ++k)
    {
       
        int biggest_col_index = k;

        double biggest_col_sq_norm = 0; // unsquared
        for (int i = k; i < cols; i++) {
            double value = m_colNormsUpdated[i];
            if (biggest_col_sq_norm < value) {
                biggest_col_sq_norm = value; // unsquare
                biggest_col_index = i;
            }
        }
        biggest_col_sq_norm = biggest_col_sq_norm * biggest_col_sq_norm; // finally squared

        //double biggest_col_sq_norm = std::pow((m_colNormsUpdated.tail(cols - k).maxCoeff(&biggest_col_index)), 2);
        //biggest_col_index += k; NO LONGER NECESSARY

        if (m_nonzero_pivots == size && biggest_col_sq_norm < threshold_helper * static_cast<double>(rows - k)) {
            m_nonzero_pivots = k;
        }

        // apply the transposition to the columns
        int current_col_index_origin = what_col_is_in_this_col_location[k];
        int big_col_index_origin = what_col_is_in_this_col_location[biggest_col_index];
        what_location_is_this_col_index[big_col_index_origin] = k; //wrong
        what_location_is_this_col_index[current_col_index_origin] = biggest_col_index;
        what_col_is_in_this_col_location[k] = big_col_index_origin;
        what_col_is_in_this_col_location[biggest_col_index] = current_col_index_origin;

        if (k != biggest_col_index) {

            for (int i = 0; i < rows; i++) {
                double temp = m_qr(i, k);
                m_qr(i, k) = m_qr(i, biggest_col_index);
                m_qr(i, biggest_col_index) = temp;
            }


            std::swap(m_colNormsUpdated[k], m_colNormsUpdated[biggest_col_index]);
            std::swap(m_colNormsDirect[k], m_colNormsDirect[biggest_col_index]);
            ++number_of_transpositions;
        }

        // generate the householder vector, store it below the diagonal
        double beta;
        makeHouseHolderInColumn(m_qr, k, m_hCoeffs[k], beta);
        // apply the householder transformation to the diagonal coefficient
        m_qr(k, k) = beta;

        // remember the maximum absolute value of diagonal coefficients
        if (abs(beta) > m_maxpivot) m_maxpivot = abs(beta);

        applyHouseHolderOnLeft(m_qr, k, m_hCoeffs[k]);

        for (int j = k + 1; j < cols; ++j) {
            if (m_colNormsUpdated[j] != 0) {
                double temp = abs(m_qr(k, j)) / m_colNormsUpdated[j];
                temp = (1 + temp) * (1 - temp);
                temp = temp < 0 ? 0 : temp;
                double temp2 = temp * std::pow((m_colNormsUpdated[j] /
                    m_colNormsDirect[j]), 2);
                if (temp2 <= norm_downdate_threshold) {

                    m_colNormsDirect[j] = internal::NormOfTailOfColumn(m_qr, j, k + 1);
                    m_colNormsUpdated[j] = m_colNormsDirect[j];
                }
                else {
                    m_colNormsUpdated[j] *= std::sqrt(temp);
                }
            }
        }
    }

    return what_location_is_this_col_index;

}





