/** 
 * hirschberg
 * @author astaroth
 * @date 10/10/20
 *
 * @brief
 */

#ifndef FA_HIRSCHBERG_HPP
#define FA_HIRSCHBERG_HPP

#include "../../types/alignment.hpp"

#include <string>
#include <vector>
#include <cmath>

#include <thread>

#include "../../utils/searcher.hpp"

namespace hirschberg {

    std::vector<long> nw_score(const std::string &a, const std::string &b, const WeightTable &weights);
    std::vector<long> nw_score_h(const std::string &a, const std::string &b, const WeightTable &weights);
    std::pair<std::string, std::string> align(const std::string &a, const std::string &b, const WeightTable &weights);
    double score(const std::string &a, const std::string &b, const WeightTable &weights);
    AlignerResults score_full(const std::string &a, const std::string &b, const WeightTable &weights);

} // namespace hirschberg


#endif //FA_HIRSCHBERG_HPP
