/** 
 * alignment
 * @author astaroth
 * @date 10/10/20
 *
 * @brief
 */

#ifndef FA_ALIGNMENT_HPP
#define FA_ALIGNMENT_HPP

#include <string>
#include <vector>

struct WeightTable {
    const long insert_cost;
    const long delete_cost;
    const long substitution_cost;
    const long match_cost;
};

struct AlignerResults {
    double score {0.};
    std::string aligned_a, aligned_b;
};

#endif //FA_ALIGNMENT_HPP
