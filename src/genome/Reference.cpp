/** 
 * Reference
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#include "Reference.hpp"

void Reference::search_proteins() {
    /// Find all ORFs, independently from their RFs
    auto orfs {find_orfs()};

    /// Find the ORF of greatest length
    std::vector<std::size_t> orf_sizes {0, 0, 0};
    orf_sizes.shrink_to_fit();  // CHECK: do not preserve memory
    for (const auto &orf: orfs)
        orf_sizes.at(orf.first % 3) += (orf.second.first - orf.first);
    // returns the index of the biggest element in `orf_sizes` == best ORF
    std::size_t index {utils::argmax(orf_sizes)};

    /// Create a map with only the CDS of the best ORF
    auto pred = [index] (const auto &p) { return p.first % 3 == index; };
    std::copy_if(orfs.begin(), orfs.end(), std::inserter(fragments, fragments.end()), pred);

    /// Save up some memory
    orfs.clear();
    //std::map< std::size_t, std::pair<std::size_t, std::string> >().swap(orfs);
}
