/** 
 * Query
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#include "Query.hpp"

void Query::search_proteins() {
    /// Find all ORFs, independently from their RFs
    auto orfs {find_orfs()};

    /// If we don't find any ORF, it is an island
    if (orfs.empty()) {
        set_whole_island();
        return;  // exit early
    }

    /*
     * We don't know much about these ORFs. We may get a different ORF from the Reference if:
     * 1) The sequence presents deletions;
     * 2) The sequence is not aligned to the reference (most, if not all, times).
     *
     * We have one way to solve this issue: test all of them against the reference proteins.
     * By doing this, we are also aligning the proteins.
     */

    /// Iterate over the ORF of the reference
    auto ref_iter {reference->fragments.begin()};
    /// Iterate over the ORF of the queries
    auto query_iter{orfs.begin()};  // no need to reset at each cycle since sequences are ordered

    /// Use block coverage indices
    std::size_t coverage_start {std::string::npos}, coverage_end {0};
    /// Track the first and the last used query_iter and ref_iter instances
    std::size_t q_first {std::string::npos}, q_last {}, r_first {}, r_last {};

    /**
     * We should not order by start but by order and length. This proves necessary to find "the best CDS to check" as
     * the smallest CDS inside an island. This is crucial, since hirschberg's algorithm will take a lot of memory
     */
    while (ref_iter != reference->fragments.end()) {
        /// Check coverage margins for the reference
        if (not (ref_iter->first > coverage_start and ref_iter->second.first < coverage_end)) {
            std::string protein{reference->sequence.substr(ref_iter->first, ref_iter->second.first - ref_iter->first + 1)};
            while (query_iter != orfs.end()) {
                /// Check coverage margins for the query
                if (not(query_iter->first > coverage_start and query_iter->second.first < coverage_end)) {
                    std::string putative{
                        sequence.substr(query_iter->first, query_iter->second.first - query_iter->first + 1)};
                    if (putative.size() == protein.size()) {
                        double score {hirschberg::score(protein, putative, weights)};
                        if (score >= similarity) {
                            /// Save the first and last proteins to be added
                            if (q_first == std::string::npos) {  // set together; one check is enough
                                q_first = query_iter->first; r_first = ref_iter->first;
                            }
                            q_last = query_iter->second.first; r_last = ref_iter->second.first;

                            /// Find gaps in the sequence, before the start and after the end
                            std::size_t gap_before{sequence.rfind(c_gap, query_iter->first)},
                                        gap_after{sequence.find(c_gap, query_iter->second.first)};
                            /// Calculate the true start and end of the protein block
                            std::size_t true_start{(gap_before != std::string::npos) ? query_iter->first : 0};
                            std::size_t true_end{(gap_after != std::string::npos) ? query_iter->second.first : sequence.size() - 1};
                            /// Update coverage so we skip already covered proteins
                            coverage_start = (true_start < coverage_start) ? true_start : coverage_start;
                            coverage_end = (true_end > coverage_end) ? true_end : coverage_end;

                            /// Store the protein block we are covering
                            fragments.insert(std::make_pair(true_start, std::make_pair(true_end, "")));
                        }
                    }
                }
                ++query_iter;
            }
        }
        ++ref_iter;
    }

    /// Calculate the difference in length
    long start_diff = static_cast<long>(r_first - q_first);
    long end_diff = static_cast<long>(((reference->sequence.size() - 1) - r_last) - ((sequence.size() - 1) - q_last));
    /// Set the padding correctly
    if (start_diff > 0)  // start
        start_pad = start_diff;
    else
        if (-start_diff > reference->start_pad)
            reference->start_pad = -start_diff;
    if (end_diff > 0)  // end
        end_pad = end_diff;
    else
        if (-end_diff > reference->end_pad)
            reference->end_pad = -end_diff;

    /// Save up some memory
    orfs.clear();

    /// Put up some flags
    if (fragments.empty())
        set_whole_island();
    else
        set_proteins_aligned();
}

void Query::search_islands() {
    /// Check the flags to know where we are at
    if (is_island) {
        /// Align the whole sequence against the whole reference sequence
        AlignerResults ar{hirschberg::score_full(reference->sequence, sequence, weights)};
        sequence = ar.aligned_b;
        /// Shrink the string
        sequence.shrink_to_fit();
    } else if (proteins_aligned) {
        /// Iterate over the proteins and check if we have an island between them
        auto protein_before {fragments.begin()}, protein_after {++fragments.begin()};
        while (protein_after != fragments.end()) {
            /// Set up some variables to be more concise
            auto pb_end {protein_before->second.first}, pa_begin {protein_after->first};

            /// Cut the sequence between the proteins
            std::string query_bw_proteins {sequence.substr(pb_end, pa_begin - pb_end + 1)};
            /// Cut the sequence from the reference
            std::string corresponding_reference {reference->sequence.substr(pb_end + start_pad,
                                                                            pa_begin - pb_end + end_pad + 1)};

            /// Align the sequence and add it to the fragments
            AlignerResults ar {hirschberg::score_full(corresponding_reference, query_bw_proteins, weights)};
            fragments.insert(std::make_pair(pb_end + 1, std::make_pair(pa_begin - 1, ar.aligned_b)));

            /// Increment the iterators
            ++protein_before; ++protein_after;
        }
    } else {
        /// In this case, we have something missing or the method has been called too early
        throw QueryNotAligned();
    }

    set_islands_aligned();
}
