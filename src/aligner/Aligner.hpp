/** 
 * Aligner
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#ifndef FA_ALIGNER_HPP
#define FA_ALIGNER_HPP


#include "../genome/Reference.hpp"
#include "../genome/Query.hpp"

#include <string>
#include <vector>
#include <thread>
#include <iostream>
#include <chrono>

#include "../file/FileManager.hpp"
//#include "../utils/threader.hpp"

class Aligner {
private:
    /// constructor
    Reference reference;
    std::vector<Query> queries;
    FileManager out_file, snps_file;
    char c_gap;
    /// !constructor
    std::vector<std::string> snp_header {"Name", "Position", "ReferenceAllele", "QueryAllele"};
    std::vector< std::vector<std::string> > mutations;
    std::size_t counter {0};  // how many query_seqs have been processed
public:
    /// constructors
    Aligner(const std::string &reference_file, const std::string &query_file, const std::string &output_file,
            const std::string & snp_file, const std::size_t &protein_threshold = 60,
            const double &similarity_threshold = 0.90, const WeightTable &weights = {-2, -2, -1, 2})
    : reference(FileManager(reference_file).read_one_shot(), protein_threshold), out_file(output_file), snps_file(snp_file)
    {
        /// Create the FileManager object
        FileManager query_manager {query_file};
        /// Assign the gap character
        c_gap = query_manager.managed_file->c_gap;
        /// Read the query sequences
        for (const auto &seq: query_manager.read_sequences()) {
            Query q {&reference, seq, c_gap, protein_threshold, similarity_threshold, weights};
            queries.push_back(q);
        }
        /// Initialize the Reference
        initialize_reference();
    }
    /// member methods
    // reference
    void initialize_reference() { reference.search_proteins(); }
    // multi threaders
    void align_genomes_bulk(const std::size_t &starting_index, const std::size_t &amount);
    void parallel_aligner(const unsigned short &threads_number);
    // finalizers
    void finalize_query(Query &query) const;
    void finalize_ref();
    // snps
    void snps_find(const std::size_t &index);
    void snps_bulk(const std::size_t &starting_index, const std::size_t &amount);
    void snps(const unsigned short &threads_number);
    // conveyor
    void perform_alignment(const unsigned short &threads_number);
};


#endif //FA_ALIGNER_HPP
