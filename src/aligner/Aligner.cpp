/** 
 * Aligner
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#include "Aligner.hpp"
#include <malloc.h>

void Aligner::align_genomes_bulk(const size_t &starting_index, const size_t &amount) {
    assert(starting_index + amount <= query_seqs.size());

    for (std::size_t i {starting_index}; i < (starting_index + amount); ++i) {
        /// Construct a Query object
        Query query {&reference, query_seqs.at(i), c_gap, protein_threshold, similarity_threshold, weights};

        /// Align the Query object
        query.search_proteins();
        query.search_islands();

        /// Finalize the Query sequence
        finalize_query(query);
        query_seqs.at(i).second = query.sequence;
        
        /// Report the amount of genomes processed
        std::cout << "\r[" << ++counter << '/' << query_seqs.size() << ']';
        std::cout.flush();
    }
}

void Aligner::parallel_aligner(const unsigned short &threads_number) {
    // CHECK: std::threads or std::async?
    std::size_t genomes_amount {query_seqs.size()};
    if (threads_number > 0) {
        unsigned short thread_n {static_cast<unsigned short>((threads_number > genomes_amount) ? genomes_amount : threads_number)};
        std::thread threading[thread_n];
        std::size_t portion_per_thread {genomes_amount / thread_n}, processed {0};
        for (unsigned short i {0}; i < thread_n; ++i) {
            if (i == thread_n - 1)
                portion_per_thread += (genomes_amount - processed - portion_per_thread);
            threading[i] = std::thread(&Aligner::align_genomes_bulk, this, processed, portion_per_thread);
            processed += portion_per_thread;
        }
        for (unsigned short i {0}; i < thread_n; ++i) {
            threading[i].join();
        }
    } else {
        align_genomes_bulk(0, genomes_amount);
    }
}

void Aligner::finalize_query(Query &query) const {
    /// This has only to be done if it is not a whole island
    if (not query.is_island) {
        /// Create now the final sequence, to which we will add fragments progressively
        std::string new_sequence {};

        /// 1. Pad the start
        new_sequence += std::string(query.start_pad, c_gap);

        /// 2. Merge `fragments` into the new sequence
        auto frag_iter{query.fragments.begin()};
        while (frag_iter != query.fragments.end()) {
            if (not frag_iter->second.second.empty()) {  // CHECK: different from .empty()?
                new_sequence += frag_iter->second.second;
            } else {
                new_sequence += query.sequence.substr(frag_iter->first, frag_iter->second.first - frag_iter->first + 1);
            }
            ++frag_iter;
        }

        /// 3. Pad the end
        new_sequence += std::string(query.end_pad, c_gap);

        /// 4. Save up some memory by clearing `fragments`
        query.fragments.clear();

        /// 5. Assign the sequence
        query.sequence.swap(new_sequence);
        query.sequence.shrink_to_fit();
    } /// If it is a whole island, we don't need to do anything - it has already been done
}

void Aligner::snps_find(const size_t &index) {
    assert(index < query_seqs.size());
    /// Point to the sequences we evaluate
    std::string *query_ptr {&query_seqs.at(index).first}, *ref_ptr {&reference.sequence};
    /// Point to the name of the query sequence
    std::string *name {&query_seqs.at(index).second};

    /// Calculate proper limits
    auto pred = [](const auto &c) { return c != '-'; };
    // CHECK VALIDITY
    std::size_t starter {static_cast<size_t>(std::distance(query_ptr->begin(),
                                                           std::find_if(query_ptr->begin(), query_ptr->end(), pred)))};
    std::size_t ender {static_cast<size_t>(std::distance(query_ptr->rbegin(),
                                                           std::find_if(query_ptr->rbegin(), query_ptr->rend(), pred)))};

    /// Initialize the snp collector
    std::vector<std::string> snp_collector;
    /// They have the same length, since they have already been aligned
    for (std::size_t i {starter}; i < ender; ++i) {
        if (query_ptr->at(i) != ref_ptr->at(i)) {
            if (query_ptr->at(i) == 'N' or ref_ptr->at(i) == 'N' or query_ptr->at(i) == '-' or ref_ptr->at(i) == '-')
                continue;
            mutations.push_back({*name, std::to_string(i), std::string(1, query_ptr->at(i)),
                                 std::string(1, ref_ptr->at(i))});
        }
    }
}

void Aligner::snps_bulk(const size_t &starting_index, const size_t &amount) {
    assert(starting_index + amount <= query_seqs.size());

    for (std::size_t i {starting_index}; i < (starting_index + amount); ++i) {
        snps_find(i);
        std::cout << '+';
    }
}

void Aligner::snps(const unsigned short &threads_number) {
    // CHECK: std::threads or std::async?
    std::size_t genomes_amount {query_seqs.size()};
    if (threads_number > 0) {
        unsigned short thread_n {static_cast<unsigned short>((threads_number > genomes_amount) ? genomes_amount : threads_number)};
        std::thread threading[thread_n];
        std::size_t portion_per_thread {genomes_amount / thread_n}, processed {0};
        for (unsigned short i {0}; i < thread_n; ++i) {
            if (i == thread_n - 1)
                portion_per_thread += (genomes_amount - processed - portion_per_thread);
            threading[i] = std::thread(&Aligner::snps_bulk, this, processed, portion_per_thread);
            processed += portion_per_thread;
        }
        for (unsigned short i {0}; i < thread_n; ++i) {
            threading[i].join();
        }
    } else {
        snps_bulk(0, genomes_amount);
    }
    std::cout << '\n';
}

void Aligner::perform_alignment(const unsigned short &threads_number) {
    /// Perform the alignment
    auto align_start {std::chrono::high_resolution_clock::now()};
    std::cout << "Aligning the sequences...\n";
    parallel_aligner(threads_number);
    auto align_end {std::chrono::high_resolution_clock::now()};
    auto align_elapsed {std::chrono::duration_cast<std::chrono::milliseconds>(align_end - align_start)};
    std::cout << "\nAligned in: " << align_elapsed.count() << " ms\n";
    /// Save up memory by clearing `reference.fragments`
    std::cout << "Clearing up data...\n";
    reference.fragments.clear();
    /// We are left with just the sequences
    // TODO: Fix all sequences (incl. reference) with reference's padding
    //  (hint: query to pad = query.start_pad - ref.start_pad)
    //std::map< std::size_t, std::pair<std::size_t, std::string> >().swap(reference.fragments);
    /// Write the output on the file
    exit(231);
    std::cout << "Writing the output sequences...\n";
    out_file.clear();
    out_file.write_results(reference, query_seqs);
    /*
    /// Find the SNPs on each query
    std::cout << "Checking for SNPs...\n";
    snps(threads_number);
    /// Write the SNPs on the file
    std::cout << "Writing the SNPs to a file...\n";
    snps_file.clear();
    snps_file.write_vsf(snp_header, mutations, '\t');
    */
}