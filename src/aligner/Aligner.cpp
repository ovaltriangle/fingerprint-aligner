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
    assert(starting_index + amount <= queries.size());

    for (std::size_t i {starting_index}; i < (starting_index + amount); ++i) {
        /// Report the amount of genomes processed
        std::cout << "\r[" << ++counter << '/' << queries.size() << ']';
        std::cout.flush();

        /// Construct a Query object
        Query *query {&queries.at(i)};

        /// Align the Query object
        query->search_proteins();
        query->search_islands();

        /// Finalize the Query sequence
        finalize_query(*query);
        queries.at(i).name = query->sequence;
    }
}

void Aligner::parallel_aligner(const unsigned short &threads_number) {
    // CHECK: std::threads or std::async?
    std::size_t genomes_amount {queries.size()};
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

void Aligner::finalize_ref() {
    /// Pad the reference's sequence
    reference.sequence = std::string(reference.start_pad, c_gap) + reference.sequence + std::string(reference.end_pad, c_gap);
}

void Aligner::snps_find(const size_t &index) {
    assert(index < queries.size());
    /// Point to the sequences we evaluate
    std::string *query_ptr {&queries.at(index).sequence}, *ref_ptr {&reference.sequence};
    /// Point to the name of the query sequence
    std::string *name {&queries.at(index).name};

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
    assert(starting_index + amount <= queries.size());

    for (std::size_t i {starting_index}; i < (starting_index + amount); ++i) {
        std::cout << "\r[" << ++counter << '/' << queries.size() << ']';
        std::cout.flush();

        snps_find(i);
    }
}

void Aligner::snps(const unsigned short &threads_number) {
    // CHECK: std::threads or std::async?
    std::size_t genomes_amount {queries.size()};
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
    counter = 0;
    parallel_aligner(threads_number);
    auto align_end {std::chrono::high_resolution_clock::now()};
    auto align_elapsed {std::chrono::duration_cast<std::chrono::milliseconds>(align_end - align_start)};
    std::cout << "\nAligned in: " << align_elapsed.count() << " ms\n";
    /// Save up memory by clearing `reference.fragments`
    std::cout << "Clearing up data...\n";
    reference.fragments.clear();
    /// We are left with just the sequences
    finalize_ref();
    /// Write the output on the file
    exit(231);
    std::cout << "Writing the output sequences...\n";
    out_file.clear();
    out_file.write_results(reference, queries);
    /*
    /// Find the SNPs on each query
    std::cout << "Checking for SNPs...\n";
    counter = 0;
    snps(threads_number);
    /// Write the SNPs on the file
    std::cout << "Writing the SNPs to a file...\n";
    snps_file.clear();
    snps_file.write_vsf(snp_header, mutations, '\t');
    */
}