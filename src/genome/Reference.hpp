/** 
 * Reference
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#ifndef FA_REFERENCE_HPP
#define FA_REFERENCE_HPP

#include "Genome.hpp"
#include "../exceptions/GenomeExceptions.hpp"

#include <utility>

#include "../utils/searcher.hpp"

class Reference : public Genome {
private:
    friend class Query;
    friend class Aligner;
    friend class FileManager;
public:
    /// constructors
    Reference(std::pair<std::string, std::string> ns, const std::size_t &protein_threshold)
    : Genome(std::move(ns), protein_threshold)
    {
        // flags for the reference genome
        set_proteins_aligned();
        set_islands_aligned();
    }
    /// abstract implementations
    void search_proteins();
    void search_islands() { throw ReferenceIslands(); };
};


#endif //FA_REFERENCE_HPP
