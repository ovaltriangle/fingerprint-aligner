/** 
 * Genome
 * @author astaroth
 * @date 10/8/20
 *
 * @brief
 */

#ifndef FA_GENOME_HPP
#define FA_GENOME_HPP

#include <string>
#include <array>
#include <utility>
#include <map>

#include <boost/algorithm/string.hpp>

#include "../utils/searcher.hpp"

class Genome {
protected:
    /// constructor
    std::string name, sequence;
    std::size_t threshold;
    bool is_rna;
    std::string start_codon;
    std::array<std::string, 3> end_codons;
    /// !constructor
    std::map< std::size_t, std::pair<std::size_t, std::string> > fragments;  // proteins + islands - aligned version
    std::size_t start_pad {0}, end_pad {0};  // useful when a query starts/ends before/after the reference
    /// flags
    bool proteins_aligned {false}, islands_aligned {false}, is_island {false};
    /// everyone needs a friend
    friend class FileManager;
    friend class Aligner;
    friend class Reference;
    friend class Query;
public:
    /// constructors
    Genome(std::pair<std::string, std::string> ns, const std::size_t &protein_threshold)
    : name(std::move(ns.first)), sequence(std::move(ns.second)), threshold(protein_threshold)
    {
        // make the sequence all uppercase to avoid issues
        boost::to_upper(sequence);
        sequence.shrink_to_fit();
        is_rna = sequence.find('U') != std::string::npos;
        set_codons();  // may be reused if we convert a RNA into a DNA and vice versa
    }
    /// member methods
    // setters
    void set_codons();
    // flags
    void set_proteins_aligned() { proteins_aligned = true; }
    void set_islands_aligned() { islands_aligned = true; }
    void set_whole_island() { is_island = true; }
    // analyzers
    std::map< std::size_t, std::pair<std::size_t, std::string> > find_orfs();
    /// abstract methods
    virtual void search_proteins() = 0;
    virtual void search_islands() = 0;
};


#endif //FA_GENOME_HPP
