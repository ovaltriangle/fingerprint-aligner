/** 
 * Query
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#ifndef FA_QUERY_HPP
#define FA_QUERY_HPP

#include "Genome.hpp"
#include "Reference.hpp"

#include "../aligner/methods/hirschberg.hpp"

#include <cstdlib>

/// MAKE THIS SUBCLASS AGAIN

class Query : public Genome {
private:
    /// constructor
    Reference *reference;
    double similarity;
    WeightTable weights;
    char c_gap;
    /// !constructor
    std::size_t s_pad_acc {0}, e_pad_acc {0};  // pad accounted for
    /// friends
    friend class Aligner;
    friend class FileManager;
public:
    /// constructors
    Query(Reference *reference, std::pair<std::string, std::string> ns, const char &c_gap, const std::size_t &protein_threshold,
          const double &similarity_threshold, const WeightTable &weights)
    : Genome(std::move(ns), protein_threshold), reference(reference), similarity(similarity_threshold), weights(weights), c_gap(c_gap)
    {}
    /// abstract implementation
    void search_proteins() override;
    void search_islands() override;
    /// fix ends
    void fix_ends();
};


#endif //FA_QUERY_HPP
