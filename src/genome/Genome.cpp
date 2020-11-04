/** 
 * Genome
 * @author astaroth
 * @date 10/8/20
 *
 * @brief
 */

#include "Genome.hpp"

void Genome::set_codons()
{
    if (is_rna) {  // RNA
        start_codon = "AUG";
        end_codons = {"UAG", "UAA", "UGA"};
    } else {  // DNA
        start_codon = "ATG";
        end_codons = {"TAG", "TAA", "TGA"};
    }
}

std::map< std::size_t, std::pair<std::size_t, std::string> > Genome::find_orfs()
{
    /// Find all starting positions
    std::set<std::size_t> starting_positions {utils::find_all(sequence, start_codon)};

    /// Find all ending positions and construct a container for them
    std::set<std::size_t> ending_position0 {utils::find_all(sequence, end_codons.at(0))};
    std::set<std::size_t> ending_position1 {utils::find_all(sequence, end_codons.at(1))};
    std::set<std::size_t> ending_position2 {utils::find_all(sequence, end_codons.at(2))};

    /// Unify them all under one vector
    std::set<std::size_t> ending_positions;
    if (not ending_position0.empty())
        ending_positions.merge(ending_position0);

    if (not ending_position1.empty())
        ending_positions.merge(ending_position1);

    if (not ending_position2.empty())
        ending_positions.merge(ending_position2);

    /// Free up some memory
    ending_position0.clear();
    ending_position1.clear();
    ending_position2.clear();

    /// Find the ORFs
    std::map< std::size_t, std::pair<std::size_t, std::string> > orfs;
    std::vector<std::size_t> seen_ends;
    for (const auto &start: starting_positions) {
        std::size_t closest {std::string::npos};
        // find the closest valid end to the start
        for (const auto &end: ending_positions) {
            if (end < start or (end - start) % 3 != 0 or (end - start) < threshold) {
                continue;
            } else {
                closest = end;
                break;
            }
        }

        // keep going if there was no closing end
        if (closest == std::string::npos)
            continue;

        // update the container
        if (not utils::contains(seen_ends, closest)) {
            orfs.insert(std::make_pair(start, std::make_pair(closest + 2, "")));
            seen_ends.push_back(closest);
        }
    }

    /// Free up some more memory
    starting_positions.clear();
    ending_positions.clear();

    /// Check if we actually have some ORFs
    if (orfs.empty()) {
        set_whole_island();
        return std::map< std::size_t, std::pair<std::size_t, std::string> >();
    }

    /// Return the found ORFs
    return orfs;  // at this point we have a list of ORFs for RF 0, 1, 2
}