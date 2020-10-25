/** 
 * searcher
 * @author astaroth
 * @date 10/10/20
 *
 * @brief
 */

#include "searcher.hpp"

namespace utils {
    /**
     * Finds all occurrences of the string `to_find` in `searched` and returns a
     * sorted collection, a set, containing all the positions in which it was found.
     *
     * @param searched The string to be searched
     * @param to_find The string to be found
     * @return set of size_t containing all positions in which the string was found
     */
    std::set<std::size_t> find_all(const std::string& searched, const std::string& to_find) {
        // initialize variables for the while loop
        std::set<std::size_t> positions;
        // find the first position of the pattern we are searching
        std::size_t position {searched.find(to_find)};
        // keep going until you can't find it
        while (position != std::string::npos) {
            positions.insert(position);  // store the position
            // find the new position wrt the old one
            position = searched.find(to_find, position + 1);
        }

        // search completed
        return positions;
    }
}  // namespace utils