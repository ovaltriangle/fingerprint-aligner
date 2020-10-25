/** 
 * searcher
 * @author astaroth
 * @date 10/10/20
 *
 * @brief
 */

#ifndef FA_SEARCHER_HPP
#define FA_SEARCHER_HPP

#include <string>
#include <vector>
#include <set>
#include <algorithm>

namespace utils {

    std::set<std::size_t> find_all(const std::string& searched, const std::string& to_find);

    /*
     * SO: 6194797
     */
    template <typename Container>
    bool contains(Container const &c, typename Container::const_reference v) {
        return std::find(c.begin(), c.end(), v) != c.end();
    }

    template <typename Container>
    std::size_t argmax(const Container &c) {
        return static_cast<size_t>(std::distance(c.begin(), std::max_element(c.begin(), c.end())));
    }

}  // namespace utils


#endif //FA_SEARCHER_HPP
