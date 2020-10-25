/** 
 * GenomeExceptions
 * @author astaroth
 * @date 10/9/20
 *
 * @brief
 */

#ifndef FA_GENOMEEXCEPTIONS_HPP
#define FA_GENOMEEXCEPTIONS_HPP

#include <exception>

/// Reference

struct ReferenceIslands : public std::exception
{
    ReferenceIslands() noexcept = default;
    ~ReferenceIslands() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "A reference genome may not contain islands.";
    }
};

struct QueryNotAligned : public std::exception
{
    QueryNotAligned() noexcept = default;
    ~QueryNotAligned() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "We cannot search for islands if the proteins are not aligned.";
    }
};

#endif //FA_GENOMEEXCEPTIONS_HPP
