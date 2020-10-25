/** 
 * SNPS
 * @author astaroth
 * @date 10/14/20
 *
 * @brief
 */

#ifndef FA_SNPS_HPP
#define FA_SNPS_HPP

#include "../File.hpp"

#include <string>
#include <utility>

class SNPS : public File
{
public:
    explicit SNPS(std::string filename)
    : File(std::move(filename), ' ', ' ', ' ', 0)
    {};
};

#endif //FA_SNPS_HPP
