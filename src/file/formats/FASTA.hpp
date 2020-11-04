/**
 * FASTA
 * @brief Header files for the well known FASTA format. Easiest to integrate.
 *
 * @author astaroth
 * @version 1.0 
 * @date 7/15/20
 */

#ifndef FP_ALIGN_FASTA_HPP
#define FP_ALIGN_FASTA_HPP

#include "../File.hpp"

#include <string>
#include <utility>

/**
 * All operations inherited by File will do just fine
 */
class FASTA : public File
{
public:
    explicit FASTA(std::string filename)
    : File(std::move(filename), '#', '-', '>', 60)
    {};
};

#endif //FP_ALIGN_FASTA_HPP
