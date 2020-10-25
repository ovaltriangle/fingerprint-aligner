/** 
 * FileManager
 * @author astaroth
 * @date 10/5/20
 *
 * @brief
 */

#ifndef FA_FILEMANAGER_HPP
#define FA_FILEMANAGER_HPP

#include "File.hpp"
#include "../exceptions/FileExceptions.hpp"

#include "formats/FASTA.hpp"
#include "formats/SNPS.hpp"

#include <memory>
#include <unordered_map>

#include "../genome/Reference.hpp"
#include "../genome/Query.hpp"

class FileManager {
private:
    std::unique_ptr<File> managed_file;
    friend class Aligner;
public:
    /// constructor
    explicit FileManager(const std::string &filename);
    /// member methods
    // read
    [[nodiscard]] std::pair<std::string, std::string> read_one_shot() const;
    [[nodiscard]] std::vector< std::pair<std::string, std::string> > read_sequences() const;
    // write
    void write_vsf(const std::vector<std::string> &cat, const std::vector< std::vector<std::string> > &data,
                                const char &sep = '\t') const;
    void write_comments(const std::vector<std::string> &comments) const;
    void write_results(const Reference &reference, const std::vector<Query> &queries) const;
    // control
    void clear() const { managed_file->clear(); };
};


#endif //FA_FILEMANAGER_HPP
