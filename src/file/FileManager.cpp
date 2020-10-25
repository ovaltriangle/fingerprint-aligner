/** 
 * FileManager
 * @author astaroth
 * @date 10/5/20
 *
 * @brief
 */

#include "FileManager.hpp"

FileManager::FileManager(const std::string &filename) {
    /// unordered_map to derive the class to be used wrt files's extension
    std::unordered_map<std::string, int> format_chooser {
            {"fasta", 1}, {"faa", 1}, {"fna", 1}, {"fa", 1},
            {"snps", 2}
    };  // mapped decisions for each format

    /// Take the extension of the files
    std::size_t dot_pos {filename.rfind('.')};
    if (dot_pos == std::string::npos)
        // cannot infer format without extension
        throw FileWithNoFormat();
    // cut the extension
    std::string ext {filename.substr(dot_pos + 1)};

    /// Choose the format & assign the pointer
    switch (format_chooser[ext]) {
        case 1:  // FASTA
            managed_file = std::make_unique<FASTA>(filename);
            break;
        case 2:
            managed_file = std::make_unique<SNPS>(filename);
            break;
        default:
            throw UnhandledFormat();  // if we don't recognize the format
    }
}

std::pair<std::string, std::string> FileManager::read_one_shot() const {
    /// Reroute to the underneath method
    long position {0};  // we have to since we modify this
    return managed_file->read_sequence(position);
}

std::vector<std::pair<std::string, std::string> > FileManager::read_sequences() const {
    /// Reroute to the underneath method
    return managed_file->read_sequences();
}

void FileManager::write_vsf(const std::vector<std::string> &cat, const std::vector< std::vector<std::string> > &data,
                            const char &sep) const {
    /// Reroute to the underneath method
    managed_file->write_vsf(cat, data, sep);
}

void FileManager::write_comments(const std::vector<std::string> &comments) const {
    /// Reroute to the underneath method
    managed_file->write_comments(comments);
}

void FileManager::write_results(const Reference &reference, const std::vector<Query> &queries) const {
    /// Write the Reference
    managed_file->write_sequence(reference.name, reference.sequence);
    /// Write all Query genomes
    for (const auto &query: queries) {
        managed_file->write_sequence(query.name, query.sequence);
    }
}
