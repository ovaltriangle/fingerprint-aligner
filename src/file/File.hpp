/** 
 * File
 * @author astaroth
 * @date 10/4/20
 *
 * @brief
 */

#ifndef FA_FILE_HPP
#define FA_FILE_HPP

#include "../exceptions/FileExceptions.hpp"

#include <string>
#include <tuple>
#include <vector>
#include <fstream>

class File {
protected:
    std::string file_name;
    char c_comm, c_gap, c_name;
    std::size_t c_per_line;
    friend class Aligner;
    friend class FileManager;
public:
    /// constructors
    File(std::string filename, const char &comment, const char &gap, const char &name, const std::size_t &chars_per_line)
    : file_name(std::move(filename)), c_comm(comment), c_gap(gap), c_name(name), c_per_line(chars_per_line)
    {};
    /// member methods
    // read
    std::pair<std::string, std::string> read_sequence(long &position) const;
    [[nodiscard]] std::vector< std::pair<std::string, std::string> > read_sequences() const;
    // write
    void write_one_line(const std::string &to_write) const;
    void write_several(const std::vector<std::string> &to_write) const;
    // users
    void write_vsf(const std::vector<std::string> &cat, const std::vector< std::vector<std::string> > &data,
                   const char &sep = '\t') const;
    void write_comments(const std::vector<std::string> &to_write) const;
    void write_sequence(const std::string &name, const std::string &sequence) const;
    // control
    void clear() const;
};


#endif //FA_FILE_HPP
