/** 
 * File
 * @author astaroth
 * @date 10/4/20
 *
 * @brief
 */

#include "File.hpp"

std::pair<std::string, std::string> File::read_sequence(long &position) const
{
    /// Open the files and check for its existence
    std::ifstream file {file_name};
    if (!file)
        throw FileNotFound();

    /// Set files position
    file.seekg(position);

    /// Read from the buffer
    std::string line {}, name {}, sequence {};
    while (getline(file, line)) {
        /// Skip comments and blank lines
        if (line.empty() or line.at(0) == c_comm) {
            continue;
        /// Preserve name and sequence
        } else if (line.at(0) == c_name) {
            if (not name.empty())
                break;
            name = line.substr(1);  // store the name
        } else {  // store the sequence
            sequence += line;
        }
    }

    /// Store the new position, then close the files
    position = file.tellg();
    if (position != -1)
        position -= static_cast<long>(line.size() + 1);  // reset cursors to the start of the line
    file.close();  // end of operations on the files

    /// If `sequence` is empty, the files has no sequences in it!
    if (sequence.empty())
        // no possibility of bad reads (see: `read_sequences`)
        throw NoSequencesInFile();

    return {name, sequence};
}

std::vector<std::pair<std::string, std::string> > File::read_sequences() const
{
    /// Create the vector object and initialize the position object
    std::vector< std::pair<std::string, std::string> > sequences;
    long position {0};  // will get modified at each iteration

    /// Read sequences and update the variables, until we reach EOF
    do {  // we suppose there is at least one sequence
        // read the sequence at position and push_back to vector
        sequences.push_back(read_sequence(position));
    } while (position != -1);  // -1 => EOF

    return sequences;
}

void File::write_one_line(const std::string &to_write) const {
    /// Open the files in append mode
    std::ofstream file {file_name, std::ios::app};
    /// Write the line to the file, then close it
    file << to_write << '\n';
    file.close();
}

void File::write_several(const std::vector<std::string> &to_write) const
{
    /// Open the files in append mode
    std::ofstream file {file_name, std::ios::app};
    /// Write each string to the files and then close it
    for (const std::string &line: to_write) {
        file << line << '\n';
    }
    file.close();
}

void File::write_vsf(const std::vector<std::string> &cat, const std::vector< std::vector<std::string> > &data,
                     const char &sep) const {
    std::vector<std::string> unraveled_data;

    std::string categories {};
    for (const auto &c: cat) {
        categories += c + sep;
    }
    unraveled_data.push_back(categories);

    for (const auto &d: data) {
        std::string line_data {};
        for (const auto &c: d) {
            line_data += c + sep;
        }
        unraveled_data.push_back(line_data);
    }

    write_several(unraveled_data);
}

void File::write_comments(const std::vector<std::string> &to_write) const {
    /// Prepare the comments
    for (const auto &comm: to_write) {
        /// Pass the "true" comment line to the method
        write_one_line(c_comm + comm);
    }
}

void File::write_sequence(const std::string &name, const std::string &sequence) const
{
    /// Create the vector element
    std::vector<std::string> to_be_written {};
    /// Push back the sequence's name
    to_be_written.push_back(c_name + name);
    /// Accumulate `c_per_line` characters and push_back them until `sequence` is exhausted
    std::size_t position {0};
    while (position < sequence.size()) {
        to_be_written.push_back(sequence.substr(position, c_per_line));
        position += c_per_line;
    }
    /// Write to files
    write_several(to_be_written);
}

void File::clear() const
{
    /// Open the files in truncation mode
    std::ofstream file {file_name, std::ios::trunc};
    /// Close the files
    file.close();
}