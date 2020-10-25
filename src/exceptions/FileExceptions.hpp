/** 
 * FileExceptions
 * @author astaroth
 * @date 10/8/20
 *
 * @brief
 */

#ifndef FA_FILEEXCEPTIONS_HPP
#define FA_FILEEXCEPTIONS_HPP

#include <exception>

/// File

struct FileNotFound : public std::exception
{
    FileNotFound() noexcept = default;
    ~FileNotFound() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "One of the supplied files was not found.";
    }
};

struct NoSequencesInFile : public std::exception
{
    NoSequencesInFile() noexcept = default;
    ~NoSequencesInFile() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "No sequences were found in one of the supplied files.";
    }
};

/// FileManager

struct FileWithNoFormat : public std::exception
{
    FileWithNoFormat() noexcept = default;
    ~FileWithNoFormat() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "One of the supplied files does not have a valid extension.";
    }
};

struct UnhandledFormat : public std::exception
{
    UnhandledFormat() noexcept = default;
    ~UnhandledFormat() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "The format of one of the supplied files is not yet supported.";
    }
};

#endif //FA_FILEEXCEPTIONS_HPP
