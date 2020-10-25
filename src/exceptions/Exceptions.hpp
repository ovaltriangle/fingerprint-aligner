/** 
 * MainExceptions
 * @author astaroth
 * @date 10/11/20
 *
 * @brief
 */

#ifndef FA_EXCEPTIONS_HPP
#define FA_EXCEPTIONS_HPP

#include <exception>

/// main

struct RQSameFile : public std::exception
{
    RQSameFile() noexcept = default;
    ~RQSameFile() override = default;
    [[nodiscard]] const char *what() const noexcept override {
        return "The Reference and Query file need to be two distinct files, not the same one.";
    }
};

#endif //FA_EXCEPTIONS_HPP
