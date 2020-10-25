/** 
 * threader
 * @author astaroth
 * @date 10/22/20
 *
 * @brief
 */

#ifndef FA_THREADER_HPP
#define FA_THREADER_HPP

#include <thread>
#include <functional>

// one day, maybe...
namespace utils {
    template <typename Container, class Class>
    void threaded_bulk(const Container &c, std::function<void (const std::size_t&, const std::size_t&)> &bulk,
                          const Class &that, const unsigned short &threads_number) {
        std::size_t amount {c.size()};
        if (threads_number > 1) {
            unsigned short thread_n {static_cast<unsigned short>((threads_number > amount) ? amount : threads_number)};
            std::thread threader[thread_n];
            std::size_t portion_per_thread {amount / thread_n}, processed {0};
            for (unsigned short i {0}; i < thread_n; ++i) {
                if (i == thread_n - 1)
                    portion_per_thread += (amount - processed - portion_per_thread);
                threader[i] = std::thread(&bulk, that, processed, portion_per_thread);
                processed += portion_per_thread;
            }
            for (unsigned short i {0}; i < thread_n; ++i)
                threader[i].join();
        } else {
            bulk(0, amount);
        }
    }
}


#endif //FA_THREADER_HPP
