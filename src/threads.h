#ifndef TSINFERR_THREADS_H
#define TSINFERR_THREADS_H

#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace tsir {

/* Run f(i) for i in [0, n) on up to num_threads threads. f must not call the
 * R API. The first exception thrown by any worker is rethrown here. */
template <typename F>
void
parallel_for(size_t n, int num_threads, F f)
{
    if (n == 0) {
        return;
    }
    size_t workers = num_threads <= 1 ? 1 : static_cast<size_t>(num_threads);
    if (workers > n) {
        workers = n;
    }
    if (workers == 1) {
        for (size_t i = 0; i < n; i++) {
            f(i);
        }
        return;
    }
    std::atomic<size_t> next(0);
    std::exception_ptr error = nullptr;
    std::mutex error_mutex;
    auto worker = [&]() {
        while (true) {
            size_t i = next.fetch_add(1);
            if (i >= n) {
                break;
            }
            try {
                f(i);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (error == nullptr) {
                    error = std::current_exception();
                }
                next.store(n);
                break;
            }
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers);
    for (size_t t = 0; t < workers; t++) {
        threads.emplace_back(worker);
    }
    for (std::thread &t : threads) {
        t.join();
    }
    if (error != nullptr) {
        std::rethrow_exception(error);
    }
}

} // namespace tsir

#endif
