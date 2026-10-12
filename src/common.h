#ifndef TSINFERR_COMMON_H
#define TSINFERR_COMMON_H

#include <Rcpp.h>

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef restrict
#define restrict __restrict__
#endif

extern "C" {
#include "tsinfer.h"
}

namespace tsir {

[[noreturn]] void stop_tsi(int err, const char *context = nullptr);
[[noreturn]] void stop_tsk(int err, const char *context = nullptr);

inline void check_tsi(int err, const char *context = nullptr)
{
    if (err != 0) {
        stop_tsi(err, context);
    }
}

inline void check_tsk(int err, const char *context = nullptr)
{
    if (err < 0) {
        stop_tsk(err, context);
    }
}

/* Table collection held by an R external pointer. */
tsk_table_collection_t *get_tables(SEXP xp);
SEXP wrap_tables(tsk_table_collection_t *tables);
SEXP new_tables_xptr(tsk_table_collection_t **out);

/* Python-compatible canonical JSON helpers (sorted keys, compact, ASCII). */
std::string json_quote(const std::string &s);
std::string json_double(double x);

int resolve_num_threads(int num_threads);

} // namespace tsir

#endif
