#include "common.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>

extern "C" void
tskr_bug_assert_failed(const char *file, int line)
{
    Rf_error("tskit internal bug detected in %s at line %d. Please report this "
             "at https://github.com/tskit-dev/tskit/issues",
        file, line);
}

namespace tsir {

void
stop_tsi(int err, const char *context)
{
    if (context != nullptr) {
        Rcpp::stop("%s: %s", context, tsi_strerror(err));
    }
    Rcpp::stop(tsi_strerror(err));
}

void
stop_tsk(int err, const char *context)
{
    if (context != nullptr) {
        Rcpp::stop("%s: tskit error %d: %s", context, err, tsk_strerror(err));
    }
    Rcpp::stop("tskit error %d: %s", err, tsk_strerror(err));
}

static void
tables_finalizer(SEXP xp)
{
    auto *tables = static_cast<tsk_table_collection_t *>(R_ExternalPtrAddr(xp));
    if (tables != nullptr) {
        tsk_table_collection_free(tables);
        std::free(tables);
        R_ClearExternalPtr(xp);
    }
}

SEXP
new_tables_xptr(tsk_table_collection_t **out)
{
    auto *tables
        = static_cast<tsk_table_collection_t *>(std::calloc(1, sizeof(tsk_table_collection_t)));
    if (tables == nullptr) {
        Rcpp::stop("Out of memory");
    }
    SEXP xp = PROTECT(
        R_MakeExternalPtr(tables, Rf_install("tsinferR_tables"), R_NilValue));
    R_RegisterCFinalizerEx(xp, tables_finalizer, TRUE);
    UNPROTECT(1);
    *out = tables;
    return xp;
}

SEXP
wrap_tables(tsk_table_collection_t *tables)
{
    SEXP xp = PROTECT(
        R_MakeExternalPtr(tables, Rf_install("tsinferR_tables"), R_NilValue));
    R_RegisterCFinalizerEx(xp, tables_finalizer, TRUE);
    UNPROTECT(1);
    return xp;
}

tsk_table_collection_t *
get_tables(SEXP xp)
{
    if (TYPEOF(xp) != EXTPTRSXP || R_ExternalPtrTag(xp) != Rf_install("tsinferR_tables")) {
        Rcpp::stop("Expected a tsinferR table collection pointer");
    }
    auto *tables = static_cast<tsk_table_collection_t *>(R_ExternalPtrAddr(xp));
    if (tables == nullptr) {
        Rcpp::stop("Tree sequence pointer is NULL (object was freed or restored "
                   "from a saved R session; use ts_dump()/ts_load() to persist "
                   "tree sequences)");
    }
    return tables;
}

std::string
json_quote(const std::string &s)
{
    std::string out = "\"";
    size_t i = 0;
    const size_t n = s.size();
    char buf[16];
    while (i < n) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        if (c < 0x80) {
            switch (c) {
                case '"':
                    out += "\\\"";
                    break;
                case '\\':
                    out += "\\\\";
                    break;
                case '\n':
                    out += "\\n";
                    break;
                case '\r':
                    out += "\\r";
                    break;
                case '\t':
                    out += "\\t";
                    break;
                case '\b':
                    out += "\\b";
                    break;
                case '\f':
                    out += "\\f";
                    break;
                default:
                    if (c < 0x20) {
                        std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                        out += buf;
                    } else {
                        out += static_cast<char>(c);
                    }
            }
            i++;
            continue;
        }
        /* Decode one UTF-8 code point and emit \uXXXX escapes like Python's
         * json.dumps(ensure_ascii=True). */
        uint32_t cp;
        size_t len;
        if ((c & 0xE0) == 0xC0) {
            cp = c & 0x1F;
            len = 2;
        } else if ((c & 0xF0) == 0xE0) {
            cp = c & 0x0F;
            len = 3;
        } else if ((c & 0xF8) == 0xF0) {
            cp = c & 0x07;
            len = 4;
        } else {
            cp = 0xFFFD;
            len = 1;
        }
        if (i + len > n) {
            cp = 0xFFFD;
            len = n - i;
        } else {
            for (size_t k = 1; k < len; k++) {
                cp = (cp << 6) | (static_cast<unsigned char>(s[i + k]) & 0x3F);
            }
        }
        if (cp >= 0x10000) {
            uint32_t v = cp - 0x10000;
            std::snprintf(buf, sizeof(buf), "\\u%04x", 0xD800 + (v >> 10));
            out += buf;
            std::snprintf(buf, sizeof(buf), "\\u%04x", 0xDC00 + (v & 0x3FF));
            out += buf;
        } else {
            std::snprintf(buf, sizeof(buf), "\\u%04x", cp);
            out += buf;
        }
        i += len;
    }
    out += "\"";
    return out;
}

std::string
json_double(double x)
{
    if (std::isnan(x)) {
        return "NaN";
    }
    if (std::isinf(x)) {
        return x > 0 ? "Infinity" : "-Infinity";
    }
    if (x == 0) {
        return std::signbit(x) ? "-0.0" : "0.0";
    }
    /* Shortest round-trip representation, formatted like Python's repr(). */
    char buf[64];
    int precision = 1;
    for (; precision <= 17; precision++) {
        std::snprintf(buf, sizeof(buf), "%.*e", precision - 1, x);
        if (std::strtod(buf, nullptr) == x) {
            break;
        }
    }
    std::string sci(buf);
    size_t epos = sci.find('e');
    std::string mantissa = sci.substr(0, epos);
    int exponent = std::atoi(sci.c_str() + epos + 1);
    bool negative = mantissa[0] == '-';
    if (negative) {
        mantissa = mantissa.substr(1);
    }
    std::string digits;
    for (char ch : mantissa) {
        if (ch != '.') {
            digits += ch;
        }
    }
    while (digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
    }
    std::string out = negative ? "-" : "";
    if (exponent >= -4 && exponent < 16) {
        int ndigits = static_cast<int>(digits.size());
        if (exponent >= 0) {
            if (ndigits <= exponent + 1) {
                out += digits + std::string(exponent + 1 - ndigits, '0') + ".0";
            } else {
                out += digits.substr(0, exponent + 1) + "." + digits.substr(exponent + 1);
            }
        } else {
            out += "0." + std::string(-exponent - 1, '0') + digits;
        }
    } else {
        out += digits.substr(0, 1);
        if (digits.size() > 1) {
            out += "." + digits.substr(1);
        }
        std::snprintf(buf, sizeof(buf), "e%c%02d", exponent < 0 ? '-' : '+',
            exponent < 0 ? -exponent : exponent);
        out += buf;
    }
    return out;
}

int
resolve_num_threads(int num_threads)
{
    if (num_threads <= 0) {
        return 1;
    }
    return num_threads;
}

} // namespace tsir
