#include "common.h"
#include "zarr.h"
#include "zarr_r.h"

using namespace Rcpp;

namespace tsir {

static DType
parse_dtype(const std::string &name)
{
    static const std::pair<const char *, DType> table[] = { { "int8", DType::Int8 },
        { "int16", DType::Int16 }, { "int32", DType::Int32 }, { "int64", DType::Int64 },
        { "uint8", DType::UInt8 }, { "uint16", DType::UInt16 }, { "uint32", DType::UInt32 },
        { "uint64", DType::UInt64 }, { "float32", DType::Float32 },
        { "float64", DType::Float64 }, { "bool", DType::Bool }, { "string", DType::VlenUtf8 },
        { "fixed_unicode", DType::FixedUnicode }, { "fixed_bytes", DType::FixedBytes } };
    for (const auto &entry : table) {
        if (name == entry.first) {
            return entry.second;
        }
    }
    stop("Unsupported Zarr data type: %s", name);
}

ZarrArray
zarr_array_from_list(List meta)
{
    ZarrArray a;
    a.path = as<std::string>(meta["path"]);
    a.zarr_format = as<int>(meta["zarr_format"]);
    NumericVector shape = meta["shape"];
    NumericVector chunks = meta["chunks"];
    if (shape.size() != chunks.size()) {
        stop("Zarr metadata: shape and chunks differ in length");
    }
    for (R_xlen_t j = 0; j < shape.size(); j++) {
        a.shape.push_back(static_cast<int64_t>(shape[j]));
        a.chunks.push_back(static_cast<int64_t>(chunks[j]));
        if (chunks[j] < 1) {
            stop("Zarr metadata: chunk sizes must be positive");
        }
    }
    a.dtype = parse_dtype(as<std::string>(meta["dtype"]));
    a.itemsize = static_cast<size_t>(as<double>(meta["itemsize"]));
    a.big_endian = as<bool>(meta["big_endian"]);
    a.serializer = as<std::string>(meta["serializer"]);
    a.key_prefix = as<std::string>(meta["key_prefix"]);
    std::string sep = as<std::string>(meta["separator"]);
    a.separator = sep.empty() ? '.' : sep[0];
    a.fill_value = as<double>(meta["fill_value"]);
    a.fill_string = as<std::string>(meta["fill_string"]);
    List codecs = meta["codecs"];
    for (R_xlen_t j = 0; j < codecs.size(); j++) {
        List c = codecs[j];
        BytesCodec bc;
        bc.id = as<std::string>(c["id"]);
        if (c.containsElementNamed("cname")) {
            bc.cname = as<std::string>(c["cname"]);
        }
        if (c.containsElementNamed("clevel")) {
            bc.clevel = as<int>(c["clevel"]);
        }
        if (c.containsElementNamed("shuffle")) {
            bc.shuffle = as<int>(c["shuffle"]);
        }
        if (c.containsElementNamed("blocksize")) {
            bc.blocksize = as<int>(c["blocksize"]);
        }
        if (c.containsElementNamed("typesize")) {
            bc.typesize = as<int>(c["typesize"]);
        }
        if (c.containsElementNamed("level")) {
            bc.level = as<int>(c["level"]);
        }
        a.codecs.push_back(bc);
    }
    return a;
}

ZarrArray *
get_zarr_array(SEXP xp)
{
    XPtr<ZarrArray> p(xp);
    return p.get();
}

} // namespace tsir

using tsir::DType;
using tsir::ZarrArray;

// [[Rcpp::export(.zarr_array)]]
SEXP
zarr_array(List meta)
{
    XPtr<ZarrArray> xp(new ZarrArray(tsir::zarr_array_from_list(meta)), true);
    return xp;
}

static SEXPTYPE
r_type_for(const ZarrArray &a)
{
    switch (a.dtype) {
        case DType::Int8:
        case DType::Int16:
        case DType::Int32:
        case DType::UInt8:
        case DType::UInt16:
            return INTSXP;
        case DType::Bool:
            return LGLSXP;
        case DType::VlenUtf8:
        case DType::FixedUnicode:
        case DType::FixedBytes:
            return STRSXP;
        default:
            return REALSXP;
    }
}

/* Calls f(chunk_index, region_lo, region_hi) for every chunk intersecting the
 * region [start, start + count). */
template <typename F>
static void
for_each_chunk(const ZarrArray &a, const std::vector<int64_t> &start,
    const std::vector<int64_t> &count, F f)
{
    size_t nd = a.ndim();
    std::vector<int64_t> first(nd), last(nd), idx(nd);
    for (size_t d = 0; d < nd; d++) {
        if (count[d] == 0) {
            return;
        }
        first[d] = start[d] / a.chunks[d];
        last[d] = (start[d] + count[d] - 1) / a.chunks[d];
    }
    idx = first;
    while (true) {
        f(idx);
        size_t d = nd;
        while (d > 0) {
            d--;
            if (idx[d] < last[d]) {
                idx[d]++;
                for (size_t e = d + 1; e < nd; e++) {
                    idx[e] = first[e];
                }
                break;
            }
            if (d == 0) {
                return;
            }
        }
        if (nd == 0) {
            return;
        }
    }
}

/* Iterate the elements in the intersection of a chunk and the region, giving
 * the chunk-local C-order offset and the R (column-major) offset in the region. */
template <typename F>
static void
for_each_element(const ZarrArray &a, const std::vector<int64_t> &chunk_idx,
    const std::vector<int64_t> &start, const std::vector<int64_t> &count, F f)
{
    size_t nd = a.ndim();
    if (nd == 0) {
        f(0, 0);
        return;
    }
    std::vector<int64_t> lo(nd), hi(nd), pos(nd);
    for (size_t d = 0; d < nd; d++) {
        int64_t c0 = chunk_idx[d] * a.chunks[d];
        lo[d] = std::max(start[d], c0);
        hi[d] = std::min(start[d] + count[d], c0 + a.chunks[d]);
        if (lo[d] >= hi[d]) {
            return;
        }
    }
    std::vector<int64_t> cstride(nd), rstride(nd);
    cstride[nd - 1] = 1;
    for (size_t d = nd - 1; d > 0; d--) {
        cstride[d - 1] = cstride[d] * a.chunks[d];
    }
    rstride[0] = 1;
    for (size_t d = 1; d < nd; d++) {
        rstride[d] = rstride[d - 1] * count[d - 1];
    }
    pos = lo;
    while (true) {
        int64_t coff = 0, roff = 0;
        for (size_t d = 0; d < nd; d++) {
            coff += (pos[d] - chunk_idx[d] * a.chunks[d]) * cstride[d];
            roff += (pos[d] - start[d]) * rstride[d];
        }
        f(static_cast<size_t>(coff), static_cast<R_xlen_t>(roff));
        size_t d = nd;
        bool done = true;
        while (d > 0) {
            d--;
            if (pos[d] + 1 < hi[d]) {
                pos[d]++;
                for (size_t e = d + 1; e < nd; e++) {
                    pos[e] = lo[e];
                }
                done = false;
                break;
            }
        }
        if (done) {
            return;
        }
    }
}

static void
check_region(const ZarrArray &a, const std::vector<int64_t> &start,
    const std::vector<int64_t> &count)
{
    if (start.size() != a.ndim() || count.size() != a.ndim()) {
        stop("region dimension mismatch");
    }
    for (size_t d = 0; d < a.ndim(); d++) {
        if (start[d] < 0 || count[d] < 0 || start[d] + count[d] > a.shape[d]) {
            stop("region out of bounds in dimension %d", static_cast<int>(d + 1));
        }
    }
}

// [[Rcpp::export(.zarr_read)]]
SEXP
zarr_read(SEXP xp, NumericVector start_, NumericVector count_)
{
    ZarrArray &a = *tsir::get_zarr_array(xp);
    std::vector<int64_t> start(start_.begin(), start_.end());
    std::vector<int64_t> count(count_.begin(), count_.end());
    check_region(a, start, count);
    R_xlen_t total = 1;
    for (int64_t c : count) {
        total *= static_cast<R_xlen_t>(c);
    }
    SEXPTYPE type = r_type_for(a);
    SEXP out = PROTECT(Rf_allocVector(type, total));
    try {
        if (type == STRSXP) {
            std::vector<std::string> buf;
            for_each_chunk(a, start, count, [&](const std::vector<int64_t> &ci) {
                tsir::read_chunk_strings(a, ci, buf);
                for_each_element(a, ci, start, count, [&](size_t coff, R_xlen_t roff) {
                    const std::string &s = buf[coff];
                    SET_STRING_ELT(out, roff,
                        Rf_mkCharLenCE(s.data(), static_cast<int>(s.size()), CE_UTF8));
                });
            });
        } else {
            std::vector<uint8_t> buf;
            for_each_chunk(a, start, count, [&](const std::vector<int64_t> &ci) {
                tsir::read_chunk(a, ci, buf);
                for_each_element(a, ci, start, count, [&](size_t coff, R_xlen_t roff) {
                    const uint8_t *p = buf.data() + coff * a.itemsize;
                    double v = tsir::element_as_double(a, p);
                    if (type == INTSXP) {
                        INTEGER(out)[roff] = static_cast<int>(v);
                    } else if (type == LGLSXP) {
                        LOGICAL(out)[roff] = v != 0;
                    } else {
                        REAL(out)[roff] = v;
                    }
                });
            });
        }
    } catch (std::exception &e) {
        UNPROTECT(1);
        stop(e.what());
    }
    if (count.size() > 1) {
        IntegerVector dims(count.begin(), count.end());
        Rf_setAttrib(out, R_DimSymbol, dims);
    }
    UNPROTECT(1);
    return out;
}

// [[Rcpp::export(.zarr_write)]]
void
zarr_write(SEXP xp, SEXP data, NumericVector start_, NumericVector count_)
{
    ZarrArray &a = *tsir::get_zarr_array(xp);
    std::vector<int64_t> start(start_.begin(), start_.end());
    std::vector<int64_t> count(count_.begin(), count_.end());
    check_region(a, start, count);
    R_xlen_t total = 1;
    for (int64_t c : count) {
        total *= static_cast<R_xlen_t>(c);
    }
    if (Rf_xlength(data) != total) {
        stop("data has %d elements but the region has %d",
            static_cast<double>(Rf_xlength(data)), static_cast<double>(total));
    }
    bool is_str = a.is_string();
    if (is_str && TYPEOF(data) != STRSXP) {
        stop("string arrays require character data");
    }
    if (!is_str && TYPEOF(data) != INTSXP && TYPEOF(data) != REALSXP
        && TYPEOF(data) != LGLSXP) {
        stop("numeric arrays require numeric or logical data");
    }
    try {
        for_each_chunk(a, start, count, [&](const std::vector<int64_t> &ci) {
            bool full = true;
            for (size_t d = 0; d < a.ndim(); d++) {
                int64_t c0 = ci[d] * a.chunks[d];
                int64_t c1 = std::min(c0 + a.chunks[d], a.shape[d]);
                if (start[d] > c0 || start[d] + count[d] < c1) {
                    full = false;
                }
            }
            if (is_str) {
                std::vector<std::string> buf;
                if (full) {
                    buf.assign(a.chunk_elements(), a.fill_string);
                } else {
                    tsir::read_chunk_strings(a, ci, buf);
                }
                for_each_element(a, ci, start, count, [&](size_t coff, R_xlen_t roff) {
                    SEXP s = STRING_ELT(data, roff);
                    buf[coff] = s == NA_STRING ? std::string() : Rf_translateCharUTF8(s);
                });
                tsir::write_chunk_strings(a, ci, buf);
            } else {
                std::vector<uint8_t> buf;
                if (full) {
                    tsir::fill_chunk(a, buf);
                } else {
                    tsir::read_chunk(a, ci, buf);
                }
                for_each_element(a, ci, start, count, [&](size_t coff, R_xlen_t roff) {
                    double v;
                    if (TYPEOF(data) == REALSXP) {
                        v = REAL(data)[roff];
                    } else {
                        int iv = TYPEOF(data) == INTSXP ? INTEGER(data)[roff]
                                                        : LOGICAL(data)[roff];
                        v = iv == NA_INTEGER ? a.fill_value : iv;
                    }
                    if (std::isnan(v) && a.dtype != DType::Float32
                        && a.dtype != DType::Float64) {
                        v = a.fill_value;
                    }
                    tsir::set_element_from_double(a, buf.data() + coff * a.itemsize, v);
                });
                tsir::write_chunk(a, ci, buf.data());
            }
        });
    } catch (std::exception &e) {
        stop(e.what());
    }
}

// [[Rcpp::export(.blosc_has_compressor)]]
bool
blosc_has_compressor_r(std::string cname)
{
    return tsir::blosc_has_compressor(cname);
}
