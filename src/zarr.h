#ifndef TSINFERR_ZARR_H
#define TSINFERR_ZARR_H

#include "codecs.h"

#include <cstdint>
#include <string>
#include <vector>

namespace tsir {

enum class DType {
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Float32,
    Float64,
    Bool,
    VlenUtf8,
    FixedUnicode,
    FixedBytes
};

/* A Zarr (v2 or v3) array stored in a local directory. Metadata is parsed in R
 * and passed in; chunk I/O and codecs are implemented here and are safe to use
 * from worker threads (no R API calls). */
struct ZarrArray {
    std::string path;
    int zarr_format = 2;
    std::vector<int64_t> shape;
    std::vector<int64_t> chunks;
    DType dtype = DType::Int8;
    size_t itemsize = 1;
    bool big_endian = false;
    std::string serializer = "raw"; /* raw, vlen-utf8, packbits */
    std::vector<BytesCodec> codecs; /* in encoding order */
    std::string key_prefix;         /* "" (v2) or "c" (v3 default encoding) */
    char separator = '.';
    double fill_value = 0;
    std::string fill_string;

    bool is_string() const;
    size_t ndim() const { return shape.size(); }
    size_t chunk_elements() const;
    int64_t grid_size(size_t dim) const;
    std::string chunk_key(const std::vector<int64_t> &index) const;
    std::string chunk_path(const std::vector<int64_t> &index) const;
};

/* Numeric/bool/fixed-width chunks: raw little-endian element bytes, always the
 * full chunk shape. Returns false (and fills with the fill value) if the chunk
 * file does not exist. */
bool read_chunk(const ZarrArray &a, const std::vector<int64_t> &index, std::vector<uint8_t> &out);
void fill_chunk(const ZarrArray &a, std::vector<uint8_t> &out);
void write_chunk(const ZarrArray &a, const std::vector<int64_t> &index, const uint8_t *data);

bool read_chunk_strings(
    const ZarrArray &a, const std::vector<int64_t> &index, std::vector<std::string> &out);
void write_chunk_strings(
    const ZarrArray &a, const std::vector<int64_t> &index, const std::vector<std::string> &data);

/* Decode a chunk of integer genotypes (any integer dtype) into int32 values. */
bool read_chunk_int32(const ZarrArray &a, const std::vector<int64_t> &index, std::vector<int32_t> &out);

double element_as_double(const ZarrArray &a, const uint8_t *p);
int64_t element_as_int64(const ZarrArray &a, const uint8_t *p);
void set_element_from_double(const ZarrArray &a, uint8_t *p, double value);

void make_dirs(const std::string &path);
bool file_exists(const std::string &path);
std::vector<uint8_t> read_file(const std::string &path);
void write_file(const std::string &path, const uint8_t *data, size_t size);

} // namespace tsir

#endif
