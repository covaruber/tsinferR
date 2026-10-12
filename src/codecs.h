#ifndef TSINFERR_CODECS_H
#define TSINFERR_CODECS_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace tsir {

/* A bytes-to-bytes compression codec (Zarr v2 "compressor" or a v3
 * bytes->bytes codec). All functions are thread-safe and never call R. */
struct BytesCodec {
    std::string id; /* blosc, zstd, zlib, gzip, lz4, crc32c */
    std::string cname = "lz4";
    int clevel = 5;
    int shuffle = 1; /* 0 none, 1 byte, 2 bit, -1 auto (numcodecs) */
    int blocksize = 0;
    int typesize = 0; /* 0 = use the array item size */
    int level = 1;
};

/* Throws std::runtime_error on failure. expected_size is a hint (0 = unknown). */
std::vector<uint8_t> codec_decode(
    const BytesCodec &codec, const std::vector<uint8_t> &in, size_t expected_size);
std::vector<uint8_t> codec_encode(
    const BytesCodec &codec, const uint8_t *data, size_t size, size_t itemsize);

bool blosc_has_compressor(const std::string &cname);

uint32_t crc32c(const uint8_t *data, size_t n);

} // namespace tsir

#endif
