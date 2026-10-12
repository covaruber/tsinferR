#include "codecs.h"

#include <array>
#include <cstring>
#include <stdexcept>

extern "C" {
#include "blosc.h"
#include "lz4.h"
#include "zstd.h"
}
#include <zlib.h>

namespace tsir {

static std::runtime_error
codec_error(const std::string &id, const std::string &msg)
{
    return std::runtime_error("codec '" + id + "': " + msg);
}

bool
blosc_has_compressor(const std::string &cname)
{
    return blosc_compname_to_compcode(cname.c_str()) >= 0;
}

static std::vector<uint8_t>
blosc_decode(const BytesCodec &codec, const std::vector<uint8_t> &in)
{
    if (in.size() < BLOSC_MIN_HEADER_LENGTH) {
        throw codec_error(codec.id, "truncated blosc buffer");
    }
    size_t nbytes = 0, cbytes = 0, blocksize = 0;
    blosc_cbuffer_sizes(in.data(), &nbytes, &cbytes, &blocksize);
    if (cbytes > in.size()) {
        throw codec_error(codec.id, "corrupt blosc header");
    }
    std::vector<uint8_t> out(nbytes);
    if (nbytes == 0) {
        return out;
    }
    int ret = blosc_decompress_ctx(in.data(), out.data(), nbytes, 1);
    if (ret < 0 || static_cast<size_t>(ret) != nbytes) {
        throw codec_error(codec.id, "decompression failed (code " + std::to_string(ret) + ")");
    }
    return out;
}

static std::vector<uint8_t>
blosc_encode(const BytesCodec &codec, const uint8_t *data, size_t size, size_t itemsize)
{
    size_t typesize = codec.typesize > 0 ? static_cast<size_t>(codec.typesize) : itemsize;
    if (typesize == 0 || typesize > BLOSC_MAX_TYPESIZE) {
        typesize = 1;
    }
    int shuffle = codec.shuffle;
    if (shuffle < 0) {
        shuffle = typesize == 1 ? BLOSC_BITSHUFFLE : BLOSC_SHUFFLE;
    }
    if (blosc_compname_to_compcode(codec.cname.c_str()) < 0) {
        throw codec_error(codec.id, "unsupported blosc compressor '" + codec.cname + "'");
    }
    std::vector<uint8_t> out(size + BLOSC_MAX_OVERHEAD);
    int ret = blosc_compress_ctx(codec.clevel, shuffle, typesize, size, data, out.data(),
        out.size(), codec.cname.c_str(), static_cast<size_t>(codec.blocksize), 1);
    if (ret <= 0) {
        throw codec_error(codec.id, "compression failed (code " + std::to_string(ret) + ")");
    }
    out.resize(static_cast<size_t>(ret));
    return out;
}

static std::vector<uint8_t>
zstd_decode(const BytesCodec &codec, const std::vector<uint8_t> &in, size_t expected)
{
    unsigned long long size = ZSTD_getFrameContentSize(in.data(), in.size());
    if (size == ZSTD_CONTENTSIZE_ERROR) {
        throw codec_error(codec.id, "not a zstd frame");
    }
    if (size != ZSTD_CONTENTSIZE_UNKNOWN) {
        std::vector<uint8_t> out(static_cast<size_t>(size));
        size_t ret = ZSTD_decompress(out.data(), out.size(), in.data(), in.size());
        if (ZSTD_isError(ret)) {
            throw codec_error(codec.id, ZSTD_getErrorName(ret));
        }
        out.resize(ret);
        return out;
    }
    ZSTD_DStream *stream = ZSTD_createDStream();
    if (stream == nullptr) {
        throw codec_error(codec.id, "out of memory");
    }
    ZSTD_initDStream(stream);
    std::vector<uint8_t> out;
    std::vector<uint8_t> buf(expected > 0 ? expected : ZSTD_DStreamOutSize());
    ZSTD_inBuffer input = { in.data(), in.size(), 0 };
    while (input.pos < input.size) {
        ZSTD_outBuffer output = { buf.data(), buf.size(), 0 };
        size_t ret = ZSTD_decompressStream(stream, &output, &input);
        if (ZSTD_isError(ret)) {
            ZSTD_freeDStream(stream);
            throw codec_error(codec.id, ZSTD_getErrorName(ret));
        }
        out.insert(out.end(), buf.data(), buf.data() + output.pos);
        if (ret == 0 && input.pos == input.size) {
            break;
        }
    }
    ZSTD_freeDStream(stream);
    return out;
}

static std::vector<uint8_t>
zstd_encode(const BytesCodec &codec, const uint8_t *data, size_t size)
{
    std::vector<uint8_t> out(ZSTD_compressBound(size));
    size_t ret = ZSTD_compress(out.data(), out.size(), data, size, codec.level);
    if (ZSTD_isError(ret)) {
        throw codec_error(codec.id, ZSTD_getErrorName(ret));
    }
    out.resize(ret);
    return out;
}

static std::vector<uint8_t>
zlib_decode(const BytesCodec &codec, const std::vector<uint8_t> &in, size_t expected)
{
    z_stream strm;
    std::memset(&strm, 0, sizeof(strm));
    /* 15 + 32: automatically detect zlib or gzip headers. */
    if (inflateInit2(&strm, 15 + 32) != Z_OK) {
        throw codec_error(codec.id, "inflateInit failed");
    }
    std::vector<uint8_t> out;
    std::vector<uint8_t> buf(expected > 0 ? expected : 65536);
    strm.next_in = const_cast<Bytef *>(in.data());
    strm.avail_in = static_cast<uInt>(in.size());
    int ret = Z_OK;
    while (ret != Z_STREAM_END) {
        strm.next_out = buf.data();
        strm.avail_out = static_cast<uInt>(buf.size());
        ret = inflate(&strm, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END) {
            inflateEnd(&strm);
            throw codec_error(codec.id, "inflate failed (code " + std::to_string(ret) + ")");
        }
        out.insert(out.end(), buf.data(), buf.data() + (buf.size() - strm.avail_out));
        if (ret == Z_OK && strm.avail_in == 0 && strm.avail_out != 0) {
            break;
        }
    }
    inflateEnd(&strm);
    return out;
}

static std::vector<uint8_t>
zlib_encode(const BytesCodec &codec, const uint8_t *data, size_t size, bool gzip)
{
    z_stream strm;
    std::memset(&strm, 0, sizeof(strm));
    int window_bits = gzip ? 15 + 16 : 15;
    if (deflateInit2(&strm, codec.level, Z_DEFLATED, window_bits, 8, Z_DEFAULT_STRATEGY)
        != Z_OK) {
        throw codec_error(codec.id, "deflateInit failed");
    }
    std::vector<uint8_t> out(deflateBound(&strm, static_cast<uLong>(size)) + 32);
    strm.next_in = const_cast<Bytef *>(data);
    strm.avail_in = static_cast<uInt>(size);
    strm.next_out = out.data();
    strm.avail_out = static_cast<uInt>(out.size());
    int ret = deflate(&strm, Z_FINISH);
    if (ret != Z_STREAM_END) {
        deflateEnd(&strm);
        throw codec_error(codec.id, "deflate failed");
    }
    out.resize(strm.total_out);
    deflateEnd(&strm);
    return out;
}

static std::vector<uint8_t>
lz4_decode(const BytesCodec &codec, const std::vector<uint8_t> &in)
{
    if (in.size() < 4) {
        throw codec_error(codec.id, "truncated buffer");
    }
    uint32_t size = static_cast<uint32_t>(in[0]) | (static_cast<uint32_t>(in[1]) << 8)
                    | (static_cast<uint32_t>(in[2]) << 16)
                    | (static_cast<uint32_t>(in[3]) << 24);
    std::vector<uint8_t> out(size);
    if (size == 0) {
        return out;
    }
    int ret = LZ4_decompress_safe(reinterpret_cast<const char *>(in.data() + 4),
        reinterpret_cast<char *>(out.data()), static_cast<int>(in.size() - 4),
        static_cast<int>(size));
    if (ret < 0 || static_cast<uint32_t>(ret) != size) {
        throw codec_error(codec.id, "decompression failed");
    }
    return out;
}

static std::vector<uint8_t>
lz4_encode(const BytesCodec &codec, const uint8_t *data, size_t size)
{
    int bound = LZ4_compressBound(static_cast<int>(size));
    std::vector<uint8_t> out(static_cast<size_t>(bound) + 4);
    out[0] = static_cast<uint8_t>(size & 0xFF);
    out[1] = static_cast<uint8_t>((size >> 8) & 0xFF);
    out[2] = static_cast<uint8_t>((size >> 16) & 0xFF);
    out[3] = static_cast<uint8_t>((size >> 24) & 0xFF);
    int ret = LZ4_compress_fast(reinterpret_cast<const char *>(data),
        reinterpret_cast<char *>(out.data() + 4), static_cast<int>(size), bound,
        codec.level > 0 ? codec.level : 1);
    if (ret <= 0 && size > 0) {
        throw codec_error(codec.id, "compression failed");
    }
    out.resize(static_cast<size_t>(ret) + 4);
    return out;
}

uint32_t
crc32c(const uint8_t *data, size_t n)
{
    static const std::array<uint32_t, 256> table = [] {
        std::array<uint32_t, 256> t{};
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t c = i;
            for (int k = 0; k < 8; k++) {
                c = (c & 1) ? (c >> 1) ^ 0x82F63B78u : c >> 1;
            }
            t[i] = c;
        }
        return t;
    }();
    uint32_t crc = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; i++) {
        crc = table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFFu;
}

std::vector<uint8_t>
codec_decode(const BytesCodec &codec, const std::vector<uint8_t> &in, size_t expected_size)
{
    if (codec.id == "blosc") {
        return blosc_decode(codec, in);
    }
    if (codec.id == "zstd") {
        return zstd_decode(codec, in, expected_size);
    }
    if (codec.id == "zlib" || codec.id == "gzip") {
        return zlib_decode(codec, in, expected_size);
    }
    if (codec.id == "lz4") {
        return lz4_decode(codec, in);
    }
    if (codec.id == "crc32c") {
        if (in.size() < 4) {
            throw codec_error(codec.id, "truncated buffer");
        }
        size_t n = in.size() - 4;
        uint32_t stored = static_cast<uint32_t>(in[n]) | (static_cast<uint32_t>(in[n + 1]) << 8)
                          | (static_cast<uint32_t>(in[n + 2]) << 16)
                          | (static_cast<uint32_t>(in[n + 3]) << 24);
        if (crc32c(in.data(), n) != stored) {
            throw codec_error(codec.id, "checksum mismatch");
        }
        return std::vector<uint8_t>(in.begin(), in.begin() + static_cast<long>(n));
    }
    throw codec_error(codec.id, "unsupported codec");
}

std::vector<uint8_t>
codec_encode(const BytesCodec &codec, const uint8_t *data, size_t size, size_t itemsize)
{
    if (codec.id == "blosc") {
        return blosc_encode(codec, data, size, itemsize);
    }
    if (codec.id == "zstd") {
        return zstd_encode(codec, data, size);
    }
    if (codec.id == "zlib") {
        return zlib_encode(codec, data, size, false);
    }
    if (codec.id == "gzip") {
        return zlib_encode(codec, data, size, true);
    }
    if (codec.id == "lz4") {
        return lz4_encode(codec, data, size);
    }
    if (codec.id == "crc32c") {
        std::vector<uint8_t> out(data, data + size);
        uint32_t c = crc32c(data, size);
        for (int k = 0; k < 4; k++) {
            out.push_back(static_cast<uint8_t>((c >> (8 * k)) & 0xFF));
        }
        return out;
    }
    throw codec_error(codec.id, "unsupported codec");
}

} // namespace tsir
