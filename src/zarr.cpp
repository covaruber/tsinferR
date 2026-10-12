#include "zarr.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <sys/stat.h>
#if defined(_WIN32)
#include <direct.h>
#endif

namespace tsir {

bool
ZarrArray::is_string() const
{
    return dtype == DType::VlenUtf8 || dtype == DType::FixedUnicode
           || dtype == DType::FixedBytes;
}

size_t
ZarrArray::chunk_elements() const
{
    size_t n = 1;
    for (int64_t c : chunks) {
        n *= static_cast<size_t>(c);
    }
    return n;
}

int64_t
ZarrArray::grid_size(size_t dim) const
{
    if (chunks[dim] <= 0) {
        return 0;
    }
    return (shape[dim] + chunks[dim] - 1) / chunks[dim];
}

std::string
ZarrArray::chunk_key(const std::vector<int64_t> &index) const
{
    std::string key = key_prefix;
    if (index.empty()) {
        return key.empty() ? "0" : key;
    }
    for (size_t j = 0; j < index.size(); j++) {
        if (!key.empty()) {
            key += separator;
        }
        key += std::to_string(index[j]);
    }
    return key;
}

std::string
ZarrArray::chunk_path(const std::vector<int64_t> &index) const
{
    return path + "/" + chunk_key(index);
}

bool
file_exists(const std::string &path)
{
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

void
make_dirs(const std::string &path)
{
    std::string partial;
    for (size_t j = 0; j <= path.size(); j++) {
        if (j == path.size() || path[j] == '/' || path[j] == '\\') {
            if (!partial.empty() && !file_exists(partial)) {
#if defined(_WIN32)
                int ret = _mkdir(partial.c_str());
#else
                int ret = mkdir(partial.c_str(), 0777);
#endif
                if (ret != 0 && errno != EEXIST) {
                    throw std::runtime_error("cannot create directory '" + partial + "'");
                }
            }
        }
        if (j < path.size()) {
            partial += path[j];
        }
    }
}

std::vector<uint8_t>
read_file(const std::string &path)
{
    FILE *f = std::fopen(path.c_str(), "rb");
    if (f == nullptr) {
        throw std::runtime_error("cannot open '" + path + "'");
    }
    std::vector<uint8_t> out;
    uint8_t buf[65536];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) {
        out.insert(out.end(), buf, buf + n);
    }
    bool failed = std::ferror(f) != 0;
    std::fclose(f);
    if (failed) {
        throw std::runtime_error("error reading '" + path + "'");
    }
    return out;
}

void
write_file(const std::string &path, const uint8_t *data, size_t size)
{
    size_t slash = path.find_last_of("/\\");
    if (slash != std::string::npos) {
        make_dirs(path.substr(0, slash));
    }
    FILE *f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
        throw std::runtime_error("cannot open '" + path + "' for writing");
    }
    size_t written = size > 0 ? std::fwrite(data, 1, size, f) : 0;
    bool failed = written != size;
    if (std::fclose(f) != 0) {
        failed = true;
    }
    if (failed) {
        throw std::runtime_error("error writing '" + path + "'");
    }
}

static void
byteswap(uint8_t *data, size_t n, size_t itemsize)
{
    if (itemsize <= 1) {
        return;
    }
    for (size_t j = 0; j < n; j++) {
        uint8_t *p = data + j * itemsize;
        for (size_t a = 0, b = itemsize - 1; a < b; a++, b--) {
            std::swap(p[a], p[b]);
        }
    }
}

static std::vector<uint8_t>
decode_bytes(const ZarrArray &a, const std::string &path, size_t expected)
{
    std::vector<uint8_t> buf = read_file(path);
    for (size_t j = a.codecs.size(); j-- > 0;) {
        buf = codec_decode(a.codecs[j], buf, j == 0 ? expected : 0);
    }
    return buf;
}

static void
encode_and_write(const ZarrArray &a, const std::string &path, std::vector<uint8_t> buf)
{
    for (const BytesCodec &c : a.codecs) {
        buf = codec_encode(c, buf.data(), buf.size(), a.itemsize);
    }
    write_file(path, buf.data(), buf.size());
}

void
set_element_from_double(const ZarrArray &a, uint8_t *p, double v)
{
    switch (a.dtype) {
        case DType::Int8: {
            int8_t x = static_cast<int8_t>(v);
            std::memcpy(p, &x, 1);
            break;
        }
        case DType::Int16: {
            int16_t x = static_cast<int16_t>(v);
            std::memcpy(p, &x, 2);
            break;
        }
        case DType::Int32: {
            int32_t x = static_cast<int32_t>(v);
            std::memcpy(p, &x, 4);
            break;
        }
        case DType::Int64: {
            int64_t x = static_cast<int64_t>(v);
            std::memcpy(p, &x, 8);
            break;
        }
        case DType::UInt8:
        case DType::Bool: {
            uint8_t x = static_cast<uint8_t>(v);
            std::memcpy(p, &x, 1);
            break;
        }
        case DType::UInt16: {
            uint16_t x = static_cast<uint16_t>(v);
            std::memcpy(p, &x, 2);
            break;
        }
        case DType::UInt32: {
            uint32_t x = static_cast<uint32_t>(v);
            std::memcpy(p, &x, 4);
            break;
        }
        case DType::UInt64: {
            uint64_t x = static_cast<uint64_t>(v);
            std::memcpy(p, &x, 8);
            break;
        }
        case DType::Float32: {
            float x = static_cast<float>(v);
            std::memcpy(p, &x, 4);
            break;
        }
        case DType::Float64:
            std::memcpy(p, &v, 8);
            break;
        default:
            throw std::runtime_error("not a numeric dtype");
    }
}

double
element_as_double(const ZarrArray &a, const uint8_t *p)
{
    switch (a.dtype) {
        case DType::Int8: {
            int8_t x;
            std::memcpy(&x, p, 1);
            return x;
        }
        case DType::Int16: {
            int16_t x;
            std::memcpy(&x, p, 2);
            return x;
        }
        case DType::Int32: {
            int32_t x;
            std::memcpy(&x, p, 4);
            return x;
        }
        case DType::Int64: {
            int64_t x;
            std::memcpy(&x, p, 8);
            return static_cast<double>(x);
        }
        case DType::UInt8:
        case DType::Bool:
            return *p;
        case DType::UInt16: {
            uint16_t x;
            std::memcpy(&x, p, 2);
            return x;
        }
        case DType::UInt32: {
            uint32_t x;
            std::memcpy(&x, p, 4);
            return x;
        }
        case DType::UInt64: {
            uint64_t x;
            std::memcpy(&x, p, 8);
            return static_cast<double>(x);
        }
        case DType::Float32: {
            float x;
            std::memcpy(&x, p, 4);
            return x;
        }
        case DType::Float64: {
            double x;
            std::memcpy(&x, p, 8);
            return x;
        }
        default:
            throw std::runtime_error("not a numeric dtype");
    }
}

int64_t
element_as_int64(const ZarrArray &a, const uint8_t *p)
{
    if (a.dtype == DType::Int64) {
        int64_t x;
        std::memcpy(&x, p, 8);
        return x;
    }
    if (a.dtype == DType::UInt64) {
        uint64_t x;
        std::memcpy(&x, p, 8);
        return static_cast<int64_t>(x);
    }
    return static_cast<int64_t>(element_as_double(a, p));
}

void
fill_chunk(const ZarrArray &a, std::vector<uint8_t> &out)
{
    size_t n = a.chunk_elements();
    out.assign(n * a.itemsize, 0);
    if (a.is_string()) {
        return;
    }
    if (a.fill_value == 0 && !std::signbit(a.fill_value)) {
        return;
    }
    std::vector<uint8_t> one(a.itemsize);
    set_element_from_double(a, one.data(), a.fill_value);
    for (size_t j = 0; j < n; j++) {
        std::memcpy(out.data() + j * a.itemsize, one.data(), a.itemsize);
    }
}

static void
unpack_bits(const std::vector<uint8_t> &in, size_t n, uint8_t *out)
{
    if (in.empty()) {
        throw std::runtime_error("packbits: empty buffer");
    }
    for (size_t j = 0; j < n; j++) {
        size_t byte = 1 + j / 8;
        if (byte >= in.size()) {
            throw std::runtime_error("packbits: truncated buffer");
        }
        out[j] = (in[byte] >> (7 - (j % 8))) & 1;
    }
}

static std::vector<uint8_t>
pack_bits(const uint8_t *in, size_t n)
{
    size_t nbytes = (n + 7) / 8;
    std::vector<uint8_t> out(1 + nbytes, 0);
    out[0] = static_cast<uint8_t>((8 - n % 8) % 8);
    for (size_t j = 0; j < n; j++) {
        if (in[j]) {
            out[1 + j / 8] |= static_cast<uint8_t>(1u << (7 - (j % 8)));
        }
    }
    return out;
}

bool
read_chunk(const ZarrArray &a, const std::vector<int64_t> &index, std::vector<uint8_t> &out)
{
    if (a.dtype == DType::VlenUtf8) {
        throw std::runtime_error("read_chunk called on a variable-length string array");
    }
    std::string path = a.chunk_path(index);
    if (!file_exists(path)) {
        fill_chunk(a, out);
        return false;
    }
    size_t n = a.chunk_elements();
    size_t nbytes = n * a.itemsize;
    if (a.serializer == "packbits") {
        std::vector<uint8_t> buf = decode_bytes(a, path, 0);
        out.assign(n, 0);
        unpack_bits(buf, n, out.data());
        return true;
    }
    out = decode_bytes(a, path, nbytes);
    if (out.size() != nbytes) {
        throw std::runtime_error("chunk '" + path + "' has " + std::to_string(out.size())
                                 + " bytes, expected " + std::to_string(nbytes));
    }
    if (a.big_endian && a.dtype != DType::FixedUnicode) {
        byteswap(out.data(), n, a.itemsize);
    } else if (a.big_endian) {
        byteswap(out.data(), nbytes / 4, 4);
    }
    return true;
}

void
write_chunk(const ZarrArray &a, const std::vector<int64_t> &index, const uint8_t *data)
{
    if (a.dtype == DType::VlenUtf8) {
        throw std::runtime_error("write_chunk called on a variable-length string array");
    }
    size_t n = a.chunk_elements();
    std::vector<uint8_t> buf;
    if (a.serializer == "packbits") {
        buf = pack_bits(data, n);
    } else {
        buf.assign(data, data + n * a.itemsize);
        if (a.big_endian) {
            byteswap(buf.data(), n, a.dtype == DType::FixedUnicode ? 4 : a.itemsize);
        }
    }
    encode_and_write(a, a.chunk_path(index), std::move(buf));
}

static std::string
ucs4_to_utf8(const uint8_t *p, size_t nchars)
{
    std::string out;
    for (size_t j = 0; j < nchars; j++) {
        uint32_t cp;
        std::memcpy(&cp, p + 4 * j, 4);
        if (cp == 0) {
            break;
        }
        if (cp < 0x80) {
            out += static_cast<char>(cp);
        } else if (cp < 0x800) {
            out += static_cast<char>(0xC0 | (cp >> 6));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        } else if (cp < 0x10000) {
            out += static_cast<char>(0xE0 | (cp >> 12));
            out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        } else {
            out += static_cast<char>(0xF0 | (cp >> 18));
            out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
            out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        }
    }
    return out;
}

bool
read_chunk_strings(
    const ZarrArray &a, const std::vector<int64_t> &index, std::vector<std::string> &out)
{
    size_t n = a.chunk_elements();
    if (a.dtype == DType::FixedUnicode || a.dtype == DType::FixedBytes) {
        std::vector<uint8_t> raw;
        bool found = read_chunk(a, index, raw);
        out.assign(n, std::string());
        if (!found) {
            std::fill(out.begin(), out.end(), a.fill_string);
            return false;
        }
        for (size_t j = 0; j < n; j++) {
            const uint8_t *p = raw.data() + j * a.itemsize;
            if (a.dtype == DType::FixedUnicode) {
                out[j] = ucs4_to_utf8(p, a.itemsize / 4);
            } else {
                size_t len = 0;
                while (len < a.itemsize && p[len] != 0) {
                    len++;
                }
                out[j] = std::string(reinterpret_cast<const char *>(p), len);
            }
        }
        return true;
    }
    std::string path = a.chunk_path(index);
    if (!file_exists(path)) {
        out.assign(n, a.fill_string);
        return false;
    }
    std::vector<uint8_t> buf = decode_bytes(a, path, 0);
    if (buf.size() < 4) {
        throw std::runtime_error("vlen-utf8: truncated chunk '" + path + "'");
    }
    uint32_t count;
    std::memcpy(&count, buf.data(), 4);
    if (count != n) {
        throw std::runtime_error("vlen-utf8: chunk '" + path + "' has "
                                 + std::to_string(count) + " items, expected "
                                 + std::to_string(n));
    }
    out.assign(n, std::string());
    size_t pos = 4;
    for (size_t j = 0; j < n; j++) {
        uint32_t len;
        if (pos + 4 > buf.size()) {
            throw std::runtime_error("vlen-utf8: truncated chunk '" + path + "'");
        }
        std::memcpy(&len, buf.data() + pos, 4);
        pos += 4;
        if (pos + len > buf.size()) {
            throw std::runtime_error("vlen-utf8: truncated chunk '" + path + "'");
        }
        out[j].assign(reinterpret_cast<const char *>(buf.data() + pos), len);
        pos += len;
    }
    return true;
}

void
write_chunk_strings(
    const ZarrArray &a, const std::vector<int64_t> &index, const std::vector<std::string> &data)
{
    size_t n = a.chunk_elements();
    if (data.size() != n) {
        throw std::runtime_error("write_chunk_strings: wrong number of items");
    }
    if (a.dtype == DType::FixedBytes || a.dtype == DType::FixedUnicode) {
        std::vector<uint8_t> raw(n * a.itemsize, 0);
        for (size_t j = 0; j < n; j++) {
            uint8_t *p = raw.data() + j * a.itemsize;
            if (a.dtype == DType::FixedBytes) {
                std::memcpy(p, data[j].data(), std::min(data[j].size(), a.itemsize));
            } else {
                size_t k = 0;
                size_t i = 0;
                const std::string &s = data[j];
                while (i < s.size() && k < a.itemsize / 4) {
                    unsigned char c = static_cast<unsigned char>(s[i]);
                    uint32_t cp;
                    size_t len;
                    if (c < 0x80) {
                        cp = c;
                        len = 1;
                    } else if ((c & 0xE0) == 0xC0) {
                        cp = c & 0x1F;
                        len = 2;
                    } else if ((c & 0xF0) == 0xE0) {
                        cp = c & 0x0F;
                        len = 3;
                    } else {
                        cp = c & 0x07;
                        len = 4;
                    }
                    for (size_t m = 1; m < len && i + m < s.size(); m++) {
                        cp = (cp << 6) | (static_cast<unsigned char>(s[i + m]) & 0x3F);
                    }
                    std::memcpy(p + 4 * k, &cp, 4);
                    i += len;
                    k++;
                }
            }
        }
        write_chunk(a, index, raw.data());
        return;
    }
    std::vector<uint8_t> buf;
    uint32_t count = static_cast<uint32_t>(n);
    buf.resize(4);
    std::memcpy(buf.data(), &count, 4);
    for (const std::string &s : data) {
        uint32_t len = static_cast<uint32_t>(s.size());
        size_t pos = buf.size();
        buf.resize(pos + 4 + len);
        std::memcpy(buf.data() + pos, &len, 4);
        std::memcpy(buf.data() + pos + 4, s.data(), len);
    }
    encode_and_write(a, a.chunk_path(index), std::move(buf));
}

bool
read_chunk_int32(const ZarrArray &a, const std::vector<int64_t> &index, std::vector<int32_t> &out)
{
    std::vector<uint8_t> raw;
    bool found = read_chunk(a, index, raw);
    size_t n = a.chunk_elements();
    out.resize(n);
    switch (a.dtype) {
        case DType::Int8:
            for (size_t j = 0; j < n; j++) {
                out[j] = static_cast<int8_t>(raw[j]);
            }
            break;
        case DType::Int16:
            for (size_t j = 0; j < n; j++) {
                int16_t x;
                std::memcpy(&x, raw.data() + 2 * j, 2);
                out[j] = x;
            }
            break;
        case DType::Int32:
            std::memcpy(out.data(), raw.data(), 4 * n);
            break;
        default:
            for (size_t j = 0; j < n; j++) {
                out[j] = static_cast<int32_t>(element_as_int64(a, raw.data() + j * a.itemsize));
            }
    }
    return found;
}

} // namespace tsir
