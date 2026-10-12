#include "common.h"
#include "threads.h"
#include "zarr_r.h"

#include <cmath>
#include <limits>
#include <map>
#include <zlib.h>

using namespace Rcpp;

namespace {

/* Line reader over plain or (b)gzipped text, via zlib. */
class LineReader {
  public:
    explicit LineReader(const std::string &path)
    {
        file_ = gzopen(path.c_str(), "rb");
        if (file_ == nullptr) {
            throw std::runtime_error("cannot open '" + path + "'");
        }
        gzbuffer(file_, 1 << 20);
    }
    ~LineReader()
    {
        if (file_ != nullptr) {
            gzclose(file_);
        }
    }
    LineReader(const LineReader &) = delete;
    LineReader &operator=(const LineReader &) = delete;

    bool next(std::string &line)
    {
        line.clear();
        char buf[65536];
        while (true) {
            if (gzgets(file_, buf, sizeof(buf)) == nullptr) {
                int err;
                const char *msg = gzerror(file_, &err);
                if (err != Z_OK && err != Z_STREAM_END) {
                    throw std::runtime_error(std::string("error reading VCF: ") + msg);
                }
                return !line.empty();
            }
            line += buf;
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                return true;
            }
        }
    }

  private:
    gzFile file_ = nullptr;
};

void
split(const std::string &s, char sep, std::vector<std::string> &out)
{
    out.clear();
    size_t start = 0;
    while (true) {
        size_t end = s.find(sep, start);
        if (end == std::string::npos) {
            out.emplace_back(s, start);
            return;
        }
        out.emplace_back(s, start, end - start);
        start = end + 1;
    }
}

/* Value of KEY=... inside a ##header=<...> line. */
std::string
header_attr(const std::string &line, const std::string &key)
{
    size_t open = line.find('<');
    if (open == std::string::npos) {
        return "";
    }
    size_t pos = open + 1;
    while (pos < line.size()) {
        size_t eq = line.find('=', pos);
        if (eq == std::string::npos) {
            return "";
        }
        std::string k = line.substr(pos, eq - pos);
        size_t vstart = eq + 1;
        size_t vend;
        if (vstart < line.size() && line[vstart] == '"') {
            vend = line.find('"', vstart + 1);
            if (vend == std::string::npos) {
                return "";
            }
            if (k == key) {
                return line.substr(vstart + 1, vend - vstart - 1);
            }
            vend++;
        } else {
            vend = line.find_first_of(",>", vstart);
            if (vend == std::string::npos) {
                vend = line.size();
            }
            if (k == key) {
                return line.substr(vstart, vend - vstart);
            }
        }
        pos = vend + 1;
    }
    return "";
}

struct GenotypeChunks {
    const tsir::ZarrArray &gt;
    const tsir::ZarrArray &phased;
    const tsir::ZarrArray &mask;
    size_t num_samples;
    int ploidy;
    size_t vchunk;
    std::vector<int8_t> calls;  /* vchunk x samples x ploidy */
    std::vector<uint8_t> phase; /* vchunk x samples */
    size_t rows = 0;
    int64_t chunk = 0;

    GenotypeChunks(const tsir::ZarrArray &g, const tsir::ZarrArray &p, const tsir::ZarrArray &m,
        size_t n, int pl, size_t vc)
        : gt(g), phased(p), mask(m), num_samples(n), ploidy(pl), vchunk(vc),
          calls(vc * n * static_cast<size_t>(pl), -1), phase(vc * n, 0)
    {
    }

    void flush(int num_threads)
    {
        if (rows == 0) {
            return;
        }
        const size_t schunk = static_cast<size_t>(gt.chunks[1]);
        const size_t nsc = (num_samples + schunk - 1) / schunk;
        const size_t pl = static_cast<size_t>(ploidy);
        tsir::parallel_for(nsc, num_threads, [&](size_t sc) {
            std::vector<int8_t> g(vchunk * schunk * pl, -1);
            std::vector<uint8_t> m(vchunk * schunk * pl, 1);
            std::vector<uint8_t> p(vchunk * schunk, 0);
            for (size_t r = 0; r < rows; r++) {
                for (size_t c = 0; c < schunk; c++) {
                    size_t s = sc * schunk + c;
                    if (s >= num_samples) {
                        break;
                    }
                    for (size_t k = 0; k < pl; k++) {
                        int8_t v = calls[(r * num_samples + s) * pl + k];
                        g[(r * schunk + c) * pl + k] = v;
                        m[(r * schunk + c) * pl + k] = v < 0 ? 1 : 0;
                    }
                    p[r * schunk + c] = phase[r * num_samples + s];
                }
            }
            tsir::write_chunk(gt, { chunk, static_cast<int64_t>(sc), 0 },
                reinterpret_cast<const uint8_t *>(g.data()));
            tsir::write_chunk(mask, { chunk, static_cast<int64_t>(sc), 0 }, m.data());
            tsir::write_chunk(phased, { chunk, static_cast<int64_t>(sc) }, p.data());
        });
        std::fill(calls.begin(), calls.end(), static_cast<int8_t>(-1));
        std::fill(phase.begin(), phase.end(), 0);
        rows = 0;
        chunk++;
    }
};

} // namespace

// [[Rcpp::export(.vcf_header)]]
List
vcf_header(std::string path)
{
    std::vector<std::string> contig_id, filter_id, filter_desc, info_id, info_type, info_number;
    std::vector<double> contig_length;
    std::vector<std::string> samples;
    int ploidy = 0;
    try {
        LineReader reader(path);
        std::string line;
        std::vector<std::string> fields;
        while (reader.next(line)) {
            if (line.rfind("##contig=", 0) == 0) {
                contig_id.push_back(header_attr(line, "ID"));
                std::string len = header_attr(line, "length");
                contig_length.push_back(len.empty() ? NA_REAL : std::atof(len.c_str()));
            } else if (line.rfind("##FILTER=", 0) == 0) {
                filter_id.push_back(header_attr(line, "ID"));
                filter_desc.push_back(header_attr(line, "Description"));
            } else if (line.rfind("##INFO=", 0) == 0) {
                info_id.push_back(header_attr(line, "ID"));
                info_type.push_back(header_attr(line, "Type"));
                info_number.push_back(header_attr(line, "Number"));
            } else if (line.rfind("#CHROM", 0) == 0) {
                split(line, '\t', fields);
                for (size_t j = 9; j < fields.size(); j++) {
                    samples.push_back(fields[j]);
                }
                /* Detect the ploidy from the genotypes of the first records. */
                std::vector<std::string> format;
                int scanned = 0;
                while (scanned < 1000 && reader.next(line)) {
                    if (line.empty() || line[0] == '#') {
                        continue;
                    }
                    scanned++;
                    split(line, '\t', fields);
                    if (fields.size() < 10) {
                        continue;
                    }
                    split(fields[8], ':', format);
                    size_t gt_index = 0;
                    while (gt_index < format.size() && format[gt_index] != "GT") {
                        gt_index++;
                    }
                    if (gt_index == format.size()) {
                        continue;
                    }
                    for (size_t j = 9; j < fields.size(); j++) {
                        split(fields[j], ':', format);
                        if (gt_index < format.size()) {
                            int p = 1;
                            for (char c : format[gt_index]) {
                                if (c == '/' || c == '|') {
                                    p++;
                                }
                            }
                            ploidy = std::max(ploidy, p);
                        }
                    }
                }
                break;
            } else if (line.empty() || line[0] != '#') {
                throw std::runtime_error("VCF header has no #CHROM line");
            }
        }
    } catch (std::exception &e) {
        stop(e.what());
    }
    return List::create(_["contig_id"] = wrap(contig_id),
        _["contig_length"] = wrap(contig_length), _["filter_id"] = wrap(filter_id),
        _["filter_description"] = wrap(filter_desc), _["info_id"] = wrap(info_id),
        _["info_type"] = wrap(info_type), _["info_number"] = wrap(info_number),
        _["samples"] = wrap(samples), _["ploidy"] = ploidy);
}

// [[Rcpp::export(.vcf_convert)]]
List
vcf_convert(std::string path, List gt_meta, List phased_meta, List mask_meta,
    int num_samples, int ploidy, std::string contig, CharacterVector info_fields,
    CharacterVector filter_ids, int num_threads)
{
    tsir::ZarrArray gt = tsir::zarr_array_from_list(gt_meta);
    tsir::ZarrArray phased = tsir::zarr_array_from_list(phased_meta);
    tsir::ZarrArray mask = tsir::zarr_array_from_list(mask_meta);
    const size_t vchunk = static_cast<size_t>(gt.chunks[0]);
    const int threads = tsir::resolve_num_threads(num_threads);
    GenotypeChunks chunks(gt, phased, mask, static_cast<size_t>(num_samples), ploidy, vchunk);

    std::vector<std::string> chrom_names;
    std::map<std::string, int> chrom_index;
    std::vector<int> v_contig;
    std::vector<double> v_pos, v_qual;
    std::vector<std::string> v_id;
    std::vector<std::vector<std::string>> v_alleles;
    std::vector<std::string> info_names(info_fields.begin(), info_fields.end());
    std::vector<std::vector<std::string>> v_info(info_names.size());
    std::map<std::string, int> filter_index;
    for (R_xlen_t j = 0; j < filter_ids.size(); j++) {
        filter_index[as<std::string>(filter_ids[j])] = static_cast<int>(j);
    }
    std::vector<std::vector<int>> v_filters;
    size_t max_ploidy_seen = 0;
    try {
        LineReader reader(path);
        std::string line;
        std::vector<std::string> fields, sub, info, kv, alts, filters;
        std::vector<int> ptr_gt;
        size_t line_no = 0;
        while (reader.next(line)) {
            line_no++;
            if (line.empty() || line[0] == '#') {
                continue;
            }
            split(line, '\t', fields);
            if (fields.size() < 8) {
                throw std::runtime_error(
                    "VCF line " + std::to_string(line_no) + " has fewer than 8 columns");
            }
            if (!contig.empty() && fields[0] != contig) {
                continue;
            }
            if (num_samples > 0 && fields.size() != static_cast<size_t>(9 + num_samples)) {
                throw std::runtime_error("VCF line " + std::to_string(line_no)
                                         + " does not have one column per sample");
            }
            auto ci = chrom_index.find(fields[0]);
            if (ci == chrom_index.end()) {
                ci = chrom_index.emplace(fields[0], static_cast<int>(chrom_names.size())).first;
                chrom_names.push_back(fields[0]);
            }
            v_contig.push_back(ci->second);
            v_pos.push_back(std::atof(fields[1].c_str()));
            v_id.push_back(fields[2] == "." ? "" : fields[2]);
            std::vector<std::string> alleles{ fields[3] };
            if (fields[4] != ".") {
                split(fields[4], ',', alts);
                alleles.insert(alleles.end(), alts.begin(), alts.end());
            }
            v_alleles.push_back(std::move(alleles));
            v_qual.push_back(fields[5] == "." ? std::numeric_limits<double>::quiet_NaN()
                                              : std::atof(fields[5].c_str()));
            std::vector<int> fv;
            if (fields[6] != ".") {
                split(fields[6], ';', filters);
                for (const std::string &f : filters) {
                    auto fi = filter_index.find(f);
                    if (fi != filter_index.end()) {
                        fv.push_back(fi->second);
                    }
                }
            }
            v_filters.push_back(std::move(fv));
            if (!info_names.empty()) {
                split(fields[7], ';', info);
                for (size_t k = 0; k < info_names.size(); k++) {
                    std::string value;
                    bool found = false;
                    for (const std::string &item : info) {
                        size_t eq = item.find('=');
                        std::string key = eq == std::string::npos ? item : item.substr(0, eq);
                        if (key == info_names[k]) {
                            value = eq == std::string::npos ? "\x01" : item.substr(eq + 1);
                            found = true;
                            break;
                        }
                    }
                    v_info[k].push_back(found ? value : std::string("\x02"));
                }
            }
            /* Genotypes. */
            size_t row = chunks.rows;
            if (num_samples > 0) {
                split(fields[8], ':', sub);
                int gt_index = -1;
                for (size_t k = 0; k < sub.size(); k++) {
                    if (sub[k] == "GT") {
                        gt_index = static_cast<int>(k);
                    }
                }
                for (int s = 0; s < num_samples; s++) {
                    const std::string &cell = fields[9 + s];
                    int8_t *dest = chunks.calls.data()
                                   + (row * static_cast<size_t>(num_samples) + s) * ploidy;
                    if (gt_index < 0) {
                        continue;
                    }
                    /* Locate the GT subfield without splitting the whole cell. */
                    size_t start = 0;
                    for (int k = 0; k < gt_index && start != std::string::npos; k++) {
                        start = cell.find(':', start);
                        if (start != std::string::npos) {
                            start++;
                        }
                    }
                    if (start == std::string::npos) {
                        continue;
                    }
                    size_t end = cell.find(':', start);
                    if (end == std::string::npos) {
                        end = cell.size();
                    }
                    int k = 0;
                    bool all_phased = true;
                    size_t pos = start;
                    while (pos < end) {
                        if (k >= ploidy) {
                            throw std::runtime_error("VCF line " + std::to_string(line_no)
                                                     + ": genotype ploidy exceeds "
                                                     + std::to_string(ploidy));
                        }
                        int value = -1;
                        if (cell[pos] == '.') {
                            pos++;
                        } else {
                            long v = 0;
                            size_t digits = 0;
                            while (pos < end && cell[pos] >= '0' && cell[pos] <= '9') {
                                v = v * 10 + (cell[pos] - '0');
                                pos++;
                                digits++;
                            }
                            if (digits == 0) {
                                throw std::runtime_error(
                                    "VCF line " + std::to_string(line_no) + ": bad genotype '"
                                    + cell.substr(start, end - start) + "'");
                            }
                            if (v > 127) {
                                throw std::runtime_error(
                                    "VCF line " + std::to_string(line_no)
                                    + ": allele indexes above 127 are not supported");
                            }
                            value = static_cast<int>(v);
                        }
                        dest[k] = static_cast<int8_t>(value);
                        k++;
                        if (pos < end) {
                            if (cell[pos] == '/') {
                                all_phased = false;
                            } else if (cell[pos] != '|') {
                                throw std::runtime_error(
                                    "VCF line " + std::to_string(line_no) + ": bad genotype '"
                                    + cell.substr(start, end - start) + "'");
                            }
                            pos++;
                        }
                    }
                    max_ploidy_seen = std::max(max_ploidy_seen, static_cast<size_t>(k));
                    /* Calls with fewer alleles than the ploidy are padded with -2. */
                    for (int p = k; p < ploidy; p++) {
                        dest[p] = -2;
                    }
                    chunks.phase[row * static_cast<size_t>(num_samples) + s]
                        = (k > 1 && all_phased) ? 1 : 0;
                }
            }
            chunks.rows++;
            if (chunks.rows == vchunk) {
                chunks.flush(threads);
                Rcpp::checkUserInterrupt();
            }
        }
        chunks.flush(threads);
    } catch (std::exception &e) {
        stop(e.what());
    }
    size_t n = v_pos.size();
    size_t max_alleles = 1;
    for (const auto &a : v_alleles) {
        max_alleles = std::max(max_alleles, a.size());
    }
    CharacterMatrix alleles(static_cast<int>(n), static_cast<int>(std::max<size_t>(max_alleles, 2)));
    std::fill(alleles.begin(), alleles.end(), "");
    for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < v_alleles[j].size(); k++) {
            alleles(static_cast<int>(j), static_cast<int>(k)) = v_alleles[j][k];
        }
    }
    LogicalMatrix filt(static_cast<int>(n), static_cast<int>(filter_ids.size()));
    for (size_t j = 0; j < n; j++) {
        for (int f : v_filters[j]) {
            filt(static_cast<int>(j), f) = TRUE;
        }
    }
    List info_out(info_names.size());
    for (size_t k = 0; k < info_names.size(); k++) {
        info_out[k] = wrap(v_info[k]);
    }
    info_out.names() = info_fields;
    return List::create(_["contig_names"] = wrap(chrom_names), _["contig"] = wrap(v_contig),
        _["position"] = wrap(v_pos), _["id"] = wrap(v_id), _["alleles"] = alleles,
        _["quality"] = wrap(v_qual), _["filters"] = filt, _["info"] = info_out,
        _["num_variants"] = static_cast<double>(n),
        _["max_ploidy"] = static_cast<double>(max_ploidy_seen));
}
