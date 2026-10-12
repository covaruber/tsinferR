#include "common.h"
#include "tables.h"

#include <cmath>
#include <ctime>

using namespace Rcpp;

namespace tsir {

static std::string
ragged_at(const char *data, const tsk_size_t *offset, tsk_size_t j)
{
    return std::string(data + offset[j], data + offset[j + 1]);
}

static CharacterVector
ragged_column(const char *data, const tsk_size_t *offset, tsk_size_t n)
{
    CharacterVector out(n);
    for (tsk_size_t j = 0; j < n; j++) {
        std::string s = ragged_at(data, offset, j);
        out[j] = Rf_mkCharLenCE(s.data(), static_cast<int>(s.size()), CE_UTF8);
    }
    return out;
}

static std::string
str_of(const char *data, tsk_size_t length)
{
    return std::string(data == nullptr ? "" : data, data == nullptr ? 0 : length);
}

void
finalise_tables(tsk_table_collection_t *tables, bool compute_mutation_parents)
{
    check_tsk(tsk_table_collection_sort(tables, nullptr, 0), "sort");
    check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    if (compute_mutation_parents) {
        check_tsk(tsk_table_collection_compute_mutation_parents(tables, 0),
            "compute_mutation_parents");
    }
    tsk_id_t ret = tsk_table_collection_check_integrity(tables, TSK_CHECK_TREES);
    if (ret < 0) {
        stop_tsk(static_cast<int>(ret), "tree sequence integrity check");
    }
}

void
add_provenance(tsk_table_collection_t *tables, const std::string &record)
{
    char timestamp[64];
    time_t now = time(nullptr);
    struct tm tm_utc;
#if defined(_WIN32)
    gmtime_s(&tm_utc, &now);
#else
    gmtime_r(&now, &tm_utc);
#endif
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
    check_tsk(tsk_provenance_table_add_row(&tables->provenances, timestamp,
                  strlen(timestamp), record.data(), record.size()),
        "add provenance");
}

} // namespace tsir

using tsir::check_tsk;
using tsir::get_tables;
using tsir::new_tables_xptr;

// [[Rcpp::export(.tc_new)]]
SEXP
tc_new(double sequence_length)
{
    tsk_table_collection_t *tables;
    SEXP xp = PROTECT(new_tables_xptr(&tables));
    check_tsk(tsk_table_collection_init(tables, 0), "init");
    tables->sequence_length = sequence_length;
    UNPROTECT(1);
    return xp;
}

// [[Rcpp::export(.tc_load)]]
SEXP
tc_load(std::string path)
{
    tsk_table_collection_t *tables;
    SEXP xp = PROTECT(new_tables_xptr(&tables));
    int err = tsk_table_collection_load(tables, path.c_str(), 0);
    if (err != 0) {
        UNPROTECT(1);
        tsir::stop_tsk(err, path.c_str());
    }
    if (!tsk_table_collection_has_index(tables, 0)) {
        check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    }
    UNPROTECT(1);
    return xp;
}

// [[Rcpp::export(.tc_dump)]]
void
tc_dump(SEXP xp, std::string path)
{
    tsk_table_collection_t *tables = get_tables(xp);
    if (!tsk_table_collection_has_index(tables, 0)) {
        check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    }
    check_tsk(tsk_table_collection_dump(tables, path.c_str(), 0), path.c_str());
}

// [[Rcpp::export(.tc_copy)]]
SEXP
tc_copy(SEXP xp)
{
    tsk_table_collection_t *src = get_tables(xp);
    tsk_table_collection_t *tables;
    SEXP out = PROTECT(new_tables_xptr(&tables));
    check_tsk(tsk_table_collection_copy(src, tables, 0), "copy");
    UNPROTECT(1);
    return out;
}

// [[Rcpp::export(.tc_equals)]]
bool
tc_equals(SEXP a, SEXP b, bool ignore_provenance, bool ignore_timestamps)
{
    tsk_flags_t options = 0;
    if (ignore_provenance) {
        options |= TSK_CMP_IGNORE_PROVENANCE;
    }
    if (ignore_timestamps) {
        options |= TSK_CMP_IGNORE_TIMESTAMPS;
    }
    return tsk_table_collection_equals(get_tables(a), get_tables(b), options);
}

// [[Rcpp::export(.tc_info)]]
List
tc_info(SEXP xp)
{
    tsk_table_collection_t *t = get_tables(xp);
    return List::create(_["sequence_length"] = t->sequence_length,
        _["num_nodes"] = static_cast<double>(t->nodes.num_rows),
        _["num_edges"] = static_cast<double>(t->edges.num_rows),
        _["num_sites"] = static_cast<double>(t->sites.num_rows),
        _["num_mutations"] = static_cast<double>(t->mutations.num_rows),
        _["num_individuals"] = static_cast<double>(t->individuals.num_rows),
        _["num_populations"] = static_cast<double>(t->populations.num_rows),
        _["num_provenances"] = static_cast<double>(t->provenances.num_rows),
        _["metadata"] = tsir::str_of(t->metadata, t->metadata_length),
        _["metadata_schema"] = tsir::str_of(t->metadata_schema, t->metadata_schema_length),
        _["time_units"] = tsir::str_of(t->time_units, t->time_units_length));
}

// [[Rcpp::export(.tc_table)]]
List
tc_table(SEXP xp, std::string name)
{
    tsk_table_collection_t *t = get_tables(xp);
    if (name == "nodes") {
        const tsk_node_table_t &n = t->nodes;
        return List::create(
            _["flags"] = NumericVector(n.flags, n.flags + n.num_rows),
            _["time"] = NumericVector(n.time, n.time + n.num_rows),
            _["population"] = IntegerVector(n.population, n.population + n.num_rows),
            _["individual"] = IntegerVector(n.individual, n.individual + n.num_rows),
            _["metadata"] = tsir::ragged_column(n.metadata, n.metadata_offset, n.num_rows),
            _["metadata_schema"] = tsir::str_of(n.metadata_schema, n.metadata_schema_length));
    }
    if (name == "edges") {
        const tsk_edge_table_t &e = t->edges;
        return List::create(_["left"] = NumericVector(e.left, e.left + e.num_rows),
            _["right"] = NumericVector(e.right, e.right + e.num_rows),
            _["parent"] = IntegerVector(e.parent, e.parent + e.num_rows),
            _["child"] = IntegerVector(e.child, e.child + e.num_rows));
    }
    if (name == "sites") {
        const tsk_site_table_t &s = t->sites;
        return List::create(
            _["position"] = NumericVector(s.position, s.position + s.num_rows),
            _["ancestral_state"] = tsir::ragged_column(
                s.ancestral_state, s.ancestral_state_offset, s.num_rows),
            _["metadata"] = tsir::ragged_column(s.metadata, s.metadata_offset, s.num_rows));
    }
    if (name == "mutations") {
        const tsk_mutation_table_t &m = t->mutations;
        return List::create(_["site"] = IntegerVector(m.site, m.site + m.num_rows),
            _["node"] = IntegerVector(m.node, m.node + m.num_rows),
            _["parent"] = IntegerVector(m.parent, m.parent + m.num_rows),
            _["time"] = NumericVector(m.time, m.time + m.num_rows),
            _["derived_state"]
            = tsir::ragged_column(m.derived_state, m.derived_state_offset, m.num_rows),
            _["metadata"] = tsir::ragged_column(m.metadata, m.metadata_offset, m.num_rows));
    }
    if (name == "individuals") {
        const tsk_individual_table_t &i = t->individuals;
        return List::create(_["flags"] = NumericVector(i.flags, i.flags + i.num_rows),
            _["metadata"] = tsir::ragged_column(i.metadata, i.metadata_offset, i.num_rows),
            _["metadata_schema"] = tsir::str_of(i.metadata_schema, i.metadata_schema_length));
    }
    if (name == "populations") {
        const tsk_population_table_t &p = t->populations;
        return List::create(
            _["metadata"] = tsir::ragged_column(p.metadata, p.metadata_offset, p.num_rows),
            _["metadata_schema"] = tsir::str_of(p.metadata_schema, p.metadata_schema_length));
    }
    if (name == "provenances") {
        const tsk_provenance_table_t &p = t->provenances;
        return List::create(
            _["timestamp"] = tsir::ragged_column(p.timestamp, p.timestamp_offset, p.num_rows),
            _["record"] = tsir::ragged_column(p.record, p.record_offset, p.num_rows));
    }
    stop("Unknown table: %s", name);
}

static std::vector<char>
pack_ragged(CharacterVector x, std::vector<tsk_size_t> &offset)
{
    std::vector<char> data;
    // tskit rejects a NULL data pointer paired with offsets.
    data.reserve(1);
    offset.assign(x.size() + 1, 0);
    for (R_xlen_t j = 0; j < x.size(); j++) {
        const char *s = Rf_translateCharUTF8(STRING_ELT(x, j));
        size_t len = strlen(s);
        data.insert(data.end(), s, s + len);
        offset[j + 1] = static_cast<tsk_size_t>(data.size());
    }
    return data;
}

// [[Rcpp::export(.tc_set_table)]]
void
tc_set_table(SEXP xp, std::string name, List cols)
{
    tsk_table_collection_t *t = get_tables(xp);
    if (name == "nodes") {
        NumericVector flags_d = cols["flags"];
        NumericVector time = cols["time"];
        IntegerVector population = cols["population"];
        IntegerVector individual = cols["individual"];
        CharacterVector md = cols["metadata"];
        std::vector<tsk_flags_t> flags(flags_d.begin(), flags_d.end());
        std::vector<tsk_size_t> off;
        std::vector<char> data = pack_ragged(md, off);
        check_tsk(tsk_node_table_set_columns(&t->nodes, flags.size(), flags.data(),
                      time.begin(), population.begin(), individual.begin(), data.data(),
                      off.data()),
            "set nodes");
    } else if (name == "edges") {
        NumericVector left = cols["left"], right = cols["right"];
        IntegerVector parent = cols["parent"], child = cols["child"];
        check_tsk(tsk_edge_table_set_columns(&t->edges, left.size(), left.begin(),
                      right.begin(), parent.begin(), child.begin(), nullptr, nullptr),
            "set edges");
    } else if (name == "sites") {
        NumericVector position = cols["position"];
        CharacterVector anc = cols["ancestral_state"];
        std::vector<tsk_size_t> off;
        std::vector<char> data = pack_ragged(anc, off);
        check_tsk(tsk_site_table_set_columns(&t->sites, position.size(), position.begin(),
                      data.data(), off.data(), nullptr, nullptr),
            "set sites");
    } else if (name == "mutations") {
        IntegerVector site = cols["site"], node = cols["node"], parent = cols["parent"];
        NumericVector time_r = cols["time"];
        CharacterVector derived = cols["derived_state"];
        std::vector<double> time(time_r.begin(), time_r.end());
        for (double &x : time) {
            if (std::isnan(x)) {
                x = TSK_UNKNOWN_TIME;
            }
        }
        std::vector<tsk_size_t> off;
        std::vector<char> data = pack_ragged(derived, off);
        check_tsk(tsk_mutation_table_set_columns(&t->mutations, site.size(), site.begin(),
                      node.begin(), parent.begin(), time.data(), data.data(), off.data(),
                      nullptr, nullptr),
            "set mutations");
    } else {
        stop("Unsupported table: %s", name);
    }
}

// [[Rcpp::export(.tc_set_metadata)]]
void
tc_set_metadata(SEXP xp, std::string schema, std::string metadata)
{
    tsk_table_collection_t *t = get_tables(xp);
    check_tsk(tsk_table_collection_set_metadata_schema(t, schema.data(), schema.size()),
        "metadata schema");
    check_tsk(tsk_table_collection_set_metadata(t, metadata.data(), metadata.size()),
        "metadata");
}

// [[Rcpp::export(.tc_finalise)]]
void
tc_finalise(SEXP xp, bool compute_mutation_parents)
{
    tsir::finalise_tables(get_tables(xp), compute_mutation_parents);
}

// [[Rcpp::export(.tc_add_provenance)]]
void
tc_add_provenance(SEXP xp, std::string record)
{
    tsir::add_provenance(get_tables(xp), record);
}

// [[Rcpp::export(.tc_num_trees)]]
double
tc_num_trees(SEXP xp)
{
    tsk_table_collection_t *t = get_tables(xp);
    tsk_treeseq_t ts;
    int err = tsk_treeseq_init(&ts, t, TSK_TS_INIT_BUILD_INDEXES);
    if (err != 0) {
        tsk_treeseq_free(&ts);
        tsir::stop_tsk(err, "tree sequence");
    }
    double n = static_cast<double>(tsk_treeseq_get_num_trees(&ts));
    tsk_treeseq_free(&ts);
    return n;
}

// [[Rcpp::export(.tc_genotype_matrix)]]
List
tc_genotype_matrix(SEXP xp)
{
    tsk_table_collection_t *t = get_tables(xp);
    tsk_treeseq_t ts;
    tsk_variant_t var;
    int err = tsk_treeseq_init(&ts, t, TSK_TS_INIT_BUILD_INDEXES);
    if (err != 0) {
        tsk_treeseq_free(&ts);
        tsir::stop_tsk(err, "tree sequence");
    }
    tsk_size_t num_samples = tsk_treeseq_get_num_samples(&ts);
    tsk_size_t num_sites = tsk_treeseq_get_num_sites(&ts);
    err = tsk_variant_init(&var, &ts, nullptr, 0, nullptr, 0);
    if (err != 0) {
        tsk_variant_free(&var);
        tsk_treeseq_free(&ts);
        tsir::stop_tsk(err, "variant");
    }
    IntegerMatrix G(static_cast<int>(num_sites), static_cast<int>(num_samples));
    List alleles(num_sites);
    for (tsk_size_t j = 0; j < num_sites; j++) {
        err = tsk_variant_decode(&var, static_cast<tsk_id_t>(j), 0);
        if (err != 0) {
            break;
        }
        for (tsk_size_t k = 0; k < num_samples; k++) {
            G(j, k) = var.genotypes[k];
        }
        CharacterVector a(var.num_alleles);
        for (tsk_size_t k = 0; k < var.num_alleles; k++) {
            a[k] = std::string(var.alleles[k], var.allele_lengths[k]);
        }
        alleles[j] = a;
    }
    IntegerVector samples(ts.samples, ts.samples + num_samples);
    tsk_variant_free(&var);
    tsk_treeseq_free(&ts);
    if (err != 0) {
        tsir::stop_tsk(err, "decode");
    }
    return List::create(_["genotypes"] = G, _["alleles"] = alleles, _["samples"] = samples);
}

// [[Rcpp::export(.tc_simplify)]]
SEXP
tc_simplify(SEXP xp, Nullable<IntegerVector> samples, bool filter_sites,
    bool filter_populations, bool filter_individuals, bool keep_unary)
{
    tsk_table_collection_t *src = get_tables(xp);
    tsk_table_collection_t *tables;
    SEXP out = PROTECT(new_tables_xptr(&tables));
    check_tsk(tsk_table_collection_copy(src, tables, 0), "copy");
    std::vector<tsk_id_t> sample_ids;
    if (samples.isNotNull()) {
        IntegerVector s(samples);
        sample_ids.assign(s.begin(), s.end());
    } else {
        for (tsk_size_t j = 0; j < tables->nodes.num_rows; j++) {
            if (tables->nodes.flags[j] & TSK_NODE_IS_SAMPLE) {
                sample_ids.push_back(static_cast<tsk_id_t>(j));
            }
        }
    }
    tsk_flags_t options = 0;
    if (filter_sites) {
        options |= TSK_SIMPLIFY_FILTER_SITES;
    }
    if (filter_populations) {
        options |= TSK_SIMPLIFY_FILTER_POPULATIONS;
    }
    if (filter_individuals) {
        options |= TSK_SIMPLIFY_FILTER_INDIVIDUALS;
    }
    if (keep_unary) {
        options |= TSK_SIMPLIFY_KEEP_UNARY;
    }
    check_tsk(tsk_table_collection_simplify(tables, sample_ids.data(),
                  static_cast<tsk_size_t>(sample_ids.size()), options, nullptr),
        "simplify");
    check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    UNPROTECT(1);
    return out;
}
