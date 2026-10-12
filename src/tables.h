#ifndef TSINFERR_TABLES_H
#define TSINFERR_TABLES_H

#include "common.h"

namespace tsir {

/* Sort, index, optionally compute mutation parents and check the tables form
 * a valid tree sequence. */
void finalise_tables(tsk_table_collection_t *tables, bool compute_mutation_parents);
void add_provenance(tsk_table_collection_t *tables, const std::string &record);

} // namespace tsir

#endif
