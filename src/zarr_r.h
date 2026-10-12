#ifndef TSINFERR_ZARR_R_H
#define TSINFERR_ZARR_R_H

#include <Rcpp.h>

#include "zarr.h"

namespace tsir {

ZarrArray zarr_array_from_list(Rcpp::List meta);
ZarrArray *get_zarr_array(SEXP xp);

} // namespace tsir

#endif
