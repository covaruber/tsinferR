#!/usr/bin/env bash
# Development helper: regenerate exports/docs, install into a local library
# and run the test suite. Usage: tools/dev.sh [test-filter]
set -euo pipefail
PKG_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LIB="${TSINFERR_DEV_LIB:-$PKG_DIR/../../rlib}"
mkdir -p "$LIB"
cd "$PKG_DIR"
Rscript -e 'Rcpp::compileAttributes("."); suppressMessages(roxygen2::roxygenise(".", roclets = c("rd", "namespace")))' > /dev/null
R CMD INSTALL --no-docs -l "$LIB" . 2>&1 | grep -E "error|warning:|DONE|ERROR" || true
cd tests
FILTER="${1:-}"
R_LIBS="$LIB" Rscript -e "library(testthat); library(tsinferR); test_dir('testthat', reporter = 'summary', filter = if (nzchar('$FILTER')) '$FILTER' else NULL)"
