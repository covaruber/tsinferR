/* Builds the bundled c-blosc library with its diagnostic output disabled, so
 * no compiled code writes to the console (CRAN policy). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define printf(...)            ((void) 0)
#define fprintf(...)           ((void) 0)
#define sprintf(buf, ...)      snprintf(buf, 256, __VA_ARGS__)

#include "blosc/blosc.c"
