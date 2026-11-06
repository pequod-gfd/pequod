#ifndef __QG3_CONFIG_H
#define __QG3_CONFIG_H

#include "debug.h"
#include "version.h"

! Enable / disable terms
#define WIND_FORCING
#undef BACKGROUND_VSHEAR
#undef BACKGROUND_HSHEAR

! CABARET options
#define CABARET
#define CABARET_LIMITER

! Centred options
#define RA_DOUBLE_FILTER
! Legacy centred options - changed only by tests
#undef CENTRED_BOTTOM_FRICTION
#define FORWARD_EULER_START

! Elliptic solver options
#define PH1D
#define PH1D_BUFFER
#undef INTERMEDIATE_TRANSPOSES
#undef PETSC
#undef SUPERLU_DIST

#endif
