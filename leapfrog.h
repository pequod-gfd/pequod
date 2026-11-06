#ifndef __QG3_LEAPFROG_H
#define __QG3_LEAPFROG_H

#ifdef CABARET
#define DOLEAPFROG
#define ENDDOLEAPFROG
#define BLEAPFROG
#define BNLEAPFROG
#define DIMLEAPFROG
#define DIMSLEAPFROG
#define TLEAPFROG
#define TSLEAPFROG
#define TNLEAPFROG
#define LLEAPFROG
#else
#define DOLEAPFROG do leapfrog = 1, 2
#define ENDDOLEAPFROG end do
#define BLEAPFROG (leapfrog)
#define BNLEAPFROG (2)
#define DIMLEAPFROG , dimension(2)
#define DIMSLEAPFROG , dimension(:)
#define TLEAPFROG , leapfrog
#define TSLEAPFROG , :
#define TNLEAPFROG , 2
#define LLEAPFROG leapfrog,
#endif

#endif
