#ifndef __QG3_DEBUG_H
#define __QG3_DEBUG_H

#ifndef NDEBUG
#define assert(cmd) if(.not. (cmd)) call pabort("Assertion failure: cmd", __FILE__, __LINE__)
#else
#define assert(cmd)
#endif
#define labort(msg) call pabort(msg, __FILE__, __LINE__)

#endif
