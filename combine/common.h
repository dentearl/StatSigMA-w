#ifndef COMMON_H_
#define COMMON_H_

#include <stdarg.h>

void verbose(char const *fmt, ...);
void debug(char const *fmt, ...);
void message(char const *type, char const *fmt, ...);

#endif // COMMON_H_
