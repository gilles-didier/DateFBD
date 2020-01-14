#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "MyError.h"


void error(const char *message, ...) {
	va_list args;

	va_start(args, message);
	vprintf(message, args);
	va_end(args);
	exit(1);
}

void warning(const char *message, ...) {
	va_list args;

	va_start(args, message);
	vprintf(message, args);
	va_end(args);
}
