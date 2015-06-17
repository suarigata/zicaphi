/***************************************************************************
 *
 * Copyright 2015 Ian Liu Rodrigues and Edson Borin
 *
 ***************************************************************************/

#include "log.h"

#include <stdio.h>
#include <stdarg.h>

void log_progress(float percentage, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	printf("[%3d%%] ", (int)(percentage * 100));
	vprintf(fmt, ap);
	printf("\033[0K\r");
	fflush(stdout);
	va_end(ap);
}
