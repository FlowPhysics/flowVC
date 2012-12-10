/*
 *  io.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "io.h"

void FatalError(char *text, ...) {
	va_list ap;
	char *p, *sval;
	int ival;
	double dval;
	
	fprintf(stderr, "\nERROR:\n");
	va_start(ap, text);
	for(p = text; *p; p++) {
		if(*p != '%') {
			putc(*p, stderr);
			continue;
		}
		switch (*++p) {
			case 'd':
				ival = va_arg(ap, int);
				fprintf(stderr, "%d", ival);
				break;
			case 'f':
				dval = va_arg(ap, double);
				fprintf(stderr, "%f", dval);
				break;
			case 's':
				for(sval = va_arg(ap, char *); *sval; sval++)
					putc(*sval, stderr);
				break;
			default:
				putc(*p, stderr);
				break;
		}
	}
	va_end(ap);
	
	fprintf(stderr, "\n");
	fflush(stderr);
	
	exit(1);
}

