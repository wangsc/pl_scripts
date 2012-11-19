#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "INLINE.h"

void greet() {
  printf("Hello, world!
");
}

MODULE = hello_pl_07fc	PACKAGE = main	

PROTOTYPES: DISABLE


void
greet ()
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	greet();
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

