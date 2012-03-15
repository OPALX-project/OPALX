/* error.C
   error handling

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: error.C 91 2007-05-06 17:47:45Z bakker $
*/

#define SVN_DATE "$Date: 2007-05-06 19:47:45 +0200 (Sun, 06 May 2007) $"


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>

#include "Algorithms/bet/error.h"


void initErrorMsg(int level,const char *fbase) {
  printf("\n\n error.cpp: initErrorMesg \n\n");
} /* initErrorMsg() */

void initErrorFilename(const char *fbase) {
  printf("\n\n error.cpp: initErrorFilename \n\n");
} /* initErrorFileName() */

void setReportLevel(int level) {
  printf("\n\n error.cpp: setReportLevel \n\n");
}

void writeError() {
  printf("\n\n error.cpp: writeError \n\n");
}

