/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

int main(int, char *[]);

#define NUMFUNCTION	4
#define EXPRLEN		3
#define NUMDIM		3
#define VECTORLEN	10
#define VARMIX		2

char Array_Indices[NUMDIM][3]	= {"I", "J", "K"};
char Function[NUMFUNCTION]	= {'+', '-', '*', '/'};
char Indent[4]			= {"  "};

int main(int argc, char *argv[])
{
   int  i, j, firstpass = 1, fnamecount = 0, timed = 0;
   int  form, function, count, dim, exprlen, strptr, tstptr, vlen;
   char filename[256], progname[256], exprstr[256], tststr[256];
   char srcString[1024], objString[1024], exeString[1024], *ptr, *typestr;
   char ObjName[NUMDIM * NUMFUNCTION * EXPRLEN][32], 
        ExeName[NUMDIM * NUMFUNCTION * EXPRLEN][32];
   FILE *OutFile;
   FILE *MakeFile;
   FILE *Run_Tests;

   if (argc >= 2)
   {

/* See if the user wants a specific type of field in the tests.
 * Handling only double, float, and int for now
 */
      if (!(strcmp("double", argv[1])))
      {
         typestr = strdup("double");
      }
      else if (!(strcmp("int", argv[1])))
      {
         typestr = strdup("int");
      }
      else if (!(strcmp("float", argv[1])))
      {
         typestr = strdup("float");
      }
      else
      {
         printf("Error!  Unknown field type %s!  Quitting\n", argv[1]);
         return -1;
      }
   }
   else
   {
/* Default type is double for now */

      typestr = strdup("double");
   }
   printf("Generating %s field tests\n", typestr);

   if (argc == 3)
      if (!(strcmp("timed", argv[2])))
      {
         timed = 1;
      }
      else
      {
         printf("Unknown third parameter!  Only acceptible parameter is 'timed'\n");
         printf("Quitting...\n");
         return -1;
      }

   srcString[0] = '\0';
   objString[0] = '\0';
   exeString[0] = '\0';

   strcat(srcString, "TSTSRC = ");
   strcat(objString, "TSTOBJ = ");
   strcat(exeString, "TSTEXE = ");
   for (dim = 1 ; dim <= NUMDIM ; dim++)
   {
      count = 0;
      vlen = (int)pow(VECTORLEN, dim);
      for (exprlen = 1 ; exprlen <= EXPRLEN ; exprlen++)
      {
         for (function = 0 ; function < NUMFUNCTION ; function++)
         {
/* Generate the file name */

            sprintf(progname, "%sTest%1.1d%2.2d", typestr, dim, count);

/* Open the file */

            sprintf(filename, "%s.C", progname);
            if ((OutFile = fopen(filename, "w")) == NULL)
            {
               fprintf(stderr, "Can't open file %s for write!\n", filename);
               return 0;
            }
            if (!firstpass) strcat(srcString, "\t");
            strcat(srcString, filename);
            if ((fnamecount + 1) < EXPRLEN * NUMDIM * NUMFUNCTION)
               strcat(srcString, " \\\n");
            else
               strcat(srcString, " \n");
            if (!firstpass) strcat(objString, "\t");
            ptr = strrchr(filename, 'C');
            *ptr = 'o';
            strcat(objString, filename);
            strcpy(ObjName[fnamecount], filename);
            if ((fnamecount + 1) < EXPRLEN * NUMDIM * NUMFUNCTION)
               strcat(objString, " \\\n");
            else
               strcat(objString, " \n");
            if (!firstpass) strcat(exeString, "\t");
            strcat(exeString, progname);
            strcpy(ExeName[fnamecount++], progname);
            if ((fnamecount) < EXPRLEN * NUMDIM * NUMFUNCTION)
               strcat(exeString, " \\\n");
            else
               strcat(exeString, " \n");
            firstpass = 0;
            if (timed)
               fprintf(OutFile, "#define BOUNDS_CHECK_DEFAULT false\n");
            else
               fprintf(OutFile, "#define BOUNDS_CHECK_DEFAULT true\n");
            fprintf(OutFile, "\n#include <iostream.h>\n");
            fprintf(OutFile, "#include \"Field.h\"\n");
            fprintf(OutFile, "#include \"Expressions.h\"\n");
            if (timed)
            {
               fprintf(OutFile, "#include \"Timer.h\"\n");
               fprintf(OutFile, "#define MILLION 1e6\n");
            }
            fprintf(OutFile, "\nint main()\n");
            fprintf(OutFile, "{\n  const unsigned Dim = %d;\n", dim);
            for (i = 0 ; i < dim ; i++)
               fprintf(OutFile,"  Index %s(%d);\n",Array_Indices[i],VECTORLEN);

            fprintf(OutFile, "  NDIndex<Dim> domain;\n");

            for (i = 0 ; i < dim ; i++)
               fprintf(OutFile,"  domain[%d] = %s;\n", i, Array_Indices[i]);

            fprintf(OutFile, "  FieldLayout<Dim> layout(domain);\n");
            fprintf(OutFile, "  Field<%s, Dim> A(layout);\n", typestr);
            fprintf(OutFile, "  Field<%s, Dim>::iterator p;\n", typestr);
            if (timed)
            {
               fprintf(OutFile, "  Timer Time;\n");
               fprintf(OutFile, "  unsigned Flops;\n");
               fprintf(OutFile, "  float Etime;\n");
            }
            fprintf(OutFile, "  %s a[%d], *ptr;\n", typestr, vlen, dim);
            fprintf(OutFile, "\n  int i;\n");
            fprintf(OutFile, "\n  ptr = a;\n");
            
            fprintf(OutFile, "\n  for (p = A.begin(), i = 0 ; p != A.end() ; ++p, ++i)");
            fprintf(OutFile, "\n  {\n    *p = i + 1;\n    *ptr++ = i + 1;\n  }\n");
            strptr = 0; tstptr = 0;
            tststr[tstptr++] = '*';
            tststr[tstptr++] = 'p';
            tststr[tstptr++] = 't';
            tststr[tstptr++] = 'r';
            tststr[tstptr++] = ' ';
            exprstr[strptr++] = 'A';
            exprstr[strptr++] = ' ';
            exprstr[strptr++] = Function[function % NUMFUNCTION];
            tststr[tstptr++] = Function[function % NUMFUNCTION];
            exprstr[strptr++] = ' ';
            tststr[tstptr++] = ' ';
            for (i = 0 ; i < exprlen ; i++)
            { 
               exprstr[strptr++] = 'A';
               exprstr[strptr++] = ' ';
               tststr[tstptr++] = '*';
               tststr[tstptr++] = 'p';
               tststr[tstptr++] = 't';
               tststr[tstptr++] = 'r';
               tststr[tstptr++] = ' ';
               if (i < (exprlen - 1)) 
               {
                  exprstr[strptr++] = Function[i % NUMFUNCTION];
                  exprstr[strptr++] = ' ';
                  tststr[tstptr++] = Function[i % NUMFUNCTION];
                  tststr[tstptr++] = ' ';
               }
            }
            exprstr[strptr] = '\0';
            tststr[tstptr] = '\0';
            if (timed)
            {
               fprintf(OutFile, "  Time.clear();\n");
               fprintf(OutFile, "  Time.start();\n");
            }
            fprintf(OutFile, "\n  assign(A, %s);\n", exprstr);
            if (timed)
            {
               fprintf(OutFile, "  Time.stop();\n");
               fprintf(OutFile, "  Etime = Time.clock_time();\n");
               fprintf(OutFile, "  Flops = %d / Etime;\n", vlen);
            }
            fprintf(OutFile, "\n  ptr = a;\n");
            fprintf(OutFile, "\n  for (i = 0 ; i < %d ; i++)", vlen);
            fprintf(OutFile, "\n  {\n    *ptr++ = %s;\n  }\n", tststr);
            fprintf(OutFile, "\n  ptr = a;\n");
            fprintf(OutFile, "\n  for (p = A.begin(), i = 0 ; p != A.end() ; ++p, ++i)\n");
            fprintf(OutFile, "  {\n    if (*ptr++ != *p)\n    {\n");
            fprintf(OutFile, "        cout << \"Test %s FAILED!\\n\";\n", progname);
            fprintf(OutFile, "        cout << \"\\nExpected : \";\n");
            fprintf(OutFile, "        for (i = 0 ; i < %d ; i++)\n", vlen);
            fprintf(OutFile, "          cout << *ptr++ << \" \" ;\n");
            fprintf(OutFile, "        cout << endl;\n\n");
            fprintf(OutFile, "        cout << \"\\nActual   : \";\n");
            fprintf(OutFile, "        cout << A << endl;\n");

            fprintf(OutFile, "        return -1;\n");
            fprintf(OutFile, "    }\n  }\n\n");
            fprintf(OutFile, "  cout << \"%s Passed. \";\n", progname);
            if (timed)
            {
               fprintf(OutFile, "  cout << \"MFlops = \" << Flops/MILLION << endl;\n", progname);
            }
            else
            {
               fprintf(OutFile, "  cout << endl;\n");
            }
             
            fprintf(OutFile, "\n  return 0;\n}\n");
            fclose(OutFile);
            count++;
         }
      }
   }
/*
 * Write out the Makefile for the generated tests.
 */

   if ((MakeFile = fopen("Makefile.tests", "w")) == NULL)
   {
      fprintf(stderr, "Can't open file 'Makefile.tests' for write!\n");
      return -1;
   }
   fprintf(MakeFile, "%s\n", srcString);
   fprintf(MakeFile, "%s\n", objString);
   fprintf(MakeFile, "%s\n", exeString);
   fprintf(MakeFile, "include ../../lib/SGI/Makefile.app\n\n");
   fprintf(MakeFile, "LINK = $(LINKER) $(LFLAGS) $@.o -o $@ $(LIBS)\n\n");
   fprintf(MakeFile, "all : $(TSTEXE)\n\n");
   fprintf(MakeFile, "\nclean : \n\trm -rf $(TSTEXE)\n\trm -rf $(TSTOBJ)\n\n");
   for (i = 0 ; i < NUMDIM * NUMFUNCTION * EXPRLEN ; i++)
   {
      fprintf(MakeFile, "%s : %s\n\t\$(LINK)\n\n", 
         ExeName[i],
         ObjName[i]);
   }
   fclose (MakeFile);

/*
 * Write out the script to run the generated tests.
 */

   if ((Run_Tests = fopen("RunTests", "w")) == NULL)
   {
      fprintf(stderr, "Can't open file 'RunTests' for write!\n");
      return -1;
   }
   fprintf(Run_Tests, "#!/bin/sh\n");
   fprintf(Run_Tests, "#\n# Self-Generated script to run tests\n#\n");
   for (i = 0 ; i < NUMDIM * NUMFUNCTION * EXPRLEN ; i++)
      fprintf(Run_Tests, "%s \n", ExeName[i]);

   fclose (Run_Tests);

   system("chmod +x RunTests");
   return 0;
}
/***************************************************************************
 * $RCSfile: Generate_Tests.c,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:07:09 $
 * IPPL_VERSION_ID: $Id: Generate_Tests.c,v 1.1.1.1 2000/11/30 21:07:09 adelmann Exp $ 
 ***************************************************************************/
