/* Code's performance monitor tools 
   See the PerfMon.h for the complete list of routines
   Example usage:
                   #include "PerfMon.h"
                   void main(){
				   ... declarations ...
				   StartMonitor();
				   ... code ...
                   OnClock("Routine1")
                   ... routine 1 code ...
                   OffClock("Routine1")
                   OnClock("Routine2")
				   ... routine 1 code ...
                   OffClock("Routine2")
                   PrintClockSummary(Filename)
                   }
   To reset any timer to 0.0 seconds        : ResetClock(char *Name)
   To print a result of a particular routine: PrintClock(char *Name, FILE *Filename)
 
                                      Zarija Lukic, November 2008
        				                   zarija@lanl.gov
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "PerfMon.h"

#define MaxClocks 50

typedef struct{const char *RName; int Count; double ElapsedTime, StartTime;} Timer;

static Timer Clock[MaxClocks+1];
static int TotClocksUsed;
static char *PerfMonStartTime;

static void PerfMonError(char *Message){
   fprintf(stderr, "ERROR in PerfMon: %s\n", Message);
	return;
}

int FindModuleIndex(const char *Name){
   int i;
   for (i=0; i<MaxClocks+1; ++i)
	   if (Clock[i].RName == Name) return(i);
	if (TotClocksUsed == MaxClocks){
		Clock[MaxClocks].RName = "Rest";
		i = MaxClocks;
		return(i);
	}
	else{
		i = TotClocksUsed+1;
		Clock[i].RName = Name;
		Clock[i].Count = 0;
		Clock[i].ElapsedTime = 0.0;
		return(i);
	}
}

void StartMonitor(){
	time_t start_t;
	start_t = time(NULL);
	TotClocksUsed = -1;
	PerfMonStartTime = asctime(localtime(&start_t));
	return;
}

void OnClock(const char *Name){
	int Index;
	Index = FindModuleIndex(Name);
	++Clock[Index].Count;
	Clock[Index].StartTime = (double)clock()/CLOCKS_PER_SEC;
	if (Index > TotClocksUsed) TotClocksUsed=Index;
	return;
}

void OffClock(const char *Name){
	int Index;
	Index = FindModuleIndex(Name);
	if (Clock[Index].Count == 0)
		PerfMonError("OffClock called for a routine which was never on!");
	else {
		Clock[Index].ElapsedTime += ((double)clock()/CLOCKS_PER_SEC - Clock[Index].StartTime);
	}
	return;
}

void ResetClock(const char *Name){
	int Index;
	Index = FindModuleIndex(Name);
	Clock[Index].ElapsedTime = 0.0;
	return;
}

long unsigned UsesOfRoutine(const char *Name){
	int Index;
	Index = FindModuleIndex(Name);
	return(Clock[Index].Count);
}

double TimeOfRoutine(const char *Name){
	int Index;
	Index = FindModuleIndex(Name);
	return(Clock[Index].ElapsedTime);
}

void PrintClock(const char *Name, FILE *Fout){
	fprintf(Fout, "Routine %s was called %lu times and used for a total "
			"time of %f seconds \n", Name, UsesOfRoutine(Name), TimeOfRoutine(Name));
	return;
}

void PrintClockSummary(FILE *Fout){
	int i;
	double total, frac;
	time_t end_t;
	end_t = time(NULL);
	total = (double)clock()/CLOCKS_PER_SEC;
	fprintf(Fout, "\n\n====================================================================\n");
	fprintf(Fout, "Initializer: code performance summary \n");
	fprintf(Fout, "             code beginning : %s", PerfMonStartTime);
	fprintf(Fout, "             code ending    : %s", asctime(localtime(&end_t)));
	fprintf(Fout, "             total run time : %G seconds \n", total);
	fprintf(Fout, "====================================================================\n");
	if (TotClocksUsed < 0)
		fprintf(Fout, "\n No part of the code was timed. \n\n");
	else{
	    fprintf(Fout, "\nRoutine            # of calls       Time [s]       Percent of total \n");
		fprintf(Fout, "--------------------------------------------------------------------\n");
    	for (i=0; i<=TotClocksUsed; ++i){
	    	frac = Clock[i].ElapsedTime*100.0/total;
		    fprintf(Fout, "%-20s %-12i   %-8.4G           %-8.4G \n", Clock[i].RName, 
					Clock[i].Count, Clock[i].ElapsedTime, frac);
		}
	}
		fprintf(Fout, "====================================================================\n\n");
	return;
}
