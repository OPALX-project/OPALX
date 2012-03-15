/* Header file for the timing module:
      Name -- name of timed routine
      Fout -- file name where results may be printed 

                     Zarija Lukic, November 2008
                           zarija@lanl.gov
*/

#ifndef PerfMon_Header_Included
#define PerfMon_Header_Included

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

#ifdef __cplusplus
extern "C" {
#endif	

   void StartMonitor();
   void OnClock(const char *Name);
   void OffClock(const char *Name);
   void ResetClock(const char *Name);
   long unsigned UsesOfRoutine(const char *Name);
   double TimeOfRoutine(const char *Name);
   void PrintClock(const char *Name, FILE *Fout);
   void PrintClockSummary(FILE *Fout);

#ifdef __cplusplus
}
#endif
      
#ifdef INITIALIZER_NAMESPACE
}
#endif

#endif
