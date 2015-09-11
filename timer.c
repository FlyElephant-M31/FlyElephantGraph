#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

struct rusage r1, r2;
// time_t tinicio, tfim;
struct timeval tinicio, tfim;

float timer()
{
   struct rusage r;
   getrusage(0, &r);
   return (float)(r.ru_utime.tv_sec+r.ru_utime.tv_usec/(float)1000000);
}

void timer1() 
{
    getrusage(RUSAGE_SELF, &r1);
    gettimeofday(&tinicio, NULL); 
}

void timer2()
{
    getrusage(RUSAGE_SELF, &r2);
    gettimeofday(&tfim, NULL);
}

double getUserTime ()
{
  double utime1, utime2;
  
  utime1 = (double) r1.ru_utime.tv_sec + (double) r1.ru_utime.tv_usec/1000000.0;
  utime2 = (double) r2.ru_utime.tv_sec + (double) r2.ru_utime.tv_usec/1000000.0;
  return (utime2 - utime1);

}

double getSystemTime ()
{
  double stime1, stime2;

  stime1 = (double) r1.ru_stime.tv_sec + (double) r1.ru_stime.tv_usec/1000000.0;
  stime2 = (double) r2.ru_stime.tv_sec + (double) r2.ru_stime.tv_usec/1000000.0;
  return (stime2 - stime1);
}

double getRealTime() 
{
  double t1 = (double)tinicio.tv_sec*1000000.0 + (double)tinicio.tv_usec;
  double t2 = (double)tfim.tv_sec*1000000.0 + (double)tfim.tv_usec;
  return (t2-t1)/1000000.0;
}

