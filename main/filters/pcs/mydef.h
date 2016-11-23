/* "$Id: mydef.h,v 1.2 1995/08/04 04:51:23 huber Exp huber $" */

#ifndef MYDEF_H
#define MYDEF_H

/* declaration of variables */

#define MINT    int
#define MLONG   long int

#ifdef USE_SINGLE
#define MFLOAT  float
#define SFO     f
#define SGO     g
#define SEO     e
#define SFO6_2   6.2f
#define SFO8_3   8.3f
#define SFO10_5  10.5f
#define SFO12_3  12.3f
#else
#define MFLOAT  double
#define SFO     lf
#define SGO     lg
#define SEO     le
#define SFO6_2   6.2lf
#define SFO8_3   8.3lf
#define SFO10_5  10.5lf
#define SFO12_3  12.3lf
#endif

#define STR(a) STRING(a)
#define STRING(a) #a
  
typedef struct RPoint {
  MFLOAT x;
  MFLOAT y;
  MFLOAT z;
} RPoint;

#endif
