/********************************************************8
*3drdllconfig.h
*author:europelee
*date:2009-7-12
************************************************************/
#ifndef MYCONFIG_H
#define MYCONFIG_H

#include <algorithm>
#include <functional>
#include <vector>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <float.h> 

#define ImHeight 128
#define ImWidth 256
#define DmSize 256
#define SWD_TWOPI	6.28318530717958647692f
#define M_PI 3.14159265358979323846
#define PI 3.141592653589793f
#define TWOSQRTPI 3.5449078083
#define SQRTTWOPI  2.5066282749
#define		MAX(a, b)       ((a) > (b) ? (a) : (b))			//!<	Returns the max value between a and b
#define		MAXMAX(a,b,c)   ((a) > (b) ? MAX (a,c) : MAX (b,c))
#define TRUE 1
#define FALSE 0
#endif
