#ifndef _LKH_H
#define _LKH_H
// undef的作用:在undef关键字后面取消以前定义的宏定义
#undef NDEBUG
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "GainType.h"

/**
 * 下面这些变量用于ReadParameters()这个函数读取输入文件
 */
char *ParameterFileName;
void ReadParameters(void);
void hello(void);

#endif