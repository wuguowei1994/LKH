#ifndef _GAINTYPE_H
#define _GAINTYPE_H


#define HAVE_LONG_LONG
/* Undefine if you don't have the long long type */
/* #undef HAVE_LONG_LONG */
/**
 * #ifdef 标识符
 * 程序段1
 * #else
 * 程序段2
 * #endif
 * 它的功能是，如果标识符已被 #define命令定义过则对程序段1进行编译；
 * 否则对程序段2进行编译。
 * 如果没有程序段2(它为空)，本格式中的#else可以没有，即可以写为：
 * #ifdef 标识符
 * 程序段
 * #endif
 */
#include <float.h>
#include <limits.h>
#ifdef HAVE_LONG_LONG
typedef long long GainType;
/**
 * 将“ifdef”改为“ifndef”。
 * 它的作用是：若标识符未被定义则编译程序段1，否则编译程序段2。
 * 这种形式与第一种形式的作用相反。
 */
#ifndef LLONG_MAX
#define LLONG_MAX 9223372036854775807LL
#endif
#ifndef LLONG_MIN
#define LLONG_MIN (-LLONG_MAX - 1LL)
#endif
#define PLUS_INFINITY LLONG_MAX
#define MINUS_INFINITY LLONG_MIN
#define GainFormat "%lld"
#define GainInputFormat "%lld"
#else
typedef double GainType;
#define PLUS_INFINITY DBL_MAX
#define MINUS_INFINITY -DBL_MAX
#define GainFormat "%0.0lf"
#define GainInputFormat "%0.0lf"
#endif

#endif
