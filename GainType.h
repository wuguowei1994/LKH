#ifndef _GAINTYPE_H
#define _GAINTYPR_H

#define HAVE_LONG_LONG

#include "LKH.h"
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
#define LLONG_MIN (-LLONG_MAX -1LL)
#endif
#define PLUS_INFINITY LLONG_MAX
#define MINUS_INFINITY LLONG_MIN
#define GainFormat "%lld"
#define GainInputFormat "%lld"
#endif

#endif