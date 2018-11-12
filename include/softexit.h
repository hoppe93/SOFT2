#ifndef _SOFTEXIT_H
#define _SOFTEXIT_H

#if defined __GNUG__ || defined __GNUC__
void softexit(int) __attribute__((noreturn));
#else
void softexit(int);
#endif

#endif/*_SOFTEXIT_H*/
