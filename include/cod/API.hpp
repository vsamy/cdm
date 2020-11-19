#pragma once

#if defined _WIN32 || defined __CYGWIN__
#define COD_DLLIMPORT __declspec(dllimport)
#define COD_DLLEXPORT __declspec(dllexport)
#define COD_DLLLOCAL
#else
#if __GNUC__ >= 4
#define COD_DLLIMPORT __attribute__((visibility("default")))
#define COD_DLLEXPORT __attribute__((visibility("default")))
#define COD_DLLLOCAL __attribute__((visibility("hidden")))
#else
#define COD_DLLIMPORT
#define COD_DLLEXPORT
#define COD_DLLLOCAL
#endif
#endif

#ifdef COD_STATIC
#define COD_DLLAPI
#define COD_LOCAL
#else
#ifdef COD_EXPORTS
#define COD_DLLAPI COD_DLLEXPORT
#else
#define COD_DLLAPI COD_DLLIMPORT
#endif
#define COD_LOCAL COD_DLLLOCAL
#endif