#pragma once

#if defined _WIN32 || defined __CYGWIN__
#define CDM_DLLIMPORT __declspec(dllimport)
#define CDM_DLLEXPORT __declspec(dllexport)
#define CDM_DLLLOCAL
#else
#if __GNUC__ >= 4
#define CDM_DLLIMPORT __attribute__((visibility("default")))
#define CDM_DLLEXPORT __attribute__((visibility("default")))
#define CDM_DLLLOCAL __attribute__((visibility("hidden")))
#else
#define CDM_DLLIMPORT
#define CDM_DLLEXPORT
#define CDM_DLLLOCAL
#endif
#endif

#ifdef CDM_STATIC
#define CDM_DLLAPI
#define CDM_LOCAL
#else
#ifdef CDM_EXPORTS
#define CDM_DLLAPI CDM_DLLEXPORT
#else
#define CDM_DLLAPI CDM_DLLIMPORT
#endif
#define CDM_LOCAL CDM_DLLLOCAL
#endif