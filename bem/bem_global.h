#ifndef BEM_GLOBAL_H
#define BEM_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(BEM_LIBRARY)
#  define BEMSHARED_EXPORT Q_DECL_EXPORT
#else
#  define BEMSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // BEM_GLOBAL_H
