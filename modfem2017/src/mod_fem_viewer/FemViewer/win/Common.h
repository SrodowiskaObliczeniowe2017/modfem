#pragma once
#ifndef Common_H_
#define Common_H_


#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>

#define _ATL_CSTRING_EXPLICIT_CONSTRUCTORS      // some CString constructors will be explicit

#include <atlbase.h>
#include <atlwin.h>


#define VERIFY  ATLVERIFY
#define ASSERT  ATLASSERT

#define V   VERIFY
#define S   SUCCEEDED

#endif /* Common_H_
*/