// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>
// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>
#include <stdio.h>
// TODO: reference additional headers your program requires here
#include "gl\gl.h"
#include "gl\glu.h"
#include "gl\glaux.h"

#pragma warning(disable: 4305) /*truncation from 'double' to 'float'*/
#pragma warning(disable: 4244) /*conversion from 'xxx' to 'yyy', possible loss of data*/
