
/*****************************************************************************\

    Copyright (c) 2011-2017

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

\*****************************************************************************/


#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <conio.h>
#include <cstdarg>
#include <mutex>
#include <set>
#include "utilities.h"

bool verboseOutput = true;

// printf is only guaranteed to handle strings up to this length.
// The Microsoft compiler handles much more than this, but
// we want this to be as portable as reasonably possible.
const int MAX_PRINTF_STRLEN = 4095;
static char formatBuffer[MAX_PRINTF_STRLEN + 1];
static std::mutex formatMutex;

#define FORMAT(fmt, buffer, bufflen)                                    \
    va_list argptr;                                                     \
    va_start(argptr, fmt);                                              \
    int num_required = vsnprintf((buffer), (bufflen), (fmt), argptr);   \
    va_end(argptr);                                                     \
    Assert(num_required < bufflen)


/***********************************************\

    Console output can slow a program down a lot.
    This is a printf-substitute that can be turned
    on or off.
 
    Use this only for stuff that you want to suppress
    when running flat-out for speed. Use plain-old
    printf for stuff that you want to always be printed.

\***********************************************/
void prn(const char* fmt, ...)
{
    if(!verboseOutput)
        return;

    std::lock_guard<std::mutex> lock{ formatMutex };
    FORMAT(fmt, formatBuffer, MAX_PRINTF_STRLEN);

    printf("%s", formatBuffer);
}

/***********************************************\

    Because printf-style formatting is just too
    nice to do without.

\***********************************************/
std::string formatString(const char* fmt, ...)
{

    std::lock_guard<std::mutex> lock{ formatMutex };
    FORMAT(fmt, formatBuffer, MAX_PRINTF_STRLEN);
    std::string result(formatBuffer);

    return result;
}

/*****************************************************************************\
    Return the size of the given file, or zero if the file can't be obtained
\*****************************************************************************/
size_t getFileSize(const char* filename) 
{
    struct stat st;
    if (stat(filename, &st) != 0) 
    {
        return 0;
    }
    return st.st_size;
}
/*****************************************************************************\
 
    Load the given file into a single contiguous string.

    This function won't return if the load operation fails. It will fail for
    the usual reasons file operations fail, or if the file contains embedded    
    zeros or characters that won't fit in a signed char.

\*****************************************************************************/
std::string loadText(const char* fname)
{
    std::string result;
    size_t size = getFileSize(fname);

    if(size > result.max_size())
    {
        failf("File %s is too large for a single string", fname);
    }

    if(size == 0)
    {
        failf("Failure loading text file %s\n", fname);
    }

    // We use "rb" because we want an exact copy, with no
    // line ending translation.
    FILE* fp = fopen(fname, "rb");
    if (!fp) 
    {
        failf("Failure loading text file %s\n", fname);
    }

    result.reserve(size);

    char c = fgetc(fp);
    int charnum = 0;
    while(!feof(fp))
    {
        if(c == 0)
        {
            failf("File contains embedded zero character in position %d\n", charnum);
        }

        if (c > 0x7E)
        {
            failf("Text contains non-printing character: (char)%d in character #%d\n", c, charnum);
        }

        result.push_back(c);
        c = fgetc(fp);
        charnum++;
    }

    return result;
}

/*****************************************************************************\
    
    A simple, quick checksum, not intended for security or UUID purposes.

\*****************************************************************************/
int checksum(const std::string& str)
{
    int result = 0;

    for(auto c : str)
    {
        result += c;
    }

    return result;
}
/*****************************************************************************\ 
 
    Fail sort of gracefully with an error message.
 
\*****************************************************************************/
void failf(const char* fmt, ...)
{   
    // Something bad has happened, therefore we won't assume we can use
    // the shared format buffer.
    const int bufsize = 255;
    char failbuf[bufsize + 1];

    FORMAT(fmt, failbuf, bufsize);

    printf("Failure: %s\n", failbuf);
    printf("Hit any key to exit..\n");
    getch();
    exit(-1);
}

/*****************************************************************************\
\*****************************************************************************/
Stopwatch::Stopwatch()
{
    clock();
    startTime = std::chrono::steady_clock::now();
}

void Stopwatch::start()
{
    startTime = std::chrono::steady_clock::now();
}

double Stopwatch::elapsedSeconds() const
{
    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> dt = now - startTime;
    double elapsed_seconds = dt.count();
    return elapsed_seconds;
}

/*****************************************************************************\
    Return the fractional difference between x and y: Abs(x-y) / Max(x,y)
\*****************************************************************************/
double fractionalDifference(double x, double y)
{
    double maxval = Max(Abs(x), Abs(y));
    if (maxval == 0) return 0;

    double delta = Abs(x - y);

    // This is not intended as a general-purpose function that may
    // be called with pathological inputs, so we're going to be
    // pragmatic and not stress about things like overflows, except 
    // for a basic sanity check.
    Assert(maxval > 1e-10);
    return delta / maxval;
}
