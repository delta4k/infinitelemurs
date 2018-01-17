#pragma once

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

// Turn of deprecation warning so that we can use posix names
// for console io functions
#pragma warning (disable : 4996)

// I really like the C idiom of testing ints and pointers as booleans,
// so tell the compiler not to warn when I do it.
#pragma warning (disable : 4800 )

#include <string>
#include <chrono>

void failf(const char* fmt, ...);
#define Assert(condition) if(!(condition)) failf("Assertion failure: %s Line %d\n", __FILE__, __LINE__)

template<typename X> void Swap(X &x0, X &x1) { X tmp; tmp = x0; x0 = x1; x1 = tmp; }

// macros that I happen to like in spite of macros not being
// C++ politically correct. They're simple and they work.

#define Between(a,b,c)   ((b)>=(a) && (b)<=(c) ? 1 : 0)
#define For(i,lim) for(int i=0;i<(lim);i++)
#define Max(a,b) ((a)>(b)?(a):(b))
#define Min(a,b) ((a)>(b)?(b):(a))
#define Abs(x) (((x)<0)?(-(x)):(x))

size_t getFileSize(const char* filename);
std::string loadText(const char* fname);

extern bool verboseOutput;
void prn(const char* fmt, ...);
std::string formatString(const char* fmt, ...);

double fractionalDifference(double x, double y);
int checksum(const std::string& str);

/*****************************************************************************\ 
 
    Timing mechanism. Not relevant to the actual simulation algorithm.
 
\*****************************************************************************/
class Stopwatch
{

private:

    std::chrono::steady_clock::time_point startTime;

public:

    Stopwatch();

    void start();   // Start accumulating, may be called multiple times

    double elapsedSeconds() const;
};
////////////////////////////////////////////////////////////////////////////////

