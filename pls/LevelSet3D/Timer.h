
/**************************************************************************
	ORIGINAL AUTHOR: 
		Emud Mokhberi (emud@ucla.edu)
	MODIFIED BY:
	
	CONTRIBUTORS:
		

-----------------------------------------------
	
 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

 This software is distributed for noncommercial use in the hope that it will 
 be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
 to anyone for the consequences of using it or for whether it serves any 
 particular purpose or works at all. No guarantee is made about the software 
 or its performance.

 You are allowed to modify the source code, add your name to the
 appropriate list above and distribute the code as long as 
 this license agreement is distributed with the code and is included at 
 the top of all header (.h) files.

 Commercial use is strictly prohibited.
***************************************************************************/

#ifndef __TIMER__
#define __TIMER__

#ifdef _WIN32
#include <windows.h>
#include "main.h"

class Timer
{
public:
	Timer();

	//in seconds
	Float GetElapsedTime();
	void Reset();
private:
	LONGLONG cur_time;

	DWORD time_count;
	LONGLONG perf_cnt;
	bool perf_flag;
	LONGLONG last_time;
	Double time_scale;

	bool QPC;
};




#else
//*****************************unix stuff****************************
#include <sys/time.h>


class Timer
{
public:
	Timer();

	void Reset();
	Float GetElapsedTime();

private:
	timeval cur_time;

};

#endif //unix

#endif // File
