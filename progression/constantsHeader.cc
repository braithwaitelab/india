#define _CRT_SECURE_NO_DEPRECATE 1	//eliminates printf deprecation warnings
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;

#include "constantsHeader.h"



static bool stripSpace(const char c) 
{
	return (c == ' ' || c == '\t');
}

static bool stripReturn(const char c) 
{
	return (c == '\n');
}

string createConstantsHeader(void)
{
	string headerString;
	string line;
	int insideComment=0;
	string temp;
	time_t rawtime;
	struct tm * timeinfo;


	ifstream conFile;

	conFile.open("constants.h", ios::in);
	if (!conFile.is_open())
	{
		conFile.clear();
		conFile.open("code/constants.h", ios::in);
	}
	if (!conFile.is_open()) 		
	{	cout << "Unable to open constants.h for file header information. \n";
		exit(0);
	}

	while (getline(conFile, line))
	{
		while (line != "\0")
		{
			if (insideComment)
			{
				if (line.compare(0,2,"*/") == 0)
				{
					line.erase(0,2);
					insideComment=0;
				}
				else line.erase(0,1);

			}

			else
			{
				if (line.compare(0,2,"//") == 0)
				{
					//The line is a comment
					break;
				}

				else if (line.compare(0,2,"/*") == 0)
				{
					insideComment=1;
				}

				else if (line.compare(0,7,"#define") == 0)
				{
					line.erase(0,7);
				}

				else 
				{
					temp = line.substr(0,1);
					headerString.append(temp);
					line.erase (0,1);
				}
			}
		}

		//Now at end of line
		headerString.append("-");

	}

	//Remove white space
	headerString.erase(std::remove_if(headerString.begin(),headerString.end(),stripSpace),headerString.end());

	temp.clear();
	temp.append(headerString);

	headerString.clear();

	time ( &rawtime);
	timeinfo = localtime (&rawtime);

	headerString.append(asctime(timeinfo));

	headerString.erase(std::remove_if(headerString.begin(),headerString.end(),stripReturn),headerString.end());

	headerString.append(" ");
	headerString.append(temp);


	headerString.append("\n\n");

	return headerString;
	

}