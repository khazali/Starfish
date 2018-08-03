#include <fstream>
#include <string.h>
#include "MIfstream.h"
#include "Globals.h"

int MIfstream::FileSearch(char *rSeek) {
	register int i;
	char str[MAX_STRING_LENGTH];

	clear();                 // clear fail and eof bits
	seekg(0, std::ios::beg); // back to the start!
	*str='\0';
	do {
		i=ReadWord(str);
		if (!strcmp(str, rSeek)) return -1;
	} while (i);

	return 0;		//Nothing found
}

int MIfstream::ReadWord(char *rWord) {
	char ch;
	register int i=0;

	*rWord='\0';
	do {
		get(ch);
		if (eof()) {
			return 0;		//Nothing has been read
		}
	} while ((ch<33) || (ch>126));

	while ((ch>32) && (ch<127)) {
		*(rWord+i)=ch;
		i++;
		get(ch);
		if (eof()) {
			*(rWord+i)='\0';
			return 1;		//Read, but end of file also encountered
		}		
	}

	*(rWord+i)='\0';
	return -1;		//Correct execution
}