#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "MIfstream.h"
#include "Globals.h"
#include <cstdio>

#define MAX_QUOTES_LENGTH 2000

void PrintQuotes(void) {
	char str[MAX_QUOTES_LENGTH];
	MIfstream QuotesFile;
	unsigned int QuotesNO, QuotesSelection, i;

	QuotesFile.open("Quotes.txt", std::ifstream::in);
	QuotesFile.ReadWord(str);
	QuotesNO=atoi(str);

	QuotesSelection=1+(rand()%QuotesNO);

	for (i=0; i<QuotesSelection; i++) {
		QuotesFile.getline(str, MAX_QUOTES_LENGTH);
	}
	QuotesFile.close();
	std::cout << std::endl << std::endl << str << std::endl << std::endl;
}