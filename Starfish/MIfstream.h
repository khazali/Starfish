#ifndef MIFSTREAM_H
#define MIFSTREAM_H
#include <fstream>

class MIfstream: public std::ifstream {
public:
	int FileSearch(char *);
	int ReadWord(char *);
};

#endif