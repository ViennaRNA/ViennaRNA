#ifndef ANCHORS_H_
#define ANCHORS_H_

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>
#include <stdarg.h>
#include <float.h>
#include "basepair.h"
#include "shape.hpp"


// global variables
 
extern int debug;
extern std::list<int> anchorlist;


// lex & parse

extern int yylex(void);
extern int yyparse(void);
extern void yyerror(char*);


// ccalc.c

extern void DumpRow(void); 
extern int GetNextChar(char *b, int maxBuffer);
extern void BeginToken(char*);
extern void PrintError(char *s, ...);
extern int parseStructure(std::string structure, std::list<int> * result);

#endif /*ANCHORS_H_*/

