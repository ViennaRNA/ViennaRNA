#include "anchors.h"

#define true 1
#define false 0

// global variables

int debug=0;
std::list<int> anchorlist;


// local variables

static const char *inputStruct;
static int eof = 0;
static int nRow = 0;
static int nBuffer = 0;
static int lBuffer = 0;
static int nTokenStart = 0;
static int nTokenLength = 0;
static int nTokenNextStart = 0;
static int lMaxBuffer = 1000;
static char *buffer;



extern
void PrintError(char *errorstring, ...) {
  static char errmsg[10000];
  va_list args;

  //int start=nTokenStart;
  //int end=start + nTokenLength - 1;
  //int i;

  // print it using variable arguments 
  va_start(args, errorstring);
  vsprintf(errmsg, errorstring, args);
  va_end(args);

  fprintf(stdout, "Error: %s\n", errmsg);
  exit(0);
}


// reads a line into the buffer

static
int getNextLine(void) {
  //int i;
  //char *p;

  nBuffer = 0;
  nTokenStart = -1;
  nTokenNextStart = 0; // was 1 once
  eof = false;

  // read a line 
  if (nRow>=1) { // we just read one line
    return 1;
  }
  strcpy(buffer, inputStruct);
  lMaxBuffer = strlen(inputStruct);

  nRow += 1;
  lBuffer = strlen(buffer);

  return 0;
}

// reads a character from input for flex
extern
int GetNextChar(char *b, int maxBuffer) {
  int frc;
  
  if (  eof  )
    return 0;
  
  // read next line if at the end of the current
  while (  nBuffer >= lBuffer  ) {
    frc = getNextLine();
    if (  frc != 0  )
      return 0;
    }

  // ok, return character
  b[0] = buffer[nBuffer];
  nBuffer += 1;

  return b[0]==0?0:1;
}


// marks the beginning of a new token

extern
void BeginToken(char *t) {
  // remember last read token
  nTokenStart = nTokenNextStart;
  nTokenLength = strlen(t);
  nTokenNextStart = nBuffer; // + 1;

  // location for bison 
  yylloc.first_line = nRow;
  yylloc.first_column = nTokenStart;
  yylloc.last_line = nRow;
  yylloc.last_column = nTokenStart + nTokenLength - 1;

}

extern
int parseStructure(std::string structure, std::list<int> * result) {
  inputStruct = structure.c_str();
  debug = 0;
  result->clear();
  anchorlist.clear();
  eof = 0;
  nRow = 0;
  nBuffer = 0;
  lBuffer = 0;
  nTokenStart = 0;
  nTokenLength = 0;
  nTokenNextStart = 0;
  lMaxBuffer = 1000;

  if (  inputStruct == NULL  ) {
    std::cerr << "No structure as input." << std::endl;
    exit(0);
  }

  buffer = (char*) malloc(lMaxBuffer);
  if (  buffer == NULL  ) {
    std::cerr << "Cannot allocate " <<lMaxBuffer << " bytes of memory\n";
    return 12;
  }

  // parse it
  if (  getNextLine() == 0  )
    yyparse();

  * result = anchorlist;
  free(buffer);

  return 0;
}
