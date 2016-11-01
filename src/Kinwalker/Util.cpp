#include "Util.h"

extern short *pair_table;
//extern void  *space(unsigned int size);
/*@exits@*/
//extern void   nrerror(const char message[]);

extern "C" {
#include "utils.h"
}

void Cout(std::string s){std::cout<<s;}

  std::string Str(double x) 
  {
    std::stringstream ss;
    ss << x;
    return ss.str();
    // return Str((float)x);
    /*
    std::stringstream ss;
    char cBuf[1000];
    sprintf(cBuf,"%12.3f",x);
    ss << cBuf;
    
    return ss.str();*/
  }


  std::string Str(float x) 
  {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }

  std::string Str(int x) 
  {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }

  std::string Str(unsigned int x) 
  {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }

std::string PrintPairTable(){
  std::string s=std::string();
  for(int i=0;i<=pair_table[0];i++){
    s+=Str(pair_table[i])+" ";
  }
  s+="\n";
  return s;
}


std::string PrintBasePair(std::pair<int,int> bp){
  return "("+Str(bp.first)+","+Str(bp.second)+")";
}


std::string PrintBasePairList(std::vector<std::pair<int,int> > list){
  std::string sText="";
  for(size_t i=0;i<list.size();i++) sText+=PrintBasePair(list[i]);
  return sText;
}


std::string
PairTableToStructure(std::vector<int> pTbl)
{
  std::string structure(pTbl.size(), '.');
  
  int i=0;
  for (std::vector<int>::const_iterator j = pTbl.begin();
       j != pTbl.end();
       j++, i++) {

    // we processed already the i < *j case
    if (i > *j) continue;

    // unpaired position
    if (*j == -1) {
      structure[*j] = '.';
      continue;
    }
    
    structure[i]  = '(';
    structure[*j] = ')';
  }

  return (structure);
}

 std::vector<std::pair<int,int> >
PairTableToBasePairList(short int * pair_table)
{
  std::vector<std::pair<int,int> > pl;

  //std::string structure(pair_table[0], '.');
  
  int i=1;
  for (int j = 1;j<= pair_table[0];j++, i++) {

    int val=pair_table[j] ;
    // we processed already the i < *j case
    if (i > val) continue;

    // unpaired position
    if (val == 0) {
      // structure[*j] = '.';
      continue;
    }
    pl.push_back(std::make_pair(i,val));
    //structure[i]  = '(';
    //structure[*j] = ')';
  }

  return pl;//(structure);
}


std::string
BasePairListToStructure(int length, std::vector<std::pair<int,int> > pl)
{
  std::string structure(length, '.');

  for (std::vector<std::pair<int,int> >::const_iterator p = pl.begin();
       p != pl.end();
       p++) {

    structure[(*p).first] = '(';
    structure[(*p).second] = ')';
  }

  return (structure);
}

std::string
BasePairListToStructure1(int length, const std::vector<std::pair<int,int> > & pl)
//BasePairListToStructure1(int length, std::vector<std::pair<int,int> >  pl)
{
  std::string structure(length, '.');

  //for(int i=0;i<pl.size();i++){
  //structure[pl[i].first-1] = '(';
  //structure[pl[i].second-1] = ')';
      for (std::vector<std::pair<int,int> >::const_iterator p = pl.begin();
       p != pl.end();
      p++) {
    structure[(*p).first-1] = '(';
    structure[(*p).second-1] = ')';
  }

  return (structure);
}

/*
void
StringToTextFile(std::string sText, std::string filePath)
{
  std::ofstream outfile(filePath.c_str());
  outfile << sText;
}
*/


 bool IntroducesPseudoKnot(const std::vector<std::pair<int,int> >& node, const std::pair<int,int>& p1){
       for(size_t i=0;i<node.size();i++){
          std::pair<int,int> p2=node[i];
          if(p1.first==p2.first || p1.first==p2.second || p1.second==p2.first || p1.second==p2.second) return true;
          if(p1.first<=p2.first && p1.second>= p2.first && p1.second<=p2.second) return true;
          if(p1.first>=p2.first && p1.first<=p2.second && p1.second>=p2.second)  return true;    
       }
       return false;
 }

bool ConformationHasPair(const std::vector<std::pair<int,int> >& node,const std::pair<int,int> & p){
  for(size_t i=0;i<node.size();i++){
    if(node[i].first==p.first && node[i].second==p.second) {
       return true;
    }
  }
  return false;
}


   bool Conflict(const std::vector<std::pair<int,int> >& node, const std::pair<int,int>& p1){
      if(ConformationHasPair(node,p1)) return true;
      if(IntroducesPseudoKnot(node,p1)) return true;
      return false;
   }      

 bool Conflict(const std::vector<std::pair<int,int> >& node1, const std::vector<std::pair<int,int> >& node2){
     for(size_t i=0;i<node2.size();i++){
       if(Conflict(node1,node2[i])) return true;
     }
     return false;

   }

//std::vector<std::vector<std::pair<int,int> > >
void
ConformationToStacks(std::vector<std::vector<std::pair<int,int> > > & stacks, std::vector<std::pair<int,int> > node,int stacksize)
{

  stacks.clear();


  // return if node is empty
  if (node.empty()) return;
  
  sort(node.begin(), node.end());

  //std::cout << "nodes" << std::endl
  //	    << PrintBasePairList(node) << std::endl;
  
  std::vector<std::pair<int,int> > stack;
  copy(node.begin(), node.begin()+1, std::inserter(stack, stack.begin()));

  //std::cout << "stack" << std::endl
  //	    << PrintBasePairList(stack) << std::endl;
  
  if (node.size()==1 && stacksize==1) {
    stacks.push_back(stack);
    return;
  }
  for(size_t i=1;i<node.size();i++){
    bool added=false;
    if(node[i].first==stack.back().first+1 &&
       node[i].second==stack.back().second-1) {
      added=true;
      stack.push_back(std::pair<int,int>(stack.back().first+1,
					 stack.back().second-1));
    }
    else if(node[i].first==stack.back().first+2 &&
	    node[i].second==stack.back().second-1) {
      added=true;
      stack.push_back(std::pair<int,int>(stack.back().first+2,
					 stack.back().second-1));
    }
    else if(node[i].first==stack.back().first+1 &&
	    node[i].second==stack.back().second-2) {
      added=true;
      stack.push_back(std::pair<int,int>(stack.back().first+1,
					 stack.back().second-2));
    }
    if(stack.size()>=(size_t)stacksize && (!added || i==node.size()-1) ) {
      stacks.push_back(stack);
      stack.clear();
      stack.push_back(node[i]);
    }
    if(!added) {
      stack.clear();
      stack.push_back(node[i]);
    }
  }
}


#ifdef WITH_DMALLOC
#define MG_space(S) calloc(1,(S))
#endif


 void MakePairTable(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   //   short *table;
   
   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   //table = (short *) space(sizeof(short)*(length+2));
   pair_table[0] = length;
   
   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(': 
	 stack[hx++]=i;
	 break;
       case ')':
	 j = stack[--hx];
	 if (hx<0) {
	   // fprintf(stderr, "%s\n", structure);
	    std::cout<<structure<<" unbalanced brackets in Michael's MakePairTable"<<std::endl;
	    //          nrerror("unbalanced brackets in Michael's MakePairTable");
	 }
	 pair_table[i]=j;
	 pair_table[j]=i;
	 break;
       default:   /* unpaired base, usually '.' */
	 pair_table[i]= 0;
	 break;
      }
   }
   if (hx!=0) {
     // fprintf(stderr, "%s\n", structure);
      std::cout<<structure<<" unbalanced brackets in Michael's MakePairTable"<<std::endl;
      //nrerror("unbalanced brackets in  Michael's MakePairTable");
   }
   free(stack);
   //return(table);
}

std::string CreateRGBString(double r,double g,double b) {
      return "["+Str(r)+" "+Str(g)+" "+Str(b)+"]";  
}



/**
  Code to create a .ps dotplot template. Based on PS_dot.c
*/

std::string PSFrontPlot(std::string sequence,std::vector<std::pair<int,int> > extrema){

  std::string sText;
       sText+=" %!PS-Adobe-3.0 EPSF-3.0\n";
       sText+=" %%Title: Maximum Matching Color Plot\n";
       sText+=" %%Creator: maxmatch.cpp\n";
       sText+=" %%BoundingBox: 66 211 518 662\n";
       sText+=" %%DocumentFonts: Helvetica\n";
       sText+=" %%Pages: 1\n";
       sText+=" %%EndComments\n";
       sText+=" \n";
       sText+=" %Options: \n";
       sText+=" % \n";
       sText+=" %Colored max matching matrixc.\n";
       sText+=" % i  j  sqrt(p(i,j)) ubox\n";
       sText+=" \n";
       sText+=" %%BeginProlog\n";
       sText+=" /DPdict 100 dict def\n";
       sText+=" DPdict begin\n";
       sText+=" /logscale false def\n";
       sText+=" \n";
       sText+=" \n";
       sText+=" \n";
       sText+=" /mylbox { % x y size [rgb] => -\n";
       sText+="    exch 4 2 roll\n";
       sText+="    len exch sub 1 add mybox\n";
       sText+=" } bind def\n";
       sText+=" \n";
       sText+=" /mybox { % [rgb] size x y box - draws box centered on x,y\n";
       sText+="    2 index 0.5 mul add            % x += 0.5\n";
       sText+="    exch 2 index 0.5 mul add exch  % x += 0.5\n";
       sText+="    newpath\n";
       sText+="    moveto\n";
       sText+="    dup neg   0 rlineto\n";
       sText+="    dup neg   0 exch rlineto\n";
       sText+="              0 rlineto\n";
       sText+="    closepath\n";
       sText+="    gsave\n";
       sText+="    aload pop\n";
       sText+="    setrgbcolor\n";
       sText+="    fill\n";
       sText+="    grestore\n";
       sText+=" } bind def\n";
       sText+=" \n";
       sText+=" \n";
       sText+=" \n";
       sText+=" \n";
       sText+=" /drawseq {\n";
       sText+=" % print sequence along all 4 sides\n";
       sText+=" [ \n";
       //       sText+=" [ [0.7 -0.3 0 ]\n"; remove bottom
       sText+="   [0.7 0.7 len add 0]\n";
       //       sText+="   [-0.3 len sub -0.4 -90]\n"; remove left
       sText+="   [-0.3 len sub 0.7 len add -90]\n";
       sText+=" ] {\n";
       sText+="    gsave\n";
       sText+="     aload pop rotate translate\n";
       sText+="     0 1 len 1 sub {\n";
       sText+="      dup 0 moveto\n";
       sText+="      sequence exch 1 getinterval\n";
       sText+="      show\n";
       sText+="     } for\n";
       sText+="    grestore\n";
       sText+="   } forall\n";
       sText+=" } bind def\n";
       sText+=" \n";
       sText+=" \n";
       sText+=" /drawgrid{\n";
       sText+="   0.01 setlinewidth\n";
       sText+="   len log 0.9 sub cvi 10 exch exp  % grid spacing\n";
       sText+="   dup 1 gt {\n";
       sText+="      dup dup 20 div dup 2 array astore exch 40 div setdash\n";
       sText+="   } { [0.3 0.7] 0.1 setdash } ifelse\n";
       sText+="   0 exch len {\n";
       sText+="      dup dup\n";
       sText+="      0 moveto\n";
       sText+="      len lineto \n";
       sText+="      dup\n";
       sText+="      len exch sub 0 exch moveto\n";
       sText+="      len exch len exch sub lineto\n";
       sText+="      stroke\n";
       sText+="   } for\n";
       sText+="   [] 0 setdash\n";
       sText+="   0.04 setlinewidth \n";
       sText+="   currentdict /cutpoint known {\n";
       sText+="     cutpoint 1 sub\n";
       sText+="     dup dup -1 moveto len 1 add lineto\n";
       sText+="     len exch sub dup\n";
       sText+="     -1 exch moveto len 1 add exch lineto\n";
       sText+="     stroke\n";
       sText+="   } if\n";
       sText+="   0.5 neg dup translate\n";
       sText+=" } bind def\n";
       sText+=" \n";
       sText+=" end\n";
       sText+=" %%EndProlog\n";
       sText+=" DPdict begin\n";
       sText+=" %delete next line to get rid of title\n";
       sText+=" %270 665 moveto /Helvetica findfont 14 scalefont setfont (dot.ps) show\n";
       sText+=" %set up the matrix\n";
       sText+=" /sequence { ("+sequence+") } def\n";
       sText+=" /len { sequence length } bind def\n";
       sText+=" 72 216 translate\n";
       sText+=" 72 6 mul len 1 add div dup scale\n";
       sText+=" /Helvetica findfont 0.95 scalefont setfont\n";
       sText+=" \n";
       sText+=" \n";
       //      sText+=" drawseq\n"; draw sequence after everything else has been drawn so it is not rendered invisible by colored squares
       sText+=" 0.5 dup translate\n";
       sText+=" % draw diagonal\n";
       sText+=" 0.04 setlinewidth\n";
       sText+=" 0 len moveto len 0 lineto stroke \n";

       sText+=" %dynamic content\n";
       sText+=" %drawgrid\n";
       sText+=" %data starts here\n";
       sText+=" % x 	y 	size 	rgb-array\n";

       //Add the actual entries and their colours for the upper half of the matrix
       for(std::vector<std::pair<int,int> >::iterator it=extrema.begin();it!=extrema.end();it++){
         for(int i=it->first; i<=it->second; i++) {
	    for (int j=i; j<=it->second; j++) {
 
              /*
              double color_depth=ratio(i,j)/max_ratio;
	      //round to two decimals
	      color_depth=.01*(int)floor(color_depth*100+.5);
	      */
	      //	      string rgb=CreateRGBString(1.0,1.0,1.0);
	      //	      string rgb=CreateRGBString(0.5,0.5,0.5);
	      //              string rgb=CreateRGBString(0.6,0.6,0.6);
	      std::string rgb=CreateRGBString(0.7,0.7,0.7);
              //string rgb=CreateRGBString(0.8,0.8,0.8);
	      sText+=Str(j)+"\t"+Str(i)+"\t"+Str(1)+"\t"+rgb+"\t"+"mylbox\n";
	      
	    }
	 }
       }
       sText+=" drawseq\n";
       //draw a colour scale in the lower half of the matrix 
       //       for(int i=0;i<10;i++) {
       //      string rgb=CreateRGBString(.1*i);
       //	      sText+=Str(i+1)+"\t"+Str(l-1)+"\t"+Str(1)+"\t"+rgb+"\t"+"mylbox\n";
       //}

       sText+=" showpage\n";
       return sText;

}
