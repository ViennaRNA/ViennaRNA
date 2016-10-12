#ifndef WIN32
#include "config.h"
#endif

#include <cmath>

#ifdef HAVE_LIBG2
#include <g2.h>
#include <g2_PS.h>
#endif

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>

#ifdef HAVE_LIBXMLPLUSPLUS
#include <libxml++/libxml++.h>
#endif

#include "matrix.h"
#include "rna_profile_alignment.h"
#include "rnafuncs.h"

#ifdef HAVE_LIBRNA

extern "C" {
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/naview.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
}

#endif

#include "forest.t.cpp"

std::ostream& operator <<(std::ostream &s,RNA_Alphabet_Profile p) {

    s.precision(2);
    s << "A " << p.p[ALPHA_PRO_BASE_A] << " ";
    s.precision(2);
    s << "C " << p.p[ALPHA_PRO_BASE_C] << " ";
    s.precision(2);
    s << "G " << p.p[ALPHA_PRO_BASE_G] << " ";
    s.precision(2);
    s << "U " << p.p[ALPHA_PRO_BASE_U] << " ";
    s.precision(2);
    s << "P " << p.p[ALPHA_PRO_BASEPAIR] << " ";
    s.precision(2);
    s << p.p[ALPHA_PRO_GAP] << " ";
    s.precision(2);
    s << p.p[ALPHA_PRO_BASE] << " ";
    s.precision(2);
    s << p.p[ALPHA_PRO_ANCHOR] << std::endl;

    //  for(unsigned int i=0;i<p.columnStr.length();i++)
    //    s << p.columnStr[i] << "\\n";

    return s;
}

/* ****************************************** */
/*    Constructor and Destructor functions    */
/* ****************************************** */

RNAProfileAlignment::RNAProfileAlignment(const std::string &baseStr, const std::string &viennaStr, const std::string &name, bool anchored)
        : ForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>(RNAFuncs::treeSize(viennaStr, false)),// just a forest of corr. size 
        name_(name),
        numStructures_(1) {
		if (anchored)
			buildForestAnchored(baseStr,viennaStr);
		else
			buildForest(baseStr,viennaStr);
    addStrName(name);
};

#ifdef HAVE_LIBRNA
RNAProfileAlignment::RNAProfileAlignment(const std::string &baseStr, const std::string &name, bool anchored, const std::string &constraint, double t)
        : ForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>(2*baseStr.length()),
        name_(name),
        numStructures_(1) {
    char *viennaStr = NULL;

    // calculate partition function for the sequence
    //do_backtrace=1;
    init_pf_fold(baseStr.length());

    //if(constraint.length()>0)
    //pf_fold((char*)baseStr.c_str(),(char*)constraint.c_str());    // expicit conversion to non-const value, but pf_fold does not alter baseStr
    //else
    pf_fold((char*)baseStr.c_str(),NULL);    // expicit conversion to non-const value, but pf_fold does not alter baseStr

    viennaStr = new char[baseStr.length()+1];
    //dangles=2;
    fold((char*)baseStr.c_str(),viennaStr);

    setSize(RNAFuncs::treeSize(std::string(viennaStr),anchored));
		bool useBPprob = true;
		if (anchored)
			buildForestAnchored(baseStr,viennaStr,useBPprob);
		else
			buildForest(baseStr,viennaStr,useBPprob);

    free_pf_arrays();
    delete[] viennaStr;

    //  hasSequence=true;
    addStrName(name);
}
#endif

RNAProfileAlignment::RNAProfileAlignment(const std::string &filename, bool anchored) {
    // read forest from file
    size_type size;
    std::ifstream s;
    char *colStr;

    s.open(filename.c_str());

    // first of all save the size
    s.read((char*)&size,sizeof(size_type));
    s.read((char*)&numStructures_,sizeof(unsigned int));

    colStr=new char[numStructures_+1];
    colStr[numStructures_]=0;   // terminate string

    // initialize structures
    //  ::ForestBase(size_);
    //  lb_=new L[f.size()];
    initialize(size,anchored);

    // save the arrays
    for (unsigned int i=0; i<size; i++) {
        label_type l;

        for (int r=0; r<RNA_ALPHABET_SIZE; r++) {
            double d;
            s.read(reinterpret_cast<char*>(&d),sizeof(double));
            l.p[r]=d;
        }

        s.read(colStr,sizeof(char)*numStructures_);

        l.columnStr=colStr;
        lb_[i]=l;
    }

    s.read(reinterpret_cast<char*>(rb_),sizeof(size_type)*size);
    s.read(reinterpret_cast<char*>(noc_),sizeof(size_type)*size);
    //  s.read((char*)sumUpCSF_,sizeof(size_type)*size);
    //  s.read((char*)rmb_,sizeof(size_type)*size);

    for (unsigned int r=0; r<numStructures_; r++) {
        char str[21];
        str[20]=0;

        s.read(str,sizeof(char)*20);
        strNames_.push_back(std::string(str));
    }


    calcSumUpCSF();
    calcRMB();
}

/*
RNAProfileForest::RNAProfileForest(const Profile_RNA_Alignment_Forest &f)
  : RNAForestBase<RNA_Alphabet_Profile>(f)
{
  numStructures_=f.numStructuresX_+f.numStructuresY_;
}*/

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

void RNAProfileAlignment::makeRepLabel(size_type node, RNA_Alphabet_Profile a,  RNA_Alphabet_Profile b, int anchor) {
    double p,q;
    double m;

    p=numStructuresX_;	// weight number of structures in profile
    q=numStructuresY_;
    m=p+q;
    p/=m;
    q/=m;

    // profile
    for (int i=0; i<RNA_ALPHABET_SIZE-1; i++)
        lb_[node].p[i]=p*a.p[i]+q*b.p[i];

		lb_[node].p[ALPHA_PRO_ANCHOR]=anchor;
		if (anchor !=0) {
		//	std::cout << "setting anchor " << anchor << " in join" << std::endl;
			setAnchor(node,anchor);
		}
    // alignment column
    lb_[node].columnStr=a.columnStr + b.columnStr;
};

void RNAProfileAlignment::makeDelLabel(size_type node, RNA_Alphabet_Profile a) {
    double p,q;
    double m;

    p=numStructuresX_;	// weight number of structures in profile
    q=1;
    m=p+q;
    p/=m;
    q/=m;

    // profile
    lb_[node]=a;
    //lb_[node].p[ALPHA_GAP]=(1+lb_[node].p[ALPHA_GAP])/2.0;

    for (int i=0; i<RNA_ALPHABET_SIZE; i++) {
        lb_[node].p[i]*=p;
    }
    lb_[node].p[ALPHA_PRO_GAP]+=q;

    // alignment column
    lb_[node].columnStr=a.columnStr + std::string(numStructuresY_,ALPHA_GAP);
};

void RNAProfileAlignment::makeInsLabel(size_type node, RNA_Alphabet_Profile b) {
    double p,q;
    double m;

    p=1;						// weight number of structures in profile
    q=numStructuresY_;
    m=p+q;
    p/=m;
    q/=m;

    // profile
    lb_[node]=b;
    //lb_[node].p[ALPHA_GAP]=(1+lb_[node].p[ALPHA_GAP])/2.0;

    for (int i=0; i<RNA_ALPHABET_SIZE; i++) {
        lb_[node].p[i]*=q;
    }
    lb_[node].p[ALPHA_PRO_GAP]+=p;

    //	  lb_[node].p[ALPHA_PRO_GAP]=(p+q*lb_[node].p[ALPHA_PRO_GAP])/2.0;

    // alignment column
    lb_[node].columnStr=std::string(numStructuresX_,ALPHA_GAP) + b.columnStr;
};

void RNAProfileAlignment::showLabel(std::ostream &s,RNA_Alphabet_Profile p) const {
    s << p;
};


void RNAProfileAlignment::makeLabel(RNA_Alphabet_Profile &p,char c) {
    int i;

    // initialize profile entries to zero - but not the anchors, which are set separately
    for (i=0; i<RNA_ALPHABET_SIZE-1; i++)
        p.p[i]=0.0;

    i=alpha2RNA_Alpha(c);
    p.p[i]=1.0;

    // if it is acgu it is also a base
    if (i<=ALPHA_PRO_BASE_U)
        p.p[ALPHA_PRO_BASE]=1.0;

    // set column string
    p.columnStr=c;
}


void RNAProfileAlignment::buildForestAnchored(const std::string &base, const std::string &anchored_vienna, bool use_bp_prob) {

		unsigned int baseStrLen = base.length();
    assert(base.length() <= anchored_vienna.length());

    unsigned long length = 0, basePairCount = 0, maxDepth = 0;
    RNAFuncs::isAnchoredViennaString(anchored_vienna, length, basePairCount, maxDepth);

    assert(length <= base.length());

    // TODO richtig via oberklassenkonstruktor machen
    unsigned long nodes = base.length()+basePairCount;
    setSize(nodes);

    rb_ = new size_type[nodes];
    noc_ = new size_type[nodes];
    sumUpCSF_ = new size_type[nodes];
    rmb_ = new size_type[nodes];
    anchors_ = new size_type[nodes];
		lb_=new RNA_Alphabet_Profile[nodes];
    memset(rb_,0,sizeof(size_type)*nodes);
    memset(noc_,0,sizeof(size_type)*nodes);
    memset(sumUpCSF_,0,sizeof(size_type)*nodes);
    memset(rmb_,0,sizeof(size_type)*nodes);
    memset(anchors_,0,sizeof(size_type)*nodes);
    //memset(lb_,0,sizeof(RNA_Alphabet_Profile)*nodes);

    unsigned int * nodeStack = new unsigned int[maxDepth+1];
    unsigned int * numChildrenStack = new unsigned int[maxDepth+1];
    unsigned int * baseposStack = new unsigned int[maxDepth+1];
    memset(nodeStack,0,sizeof(unsigned int)*maxDepth+1);
    memset(numChildrenStack,0,sizeof(unsigned int)*maxDepth+1);
    memset(baseposStack,0,sizeof(unsigned int)*maxDepth+1);

    // keep track of the number of built nodes
    unsigned long node = 0;
    unsigned long stackPtr = 0;

    //std::cout << "reading anchored_vienna string " << anchored_vienna << std::endl;
    int anchor = 1;

    // fill Forest structure
    for (unsigned int i = 0, j = 0; i< base.length() && j < anchored_vienna.length();j++) {
        assert(node <= nodes);
        switch (anchored_vienna[j]) {
            case 'a': {
                // record read anchor
                setAnchor(node,anchor);
								lb_[node].p[ALPHA_PRO_ANCHOR]=anchor;
								//if (anchor !=0) {
								//   std::cout << "setting anchor " << anchor << "in build" << std::endl;
								//}
								anchor++;

                //std::cout << "node " << node << " got anchor " << getAnchor(node) << std::endl;
                break;
            } 
						case '.': {
		            // set label
		            if (baseStrLen)
		                makeLabel(lb_[node],base[i++]);
		            else
		                makeLabel(lb_[node],'B');
		
		            // set right brother
		            if (node==size()-1)
		                setRightBrotherIndex(node,0);
		            else
		                setRightBrotherIndex(node,node+1);
		
		            // set num children
		            setNumChildren(node,0);
		
		            // increase stack values
		            numChildrenStack[stackPtr]++;
		
		            node++;
		            break;
						}
		        case '(':
		            // set label
		            makeLabel(lb_[node],'P');
		
		            // increase stack values
		            numChildrenStack[stackPtr]++;
		            baseposStack[stackPtr]=i+1;
		
		            // push
		            stackPtr++;
		            nodeStack[stackPtr]=node;
		            numChildrenStack[stackPtr]=1;
		
		            node++;
		
		            // set label
		            if (baseStrLen)
		                makeLabel(lb_[node],base[i++]);
		            else
		                makeLabel(lb_[node],'B');
		
		            // set right brother
		            setRightBrotherIndex(node,node+1);
		
		            // set num children
		            setNumChildren(node,0);
		
		            node++;
		            break;
		        case ')':
		            // set label
		            if (baseStrLen)
		                makeLabel(lb_[node],base[i++]);
		            else
		                makeLabel(lb_[node],'B');
		
		            // set right brother
		            setRightBrotherIndex(node,0);
		
		            // set num children
		            setNumChildren(node,0);
		
		            // pop
		            if (node==size()-1)
		                setRightBrotherIndex(nodeStack[stackPtr],0);
		            else
		                setRightBrotherIndex(nodeStack[stackPtr],node+1);
		
		            setNumChildren(nodeStack[stackPtr],numChildrenStack[stackPtr]+1);
		            stackPtr--;

#ifdef HAVE_LIBRNA
            // set basepair probability
            if (use_bp_prob) {
                lb_[nodeStack[stackPtr]].p[ALPHA_PRO_BASEPAIR]=pr[iindx[baseposStack[stackPtr]]-(i+1)];
                lb_[nodeStack[stackPtr]].p[ALPHA_PRO_GAP]=1-lb_[nodeStack[stackPtr]].p[ALPHA_PRO_BASEPAIR];
            }
#endif
            node++;
            break;
        }
    }

    delete[] nodeStack;
    delete[] numChildrenStack;
    delete[] baseposStack;
    calcSumUpCSF();
    calcRMB();
    assert (noc_[size_-1]==0);
}


void RNAProfileAlignment::buildForest(const std::string &baseStr, const std::string &viennaStr, bool use_bp_prob) {

    assert(RNAFuncs::isRNAString(baseStr));

    unsigned long basePairCount = 0, maxDepth = 0;
    RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);

    unsigned int baseStrLen = baseStr.length();
    unsigned int viennaStrLen = viennaStr.length();

		// TODO unneeded
    if (baseStrLen)
        hasSequence_ = true;
    else
        hasSequence_ = false;

    // check if base string and vienna string have the same length
    //  if(baseStr.length() > 0)
    //    if(baseStr.length() != viennaStrLen)
    //      throw RNAForestExceptionInput(RNAForestExceptionInput::Error_BaseStringAndViennaStringIncompatible);


    unsigned int * nodeStack = new unsigned int[maxDepth+1];
    unsigned int * numChildrenStack = new unsigned int[maxDepth+1];
    unsigned int * baseposStack = new unsigned int[maxDepth+1];
    memset(nodeStack,0,sizeof(unsigned int)*maxDepth+1);
    memset(numChildrenStack,0,sizeof(unsigned int)*maxDepth+1);
    memset(baseposStack,0,sizeof(unsigned int)*maxDepth+1);

    // keep track of the number of built nodes
    unsigned long node = 0;
    unsigned long stackPtr = 0;

    // fill Forest structure
    for (unsigned int i=0; i<viennaStrLen; i++) {
        switch (viennaStr[i]) {
        case '.':
            // set label
            if (baseStrLen)
                makeLabel(lb_[node],baseStr[i]);
            else
                makeLabel(lb_[node],'B');

            // set right brother
            if (node==size()-1)
                setRightBrotherIndex(node,0);
            else
                setRightBrotherIndex(node,node+1);

            // set num children
            setNumChildren(node,0);

            // increase stack values
            numChildrenStack[stackPtr]++;

            node++;
            break;

        case '(':
            // set label
            makeLabel(lb_[node],'P');

            // increase stack values
            numChildrenStack[stackPtr]++;
            baseposStack[stackPtr]=i+1;

            // push
            stackPtr++;
            nodeStack[stackPtr]=node;
            numChildrenStack[stackPtr]=1;

            node++;

            // set label
            if (baseStrLen)
                makeLabel(lb_[node],baseStr[i]);
            else
                makeLabel(lb_[node],'B');

            // set right brother
            setRightBrotherIndex(node,node+1);

            // set num children
            setNumChildren(node,0);

            node++;
            break;
        case ')':
            // set label
            if (baseStrLen)
                makeLabel(lb_[node],baseStr[i]);
            else
                makeLabel(lb_[node],'B');

            // set right brother
            setRightBrotherIndex(node,0);

            // set num children
            setNumChildren(node,0);

            // pop
            if (node==size()-1)
                setRightBrotherIndex(nodeStack[stackPtr],0);
            else
                setRightBrotherIndex(nodeStack[stackPtr],node+1);

            setNumChildren(nodeStack[stackPtr],numChildrenStack[stackPtr]+1);
            stackPtr--;

#ifdef HAVE_LIBRNA
            // set basepair probability
            if (use_bp_prob) {
                lb_[nodeStack[stackPtr]].p[ALPHA_PRO_BASEPAIR]=pr[iindx[baseposStack[stackPtr]]-(i+1)];
                lb_[nodeStack[stackPtr]].p[ALPHA_PRO_GAP]=1-lb_[nodeStack[stackPtr]].p[ALPHA_PRO_BASEPAIR];
            }
#endif
            node++;
            break;
        }
    }
    delete[] nodeStack;
    delete[] numChildrenStack;
    delete[] baseposStack;
    calcSumUpCSF();
    calcRMB();
    assert (noc_[size_-1]==0);
}

void RNAProfileAlignment::getStructureAlignmentFromCSF(std::string &s, std::vector<double> &pairprob, double t,size_type i, size_type j) const {
    size_type h;
    double bestPairScore;
    //  vector<unsigned int> leftPairList;
    //vector<unsigned int> rightPairList;
    //unsigned int bestLeftIndex,bestRightIndex;
    //unsigned int lastLeftIndex,lastRightIndex;

    //  QWATCH(i);
    //  QWATCH(j);

    if (j==0)
        return;

    if (isPair(i)) {
        // TRACE(DBG_GET_PROFILE_STRUCTURE,"Profile_RNA_Alignment_Forest::getStructureAlignmentFromCSF","basepair");
        // WATCH(DBG_GET_PROFILE_STRUCTURE,"Profile_RNA_Alignment_Forest::getStructureAlignmentFromCSF",lb_[i].p[ALPHA_BASEPAIR]);

        bestPairScore=bestPairs(i);

        // backtrace best pairs
        if (bestPairScore == bestPairs(i+1) + bestPairs(getRightmostBrotherIndex(i+1)) || lb_[i].p[ALPHA_PRO_BASEPAIR] < t) {
            // i pairs not
            //	std::cout << "unpaired" << std::endl;
            getStructureAlignmentFromCSF(s,pairprob,t,i+1,noc(i));
            //	std::cout << "back to:" << std::endl;
            //	QWATCH(i);
            //	QWATCH(j);
        } else {
            //	std::cout << "paired" << std::endl;
            // i pairs
            s += '(';
            pairprob.push_back(lb_[i].p[ALPHA_PRO_BASEPAIR]);

            // left path - righthand best pairs
            h=i+1;
            while (h < size() && isPair(h)) {
                //	  std::cout << "left" << std::endl;
                //	  QWATCH(h);

                assert((int)noc(h)-1>=0);
                getStructureAlignmentFromCSF(s,pairprob,t,rb(h+1),noc(h)-1);
                h=h+1;
            }


            assert((int)noc(i)-2>=0);
            getStructureAlignmentFromCSF(s,pairprob,t,rb(i+1),noc(i)-2);
            //	std::cout << "back to:" << std::endl;
            //	QWATCH(i);
            //	QWATCH(j);

            // right path - lefthand best pairs
            h=getRightmostBrotherIndex(i+1);
            while (h < size() && isPair(h)) {
                //	  std::cout << "right" << std::endl;
                //	  QWATCH(h);

                assert((int)noc(h)-1>=0);
                getStructureAlignmentFromCSF(s,pairprob,t,h+1,noc(h)-1);
                //h=h+1;
                h=getRightmostBrotherIndex(h+1);
            }


            s += ')';
            pairprob.push_back(lb_[i].p[ALPHA_PRO_BASEPAIR]);
        }
    } else {
        s+= '.';
        pairprob.push_back(lb_[i].p[ALPHA_PRO_BASE]);
    }

    // right forest
    getStructureAlignmentFromCSF(s,pairprob,t,rb(i),j-1);
}

double RNAProfileAlignment::bestPairs(size_type node) const {
    size_type i=node;
    double d1,d2;

    WATCH(DBG_GET_PROFILE_STRUCTURE,"RNAProfileForest::getStructureAlignment",node);

    if (isBase(node))
        return 0;
    else {
        // node pairs not
        d1=bestPairs(node+1) + bestPairs(getRightmostBrotherIndex(node+1));

        // node pairs
        d2=label(node).p[ALPHA_PRO_BASEPAIR];

        // left path - righthand best pairs
        i=node+1;
        while (i < size() && isPair(i)) {
            d2+=bestPairs(getRightmostBrotherIndex(i+1));
            i=i+1;
        }

        // right path - lefthand best pairs
        i=getRightmostBrotherIndex(node+1);
        while (isPair(i) && i < size()) {
            d2+=bestPairs(i+1);
            i=getRightmostBrotherIndex(i+1);
        }

        return std::max(d1,d2);
    }
}

#ifdef HAVE_LIBG2
void RNAProfileAlignment::drawBaseCircles(int device_id,const BaseProbs &bp,double center_x,double center_y) const {
    const double box_size=6.0;
    const double max_radius=3.0;

    double xpos,ypos;
    int color;
    //  double dashes=0.5;

    // draw the base probabilities as circles on the edges of the square and the gap probability
    // as the center square

    // upper left corner = a
    xpos=center_x-box_size/2.0;
    ypos=center_y+box_size/2.0;

    color=g2_ink(device_id,1,0,0);
    g2_pen(device_id,color);
    g2_filled_circle(device_id,xpos,ypos,max_radius*bp.a);

    // upper right corner = c
    xpos=center_x+box_size/2.0;
    ypos=center_y+box_size/2.0;

    color=g2_ink(device_id,0,1,0);
    g2_pen(device_id,color);
    g2_filled_circle(device_id,xpos,ypos,max_radius*bp.c);

    // lower right corner = g
    xpos=center_x+box_size/2.0;
    ypos=center_y-box_size/2.0;

    color=g2_ink(device_id,0,0,1);
    g2_pen(device_id,color);
    g2_filled_circle(device_id,xpos,ypos,max_radius*bp.g);

    // lower right corner = u
    xpos=center_x-box_size/2.0;
    ypos=center_y-box_size/2.0;

    color=g2_ink(device_id,1,0,1);
    g2_pen(device_id,color);
    g2_filled_circle(device_id,xpos,ypos,max_radius*bp.u);

    // gap in the center
    color=g2_ink(device_id,0,0,0);
    g2_pen(device_id,color);
    g2_filled_circle(device_id,center_x,center_y,max_radius*bp.gap);

    // draw rectangle for orientation
    //g2_set_dash(device_id,1,&dashes);
    color=g2_ink(device_id,0.75,0.75,0.75);
    g2_pen(device_id,color);
    g2_rectangle(device_id,center_x-box_size/2,center_y-box_size/2,center_x+box_size/2,center_y+box_size/2);
}
#endif

inline double RNAProfileAlignment::getMlBaseFreq(const BaseProbs &bp) const {
    return std::max(std::max(std::max(bp.a,bp.c),bp.g),bp.u);
}

char RNAProfileAlignment::getMlBase(const BaseProbs &bp) const {
    double p = getMlBaseFreq(bp);
    if (bp.a==p)
        return ALPHA_BASE_A;
    if (bp.c==p)
        return ALPHA_BASE_C;
    if (bp.g==p)
        return ALPHA_BASE_G;
    if (bp.u==p)
        return ALPHA_BASE_U;

    assert(false); // you should never get so far ...
    return 'x';
}

void RNAProfileAlignment::getSeqAli(std::string &seq,unsigned int row, unsigned int i, unsigned int j) const {
    if (i==0 && j==getMaxLength(0))
        seq="";

    if (j==0)
        return;

    // basepair=internal node
    if (isPair(i)) {
        getSeqAli(seq,row,i+1,noc(i));
        getSeqAli(seq,row,rb(i),j-1);
    } else {
        // base=leaf
        seq+=lb_[i].columnStr[row];
        getSeqAli(seq,row,rb(i),j-1);
    }
}

void RNAProfileAlignment::getStructAli(std::string &s,unsigned int row) const {
    s="";

    std::map<unsigned int,unsigned int> pairs;
    makePairTable(pairs, row);

    // iterate through leaves nodes and use information of pairs
    for (unsigned int i=0; i<size_; i++) {
        RNA_Alphabet c;

        c=lb_[i].columnStr[row];

        if (isBase(i)) {
            if (c != '-') {
                if (pairs.find(i) != pairs.end()) {	// is base paired ?
                    unsigned int j=pairs[i];
                    if (i<j)
                        s+='(';
                    else
                        s+=')';
                } else
                    s+='.';
            } else
                s+='-';
        }
    }
}

void RNAProfileAlignment::makePairTable(std::map<unsigned int,unsigned int> &pairs, unsigned int row) const {
    std::pair<int,int> *baseIndex;		// left and right base
    RNA_Alphabet c;

    assert(row<numStructures_);
    baseIndex=new std::pair<int,int>[size_];

    // initialize pairs, all gaps
    for (int i=size()-1; i>=0; i--) {
        baseIndex[i].first=-1;
        baseIndex[i].second=-1;
    }

    for (int i=size()-1; i>=0; i--) {
        c=lb_[i].columnStr[row];

        if (isLeaf(i)) {
            if (c!=ALPHA_GAP) {
                baseIndex[i].first=i;
                baseIndex[i].second=i;
            }
        } else {
            // internal node
            // leftmost and rightmost base
            bool lmBaseFound=false;
            for (size_type r=0,h=i+1; r<getMaxLength(i+1); r++,h=rb(h)) {
                // leftmost base
                if (!lmBaseFound && baseIndex[h].first != -1) {
                    baseIndex[i].first=baseIndex[h].first;
                    lmBaseFound=true;
                }

                // rightmost base
                if (baseIndex[h].second != -1) {
                    baseIndex[i].second=baseIndex[h].second;
                }
            }

            // report pairing bases if P node
            if (c==ALPHA_BASEPAIR) {
                assert(baseIndex[i].first != -1);
                assert(baseIndex[i].second != -1);
                assert(baseIndex[i].first < baseIndex[i].second);

                pairs[baseIndex[i].first]=baseIndex[i].second;
                pairs[baseIndex[i].second]=baseIndex[i].first;
            }
        }
    }

    delete[] baseIndex;
}

void RNAProfileAlignment::filterConsensus(std::string &structure, std::vector<double> &pairprob, std::vector<BaseProbs> &baseprobs, double minFreq) const {
    std::vector<double>::iterator itPP;
    //	std::vector<BaseProbs>::iterator itBP;
    std::string::iterator it;

    for (itPP=pairprob.begin(),it=structure.begin(); itPP!=pairprob.end(); itPP++,it++) {
        if (*itPP<minFreq) {
            *itPP=0.0;
            *it='.';
        }
    }

    /*
    for(it=structure.begin(),itPP=pairprob.begin(),itBP=baseprobs.begin();it!=structure.end();)
    {
    	if(itBP->base<minFreq)
    	{
    		structure.erase(it);
    		baseprobs.erase(itBP);
    		pairprob.erase(itPP);

    		it=structure.begin();
    		itPP=pairprob.begin();
    		itBP=baseprobs.begin();
    	}
    	else
    	{
    		it++;
    		itPP++;
    		itBP++;
    	}
    }
    */
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

void RNAProfileAlignment::getSequenceAlignment(std::vector<BaseProbs> &baseprob) const {
    BaseProbs bp;

    // generate base strings
    for (size_type i=0; i<size(); i++) {
        if (isBase(i)) {
            bp.a=label(i).p[ALPHA_PRO_BASE_A];
            bp.c=label(i).p[ALPHA_PRO_BASE_C];
            bp.g=label(i).p[ALPHA_PRO_BASE_G];
            bp.u=label(i).p[ALPHA_PRO_BASE_U];
            bp.gap=label(i).p[ALPHA_PRO_GAP];
            bp.base=bp.a+bp.c+bp.g+bp.u;

            baseprob.push_back(bp);
        }
    }
}

void RNAProfileAlignment::getStructureAlignment(double t,std::string &s, std::vector<double> &pairprob) const {
    WATCH(DBG_GET_PROFILE_STRUCTURE,"Profile_RNA_Alignment_Forest::getStructureAlignment",size());

    s="";
    getStructureAlignmentFromCSF(s,pairprob,t,0,getMaxLength(0));
}

#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
void RNAProfileAlignment::squigglePlot(const std::string &filename, SquigglePlotOptions &options) const {
    const double base_fontsize=8;
    const unsigned int nugrey_colors_=100;
    const double min_grey_color=1.0;

    std::string seq,structure;
    std::string base,structname;
    float *X,*Y,min_X=0,max_X=0,min_Y=0,max_Y=0;
    unsigned int i;
    short *pair_table;
    int id_PS,id;
    int ps_grey_colors[nugrey_colors_];
    int ps_color_red;
    int ps_color_black;
    double xpos,ypos;

    std::vector<double> pairprob;
    std::vector<BaseProbs> baseprobs;

    getStructureAlignment(options.minPairProb,structure,pairprob);
    getSequenceAlignment(baseprobs);

    //  filterConsensus(structure,pairprob,baseprobs,0.5);

    //assert(baseprobs.size() == structure.size());
    if (baseprobs.size() != structure.size())
        std::cerr <<  "Error in resolving consensus structure!" << std::endl;

    X = new float[structure.size()];
    Y = new float[structure.size()];

    pair_table = make_pair_table(structure.c_str());
    i = naview_xy_coordinates(pair_table, X, Y);
    if (i!=structure.size())
        std::cerr << "strange things happening in squigglePlot ..." << std::endl;

    // calculate image dimesions
    for (i=0; i<structure.size(); i++) {
        min_X=std::min(min_X,X[i]);
        max_X=std::max(max_X,X[i]);
        min_Y=std::min(min_Y,Y[i]);
        max_Y=std::max(max_Y,Y[i]);
    }

    //  id_PS  = g2_open_PS("ali.ps", g2_A4, g2_PS_port);
    id_PS  = g2_open_EPSF((char*)filename.c_str());
    id     = g2_open_vd();
    g2_attach(id,id_PS);

    //  std::cout << "min_X: " << min_X <<",max_X: " << max_X << ",min_Y: " << min_Y << "max_Y: " << max_Y << std::endl;
    g2_set_coordinate_system(id_PS,595/2.0,842/2.0,0.5,0.5);
    g2_set_line_width(id,0.2);


    // set colors
    double intv=min_grey_color/(double)nugrey_colors_;
    for (i=0; i<nugrey_colors_; i++) {
        double grey_color=min_grey_color-i*intv;
        ps_grey_colors[i]=g2_ink(id_PS,grey_color,grey_color,grey_color);
    }

    ps_color_black=g2_ink(id_PS,0,0,0);
    if (options.greyColors)
        ps_color_red=g2_ink(id_PS,0,0,0);
    else
        ps_color_red=g2_ink(id_PS,1,0,0);

    // draw sequence
    g2_set_font_size(id,base_fontsize);
    for (i=0; i<structure.size(); i++) {

        if (options.mostLikelySequence) {
            double p=getMlBaseFreq(baseprobs[i]);

            //base color
            if (p==1)
                g2_pen(id,ps_color_red);
            else
                g2_pen(id,ps_grey_colors[(int)floor(p*nugrey_colors_-1)]);

            base=getMlBase(baseprobs[i]);

            xpos=X[i]-base.length()*base_fontsize/2.0;
            ypos=Y[i]-4;
            g2_string(id,xpos,ypos,(char*)base.c_str());
        } else {
            drawBaseCircles(id_PS,baseprobs[i],X[i],Y[i]);
        }

        // connection to next base
        if (i<structure.size()-1) {
            if ((1-baseprobs[i].gap)*(1-baseprobs[i+1].gap)==1)
                g2_pen(id,ps_color_red);
            else
                g2_pen(id,ps_grey_colors[(int)floor((1-baseprobs[i].gap)*(1-baseprobs[i+1].gap)*nugrey_colors_-1)]);

            g2_line(id,X[i],Y[i],X[i+1],Y[i+1]);
        }
    }

    // draw pairings
    // !!! pair_table indexing begins at 1 !!!
    for (i=0; i<structure.size(); i++) {
        if ((unsigned short)pair_table[i+1]>i+1) {
            // pairs in both structures
            if (pairprob[i]==1)
                g2_pen(id,ps_color_red);
            else
                g2_pen(id,ps_grey_colors[(int)floor(pairprob[i]*nugrey_colors_-1)]);

            g2_line(id,X[i],Y[i],X[pair_table[i+1]-1],Y[pair_table[i+1]-1]);
        }
    }

    g2_flush(id);
    g2_close(id);

    free(pair_table);
    DELETE(X);
    DELETE(Y);
}
#endif
#endif

void RNAProfileAlignment::printSeqAli() const {
    unsigned int i,l;
    std::vector<std::string> seqs;
    std::string seq,info;

    // get alignment rows
    for (i=0; i<numStructures_; i++) {

        getSeqAli(seq,i,0,getMaxLength(0));
        seqs.push_back(seq);
    }

    l=seq.length();

    // sequence
    //	TODO if(hasSequence)
    //	{
    // calculate info line
    for (i=0; i<l; i++) {
        bool equal=true;

        for (unsigned int r=1; r<numStructures_; r++) {
            if ((seqs[r])[i]!=(seqs[r-1])[i]) {
                equal=false;
                break;
            }
        }

        info += equal ? "*" : " ";
    }

    // print it
    // sequences
    for (i=0; i<l; i+=55) {
        for (unsigned int r=0; r<numStructures_; r++)
            std::cout << std::setw(20) << std::setfill(' ') << std::left << strNames_[r].substr(0,20) << std::setw(5) << " " << RNAFuncs::UpperCase(seqs[r].substr(i,55)) << std::endl;

        std::cout << std::setw(25) << " " << info.substr(i,55) << std::endl;
        std::cout << std::endl;
    }
    //	}
}

std::vector<std::pair<std::string,std::string> > RNAProfileAlignment::getSeqAli() const {
    unsigned int i;
    std::vector<std::pair<std::string,std::string> > seqAli;
    std::pair<std::string,std::string> p;
    std::string seq;

    // get alignment rows
    for (i=0; i<numStructures_; i++) {
        getSeqAli(seq,i,0,getMaxLength(0));
        p.first = strNames_[i];
        p.second = seq;
        seqAli.push_back(p);
    }

    return seqAli;
}

void RNAProfileAlignment::printStrAli() const {
    unsigned int i,l;
    std::vector<std::string> strs;
    std::string str,info;

    // get alignment rows
    for (i=0; i<numStructures_; i++) {
        getStructAli(str,i);

        strs.push_back(str);
    }

    l=str.length();


    for (i=0; i<l; i++) {
        bool equal=true;

        for (unsigned int r=1; r<numStructures_; r++) {
            if ((strs[r])[i]!=(strs[r-1])[i]) {
                equal=false;
                break;
            }
        }

        info += equal ? "*" : " ";
    }

    // structures
    for (i=0; i<l; i+=55) {
        for (unsigned int r=0; r<numStructures_; r++)
            std::cout << std::setw(20) << std::setfill(' ') << std::left << strNames_[r].substr(0,20) << std::setw(5) << " " << strs[r].substr(i,55) << std::endl;

        std::cout << std::setw(25) << " " <<  info.substr(i,55) << std::endl;
        std::cout << std::endl;
    }


}

std::vector<std::string> RNAProfileAlignment::getStrAli() const {
    unsigned int i;
    std::vector<std::string> strs;
    std::string str;

    // get alignment rows
    for (i=0; i<numStructures_; i++) {
        getStructAli(str,i);
        strs.push_back(str);
    }

    return strs;
}

std::string RNAProfileAlignment::getConsSeq() const {
    std::string seq;
    std::vector<BaseProbs> baseprobs;
    getSequenceAlignment(baseprobs);

    // build sequence
    std::vector<BaseProbs>::const_iterator it;
    for (it=baseprobs.begin(); it!=baseprobs.end(); it++) {
        seq += getMlBase(*it);
    }

    return seq;
}

std::string RNAProfileAlignment::getConsStr(double minPairProb) const {
    std::string str;
    std::vector<double> pairprob;
    getStructureAlignment(minPairProb,str,pairprob);

    return str;
}

std::vector<double> RNAProfileAlignment::getBaseProb() const {
    std::vector<BaseProbs> baseprobs;
    std::vector<double> bestBaseprob;

    getSequenceAlignment(baseprobs);
    // build sequence
    std::vector<BaseProbs>::const_iterator it;
    for (it=baseprobs.begin(); it!=baseprobs.end(); it++) {
        bestBaseprob.push_back(getMlBaseFreq(*it));
    }

    return bestBaseprob;
}

//std::vector<pair<pair<int,int>,double> > RNAProfileAlignment::getPairProb(double minPairProb)
void RNAProfileAlignment::getPairProb(double &minPairProb, std::vector<std::pair<std::pair<int,int>,double> > &pairprobs) {
    //std::vector<pair<pair<int,int>,double> > structureProbs;
    std::vector<double> pairprob;
    std::string str;

    getStructureAlignment(minPairProb,str,pairprob);

    // store the position of open brackets in an integer vector
    std::vector<int> openBrackets;
    // iterator for structureProbs
    std::pair<int,int> tmpPair;
    std::pair<std::pair<int,int>,double> tmpPair2;

    // iterate through structure str
    for (unsigned int i=0; i<=str.length(); i++) {
        if (str[i]=='(') {
            openBrackets.push_back(i);
        } else if (str[i]==')') {
            // insert the int value of the corresponding opening bracket into structureProbs
            // insert position of the actual closing bracket
            // insert probability according to the actual bases
            tmpPair = std::pair<int,int> (openBrackets.back(),i);
            tmpPair2 = std::pair<std::pair<int,int>,double> (tmpPair,pairprob[i]);
            pairprobs.push_back(tmpPair2);
            // delete position of the last opening bracket
            openBrackets.pop_back();
        }
    }
}

void RNAProfileAlignment::printFastaAli(bool noStructure) const {
    std::string str,seq;

    // get alignment rows
    for (unsigned int i=0; i<numStructures_; i++) {
        getStructAli(str,i);
        getSeqAli(seq,i,0,getMaxLength(0));

        std::cout << ">" << strNames_[i] << std::endl;
        std::cout << seq << std::endl;
        if (!noStructure)
            std::cout << str << std::endl;
    }
}

void RNAProfileAlignment::printConsensus(double minPairProb) const {
    const int mountain_high=10;
    int j;

    std::string seq,structure;

    std::vector<BaseProbs> baseprobs;
    std::vector<double> pairprob;

    getSequenceAlignment(baseprobs);
    getStructureAlignment(minPairProb,structure,pairprob);

    // build sequence
    std::vector<BaseProbs>::const_iterator it;
    for (it=baseprobs.begin(); it!=baseprobs.end(); it++) {
        seq += getMlBase(*it);
    }

//	assert(seq.size() == structure.size());
    if (seq.size() != structure.size())
        std::cerr <<  "Error in resolving consensus structure!" << std::endl;

    for (unsigned int i=0; i<seq.length(); i+=55) {
        // show structure frequency mountain
        for (j=0; j<mountain_high; j++) {
            std::cout << std::setw(20) << std::setfill(' ') << " ";
            std::cout << std::setw(3) << std::right << 100 - (100 / mountain_high) * j << "% ";

            for (int k=0; k<55; k++) {
                if (i+k<seq.length()) {
                    if (getMlBaseFreq(baseprobs[i+k])>=1-(1/(double)mountain_high)*(double)j)
                        std::cout << "*";
                    else
                        std::cout << " ";
                }
            }
            std::cout << std::endl;
        }

        std::cout << std::setw(25) << " " << RNAFuncs::UpperCase(seq.substr(i,55)) << std::endl;
        std::cout << std::setw(25) << " " << structure.substr(i,55) << std::endl;

        // show structure frequency mountain
        for (j=0; j<mountain_high; j++) {
            std::cout << std::setw(20) << std::setfill(' ') << " ";
            std::cout << std::setw(3) << std::right << (100 / mountain_high) * (j+1) << "% ";

            for (int k=0; k<55; k++) {
                if (i+k<seq.length()) {
                    if (pairprob[i+k]>=(1/(double)mountain_high)*(double)j)
                        std::cout << "*";
                    else
                        std::cout << " ";
                }
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
    }
}

void RNAProfileAlignment::addStrNames(const std::vector<std::string>& strNames) {
    std::vector<std::string>::const_iterator it;
    for (it=strNames.begin(); it!=strNames.end(); it++)
        strNames_.push_back(*it);
}

void RNAProfileAlignment::save(const std::string &filename) {
    std::ofstream s(filename.c_str());

    // save the pforest to stream in binary format

    // first of all save the size
    s.write(reinterpret_cast<const char*>(&size_),sizeof(size_type));
    s.write(reinterpret_cast<const char*>(&numStructures_),sizeof(unsigned int));

    // save the arrays
    for (unsigned int i=0; i<size_; i++) {
        for (int r=0; r<RNA_ALPHABET_SIZE; r++)
            s.write(reinterpret_cast<const char*>(&lb_[i].p[r]),sizeof(double));

        s.write(reinterpret_cast<const char*>(lb_[i].columnStr.c_str()),sizeof(char)*numStructures_);
    }


    //  for(unsigned int i=0;i<size_;i++)
    //  {
    //    s.write(reinterpret_cast<char*>(&lb_[i]),sizeof(label_type)*size_);
    //  }

    //  s.write(reinterpret_cast<char*>(lb_),sizeof(label_type)*size_);
    s.write(reinterpret_cast<const char*>(rb_),sizeof(size_type)*size_);
    s.write(reinterpret_cast<const char*>(noc_),sizeof(size_type)*size_);

    for (unsigned int r=0; r<numStructures_; r++) {
        std::ostringstream ss;

        ss << std::setw(20) << std::setfill(' ') << std::left <<  strNames_[r];
        s.write(reinterpret_cast<const char*>(ss.str().c_str()),sizeof(char)*20);
    }


    //  s.write(reinterpret_cast<char*>(sumUpCSF_),sizeof(size_type)*size_);
    //  s.write(reinterpret_cast<char*>(rmb_),sizeof(size_type)*size_);

    //  s.write(name_;
    //unsigned int numStructures_;
    //unsigned int numStructuresX_;
    //unsigned int numStructuresY_;
    //std::vector<std::string> strNames_;
    //bool hasSequence;
}
