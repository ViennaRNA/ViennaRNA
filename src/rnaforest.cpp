#ifndef WIN32
#include "config.h"
#endif

#ifndef HAVE_LIBRNA
#undef HAVE_LIBG2
#endif

#include <string.h>

#include "rnaforest.h"
#include "rnafuncs.h"

#include "forest.t.cpp"

/* ****************************************** */
/*    Constructor and Destruktor functions    */
/* ****************************************** */

RNAForest::RNAForest(const std::string &baseStr, const std::string &viennaStr, const std::string &name, bool anchored)
        : Forest<RNA_Alphabet>(RNAFuncs::treeSize(viennaStr, anchored)),
        baseStr_(baseStr),
        viennaStr_(viennaStr) {
    if (name.empty())
        name_="unknown";
    else
        name_=name;

		if (anchored)
			buildForestAnchored(baseStr,viennaStr);
		else
			buildForest(baseStr,viennaStr);
};

RNAForest::~RNAForest() {
	// TODO
		DELETE_ARRAY(anchors_);
};

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

inline void RNAForest::makeLabel(RNA_Alphabet &a,char c) {
    a=c;
}

inline void RNAForest::showLabel(std::ostream &s,RNA_Alphabet a) const {
    s << a;
}

void RNAForest::buildForest(const std::string &baseStr, const std::string &viennaStr) {
    unsigned long basePairCount=0,maxDepth=0,stackPtr;
    unsigned int baseStrLen,viennaStrLen,node,*nodeStack, *numChildrenStack;

    assert(RNAFuncs::isRNAString(baseStr));

    RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);

    baseStrLen=baseStr.length();
    viennaStrLen=viennaStr.length();

    // check if base string and vienna string have the same length
    //  if(baseStr.length() > 0)
    //    if(baseStr.length() != viennaStrLen)
    //      throw RNAForestExceptionInput(RNAForestExceptionInput::Error_BaseStringAndViennaStringIncompatible);

    nodeStack=new unsigned int[maxDepth+1];
    numChildrenStack=new unsigned int[maxDepth+1];
    memset(nodeStack,0,sizeof(unsigned int)*maxDepth+1);
    memset(numChildrenStack,0,sizeof(unsigned int)*maxDepth+1);

    // fill Forest structure
    stackPtr=0;
    node=0;
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

            node++;
            break;
        }
    }

    delete[] nodeStack;
    delete[] numChildrenStack;
    calcSumUpCSF();
    calcRMB();

    assert (noc_[size_-1]==0);
}

// baut einen Baum mit anchors aus einem string mit anchoring-symbolen
// TODO wir brauchen auch die Laenge ohne Anchors, damit wir den Baum auf die richtige Groesse bringen koennen
void RNAForest::buildForestAnchored(const std::string &base, const std::string &anchored_vienna) {
    unsigned long length = 0, basePairCount = 0, maxDepth = 0;

    assert(base.length() <= anchored_vienna.length());
    // TODO assert(RNAFuncs::isAnchoredRNAString(base));
    RNAFuncs::isAnchoredViennaString(anchored_vienna, length, basePairCount, maxDepth);
    assert(length <= base.length());

    unsigned int *nodeStack, *numChildrenStack;
    nodeStack=new unsigned int[maxDepth+1];
    numChildrenStack=new unsigned int[maxDepth+1];
    memset(nodeStack,0,sizeof(unsigned int)*maxDepth+1);
    memset(numChildrenStack,0,sizeof(unsigned int)*maxDepth+1);


    unsigned long nodes = base.length()+basePairCount;
    // TODO richtig via oberklassenkonstruktor machen
    setSize(nodes);

    rb_ = new size_type[nodes];
    noc_ = new size_type[nodes];
    sumUpCSF_ = new size_type[nodes];
    rmb_ = new size_type[nodes];
    anchors_ = new size_type[nodes];
		lb_=new RNA_Alphabet[nodes];
    memset(rb_,0,sizeof(size_type)*nodes);
    memset(noc_,0,sizeof(size_type)*nodes);
    memset(sumUpCSF_,0,sizeof(size_type)*nodes);
    memset(rmb_,0,sizeof(size_type)*nodes);
    memset(anchors_,0,sizeof(size_type)*nodes);
    memset(lb_,0,sizeof(RNA_Alphabet)*nodes);

    // keep track of the number of built nodes
    unsigned long node = 0;
    unsigned long stackPtr = 0;

    //std::cout << "reading anchored_vienna string " << anchored_vienna << std::endl;
    int anchor = 1;
    // read anchored_vienna string char by char
    for (unsigned int i = 0, j = 0; i< base.length() && j < anchored_vienna.length();j++) {
        assert(node <= nodes);
        switch (anchored_vienna[j]) {
            case 'a': {
                // record read anchor
                setAnchor(node,anchor++);
                //std::cout << "node " << node << " got anchor " << getAnchor(node) << std::endl;
                break;
            }
            case '.': {
                // on reading a dot, just add a base
                // set label
                if (anchored_vienna.length()) // no basestring
                    lb_[node] = base[i++];
                else
                    lb_[node] = 'B';

                // set right brother (preorder index), 0 if struct ends with
                // this base
                if (node == size() - 1)
                    setRightBrotherIndex(node, 0);
                else
                    setRightBrotherIndex(node, node + 1);

                // this node as a base has no children
                setNumChildren(node, 0);

                // it is a child for its parent
                numChildrenStack[stackPtr]++;

                node++; // created a node -> increase node number
                break;
            }
            case '(': {
                // on reading a opening brace, build pair node,
                // descend, build opening base node
                // set label: pair node
                lb_[node] = 'P';

                // increase stack value: p node is a new child of its parent
                numChildrenStack[stackPtr]++;

                // descend - push node stacks
                stackPtr++;
                // store parent node number on node stack
                nodeStack[stackPtr] = node;
                //ancestorStack.push_back(node);
                // store noc on noc stack
                numChildrenStack[stackPtr] = 1;

                node++; // created a node -> increase node number

                // set label for opening base in basepair
                if (anchored_vienna.length()) // no basestring
                    lb_[node] = base[i++];
                else
                    lb_[node] = 'B';

                // set right brother - for opening base there must be one
                setRightBrotherIndex(node, node + 1);

                // this node, as an opening base, has no children
                setNumChildren(node, 0);

                node++; // created a node -> increase node number
                break;
            }
            case ')': {
                //3. on reading a closing brace, build closing base node
                // ascend
                // set label
                if (anchored_vienna.length()) // no basestring
                    lb_[node] = base[i++];
                else
                    lb_[node] = 'B';

                // set right brother - for closing base there is none
                setRightBrotherIndex(node, 0);

                // closing base has no children
                setNumChildren(node, 0);

                // ascend - pop node stack
                if (node == size() - 1) // sequence ends with closing brace
                    setRightBrotherIndex(nodeStack[stackPtr], 0);
                else
                    setRightBrotherIndex(nodeStack[stackPtr], node + 1);

                // add as a child to its parent
                setNumChildren(nodeStack[stackPtr],
                        numChildrenStack[stackPtr] + 1);
                // ascend
                stackPtr--;
                //ancestorStack.pop_back();

                node++;
                break;
              }
        }
    }
    calcSumUpCSF();
    calcRMB();

    delete[] nodeStack;
    delete[] numChildrenStack;
    assert(noc_[size_ - 1] == 0); // the rightmost node is a leaf

    //printMembers();
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
void RNAForest::plot2d(const std::string &filename_prefix, const std::vector<std::pair<unsigned int,unsigned int> > &regions, const RNAFuncs::SquigglePlotOptions &sqOptions) const {
    RNAFuncs::drawRNAStructure(baseStr_,viennaStr_,filename_prefix,name_,regions,sqOptions);
}
#endif
#endif
