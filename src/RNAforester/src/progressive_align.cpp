#include <algorithm>
#include <fstream>
#include "alignment.h"
#include "progressive_align.h"
#include "alignment.t.cpp"


// !!! this operator is defined as > !!!
bool operator < (std::pair<double,RNAProfileAlignment*> &l, std::pair<double,RNAProfileAlignment*> &r) {
    if (l.second->getNumStructures() < r.second->getNumStructures())
        return false;
    else
        return true;
}

void progressiveAlign(std::vector<RNAProfileAlignment*> &inputList, 
											std::vector<std::pair<double,RNAProfileAlignment*> > &resultList, const Score& score, const Options &options, bool anchored) {

		Algebra<double,RNA_Alphabet_Profile> *alg = NULL;
		AlgebraAffine<double,RNA_Alphabet_Profile> *alg_affine = NULL;

    if (options.has(Options::Affine)) {
			if (options.has(Options::CalculateDistance)) {
        alg_affine = new AffineDoubleDistProfileAlgebra(score);
			}
			else {
        alg_affine = new AffineDoubleSimiProfileAlgebra(score);
			}
		}
		else {
			// distance or similarity
			if (options.has(Options::CalculateDistance))
        alg = new DoubleDistProfileAlgebra(score);
			else
        alg = new DoubleSimiProfileAlgebra(score);
		}

    std::cout << "*** Calculation ***" << std::endl << std::endl;

    // create inputMapProfile to access a profile by an index value
    RNAProfileAliMapType inputMapProfile;
    std::vector<RNAProfileAlignment*>::const_iterator inpIt;
    long i = 1;
    for (inpIt=inputList.begin(); inpIt!=inputList.end(); inpIt++) {
        inputMapProfile[i]=*inpIt;
        i++;
    }
    inputList.clear();

    // create matrix for all against all comparison
		double bestScore;
		if (options.has(Options::Affine))
			bestScore = alg_affine->worst_score();
		else
			bestScore = alg->worst_score();
    Matrix<double> *score_mtrx = new Matrix<double>(inputMapProfile.size(),inputMapProfile.size());

    // set threshold for the clustering algorithm
    double threshold = 0; 
    if (options.has(Options::CalculateDistance))
        options.get(Options::ClusterThreshold, threshold, 20.0);
    else
        options.get(Options::ClusterThreshold, threshold, 0.7);
    std::cout << "clustering threshold is: " << threshold << std::endl;

    // set cutoff value for clustering
		double cutoff = 0;
    if (options.has(Options::CalculateDistance))
        options.get(Options::ClusterJoinCutoff, cutoff, 100.0);
    else
        options.get(Options::ClusterJoinCutoff, cutoff, 0.0);
    std::cout << "join clusters cutoff is: " << cutoff << std::endl << std::endl;

		bool topdown = options.has(Options::Topdown);
		//bool anchored = options.has(Options::Anchoring);
    bool local = options.has(Options::LocalSimilarity);
		bool printBT = options.has(Options::Backtrace);

    // generate dot file
		std::string clusterfilename = options.generateFilename(Options::Help,"_cluster.dot", "cluster.dot");  // use Help as dummy
    std::ofstream s;
    s.open(clusterfilename.c_str());
    s << "digraph forest" << std::endl << "{" << std::endl;

		// generate nodes for the input forests
	  RNAProfileAliMapType::iterator it;
    for (it = inputMapProfile.begin(); it!=inputMapProfile.end(); it++) {
        RNAProfileAlignment * f = it->second;

        s << "\"" << f->getName() << "\"" << "[label=\"" << f->getName() << "\"]" << std::endl;
        s << "\"" << f->getName() << "\"" << "[label=\"" << f->getName() << "\"]" << std::endl;
    }

		// compute all pairwise alignment scores
    // !! NOTE !! iterating through the map is ordered by key number
    // as i only calculate a triangle matrix this is a prerequisite
		long x = 0, y = 0;
		RNAProfileAlignment *f1 = NULL, *f2 = NULL;
    std::cout << "Computing all pairwise similarities" << std::endl;

    RNAProfileAliMapType::iterator it2;
    for (it=inputMapProfile.begin(); it!=inputMapProfile.end(); it++) {
        x = it->first;
        f1 = it->second;
        for (it2=inputMapProfile.begin(); it2->first<it->first; it2++) {
            y = it2->first;
            f2 = it2->second;

						if (options.has(Options::Affine)) {
							AlignmentAffine<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> * ali = new AlignmentAffine<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>(f1,f2,*alg_affine,topdown,anchored,local,printBT);
	            if (local)
		              score_mtrx->setAt(x-1,y-1,ali->getLocalOptimum());
			        else
				          score_mtrx->setAt(x-1,y-1,ali->getGlobalOptimumRelative());

					    std::cout << x << "," << y << ": " << score_mtrx->getAt(x-1,y-1) << std::endl;
							delete ali;
						}
						else {
							AlignmentLinear<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> * ali = new AlignmentLinear<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>(f1,f2,*alg,topdown,anchored,local,printBT);
							if (local)
								  score_mtrx->setAt(x-1,y-1,ali->getLocalOptimum());
							else
								  score_mtrx->setAt(x-1,y-1,ali->getGlobalOptimumRelative());

							std::cout << x << "," << y << ": " << score_mtrx->getAt(x-1,y-1) << std::endl;
							delete ali;
						}
        }
    }
    std::cout << std::endl;

    std::vector<RNAProfileAliKeyPairType> inputListMult;
    long bestx = 0, besty = 0;  // keys of stuctures in inputMapProfile
    RNAProfileAlignment *f = NULL;
    int level = 1;


    while (inputMapProfile.size()>1) {
        // find the best score of all pairwise alignments
				if (options.has(Options::Affine))
					bestScore = alg_affine->worst_score();
				else
					bestScore = alg->worst_score();
        for (it=inputMapProfile.begin(); it!=inputMapProfile.end(); it++) {
            x = it->first;
            for (it2=inputMapProfile.begin(); it2->first < it->first; it2++) {
                double old_bestScore = bestScore;

                y = it2->first;

								if (options.has(Options::Affine))
									bestScore = alg_affine->choice(bestScore,score_mtrx->getAt(x-1,y-1));
								else
									bestScore = alg->choice(bestScore,score_mtrx->getAt(x-1,y-1));

                if (bestScore != old_bestScore) {
                    bestx = it->first;
                    besty = it2->first;
                    old_bestScore = bestScore;
                }
            }
        }

        std::cout << "joining alignments:" << std::endl;

        // test, if score is better than the threshold value
        // if threshold is set generate a vector of best pairs within threshold
				bool scoreBetterThanThreshold = false;
				if (options.has(Options::Affine))
					scoreBetterThanThreshold = (alg_affine->choice(bestScore,threshold)!= threshold);
				else
					scoreBetterThanThreshold = (alg->choice(bestScore,threshold)!= threshold);
        if (scoreBetterThanThreshold) {

            int maximize;

            Graph graph;
						if (options.has(Options::Affine))
							graph = makePairsGraph(inputMapProfile,alg_affine,score_mtrx,threshold);
						else
							graph = makePairsGraph(inputMapProfile,alg,score_mtrx,threshold);
            if (options.has(Options::CalculateDistance))
                maximize = 0;
            else
                maximize = 1;

            int *mate;
            mate = Weighted_Match(graph, 1, maximize);
            for (x=1; x<=score_mtrx->xDim(); x++) { // !! begins at 1 !!
                if (x<mate[x]) {
                    // if it is a best pair put it in the align vector
                    f1 = inputMapProfile[x];
                    inputMapProfile.erase(x);
                    inputListMult.push_back(std::make_pair(x,f1));
                    f2 = inputMapProfile[mate[x]];
                    inputMapProfile.erase(mate[x]);
                    inputListMult.push_back(std::make_pair(mate[x],f2));
                }
            }
            free(mate);
            free(graph);
        } else {
            // if there us no pair below the threshold
            // combine those two profile forests, that produced the best score
            f1 = inputMapProfile[bestx];
            f2 = inputMapProfile[besty];
            inputMapProfile.erase(bestx);
            inputMapProfile.erase(besty);
            inputListMult.push_back(std::make_pair(bestx,f1));
            inputListMult.push_back(std::make_pair(besty,f2));
        }

        // align all forests in the align vector.
        while (inputListMult.size()>1) {
            std::string aliName;

            x = inputListMult.front().first;
            f1 = inputListMult.front().second;
            inputListMult.erase(inputListMult.begin());
            y = inputListMult.front().first;
            f2 = inputListMult.front().second;
            inputListMult.erase(inputListMult.begin());

            // compute alignment again
            Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> * bestali = NULL;
						if (options.has(Options::Affine))
							bestali = new AlignmentAffine<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>(f1,f2,*alg_affine,topdown,anchored,local,printBT);
						else
							bestali = new AlignmentLinear<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>(f1,f2,*alg,topdown,anchored,local,printBT);
            if (local)
                bestScore = bestali->getLocalOptimum();
            else
                bestScore = bestali->getGlobalOptimumRelative();

            // test, if score is worse than the cutoff value
						bool scoreWorseThanCutoff = false;
						if (options.has(Options::Affine))
							scoreWorseThanCutoff = (alg_affine->choice(bestScore,cutoff) == cutoff);
						else
							scoreWorseThanCutoff = (alg->choice(bestScore,cutoff) == cutoff);

            if (scoreWorseThanCutoff) {
                // copy involved forests to the result vector
                std::cout << x << "," << y << ": alignment is below cutoff." << std::endl;
								std::cout << f1->getNumStructures() << "  " << f2->getNumStructures() << std::endl;
								std::cout << resultList.size() << std::endl;
								if (options.has(Options::Affine)) {
									if (f1->getNumStructures()>1)
                    resultList.push_back(std::make_pair(f1->maxScore(*alg_affine),f1));
									if (f2->getNumStructures()>1)
                    resultList.push_back(std::make_pair(f2->maxScore(*alg_affine),f2));
								}
								else {
									if (f1->getNumStructures()>1)
                    resultList.push_back(std::make_pair(f1->maxScore(*alg),f1));
									if (f2->getNumStructures()>1)
                    resultList.push_back(std::make_pair(f2->maxScore(*alg),f2));
								}
            } 
						else {
                // calculate optimal alignment and add it to inputMapProfile

								// the anchors are part of the profile alignment
                f = new RNAProfileAlignment(f1->getNumStructures(),f2->getNumStructures());

                if (local) {
                    unsigned int xbasepos,ybasepos;
                    bestali->getOptLocalAlignment(*f,xbasepos,ybasepos);
                } else {
                    bestali->getOptGlobalAlignment(*f);
                }

                f->addStrNames(f1->getStrNames());
                f->addStrNames(f2->getStrNames());
                aliName = f1->getName() + "." +f2->getName();
                f->setName(aliName);

                // for debug purposes
                //std::string dotfilename;
                //dotfilename=f->getName() + std::string(".dot");
                //				std::ofstream s("test.dot");
                //				f->printDot(s);

                // generate nodes in cluster file (dot format)
                s << "\"" << f->getName() << "\"" << "[shape=\"diamond\" label=\"" << bestScore << "\"]" << std::endl;
                s << "\"" << f->getName() << "\"" << "-> {\"" <<  f1->getName() << "\" \"" << f2->getName() << "\"}";

								long joinedClusterNumber = std::min(x, y);
                std::cout << x << "," << y << ": " << bestScore << " -> " << joinedClusterNumber << std::endl;

								// TODO problem when deleting (the anchors?)
                //delete f1;
                //delete f2;

                // calculate distance to all forests in the vector
                std::cout << "Calculate similarities to other clusters" << std::endl;
                f1 = f;
                // x remains x !!
								x = joinedClusterNumber;
                for (it=inputMapProfile.begin(); it!=inputMapProfile.end(); it++) {
                    y = it->first;
                    f2 = it->second;

										Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> * ali = NULL;
										if (options.has(Options::Affine))
											ali = new AlignmentAffine<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>(f1,f2,*alg_affine,topdown,anchored,local,printBT);
										else
											ali = new AlignmentLinear<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>(f1,f2,*alg,topdown,anchored,local,printBT);
                    if (local) {
                        score_mtrx->setAt(std::min(x-1,y-1),std::max(x-1,y-1),ali->getLocalOptimum());  // min - max = fill the upper triangle
                        std::cout << std::min(x,y) << "," << std::max(x,y) << ": " << ali->getLocalOptimum() <<  std::endl;
                    } 
										else {
                        score_mtrx->setAt(std::min(x-1,y-1),std::max(x-1,y-1),ali->getGlobalOptimumRelative());  // min - max = fill the upper triangle
                        std::cout << std::min(x,y) << "," << std::max(x,y) << ": " << ali->getGlobalOptimumRelative() <<  std::endl;
                    }
										delete ali;
                }
                std::cout << std::endl;

                // ... and append it to the vector
                inputMapProfile.insert(std::make_pair(x,f));
            }
						delete bestali;
        }

        level++;
    }

    assert(inputMapProfile.size()<2);

    // copy last profile to resultList
    if (inputMapProfile.size()==1) {
        f=inputMapProfile.begin()->second;
				if (options.has(Options::Affine))
					resultList.push_back(std::make_pair(f->maxScore(*alg_affine),f));
				else
					resultList.push_back(std::make_pair(f->maxScore(*alg),f));
        inputMapProfile.clear();
    }

    // end of dot file
    s << "}" << std::endl;

    // sort result vector
    //	std::sort(resultList.begin(),resultList.end());

		if (options.has(Options::Affine)) 
			delete alg_affine;
		else 
			delete alg;

    delete score_mtrx;
}

Graph makePairsGraph(const RNAProfileAliMapType &inputMapProfile, const Algebra<double,RNA_Alphabet_Profile> *alg, const Matrix<double> *score_mtrx, double threshold) {
    Graph graph;
    RNAProfileAliMapType::const_iterator it,it2;
    RNAProfileAlignment *f1=NULL,*f2=NULL;

    graph = NewGraph(score_mtrx->xDim());

    for (int i=0; i<score_mtrx->xDim(); i++) {
        Xcoord(graph,i+1) = 0;
        Ycoord(graph,i+1) = 0;
    }

    for (it=inputMapProfile.begin(); it!=inputMapProfile.end(); it++) {
        f1=it->second;
        for (it2=inputMapProfile.begin(); it2->first<it->first; it2++) {
            double score;
            f2=it2->second;

            score=score_mtrx->getAt(it->first-1,it2->first-1);
            if (alg->choice(score,threshold) != threshold) { // is it better than the threshold ?
                AddEdge (graph,it->first,it2->first,(int)(score*100.0));
            }
        }
    }


    WriteGraph (graph,(char*)"test.out");
    return graph;
}
