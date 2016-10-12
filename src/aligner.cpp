#include <unistd.h>

#include "aligner.h"
#include "config.h"

void printCluster(std::pair<double, RNAProfileAlignment*> cluster,
                  unsigned int clusterNr, double minPairProb, const Options & options) {
    std::cout << "RNA Structure Cluster Nr: " << clusterNr << std::endl;
    std::cout << "Score: " << cluster.first << std::endl;
    std::cout << "Members: " << cluster.second->getNumStructures() << std::endl << std::endl;
    if (options.has(Options::FastaOutput)) {
        cluster.second->printFastaAli(false);
        std::cout << std::endl;
    } else {

			if (not options.has(Options::ShowOnlyScore)){
        cluster.second->printSeqAli();

#ifdef HAVE_LIBRNA
        // print alignment
        cluster.second->printStrAli();
        std::cout << std::endl;
#endif
			}
    }

    // save profile
    if (options.has(Options::SaveProfile)) {
        std::string filename;
        filename = options.generateFilename(Options::SaveProfile,".pro", "rna.pro",clusterNr);
        cluster.second->save(filename);
    }


		if (not options.has(Options::ShowOnlyScore)){
			// print consensus structure
			std::cout << "Consensus sequence/structure:" << std::endl;
			cluster.second->printConsensus(minPairProb);
		}
    std::cout << std::endl;

#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
    // generate squiggle plots
    if (options.has(Options::MakeSquigglePlot)) {
        RNAProfileAlignment::SquigglePlotOptions sqOptions;
        std::string filename;

        // plot showing full base information
        filename = options.generateFilename(Options::MakeSquigglePlot,".ps", "rnaprofile.ps",clusterNr);
        sqOptions.greyColors = options.has(Options::SquiggleGreyColors);
        sqOptions.minPairProb = minPairProb;
        sqOptions.mostLikelySequence = false;
        cluster.second->squigglePlot(filename,sqOptions);

        // plot showing consensus sequence
        filename = options.generateFilename(Options::MakeSquigglePlot,"_cons.ps", "rnaprofile_cs.ps",clusterNr);
        sqOptions.mostLikelySequence = true;
        cluster.second->squigglePlot(filename,sqOptions);
    }
#endif
#endif
}



void alignMultiple(std::vector<RNAProfileAlignment*> &inputListMult, Score &score,const Options &options,bool anchoring,RNAFuncs::AddXmlInfos &xmlInfos) {

    std::vector<std::pair<double,RNAProfileAlignment*> > resultList;
    progressiveAlign(inputListMult, resultList, score, options, anchoring);

    double minPairProb;
    options.get(Options::ConsensusMinPairProb, minPairProb, 0.5);

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "*** Results ***" << std::endl << std::endl;
    std::cout << "Minimum basepair probability for consensus structure (-cmin): " << minPairProb << std::endl << std::endl;

    unsigned int clusterNr = 1;
    std::vector<std::pair<double,RNAProfileAlignment*> >::const_iterator cluster;
    for (cluster=resultList.begin(); cluster!=resultList.end(); cluster++) {
        printCluster(*cluster, clusterNr, minPairProb, options);
        clusterNr++;
    }
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    // generate xml
    std::string outputFile;
    if (options.has(Options::XmlOutputFile)) {
        options.get(Options::XmlOutputFile,outputFile,std::string(""));
    } else {
        outputFile = "result.xml";
    }

    if (options.has(Options::GenerateXML)) {
        RNAFuncs::printMAliXML(resultList,options,minPairProb,xmlInfos,outputFile);
    }
#endif
#endif

    /* TODO das geht alle noch nicht, geht aber in pairwise
    if(!options.has(Options::ShowOnlyScore))
    {
    	std::cout << "global optimal score: ";
    }
    std::cout << optScore << std::endl;

    vector<RNAProfileAlignment*>::iterator it=inputMapProfile.begin();
    RNAProfileAlignment *f=it;
    if(!options.has(Options::ShowOnlyScore))
    {
    	// generate dot file
    	makeDotFileAli(*f,options);

    	f->print();
    	std::cout << std::endl;
    	f->printConsensus();
    }
    */

    // save profile alignment to binary file
    /*	      if(options.has(Options::SaveMultipleAliFile))
    {
    	std::string filename;

    filename=generateFilename(options,Options::SaveMultipleAliFile,".mta", "unknown.dot");
    	std::ofstream s(filename.c_str());
    f1->save(s);
    }
    */

    /*
    // generate squiggle plots
    if(options.has(Options::MakeSquigglePlot))
    {
    	RNAProfileAlignment::SquigglePlotOptions sqOptions;
    		std::string filename;

    	// plot showing full base information
    	filename=options.generateFilename(Options::MakeSquigglePlot,".ps", "rnaprofile.ps");
    	sqOptions.greyColors=options.has(Options::SquiggleGreyColors);
    	sqOptions.mostLikelySequence=false;
    	f->squigglePlot(filename,sqOptions);

    	// plot showing consensus sequence
    	filename=options.generateFilename(Options::MakeSquigglePlot,"_cons.ps", "rnaprofile_cons.ps");
    	sqOptions.greyColors=options.has(Options::SquiggleGreyColors);
    	sqOptions.mostLikelySequence=true;
    	f->squigglePlot(filename,sqOptions);
    }
    */


}

//void computeSpaceTimeInfo(const RNAForest* f1, const RNAForest* f2, const Score &score, const Options &options) {
//
//    tms tmsStart, tmsEnd;
//    Algebra<double,RNA_Alphabet> *algGlobClassic = new DoubleDist_Algebra(score);
//
//    if (!options.has(Options::ShowOnlyScore)) {
//        std::cout << "F1_NUMNODES" << ";";
//        std::cout << "F2_NUMNODES" << ";";
//        std::cout << "F1_DEGREE" << ";";
//        std::cout << "F2_DEGREE" << ";";
//        std::cout << "F1_LEAVES" << ";";
//        std::cout << "F2_LEAVES" << ";";
//        std::cout << "F1_DEPTH" << ";";
//        std::cout << "F2_DEPTH" << ";";
//        std::cout << "F1_NUMCSFS" << ";";
//        std::cout << "F2_NUMCSFS" << ";";
//        std::cout << "TABLE_SIZE_4D" << ";";
//        std::cout << "TABLE_SIZE_2D" << ";";
//        std::cout << "GLOBALI_TIME" << ";";
//        std::cout << "GLOBALI_TIME_SPEEDUP" << ";";
//        std::cout << "LOCALALI_TIME" << ";";
//        std::cout << "LOCALALI_TIME_SPEEDUP" << ";";
//        std::cout << std::endl;
//    }
//
//    std::cout << f1->size() << "\t";
//    std::cout << f2->size() << "\t";
//    std::cout << f1->maxDegree() << "\t";
//    std::cout << f2->maxDegree() << "\t";
//    std::cout << f1->numLeaves() << "\t";
//    std::cout << f2->numLeaves() << "\t";
//    std::cout << f1->maxDepth() << "\t";
//    std::cout << f2->maxDepth() << "\t";
//    std::cout << f1->getNumCSFs() << "\t";
//    std::cout << f2->getNumCSFs() << "\t";
//    std::cout << f1->size()*f2->size()*f1->maxDegree()*f1->maxDegree() << "\t";
//    std::cout << f1->getNumCSFs()*f2->getNumCSFs() << "\t";
//
//		bool topdown = false;
//		bool anchored = options.has(Options::Anchoring);
//    // global alignment
//    {
//        times(&tmsStart);
//        AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,topdown,anchored,false,false);
//        times(&tmsEnd);
//        std::cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t";
//    }
//    // global alignment speedup
//    {
//        times(&tmsStart);
//        AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,topdown,anchored,false,true);
//        times(&tmsEnd);
//        std::cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t";
//    }
//    // local alignment
//    {
//        times(&tmsStart);
//        AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,topdown,anchored,true,false);
//        times(&tmsEnd);
//        std::cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t";
//    }
//    // local alignment speedup
//    {
//        times(&tmsStart);
//        AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,topdown,anchored,true,true);
//        times(&tmsEnd);
//        std::cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t";
//    }
//    std::cout << std::endl;
//}

void computeSpaceTimeInfo(const RNAForest* f1, const RNAForest* f2, const Score &score, const Options &options) {

//    tms tmsStart, tmsEnd;
//    Algebra<double,RNA_Alphabet> *algGlobClassic = new DoubleDist_Algebra(score);
    if (!options.has(Options::ShowOnlyScore)) {
        std::cout << "F1_NUMNODES" << ";";
        std::cout << "F2_NUMNODES" << ";";
        std::cout << "F1_DEGREE" << ";";
        std::cout << "F2_DEGREE" << ";";
        std::cout << "F1_LEAVES" << ";";
        std::cout << "F2_LEAVES" << ";";
        std::cout << "F1_DEPTH" << ";";
        std::cout << "F2_DEPTH" << ";";
        std::cout << "F1_NUMCSFS" << ";";
        std::cout << "F2_NUMCSFS" << ";";
        std::cout << "TABLE_SIZE_4D" << ";";
        std::cout << "TABLE_SIZE_2D" << ";";

				for (int topdown=0; topdown<=1; topdown++) {
					for (int anchored=0; anchored<=1; anchored++) {
						for (int local=0; local<=1; local++) {
							for (int affine=0; affine <=1; affine++) {
								if ((anchored && !options.has(Options::Anchoring))||(anchored && !topdown))
									continue;
								if (affine) 
									std::cout << "AFFINE_";
								else
									std::cout << "LINEAR_";
								if (topdown)
									std::cout << "TOPDOWN_";
								else
									std::cout << "BOTTOM_UP_";
								if (anchored)
									std::cout << "ANCHORED_";
								else
									std::cout << "NONANCHORED_";
								if (local)
									std::cout << "LOCAL" << ";";
								else
									std::cout << "GLOBAL" << ";";
							}
						}
					}
				}

        std::cout << std::endl;
    }
    tms tmsStart, tmsEnd;
    Algebra<double,RNA_Alphabet> *alg = new DoubleDist_Algebra(score);
    AlgebraAffine<double,RNA_Alphabet> *alg_affine = new AffineDoubleDist_Algebra(score);

    std::cout << f1->size() << "\t";
    std::cout << f2->size() << "\t";
    std::cout << f1->maxDegree() << "\t";
    std::cout << f2->maxDegree() << "\t";
    std::cout << f1->numLeaves() << "\t";
    std::cout << f2->numLeaves() << "\t";
    std::cout << f1->maxDepth() << "\t";
    std::cout << f2->maxDepth() << "\t";
    std::cout << f1->getNumCSFs() << "\t";
    std::cout << f2->getNumCSFs() << "\t";
    std::cout << f1->size()*f2->size()*f1->maxDegree()*f1->maxDegree() << "\t";
    std::cout << f1->getNumCSFs()*f2->getNumCSFs() << "\t";

		bool printBT = true; // never enabled
		bool speedup = true; // always enabled

		for (int topdown=0; topdown<=1; topdown++) {
			for (int anchored=0; anchored<=1; anchored++) {
				for (int local=0; local<=1; local++) {
					for (int affine=0; affine <=1; affine++) {

						// TODO man muesste eig abpruefen, ob der input fuer's anchoring ausreicht, nicht, ob es option ist
						if ((anchored && !options.has(Options::Anchoring))||(anchored && !topdown))
							continue;

						std::cout << "\t" << std::flush;
						//std::cout << std::endl;
						//if (affine)
						//	std::cout << "affine ";
						//else
						//	std::cout << "linear ";
						//if (topdown)
						//	std::cout << "topdown ";
						//else
						//	std::cout << "bottomup ";
						//if (anchored)
						//	std::cout << "anchored ";
						//else 
						//	std::cout << "nonanchored ";
						//if (local)
						//	std::cout << "local ";
						//else
						//	std::cout << "global ";
						//std::cout << std::endl;

						if (affine) {
	        		times(&tmsStart);
	        		AlignmentAffine<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*alg_affine,topdown,anchored,local,!printBT,speedup);
		     		 	times(&tmsEnd);
			  		  std::cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t" << std::flush;
						}
						else {
				  		times(&tmsStart);
							AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*alg,topdown,anchored,local,!printBT,speedup);
						 	times(&tmsEnd);
						  std::cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t" << std::flush;
						}

					}
				}
			}
		}
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}

void alignPairwise(std::vector<RNAForest*> &inputListPW, Score &score, const Options &options, bool anchored, RNAFuncs::AddXmlInfos &xmlInfos) {

#ifdef HAVE_LIBG2
    RNAFuncs::SquigglePlotOptions sqOptions;
    // generate squiggle plot
    if (options.has(Options::MakeSquigglePlot)) {
        // get sq options
        sqOptions.hideBaseNumbers = options.has(Options::SquiggleHideBaseNumbers);
        options.get(Options::SquiggleBaseNumberInterval,sqOptions.baseNumInterval,(unsigned int)20);
        sqOptions.greyColors = options.has(Options::SquiggleGreyColors);
        options.get(Options::SquiggleScaleFactor,sqOptions.scale,1.0);
        sqOptions.generateFIG = options.has(Options::SquiggleGenerateFIG);
#ifdef HAVE_LIBGD
        sqOptions.generatePNG = options.has(Options::SquiggleGeneratePNG);
        sqOptions.generateJPG = options.has(Options::SquiggleGenerateJPG);
#endif

    }
#endif

    RNAForest *f1 = inputListPW.front();
    RNAForest *f2 = inputListPW.back();

    if (options.has(Options::SpaceTimeInfo)) {
        computeSpaceTimeInfo(f1, f2, score, options);
        exit(EXIT_SUCCESS);
    }

    // for global we compute less entries in the matrix
		bool topdown = options.has(Options::Topdown);
		//bool anchored = options.has(Options::Anchoring);
    bool local = options.has(Options::LocalSimilarity);
		bool printBT = options.has(Options::Backtrace);

		RNA_Algebra<double,RNA_Alphabet> *alg = NULL;
    RNA_AlgebraAffine<double,RNA_Alphabet> * alg_affine = NULL;
		Alignment<double,RNA_Alphabet,RNA_AlphaPair> * ali = NULL;

		if (options.has(Options::Affine)) {
        if (options.has(Options::CalculateDistance))
            alg_affine = new AffineDoubleDistRNA_Algebra(score);
        else if (options.has(Options::RIBOSUMScore))
            alg_affine = new AffineRIBOSUM8560(score);
        else
            alg_affine = new AffineDoubleSimiRNA_Algebra(score);
        ali = new AlignmentAffine<double,RNA_Alphabet,RNA_AlphaPair>(f1,f2,*alg_affine,topdown,anchored,local,printBT);
		}
		else {
        if (options.has(Options::CalculateDistance))
            alg = new DoubleDistRNA_Algebra(score);
        else if (options.has(Options::RIBOSUMScore))
            alg = new RIBOSUM8560(score);
        else
            alg = new DoubleSimiRNA_Algebra(score);

        ali = new AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair>(f1,f2,*alg,topdown,anchored,local,printBT);
		}
			
 		if (options.has(Options::Tables)) { // TODO mit stringstream zusammenbauen
 			if (options.has(Options::Topdown)) {
 				if (options.has(Options::Affine))
 					ali->printToFile(std::string("tables_affine_topdown.txt"));
 				else
 					ali->printToFile(std::string("tables_linear_topdown.txt"));
 			}
 			else {
 				if (options.has(Options::Affine))
 					ali->printToFile(std::string("tables_affine.txt"));
 				else
 					ali->printToFile(std::string("tables_linear.txt"));
 			}
 		}

    RNA_Alignment fali;
    fali.setStructureNames(f1->getName(),f2->getName());

#if defined(HAVE_LIBRNA) && defined(HAVE_LIBXMLPLUSPLUS) && defined(HAVE_LIBXML2)
    double optScore = getAndPrintOptimum(ali, options);
#else
    getAndPrintOptimum(ali, options);
#endif

            std::vector<std::pair<unsigned int,unsigned int> > xsubopts;
            std::vector<std::pair<unsigned int,unsigned int> > ysubopts;
    if (!options.has(Options::ShowOnlyScore)) {

            // read options for local simi
            int suboptPercent;
            options.get(Options::LocalSubopts,suboptPercent,100);
            // for local simi
            unsigned int xbasepos,ybasepos,xlen,ylen;
            if (options.has(Options::SmallInLarge)) {
                ali->getOptSILAlignment(fali,ybasepos);
                std::cout << "starting at position: " << ybasepos << std::endl << std::endl;
            } 
						else if (options.has(Options::LocalSimilarity)) {
                ali->resetOptLocalAlignment(suboptPercent);
                ali->getOptLocalAlignment(fali,xbasepos,ybasepos);
                std::cout << "starting at positions: " << xbasepos << "," << ybasepos << std::endl << std::endl;
                xmlInfos.xbasepos = xbasepos;
                xmlInfos.ybasepos = ybasepos;
            } 
						else
                ali->getOptGlobalAlignment(fali);

            // generate dot file
            makeDotFileAli(fali,options);

            // print alignment
            //printRNAAlignment(RNA_Alignment fali, Options options);
            std::string seq1,seq2,str1,str2;
            fali.getSequenceAlignments(seq1,seq2);
            fali.getStructureAlignment(str1,true);
            fali.getStructureAlignment(str2,false);

            if (options.has(Options::FastaOutput)) {
                std::cout << fali.getStructureNameX() << std::endl;
                std::cout << seq1 << std::endl;
                std::cout << str1 << std::endl;
                std::cout << fali.getStructureNameY() << std::endl;
                std::cout << seq2 << std::endl;
                std::cout << str2 << std::endl;
                std::cout << std::endl;
            } 
						else {
                RNAFuncs::printAli(fali.getStructureNameX(),fali.getStructureNameY(),seq1,seq2,str1,str2);
            }

            xlen = seq1.size();
            ylen = seq2.size();
            if (options.has(Options::LocalSimilarity)) {
                xsubopts.push_back(std::pair<unsigned int,unsigned int>(xbasepos,xlen));
                ysubopts.push_back(std::pair<unsigned int,unsigned int>(ybasepos,ylen));
            }

            // for squiggle plot and subopts squiggle plot
            unsigned int count = 1;

#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
            char s[8];
            if (options.has(Options::MakeSquigglePlot)) {
                sprintf(s,"%d",count);
                fali.squigglePlot(s,sqOptions);
            }
#endif
#endif

            // suboptimal alignments
            if (options.has(Options::LocalSubopts)) {
                while (ali->nextLocalSuboptimum()) {
                    count++;
                    std::cout << "local optimal score: ";
                    std::cout << ali->getLocalOptimum() << std::endl;
										RNA_Alignment fali;
										unsigned int xbasepos = 10, ybasepos = 10;
                    ali->getOptLocalAlignment(fali,xbasepos,ybasepos);
                    std::cout << "starting at positions: " << xbasepos << "," << ybasepos << std::endl << std::endl;

                    // ab hier und das pushback herausnehmen - fasta option abfragen (print subopt)
                    // print alignment
                    fali.getSequenceAlignments(seq1,seq2);
                    fali.getStructureAlignment(str1,true);
                    fali.getStructureAlignment(str2,false);
                    RNAFuncs::printAli(fali.getStructureNameX(),fali.getStructureNameY(),seq1,seq2,str1,str2);
                    xlen = seq1.size();
                    ylen = seq2.size();
                    xsubopts.push_back(std::pair<unsigned int,unsigned int>(xbasepos,xlen));
                    ysubopts.push_back(std::pair<unsigned int,unsigned int>(ybasepos,ylen));

#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
                    if (options.has(Options::MakeSquigglePlot)) {
                        sprintf(s,"%d",count);
                        fali.squigglePlot(s,sqOptions);
                    }
#endif
#endif
                }
            }
        }

        // generateXML(Options & options, RNA_Alignment & fali, seq1, seq2, str1, str2, optScore, xmlInfos)
        // NOTE der folgende Teil druckt bei den subopts nur das optimale
#ifdef HAVE_LIBRNA
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
        // generate xml
        std::string outputFile;
        if (options.has(Options::XmlOutputFile)) {
            options.get(Options::XmlOutputFile,outputFile,std::string(""));
        } else {
            outputFile = "result.xml";
        }
        if (options.has(Options::GenerateXML)) {
            RNAFuncs::printPAliXML(fali.getStructureNameX(),fali.getStructureNameY(),seq1,seq2,str1,str2,optScore,options,xmlInfos,outputFile);
        }
#endif
#endif
#endif

        // der plot2d bekommt alle opt und subopt. alis zum ausdruck
#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
        // generate squiggle plot
        if (options.has(Options::MakeSquigglePlot)) {
            f1->plot2d("x_str", xsubopts, sqOptions);
            f2->plot2d("y_str", ysubopts, sqOptions);
        }
#endif
#endif

        // clear input vector
        //std::vector<RNAForest*>::const_iterator it;
        //for (it = inputListPW.begin(); it!=inputListPW.end(); it++)
        //    delete *it;
        inputListPW.clear();
				if (options.has(Options::Affine)) 
					delete alg_affine;
				else
					delete alg;
}



void editPairwise(std::vector<RNAForestSZ*> &inputListSZ,Score &score,Options &options, bool anchored) {
    //  timeb t1,t2;
    SZAlgebra<double,RNA_Alphabet> *alg;

//  if(options.has(Options::CalculateDistance))
    //alg=new DoubleDistSZAlgebra(score);
    alg=new IntDistSZAlgebra(score);
//  else
//    alg=new DoubleSimiSZAlgebra(score);

    RNAForestSZ *f1=inputListSZ.front();
    RNAForestSZ *f2=inputListSZ.back();

    //  ftime(&t1);
    Mapping<double,RNA_Alphabet> mapping(f1,f2,*alg);
    //  ftime(&t2);

    std::cout << "Global optimum: " << mapping.getGlobalOptimum() << std::endl;
    //  std::cout << "Calculation Time ms: " << (t2.time*1000+t2.millitm) - (t1.time*1000+t1.millitm) << std::endl;
}


void alignPairwiseSimple(std::vector<RNAForest*> &inputListPW,Score &score,Options &options, bool anchored) {
    //  timeb t1,t2;
    Algebra<double,RNA_Alphabet> *alg;

    //  if(options.has(Options::CalculateDistance))
    alg = new DoubleDist_Algebra(score);
//  else
//    alg=new ScoreAlgebraSimple(score);

    RNAForest *f1 = inputListPW.front();
    RNAForest *f2 = inputListPW.back();

    //  ftime(&t1);
		//  no local, no print backtrace, but speedup
    AlignmentLinear<double,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*alg,false,false,false,true);
    //  ftime(&t2);

    std::cout << "Global optimum: " << ali.getGlobalOptimum() << std::endl;
    //  std::cout << "Calculation Time ms: " << (t2.time*1000+t2.millitm) - (t1.time*1000+t1.millitm) << std::endl;
}

double getAndPrintOptimum(Alignment<double,RNA_Alphabet,RNA_AlphaPair>* ali, const Options & options) {
    if (!options.has(Options::ShowOnlyScore)) {
        if (options.has(Options::SmallInLarge))
            std::cout << "small-in-large optimal score: ";
        else if (options.has(Options::LocalSimilarity))
            std::cout << "local optimal score: ";
        else
            std::cout << "global optimal score: ";
    }

    if (options.has(Options::SmallInLarge)) {
        std::cout << ali->getSILOptimum() << std::endl;
        return ali->getSILOptimum();
    } 
		else if (options.has(Options::LocalSimilarity)) {
        std::cout << ali->getLocalOptimum() << std::endl;
        return ali->getLocalOptimum();
    } 
		else {
        std::cout << ali->getGlobalOptimum() << std::endl;
        if (options.has(Options::RelativeScore))
            std::cout << ali->getGlobalOptimumRelative() << std::endl;
        return ali->getGlobalOptimum();
    }
}


