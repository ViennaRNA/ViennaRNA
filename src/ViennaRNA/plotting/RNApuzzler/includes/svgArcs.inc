#ifndef RNAPUZZLER_SVG_ARCS
#define RNAPUZZLER_SVG_ARCS

/*
 *      RNApuzzler Calculation of arc data for svg files
 *
 *      ViennaRNA package
 */

#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"

#include "vector_math.inc"
#include "../headers/tBaseInformation_struct.h"


//Probably would be cleverer to move this method elsewhere and delete this file. It was planned that this File would have more methods but the transformation of the coords for the SVG Arc was easier than expected.
PRIVATE void transformPSArcsToSVG(int n, double * oldArcs, double ** newArcs){
	*newArcs = (double*)vrna_alloc(n*2*sizeof(double));
	for(int i=0;i<n;i++){
		//Arc exists
		if(oldArcs[6*i + 2] > 0){	//radius
			(*newArcs)[2*i] = oldArcs[6*i + 2];
			(*newArcs)[2*i +1] = oldArcs[6*i + 5];  //goClockwise
		} else{
			(*newArcs)[2*i] = -1;
			(*newArcs)[2*i+1] = -1;
		}
	}
	
	
}

#endif