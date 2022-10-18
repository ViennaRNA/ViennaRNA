/**********************************************/
/* BEGIN interface for distorted_sampling functions   */
/**********************************************/

%apply double *OUTPUT { double *result };
void computeInitialDistortion(vrna_fold_compound_t *vc, const char *s1, const char *s2, double *OUTPUT, double *OUTPUT);

void computeDistortion(vrna_fold_compound_t *vc, const char *s0, const char *s1, const char *s2, double *OUTPUT, double *OUTPUT);

%typemap(out) gridLandscapeT* estimate_landscape, gridLandscapeT* convertGrid_toList {
	//get number of cells
	int numberOfCells = 0;
	for(int i = 0; i < $1->size1;i++){
		for(int j = 0; j < $1->size2;j++){
			if($1->landscape[i][j].num_structs > 0){
				numberOfCells++;
			}
		}
	}
	
	int cellIndex = 0;
	PyObject* res = PyList_New(numberOfCells);
	for(int i = 0; i < $1->size1; i++){
		for(int j = 0; j < $1->size2;j++){
			if($1->landscape[i][j].num_structs > 0){
				PyObject* cellList = PyTuple_New(4);
				
				PyObject* pK = PyInt_FromLong($1->landscape[i][j].k);
				PyObject* pL = PyInt_FromLong($1->landscape[i][j].l);
				
				double num = $1->landscape[i][j].mfe;
				PyObject* pMFE = PyFloat_FromDouble(num);
				
				PyObject* structureList = PyList_New($1->landscape[i][j].num_structs);
				for(int k = 0; k < $1->landscape[i][j].num_structs; k++){
					PyObject* pStr = PyString_FromString($1->landscape[i][j].structures[k]);
					PyList_SetItem(structureList,k,pStr);
					free($1->landscape[i][j].structures[k]);
				}
				
				PyTuple_SET_ITEM(cellList,0,pK);
				PyTuple_SET_ITEM(cellList,1,pL);
				PyTuple_SET_ITEM(cellList,2,pMFE);
				PyTuple_SET_ITEM(cellList,3,structureList);
				PyList_SetItem(res,cellIndex,cellList);
				cellIndex++;
			}
			free($1->landscape[i][j].structures);
		}
		free($1->landscape[i]);
	}
	free($1->landscape);
	free($1);
	$result=res;
}



%{
#include <vector>
#include <string>
%}


%extend gridLandscapeT {
 void addStructure(char * structure){
 	addStructure($self, structure);
 }
};

%inline %{
#include <ViennaRNA/params.h>

void rescaleEnergy(vrna_fold_compound_t *vc, double rescale){
	vrna_exp_params_rescale(vc, &rescale);
}
%}


//define only the conversion method and map the output via %typemap.
%inline %{
gridLandscapeT* convertGrid_toList(gridLandscapeT* grid) {
    return grid;
}
%}

%inline %{
 void addStructuresToGrid(gridLandscapeT* grid, std::vector<std::string> structures) {
    for(int i=0; i < structures.size();i++){
 		addStructure(grid,(char *)structures[i].c_str());
 	}
}
%}


%include "../src/distorted_sampling.h"
