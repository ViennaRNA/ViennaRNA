/**********************************************/
/* BEGIN interface for distorted_sampling functions   */
/**********************************************/


// Map a Python sequence into any sized C double array
%typemap(in) double* {
 /* Check if is a list */
  if (PyList_Check($input)) {
  	int size = PyList_Size($input);
    int i = 0;
    $1 = (double *) malloc((size)*sizeof(double));
    for (i = 0; i < size; i++) {
       $1[i] = PyFloat_AsDouble(PyList_GetItem($input,i));
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $1[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($1);
        return NULL;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

%typemap(out) gridLandscapeT* estimate_landscapeMD {
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


%inline %{
	std::vector<double> computeDistortions(vrna_fold_compound_t* fc,const char **structures, size_t numberOfStructures){
	 	double * tmp = rxp_computeDistortionsWithMFE(fc, structures, numberOfStructures);
	 	std::vector<double> result;
	    for(int i=0; i < numberOfStructures;i++){
	 		result.push_back(tmp[i]);
	 	}
	 	return result;
    }
    
	std::vector<double> computeDistortionsMaxDist(vrna_fold_compound_t* fc,const char **structures, size_t numberOfStructures,double *distances){
	 	double * tmp = rxp_computeDistortionsWRTMaxDistance(fc, structures, numberOfStructures,distances);
	 	std::vector<double> result;
	    for(int i=0; i < numberOfStructures;i++){
	 		result.push_back(tmp[i]);
	 	}
	 	return result;
    }
%}


%include "../src/distorted_samplingMD.h"
%include "../src/dist_class_sc.h"
