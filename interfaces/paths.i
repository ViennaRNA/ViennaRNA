/**********************************************/
/* BEGIN interface for PathFinder functions   */
/**********************************************/


%typemap(out) vrna_path_t* convertVRNA_PATH_toList {
	int len = 0;
	while($1[len].s) len++;
	
	PyObject* res = PyList_New(len);
	for(int i = 0; i < len; i++){
		float e = $1[i].en;
		char* s = $1[i].s;
		PyObject* energy = PyFloat_FromDouble(e);
		PyObject* structure = PyString_FromString(s);
		PyObject* tupel = PyTuple_New(2);
		PyTuple_SET_ITEM(tupel,0,structure);
		PyTuple_SET_ITEM(tupel,1,energy); 
		PyList_SetItem(res,i,tupel);
	}
    $result = res;
    
}

//define only the conversion method and map the output via %typemap.
%inline %{
vrna_path_t* convertVRNA_PATH_toList(vrna_path_t* path) {
    return path;
}
%}

%include "../src/PathFinder.h"
