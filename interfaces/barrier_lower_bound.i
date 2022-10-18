/**********************************************/
/* BEGIN interface for barrier_lowerbound functions   */
/**********************************************/

%typemap(out) float barrier_estimate_2D %{
	$result = PyFloat_FromDouble($1);
%}

%include "../src/barrier_lower_bound.h"
