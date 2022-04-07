/*#######################################*/
/* Typemap constraints and checks        */
/*#######################################*/


/**************************************************/
/* Input value type check for pair_table-like     */
/* var_array<short> types. It will be applied     */
/* to all occurences of var_array<short> derived  */
/* types with variable name 'pt', e.g.            */
/* foo(var_array<short> &pt);                     */
/**************************************************/

%typemap(check) var_array<short> pair_table_like {
   if ((!((*$1).type & VAR_ARRAY_LINEAR)) ||
       (!((*$1).type & VAR_ARRAY_ONE_BASED)) ||
       ((size_t)((*$1).data[0]) != (*$1).length ))
   {
       SWIG_exception(SWIG_ValueError,"Expected var_array<short> with pair_table properties, i.e. data[0] == length, type = VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED.");
   }
}

%apply var_array<short> pair_table_like { var_array<short> *pt, var_array<short> *const pt, var_array<short> &pt, var_array<short> const &pt };

