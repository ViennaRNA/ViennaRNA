/**********************************************/
/* BEGIN interface for Heat Capacity          */
/* computation                                */
/**********************************************/

%{

extern "C" {
  typedef struct {
    float temperature;    /**< @brief   The temperature in &deg;C */
    float heat_capacity;  /**< @brief   The specific heat at this temperature in Kcal/(Mol * K) */
  } heat_capacity_result;
}

%}

typedef struct {
  float temperature;    /**< @brief   The temperature in &deg;C */
  float heat_capacity;  /**< @brief   The specific heat at this temperature in Kcal/(Mol * K) */
} heat_capacity_result;


#ifdef SWIGPYTHON
%extend heat_capacity_result {

  std::string
  __str__()
  {
    std::ostringstream out;
    out << "{ temperature: \"" << $self->temperature << "\"";
    out << ", heat_capacity: " << $self->heat_capacity;
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}

}
#endif

namespace std {

%template(HeatCapacityVector) std::vector<heat_capacity_result>;

};

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") heat_capacity;
%feature("kwargs") heat_capacity;
#endif

  std::vector<heat_capacity_result>
  heat_capacity(float         T_min       = 0.,
                float         T_max       = 100.,
                float         T_increment = 1.,
                unsigned int  mpoints     = 2U)
  {
    vrna_heat_capacity_t              *result_c;
    std::vector<heat_capacity_result> result;

    result_c = vrna_heat_capacity($self, T_min, T_max, T_increment, mpoints);

    if (result_c) {
      for (size_t i = 0; result_c[i].temperature >= T_min; i++) {
        heat_capacity_result r;
        r.temperature = result_c[i].temperature;
        r.heat_capacity = result_c[i].heat_capacity;
        result.push_back(r);
      }
    }

    free(result_c);

    return result;
  }
}

%rename (heat_capacity) my_heat_capacity;

#ifdef SWIGPYTHON
%feature("autodoc") heat_capacity;
%feature("kwargs") heat_capacity;
#endif

%{
  std::vector<heat_capacity_result>
  my_heat_capacity(std::string   sequence,
                   float         T_min        = 0.,
                   float         T_max        = 100.,
                   float         T_increment  = 1.,
                   unsigned int  mpoints      = 2U)
  {
    vrna_heat_capacity_t              *result_c;
    std::vector<heat_capacity_result> result;

    result_c = vrna_heat_capacity_simple(sequence.c_str(), T_min, T_max, T_increment, mpoints);

    if (result_c) {
      for (size_t i = 0; result_c[i].temperature >= T_min; i++) {
        heat_capacity_result r;
        r.temperature = result_c[i].temperature;
        r.heat_capacity = result_c[i].heat_capacity;
        result.push_back(r);
      }
    }

    free(result_c);

    return result;
  }
%}

std::vector<heat_capacity_result>
my_heat_capacity(std::string   sequence,
                 float         T_min        = 0.,
                 float         T_max        = 100.,
                 float         T_increment  = 1.,
                 unsigned int  mpoints      = 2U);

%include  <ViennaRNA/heat_capacity.h>
