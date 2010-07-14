#ifndef __VIENNA_RNA_PACKAGE_FOLD_H__
#define __VIENNA_RNA_PACKAGE_FOLD_H__

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/** \file fold.h **/

/** if nonzero use logarithmic ML energy in energy_of_struct **/
extern  int logML;
/** do ML decomposition uniquely (for subopt) **/
extern  int uniq_ML;
/** set to first pos of second seq for cofolding **/
extern  int cut_point;
/** verbose info from energy_of_struct **/
extern  int eos_debug;

/**
*** Compute minimum free energy and an appropriate secondary
*** structure of an RNA sequence
***
*** \see            circfold()
*** \param sequence RNA sequence
*** \param structure a pointer to the character array the
***        secondary structure in dot-bracket notation will be written to
*** \returns the minimum free energy (MFE) in kcal/mol
**/
float fold(const char *sequence, char *structure);

/**
*** Calculate the free energy of an already folded RNA
***
*** \note This function is not entirely threadsafe! Depending on the state of the global
*** variable \ref eos_debug it prints energy information to stdout or not...\n
*** Better use \ref energy_of_structure() instead
***
*** \see              energy_of_structure, energy_of_circ_struct(), energy_of_struct_pt()
*** \param string     RNA sequence
*** \param structure  secondary structure in dot-bracket notation
*** \returns          the free energy of the input structure given the input sequence in kcal/mol
**/
float energy_of_struct(const char *string, const char *structure);

/**
*** Calculate the free energy of an already folded RNA
***
*** If verbosity level is set to a value >0, energies of structure elements are printed to stdout
***
*** \see              energy_of_circ_structure(), energy_of_structure_pt()
*** \param string     RNA sequence
*** \param structure  secondary structure in dot-bracket notation
*** \param verbosity_level a flag to turn verbose output on/off
*** \returns          the free energy of the input structure given the input sequence in kcal/mol
**/
float energy_of_structure(const char *string, const char *structure, int verbosity_level);

/**
*** Calculate the free energy of an already folded RNA
***
*** \note This function is not entirely threadsafe! Depending on the state of the global
*** variable \ref eos_debug it prints energy information to stdout or not...\n
*** Better use \ref energy_of_structure_pt() instead
***
*** \see              make_pair_table(), energy_of_structure()
*** \param string     RNA sequence
*** \param ptable     the pair table of the secondary structure
*** \param s          encoded RNA sequence
*** \param s1         encoded RNA sequence
*** \returns          the free energy of the input structure given the input sequence in 10kcal/mol
**/
int   energy_of_struct_pt(const char *string, short *ptable, short *s, short *s1);

/**
*** Calculate the free energy of an already folded RNA
***
*** If verbosity level is set to a value >0, energies of structure elements are printed to stdout
***
*** \see              make_pair_table(), energy_of_struct()
*** \param string     RNA sequence
*** \param ptable     the pair table of the secondary structure
*** \param s          encoded RNA sequence
*** \param s1         encoded RNA sequence
*** \param verbosity_level a flag to turn verbose output on/off
*** \returns          the free energy of the input structure given the input sequence in 10kcal/mol
**/
int   energy_of_structure_pt(const char *string, short *ptable, short *s, short *s1, int verbosity_level);

/**
*** free arrays for mfe folding
**/
void  free_arrays(void);
/**
*** Allocate arrays for folding\n
*** \deprecated {This function is deprecated and will be removed soon!}
***
**/
void  DEPRECATED(initialize_fold(int length));
/**
*** recalculate energy parameters
**/
void  update_fold_params(void);
/**
***
**/
char  *backtrack_fold_from_pair(char *sequence, int i, int j);
/**
***
**/
int   loop_energy(short * ptable, short *s, short *s1, int i);
/**
***
**/
void  export_fold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);
/**
*** Compute minimum free energy and an appropriate secondary
*** structure of an RNA sequence assuming it to be circular instead of linear
***
*** \see            fold()
*** \param sequence RNA sequence
*** \param structure a pointer to the character array the
***        secondary structure in dot-bracket notation will be written to
*** \returns the minimum free energy (MFE) in kcal/mol
**/
float circfold(const char *string, char *structure);
/**
*** Calculate the free energy of an already folded  circular RNA
***
*** \note This function is not entirely threadsafe! Depending on the state of the global
*** variable \ref eos_debug it prints energy information to stdout or not...\n
*** Better use \ref energy_of_circ_structure() instead
***
*** \see              energy_of_circ_structure(), energy_of_struct(), energy_of_struct_pt()
*** \param string     RNA sequence
*** \param structure  secondary structure in dot-bracket notation
*** \returns          the free energy of the input structure given the input sequence in kcal/mol
**/
float energy_of_circ_struct(const char *string, const char *structure);
/**
*** Calculate the free energy of an already folded  circular RNA
***
*** If verbosity level is set to a value >0, energies of structure elements are printed to stdout
*** \see              energy_of_struct(), energy_of_struct_pt()
*** \param string     RNA sequence
*** \param structure  secondary structure in dot-bracket notation
*** \param verbosity_level a flag to turn verbose output on/off
*** \returns          the free energy of the input structure given the input sequence in kcal/mol
**/
float energy_of_circ_structure(const char *string, const char *structure, int verbosity_level);

/**
***
**/
void  export_circfold_arrays(int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p, int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);


/* finally moved the loop energy function declarations to this header...  */
/* BUT: The functions only exist for backward compatibility reasons!      */
/* You better include "loop_energies.h" and call the functions:           */
/* E_Hairpin() and E_IntLoop() which are (almost) threadsafe as they get  */
/* a pointer to the energy parameter datastructure as additional argument */

/**
*** \deprecated {This function is deprecated and will be removed soon.
*** Use \func E_IntLoop() instead!}
***/
int   DEPRECATED(LoopEnergy(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1));
/**
*** \deprecated {This function is deprecated and will be removed soon.
*** Use \func E_Hairpin() instead!}
***/
int   DEPRECATED(HairpinE(int size, int type, int si1, int sj1, const char *string));

#endif
