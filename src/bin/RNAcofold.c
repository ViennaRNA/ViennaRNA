/* Last changed Time-stamp: <2006-04-05 12:58:49 ivo> */
/*
                  c Ivo L Hofacker, Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func_co.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "RNAcofold_cmdl.h"

/*@unused@*/
PRIVATE char rcsid[] = "$Id: RNAcofold.c,v 1.7 2006/05/10 15:14:27 ivo Exp $";

PRIVATE vrna_dimer_pf_t do_partfunc(char *string, int length, int Switch, plist **tpr, plist **mf, vrna_exp_param_t *parameters);
PRIVATE double *read_concentrations(FILE *fp);
PRIVATE void do_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconces, vrna_exp_param_t *parameters);

PRIVATE double bppmThreshold;

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  struct        RNAcofold_args_info args_info;
  unsigned int  input_type;
  char          *string, *input_string;
  char    *constraints_file, *structure, *cstruc, *rec_sequence, *orig_sequence, *rec_id, **rec_rest;
  char    fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH];
  char    *ParamFile;
  char    *ns_bases, *c;
  char    *Concfile;
  int     i, length, l, sym, r, cl;
  double  min_en;
  double  kT, sfact, betaScale;
  int     pf, istty;
  int     noconv, noPS, enforceConstraints;
  int     doT;    /*compute dimere free energies etc.*/
  int     doC;    /*toggle to compute concentrations*/
  int     doQ;    /*toggle to compute prob of base being paired*/
  int     cofi;   /*toggle concentrations stdin / file*/
  plist   *prAB;
  plist   *prAA;   /*pair probabilities of AA dimer*/
  plist   *prBB;
  plist   *prA;
  plist   *prB;
  plist   *mfAB;
  plist   *mfAA;   /*pair mfobabilities of AA dimer*/
  plist   *mfBB;
  plist   *mfA;
  plist   *mfB;
  double  *ConcAandB;
  unsigned int    rec_type, read_opt;
  vrna_param_t      *P;
  vrna_exp_param_t  *pf_parameters;
  vrna_md_t         md;


  /*
  #############################################
  # init variables and parameter options
  #############################################
  */
  dangles       = 2;
  sfact         = 1.07;
  bppmThreshold = 1e-5;
  noconv        = 0;
  noPS          = 0;
  do_backtrack  = 1;
  pf            = 0;
  doT           = 0;
  doC           = 0;
  doQ           = 0;
  cofi          = 0;
  betaScale     = 1.;
  gquad         = 0;
  ParamFile     = NULL;
  pf_parameters = NULL;
  string        = NULL;
  Concfile      = NULL;
  structure     = NULL;
  cstruc        = NULL;
  ns_bases      = NULL;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  constraints_file = NULL;
  enforceConstraints = 0;

  set_model_details(&md);
  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(RNAcofold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given){
    fold_constrained=1;
    if(args_info.constraint_arg[0] != '\0')
      constraints_file = strdup(args_info.constraint_arg);
  }
  /* enforce base pairs given in constraint string rather than weak enforce */
  if(args_info.enforceConstraint_given)
    enforceConstraints = 1;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
     md.dangles = dangles = args_info.dangles_arg;
  }
  /* do not allow weak pairs */
  if(args_info.noLP_given)
    md.noLP = noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)
    md.noGU = noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;
  /* gquadruplex support */
  if(args_info.gquad_given)
    md.gquad = gquad = 1;
  /* enforce canonical base pairs in any case? */
  if(args_info.canonicalBPonly_given)
    md.canonicalBPonly = canonicalBPonly = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)
    noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;
  /*  */
  if(args_info.noPS_given)
    noPS = 1;
  /* take another energy parameter set */
  if(args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)
    md.sfact = sfact = args_info.pfScale_arg;

  if(args_info.all_pf_given)
    doT = pf = 1;
  /* concentrations from stdin */
  if(args_info.concentrations_given)
    doC = doT = pf = 1;
  /* set the bppm threshold for the dotplot */
  if(args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0.,args_info.bppmThreshold_arg));
  /* concentrations in file */
  if(args_info.betaScale_given)
    md.betaScale = betaScale = args_info.betaScale_arg;
  if(args_info.concfile_given){
    Concfile = strdup(args_info.concfile_arg);
    doC = cofi = doT = pf = 1;
  }
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != -1)
      md.compute_bpp = do_backtrack = args_info.partfunc_arg;
  }
  /* free allocated memory of command line data structure */
  RNAcofold_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  if(pf && gquad){
    vrna_message_error("G-Quadruplex support is currently not available for partition function computations");
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards = vrna_alloc(33);
    c=ns_bases;
    i=sym=0;
    if (*c=='-') {
      sym=1; c++;
    }
    while (*c!='\0') {
      if (*c!=',') {
        nonstandards[i++]=*c++;
        nonstandards[i++]=*c;
        if ((sym)&&(*c!=*(c-1))) {
          nonstandards[i++]=*c;
          nonstandards[i++]=*(c-1);
        }
      }
      c++;
    }
  }
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));

  /* get energy parameters */
  P = vrna_params(&md);

  /* print user help if we get input from tty */
  if(istty){
    printf("Use '&' to connect 2 sequences that shall form a complex.\n");
    if(fold_constrained){
      vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK);
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n");
    }
    else vrna_message_input_seq_simple();
  }

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if(istty)             read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  if(!fold_constrained) read_opt |= VRNA_INPUT_NO_REST;

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    if(rec_id){
      if(!istty) printf("%s\n", rec_id);
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) vrna_seq_toRNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound_t *vc = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE |  VRNA_OPTION_HYBRID | ((pf) ? VRNA_OPTION_PF : 0));
    length    = vc->length;
    structure = (char *) vrna_alloc((unsigned) length+1);

    /* parse the rest of the current dataset to obtain a structure constraint */
    if(fold_constrained){
      if(constraints_file){
        vrna_constraints_add(vc, constraints_file, VRNA_CONSTRAINT_FILE | VRNA_CONSTRAINT_SOFT_MFE | ((pf) ? VRNA_CONSTRAINT_SOFT_PF : 0));
      } else {
        cstruc = NULL;
        cstruc = NULL;
        int cp = -1;
        unsigned int coptions = (rec_id) ? VRNA_CONSTRAINT_MULTILINE : 0;
        coptions |= VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK | VRNA_CONSTRAINT_DB_PIPE;
        vrna_extract_record_rest_constraint(&cstruc, (const char **)rec_rest, coptions);
        cstruc = vrna_cut_point_remove(cstruc, &cp);
        if(vc->cutpoint != cp){
          fprintf(stderr,"cut_point = %d cut = %d\n", vc->cutpoint, cp);
          vrna_message_error("Sequence and Structure have different cut points.");
        }

        cl = (cstruc) ? (int)strlen(cstruc) : 0;

        if(cl == 0)           vrna_message_warning("structure constraint is missing");
        else if(cl < length)  vrna_message_warning("structure constraint is shorter than sequence");
        else if(cl > length)  vrna_message_error("structure constraint is too long");

        if(cstruc){
          strncpy(structure, cstruc, sizeof(char)*(cl+1));

          unsigned int constraint_options = 0;
          constraint_options |= VRNA_CONSTRAINT_DB
                                | VRNA_CONSTRAINT_DB_PIPE
                                | VRNA_CONSTRAINT_DB_DOT
                                | VRNA_CONSTRAINT_DB_X
                                | VRNA_CONSTRAINT_DB_ANG_BRACK
                                | VRNA_CONSTRAINT_DB_RND_BRACK;

          if(enforceConstraints)
            constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;
          vrna_constraints_add(vc, (const char *)structure, constraint_options);
        }
      }
    }

    if(istty){
      if (cut_point == -1)
        printf("length = %d\n", length);
      else
        printf("length1 = %d\nlength2 = %d\n", cut_point-1, length-cut_point+1);
    }

    if (doC) {
      FILE *fp;
      if (cofi) { /* read from file */
        fp = fopen(Concfile, "r");
        if (fp==NULL) {
          fprintf(stderr, "could not open concentration file %s", Concfile);
          vrna_message_error("\n");
        }
        ConcAandB = read_concentrations(fp);
        fclose(fp);
      } else {
        printf("Please enter concentrations [mol/l]\n format: ConcA ConcB\n return to end\n");
        ConcAandB = read_concentrations(stdin);
      }
    }
    /*
    ########################################################
    # begin actual computations
    ########################################################
    */

    /* compute mfe of AB dimer */
    min_en  = vrna_mfe_dimer(vc, structure);
    mfAB    = vrna_plist(structure, 0.95);

    {
      char *pstring, *pstruct;
      pstring = strdup(orig_sequence);
      pstruct = vrna_cut_point_insert(structure, vc->cutpoint);

      printf("%s\n%s", pstring, pstruct);
      if (istty)
        printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
      else
        printf(" (%6.2f)\n", min_en);

      (void) fflush(stdout);

      if (!noPS) {
        char annot[512] = "";
        if (fname[0]!='\0') {
          strcpy(ffname, fname);
          strcat(ffname, "_ss.ps");
        } else {
          strcpy(ffname, "rna.ps");
        }
        if (vc->cutpoint >= 0)
          sprintf(annot,
                  "1 %d 9  0 0.9 0.2 omark\n%d %d 9  1 0.1 0.2 omark\n",
                  vc->cutpoint-1, vc->cutpoint+1, length+1);
        if (!noPS)
          (void) vrna_file_PS_rnaplot_a(pstring, pstruct, ffname, annot, NULL, &md);
      }
      free(pstring);
      free(pstruct);
    }

    if (length>2000)
      vrna_mx_mfe_free(vc);

    /* compute partition function */
    if (pf) {
      vrna_dimer_pf_t AB, AA, BB;
      if (dangles==1){
        vc->params->model_details.dangles = dangles = 2;   /* recompute with dangles as in pf_fold() */
        min_en = vrna_eval_structure(vc, structure);
        vc->params->model_details.dangles = dangles = 1;
      }

      vrna_exp_params_rescale(vc, &min_en);
      kT = vc->exp_params->kT/1000.;

      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      /* do we need to add hard constraints? */
      if (cstruc!=NULL) strncpy(structure, cstruc, length+1);

      /* compute partition function */
      AB = vrna_pf_dimer(vc, structure);

      if (do_backtrack) {
        char *costruc;
        costruc = vrna_cut_point_insert(structure, vc->cutpoint);
        printf("%s", costruc);

        if (!istty) printf(" [%6.2f]\n", AB.FAB);
        else printf("\n");/*8.6.04*/
        free(costruc);
        prAB = vrna_plist_from_probs(vc, bppmThreshold);
      }

      if ((istty)||(!do_backtrack))
        printf(" free energy of ensemble = %6.2f kcal/mol\n", AB.FAB);

      printf(" frequency of mfe structure in ensemble %g",
             exp((AB.FAB-min_en)/kT));

      printf(" , delta G binding=%6.2f\n", AB.FcAB - AB.FA - AB.FB);

      /* if (doQ) make_probsum(length,fname); */ /*compute prob of base paired*/
      /* free_co_arrays(); */
      if (doT) { /* cofold of all dimers, monomers */
        int Blength, Alength;
        char  *Astring, *Bstring, *orig_Astring, *orig_Bstring;
        char *Newstring;
        char Newname[30];
        char comment[80];
        if (vc->cutpoint <= 0) {
          printf("Sorry, i cannot do that with only one molecule, please give me two or leave it\n");
          free(mfAB);
          free(prAB);
          continue;
        }
        if (dangles==1) dangles=2;
        Alength = vc->cutpoint - 1;           /* length of first molecule */
        Blength = length - vc->cutpoint + 1;  /* length of 2nd molecule   */

        Astring = (char *)vrna_alloc(sizeof(char)*(Alength+1));/*Sequence of first molecule*/
        Bstring = (char *)vrna_alloc(sizeof(char)*(Blength+1));/*Sequence of second molecule*/
        strncat(Astring,rec_sequence,Alength);
        strncat(Bstring,rec_sequence+Alength+1,Blength);

        orig_Astring=(char *)vrna_alloc(sizeof(char)*(Alength+1));/*Sequence of first molecule*/
        orig_Bstring=(char *)vrna_alloc(sizeof(char)*(Blength+1));/*Sequence of second molecule*/
        strncat(orig_Astring,orig_sequence,Alength);
        strncat(orig_Bstring,orig_sequence+Alength+1,Blength);

        /* compute AA dimer */
        AA=do_partfunc(Astring, Alength, 2, &prAA, &mfAA, pf_parameters);
        /* compute BB dimer */
        BB=do_partfunc(Bstring, Blength, 2, &prBB, &mfBB, pf_parameters);
        /*free_co_pf_arrays();*/

        /* compute A monomer */
        do_partfunc(Astring, Alength, 1, &prA, &mfA, pf_parameters);

        /* compute B monomer */
        do_partfunc(Bstring, Blength, 1, &prB, &mfB, pf_parameters);

        if(do_backtrack){
          vrna_pf_dimer_probs(AB.F0AB, AB.FA, AB.FB, prAB, prA, prB, Alength, pf_parameters);
          vrna_pf_dimer_probs(AA.F0AB, AA.FA, AA.FA, prAA, prA, prA, Alength, pf_parameters);
          vrna_pf_dimer_probs(BB.F0AB, BB.FA, BB.FA, prBB, prA, prB, Blength, pf_parameters);
        }
        printf("Free Energies:\nAB\t\tAA\t\tBB\t\tA\t\tB\n%.6f\t%6f\t%6f\t%6f\t%6f\n",
               AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB);

        if (doC) {
          vrna_pf_dimer_concentrations(AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB, ConcAandB, pf_parameters);
          free(ConcAandB);/*freeen*/
        }

        if (fname[0]!='\0') {
          strcpy(ffname, fname);
          strcat(ffname, "_dp5.ps");
        } else strcpy(ffname, "dot5.ps");
        /*output of the 5 dot plots*/

        if(do_backtrack){

          /*AB dot_plot*/

          /*write Free Energy into comment*/
          sprintf(comment,"\n%%Heterodimer AB FreeEnergy= %.9f\n", AB.FcAB);
          /*reset cut_point*/
          cut_point=Alength+1;

          /*write New name*/
          strcpy(Newname,"AB");
          strcat(Newname,ffname);
          int cp = -1;
          char *coseq = vrna_cut_point_remove(orig_sequence, &cp);
          (void)PS_dot_plot_list(coseq, Newname, prAB, mfAB, comment);
          free(coseq);

          /*AA dot_plot*/
          sprintf(comment,"\n%%Homodimer AA FreeEnergy= %.9f\n",AA.FcAB);
          /*write New name*/
          strcpy(Newname,"AA");
          strcat(Newname,ffname);
          /*write AA sequence*/
          Newstring=(char*)vrna_alloc((2*Alength+1)*sizeof(char));
          strcpy(Newstring,orig_Astring);
          strcat(Newstring,orig_Astring);
          (void)PS_dot_plot_list(Newstring, Newname, prAA, mfAA, comment);
          free(Newstring);

          /*BB dot_plot*/
          sprintf(comment,"\n%%Homodimer BB FreeEnergy= %.9f\n",BB.FcAB);
          /*write New name*/
          strcpy(Newname,"BB");
          strcat(Newname,ffname);
          /*write BB sequence*/
          Newstring=(char*)vrna_alloc((2*Blength+1)*sizeof(char));
          strcpy(Newstring,orig_Bstring);
          strcat(Newstring,orig_Bstring);
          /*reset cut_point*/
          cut_point=Blength+1;
          (void)PS_dot_plot_list(Newstring, Newname, prBB, mfBB, comment);
          free(Newstring);

          /*A dot plot*/
          /*reset cut_point*/
          cut_point=-1;
          sprintf(comment,"\n%%Monomer A FreeEnergy= %.9f\n",AB.FA);
          /*write New name*/
          strcpy(Newname,"A");
          strcat(Newname,ffname);
          /*write BB sequence*/
          (void)PS_dot_plot_list(orig_Astring, Newname, prA, mfA, comment);

          /*B monomer dot plot*/
          sprintf(comment,"\n%%Monomer B FreeEnergy= %.9f\n",AB.FB);
          /*write New name*/
          strcpy(Newname,"B");
          strcat(Newname,ffname);
          /*write BB sequence*/
          (void)PS_dot_plot_list(orig_Bstring, Newname, prB, mfB, comment);
        }

        free(Astring); free(Bstring); free(orig_Astring); free(orig_Bstring);
        free(prAB); free(prAA); free(prBB); free(prA); free(prB);
        free(mfAB); free(mfAA); free(mfBB); free(mfA); free(mfB);

      } /*end if(doT)*/

      free(pf_parameters);
    }/*end if(pf)*/


    if (do_backtrack) {
      if (fname[0]!='\0') {
        strcpy(ffname, fname);
        strcat(ffname, "_dp.ps");
      } else strcpy(ffname, "dot.ps");

      if (!doT) {
        if (pf) {
          int cp;
          char *seq = vrna_cut_point_remove(rec_sequence, &cp);
          (void) vrna_plot_dp_PS_list(seq, cp, ffname, prAB, mfAB, "doof");
          free(prAB);
          free(seq);
        }
        free(mfAB);
      }
    }
    if (!doT)
      vrna_mx_pf_free(vc);

    (void) fflush(stdout);
    
    /* clean up */
    free(cstruc);
    free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    free(structure);
    /* free the rest of current dataset */
    if(rec_rest){
      for(i=0;rec_rest[i];i++) free(rec_rest[i]);
      free(rec_rest);
    }
    rec_id = rec_sequence = orig_sequence = structure = cstruc = NULL;
    rec_rest = NULL;

    /* print user help for the next round if we get input from tty */
    if(istty){
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      if(fold_constrained){
        vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK);
        vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n");
      }
      else vrna_message_input_seq_simple();
    }
  }
  return EXIT_SUCCESS;
}

PRIVATE vrna_dimer_pf_t
do_partfunc(char *string,
            int length,
            int Switch,
            plist **tpr,
            plist **mfpl,
            vrna_exp_param_t *parameters){

  /*compute mfe and partition function of dimer or monomer*/
  char *Newstring;
  char *tempstruc;
  double min_en;
  double sfact=1.07;
  double kT;
  vrna_exp_param_t *par;
  vrna_dimer_pf_t X;
  vrna_fold_compound_t *vc;
  kT = parameters->kT/1000.;
  switch (Switch){
    case 1:   /* monomer */
              tempstruc = (char *) vrna_alloc((unsigned)length+1);
              //parameters->model_details.min_loop_size = TURN; /* we need min_loop_size of 0 to correct for Q_AB */
              vc = vrna_fold_compound(string, &(parameters->model_details), VRNA_OPTION_MFE | VRNA_OPTION_PF);
              min_en = vrna_mfe(vc, tempstruc);
              *mfpl = vrna_plist(tempstruc, 0.95);
              vrna_mx_mfe_free(vc);

              par = vrna_exp_params_copy(parameters);
              par->pf_scale = exp(-(sfact*min_en)/kT/(length));
              vrna_exp_params_subst(vc, par);
              X = vrna_pf_dimer(vc, tempstruc);
              if(*tpr){
                *tpr = vrna_plist_from_probs(vc, bppmThreshold);
              }
              vrna_fold_compound_free(vc);
              free(tempstruc);
              free(par);
              parameters->model_details.min_loop_size = 0;
              break;

    case 2:   /* dimer */
              tempstruc = (char *) vrna_alloc((unsigned)length*2+2);
              Newstring = (char *)vrna_alloc(sizeof(char)*(length*2+2));
              strcat(Newstring, string); strcat(Newstring, "&"); strcat(Newstring, string);
              parameters->model_details.min_loop_size = 0;
              vc = vrna_fold_compound(Newstring, &(parameters->model_details), VRNA_OPTION_MFE | VRNA_OPTION_PF | VRNA_OPTION_HYBRID);
              min_en = vrna_mfe_dimer(vc, tempstruc);
              *mfpl = vrna_plist(tempstruc, 0.95);
              vrna_mx_mfe_free(vc);

              par = vrna_exp_params_copy(parameters);
              par->pf_scale = exp(-(sfact*min_en)/kT/(2*length));
              vrna_exp_params_subst(vc, par);
              X = vrna_pf_dimer(vc, tempstruc);
              if(*tpr){
                *tpr = vrna_plist_from_probs(vc, bppmThreshold);
              }
              vrna_fold_compound_free(vc);

              free(Newstring);
              free(tempstruc);
              free(par);
              break;

    default:  printf("Error in get_partfunc\n, computing neither mono- nor dimere!\n");
              exit (42);
  }

  return X;
}


PRIVATE void
do_concentrations(double FEAB,
                  double FEAA,
                  double FEBB,
                  double FEA,
                  double FEB,
                  double *startconc,
                  vrna_exp_param_t *parameters){

  /* compute and print concentrations out of free energies, calls get_concentrations */
  vrna_dimer_conc_t *result;
  int i, n;

  result=vrna_pf_dimer_concentrations(FEAB, FEAA, FEBB, FEA, FEB, startconc, parameters);

  printf("Initial concentrations\t\trelative Equilibrium concentrations\n");
  printf("A\t\t B\t\t AB\t\t AA\t\t BB\t\t A\t\t B\n");
  for (n=0; (startconc[2*n]>0) || (startconc[2*n+1]>0); n++); /* count */
  for (i=0; i<n ;i++) {
    double tot = result[i].A0+result[i].B0;
    printf("%-10g\t%-10g\t%.5f \t%.5f \t%.5f \t%.5f \t%.5f\n", result[i].A0, result[i].B0,
           result[i].ABc/tot, result[i].AAc/tot, result[i].BBc/tot ,result[i].Ac/tot, result[i].Bc/tot);
  }
  free(result);
}


PRIVATE double *read_concentrations(FILE *fp) {
  /* reads concentrations, returns list of double, -1. marks end */
  char *line;
  double *startc;
  int i=0, n=2;

  startc = (double *) vrna_alloc((2*n+1)*sizeof(double));

  while ((line=get_line(fp))!=NULL) {
    int c;
    if (i+4>=2*n) {
      n*=2;
      startc=(double *)vrna_realloc(startc,(2*n+1)*sizeof(double));
    }
    c = sscanf(line,"%lf %lf", &startc[i], &startc[i+1]);
    free(line);
    if (c<2) break;
    i+=2;
  }
  startc[i]=startc[i+1]=0;
  return startc;
}
