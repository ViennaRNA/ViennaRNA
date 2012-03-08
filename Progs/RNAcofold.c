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
#include "PS_dot.h"
#include "cofold.h"
#include "fold.h"
#include "part_func_co.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"
#include "read_epars.h"
#include "params.h"
#include "RNAcofold_cmdl.h"

/*@unused@*/
PRIVATE char rcsid[] = "$Id: RNAcofold.c,v 1.7 2006/05/10 15:14:27 ivo Exp $";

PRIVATE char *costring(char *string);
PRIVATE char *tokenize(char *line);
PRIVATE cofoldF do_partfunc(char *string, int length, int Switch, struct plist **tpr, struct plist **mf, pf_paramT *parameters);
PRIVATE double *read_concentrations(FILE *fp);
PRIVATE void do_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconces);

PRIVATE double bppmThreshold;

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  struct        RNAcofold_args_info args_info;
  unsigned int  input_type;
  char          *string, *input_string;
  char    *structure, *cstruc, *rec_sequence, *orig_sequence, *rec_id, **rec_rest;
  char    fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH];
  char    *ParamFile;
  char    *ns_bases, *c;
  char    *Concfile;
  int     i, length, l, sym, r, cl;
  double  min_en;
  double  kT, sfact, betaScale;
  int     pf, istty;
  int     noconv, noPS;
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
  pf_paramT       *pf_parameters;
  model_detailsT  md;

  set_model_details(&md);

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

  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(RNAcofold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)            temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given)      fold_constrained=1;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)         md.special_hp = tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)         md.dangles = dangles = args_info.dangles_arg;
  /* do not allow weak pairs */
  if(args_info.noLP_given)            md.noLP = noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)            md.noGU = noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given)     md.noGUclosure = no_closingGU = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)          noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given)     energy_set = args_info.energyModel_arg;
  /*  */
  if(args_info.noPS_given)            noPS = 1;
  /* take another energy parameter set */
  if(args_info.paramFile_given)       ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)             ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)         sfact = args_info.pfScale_arg;

  if(args_info.all_pf_given)          doT = pf = 1;
  /* concentrations from stdin */
  if(args_info.concentrations_given)  doC = doT = pf = 1;
  /* set the bppm threshold for the dotplot */
  if(args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0.,args_info.bppmThreshold_arg));
  /* concentrations in file */
  if(args_info.betaScale_given)       betaScale = args_info.betaScale_arg;
  if(args_info.concfile_given){
    Concfile = strdup(args_info.concfile_arg);
    doC = cofi = doT = pf = 1;
  }
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != -1)
      do_backtrack = args_info.partfunc_arg;
  }
  /* free allocated memory of command line data structure */
  RNAcofold_cmdline_parser_free (&args_info);


  /*
  #############################################
  # begin initializing
  #############################################
  */
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards = space(33);
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

  /* print user help if we get input from tty */
  if(istty){
    printf("Use '&' to connect 2 sequences that shall form a complex.\n");
    if(fold_constrained){
      print_tty_constraint(VRNA_CONSTRAINT_DOT | VRNA_CONSTRAINT_X | VRNA_CONSTRAINT_ANG_BRACK | VRNA_CONSTRAINT_RND_BRACK);
      print_tty_input_seq_str("Input sequence (upper or lower case) followed by structure constraint\n");
    }
    else print_tty_input_seq();
  }

  /* set options we wanna pass to read_record */
  if(istty)             read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  if(!fold_constrained) read_opt |= VRNA_INPUT_NO_REST;

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt))
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

    cut_point = -1;

    rec_sequence  = tokenize(rec_sequence); /* frees input_string and sets cut_point */
    length    = (int) strlen(rec_sequence);
    structure = (char *) space((unsigned) length+1);

    /* parse the rest of the current dataset to obtain a structure constraint */
    if(fold_constrained){
      cstruc = NULL;
      int cp = cut_point;
      unsigned int coptions = (rec_id) ? VRNA_CONSTRAINT_MULTILINE : 0;
      coptions |= VRNA_CONSTRAINT_DOT | VRNA_CONSTRAINT_X | VRNA_CONSTRAINT_ANG_BRACK | VRNA_CONSTRAINT_RND_BRACK;
      getConstraint(&cstruc, (const char **)rec_rest, coptions);
      cstruc = tokenize(cstruc);
      if(cut_point != cp) nrerror("cut point in sequence and structure constraint differs");
      cl = (cstruc) ? (int)strlen(cstruc) : 0;

      if(cl == 0)           warn_user("structure constraint is missing");
      else if(cl < length)  warn_user("structure constraint is shorter than sequence");
      else if(cl > length)  nrerror("structure constraint is too long");

      if(cstruc) strncpy(structure, cstruc, sizeof(char)*(cl+1));
    }

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) str_DNA2RNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    str_uppercase(rec_sequence);

    if(istty){
      if (cut_point == -1)
        printf("length = %d\n", length);
      else
        printf("length1 = %d\nlength2 = %d\n", cut_point-1, length-cut_point+1);
    }

    /*
    ########################################################
    # begin actual computations
    ########################################################
    */

    if (doC) {
      FILE *fp;
      if (cofi) { /* read from file */
        fp = fopen(Concfile, "r");
        if (fp==NULL) {
          fprintf(stderr, "could not open concentration file %s", Concfile);
          nrerror("\n");
        }
        ConcAandB = read_concentrations(fp);
        fclose(fp);
      } else {
        printf("Please enter concentrations [mol/l]\n format: ConcA ConcB\n return to end\n");
        ConcAandB = read_concentrations(stdin);
      }
    }
    /*compute mfe of AB dimer*/
    min_en = cofold(rec_sequence, structure);
    assign_plist_from_db(&mfAB, structure, 0.95);

    {
      char *pstring, *pstruct;
      if (cut_point == -1) {
        pstring = strdup(orig_sequence);
        pstruct = strdup(structure);
      } else {
        pstring = costring(orig_sequence);
        pstruct = costring(structure);
      }
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
        if (cut_point >= 0)
          sprintf(annot,
                  "1 %d 9  0 0.9 0.2 omark\n%d %d 9  1 0.1 0.2 omark\n",
                  cut_point-1, cut_point+1, length+1);
        (void) PS_rna_plot_a(pstring, pstruct, ffname, annot, NULL);
      }
      free(pstring);
      free(pstruct);
    }

    if (length>2000)  free_co_arrays();

    /*compute partition function*/
    if (pf) {
      cofoldF AB, AA, BB;
      FLT_OR_DBL *probs;
      if (dangles==1) {
        dangles=2;   /* recompute with dangles as in pf_fold() */
        min_en = energy_of_structure(rec_sequence, structure, 0);
        dangles=1;
      }

      kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      pf_parameters = get_boltzmann_factors(temperature, betaScale, md, pf_scale);

      if (cstruc!=NULL)
        strncpy(structure, cstruc, length+1);
      AB = co_pf_fold_par(rec_sequence, structure, pf_parameters, do_backtrack, fold_constrained);

      if (do_backtrack) {
        char *costruc;
        costruc = (char *) space(sizeof(char)*(strlen(structure)+2));
        if (cut_point<0) printf("%s", structure);
        else {
          strncpy(costruc, structure, cut_point-1);
          strcat(costruc, "&");
          strcat(costruc, structure+cut_point-1);
          printf("%s", costruc);
        }
        if (!istty) printf(" [%6.2f]\n", AB.FAB);
        else printf("\n");/*8.6.04*/
      }
      if ((istty)||(!do_backtrack))
        printf(" free energy of ensemble = %6.2f kcal/mol\n", AB.FAB);
      printf(" frequency of mfe structure in ensemble %g",
             exp((AB.FAB-min_en)/kT));

      printf(" , delta G binding=%6.2f\n", AB.FcAB - AB.FA - AB.FB);

      probs = export_co_bppm();
      assign_plist_from_pr(&prAB, probs, length, bppmThreshold);

      /* if (doQ) make_probsum(length,fname); */ /*compute prob of base paired*/
      /* free_co_arrays(); */
      if (doT) { /* cofold of all dimers, monomers */
        int Blength, Alength;
        char  *Astring, *Bstring, *orig_Astring, *orig_Bstring;
        char *Newstring;
        char Newname[30];
        char comment[80];
        if (cut_point<0) {
          printf("Sorry, i cannot do that with only one molecule, please give me two or leave it\n");
          free(mfAB);
          free(prAB);
          continue;
        }
        if (dangles==1) dangles=2;
        Alength=cut_point-1;        /*length of first molecule*/
        Blength=length-cut_point+1; /*length of 2nd molecule*/

        Astring=(char *)space(sizeof(char)*(Alength+1));/*Sequence of first molecule*/
        Bstring=(char *)space(sizeof(char)*(Blength+1));/*Sequence of second molecule*/
        strncat(Astring,rec_sequence,Alength);
        strncat(Bstring,rec_sequence+Alength,Blength);

        orig_Astring=(char *)space(sizeof(char)*(Alength+1));/*Sequence of first molecule*/
        orig_Bstring=(char *)space(sizeof(char)*(Blength+1));/*Sequence of second molecule*/
        strncat(orig_Astring,orig_sequence,Alength);
        strncat(orig_Bstring,orig_sequence+Alength,Blength);

        /* compute AA dimer */
        AA=do_partfunc(Astring, Alength, 2, &prAA, &mfAA, pf_parameters);
        /* compute BB dimer */
        BB=do_partfunc(Bstring, Blength, 2, &prBB, &mfBB, pf_parameters);
        /*free_co_pf_arrays();*/

        /* compute A monomer */
        do_partfunc(Astring, Alength, 1, &prA, &mfA, pf_parameters);

        /* compute B monomer */
        do_partfunc(Bstring, Blength, 1, &prB, &mfB, pf_parameters);

        compute_probabilities(AB.F0AB, AB.FA, AB.FB, prAB, prA, prB, Alength);
        compute_probabilities(AA.F0AB, AA.FA, AA.FA, prAA, prA, prA, Alength);
        compute_probabilities(BB.F0AB, BB.FA, BB.FA, prBB, prA, prB, Blength);
        printf("Free Energies:\nAB\t\tAA\t\tBB\t\tA\t\tB\n%.6f\t%6f\t%6f\t%6f\t%6f\n",
               AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB);

        if (doC) {
          do_concentrations(AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB, ConcAandB);
          free(ConcAandB);/*freeen*/
        }

        if (fname[0]!='\0') {
          strcpy(ffname, fname);
          strcat(ffname, "_dp5.ps");
        } else strcpy(ffname, "dot5.ps");
        /*output of the 5 dot plots*/

        /*AB dot_plot*/
        /*write Free Energy into comment*/
        sprintf(comment,"\n%%Heterodimer AB FreeEnergy= %.9f\n", AB.FcAB);
        /*reset cut_point*/
        cut_point=Alength+1;
        /*write New name*/
        strcpy(Newname,"AB");
        strcat(Newname,ffname);
        (void)PS_dot_plot_list(orig_sequence, Newname, prAB, mfAB, comment);

        /*AA dot_plot*/
        sprintf(comment,"\n%%Homodimer AA FreeEnergy= %.9f\n",AA.FcAB);
        /*write New name*/
        strcpy(Newname,"AA");
        strcat(Newname,ffname);
        /*write AA sequence*/
        Newstring=(char*)space((2*Alength+1)*sizeof(char));
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
        Newstring=(char*)space((2*Blength+1)*sizeof(char));
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
        if (pf) {          (void) PS_dot_plot_list(rec_sequence, ffname, prAB, mfAB, "doof");
        free(prAB);}
        free(mfAB);
      }
    }
    if (!doT) free_co_pf_arrays();

    (void) fflush(stdout);
    
    /* clean up */
    if(cstruc) free(cstruc);
    if(rec_id) free(rec_id);
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
        print_tty_constraint(VRNA_CONSTRAINT_DOT | VRNA_CONSTRAINT_X | VRNA_CONSTRAINT_ANG_BRACK | VRNA_CONSTRAINT_RND_BRACK);
        print_tty_input_seq_str("Input sequence (upper or lower case) followed by structure constraint\n");
      }
      else print_tty_input_seq();
    }
  }
  return EXIT_SUCCESS;
}

PRIVATE char *tokenize(char *line)
{
  char *pos, *copy;
  int cut = -1;

  copy = (char *) space(strlen(line)+1);
  (void) sscanf(line, "%s", copy);
  pos = strchr(copy, '&');
  if (pos) {
    cut = (int) (pos-copy)+1;
    if (cut >= strlen(copy)) cut = -1;
    if (strchr(pos+1, '&')) nrerror("more than one cut-point in input");
    for (;*pos;pos++) *pos = *(pos+1); /* splice out the & */
  }
  if (cut > -1) {
    if (cut_point==-1) cut_point = cut;
    else if (cut_point != cut) {
      fprintf(stderr,"cut_point = %d cut = %d\n", cut_point, cut);
      nrerror("Sequence and Structure have different cut points.");
    }
  }
  free(line);
  return copy;
}

PRIVATE char *costring(char *string)
{
  char *ctmp;
  int len;

  len = strlen(string);
  ctmp = (char *)space((len+2) * sizeof(char));
  /* first sequence */
  (void) strncpy(ctmp, string, cut_point-1);
  /* spacer */
  ctmp[cut_point-1] = '&';
  /* second sequence */
  (void) strcat(ctmp, string+cut_point-1);

  return ctmp;
}

PRIVATE cofoldF do_partfunc(char *string, int length, int Switch, struct plist **tpr, struct plist **mfpl, pf_paramT *parameters){
  /*compute mfe and partition function of dimere or  monomer*/
  char *Newstring;
  char *tempstruc;
  double min_en;
  double sfact=1.07;
  double kT;
  pf_paramT *par;
  FLT_OR_DBL *probs;
  cofoldF X;
  kT = parameters->kT/1000.;
  switch (Switch)
    {
    case 1: /*monomer*/
      cut_point=-1;
      tempstruc = (char *) space((unsigned)length+1);
      min_en = fold(string, tempstruc);
      assign_plist_from_db(mfpl, tempstruc, 0.95);
      free_arrays();
      /*En=pf_fold(string, tempstruc);*/
      /* init_co_pf_fold(length); <- obsolete */
      par = get_boltzmann_factor_copy(parameters);
      par->pf_scale = exp(-(sfact*min_en)/kT/(length));
      X=co_pf_fold_par(string, tempstruc, par, 1, fold_constrained);
      probs = export_co_bppm();
      assign_plist_from_pr(tpr, probs, length, bppmThreshold);
      free_co_pf_arrays();
      free(tempstruc);
      free(par);
      break;

    case 2: /*dimer*/
      Newstring=(char *)space(sizeof(char)*(length*2+1));
      strcat(Newstring,string);
      strcat(Newstring,string);
      cut_point=length+1;
      tempstruc = (char *) space((unsigned)length*2+1);
      min_en = cofold(Newstring, tempstruc);
      assign_plist_from_db(mfpl, tempstruc, 0.95);
      free_co_arrays();
      /* init_co_pf_fold(2*length); <- obsolete */
      par = get_boltzmann_factor_copy(parameters);
      par->pf_scale =exp(-(sfact*min_en)/kT/(2*length));
      X=co_pf_fold_par(Newstring, tempstruc, par, 1, fold_constrained);
      probs = export_co_bppm();
      assign_plist_from_pr(tpr, probs, 2*length, bppmThreshold);
      free_co_pf_arrays();
      free(Newstring);
      free(tempstruc);
      free(par);
      break;

    default:
      printf("Error in get_partfunc\n, computing neither mono- nor dimere!\n");
      exit (42);

    }
  return X;
}


PRIVATE void do_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc) {
  /* compute and print concentrations out of free energies, calls get_concentrations */
  struct ConcEnt *result;
  int i, n;

  result=get_concentrations(FEAB, FEAA, FEBB, FEA, FEB, startconc);

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

  startc = (double *) space((2*n+1)*sizeof(double));

  while ((line=get_line(fp))!=NULL) {
    int c;
    if (i+4>=2*n) {
      n*=2;
      startc=(double *)xrealloc(startc,(2*n+1)*sizeof(double));
    }
    c = sscanf(line,"%lf %lf", &startc[i], &startc[i+1]);
    free(line);
    if (c<2) break;
    i+=2;
  }
  startc[i]=startc[i+1]=0;
  return startc;
}
