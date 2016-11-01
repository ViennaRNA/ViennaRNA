/*
  Last changed Time-stamp: <2015-04-08 15:37:56 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: baum.c,v 1.9 2008/05/21 10:15:45 ivo Exp $
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if HAVE_LIBRNA_API3
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/cofold.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/pair_mat.h>
#else
#include <fold_vars.h>
#include <fold.h>
#include <cofold.h>
#include <utils.h>
#include <pair_mat.h>
#endif

#include "nachbar.h"
#include "globals.h"

#define MYTURN 1
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))
#define ORDER(x,y) if ((x)->nummer>(y)->nummer) {tempb=x; x=y; y=tempb;}

/* item of structure ringlist */
typedef struct _baum {
  int nummer; /* number of base in sequence */
  char typ;   /* 'r' virtualroot, 'p' or 'q' paired, 'u' unpaired */
  unsigned short base; /* 0<->unknown, 1<->A, 2<->C, 3<->G, 4<->U */
  int loop_energy;
  struct _baum *up;
  struct _baum *next;
  struct _baum *prev;
  struct _baum *down;
} baum;

static char UNUSED rcsid[]="$Id: baum.c,v 1.9 2008/05/21 10:15:45 ivo Exp $";
static short *pairList = NULL;
static short *typeList = NULL;
static short *aliasList = NULL;
static baum *rl = NULL;         /* ringlist */
static baum *wurzl = NULL;      /* virtualroot of ringlist-tree */
static char **ptype = NULL;

static int comp_struc(const void *A, const void *B);
/* PUBLIC FUNCTIONES */
void ini_or_reset_rl (void);
void move_it (void);
void update_tree (int i, int j);
void clean_up_rl (void);

/* PRIVATE FUNCTIONES */
static void ini_ringlist(void);
static void reset_ringlist(void);
static void struc2tree (char *struc);
static void close_bp_en (baum *i, baum *j);
static void close_bp (baum *i, baum *j);
static void open_bp (baum *i);
static void open_bp_en (baum *i);
static void inb (baum *root);
static void inb_nolp (baum *root);
static void dnb (baum *rli);
static void dnb_nolp (baum *rli);
static void fnb (baum *rli);
static void make_ptypes(const short *S);
/* debugging tool(s) */
#if 0
static void rl_status(void);
#endif

/* convert structure in bracked-dot-notation to a ringlist-tree */
static void struc2tree(char *struc) {
  char* struc_copy;
  int ipos, jpos, balance = 0;
  baum *rli, *rlj;

  struc_copy = (char *)calloc(GSV.len+1, sizeof(char));
  assert(struc_copy);
  strcpy(struc_copy,struc);

  for (ipos = 0; ipos < GSV.len; ipos++) {
    if (struc_copy[ipos] == ')') {
      jpos = ipos;
      struc_copy[ipos] = '.';
      balance++;
      while (struc_copy[--ipos] != '(');
      struc_copy[ipos] = '.';
      balance--;
      rli = &rl[ipos];
      rlj = &rl[jpos];
      close_bp(rli, rlj);
    }
  }

  if (balance) {
    fprintf(stderr,
	    "struc2tree(): start structure is not balanced !\n%s\n%s\n",
	    GAV.farbe, struc);
    exit(1);
  }

#if HAVE_LIBRNA_API3
  GSV.currE = GSV.startE = (float)vrna_eval_structure_pt(GAV.vc, pairList) / 100.0;
#else
  GSV.currE = GSV.startE =
    (float )energy_of_struct_pt_par(GAV.farbe, pairList, typeList,
				    aliasList, GAV.params, 0) / 100.0;
#endif
  {
    int i;
    for(i = 0; i < GSV.len; i++) {
      if (pairList[i+1]>i+1)
#if HAVE_LIBRNA_API3
        rl[i].loop_energy = vrna_eval_loop_pt(GAV.vc, i+1, pairList);
#else
	rl[i].loop_energy = loop_energy(pairList, typeList, aliasList,i+1);
#endif
    }
#if HAVE_LIBRNA_API3
    wurzl->loop_energy = vrna_eval_loop_pt(GAV.vc, 0, pairList);
#else
    wurzl->loop_energy = loop_energy(pairList, typeList, aliasList,0);
#endif
  }

  free(struc_copy);
}

/**/
static void ini_ringlist(void) {
  int i;

  /* needed by function energy_of_struct_pt() from Vienna-RNA-1.4 */
  pairList = (short *)calloc(GSV.len + 2, sizeof(short));
  assert(pairList != NULL);
  typeList = (short *)calloc(GSV.len + 2, sizeof(short));
  assert(typeList != NULL);
  aliasList = (short *)calloc(GSV.len + 2, sizeof(short));
  assert(aliasList != NULL);
  pairList[0] = typeList[0] = aliasList[0] = GSV.len;
  ptype =  (char **)calloc(GSV.len + 2, sizeof(char *));
  assert(ptype != NULL);
  for (i=0; i<=GSV.len; i++) {
    ptype[i] =   (char*)calloc(GSV.len + 2, sizeof(char));
    assert(ptype[i] != NULL);
  }

  /* allocate virtual root */
  wurzl = (baum *)calloc(1, sizeof(baum));
  assert(wurzl != NULL);
  /* allocate ringList */
  rl = (baum *)calloc(GSV.len+1, sizeof(baum));
  assert(rl != NULL);
  /* allocate PostOrderList */

  /* initialize virtualroot */
  wurzl->typ = 'r';
  wurzl->nummer = -1;
  /* connect virtualroot to ringlist-tree in down direction */
  wurzl->down = &rl[GSV.len];
  /* initialize post-order list */

  make_pair_matrix();

  /* initialize rest of ringlist-tree */
  for(i = 0; i < GSV.len; i++) {
    int c;
    GAV.currform[i] = '.';
    GAV.prevform[i] = 'x';
    pairList[i+1] = 0;
    rl[i].typ = 'u';
    /* decode base to numeric value */
    c = encode_char(GAV.farbe[i]);
    rl[i].base = typeList[i+1] = c;
    aliasList[i+1] = alias[typeList[i+1]];
    /* astablish links for node of the ringlist-tree */
    rl[i].nummer = i;
    rl[i].next = &rl[i+1];
    rl[i].prev = ((i == 0) ? &rl[GSV.len] : &rl[i-1]);
    rl[i].up = rl[i].down = NULL;
  }
  GAV.currform[GSV.len] =   GAV.prevform[GSV.len] = '\0';
  make_ptypes(aliasList);

  rl[i].nummer = i;
  rl[i].base = 0;
  /* make ringlist circular in next, prev direction */
  rl[i].next = &rl[0];
  rl[i].prev = &rl[i-1];
  /* make virtual basepair for virtualroot */
  rl[i].up = wurzl;
  rl[i].typ = 'x';

}

/**/
void ini_or_reset_rl(void) {

  /* if there is no ringList-tree make a new one */
  if (wurzl == NULL) {
    ini_ringlist();

    /* start structure */
    struc2tree(GAV.startform);
#if HAVE_LIBRNA_API3
    GSV.currE = GSV.startE = vrna_eval_structure(GAV.vc, GAV.startform);
#else
    GSV.currE = GSV.startE = energy_of_structure(GAV.farbe, GAV.startform, 0);
#endif

    /* stop structure(s) */
    if ( GTV.stop )  {
      int i;

      qsort(GAV.stopform, GSV.maxS, sizeof(char *), comp_struc);
#if HAVE_LIBRNA_API3
      /*
        note that we need to hack the full length into GAV.vc again,
        in case it was shortened due to chain growth simulation
      */
      unsigned int n, tmp_n;
      n     = strlen(GAV.farbe_full);
      tmp_n = GAV.vc->length;
      GAV.vc->length = n;
      for (i = 0; i< GSV.maxS; i++)
        GAV.sE[i] = vrna_eval_structure(GAV.vc, GAV.stopform[i]);
      GAV.vc->length = tmp_n;
#else
      for (i = 0; i< GSV.maxS; i++)
	GAV.sE[i] = energy_of_structure(GAV.farbe_full, GAV.stopform[i], 0);
#endif
    }
    else {
#if HAVE_LIBRNA_API3
      /* fold sequence to get Minimum free energy structure (Mfe) */
      /*
        note that we need to hack the full length into GAV.vc again,
        in case it was shortened due to chain growth simulation
      */
      unsigned int n, tmp_n;
      n     = strlen(GAV.farbe_full);
      tmp_n = GAV.vc->length;
      GAV.vc->length = n;
      GAV.sE[0] = vrna_mfe_dimer(GAV.vc, GAV.stopform[0]);
      vrna_mx_mfe_free(GAV.vc);
      /* revaluate energy of Mfe (maye differ if --logML=logarthmic */
      GAV.sE[0] = vrna_eval_structure(GAV.vc, GAV.stopform[0]);
      GAV.vc->length = tmp_n;
#else
      if(GTV.noLP)
	noLonelyPairs=1;
      initialize_cofold(GSV.len);
      /* fold sequence to get Minimum free energy structure (Mfe) */
      GAV.sE[0] = cofold(GAV.farbe_full, GAV.stopform[0]);
      free_arrays();
      /* revaluate energy of Mfe (maye differ if --logML=logarthmic */
      GAV.sE[0] = energy_of_structure(GAV.farbe_full, GAV.stopform[0], 0);
#endif
    }
    GSV.stopE = GAV.sE[0];
    ini_nbList(strlen(GAV.farbe_full)*strlen(GAV.farbe_full));
  }
  else {
    /* reset ringlist-tree to start conditions */
    reset_ringlist();
    if(GTV.start) struc2tree(GAV.startform);
    else {
      GSV.currE = GSV.startE;
    }
  }
}

/**/
static void reset_ringlist(void) {
  int i;

  for(i = 0; i < GSV.len; i++) {
    GAV.currform[i] = '.';
    GAV.prevform[i] = 'x';
    pairList[i+1] = 0;
    rl[i].typ = 'u';
    rl[i].next = &rl[i + 1];
    rl[i].prev = ((i == 0) ? &rl[GSV.len] : &rl[i - 1]);
    rl[i].up = rl[i].down = NULL;
  }
  rl[i].next = &rl[0];
  rl[i].prev = &rl[i-1];
  rl[i].up = wurzl;
}

/* update ringlist-tree */
void update_tree(int i, int j) {

  baum *rli, *rlj, *tempb;

  if ( abs(i) < GSV.len) { /* >> single basepair move */
    if ((i > 0) && (j > 0)) { /* insert */
      rli = &rl[i-1];
      rlj = &rl[j-1];
      close_bp_en(rli, rlj);
    }
    else if ((i < 0)&&(j < 0)) { /* delete */
      i = -i;
      rli = &rl[i-1];
      open_bp_en(rli);
    }
    else { /* shift */
      if (i > 0) { /* i remains the same, j shifts */
	j=-j;
	rli=&rl[i-1];
	rlj=&rl[j-1];
	open_bp_en(rli);
	ORDER(rli, rlj);
	close_bp_en(rli, rlj);
      }
      else { /* j remains the same, i shifts */
	baum *old_rli;
	i = -i;
	rli = &rl[i-1];
	rlj = &rl[j-1];
	old_rli = rlj->up;
	open_bp_en(old_rli);
	ORDER(rli, rlj);
	close_bp_en(rli, rlj);
      }
    }
  } /* << single basepair move */
  else { /* >> double basepair move */
    if ((i > 0) && (j > 0)) { /* insert */
      rli = &rl[i-GSV.len-2];
      rlj = &rl[j-GSV.len-2];
      close_bp_en(rli->next, rlj->prev);
      close_bp_en(rli, rlj);
    }
    else if ((i < 0)&&(j < 0)) { /* delete */
      i = -i;
      rli = &rl[i-GSV.len-2];
      open_bp_en(rli);
      open_bp_en(rli->next);
    }
  } /* << double basepair move */

}

/* open a particular base pair */
void open_bp(baum *i) {

  baum *in; /* points to i->next */

  /* change string representation */
  GAV.currform[i->nummer] = '.';
  GAV.currform[i->down->nummer] = '.';

  /* change pairtable representation */
  pairList[1 + i->nummer] = 0;
  pairList[1 + i->down->nummer] = 0;

  /* change tree representation */
  in = i->next;
  i->typ = 'u';
  i->down->typ = 'u';
  i->next = i->down->next;
  i->next->prev = i;
  in->prev = i->down;
  i->down->next = in;
  i->down = in->prev->up = NULL;
}

/* close a particular base pair */
void close_bp (baum *i, baum *j) {

  baum *jn; /* points to j->next */

  /* change string representation */
  GAV.currform[i->nummer] = '(';
  GAV.currform[j->nummer] = ')';

  /* change pairtable representation */
  pairList[1 + i->nummer] = 1+ j->nummer;
  pairList[1 + j->nummer] = 1 + i->nummer;

  /* change tree representation */
  jn = j->next;
  i->typ = 'p';
  j->typ = 'q';
  i->down = j;
  j->up = i;
  i->next->prev = j;
  j->next->prev = i;
  j->next = i->next;
  i->next = jn;
}

# if 0
/* for a given tree, generate postorder-list */
static void make_poList (baum *root) {

  baum *stop, *rli;

  if (!root) root = wurzl;
  stop = root->down;

  /* foreach base in ringlist ... */
  for (rli = stop->next; rli != stop; rli = rli->next) {
    /* ... explore subtee if bp found */
    if (rli->typ == 'p') {
      /*  fprintf(stderr, "%d >%d<\n", poListop, rli->nummer); */
      poList[poListop++] = rli;
      if ( poListop > GSV.len+1 ) {
	fprintf(stderr, "Something went wrong in make_poList()\n");
	exit(1);
      }
      make_poList(rli);
    }
  }
  return;
}
#endif

/* for a given ringlist, generate all structures
   with one additional basepair */
static void inb(baum *root) {

  int EoT;
  int E_old, E_new_in, E_new_out;
  baum *stop,*rli,*rlj;

  E_old = root->loop_energy;
  stop=root->down;
  /* loop ringlist over all possible i positions */
  for(rli=stop->next;rli!=stop;rli=rli->next){
    /* potential i-position is already paired */
    if(rli->typ=='p') continue;
    /* loop ringlist over all possible j positions */
    for(rlj=rli->next;rlj!=stop;rlj=rlj->next){
      /* base pair must enclose at least 3 bases */
      if(rlj->nummer - rli->nummer < MYTURN) continue;
      /* potential j-position is already paired */
      if(rlj->typ=='p') continue;
      /* if i-j can form a base pair ... */
      if(ptype[rli->nummer][rlj->nummer]){
	/* close the base bair and ... */
	close_bp(rli,rlj);
#if HAVE_LIBRNA_API3
        E_new_in  = vrna_eval_loop_pt(GAV.vc, rli->nummer+1, pairList);
        E_new_out = vrna_eval_loop_pt(GAV.vc, root->nummer+1, pairList);
#else
	E_new_in  = loop_energy(pairList, typeList, aliasList,rli->nummer+1);
	E_new_out = loop_energy(pairList, typeList, aliasList,root->nummer+1);
#endif
	/* ... evaluate energy of the structure */
	EoT = (int) (GSV.currE*100 + ((GSV.currE<0)?-0.4:0.4)) +  E_new_in + E_new_out - E_old ;
	/* assert(EoT ==  energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params)); */
	/* open the base pair again... */
	open_bp(rli);
	/* ... and put the move and the enegy
	   of the structure into the neighbour list */
	update_nbList(1 + rli->nummer, 1 + rlj->nummer, EoT);
      }
    }
  }
}

/* for a given ringlist, generate all structures (canonical)
   with one additional base pair (BUT WITHOUT ISOLATED BASE PAIRS) */
static void inb_nolp(baum *root) {

  int EoT = 0;
  baum *stop, *rli, *rlj;

  stop = root->down;
    /* loop ringlist over all possible i positions */
  for (rli=stop->next;rli!=stop;rli=rli->next) {
    /* potential i-position is already paired */
    if (rli->typ=='p') continue;
    /* loop ringlist over all possible j positions */
    for (rlj=rli->next;rlj!=stop;rlj=rlj->next) {
      /* base pair must enclose at least 3 bases */
      if (rlj->nummer - rli->nummer < MYTURN) continue;
      /* potential j-position is already paired */
      if (rlj->typ=='p') continue;
      /* if i-j can form a base pair ... */
      if (ptype[rli->nummer][rlj->nummer]) {
	/* ... and extends a helix ... */
	if (((rli->prev==stop && rlj->next==stop) && stop->typ != 'x') ||
	    (rli->next == rlj->prev)) {
	  /* ... close the base bair and ... */
	  close_bp(rli,rlj);
	  /* ... evaluate energy of the structure */
#if HAVE_LIBRNA_API3
	  EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
	  EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
	  /* open the base pair again... */
	  open_bp(rli);
	  /* ... and put the move and the enegy
	     of the structure into the neighbour list */
	  update_nbList(1 + rli->nummer, 1 + rlj->nummer, EoT);
	}
	/* if double insertion is possible ... */
	else if ((rlj->nummer - rli->nummer >= MYTURN+2)&&
		 (rli->next->typ != 'p' && rlj->prev->typ != 'p') &&
		 (rli->next->next != rlj->prev->prev) &&
		 (ptype[rli->next->nummer][rlj->prev->nummer])) {
	  /* close the two base bair and ... */
	  close_bp(rli->next, rlj->prev);
	  close_bp(rli, rlj);
	  /* ... evaluate energy of the structure */
#if HAVE_LIBRNA_API3
	  EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
	  EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
	  /* open the two base pair again ... */
	  open_bp(rli);
	  open_bp(rli->next);
	  /* ... and put the move and the enegy
	     of the structure into the neighbour list */
	  update_nbList(1+rli->nummer+GSV.len+1, 1+rlj->nummer+GSV.len+1, EoT);
	}
      }
    }
  }
}

/* for a given ringlist, generate all structures
 with one less base pair */
static void dnb(baum *rli){

  int EoT, E_old_in, E_old_out, E_new;

  baum *rlj, *r;

  rlj=rli->down;
  open_bp(rli);
  /* ... evaluate energy of the structure */

  for (r=rli->next; r->up==NULL; r=r->next);
  E_old_in = rli->loop_energy;
  E_old_out = r->up->loop_energy;
#if HAVE_LIBRNA_API3
  E_new = vrna_eval_loop_pt(GAV.vc, r->up->nummer+1, pairList);
#else
  E_new = loop_energy(pairList,typeList,aliasList,r->up->nummer+1);
#endif
  EoT = (int) (GSV.currE*100 + ((GSV.currE<0)?-0.4:0.4)) -
    E_old_in - E_old_out + E_new;

  /* assert(EoT== energy_of_struct_pt(GAV.farbe, pairList, typeList, aliasList));*/
  close_bp(rli,rlj);
  update_nbList(-(1 + rli->nummer), -(1 + rlj->nummer), EoT);
}

/* for a given ringlist, generate all structures (canonical)
 with one less base pair (BUT WITHOUT ISOLATED BASE PAIRS) */
static void dnb_nolp(baum *rli) {

  int EoT = 0;
  baum *rlj;
  baum *rlin = NULL; /* pointers to following pair in helix, if any */
  baum *rljn = NULL;
  baum *rlip = NULL; /* pointers to preceding pair in helix, if any */
  baum *rljp = NULL;

  rlj = rli->down;

  /* immediate interior base pair ? */
  if (rlj->next == rlj->prev) {
    rlin = rlj->next;
    rljn = rlin->down;
  }

  /* immediate exterior base pair and not virtualroot ? */
  if (rli->prev == rli->next && rli->next->typ != 'x') {
    rlip = rli->next->up;
    rljp = rli->next;
  }

  /* double delete ? */
  if (rlip==NULL && rlin && rljn->next != rljn->prev ) {
    /* open the two base pairs ... */
    open_bp(rli);
    open_bp(rlin);
    /* ... evaluate energy of the structure ... */
#if HAVE_LIBRNA_API3
    EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
    EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
    /* ... and put the move and the enegy
       of the structure into the neighbour list ... */
    update_nbList(-(1+rli->nummer+GSV.len+1),-(1+rlj->nummer+GSV.len+1), EoT);
    /* ... and close the two base pairs again */
    close_bp(rlin, rljn);
    close_bp(rli, rlj);
  } else { /* single delete */
    /* the following will work only if boolean expr are shortcicuited */
    if (rlip==NULL || (rlip->prev == rlip->next && rlip->prev->typ != 'x'))
      if (rlin ==NULL || (rljn->next == rljn->prev)) {
	/* open the base pair ... */
	open_bp(rli);
	/* ... evaluate energy of the structure ... */
#if HAVE_LIBRNA_API3
	EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
	EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
	/* ... and put the move and the enegy
	   of the structure into the neighbour list ... */
	update_nbList(-(1 + rli->nummer),-(1 + rlj->nummer), EoT);
	/* and close the base pair again */
	close_bp(rli, rlj);
      }
  }
}

/* for a given ringlist, generate all structures
 with one shifted base pair */
static void fnb(baum *rli) {

  int EoT = 0, x;
  baum *rlj, *stop, *help_rli, *help_rlj;

  stop = rli->down;

  /* examin interior loop of bp(ij); (.......)
     i of j move                      ->   <- */
  for (rlj = stop->next; rlj != stop; rlj = rlj->next) {
    /* prevent shifting to paired position */
    if ((rlj->typ=='p')||(rlj->typ=='q')) continue;
    /* j-position of base pair shifts to k position (ij)->(ik) i<k<j */
    if ( (rlj->nummer-rli->nummer >= MYTURN)
	 && (ptype[rli->nummer][rlj->nummer]) ) {
      /* open original basepair */
      open_bp(rli);
      /* close shifted version of original basepair */
      close_bp(rli, rlj);
      /* evaluate energy of the structure */
#if HAVE_LIBRNA_API3
      EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
      EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
      /* put the move and the enegy of the structure into the neighbour list */
      update_nbList(1+rli->nummer, -(1+rlj->nummer), EoT);
      /* open shifted basepair */
      open_bp(rli);
      /* restore original basepair */
      close_bp(rli, stop);
    }
    /* i-position of base pair shifts to position k (ij)->(kj) i<k<j */
    if ( (stop->nummer-rlj->nummer >= MYTURN)
	 && (ptype[stop->nummer][rlj->nummer]) ) {
      /* open original basepair */
      open_bp(rli);
      /* close shifted version of original basepair */
      close_bp(rlj, stop);
      /* evaluate energy of the structure */
#if HAVE_LIBRNA_API3
      EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
      EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
      /* put the move and the enegy of the structure into the neighbour list */
      update_nbList(-(1 + rlj->nummer), 1 + stop->nummer, EoT);
      /* open shifted basepair */
      open_bp(rlj);
      /* restore original basepair */
      close_bp(rli, stop);
    }
  }
  /* examin exterior loop of bp(ij);   (.......)
     i or j moves                    <-         -> */
  for (rlj=rli->next;rlj!=rli;rlj=rlj->next) {
    if ((rlj->typ=='p') || (rlj->typ=='q' ) || (rlj->typ=='x')) continue;
    x=rlj->nummer-rli->nummer;
    if (x<0) x=-x;
    /* j-position of base pair shifts to position k */
    if ((x >= MYTURN) && (ptype[rli->nummer][rlj->nummer])) {
      if (rli->nummer<rlj->nummer) {
	help_rli=rli;
	help_rlj=rlj;
      }
      else {
	help_rli=rlj;
	help_rlj=rli;
      }
      /* open original basepair */
      open_bp(rli);
      /* close shifted version of original basepair */
      close_bp(help_rli,help_rlj);
      /* evaluate energy of the structure */
#if HAVE_LIBRNA_API3
      EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
      EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
      /* put the move and the enegy of the structure into the neighbour list */
      update_nbList(1 + rli->nummer, -(1 + rlj->nummer), EoT);
      /* open shifted base pair */
      open_bp(help_rli);
      /* restore original basepair */
      close_bp(rli,stop);
    }
    x = rlj->nummer-stop->nummer;
    if (x < 0) x = -x;
    /* i-position of base pair shifts to position k */
    if ((x >= MYTURN) && (ptype[stop->nummer][rlj->nummer])) {
      if (stop->nummer < rlj->nummer) {
	help_rli = stop;
	help_rlj = rlj;
      }
      else {
	help_rli = rlj;
	help_rlj = stop;
      }
      /* open original basepair */
      open_bp(rli);
       /* close shifted version of original basepair */
      close_bp(help_rli, help_rlj);
      /* evaluate energy of the structure */
#if HAVE_LIBRNA_API3
      EoT = vrna_eval_structure_pt(GAV.vc, pairList);
#else
      EoT = energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0);
#endif
      /* put the move and the enegy of the structure into the neighbour list */
      update_nbList(-(1 + rlj->nummer), 1 + stop->nummer, EoT);
      /* open shifted basepair */
      open_bp(help_rli);
      /* restore original basepair */
      close_bp(rli,stop);
    }
  }
}

/* for a given tree (structure),
   generate all neighbours according to moveset */
void move_it (void) {
  int i;
  
#if HAVE_LIBRNA_API3
  GSV.currE = (float)vrna_eval_structure_pt(GAV.vc, pairList)/100.;
#else
  GSV.currE =
    energy_of_struct_pt_par(GAV.farbe, pairList, typeList, aliasList, GAV.params, 0)/100.;
#endif
  
  if ( GTV.noLP ) { /* canonical neighbours only */
    inb_nolp(wurzl);
    for (i = 0; i < GSV.len; i++) {
      
      if (pairList[i+1]>i+1) {
	inb_nolp(rl+i);      /* insert pair neighbours */
	dnb_nolp(rl+i);  /* delete pair neighbour */
      }
    }
  }
  else { /* all neighbours */
    inb(wurzl);
    for (i = 0; i < GSV.len; i++) {
      
      if (pairList[i+1]>i+1) {
	inb(rl+i); 	 /* insert pair neighbours */
	dnb(rl+i);  /* delete pair neighbour */
	if ( GTV.noShift == 0 ) fnb(rl+i);
      }
    }
  }
}


/**/
void clean_up_rl(void) {
  int i;
  free(pairList); pairList=NULL;
  free(typeList); typeList = NULL;
  free(aliasList); aliasList = NULL;
  free(rl); rl=NULL;
  free(wurzl);  wurzl=NULL;
  for (i=0; i<=GSV.len; i++)
    free(ptype[i]);
  free(ptype);
  ptype=NULL;
}

/**/
static int comp_struc(const void *A, const void *B) {
  int aE, bE;
  aE = (int)(100 * energy_of_structure(GAV.farbe_full, ((char **)A)[0], 0));
  bE = (int)(100 * energy_of_structure(GAV.farbe_full, ((char **)B)[0], 0));
  return (aE-bE);
}

#if 0
/**/
static void rl_status(void) {

  int i;

  printf("\n%s\n%s\n", GAV.farbe, GAV.currform);
  for (i=0; i <= GSV.len; i++) {
    printf("%2d %c %c %2d %2d %2d %2d\n",
	   rl[i].nummer,
	   i == GSV.len ? 'X': GAV.farbe[i],
	   rl[i].typ,
	   rl[i].up==NULL?0:(rl[i].up)->nummer,
	   rl[i].down==NULL?0:(rl[i].down)->nummer,
	   (rl[i].prev)->nummer,
	   (rl[i].next)->nummer);
  }
  printf("---\n");
}
#endif

#define TURN 3
static void make_ptypes(const short *S) {
  int n,i,j,k,l;
  n=S[0];
  for (k=1; k<n; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+l+TURN;
      if (!SAME_STRAND(i,j)) j=cut_point;
      if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
	if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
	if (noLonelyPairs && (!otype) && (!ntype))
	  type = 0; /* i.j can only form isolated pairs */
	ptype[i-1][j-1] = ptype[j-1][i-1] = (char) type;
	otype =  type;
	type  = ntype;
	i--; j++;
      }
    }
}

static void close_bp_en (baum *i, baum *j) {
  /* close bp and update energy */
  baum *r;
  close_bp(i,j);

#if HAVE_LIBRNA_API3
  i->loop_energy = vrna_eval_loop_pt(GAV.vc, i->nummer+1, pairList);
#else
  i->loop_energy = loop_energy(pairList,typeList,aliasList,i->nummer+1);
#endif

  for (r=i->next; r->up==NULL; r=r->next);

#if HAVE_LIBRNA_API3
  r->up->loop_energy = vrna_eval_loop_pt(GAV.vc, r->up->nummer+1, pairList);
#else
  r->up->loop_energy = loop_energy(pairList,typeList,aliasList,r->up->nummer+1);
#endif
};

static void open_bp_en (baum *i) {
  /* open bp and update energy */
  baum *r;
  i->loop_energy=0;
  open_bp(i);
  for (r=i->next; r->up==NULL; r=r->next);
#if HAVE_LIBRNA_API3
  r->up->loop_energy = vrna_eval_loop_pt(GAV.vc, r->up->nummer+1, pairList);
#else
  r->up->loop_energy = loop_energy(pairList,typeList,aliasList,r->up->nummer+1);
#endif
};
