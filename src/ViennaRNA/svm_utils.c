
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "svm_utils.h"

#include "ViennaRNA/model_avg.inc"  /* defines avg_model_string */
#include "ViennaRNA/model_sd.inc"   /* defines sd_model_string */



PRIVATE struct svm_model *avg_model;
PRIVATE struct svm_model *sd_model;

PRIVATE void    freeFields(char** fields);
PRIVATE char**  splitFields(char* string);
PRIVATE char**  splitLines(char* string);

PUBLIC float get_z(char *sequence, double energy) {
  double average_free_energy;
  double sd_free_energy;
  float my_z = 0.;
  int info_avg;
  make_pair_matrix();
  short *S      = encode_sequence(sequence, 0);
  unsigned int   length  = strlen(sequence);
  int   *AUGC   = get_seq_composition(S, 1, length, length);
  avg_model     = svm_load_model_string(avg_model_string);
  sd_model      = svm_load_model_string(sd_model_string);
  average_free_energy = avg_regression(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4], avg_model, &info_avg);

  if(info_avg == 0){
    double difference = (energy/* /100*/) - average_free_energy;
    sd_free_energy    = sd_regression(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4], sd_model);
    my_z              = difference / sd_free_energy;
  }
  else{
    vrna_message_warning("sequence out of bounds");
#if 0
    my_z = shuffle_score(sequence, energy);
#endif
  }
  free(AUGC);
  free(S);
  svm_free_model_content(avg_model);
  svm_free_model_content(sd_model);
  return my_z;
}

PUBLIC int *get_seq_composition(short *S, unsigned int start, unsigned int stop, unsigned int length){
  unsigned int i;
  int *ret = (int *)vrna_alloc(sizeof(int) * 6);

  for (i=MAX2(start, 1); i <= MIN2(stop, length); i++){
    if(S[i] > 4)  ret[0]++;
    else          ret[S[i]]++;
  }
  ret[5] = -1; /* indicate last entry */
  return ret;
}

PUBLIC double sd_regression(int N, int A, int C, int G, int T,  struct svm_model *sd_model){
  double sd_free_energy = 0.0;
  int length = A + C + G + T + N;
  double GC_content  = (double) (G + C)/length;
  double AT_ratio    = (double) A/(A+T);
  double CG_ratio    = (double) C/(C+G);
  double norm_length = (double) (length-50)/350.0;
  struct svm_node node_mono[5];

  node_mono[0].index = 1; node_mono[0].value = GC_content;
  node_mono[1].index = 2; node_mono[1].value = AT_ratio;
  node_mono[2].index = 3; node_mono[2].value = CG_ratio;
  node_mono[3].index = 4; node_mono[3].value = norm_length;
  node_mono[4].index =-1;

  sd_free_energy = svm_predict(sd_model,node_mono);

  sd_free_energy = (double) sd_free_energy * sqrt(length);

  return sd_free_energy;
}

PUBLIC double avg_regression(int N, int A, int C, int G, int T, struct svm_model *avg_model, int *info ){
  double average_free_energy = 0.0;

  int length = A + C + G + T + N;
  double N_fraction = (double) N/length;
  double GC_content = (double) (G + C)/length;
  double AT_ratio   = (double) A/(A+T);
  double CG_ratio   = (double) C/(C+G);

  double norm_length = (double) (length-50)/350.0;

  struct svm_node node_mono[5];
  *info = 0;
  if ( length < 50 || length > 400 ) {
    *info = 1;
    return 0.0;
  }
  if ( N_fraction > 0.05 ) {
    *info = 2;
    return 0.0;
  }
  if ( GC_content < 0.20 || GC_content > 0.80 ) {
    *info = 3;
    return 0.0;
  }
  if ( AT_ratio < 0.20 || AT_ratio > 0.80 ) {
    *info = 4;
    return 0.0;
  }
  if ( CG_ratio < 0.20 || CG_ratio > 0.80 ) {
    *info = 5;
    return 0.0;
  }

  node_mono[0].index = 1; node_mono[0].value = GC_content;
  node_mono[1].index = 2; node_mono[1].value = AT_ratio;
  node_mono[2].index = 3; node_mono[2].value = CG_ratio;
  node_mono[3].index = 4; node_mono[3].value = norm_length;
  node_mono[4].index =-1;

  average_free_energy = svm_predict(avg_model,node_mono);

  average_free_energy = (double) average_free_energy * length;

  return average_free_energy;
}

PUBLIC double minimal_sd(int N, int A, int C, int G, int T ){
  int length = A + C + G + T + N;
  if ( length <  60 ) return 0.450324;
  if ( length <  70 ) return 0.749771;
  if ( length <  80 ) return 1.029421;
  if ( length <  90 ) return 1.027517;
  if ( length <  100 ) return 1.347283;
  if ( length <  120 ) return 1.112086;
  if ( length <  150 ) return 1.574339;
  if ( length <  170 ) return 1.779043;
  if ( length <  200 ) return 1.922908;
  if ( length <  250 ) return 2.226856;
  if ( length <  300 ) return 2.349300;
  if ( length <  350 ) return 2.589703;
  if ( length <  400 ) return 2.791215;

  return 0.450324;
}

PUBLIC struct svm_model  *svm_load_model_string(char *modelString){

  /* redefinition from svm.cpp */
  char *svm_type_table[]={"c_svc","nu_svc","one_class","epsilon_svr","nu_svr",NULL};
  char *kernel_type_table[]={"linear","polynomial","rbf","sigmoid",NULL};

  struct svm_model *model;
  char **lines, **fields;
  int i,j,k,l,m;
  char *key, *value, *field;
  char c;
  int dataStart, elements;
  int isColon;
  struct svm_node *x_space=NULL;

  model = (struct svm_model*)vrna_alloc(sizeof(struct svm_model));

  model->rho = NULL;
  model->probA = NULL;
  model->probB = NULL;
  model->label = NULL;
  model->nSV = NULL;


  /* Read header until support vectors start */
  lines=splitLines(modelString);
  i=0;
  while (lines[i] && (strcmp(lines[i],"SV")!=0)){
        fields=splitFields(lines[i]);

        key=fields[0];

        if(strcmp(key,"svm_type")==0){
          value=fields[1];
          for(j=0;svm_type_table[j];j++){
                if(strcmp(svm_type_table[j],value)==0){
                  model->param.svm_type=j;
                  break;
                }
          }
          if(svm_type_table[i] == NULL){
                vrna_message_warning("unknown svm type.");
                free(model->rho);
                free(model->label);
                free(model->nSV);
                free(model);
                return NULL;
          }
        } else

        if(strcmp(key,"kernel_type")==0){
          value=fields[1];
          for(j=0;kernel_type_table[j];j++){
                if(strcmp(kernel_type_table[j],value)==0){
                  model->param.kernel_type=j;
                  break;
                }
          }
          if(kernel_type_table[i] == NULL){
                vrna_message_warning("unknown kernel type.");
                free(model->rho);
                free(model->label);
                free(model->nSV);
                free(model);
                return NULL;
          }
        } else

        if (strcmp(key,"gamma")==0){
          value=fields[1];
          sscanf(value,"%lf",&model->param.gamma);
        }

        if (strcmp(key,"degree")==0){
          value=fields[1];
          sscanf(value,"%d",&model->param.degree);
        } else

        if (strcmp(key,"coef0")==0){
          value=fields[1];
          sscanf(value,"%lf",&model->param.coef0);
        } else
        if (strcmp(key,"nr_class")==0){
          value=fields[1];
          sscanf(value,"%d",&model->nr_class);
        } else
        if (strcmp(key,"total_sv")==0){
          value=fields[1];
          sscanf(value,"%d",&model->l);
        } else

        if (strcmp(key,"rho")==0){
          int n = model->nr_class * (model->nr_class-1)/2;
          model->rho = (double*)vrna_alloc(sizeof(double)*n);
          for(j=0;j<n;j++){
                sscanf(fields[j+1],"%lf",&model->rho[j]);
          }
        } else

        if (strcmp(key,"nr_sv")==0){
          int n = model->nr_class;
          model->nSV = (int*)vrna_alloc(sizeof(int)*n);
          for(j=0;j<n;j++){
                sscanf(fields[j+1],"%d",&model->nSV[j]);
          }
        } else

        if (strcmp(key,"label")==0){
          int n = model->nr_class;
          model->label = (int*)vrna_alloc(sizeof(int)*n);
          for(j=0;j<n;j++){
                sscanf(fields[j+1],"%d",&model->label[j]);
          }
        } else

        if (strcmp(key,"probA")==0){
          int n = model->nr_class * (model->nr_class-1)/2;
          model->probA = (double*)vrna_alloc(sizeof(double)*n);
          for(j=0;j<n;j++){
                sscanf(fields[j+1],"%lf",&model->probA[j]);
          }
        } else

        if (strcmp(key,"probB")==0){
          int n = model->nr_class * (model->nr_class-1)/2;
          model->probB = (double*)vrna_alloc(sizeof(double)*n);
          for(j=0;j<n;j++){
                sscanf(fields[j+1],"%lf",&model->probB[j]);
          }
        }
        i++;
        freeFields(fields);
  }

  dataStart=i+1;
  elements=0;

  /* Count number of nodes (by counting colons) in advance to allocate
         memory in one block */
  while (lines[i]!=NULL){
        j=0;
        while ((c=lines[i][j])!='\0'){
          if (c==':'){
                elements++;
          }
          j++;
        }
        elements++;
        i++;
  }

  /* allocate memory for SVs and coefficients */
  m = model->nr_class - 1;
  l = model->l;
  model->sv_coef = (double**)vrna_alloc(sizeof(double*)*m);
  for(i=0;i<m;i++){
        model->sv_coef[i] = (double*)vrna_alloc(sizeof(double)*l);
  }
  model->SV = (struct svm_node**)vrna_alloc(sizeof(struct svm_node*)*l);

  if(l>0){
    x_space = (struct svm_node*)vrna_alloc(sizeof(struct svm_node)*(elements));
  }


  /* parse support vector data */
  j=0;
  for(i=0;i<l;i++){
        fields=splitFields(lines[dataStart+i]);
        model->SV[i] = &x_space[j];
        k=0;
        while ((field=fields[k])!=NULL){
          if (k<m){
            sscanf(fields[k],"%lf",&model->sv_coef[k][i]);
          } else {
            sscanf(fields[k],"%d:%lf",&(x_space[j].index),&(x_space[j].value));
            j++;
          }
          k++;
        }
        x_space[j++].index = -1;
        freeFields(fields);
  }
  freeFields(lines);

  model->free_sv = 1;

  return(model);
}

PRIVATE char **splitFields(char* string){

  char c;
  char* currField;
  char** output=NULL;
  int* seps;
  int nSep;
  int nField=0;
  int i=0;

  if (strlen(string)==0 || string==NULL){
        return NULL;
  }

  /* First find all characters which are whitespaces and store the
         positions in the array seps */

  seps=(int *)vrna_alloc(sizeof(int));
  seps[0]=-1;
  nSep=1;

  while ((c=string[i])!='\0' && (c!='\n')){
        if (isspace(c)){
          seps=(int*)vrna_realloc(seps,sizeof(int)*(nSep+1));
          seps[nSep++]=i;
        }
        i++;
  }

  seps=(int*)vrna_realloc(seps,sizeof(int)*(nSep+1));
  seps[nSep]=strlen(string);


  /* Then go through all intervals in between of two whitespaces (or
         end or start of string) and store the fields in the array
         "output"; if there are two adjacent whitespaces this is ignored
         resulting in a behaviour like "split /\s+/" in perl */

  for (i=0;i<nSep;i++){

        int start=seps[i];
        int stop=seps[i+1];
        int length=(stop-start);
        int notSpace,j;


        currField=(char *)vrna_alloc(sizeof(char)*(length+1));
        strncpy(currField,string+start+1,length-1);
        currField[length]='\0';

        /* check if field is not only whitespace */
        notSpace=0;
        j=0;
        while ((c=currField[j])!='\0'){
          if (!isspace(c)){
                notSpace=1;
                break;
          }
        }

        if (notSpace){
          output=(char**)vrna_realloc(output,sizeof(char**)*(nField+1));
          output[nField++]=currField;
          currField=NULL;
        } else {
          free(currField);
          currField=NULL;
        }

        /* printf("%s|\n",output[nField-1]); */
  }

  if (nField==0){
        return NULL;
  }


  output=(char**)vrna_realloc(output,sizeof(char**)*(nField+1));
  output[nField]=NULL;

  free(seps);
  return output;

}

PRIVATE char **splitLines(char* string){

  char c;
  char* currLine=NULL;
  char** output=NULL;
  int i=0;
  int currLength=0;
  int lineN=0;

  while ((c=string[i])!='\0'){

        if (c=='\n'){
          output=(char**)vrna_realloc(output,sizeof(char**)*(lineN+1));
          currLine=(char*)vrna_realloc(currLine,sizeof(char)*(currLength+1));
          currLine[currLength]='\0';
          output[lineN]=currLine;
          currLength=0;
          currLine=NULL;
          lineN++;
        } else {

          currLine=(char*)vrna_realloc(currLine,sizeof(char)*(currLength+1));
          currLine[currLength]=c;
          currLength++;
        }
        i++;
  }

  output=(char**)vrna_realloc(output,sizeof(char**)*(lineN+1));
  output[lineN]=NULL;

  return output;

}

/*  for both splitLines and splitFields */
void freeFields(char** fields){

  int i=0;
  while (fields[i]!=NULL){
        free(fields[i++]);
  }
  free(fields);
}
