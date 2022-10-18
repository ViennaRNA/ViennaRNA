#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/utils/basic.h>

#include "meshpoint.h"

int
insert_meshpoint_hidden(char            *s,
                        float           en,
                        float           structureEnergy,
                        meshpoint_list  *list,
                        int             maxEntries,
                        int             with_struct_en)
{
  if (list->count >= maxEntries) {
    /* do we insert the new meshpoint or do we skip it? */
    if (en > list->worstEnergy) {
      return 0;
    } else if (en == list->worstEnergy) {
      if (with_struct_en) {
        if (structureEnergy >= list->worstStructureEnergy) {
          return 0;
        }
        /* however in case we also take structural energy into account and ('en == list->worstEnergy' && 'structureEnergy < list->worstStructureEnergy') we insert this meshpoint */
        else {
          meshpoint *cur  = list->first;
          meshpoint *prev = NULL;
          for (cur = list->first; cur != NULL; prev = cur, cur = cur->next) {
            if (cur->en != en) {
              continue;
            } else if (cur->struct_en <= structureEnergy) {
              continue;
            }
            /* now 'cur' should point to the meshpoint in which's front we want to insert the new meshpoint */
            /* so we allocate memory for our new meshpoint and connect it into the list */
            else {
              meshpoint *tmp    = cur;
              meshpoint *newMP  = (meshpoint *)vrna_alloc(sizeof(meshpoint));
              if (prev == NULL) {
                list->first = newMP;
                cur         = list->first;
              } else {
                prev->next  = newMP;
                cur         = prev->next;
              }

              cur->s          = strdup(s);
              cur->en         = en;
              cur->struct_en  = structureEnergy;
              cur->next       = tmp;
              break;
            }
          }
          /* now after we've inserted the new meshpoint we have to delete the last entry in the list to not exceed the maximum entries */
          while (cur->next != NULL) {
            prev  = cur;
            cur   = cur->next;
          }
          free(cur->s);
          free(cur);
          prev->next                  = NULL;
          list->worstEnergy           = prev->en;
          list->worstStructureEnergy  = prev->struct_en;
        }
      } else {
        return 0;
      }
    } else {
      /* we search for a good place to insert and throw the last meshpoint in list into the void */
      meshpoint *cur  = list->first;
      meshpoint *prev = NULL;
      for (cur = list->first; cur != NULL; prev = cur, cur = cur->next) {
        /* if the current energy is better to the one we want to insert,
         *  we have to go further...
         */
        if (cur->en < en) {
          continue;
        } else if (cur->en == en && !with_struct_en) {
          continue;
        } else if (cur->en == en && (cur->struct_en <= structureEnergy)) {
          continue;
        } else {
          /* else we just insert our new meshpoint right before the next worse meshpoint */
          meshpoint *tmp    = cur;
          meshpoint *newMP  = (meshpoint *)vrna_alloc(sizeof(meshpoint));
          /* if we are at the beginning of our list, we have to do smth. different... */
          if (prev == NULL) {
            list->first = newMP;
            cur         = list->first;
          } else {
            prev->next  = newMP;
            cur         = prev->next;
          }

          cur->s  = strdup(s);
          cur->en = en;
          if (with_struct_en)
            cur->struct_en = structureEnergy;

          cur->next = tmp;
          break;
        }
      }
      /* after we have inserted a new meshpoint into the full list, we have to delete the last element in the list,
       *  as it is one too much. We also have to update the lists worst energy, while the list counter remains as it
       *  is
       */
      while (cur->next != NULL) {
        prev  = cur;
        cur   = cur->next;
      }
      /* now we've reached the end of the list...
       *  if we want to delete the last element in the list, we do not forget to free its pair_table before freeing
       *  itself
       */
      free(cur->s);
      free(cur);
      prev->next        = NULL;
      list->worstEnergy = prev->en;
      if (with_struct_en)
        list->worstStructureEnergy = prev->struct_en;
    }
  } else {
    /* so we have enough space to insert the new meshpoint in our list...
     *  we search for a good place to insert
     */
    meshpoint *cur  = NULL;
    meshpoint *prev = NULL;
    for (cur = list->first; cur != NULL; prev = cur, cur = cur->next) {
      /* if the current energy is better or equal to the one we want to insert,
       *  we have to go further...
       */
      if (cur->en < en) {
        continue;
      } else if (cur->en == en && !with_struct_en) {
        continue;
      } else if (cur->en == en && (cur->struct_en <= structureEnergy)) {
        continue;
      } else {
        /* else we just insert our new meshpoint right before the next worse meshpoint */
        meshpoint *tmp    = cur;
        meshpoint *newMP  = (meshpoint *)vrna_alloc(sizeof(meshpoint));
        /* if we are at the beginning of our list, we have to do smth. different... */
        if (prev == NULL) {
          list->first = newMP;
          cur         = list->first;
        } else {
          prev->next  = newMP;
          cur         = prev->next;
        }

        cur->s  = strdup(s);
        cur->en = en;
        if (with_struct_en)
          cur->struct_en = structureEnergy;

        cur->next = tmp;
        break;
      }
    }
    if (cur == NULL) {
      if (prev == NULL) {
        /* this is the very first element in our list, so we have to do the memory
         *  allocating stuff now, as it was not done in the loop above
         */
        list->first = (meshpoint *)vrna_alloc(sizeof(meshpoint));
        cur         = list->first;
        cur->s      = strdup(s);
        cur->en     = en;
        if (with_struct_en)
          cur->struct_en = structureEnergy;

        cur->next = NULL;
      } else {
        /* we insert our new meshpoint at the end of the list */
        prev->next  = (meshpoint *)vrna_alloc(sizeof(meshpoint));
        cur         = prev->next;
        cur->s      = strdup(s);
        cur->en     = en;
        if (with_struct_en)
          cur->struct_en = structureEnergy;

        cur->next = NULL;
      }

      list->worstEnergy = en;
      if (with_struct_en)
        list->worstStructureEnergy = structureEnergy;
    }

    /* and we increment the list counter.... */
    list->count++;
  }

  return 1;
}


int
insert_meshpoint(char           *s,
                 float          en,
                 meshpoint_list *list,
                 int            maxEntries)
{
  return insert_meshpoint_hidden(s, en, 0., list, maxEntries, 0);
}


int
insert_meshpoint_with_struct_energy(char            *s,
                                    float           en,
                                    float           structEn,
                                    meshpoint_list  *list,
                                    int             maxEntries)
{
  return insert_meshpoint_hidden(s, en, structEn, list, maxEntries, 0);
}


void
clear_meshpoints(meshpoint_list *list)
{
  meshpoint *cur  = NULL;
  meshpoint *prev = NULL;

  for (cur = list->first; cur != NULL; cur = cur->next) {
    free(cur->s);
    if (prev != NULL)
      free(prev);

    prev = cur;
  }
  if (prev != NULL)
    free(prev);

  list->count                 = 0;
  list->worstEnergy           = 10000.;
  list->worstStructureEnergy  = 10000.;
}


void
init_meshpoint_list(meshpoint_list *list)
{
  list->count                 = 0;
  list->worstEnergy           = 10000.;
  list->worstStructureEnergy  = 10000.;
  list->first                 = NULL;
}


void
init_structure_queue(structure_queue *queue)
{
  queue->count  = 0;
  queue->first  = NULL;
  queue->last   = NULL;
}


void
clear_structure_queue(structure_queue *queue)
{
  meshpoint *cur = NULL;

  for (cur = queue->first; cur != NULL; cur = queue->first) {
    queue->first = cur->next;
    if (is_in_queue(cur->s, queue) == NULL)
      free(cur->neighbor_list);

    free(cur->s);
    free(cur);
  }
  queue->count  = 0;
  queue->first  = NULL;
  queue->last   = NULL;
}


void
insert_structure_in_queue(structure_queue *queue,
                          char            *s,
                          float           en,
                          neighbor        *neighbor_list,
                          int             neighbor_cnt,
                          float           bestNeighborEn,
                          int             bestNeighborIdx,
                          int             maxStructs)
{
  //fprintf(stdout, "inserting %s\n", s);
  meshpoint *newMP = (meshpoint *)vrna_alloc(sizeof(meshpoint));

  newMP->s                = strdup(s);
  newMP->en               = en;
  newMP->neighbor_cnt     = neighbor_cnt;
  newMP->bestNeighborEn   = bestNeighborEn;
  newMP->bestNeighborIdx  = bestNeighborIdx;
  newMP->neighbor_list    = neighbor_list;
  newMP->next             = NULL;
  if (queue->count < maxStructs) {
    if (queue->first != NULL) {
      (queue->last)->next = newMP;
      queue->last         = newMP;
    } else {
      queue->first  = newMP;
      queue->last   = newMP;
    }

    queue->count++;
  } else {
    /* now we have to remove the first elelemt from the queue and insert
     *  the new one
     */
    meshpoint *tmp = queue->first;
    queue->first = tmp->next;

    /* insert new state */
    (queue->last)->next = newMP;
    queue->last         = newMP;
    /* completely delete the removed element */
    /* only remove the neighbor list, if there is no other
     *  state in the queue using this neighbor list
     */
    if (is_in_queue(tmp->s, queue) == NULL)
      free(tmp->neighbor_list);

    free(tmp->s);
    free(tmp);
  }
}


meshpoint *
is_in_queue(char            *s,
            structure_queue *queue)
{
  meshpoint *cur = NULL;

  if (queue == NULL)
    return NULL;

  for (cur = queue->first; cur != NULL; cur = cur->next)
    if (!strcmp(s, cur->s))
      return cur;

  return NULL;
}


int
sort_neighbors_by_energy_asc(const void *p1,
                             const void *p2)
{
  if (((neighbor *)p1)->en > ((neighbor *)p2)->en)
    return 1;
  else if (((neighbor *)p1)->en < ((neighbor *)p2)->en)
    return -1;

  return 0;
}
