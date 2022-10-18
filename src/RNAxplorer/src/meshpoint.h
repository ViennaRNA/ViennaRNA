#ifndef __MESHPOINT_LIST__
#define __MESHPOINT_LIST__

typedef struct neighbor {
  int   i;
  int   j;
  float en;
} neighbor;

typedef struct meshpoint {
  char              *s;
  float             en;
  float             struct_en;
  int               neighbor_cnt;
  float             bestNeighborEn;
  int               bestNeighborIdx;
  struct neighbor   *neighbor_list;
  struct meshpoint  *next;
} meshpoint;

typedef struct meshpoint_list {
  int       count;
  float     worstEnergy;
  float     worstStructureEnergy;
  meshpoint *first;
} meshpoint_list;

typedef struct structure_queue {
  int       count;
  meshpoint *first;
  meshpoint *last;
} structure_queue;

int insert_meshpoint(char           *s,
                     float          en,
                     meshpoint_list *list,
                     int            maxEntries);


int insert_meshpoint_with_struct_energy(char            *s,
                                        float           en,
                                        float           struct_en,
                                        meshpoint_list  *list,
                                        int             maxEntries);


void clear_meshpoints(meshpoint_list *list);


void init_meshpoint_list(meshpoint_list *list);


void init_structure_queue(structure_queue *queue);


void clear_structure_queue(structure_queue *queue);


void insert_structure_in_queue(structure_queue  *queue,
                               char             *s,
                               float            en,
                               neighbor         *neighbor_list,
                               int              neighbor_cnt,
                               float            bestNeighborEn,
                               int              bestNeighborIdx,
                               int              maxStructs);


meshpoint *is_in_queue(char             *s,
                       structure_queue  *queue);


int sort_neighbors_by_energy_asc(const void *p1,
                                 const void *p2);


#endif
