#ifdef _WIN32

int
sanitize_input(const char *string)
{
  return 0; /* always successful if MS Windows */
}

#else

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <termios.h>
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/layouts.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/special_const.h"
#include "ViennaRNA/io/sanitize.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

typedef struct {
  float x;
  float y;
  int   ml;
  float vel[2];
} sectf;

typedef struct {
  int   i;
  int   j;
  int   ml;
} secti;

typedef struct {
  char  *sequence;
  char  *structure;
  float en;
  int   length;
  secti  *struct_list;
  secti  *reactive;
  int   reactive_n;
  secti  *probe;
  int   probe_n;
  sectf *fragments;
  int   fragments_n;
} bt_data;



PRIVATE INLINE int
n2c(char n)
{
  switch (n) {
    case 'A':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 3;
    case 'U':
      return 4;
  }
  return 0;
}



PRIVATE INLINE char
c2n(int c)
{
  switch (c) {
    case 1:
      return 'A';
    case 2:
      return 'C';
    case 3:
      return 'G';
    case 4:
      return 'U';
  }
  return 'N';
}


PRIVATE INLINE void
printCharColXY(const char c,
               int        cc,
               int        x,
               int        y)
{
  printf("\033[%d;%df\033[1;%d;40m%c\033[0m", y, x, cc, c);
}


PRIVATE INLINE void
removeFromGrid(secti *list)
{

  for (; (*list).ml >= 0; list++)
    if (((*list).j < 21) && ((*list).j >= 0))
      printf("\033[%d;%df\033[1;37;40m%c\033[0m", 23 - (*list).j, (*list).i + 1, bg);

  fflush(stdout);
}


PRIVATE INLINE void
putOnGridChar(secti        *list,
              const char  c)
{
  for (; (*list).ml >= 0; list++)
    if (((*list).j < 21) && ((*list).j >= 0))
      printCharColXY(c, (*list).ml + 30, (*list).i + 1, 23 - (*list).j);

  fflush(stdout);
}


PRIVATE INLINE void
putOnGridSeq(secti *list,
             char *seq)
{
  for (; (*list).ml >= 0; list++, seq++)
    if (((*list).j < 21) && ((*list).j >= 0))
      printCharColXY(*seq, (*list).ml + 30, (*list).i + 1, 23 - (*list).j);

  fflush(stdout);
}


PRIVATE INLINE void
getBoundary(secti  *list,
            int   position[])
{
  position[1] = 0;
  position[0] = position[2] = INF;
  for (; (*list).ml >= 0; list++) {
    position[0] = MIN2(position[0], (*list).i);
    position[1] = MAX2(position[1], (*list).i);
    position[2] = MIN2(position[2], (*list).j);
  }
}



PRIVATE INLINE void
placeInjection(unsigned int hpos)
{
  const char *c = &injector[0];

  for (; (*c); c++)
    printf("\033[23;%df\033[1;37;40m%c\033[00m", (int)(hpos + (c - &injector[0]) + 1), *c);
}


PRIVATE INLINE void
clearOutput(void)
{
  int i;

  printf("\033[0;0f\033[1;31;40m\033[J%c", wall);
  for (i = 1; i < 1896; i++)
    !(i % 79) ? printf("%c\n%c", wall, wall) : printf("%c", bg);
  printf("\033[24;0f\033[1;31;40m");
  for (i = 0; i < 80; i++)
    printf("%c", wall);
  printf("\033[00m\033[0;0f");
  fflush(stdout);
}


PRIVATE INLINE sectf *
initFragments(void)
{
  sectf *list = (sectf *)vrna_alloc(sizeof(sectf));

  list[0].ml = -1;
  return list;
}


PRIVATE INLINE secti *
initList(void)
{
  secti *list = (secti *)vrna_alloc(sizeof(secti));

  list[0].ml = -1;
  return list;
}


PRIVATE INLINE void
addFragments(bt_data  *d,
             int      i,
             int      j)
{
  int   cnt;
  float r     = vrna_urn();
  int   *size = &(d->fragments_n);
  sectf *list = d->fragments;

  for (cnt = 0; cnt < 1 + (int)(4. * r); cnt++) {
    float x = vrna_urn();
    float y = vrna_urn();
    list                      = (sectf *)vrna_realloc(list, sizeof(sectf) * (++(*size) + 1));
    list[(*size) - 1].ml      = *size;
    list[(*size) - 1].x       = i;
    list[(*size) - 1].y       = j;
    list[(*size) - 1].vel[0]  = 3. - (6 * x);
    list[(*size) - 1].vel[1]  = 1. + (1.5 * y);
    list[(*size)].ml          = -1;
  }
  d->fragments = list;
}


PRIVATE INLINE void
updateFragments(bt_data *d)
{
  float gravity = 0.3;
  sectf *bla    = d->fragments;
  int   *size   = &(d->fragments_n);
  sectf *list   = d->fragments;

  while ((*list).ml >= 0) {
    printf("\033[%d;%df\033[1;37;40m%c\033[0m", 23 - (int)((*list).y), (int)((*list).x + 1), bg);

    (*list).x += (*list).vel[0];
    (*list).y += (*list).vel[1];

    (*list).vel[1] -= gravity;

    if (((*list).y < 0.) || ((*list).x < 1.) || ((*list).x > 78.)) {
      memmove(list, list + 1, ((*size) - (list - bla)) * sizeof(sectf));
      bla = (sectf *)vrna_realloc(bla, sizeof(sectf) * (*size));
      (*size)--;
      list = bla;
      continue;
    }

    printf("\033[%d;%df\033[1;33;40m'\033[0m", 23 - (int)((*list).y), (int)((*list).x + 1));
    fflush(stdout);
    list++;
  }
  d->fragments = bla;
}


PRIVATE INLINE int
insideBB(int  i,
         int  j,
         int  bb[])
{
  return (j >= bb[2]) && (i >= bb[0]) && (i <= bb[1]);
}



PRIVATE INLINE bt_data *
initStacks(bt_data  *d,
           int      stage)
{
  short *pl     = NULL;
  int   length  = 40, maxsize = 30, n;
  float xmax, xmin, ymax, ymin, *x, *y, xyscaling;
  secti  *list = NULL;

  if (d) {
    free(d->reactive);
    free(d->probe);
    free(d->sequence);
    free(d->structure);
    free(d->struct_list);
    free(d->fragments);
  } else {
    d = (bt_data *)vrna_alloc(sizeof(bt_data));
  }

  d->length       = length;
  d->reactive     = initList();
  d->reactive_n   = 0;
  d->probe        = initList();
  d->probe_n      = 0;
  d->fragments    = initFragments();
  d->fragments_n  = 0;
  d->sequence     = vrna_random_string(length, "AUCG");
  d->structure    = (char *)vrna_alloc((length + 1) * sizeof(char));
  d->en           = vrna_fold(d->sequence, d->structure);
  pl              = vrna_ptable(d->structure);

  x     = y = NULL;
  list  = (secti *)vrna_alloc((length + 1) * sizeof(secti));
  xmax  = ymax = 0;
  xmin  = ymin = INF;
  if (length != vrna_plot_coords_pt(pl, &x, &y, VRNA_PLOT_TYPE_DEFAULT)) {
    vrna_log_error("make_coord_list: error in number of coordinates");
    exit(1);
  }

  for (n = 0; n < length; n++) {
    xmin  = MIN2(xmin, x[n]);
    xmax  = MAX2(xmax, x[n]);
    ymin  = MIN2(ymin, y[n]);
    ymax  = MAX2(ymax, y[n]);
  }
  xyscaling = (float)(MAX2(xmax - xmin, ymax - ymin)) / (float)maxsize;
  for (n = 0; n < length; n++) {
    list[n].i   = (int)((x[n] - xmin) / xyscaling) + 5;
    list[n].j   = (int)((y[n] - ymin) / xyscaling) + (15 - stage);
    list[n].ml  = n2c(d->sequence[n]);
  }
  free(x);
  free(y);
  list[length].ml = -1;
  free(pl);
  d->struct_list = list;

  return d;
}


PRIVATE INLINE void
printHeader(bt_data *d,
            int     stage,
            float   score,
            int     attempts,
            int     probe_max)
{
  secti *t;


  printf("\033[0;0f\033[0;37;40m\033[K\033[0;0f");
  printf(head1l, (stage), (attempts), (d->probe_n), (probe_max));

  if (d->length < 40)
    printf("\033[%dC", 40 - d->length);

  t = (d->struct_list);
  while ((*t).ml >= 0) {
    printf("\033[1;%d;40m%c", (*t).ml + 30, c2n((*t).ml));
    t++;
  }
  printf(" \033[1;31;40m##\033[0;37;40m");

  printf("\033[2;0f\033[0;37;40m\033[K\033[2;0f");
  printf(head2l, (score), (d->en));

  if (d->length < 40)
    printf("\033[%dC", 40 - d->length);

  printf("\033[0;37;40m%s \033[1;31;40m##\033[0;37;40m", (d->structure));

  fflush(stdout);
}


PRIVATE INLINE int
reactiveState(float r)
{
  return (r < 0.1) ? 5 : (r < 0.25) ? 6 : 7;
}


int
sanitize_input(const char *string)
{
  char            state, last_state;
  int             ticks, init_ticks, hpos,
                  cnt2, direction,
                  stage, attempts, probe_max;
  int             inv_bb[3] = { 0 };
  float           score;
  struct termios  neu, alt;
  secti            *ts1, *ts2, *tmp, *t;
  bt_data         *d = NULL;

  if ((strcmp(string, inv)) || (!isatty(fileno(stdin))) || (!isatty(fileno(stdout))))
    return 0;

  hpos      = 50;
  cnt2      = 0;
  inv_bb[2] = 10;
  direction = 1;
  stage     = 1;
  attempts  = 3;
  probe_max = 3;
  score     = 0;

  init_ticks  = ticks = 200;
  last_state  = state = 0;

  tcgetattr(0, &alt);
  tcgetattr(0, &neu);
  neu.c_lflag     &= ~(ICANON | ECHO);
  neu.c_cc[VMIN]  = 0;
  neu.c_cc[VTIME] = 0;
  tcsetattr(0, TCSANOW, &neu);
  printf("%s", start);
  vrna_init_rand();
  while (1) {
    usleep(1000);
    switch (state) {

      case 0:
        if (last_state != 0) {
          last_state = 0;
          printf("%s", start);
        }

        break;


      case 1:
        break;


      case 2:
        if (!cnt2) {
          updateFragments(d);
          placeInjection(hpos);
          removeFromGrid(d->struct_list);
          removeFromGrid(d->probe);
          t = (d->probe);
          while ((*d->probe).ml >= 0) {
            if ((*(d->probe)).j == (20)) {
              memmove(((d->probe)), ((d->probe)) + 1,
                      (((d->probe_n)) - (((d->probe)) - (t))) * sizeof(secti));
              (t)           = (secti *)vrna_realloc((t), sizeof(secti) * ((d->probe_n)));
              ((d->probe))  = (t);
              ((d->probe_n))--;
              continue;
            }

            (*(d->probe)).j += (1);
            (d->probe)++;
          }
          (d->probe) = t;
          putOnGridChar(d->probe, probe);
          removeFromGrid(d->reactive);

          t = (d->reactive);
          while ((*d->reactive).ml >= 0) {
            if ((*(d->reactive)).j == (0)) {
              memmove(((d->reactive)), ((d->reactive)) + 1,
                      (((d->reactive_n)) - (((d->reactive)) - (t))) * sizeof(secti));
              (t)             = (secti *)vrna_realloc((t), sizeof(secti) * ((d->reactive_n)));
              ((d->reactive)) = (t);
              ((d->reactive_n))--;
              continue;
            }

            (*(d->reactive)).j += (-1);
            (d->reactive)++;
          }
          (d->reactive) = t;

          putOnGridChar(d->reactive, star);
          printHeader(d, stage, score, attempts, probe_max);

          int down = 0;
          if (inv_bb[1] >= 78) {
            direction = -1;
            down--;
          } else if (inv_bb[0] <= 1) {
            direction = 1;
            down--;
          }


          tmp = (d->struct_list);
          for (; (*tmp).ml >= 0; tmp++) {
            (*tmp).i  += (direction);
            (*tmp).j  += (down);
          }
          putOnGridSeq(d->struct_list, d->sequence);
          getBoundary(d->struct_list, inv_bb);

          char  *ts = d->sequence;
          char  hit;
          ts1 = d->probe;
          while ((*(d->probe)).ml >= 0) {
            int i, j, p, q;
            hit = 0;
            i   = (*(d->probe)).i;
            j   = (*(d->probe)).j;
            if (insideBB(i, j, inv_bb)) {
              ts2 = d->struct_list;
              ts  = d->sequence;
              while ((*(d->struct_list)).ml >= 0) {
                p = (*(d->struct_list)).i;
                q = (*(d->struct_list)).j;
                if ((i == p) && (j == q)) {
                  float r = vrna_urn();

                  (d->reactive) = (secti *)vrna_realloc((d->reactive),
                                                       sizeof(secti) *
                                                       ((d->reactive_n) + 2));
                  (d->reactive)[(d->reactive_n)].i    = i;
                  (d->reactive)[(d->reactive_n)].j    = j;
                  (d->reactive)[(d->reactive_n)].ml   = reactiveState(r);
                  (d->reactive)[++(d->reactive_n)].ml = -1;

                  printf("\033[%d;%df\033[1;%dm%c\033[00m", 23 - j, i + 1,
                         reactiveState(r) + 30, star);

                  memmove((d->probe), (d->probe) + 1,
                          ((d->probe_n) - ((d->probe) - (ts1))) * sizeof(secti));
                  (ts1)       = (secti *)vrna_realloc((ts1), sizeof(secti) * (d->probe_n));
                  (d->probe)  = (ts1);
                  (d->probe_n)--;

                  memmove(ts, ts + 1, ((d->length) - ((ts) - (d->sequence))) * sizeof(char));
                  d->sequence = (char *)vrna_realloc((d->sequence), sizeof(char) * (d->length));
                  ts          = d->sequence;

                  memmove((d->struct_list), (d->struct_list) + 1,
                          ((d->length) - ((d->struct_list) - (ts2))) * sizeof(secti));
                  (ts2)             = (secti *)vrna_realloc((ts2), sizeof(secti) * (d->length));
                  (d->struct_list)  = (ts2);
                  (d->length)--;


                  d->structure =
                    (char *)vrna_realloc(d->structure, (d->length + 1) * sizeof(char));
                  d->structure[d->length] = '\0';
                  float new_en = (d->length > TURN) ? vrna_fold(d->sequence, d->structure) : 0.;
                  score += d->en - new_en;
                  d->en = new_en;
                  hit   = 1;
                  if (ticks > 10)
                    ticks--;

                  clearOutput();
                  getBoundary(d->struct_list, inv_bb);
                  putOnGridSeq(d->struct_list, d->sequence);
                  printHeader(d, stage, score, attempts, probe_max);
                  break;
                }

                d->struct_list++;
                ts++;
              }
              d->struct_list = ts2;
            }

            if (!hit)
              d->probe++;
          }
          d->probe = ts1;
          fflush(stdout);
          getBoundary(d->struct_list, inv_bb);

          ts1 = d->reactive;
          while ((*(d->reactive)).ml >= 0) {
            int i = (*(d->reactive)).i;
            int j = (*(d->reactive)).j;
            if ((j == 0) && (i >= hpos) && (i < (hpos + injector_len))) {
              switch ((*(d->reactive)).ml) {
                case 5:
                  attempts++;
                  break;
                case 6:
                  probe_max++;
                  break;
                default:
                  attempts--;
                  {
                    const char  *c;
                    int         ccc, cc[] = {
                      1, 3, 7
                    };
                    for (ccc = 0, c = flash; *c; c++, ccc++) {
                      printf("\033[23;%df\033[1;%d;40m%c\033[00m\033[0;0f",
                             i + 1,
                             30 + cc[ccc % 3],
                             *c);
                      fflush(stdout);
                      usleep(50000);
                    }
                    printf("\033[23;%df\033[1;37;40m\x20\033[00m\033[0;0f", i + 1);
                    fflush(stdout);
                  }
                  addFragments(d, i, j);
              }

              memmove((d->reactive),
                      (d->reactive) + 1,
                      ((d->reactive_n) - ((d->reactive) - (ts1))) * sizeof(secti));
              (ts1)         = (secti *)vrna_realloc((ts1), sizeof(secti) * (d->reactive_n));
              (d->reactive) = (ts1);
              (d->reactive_n)--;

              placeInjection(hpos);
              printHeader(d, stage, score, attempts, probe_max);
              continue;
            }

            d->reactive++;
          }
          d->reactive = ts1;
          if ((inv_bb[2] < 1) || (attempts <= 0)) {
            state = 4;
            break;
          }

          if (d->length == 0) {
            if (stage == 10) {
              state = 6;
              break;
            }

            clearOutput();
            stage++;
            ticks = init_ticks - (10 * stage);
            printf(lvlstr, stage);
            fflush(stdout);
            usleep(1000000);
            d = initStacks(d, stage);
            putOnGridSeq(d->struct_list, d->sequence);
            getBoundary(d->struct_list, inv_bb);
            clearOutput();
            printHeader(d, stage, score, attempts, probe_max);
          }

          fflush(stdout);
        }

        cnt2++;
        cnt2 %= ticks;
        fflush(stdout);
        break;


      case 3:
        break;


      case 4:
        last_state = 4;
        printf("%s", end);
        fflush(stdout);
        usleep(3000000);
        stage = 1;
        state = 0;
        break;


      case 5:
        break;


      case 6:
        last_state = 6;
        printf("%s", success);
        fflush(stdout);
        usleep(7000000);
        stage = 1;
        state = 0;
        break;

      default:
        break;
    }


    char c = 0;
    if (read(0, &c, sizeof(char)) > 0) {
      switch (c) {

        case 'q':
          goto runnaway;
          break;

        case 27:
          if (read(0, &c, sizeof(char)) > 0) {

            if (c == EOF)
              goto runnaway;


            if (c == 91) {
              if (read(0, &c, sizeof(char)) > 0) {
                if (state == 2) {

                  if (c == 68) {
                    if (hpos > 1)
                      hpos--;

                    printf("\033[23;%df\033[1;37;40m%c\033[0m", hpos + injector_len + 1, bg);
                  }

                  else if (c == 67) {
                    if (hpos < (79 - injector_len))
                      hpos++;

                    printf("\033[23;%df\033[1;37;40m%c\033[0m", hpos, bg);
                  }

                  else if (c == 65) {
                    if (d->probe_n < probe_max) {
                      (d->probe) =
                        (secti *)vrna_realloc((d->probe), sizeof(secti) * ((d->probe_n) + 2));
                      (d->probe)[(d->probe_n)].i    = (hpos + injector_len / 2);
                      (d->probe)[(d->probe_n)].j    = 0;
                      (d->probe)[(d->probe_n)].ml   = probe;
                      (d->probe)[++(d->probe_n)].ml = -1;
                    }
                  }


                  placeInjection(hpos);
                  printHeader(d, stage, score, attempts, probe_max);
                }
              }
            }
          }

          break;


        case  'p':
          if (state == 3) {
            state = last_state;
          } else if (state > 1) {
            last_state  = state;
            state       = 3;
          }

          break;


        case 'b':
          if (state == 5) {
            if (last_state == 2) {
              clearOutput();
              placeInjection(hpos);
              putOnGridSeq(d->struct_list, d->sequence);
              getBoundary(d->struct_list, inv_bb);
            }

            state = last_state;
          } else {
            last_state  = state;
            state       = 5;
          }

          break;


        case  'n':
          stage     = 1;
          attempts  = 3;
          probe_max = 3;
          d         = initStacks(d, stage);
          putOnGridSeq(d->struct_list, d->sequence);
          getBoundary(d->struct_list, inv_bb);
          clearOutput();
          placeInjection(hpos);
          printHeader(d, stage, score, attempts, probe_max);
          last_state  = state;
          state       = 2;
          break;


        case  'h':
          if (state == 1) {
            if (last_state == 2) {
              clearOutput();
              placeInjection(hpos);
              putOnGridSeq(d->struct_list, d->sequence);
              getBoundary(d->struct_list, inv_bb);
            }

            state = last_state;
          } else {
            printf("%s", start);
            last_state  = state;
            state       = 1;
          }

          break;


        case 's':
          if (probe_max < 50)
            probe_max++;

          break;
        case 'S':
          if (probe_max > 1)
            probe_max--;

          break;
        default:
          break;
      }
    }

    if (d) {

      ts1 = d->reactive;
      while ((*(d->reactive)).ml >= 0) {
        int i = (*(d->reactive)).i;
        int j = (*(d->reactive)).j;
        if ((j == 0) && (i >= hpos) && (i < (hpos + 5))) {
          switch ((*(d->reactive)).ml) {
            case 5:
              attempts++;
              break;
            case 6:
              probe_max++;
              break;
            default:
              attempts--;
          }
          memmove((d->reactive),
                  (d->reactive) + 1,
                  ((d->reactive_n) - ((d->reactive) - (ts1))) * sizeof(secti));
          (ts1)         = (secti *)vrna_realloc((ts1), sizeof(secti) * (d->reactive_n));
          (d->reactive) = (ts1);
          (d->reactive_n)--;
          printHeader(d, stage, score, attempts, probe_max);
          continue;
        }

        d->reactive++;
      }
      d->reactive = ts1;
      fflush(stdout);
    }

    printf("\033[0;0f");
    fflush(stdout);
  }
runnaway:
  tcsetattr(0, TCSANOW, &alt);
  printf("\033[0m\033[2J\033[0;0f");
  fflush(stdout);
  exit(0);
}

#endif
