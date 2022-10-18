#/usr/bin/env python2
#

from __future__ import print_function
import sys
import math
import argparse

import RNA
import RNAxplorer

from pipeline.clusteralgorithms.diana import DIANA


granularity               = 1000
num_samples               = 10000
kt_fact                   = 1
do_clustering             = False
fake_2D_file              = False
# switch to select penalization of base pairs instead of computing complicated structure penalties
penalize_base_pairs       = True
verbose                   = False
debug                     = False
lmin_file                 = "local_minima.txt"
TwoD_file                 = "local_minima.2D.out"
nonredundant_sample_file  = False
nonredundant              = False
explore_two_neighborhood  = False
ediff_penalty             = False
mu                        = 0.1
post_filter_two           = False

# this is SV11
sequence = "GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA"
structure1="(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..)))."
structure2="(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....))))))))....."

"""
Error print
"""
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def RNAlocmin_output(sequence, m, filename = None):
    if filename:
        f = open(filename, 'w')
    else:
        f = sys.stdout

    f.write("     %s\n" % sequence)
    for i,s in enumerate(sorted(m.keys(), key=lambda x: m[x]['energy'])):
        f.write("%4d %s %6.2f %6d\n" % (i, s, m[s]['energy'], m[s]['count']))

    if filename:
        f.close()


def RNA2Dfold_output(sequence, s_mfe, mfe, s1, s2, m, filename = None):
    n         = len(sequence)
    distances = [ [ None for j in range(0, n) ] for i in range(0, n) ];

    for s in m.keys():
        d1 = RNA.bp_distance(s1, s)
        d2 = RNA.bp_distance(s2, s)
        if not distances[d1][d2] or m[s]['energy'] < distances[d1][d2]['e']:
            distances[d1][d2] = { 's': s, 'e': m[s]['energy'] }

    if filename:
        f = open(filename, 'w')
    else:
        f = sys.stdout

    f.write("%s\n%s (%6.2f)\n%s\n%s\n\n\n" % (sequence, s_mfe, mfe, s1, s2))

    for i in range(0, n):
        for j in range(0, n):
            if distances[i][j] != None:
                f.write("%d\t%d\ta\tb\tc\t%6.2f\td\t%s\n" % (i, j, distances[i][j]['e'], distances[i][j]['s']))

    if filename:
        f.close()


def mfeStructure(seq, cluster):
    mfe = sys.float_info.max
    mfe_s = ""
    for s in cluster:
        e = RNA.energy_of_struct(seq, s)
        if e < mfe:
            mfe = e
            mfe_s = s
    return mfe_s


def selectRepresentatives(seq, clusters):
    """
    extract mfe structures from clusters
    @param seq - string - the RNA sequence.
    @return list - structures in dot-bracket notation
    """
    representatives = []
    for c in clusters:
        # select mfe or TODO: centroid structure.
        s = mfeStructure(seq, c)
        representatives.append(s)
    return representatives


"""
Start argument parsing
"""
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbose", help="Be verbose", action="store_true")
parser.add_argument("-d", "--debug", help="Be even more verbose", action="store_true")
parser.add_argument("-s", "--sequence", type=str, help="Input sequence")
parser.add_argument("--struc1", type=str, help="Input structure 1")
parser.add_argument("--struc2", type=str, help="Input structure 2")
parser.add_argument("--granularity", type=int, help="Granularity, i.e. number of samples after which distortion checks are performed")
parser.add_argument("-n", "--num-samples", type=int, help="Number of samples per iteration")
parser.add_argument("-f", "--exploration-factor", type=float, help="Exploration factor")
parser.add_argument("--min-exploration-percent", type=float, help="Minimum exploration percentage before adding new repelled structures")
parser.add_argument("-c", "--cluster", help="Cluster resulting local minima to reduce effective size", action="store_true")
parser.add_argument("--lmin-file", type=str, help="Output filename for local minima")
parser.add_argument("--TwoD-file", type=str, help="Output filename for pseudo-2D file")
parser.add_argument("--nonred", help="Do sampling with non-redundant pbacktrack", action="store_true")
parser.add_argument("--nonred-file", type=str, help="Input filename for nonredundant samples")
parser.add_argument("-2", "--explore-two-neighborhood", help="Explore 2-Neighborhood of local minima, i.e. eliminate shallow minima", action="store_true")
parser.add_argument("--post-filter-two", help="Post processing Filter local minima according to 2-Neighborhood, i.e. eliminate shallow minima", action="store_true")
parser.add_argument("-e", "--ediff-penalty", help="Use energy difference instead of kT for penalty", action="store_true")
parser.add_argument("-m", "--mu", type=float, help="proportion factor used to decide whether sampling round was sufficient")

args = parser.parse_args()

if args.verbose:
    verbose = True

if args.debug:
    debug = True

if args.sequence:
    sequence = args.sequence

if args.struc1:
    structure1 = args.struc1

if args.struc2:
    structure2 = args.struc2

if args.granularity:
    granularity = args.granularity

if args.num_samples:
    num_samples = args.num_samples

if args.cluster:
    do_clustering = True

if args.exploration_factor:
    kt_fact = args.exploration_factor

if args.min_exploration_percent:
    min_explore_min_percent = args.min_exploration_percent

if args.lmin_file:
    lmin_file = args.lmin_file

if args.TwoD_file:
    TwoD_file = args.TwoD_file
    fake_2D_file = True

if args.nonred_file:
    nonredundant_sample_file = args.nonred_file

if args.nonred:
    nonredundant = True

if args.explore_two_neighborhood:
    explore_two_neighborhood = True

if args.ediff_penalty:
    ediff_penalty = True

if args.mu:
    mu  = args.mu

if args.post_filter_two:
    post_filter_two = True


sc_data = {
  'base_pairs': {},
  'weights': {},
}


def store_basepair_sc(fc, data, structure, weight, distance_based = False):
    if distance_based:
        return
    else:
        pt = RNA.ptable(structure)
        # count number of pairs in structure to repell
        cnt = 0
        for i in range(1, len(pt)):
            if pt[i] > i:
                cnt = cnt + 1

        if cnt > 0:
            weight = weight / cnt

        # add repulsion
        for i in range(1, len(pt)):
            if pt[i] > i:
                key = (i, pt[i])
                if key not in data['base_pairs']:
                    data['base_pairs'][key] = 1
                    data['weights'][key] = weight
                else:
                    data['base_pairs'][key] = data['base_pairs'][key] + 1
                    data['weights'][key] = data['weights'][key] + weight

        # remove previous soft constraints
        fc.sc_remove()

        # add latest penalties for base pairs
        for k in data['weights'].keys():
            i = k[0]
            j = k[1]
            fc.sc_add_bp(i, j, data['weights'][k])


def move_apply(structure_in, move):
    #s_length = structure_in[0]
    output_structure = list(structure_in)
    if(move.pos_5 > 0 and move.pos_3 > 0):
        output_structure[move.pos_5] = move.pos_3
        output_structure[move.pos_3] = move.pos_5
    elif(move.pos_5 < 0 and move.pos_3 < 0):
        output_structure[-move.pos_5] = 0
        output_structure[-move.pos_3] = 0
    else:
        print("Error: shfit moves are not supported")
    res = output_structure
    return res


def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))


def detect_local_minimum_two(fc, structure):
    """
    Take a local minimum and detect nearby local minimum via extended gradient walks.
    If a lower energy structure within radius 2 is detected (lower than the local minimum
    of a normal gradient walk), then another gradient walk is applied for the lower structure.
    """
    deepest_neighbor      = None
    deepest_neighbor_ddG  = 1
    pt                    = RNA.IntVector(RNA.ptable(structure))
    neigh                 = fc.neighbors(pt, RNA.MOVESET_DELETION | RNA.MOVESET_INSERTION)

    for nb in neigh:
        dG_nb   = fc.eval_move_pt(pt, nb.pos_5, nb.pos_3)
        pt_nb   = RNA.IntVector(move_apply(pt, nb))
        neigh_2 = fc.neighbors(pt_nb, RNA.MOVESET_DELETION | RNA.MOVESET_INSERTION)

        for nb_2 in neigh_2:
            dG_nb2  = fc.eval_move_pt(pt_nb, nb_2.pos_5, nb_2.pos_3)
            ddG     = dG_nb + dG_nb2
            if ddG < 0 and ddG < deepest_neighbor_ddG:
                deepest_neighbor_ddG  = ddG
                deepest_neighbor      = move_apply(pt_nb, nb_2)

    if deepest_neighbor:
        deepest_neighbor = detect_local_minimum(fc, RNA.db_from_ptable(deepest_neighbor))
        return detect_local_minimum_two(fc, deepest_neighbor)
    else:
        return structure


def reduce_lm_two_neighborhood(fc, lm, verbose = False):
    cnt       = 1;
    cnt_max   = len(lm)
    lm_remove = list()
    lm_novel  = dict()

    if verbose:
        sys.stderr.write("Applying 2-Neighborhood Filter...")
        sys.stderr.flush()

    for s in lm:
        if verbose:
            sys.stderr.write("\rApplying 2-Neighborhood Filter...%6d / %6d" % (cnt, cnt_max))
            sys.stderr.flush()

        ss = detect_local_minimum_two(fc, s)
        if s != ss:
            # store structure for removal
            lm_remove.append(s)
            if ss not in lm:
                # store structure for novel insertion
                if ss not in lm_novel:
                    lm_novel[ss] = { 'count' : lm[s]['count'], 'energy' : fc.eval_structure(ss) }
                else:
                    lm_novel[ss]['count'] = lm_novel[ss]['count'] + lm[s]['count']
            else:
                # update current minima list
                lm[ss]['count'] = lm[ss]['count'] + lm[s]['count']

        cnt = cnt + 1

    # remove obsolete local minima
    for a in lm_remove:
        del lm[a]

    # add newly detected local minima
    lm.update(lm_novel)

    if verbose:
        sys.stderr.write("\rApplying 2-Neighborhood Filter...done             \n")


def generate_samples(fc, number, non_redundant=False):
    samples = list()

    if non_redundant:
        samples = fc.pbacktrack_nr(number)
    else:
        for i in range(0, number):
            samples.append(fc.pbacktrack())

    return samples


"""
Do main stuff
"""
# init random number generator in RNAlib
RNA.init_rand()

# prepare RNAlib fold_compound
md = RNA.md()
md.uniq_ML = 1
md.compute_bpp = 0

kT = RNA.exp_param(md).kT

fc = RNA.fold_compound(sequence, md)
fc_base = RNA.fold_compound(sequence, md)

# compute MFE structure
(ss, mfe) = fc.mfe()
fc.pf()


minima = dict()
sample_list = []

pending_lm = dict()
num_sc = 1

num_iter      = int(math.ceil(float(num_samples) / float(granularity)))
samples_left  = num_samples

for it in range(0, num_iter):
    # determine number of samples for this round
    # usually, this is 'granularity'
    if samples_left < granularity:
        current_num_samples = samples_left
    else:
        current_num_samples = granularity

    samples_left = samples_left - current_num_samples

    # generate samples through stocastic backtracing
    sample_set = generate_samples(fc, current_num_samples, nonredundant)

    sys.stderr.write("\rsamples so far: %6d / %6d" % (num_samples - samples_left, num_samples))
    sys.stderr.flush()

    # store samples of this round to global list of samples
    sample_list = sample_list + sample_set

    current_lm = dict()

    # go through list of sampled structures and determine corresponding local minima
    for s in sample_set:
        ss = detect_local_minimum(fc_base, s)

        if ss not in current_lm:
            current_lm[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
        else:
            current_lm[ss]['count'] = current_lm[ss]['count'] + 1

    # explore the 2-neighborhood of current local minima to reduce total number of local minima
    if explore_two_neighborhood:
        reduce_lm_two_neighborhood(fc_base, current_lm, verbose)

    # transfer local minima obtained in this iteration to list of pending local minima
    for ss in current_lm:
        if ss not in pending_lm:
            pending_lm[ss] = current_lm[ss]
        else:
            pending_lm[ss]['count'] = pending_lm[ss]['count'] + current_lm[ss]['count']

    del current_lm

    if it < num_iter - 1:
        # find out which local minima we've seen the most in this sampling round
        struct_cnt_max = max(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['count']))
        struct_cnt_min = min(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['count']))
        struct_en_max = max(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['energy']))
        struct_en_min = min(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['energy']))

        cnt_once  = 0
        cnt_other = 0

        for key, value in sorted(pending_lm.iteritems(), key=lambda (k, v): v['energy']):
            if value['count'] == 1:
                cnt_once = cnt_once + 1
            else:
                cnt_other = cnt_other + 1

            # check whether we've seen other local minim way more often than those we've seen just once
            if cnt_other > 0 and cnt_once > 0:
                if cnt_once <= (mu * cnt_other):
                    if ediff_penalty:
                        repell_en = kt_fact * (value['energy'] - pending_lm[struct_en_min]['energy'])
                    else:
                        repell_en = kt_fact * kT / 1000.

                    store_basepair_sc(fc, sc_data, struct_cnt_max, repell_en)

                    fc.pf()

                    for cmk in pending_lm.keys():
                        if cmk not in minima:
                            minima[cmk] = pending_lm[cmk]
                        else:
                            minima[cmk]['count'] = minima[cmk]['count'] + pending_lm[cmk]['count']

                    pending_lm = dict()

                    break
    else:
        for cmk in pending_lm.keys():
            if cmk not in minima:
                minima[cmk] = pending_lm[cmk]
            else:
                minima[cmk]['count'] = minima[cmk]['count'] + pending_lm[cmk]['count']

eprint(" ... done")

if post_filter_two:
    reduce_lm_two_neighborhood(fc_base, minima, verbose)


RNAlocmin_output(sequence, minima, lmin_file)


sample_file=""
rind = lmin_file.rfind(".")
if rind >= 0 :
    sample_file = lmin_file[:rind] + ".samples"
else:
    sample_file = lmin_file[:rind] + ".samples"
f = open(sample_file, 'w')
f.write("     %s\n" % sequence)
for s in sample_list:
    f.write(s+"\n")
f.close()


if fake_2D_file:
    (ss, mfe) = fc_base.mfe()
    RNA2Dfold_output(sequence, ss, mfe, structure1, structure2, minima, TwoD_file)


# read a list of sample structures and produce list of local minima for it
if nonredundant_sample_file:
    lmin_nonred_file = "local_minima_nonred.txt"
    nonredundant_samples = []
    with open(nonredundant_sample_file) as f:
        nonredundant_samples = f.readlines()

    nonredundant_samples = [x.strip() for x in nonredundant_samples]

    nonredundant_samples.pop(0)

    nonredundant_minima = dict()

    for s in nonredundant_samples:
        pt = RNA.IntVector(RNA.ptable(s))
        fc_base.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
        ss = RNA.db_from_ptable(list(pt))
        if ss not in nonredundant_minima:
             nonredundant_minima[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
             new_minima = new_minima + 1
        else:
             nonredundant_minima[ss]['count'] = nonredundant_minima[ss]['count'] + 1

    f = open(lmin_nonred_file, 'w')
    f.write("     %s\n" % sequence)
    for i,s in enumerate(sorted(nonredundant_minima.keys(), key=lambda x: nonredundant_minima[x]['energy'])):
        f.write("%4d %s %6.2f %6d\n" % (i, s, nonredundant_minima[s]['energy'], nonredundant_minima[s]['count']))
    f.close()

    if fake_2D_file:
        distances = [ [ None for j in range(0, 200) ] for i in range(0, 200) ];

        for s in nonredundant_minima.keys():
            d1 = RNA.bp_distance(structure1, s)
            d2 = RNA.bp_distance(structure2, s)
            if not distances[d1][d2] or nonredundant_minima[s]['energy'] < distances[d1][d2]:
                distances[d1][d2] = nonredundant_minima[s]['energy']

        f = open("sv11_fake_nonred.2D.out", 'w')
        f.write("%s\n%s (%6.2f)\n%s\n%s\n\n\n" % (sequence, structure1, mfe, structure1, structure2))
        for i in range(0, 200):
            for j in range(0, 200):
                if distances[i][j] != None:
                    f.write("%d\t%d\ta\tb\tc\t%6.2f\n" % (i, j, distances[i][j]))
        f.close()






if do_clustering:
    initial_minima = [ s for s in minima.keys() ]

    maxDiameterThreshold = 0
    maxAverageDiameterThreshold = 6

    clusters = DIANA.doClustering(initial_minima, maxDiameterThreshold, maxAverageDiameterThreshold)
    #DIANA.printClusters(clusters)

    if verbose:
        eprint("done clustering %d structures into %d clusters" % (len(initial_minima), len(clusters)))

    representatives = selectRepresentatives(sequence, clusters)

    # create a minima-like dict for output printing
    red_minima = dict()
    for r in representatives:
        red_minima[r] = {'count': 1, 'energy' : fc_base.eval_structure(r)}

    # produce RNAlocmin - like output
    #RNAlocmin_output(sequence, minima)
    RNAlocmin_output(sequence, red_minima)

