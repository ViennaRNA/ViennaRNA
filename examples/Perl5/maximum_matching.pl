use strict;
use warnings;
use Data::Dumper;
use RNA;

my $seq1 = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCAA";

# Turn-off dangles globally
$RNA::dangles = 0;

# Data structure that will be passed to our MaximumMatching() callback with two components:
# 1. a 'dummy' fold_compound to evaluate loop energies w/o constraints, 2. a fresh set of energy parameters
my %mm_data = ( 'dummy' => new RNA::fold_compound($seq1), 'params' => new RNA::param() );

# Nearest Neighbor Parameter reversal functions
my %revert_NN = ( 
    RNA::DECOMP_PAIR_HP => sub { my ($i, $j, $k, $l, $f, $p) = @_; return - $f->eval_hp_loop($i, $j) - 100;},
    RNA::DECOMP_PAIR_IL => sub { my ($i, $j, $k, $l, $f, $p) = @_; return - $f->eval_int_loop($i, $j, $k, $l) - 100},
    RNA::DECOMP_PAIR_ML => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - $p->{MLclosing} - $p->{MLintern}[0] - ($j - $i - $k + $l - 2) * $p->{MLbase} - 100},
    RNA::DECOMP_ML_ML_STEM => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - $p->{MLintern}[0] - ($l - $k - 1) * $p->{MLbase}},
    RNA::DECOMP_ML_ML_ML => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  0},
    RNA::DECOMP_ML_STEM => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - $p->{MLintern}[0] - ($j - $i - $k + $l) * $p->{MLbase}},
    RNA::DECOMP_ML_ML => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - ($j - $i - $k + $l) * $p->{MLbase}},
    RNA::DECOMP_ML_UP => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - ($j - $i + 1) * $p->{MLbase}},
    RNA::DECOMP_EXT_STEM => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - $f->E_ext_loop($k, $l)},
    RNA::DECOMP_EXT_EXT => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  0},
    RNA::DECOMP_EXT_STEM_EXT => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - $f->E_ext_loop($i, $k)},
    RNA::DECOMP_EXT_EXT_STEM => sub { my ($i, $j, $k, $l, $f, $p) = @_; return : - $f->E_ext_loop($l, $j)},
    RNA::DECOMP_EXT_EXT_STEM1 => sub { my ($i, $j, $k, $l, $f, $p) = @_; return  - $f->E_ext_loop($l, $j - 1)},
);

# Maximum Matching callback function (will be called by RNAlib in each decomposition step)
sub MaximumMatching {
  my ($i, $j, $k, $l, $d, $data) = @_;
  return $revert_NN{$d}->($i, $j, $k, $l, $data->{'dummy'}, $data->{'params'}) if defined $revert_NN{$d};
  return 0;
}

# Create a 'fold_compound' for our sequence
my $fc = new RNA::fold_compound($seq1);

# Add maximum matching soft-constraints
$fc->sc_add_f(\&MaximumMatching);
$fc->sc_add_data(\%mm_data, undef);

# Call MFE algorithm
my ($s, $mm) = $fc->mfe();

# print result
printf("%s\n%s (MM: %d)\n", $seq1, $s, - $mm);

