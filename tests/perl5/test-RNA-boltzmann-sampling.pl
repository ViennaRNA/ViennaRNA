use strict;
use warnings;
use Test::More tests => 32;

use RNA;

my $sequence = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUACUUAUACCGCCUGUGCGGACUACUAUCCUGACCACAUAGU";

sub store_structure {
  my $s     = shift;
  my $data  = shift;
  push @{$data}, $s if defined($s);
}


sub prepare_fc {
  my $md = new RNA::md();
  $md->{uniq_ML} = 1;
  my $fc = new RNA::fold_compound($sequence, $md);

  my ($ss, $mfe) = $fc->mfe();

  $fc->exp_params_rescale($mfe);
  $fc->pf();

  return $fc;
}


sub uniq {
  my %seen;
  grep !$seen{$_}++, @_;
}


####################
# Begin actual tests
####################
my $fc          = prepare_fc();
my $failure     = 0;
my $num_samples = 500;
my $iterations  = 15;
my $s;
my $i;
my @ss;
my @sss;
my $d;

print "test_pbacktrack5    (single sub-structure, a.k.a. pbacktrack5)\n";
$s = $fc->pbacktrack5(10);
ok(length($s) == 10);

$s = $fc->pbacktrack5(50);
ok(length($s) == 50);


print "test_pbacktrack5    (multiple sub-structures, a.k.a. pbacktrack5_num)\n";
$s = $fc->pbacktrack5(20, 10);

ok(scalar(@{$s}) == 20);

$failure = 0;
foreach my $s (@{$s}) {
  if (length($s) != 10) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);

$s = $fc->pbacktrack5(100, 50);
ok(scalar(@{$s}) == 100);

$failure = 0;
foreach my $s (@{$s}) {
  if (length($s) != 50) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);


print "test_pbacktrack     (single structure, a.k.a. pbacktrack)\n";
$s = $fc->pbacktrack();
ok(length($s) == length($sequence));


print "test_pbacktrack     (multiple structures, a.k.a. pbacktrack_num)\n";
$s = $fc->pbacktrack(100);
ok(scalar(@{$s}) == 100);

$failure = 0;
foreach my $s (@{$s}) {
  if (length($s) != length($sequence)) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);


print "test_pbacktrack_nr\n";
$s = $fc->pbacktrack(100, RNA::PBACKTRACK_NON_REDUNDANT);
ok(scalar(@{$s}) == 100);

$failure = 0;
foreach my $s (@{$s}) {
  if (length($s) != length($sequence)) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);

# check for uniqueness, i.e. no duplicates
@ss = uniq(@{$s});
ok(scalar(@{$s}) == scalar(@ss));


print "test_pbacktrack_nr_resume\n";
$d = undef;
@ss = ();

foreach my $i (1..$iterations) {
  ($d, $s) = $fc->pbacktrack($num_samples, $d, RNA::PBACKTRACK_NON_REDUNDANT);
  push(@ss, @{$s});
}

ok(scalar(@ss) == ($iterations * $num_samples));

$failure = 0;
foreach my $s (@ss) {
  if (length($s) != length($sequence)) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);

# check for uniqueness, i.e. no duplicates
@sss = uniq(@ss);
ok(scalar(@ss) == scalar(@sss));


print "test_pbacktrack5_cb (multiple sub-structures, a.k.a. pbacktrack5_cb)\n";
@ss = ();
$i  = $fc->pbacktrack5(20, 10, \&store_structure, \@ss);
ok($i == 20);
ok(scalar(@ss) == 20);

$failure = 0;
foreach my $s (@ss) {
  if (length($s) != 10) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);

@ss = ();
$i  = $fc->pbacktrack5(100, 50, \&store_structure, \@ss);
ok($i == 100);
ok(scalar(@ss) == 100);

$failure = 0;
foreach my $s (@ss) {
  if (length($s) != 50) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);


print "test_pbacktrack_cb  (multiple structures, a.k.a. pbacktrack_cb)\n";
@ss = ();
$i  = $fc->pbacktrack(100, \&store_structure, \@ss);
ok($i == 100);
ok(scalar(@ss) == 100);

$failure = 0;
foreach my $s (@ss) {
  if (length($s) != length($sequence)) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);


print "test_pbacktrack_nr_cb\n";
@ss = ();
$i  = $fc->pbacktrack(100, \&store_structure, \@ss, RNA::PBACKTRACK_NON_REDUNDANT);
ok($i == 100);
ok(scalar(@ss) == 100);

$failure = 0;
foreach my $s (@ss) {
  if (length($s) != length($sequence)) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);

# check for uniqueness, i.e. no duplicates
@sss = uniq(@ss);
ok(scalar(@ss) == scalar(@sss));


print "test_pbacktrack_nr_resume_cb\n";
$d  = undef;
@ss = ();

$failure = 0;
foreach my $i (1..$iterations) {
  ($d, $i) = $fc->pbacktrack($num_samples, \&store_structure, \@ss, $d, RNA::PBACKTRACK_NON_REDUNDANT);
  $failure = 1, last if $i != $num_samples;
}
ok($failure == 0);

ok(scalar(@ss) == ($iterations * $num_samples));

$failure = 0;
foreach my $s (@ss) {
  if (length($s) != length($sequence)) {
    $failure = 1;
    last;
  }
}
ok($failure == 0);

# check for uniqueness, i.e. no duplicates
@sss = uniq(@ss);
ok(scalar(@ss) == scalar(@sss));
