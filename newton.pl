#################
### newton.pl ###
#################

# Compute the Newton polytope of f_t,d,n using the Polymake software package,
# and count the number of vertices.

# Takes as input the file produced by computeSupport.m2.

use application 'polytope';

open MACAULAY_OUTPUT, "computeSupport.out";

my @points = ();

while(<MACAULAY_OUTPUT>) {
  chomp;
  s/\}//;
  s/\{//;
  push @points, [1, split ", "];
}

print scalar @points, "\n";
print scalar @{$points[0]}, "\n";

my $newton_polytope = new Polytope(POINTS=>\@points);
my $vertices = $newton_polytope->VERTICES;
print $vertices, "\n";
print "Number of vertices: ", $vertices->rows(), "\n";
print "Number of lattice points: ", $newton_polytope->N_LATTICE_POINTS, "\n";

