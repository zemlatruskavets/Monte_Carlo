#!/usr/bin/perl

# Usage: Initial <N> <L> <n>
#
# Create an XYZ file with N particles in a cubic box of length L with reduced density n

$N  = shift(@ARGV);
$L  = shift(@ARGV);
$n  = shift(@ARGV);

use constant PI => 4 * atan2(1, 1);

$N_C = $N ** (1 / 3);
$D = ((6 * $n) / ($N * PI)) ** (1 / 3);

print "N = $N\nCube root = $N_C\nBox length = $L\nDiameter = $D\n\n";

for ($i = 0; $i < $N_C; $i++)
{
    $x[$i] = $D * $i;
    for ($j = 0; $j < $N_C; $j++)
    {
        $y[$j] = $D * $j;
        for ($k = 0; $k < $N_C; $k++)
        {
            $z[$k] = $D * $k;
            print "$x[$i] $y[$j] $z[$k]\n";
        }
    }

}

