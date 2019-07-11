#!/usr/bin/perl

# Usage: Initial_random <N> <L> <n>
#
# Create an XYZ file with N particles in a cubic box of length L with reduced density n

$N  = shift(@ARGV);
$L  = shift(@ARGV);
$n  = shift(@ARGV);

use constant PI => 4 * atan2(1, 1);

$N_C = $N ** (1 / 3);
$D = ((6 * $n) / ($N * PI)) ** (1 / 3);

#print "N = $N\nCube root = $N_C\nBox length = $L\nDiameter = $D\n\n";


for ($i = 0; $i < $N; $i++)
{
    do
    {
        do
        {
            $x[$i] = $L * (rand() * 1) + ($L / 4.);
            $y[$i] = $L * (rand() * 1) + ($L / 4.);
            $z[$i] = $L * (rand() * 1) + ($L / 4.);
        }
        while ($x[$i] * $x[$i] + $y[$i] * $y[$i] + $z[$i] * $z[$i] == $L * $L);
        
        $accept = 1;
        for ($j = 0; $j < $i - 1; $j++)
        {
            $dx = $x[$i] - $x[$j];
            $dy = $y[$i] - $y[$j];
            $dz = $z[$i] - $z[$j];
            $r2 = $dx * $dx + $dy * $dy + $dz * $dz;
            if ($r2 < ($L * $L)) {$accept = 0; break;}
        }
    } while ($accept == 0);
    
    print "$x[$i] $y[$i] $z[$i]\n";
}
