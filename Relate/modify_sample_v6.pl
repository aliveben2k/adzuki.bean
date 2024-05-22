#!/usr/local/bin/perl

chomp(@ARGV);
my $o_path = $ARGV[0];
my $ran = $ARGV[1]; 
my $num = $ARGV[2];
#modify *.sample file
open(SAM,"<$o_path\/inputs\/$ran\_input_chr$num.sample") || die "Cannot open $o_path\/inputs\/$ran\_input_chr$num.sample: $!\n";
my @lines = <SAM>;
close(SAM);
chomp(@lines);
open(OUT, ">$o_path\/inputs\/$ran\_input_chr$num.sample") || die "Cannot write $o_path\/inputs\/$ran\_input_chr$num.sample: $!\n";
print OUT "$lines[0]\n$lines[1]\n";
foreach my $i (2..$#lines){
	my @eles = split(/\s+|\t+/, $lines[$i]);
	$eles[1] = "NA";
	print OUT join("\t", @eles), "\n";
}
close(OUT);