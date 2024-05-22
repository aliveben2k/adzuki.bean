#!/usr/local/bin/perl

chomp(@ARGV);
my $o_path = $ARGV[0];
my $ran = $ARGV[1]; 
my $num = $ARGV[2];
#modify *.haps file
open(HAP,"<$o_path\/inputs\/$ran\_input_chr$num.haps") || die "Cannot open $o_path\/inputs\/$ran\_input_chr$num.haps: $!\n";
my @lines = <HAP>;
close(HAP);
chomp(@lines);
open(OUT, ">$o_path\/inputs\/$ran\_input_chr$num.haps") || die "Cannot write $o_path\/inputs\/$ran\_input_chr$num.haps: $!\n";
foreach my $i (0..$#lines){
	my @eles = split(/\s+|\t+/, $lines[$i]);
	foreach my $j (0..4){
		print OUT "$eles[$j] ";
	}
	for (my $k=5; $k<=$#eles; $k+=2){
		if ($k < $#eles-1){
			print OUT "$eles[$k] ";
		}
		else {
			print OUT "$eles[$k]\n";
		}
	}
}
close(OUT);