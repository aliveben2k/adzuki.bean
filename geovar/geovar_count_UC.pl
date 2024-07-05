#!/usr/bin/perl
use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
my $title;
if ($ARGV[0] == -1){
	print "Usage: perl geovar_count.pl frequency_table.csv\n";
	exit;
}
if (-e $ARGV[0]){}
else {
	print "Usage: perl geovar_count.pl frequency_table.csv\n";
	exit;
}

open(OUTPUT, ">$ARGV[0].cnt") || die BOLD "Cannot write $ARGV[0].cnt: $!", RESET, "\n";
my $U = 0; my @types; my $total = 0;
open(INPUT, "<$ARGV[0]") || die BOLD "Cannot open $ARGV[0]: $!", RESET, "\n";
my @samples; my @header;
while (my $line = <INPUT>){
	chomp($line);
	if ($line =~ /^CHR/){
		print OUTPUT "CNT\tCFV\t";
		@header = split(/\s/, $line);
		for (my $i=6; $i<=$#header; $i++){
			if ($i<$#header){
				print OUTPUT "$header[$i]\t";
			}
			else {
				print OUTPUT "$header[$i]\n";
			}
			if ($i >= 6){
				push(@samples, $header[$i]);
			}
		}	
		
	}
	else {
		$total++;
		my $tmp_type;
		my @values = split(/\s/, $line);
		for (my $j=6; $j<=$#header; $j++){
			if ($values[$j] == 0){
				$tmp_type .= "U ";
			}
			elsif ($values[$j] > 0 && $values[$j] <= 1){
				$tmp_type .= "C ";
			}
			else {
				$tmp_type .= "- ";
			}
		}
		$tmp_type =~ s/\s$//;
		if ($tmp_type !~ /C/){
			$U++;
			next;
		}
		if (@types){
			my $true = 0;
			for (my $k=0; $k<=$#types; $k++){
				if ($types[$k] eq $tmp_type){
					$counts[$k]++;
					$true = 1;
				}
			}
			if ($true == 0){
				push(@types, $tmp_type);
				push(@counts, 1);
			}
		}
		else {
			@types = $tmp_type;
			@counts = 1;
		}
	}
}
close(INPUT);
my @sort_idx = sort {@counts[$b] <=> @counts[$a]} 0..$#counts;
my @stypes = @types[@sort_idx];
my @scounts = @counts[@sort_idx];
$total = $total - $U;
for (my $l=0; $l<=$#stypes; $l++){
	$stypes[$l] =~ s/\s/\t/g;
	my $cfv = $scounts[$l]/$total;
	print OUTPUT "$scounts[$l]\t$cfv\t$stypes[$l]\n";
}
close(OUTPUT);

open(OUT2, ">$ARGV[0].Rplot.txt") || die BOLD "Cannot write $ARGV[0].Rplot.txt: $!", RESET, "\n";
print OUT2 "X\tY\tCNT\tFV\tCFV\tGroup\n";
open(OUT3, ">$ARGV[0].Rplot_private.txt") || die BOLD "Cannot write $ARGV[0].Rplot_private.txt: $!", RESET, "\n";
print OUT3 "X\tY\tCNT\tFV\tGroup\n";
my $cfv = 0;
for (my $n=0; $n<=$#stypes; $n++){
	my $Ccount=0; my $Mcount=0;
	my @privates;
	my @blocks = split(/\t/, $stypes[$n]);
	my $fv = $scounts[$n]/$total * 100;
	$cfv += $fv;
	my $y_axis = $n + 1;
	for (my $o=0; $o<=$#blocks; $o++){
		print OUT2 "$samples[$o]\t$y_axis\t$scounts[$n]\t";
		push(@privates, "$samples[$o]\t$y_axis\t$scounts[$n]\t");
		printf OUT2 "%.4f\t", $fv;
		push(@privates, "$fv\t");
		printf OUT2 "%.2f\t", $cfv;
		if ($blocks[$o] =~ /U/){
			print OUT2 "U\n";
			push(@privates, "U\n");
		}
		elsif ($blocks[$o] =~ /C/){
			print OUT2 "C\n";
			push(@privates, "C\n");
			$Ccount++;
		}
		else {
			print OUT2 "-\n";
			$Mcount++;
		}
	}
	my $private = $Ccount;
	if ($private == 1 && $Mcount == 0){
		print OUT3 @privates;
	}
}
close(OUT2);
close(OUT3);
print "Total counts: $total\n";
my $type_n = $#stypes+1;
print "Uni-type number: $type_n\n";
