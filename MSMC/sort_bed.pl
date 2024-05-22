#!/usr/bin/perl

use Term::ANSIColor qw(:constants);
chomp(@ARGV);
my $list = $ARGV[0];
if ($list !~ /bed/){
	print "Only accept a \"\*.bed\(\.gz\)\" file\n";
	exit;
}
if ($list !~ /gz$/){
	open(LIST, "<$list") || die BOLD "Cannot open $list: $!", RESET, "\n";
	$list = $list."\.gz";
}
else {
	open(LIST, "-|", "gzip -dc $list") || die BOLD "Cannot open $list: $!", RESET, "\n";
}
my @lists = <LIST>;
chomp(@lists);
close(LIST);

my $out = $list;
$out =~ s/bed/sorted.bed/;
open (OUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";

my $chr; my $last_end = -1;
foreach (@lists){
	@eles = split(/\t/, $_);
	chomp(@eles);
	if ($chr){}
	else {
		$chr = $eles[0];
	}
	if ($eles[0] ne $chr){
		print OUT "$last_end\n";
		$chr = $eles[0];
		$last_end = -1;
	}
	if ($eles[0] == $chr){
		if ($eles[2] < $eles[1]){
			next;
		}
		if ($eles[1] < $last_end && $eles[2] <= $last_end){
			next;
		}
		elsif ($eles[1] > $last_end){
			if ($last_end != -1){
				print OUT "$last_end\n";
			}	
			print OUT "$eles[0]\t$eles[1]\t";
			$last_end = $eles[2];
		}
		elsif ($eles[2] >= $last_end){
			$last_end = $eles[2];
		}
		else {
			print "exception found!\n";
		}
	}
}
print OUT "$last_end";
close(OUT);
close(LIST);