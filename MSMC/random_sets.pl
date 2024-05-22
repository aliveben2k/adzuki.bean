#!/usr/bin/perl
use Term::ANSIColor qw(:constants);

my $syn_num = 2;
my $list = $ARGV[0];
my $repeats = $ARGV[1];
my $syn = $ARGV[2];
if ($ARGV[3]){
	$syn_num = $ARGV[3];
}

open(IN, "<$list") || die "Cannot open $list: $!\n";
my @lists = <IN>;
chomp(@lists);
close(IN);
shift(@lists);

my %pop;
foreach my $i (0..$#lists){
	my @eles = split(/\t|\s+/, $lists[$i]);
	if (@{$pop{$eles[1]}}){
		push(@{$pop{$eles[1]}}, $eles[0]);
	}
	else {
		@{$pop{$eles[1]}} = $eles[0];
	}
}

foreach my $l (1..$repeats){
	my $out = $list;
	$out =~ s/txt$|list$/$l.txt/;
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
	foreach my $key (keys %pop){
		my @samples = @{$pop{$key}};
		my $sample_num;
		if ($syn == 0){
			$sample_num = $syn_num;
		}
		else {
			$sample_num = $syn_num*2;
		}
		if (scalar(@samples) < $sample_num){
			die "Sample number is not enought for $pop.\n";
		}
		my @picked_samples;
		foreach my $j (0..$sample_num-1){
			my $idx = int(rand($#samples));
			push(@picked_samples, $samples[$idx]);
			@samples = grep {!/\b$samples[$idx]\b/} @samples;
		}
		if ($syn == 0){
			foreach my $sample (@picked_samples){
				print OUT "$sample\t$key\n";
			}
		}
		else {
			my $syn_cnt = 1;
			for (my $k=0; $k<=$#picked_samples; $k+=2){
				print OUT "$picked_samples[$k]\t$picked_samples[$k+1]\t$key\_$syn_cnt\t$key\n";
				$syn_cnt++;
			}
		}
	}
	close(OUT);
}
