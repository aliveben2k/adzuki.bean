#!/usr/bin/perl
use Term::ANSIColor qw(:constants);

my $hap = "NA";
my $list = $ARGV[0];
my $repeats = $ARGV[1];
my $sample_num = $ARGV[2];
my $hap = $ARGV[3];
if ($ARGV[3] eq "hap"){
	$hap = 1;
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
	$out =~ s/poplabels$/$l.poplabels/;
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
    print OUT "sample population group sex\n";
	foreach my $key (keys %pop){
		my @samples = @{$pop{$key}};
		if (scalar(@samples) < $sample_num){
			die "Sample number is not enought for $pop.\n";
		}
		my @picked_samples;
		foreach my $j (0..$sample_num-1){
			my $idx = int(rand($#samples));
			push(@picked_samples, $samples[$idx]);
			@samples = grep {!/\b$samples[$idx]\b/} @samples;
		}
		foreach my $sample (@picked_samples){
			print OUT "$sample $key $key $hap\n";
		}
	}
	close(OUT);
}
