#!/usr/bin/perl
use Term::ANSIColor qw(:constants);

my $list = $ARGV[0];
my $repeats = $ARGV[1];

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

#count the least number of population samples, use half of them as the minimum number of the samples within a population in each run
my $min_num;
foreach my $key (keys %pop){
    my @samples = @{$pop{$key}};
    unless ($min_num){
        $min_num = scalar(@samples);
        next;
    }
    if (scalar(@samples) < $min_num){
        $min_num = scalar(@samples);
    }
}
$min_num = $min_num/2;
$min_num++ if ($min_num - int($min_num) >= 0.5);
$min_num = int($min_num);
$min_num = 1 if ($min_num < 1);

foreach my $l (1..$repeats){
	my $out = $list;
	$out =~ s/txt$|list$/$l.txt/;
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
    print OUT "sample pop\n";
	foreach my $key (keys %pop){
		my @samples = @{$pop{$key}};
		my @picked_samples;
		foreach my $j (0..$min_num-1){
			my $idx = int(rand($#samples));
			push(@picked_samples, $samples[$idx]);
			@samples = grep {!/\b$samples[$idx]\b/} @samples;
		}
		foreach my $sample (@picked_samples){
			print OUT "$sample $key\n";
		}
	}
	close(OUT);
}
