#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
#use PerlIO::gzip;
my $time = scalar localtime();

chomp(@ARGV);
my $input = $ARGV[0];
my $dp = 3;
my $out = $input;
if ($out =~ /raw\.vcf/){
	$out =~ s/raw\.vcf/GTonly\.vcf/;
}
elsif ($out =~ /filtered\.vcf/){
    $out =~ s/filtered\.vcf/GTonly\.vcf/;
}
else {
	$out =~ s/\.vcf/\.GTonly\.vcf/;
}

print "\[$time\]\: Filtering start...\n";
if ($input =~ /\.gz$/){
	open(INPUT, "-|", "gzip -dc $input") || die BOLD "Cannot open $input: $!", RESET, "\n";
}
else {
	open (INPUT, "<$input") || die BOLD "Cannot open $input: $!", RESET, "\n";
	$out = $out."\.gz";
}
open (OUTPUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
my $num = 1; 
my @sample_two_alleles;
my $total = 0;
while (my $line = <INPUT>){
	my @ALT_alleles_noAster; my $which_allele_aster;
	chomp $line;
	if ($line =~ /^#/){
		print OUTPUT "$line\n";
	}
	else {
		$total += 1;
		my @line_vec = split('\t', $line);
		my @ALT_alleles = split(',', $line_vec[4]);
		for (my $i=0; $i<=$#ALT_alleles; $i++){
			if ($ALT_alleles[$i] !~ /\*/){
				push(@ALT_alleles_noAster, $ALT_alleles[$i]);
			}
			else {
				$which_allele_aster = $i+1;
			}
		}
		if ($#ALT_alleles_noAster == -1){
            next;
		}
		my $ALT_join = join(',', @ALT_alleles_noAster);
		if ($line !~ /\*/){
			$which_allele_aster = 500;
		}
		for (my $j=0; $j<=8; $j++){
			if ($j == 4){
				$line_vec[$j] = $ALT_join;
			}
			if ($j == 8){
				$line_vec[$j] = "GT";
			}
			if ($j == 0){
                if ($line_vec[4] ne "."){
                    print OUTPUT "$line_vec[$j]";
				}
			}
			else {
                if ($line_vec[4] ne "."){
                    print OUTPUT "\t$line_vec[$j]";
				}
			}
		}
		for (my $l=9; $l<=$#line_vec; $l++){
			my @sample_vec;
			@sample_vec = split(':', $line_vec[$l]);
			@sample_two_alleles = split(/\/|\|/, $sample_vec[0]);
			if (scalar(@sample_two_alleles) >= 2){
				$sample_vec[0] = "$sample_two_alleles[0]\|$sample_two_alleles[1]";
			}
			elsif ($sample_two_alleles[0] eq '.'){
				$sample_vec[0] = ".\|.";
			}
			else {
				$sample_vec[0] = ".\|.";
			}
            if ($line_vec[4] ne "."){
                print OUTPUT "\t$sample_vec[0]";
			}
		}
		if ($line_vec[4] ne "."){
            print OUTPUT "\n";
		}
	}
}
close(INPUT);
close(OUTPUT);
$time = scalar localtime();
print "\[$time\]\: Filter_vcf done.\n";
