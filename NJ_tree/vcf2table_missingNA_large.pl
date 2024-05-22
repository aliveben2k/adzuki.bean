#!/usr/bin/perl
# perl this_code.pl INPUT OUTPUT
# Only accepts bi-allelic variants with alleles 0 and 1

chomp(@ARGV);
my $input_file = $ARGV[0];
my $output_file = $ARGV[1];
my $pop_file = $ARGV[2];
my $number_of_title_or_comment_lines;

if ($input_file =~ /gz$/){
	$number_of_title_or_comment_lines = `gzip -cd $input_file \| grep -c "^#"`;	# How many lines start with #? Also the last line is the title line
}
else {
	$number_of_title_or_comment_lines = `grep -c "^#" $input_file`;
}
my $how_many_fields_before_data = 9;	# For each line, how many fields (with other informations) are there before actually goes into the genotype data

if ($input_file =~ /gz$/){
	open (INPUT, "-|", "gzip -cd $input_file") || die "\[$time\]\: Cannot open $input: $!\n";
}
else {
	open (INPUT,"<$input_file") || die "Open input file error~!\n";
}
# First ? lines are info. Ignore.
my $title_line_raw;
for (my $i = 1; $i <= $number_of_title_or_comment_lines; $i ++) {
	my $temp = <INPUT>;	# Just want to call INPUT here to let those lines being skipped
	if ($i == $number_of_title_or_comment_lines) {	# Last line is title line
		$title_line_raw = $temp;
	}
}

#if population file is provided, generate output file by population frequency
%pop_freq; my @uni_pop;
if ($pop_file ne ""){
	open(POP, "<$pop_file") || die "Cannot open $pop_file: $!\n";
	my @pop_list = <POP>;
	chomp(@pop_list);
	if ($pop_list[0] =~ /^ind\B|^id|^taxa|^samp\B/i){
		shift(@pop_list);
	}
	my @ids; my @pop_names;
	foreach (@pop_list){
		my @eles = split(/\t/, $_);
		push(@ids, $eles[0]);
		push(@pop_names, $eles[1]);
	}
	@uni_pop = do { my %seen; grep { !$seen{$_}++ } @pop_names };
	@uni_pop = sort @uni_pop;
	foreach (@uni_pop){ #define samples in the populations
		foreach my $i (0..$#pop_names){
			if ($pop_names[$i] eq $_){
				if (@{$pop_freq{$_}}){
					push(@{$pop_freq{$_}}, $ids[$i]);
				}
				else {
					@{$pop_freq{$_}} = $ids[$i];
				}
			}
		}
	}
}

# For @title_line and every following line, element 0 is chromosome, 1 is position, 3 is ref allele, 4 is alt allele, starting from element $how_many_fields_before_data onwords are 1135 data
$title_line_raw =~ s/^.//;	# Cut off the first character #
my @title_line = split(/\t/,$title_line_raw);
my @accession = @title_line;
for (my $i = 0; $i <= ($how_many_fields_before_data - 1); $i ++) {	# Get rid of non-accession fields
	shift(@accession);
}
chomp(@accession);

if ($pop_file ne ""){
	foreach my $key (sort keys %pop_freq){
		my @current_pop_ids = @{$pop_freq{$key}};
		foreach my $k (0..$#current_pop_ids){
			my $exist = 0;
			foreach my $j (0..$#accession){
				if ($current_pop_ids[$k] eq $accession[$j]){
					$exist = 1;
					$pop_freq{$key}[$k] = $j;
					last;
				}
			}
			if ($exist == 0){
				$pop_freq{$key}[$k] = -1;
			}
		}
		@{$pop_freq{$key}} = grep { $_ >= 0 } @{$pop_freq{$key}};
	}
	#debug:
	#foreach my $key (sort keys %pop_freq){
	#	print "$key:\n";
	#	my @current_pop_ids = @{$pop_freq{$key}};
	#	print join(" ", @current_pop_ids), "\n";
	#}
}

unless ($output_file =~ /gz$/){
	$output_file .= '.gz';
}
my $cnt = 1;
if ($output_file =~ /\.txt/){
	$output_file =~ s/\.txt/.$cnt.tmp.txt/;
}
else {
	$output_file =~ s/\.gz$/.$cnt.tmp.txt.gz/;
}
open (OUTPUT, "|-", "gzip \> $output_file") || die "Cannot write $output_file: $!\n";
print OUTPUT "SNP\tchr\tpos\tref\talt\t";
if ($pop_file ne ""){
	print OUTPUT join("\t", @uni_pop), "\n";
}
else {
	print OUTPUT join("\t", @accession), "\n";
}

my $line_cnt = 0;
while (<INPUT>) {
	# The output is for R to read, so I don't have to bother transforming things into nucleotide
	chomp $_;
	my @this_line = split(/\t/,$_);
	my $this_chr = $this_line[0];
	my $pos = $this_line[1];
	my $ref = $this_line[3];
	my $alt = $this_line[4];
	my $header = "$this_chr-$pos\t$this_chr\t$pos\t$ref\t$alt";

	my @genotypes = @this_line;
	for (my $i = 0; $i <= ($how_many_fields_before_data - 1); $i ++) {	# Get rid of non-accession fields
		shift(@genotypes);
	}
	my @geno_sum;
	for (my $i = 0; $i <= $#genotypes; $i ++) {
		my $allele_one = substr($genotypes[$i],0,1);	# The first character is the 1st allele num of each accession
		my $allele_two = substr($genotypes[$i],2,1);	# The third character is the 2nd allele num of each accession
		if (($allele_one eq '.') || ($allele_two eq '.')) {
			push(@geno_sum, "NA");
		} elsif ($allele_one == $allele_two) {	# Homozygous
			push(@geno_sum, "$allele_one");
		} else {	# Heterozygous
			push(@geno_sum, "0.5");
		}
	}
	my @unique = do { my %seen; grep { !$seen{$_}++ } @geno_sum };
	@unique = do { grep { !/NA/ } @unique };
	if (scalar(@unique) > 1){
		print OUTPUT "$header\t";
		if ($pop_file ne ""){
			my @pop_freq;
			foreach my $key (sort keys %pop_freq){ #point to the population level
				my @current_pop_index = @{$pop_freq{$key}};
				my $current_pop_freq = 0;
				my $sample_cnt = 0;
				foreach my $idx (@current_pop_index){ #get the population sample value for calculation
					if ($geno_sum[$idx] !~ /[^0-9\.]/){
						$current_pop_freq += $geno_sum[$idx];
						$sample_cnt++;
					}
				}
				if ($sample_cnt > 0){
					$current_pop_freq = sprintf("%.3f", $current_pop_freq/$sample_cnt);
				}
				else {
					$current_pop_freq = "NA";
				}
				push(@pop_freq, $current_pop_freq);
			}
			print OUTPUT join("\t", @pop_freq), "\n";
		}
		else {
			print OUTPUT join("\t", @geno_sum), "\n";
		}
	}
	$line_cnt++;
	if ($line_cnt == 10002){
		close(OUTPUT);
		$cnt++;
		my $pre_cnt = $cnt - 1;
		if ($output_file =~ /\.tmp.txt/){
			$output_file =~ s/\.$pre_cnt\.tmp.txt/.$cnt.tmp.txt/;
		}
		open (OUTPUT, "|-", "gzip \> $output_file") || die "Cannot write $output_file: $!\n";
		print OUTPUT "SNP\tchr\tpos\tref\talt\t";
		if ($pop_file ne ""){
			print OUTPUT join("\t", @uni_pop), "\n";
		}
		else {
			print OUTPUT join("\t", @accession), "\n";
		}
		$line_cnt = 0;
	}
}
close (INPUT);
close (OUTPUT);
