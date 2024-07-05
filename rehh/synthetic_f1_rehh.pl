#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
my $time = scalar localtime();

chomp(@ARGV);
my $input = $ARGV[0];
my $syn_list = $ARGV[1];

#my $in_path;
my $in_name;
if ($input =~ /\//){
my @in_paths = split(/\//, $input);
	$in_name = pop(@in_paths);
	#$in_path = join("\/", @in_paths);
}
else {
	#$in_path = ".";
	$in_name = $input;
}

my $out_path;
if ($syn_list =~ /\//){
	my @out_paths = split(/\//, $syn_list);
	pop(@out_paths);
	$out_path = join("\/", @out_paths);
}
else {
	$out_path = ".";
}


open(SYNLIST, "<$syn_list") || die BOLD "Cannot open $syn_list: $!", RESET, "\n";
my @lists = <SYNLIST>;
chomp(@lists);
close(SYNLIST);

my $out = "$out_path\/$in_name";
$out =~ s/\.vcf\.gz$|\.vcf$/.syn\.vcf.gz/;


if ($input =~ /\.gz$/){
	open(INPUT, "-|", "gzip -dc $input") || die BOLD "Cannot open $input: $!", RESET, "\n";
}
else {
	open (INPUT, "<$input") || die BOLD "Cannot open $input: $!", RESET, "\n";
	$out = $out."\.gz";
}
print "$out\n";
open (OUTPUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";

my @sample_line; my @samples; my $pop_sample = 0; my @pos; my @name;
while (my $line = <INPUT>){
	$fixed = 0;
	chomp $line;
	if ($line =~ /^#/){
		if ($line =~ /^#CHROM/){
			@sample_line = split(/\t/, $line);
			for (my $i=0; $i<=8; $i++){
				print OUTPUT "$sample_line[$i]\t";
			}
			my $list_tmp;
			foreach (@lists){
				for (my $x=0; $x<=$#sample_line; $x++){
					if ($_ =~ /\b$sample_line[$x]\b/){
						push(@pos, $x);
					}
				}
				my @list_tmps = split(/\t/, $_);
				$list_tmp .= "$list_tmps[2]\t";
			}
			$list_tmp =~ s/\t$//;
			print OUTPUT "$list_tmp\n";
		}
		else {
			print OUTPUT "$line\n";
		}
	}
	else {
		my @line_vec = split('\t', $line);
		for (my $j=0; $j<=8; $j++){
				print OUTPUT "$line_vec[$j]\t";				
		}
		for (my $k=0; $k<=$#pos; $k+=2){
			my $fst_pos = int($pos[$k]);
			my $snd_pos = int($pos[$k+1]);
			my @fst_sample = split(/\/|\|/, $line_vec[$fst_pos]);
			my @snd_sample = split(/\/|\|/, $line_vec[$snd_pos]);
			if ($k < $#pos-1){
				print OUTPUT "$fst_sample[0]\/$snd_sample[0]\t";
                #print "$fst_sample[0]\/$snd_sample[0]\t";
			}
			else {
				print OUTPUT "$fst_sample[0]\/$snd_sample[0]\n";
                #print "$fst_sample[0]\/$snd_sample[0]\n";
			}
		}
	}
}
close(INPUT);
close(OUTPUT);

