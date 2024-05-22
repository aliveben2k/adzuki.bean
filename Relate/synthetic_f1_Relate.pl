#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
my $time = scalar localtime();

chomp(@ARGV);
my $input = $ARGV[0];
my $syn_list = $ARGV[1];
my $missing;
if ($ARGV[2] eq "--no_missing"){
	$missing = 0;
}
else {
	$missing = 1;
}

open(SYNLIST, "<$syn_list") || die BOLD "Cannot open $syn_list: $!", RESET, "\n";
my @lists = <SYNLIST>;
chomp(@lists);
close(SYNLIST);

my $out = $input;
if ($out =~ /\.vcf$/){
	$out =~ s/\.vcf$/.syn\.vcf/;
}
else {
	$out =~ s/\.vcf\.gz$/\.syn\.vcf\.gz/;
}

my $out2 = $syn_list;
$out2 =~ s/.list$//;
$out2 =~ s/.txt$//;
$out2 .= ".poplabels";

if ($input =~ /\.gz$/){
	open(INPUT, "-|", "gzip -dc $input") || die BOLD "Cannot open $input: $!", RESET, "\n";
}
else {
	open (INPUT, "<$input") || die BOLD "Cannot open $input: $!", RESET, "\n";
	$out = $out."\.gz";
}
open (OUTPUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
my $exist;
unless (-e $out2){
	open (POP, ">$out2") || die BOLD "Cannot write $out2: $!", RESET, "\n";
	print POP "ID POP GROUP SEX\n";
	$exist = 0;
}	
my @sample_line; my @samples; my $pop_sample = 0; my @pos; my @name;
while (my $line = <INPUT>){
	my @out_eles;
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
					if ($_ =~ /$sample_line[$x]/){
						push(@pos, $x);
					}
				}
				my @list_tmps = split(/\t/, $_);
				$list_tmp .= "$list_tmps[2]\t";
				unless ($exist){
					print POP "$list_tmps[2] $list_tmps[3] $list_tmps[3] NA\n";
				}
			}
			unless ($exist){
				close(POP);
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
		for (my $j=0; $j<=7; $j++){	
				push(@out_eles, "$line_vec[$j]");
		}
		push(@out_eles, "GT"); 
		for (my $k=0; $k<=$#pos; $k+=2){
			my $fst_pos = int($pos[$k]);
			my $snd_pos = int($pos[$k+1]);
			my @fst_sample = split(/\/|\|/, $line_vec[$fst_pos]);
			my @snd_sample = split(/\/|\|/, $line_vec[$snd_pos]);
			push(@out_eles, "$fst_sample[0]\|$snd_sample[0]");
		}
		my $dot = 0;
		if ($missing == 0){
			foreach my $l (9..$#out_eles){
				if ($out_eles[$l] =~ /\./){
					$dot = 1;
					last;
				}
			}
		}
		if ($dot == 0){
			print OUTPUT join("\t", @out_eles), "\n";
		}
	}
}
close(INPUT);
close(OUTPUT);

