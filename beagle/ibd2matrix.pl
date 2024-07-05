#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);
#my $r_env = '/home/hpc/crlee/miniconda3/envs/R-3.5/bin';
#my $Relate = "\$Relate";

chomp(@ARGV);
print "The script is written by Ben Chien. Dec. 2021.\n";
print "Input command line:\n";
print "perl ibd2matrix\.pl @ARGV\n";

if ($#ARGV == -1){
	&usage;
	exit;
}

my @ibds; my $path; my $ibd; my @samples;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-i"){
		if (-d $ARGV[$i+1]){
			if ($ARGV[$i+1] =~ /\/$/){
				$ARGV[$i+1] =~ s/\/$//;
			}
			$path = $ARGV[$i+1];
			if (-e "$path\/all.ibd.gz"){
				system("rm $path\/all.ibd.gz");
			}
			if (-e "$path\/all.ibd.gz.matrix"){
				system("rm $path\/all.ibd.gz.matrix");
			}
			@ibds = <$ARGV[$i+1]\/*.ibd.gz>;
			unless (@ibds){
				die "Cannot find any ibd file.\n";
			}
		}
		elsif (-e $ARGV[$i+1]){
			$ibd = $ARGV[$i+1];
		}
		else {
			die "Cannot find any ibd file.\n";
		}
	}
	if ($ARGV[$i] eq "\-vcf"){
		if (-e $ARGV[$i+1]){
			@samples = &sample_name($ARGV[$i+1]);
		}
		else {
			die "Cannot find the vcf file.\n";
		}
	}
	if ($ARGV[$i] eq "\-list"){
		if (-e $ARGV[$i+1]){
			open(LIST, "<$ARGV[$i+1]") || die "Cannot find the list.\n";
			@samples = <LIST>;
			chomp(@samples);
			close(LIST);
		}
	}
}
if (@ibds){
	system("cat @ibds \> $path\/all.ibd.gz");
	$ibd = "$path\/all.ibd.gz";
}
unless (@samples){
	die "Cannot find the sample list. Please provide correct vcf file.\n";
}

open(INPUT, "-|", "gzip -dc $ibd") || die BOLD "Cannot open $ibd: $!", RESET, "\n";
my $out = "$ibd.matrix";
my %main_hash;
foreach $i (0..$#samples){
	my @tmp;
	foreach (@samples){
		push(@tmp, 0);
	}
	@{$main_hash{$i}} = @tmp;
}
while (my $line=<INPUT>){
	chomp($line);
	$line =~ s/[\x0A\x0D]//g;
	my @eles = split(/\t+|\s+/, $line);
	for my $j (0..$#samples){
		$samples[$j] =~ s/[\x0A\x0D]//g;
		if ($eles[0] eq $samples[$j]){
			for my $k (0..$#samples){
				$samples[$k] =~ s/[\x0A\x0D]//g;
				if ($eles[2] eq $samples[$k]){
					${$main_hash{$j}}[$k] += $eles[-1];
					${$main_hash{$k}}[$j] += $eles[-1];
				}
			}
		}
	}
}
for my $i (0..$#samples){
	for my $j (0..$#samples){
		if ($i == $j){
			${$main_hash{$i}}[$j] = 'NA';
		}
	}
}

open(OUT, ">$out") || die "Cannot write $out: $!.\n";
print OUT "\t", join("\t", @samples), "\n";

my @finals;
foreach my $key (keys %main_hash){
	my @values = @{$main_hash{$key}};
	my $final = "$samples[$key]\t";
	$final .= join("\t", @values);
	push(@finals, $final);
}
foreach (@samples){
	foreach my $final (@finals){
		if ($final =~ /^$_/){
			print OUT "$final\n";
		}
	}
}

sub usage {
	print BOLD "Usage: perl ibd2matrix.pl -i ibd_file\(s\) -vcf\|-list VCF_FILE\|LIST_FILE\n", RESET;
	return;
}

sub sample_name {
	my $file = shift;
	my @content; my @line;
	if (-e $file){}
	else {
		return 2;
	}
	if ($file =~ /\.vcf\.gz/){
		@content = `gzip \-cd $file \| head \-n 10000`;
	}
	elsif ($file =~ /\.vcf$/){
		@content = `head -n 10000 $file`;
	}
	else {
		return 2;
	}
	foreach (@content){
		if ($_ =~ /\#CHROM/){
			@line = split(/\t/, $_);
			for (my $i=9; $i<=$#line; $i++){
                push(@samples, $line[$i]);
			}
			last;
		}
	}
	chomp(@samples);
	return (@samples);
} #get sample name