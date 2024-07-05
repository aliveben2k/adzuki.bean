#!/usr/bin/perl
use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;	
}
my $path; my @files; my $vcf;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "-d"){
		if (-d $ARGV[$i+1]){
			$path = $ARGV[$i+1];
			if ($path =~ /\/$/){
				$path =~ s/\/$//;
			}
			@files = <$path\/*.vcf.gz>;
			@files = grep{!/\.m\./} @files;
		}
		elsif (-e $ARGV[$i+1]){
			@files = $ARGV[$i+1];
			if ($ARGV[$i+1] =~ /\//){
				@tmp = split(/\//, $ARGV[$i+1]);
				pop(@tmp);
				$path = join("\/", @tmp);
			}
			else {
				$path = ".";
			}
		}		
	}
	if ($ARGV[$i] eq "-v"){
		if (-e "$ARGV[$i+1]"){
			$vcf = $ARGV[$i+1];
		}
	}
}
die "Cannot find -v or -d\n" unless ($path && $vcf);

my $ID_info;
open(INPUT, "-|", "gzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
while (my $line = <INPUT>){
	chomp($line);
	if ($line =~ /^\#/){
		if ($line =~ /^\#\#contig\=\<ID\=/){
			$ID_info .= "$line\n";
#			push(@ID_info, $line);
		}
	}
	else {
		last;
	}
}
close(INPUT);

foreach (@files){
	open(INPUT, "-|", "gzip -dc $_") || die BOLD "Cannot open $_: $!", RESET, "\n";
	my $out = $_;
	$out =~ s/\.vcf\.gz$/\.m\.vcf\.gz/;
	open (OUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
	while (my $line = <INPUT>){
		chomp($line);
		if ($line =~ /^\#\#/){
			unless ($line =~ /^\#\#contig\=\<ID\=/){
				print OUT "$line\n";
			}
		}
		elsif ($line =~ /^\#CHROM/){
			print OUT "$ID_info";
			print OUT "$line\n";
		}
		else {
			print OUT "$line\n";
		}
	}
	close(INPUT);
	close(OUT);
}
sub usage{
	print "Usage: perl add_contig_info.pl -d VCF_PATH -v ORINGINAL_VCF\n";
}