#!/usr/bin/perl
#usage: perl split_fasta.pl MASKED_FASTA_FILE OUTPUT_FOLDER [PREFIX]

use Term::ANSIColor qw(:constants);
chomp(@ARGV);

my $file = $ARGV[0];
my $out = $ARGV[1];
my $pre = $ARGV[2];
my $err;
unless ($file){
    die "Cannot find the input fasta file.\n";
}
unless ($out){
    die "Output folder is not defined.\n";
}
if (-d $out){
	$out = "$out";
}
else{
    $err = `mkdir $out 2>&1`;
    if ($err){
        die "$!\n";
    }
    $out = "$out";
}
if ($out =~ /\/$/){
    $out =~ s/\/$//;
}

open (IN, "<$file") || die ("Cannot open $file\n");
my $cnt = 1; my $exist = 1;
while (my $line = <IN>){
    chomp($line);
    if ($line =~ /\>/){
        if ($cnt != 1){
            close(OUT);
        }
        if ($pre){
        	if ($line !~ /^\>$pre/){
        		$exist = 0;
        		next;
        	}
        	else {
        		$exist = 1;
        	}
        }
        open(OUT, ">$out\/$cnt.masked.fasta") || die "Cannot write $out\/$cnt.masked.fasta: $!\n";
        print OUT "\>$cnt\n";
        $cnt++;
    }
    else {
    	if ($exist == 0){
    		next;
    	}
        $line =~ s/[^N]/P/ig;
        print OUT "$line\n";
    }
}
close(OUT);
close(IN);