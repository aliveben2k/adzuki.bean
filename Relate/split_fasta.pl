#!/usr/local/bin/perl
#usage: perl split_fasta.pl FASTA_FILE OUTPUT_FOLDER [PREFIX]

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
$cnt = 1; $exist = 1;
while (my $line = <IN>){
    chomp($line);
    if ($line =~ /\>/){
        if ($cnt != 1){
            close(OUT);
        }
        if ($pre){
        	if ($line !~ /^\>$pre/i){
        		$exist = 0;
        		next;
        	}
        	else {
        		$exist = 1;
        	}
        }
        open(OUT, ">$out\/$cnt.fasta") || die "Cannot write $out\/$cnt.fasta: $!\n";
        print OUT "\>$cnt\n";
        $cnt++;
    }
    else {
    	if ($exist == 0){
    		next;
    	}
        print OUT "$line\n";
    }
}
close(OUT);
close(IN);
