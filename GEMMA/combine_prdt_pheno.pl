#!/usr/bin/perl

chomp(@ARGV);
if ($#ARGV == -1){
    print "Usage: perl combine_prdt_pheno.pl *.pheno_file prefix_of_*.X.prdt.txt\n";
    print "* is prefix in your file name. X is phenotype index in the file name.\n";
	die "This script is a part of gemma2.pl.\n";
}

my $ori = $ARGV[0];
my $prefix = $ARGV[1];
open (ORI, "<$ori") || die "Cannot open $ori: $!\n";
my @ori_lines = <ORI>;
close(ORI);
chomp(@ori_lines);
my $prdt_out = $ori;
$prdt_out =~ s/pheno$/prdt.pheno/;
shift(@ARGV);
open (OUT, ">$prdt_out") || die "Cannot write $prdt_out: $!\n";
print "Combining *.prdt.txt files\n";
foreach my $i (0..$#ori_lines){ # $i is sample (row) index
	@line_eles = split(/\t+|\s+|\,/, $ori_lines[$i]);
	foreach my $j (0..$#line_eles){ # $j is phenotype (file & column) index
		if ($line_eles[$j] eq "NA"){ #column index start from 0
            my $num = $j+1; #file index start from 1
            if (-e "$prefix.$num.prdt.txt"){
                open(PHE, "<$prefix.$num.prdt.txt") || die "Cannot find $prefix.$num.prdt.txt: $!\n";
                @prdt_vals = <PHE>;
                close(PHE);
                chomp(@prdt_vals);
#                print "$line_eles[$j]\t$prdt_vals[$i]\n";
                $line_eles[$j] = $prdt_vals[$i]; #replace original position $i(row), $j(column) with prdt_val $i(row)
			}
		}
	}
	$ori_lines[$i] = join("\t", @line_eles);
	print OUT "$ori_lines[$i]\n";
}
close(OUT);
print "Finished.\n";
