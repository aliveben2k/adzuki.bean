#!/usr/local/bin/perl
#usage: perl split_gmap.pl GENETIC_MAP_FILE OUTPUT_FOLDER [PREFIX]

chomp(@ARGV);

my $file = $ARGV[0];
my $out = $ARGV[1];
my $pre = $ARGV[2];
unless ($file){
    die "Cannot find the input file.\n";
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
open(GMAP, "<$file") || die "Cannot open $file: $!\n";
my @content = <GMAP>;
chomp(@content);
close(GMAP);

my $chr; my $cnt = 1;
foreach (@content){
	$_ =~ s/[\x0A\x0D]//g;
	@line_eles = split(/\t|\s+/, $_);
	if ($pre){
		if ($line_eles[0] !~ /^$pre/){
			next;
		}
	}
	$chr = $line_eles[0] unless $chr;
	if ($cnt == 1){
		open(OUT, ">$out\/$cnt.gmap.map") || die "Cannot write $out\/$cnt.gmap.map: $!\n";
		print OUT "pos COMBINED_rate Genetic_Map\n";
	}
	if ($line_eles[0] ne $chr){
		close(OUT);
		$cnt++;
		open(OUT, ">$out\/$cnt.gmap.map") || die "Cannot write $out\/$cnt.gmap.map: $!\n";
		print OUT "pos COMBINED_rate Genetic_Map\n";
		$chr = $line_eles[0];
	}
	shift(@line_eles);
	print OUT join(' ', @line_eles), "\n";
}
close(OUT);

	
