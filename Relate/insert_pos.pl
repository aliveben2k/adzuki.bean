#!/usr/bin/perl

chomp(@ARGV);
if ($#ARGV == -1){
	print "Usage: perl insert_pos.pl [-i insertion_file -d deletion_file] -r reference_file\n";
	exit;	
}

my $input; my $ref; my $del;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-i"){
		$input = $ARGV[$i+1];
		unless (-e $input){
			die "Cannot find the input file.\n";
		}
	}
	if ($ARGV[$i] eq "\-d"){
		$del = $ARGV[$i+1];
		unless (-e $del){
			die "Cannot find the input file.\n";
		}
	}
	if ($ARGV[$i] eq "\-r"){
		$ref = $ARGV[$i+1];
		unless (-e $ref){
			die "Cannot find the reference file.\n";
		}
	}
}

unless (($input || $del) && $ref){
	die "Not enough arguments provided.\n";
}

my $hap_fmt = 0;
if ($ref =~ /vcf\.gz$/) {
	open(REF, "-|", "gzip -dc $ref") || die "Cannot open $ref: $!\n";
}
elsif ($ref =~ /vcf$/){
	open(REF, "<$ref") || die "Cannot open $ref: $!\n";
}
elsif ($ref =~ /haps\.gz$/) {
	open(REF, "-|", "gzip -dc $ref") || die "Cannot open $ref: $!\n";
	$hap_fmt = 1;
}
elsif ($ref =~ /haps$/){
	open(REF, "<$ref") || die "Cannot open $ref: $!\n";
	$hap_fmt = 1;
}
print "debug: hap_fmt $hap_fmt\n";
my @chr; my @pos; my @inputs;
if ($input){
	if ($input =~ /gz$/){
		open(INPUT, "-|", "gzip -dc $input") || die "Cannot open $input: $!\n";
	}
	else {
		open(INPUT, "<$input") || die "Cannot open $input: $!\n";
	}
	@inputs = <INPUT>;
	chomp(@inputs);
	close(INPUT);
	foreach (@inputs){
		my @tmp = split(/\s+|\t+/, $_);
		push(@chr, $tmp[0]);
		if ($hap_fmt == 0){
			push(@pos, $tmp[1]);
		}
		else {
			push(@pos, $tmp[2]);
		}
	}
}
my @chr_del; my @pos_del_start; my @pos_del_end;
if ($del){
	open(INPUT, "<$del") || die "Cannot open $del: $!\n";
	my @dels = <INPUT>;
	chomp(@dels);
	close(INPUT);
	foreach (@dels){
		my @tmp = split(/\:|\-/, $_);
		push(@chr_del, $tmp[0]);
		push(@pos_del_start, $tmp[1]);
		push(@pos_del_end, $tmp[2]);
	}
	print "debug: @chr_del\n@pos_del_start\n@pos_del_end\n";
}

my $out;
$out = $ref;
if ($out =~ /vcf\.gz$/) {
	$out =~ s/vcf\.gz$/m.vcf\.gz/;
	open(OUT, "|-", "bgzip -c \> $out") || die "Cannot write $out: $!\n";
}
elsif ($out =~ /vcf$/){
	$out =~ s/vcf$/m.vcf/;
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
}
elsif ($out =~ /haps\.gz$/) {
	$out =~ s/haps\.gz$/m.haps\.gz/;
	open(OUT, "|-", "gzip -c \> $out") || die "Cannot write $out: $!\n";
}
elsif ($out =~ /haps$/){
	$out =~ s/haps$/m.haps/;
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
}

my $last_chr; my $last_pos = 0;
while (my $line = <REF>){
	chomp($line);
	if ($line =~ /^\#/){
		print OUT "$line\n";
		next;
	}
	else {
		my $replace = 0;
		my @tmp = split(/\s+|\t+/, $line); #reference
		$last_chr = $tmp[0] unless ($last_chr);
		if ($last_chr ne $tmp[0]){
			$last_chr = $tmp[0];
			$last_pos = 0;
		}
		if ($del){
			foreach my $j (0..$#chr_del){
				if ($tmp[0] eq $chr_del[$j]){
					if ($hap_fmt == 0){
						 if ($tmp[1] > $pos_del_start[$j] && $tmp[1] < $pos_del_end[$j]){
							#print "debug: vcf $tmp[1]\n";
							$replace = 1;
						 }
					}
					else {
						 if ($tmp[2] > $pos_del_start[$j] && $tmp[2] < $pos_del_end[$j]){
							#print "debug: hap $tmp[2]\n";
							$replace = 1;
						 }
					}
				}
			}
		}
		if ($input){
			foreach my $i (0..$#chr){
				#print "debug: $tmp[0] $chr[$i]\n"; 
				if ($tmp[0] eq $chr[$i]){
					if ($hap_fmt == 0){
						#print "debug: $tmp[1] $pos[$i] $last_pos\n"; 
						if ($tmp[1] < $pos[$i]){
							$last_pos = $tmp[1];
							next;
						}
						elsif ($tmp[1] > $pos[$i] && $pos[$i] > $last_pos) {
							print OUT "$inputs[$i]\n";
							#print "debug: insert 1\n";
							$last_pos = $pos[$i];
						}
						elsif ($tmp[1] == $pos[$i]) {
							print OUT "$inputs[$i]\n";
							#print "debug: insert 2\n";
							$last_pos = $pos[$i];
							$replace = 1;
						}
					}
					else {
						if ($tmp[2] < $pos[$i]){
							$last_pos = $tmp[2];
							next;
						}
						elsif ($tmp[2] > $pos[$i] && $pos[$i] > $last_pos) {
							print OUT "$inputs[$i]\n";
							#print "debug: insert 3\n";
							$last_pos = $pos[$i];
						}
						elsif ($tmp[2] == $pos[$i]) {
							print OUT "$inputs[$i]\n";
							#print "debug: insert 4\n";
							#print "debug: $inputs[$i]\n";
							$last_pos = $pos[$i];
							$replace = 1;
						}
					}
				}
			}
		}
		if ($replace == 0){
			print OUT "$line\n";
		}
	}
}
close(OUT);
close(REF);

