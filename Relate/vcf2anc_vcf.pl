#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
print "The script is written by Ben Chien. Jan. 2023.\n";
print "Usage: perl vcf2anc_vcf.pl -vcf TARGET_VCF -aid ANCESTRAL_ID\|A_LIST_FILE [-o OUTPUT_FILE_NAME] [-list A_LIST_FILE] [-bi] [-rchr RENAME_FILE] [-keep] [-hap] [-nm] [-sf] [-thap]\n";
print "-aid: an ancestral id. Multiple samples could be seperated by comma, or listed in a file \(one sample per line\).\n";
print "-o: output file name without extension.\n";
print "-list: only keep the listed samples.\n";
print "-bi: only keep bi-allele SNPs. Default: false\n";
print "-rchr: rename CHR column. Default: false\n";
print "-keep: also keep the ancestral genotype. Default: false\n";
print "-hap: only use the first allele \(haplotype\). Default: false\n";
print "-nm: no missing, the missing haploid will be imputed as the other available haploid or the ancient one.\n";
print "-sf: outputting Shapeit format instead of vcf format.\n";
print "-thap: outputting rehh\:thap format instead of vcf format.\n";
print "The output file is the same as the input vcf, and the name is PREFIX.anc.vcf.gz\n";
print "This is a local script, not a server script.\n\n";
print "Input command line:\n";
print "perl vcf2anc_vcf\.pl @ARGV\n\n";

my $vcf; my @aids; my $hap = 0; my $kaid = 0; my $list; my @lists; my $bi = 0; my @rm_lists; my $nm = 0; my $sf = 0; my $out; my $thap = 0;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-vcf"){
		if (-e $ARGV[$i+1]){
            $vcf = $ARGV[$i+1];
            if ($ARGV[$i+1] =~ /\//){
                my @tmps = split(/\//, $ARGV[$i+1]);
                pop(@tmps);
                $path = join("\/", @tmps);
            }
        }
		else {
            die "Cannot find the file\(s\).\n";
		}
	}
	if ($ARGV[$i] eq "\-aid"){
		if (-e $ARGV[$i+1]){
			open(ALIST, "<$ARGV[$i+1]") || die BOLD "Cannot open aid list, $ARGV[$i+1]: $!", RESET, "\n";
			@aids = <ALIST>;
			chomp(@aids);
			close(ALIST);
            foreach my $j (0..$#aids){
            	my @sp_lines = split(/\t+|\s+/, $aids[$j]);
            	$aids[$j] = $sp_lines[0];
            }	  
		}
		elsif ($ARGV[$i+1] =~ /\,/){
			@aids = split(/\,/, $ARGV[$i+1]);
		}
		else {
			push(@aids, $ARGV[$i+1]);
		}
	}
	if ($ARGV[$i] eq "\-keep"){ #keep ancestral samples
        $kaid = 1;
	}
	if ($ARGV[$i] eq "\-hap"){
        $hap = 1;
	}
	if ($ARGV[$i] eq "\-bi"){
        $bi = 1;
	}
	if ($ARGV[$i] eq "\-nm"){ #no missing
        $nm = 1;
	}
	if ($ARGV[$i] eq "\-rchr"){
		if (-e $ARGV[$i+1]){
			open(RE, "<$ARGV[$i+1]") || die BOLD "Cannot open the rename list, $ARGV[$i+1]: $!", RESET, "\n";
			@rm_lists = <RE>;
			chomp(@rm_lists);
			close(RE);
        }
        else {
        	die BOLD "Cannot find the rename list.", RESET, "\n";
        }
	}
	if ($ARGV[$i] eq "\-list"){
        if (-e $ARGV[$i+1]){
            open (LIST, "<$ARGV[$i+1]") || die BOLD "Cannot open sample list, $ARGV[$i+1]: $!", RESET, "\n";
            @lists = <LIST>;
            chomp(@lists);
            if ($ARGV[$i+1] =~ /poplabels/){
            	shift(@lists);
            }
            close(LIST);
            foreach my $j (0..$#lists){
            	my @sp_lines = split(/\t+|\s+/, $lists[$j]);
            	$lists[$j] = $sp_lines[0];
            }
            $list = $ARGV[$i+1];
        }
	}
	if ($ARGV[$i] eq "\-sf"){ #shapeit format for Relate
        $sf = 1;
	}
	if ($ARGV[$i] eq "\-thap"){ #thap format: for rehh
        $thap = 1;
	}	
	if ($ARGV[$i] eq "\-o"){
        $out = $ARGV[$i+1];
	}
}

unless (-e $vcf){
	die "Cannot find target vcf: $vcf.\n";
}

#processing vcf
if ($vcf =~ /\.gz$/){
    open(INPUT, "-|", "gzip -dc $vcf") || die BOLD "Cannot open vcf $vcf: $!", RESET, "\n";
}
else {
    open (INPUT, "<$vcf") || die BOLD "Cannot open vcf $vcf: $!", RESET, "\n";
}
print "Start processing vcf...\n";

unless ($out){
	$out = $vcf;
	$out =~ s/.gz$//;
	if ($list){
		if ($list =~ /\//){
			my @tmp = split(/\//, $list);
			$list = $tmp[-1];
		}
		$list =~ s/\.txt$|\.list$//;
		$out =~ s/\.vcf$/.$list.vcf/;
	}
	$out =~ s/\.vcf$/.anc.vcf.gz/;
	if ($sf == 1){
		$out =~ s/vcf\.gz$/haps/;
	}
	elsif ($thap == 1){
		$out =~ s/vcf\.gz$/thap/;
	}
}
else {
	if ($sf == 0){ #vcf
		$out .= ".vcf.gz";
	}
	elsif ($thap == 1){
		$out .= ".thap";
	}
	else { #haps format
		$out .= ".haps";
	}
}
if ($sf == 1){
	my $out2 = $out;
	$out2 =~ s/\.haps$/.sample/;
	open(OUT, ">$out") || die BOLD "Cannot write $out: $!", RESET, "\n";
	open(SAMPLE, ">$out2") || die BOLD "Cannot write $out2: $!", RESET, "\n";
}
elsif ($thap == 1){
	my $out2 = $out;
	$out2 =~ s/\.thap$/.map/;
	open(OUT, ">$out") || die BOLD "Cannot write $out: $!", RESET, "\n";
	open(SAMPLE, ">$out2") || die BOLD "Cannot write $out2: $!", RESET, "\n";
}
else {
	open(OUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
}
#print "debug\n";
my $cnt = 1; my $ck_aid = 0; my @samples; my @index; my $indexes;
while (my $line = <INPUT>){
	chomp($line);
	if ($line =~ /^\#/){
		if ($line =~ /^\#CHROM/){
			my @ids = split(/\t/, $line);
			foreach my $i (0..$#ids){
				if ($i < 9){
					if ($sf == 0 && $thap == 0){
						print OUT "$ids[$i]\t";
					}
				}
				if ($i >= 9){
					foreach my $j (0..$#aids){
						if ($ids[$i] eq $aids[$j]){
							$aids[$j] = $i;
							if ($ck_aid == 0){
								$ck_aid = 1;
							}
						}
					}
				}
			}
			my @l_samples;
			foreach my $k (9..$#ids){
				if ($kaid == 0){
					if (join(' ', @ids[@aids]) !~ /\b$ids[$k]\b/){
						if ($list){
							if (join(' ', @lists) =~ /\b$ids[$k]\b/){
								push(@index, $k);
								push(@l_samples, $ids[$k]);
							}
						}
						else {
							push(@index, $k);
							push(@l_samples, $ids[$k]);
						}
					}
                }
                else {
					if ($list){
						if (join(' ', @lists) =~ /\b$ids[$k]\b/){
							push(@index, $k);
							push(@l_samples, $ids[$k]);
						}
					}
					else {
						push(@index, $k);
						push(@l_samples, $ids[$k]);
					}
                }
			}
            if ($list){
                $indexes = join(" ", @index);
            }
            if ($sf == 0){ #output vcf
				print OUT join("\t", @l_samples), "\n";
			}
			else { #output haps and sample
				print SAMPLE "ID_1 ID_2 missing\n0 0 0\n";
				if ($hap == 0){
					foreach my $l (0..$#l_samples){
						print SAMPLE "$l_samples[$l] $l_samples[$l] 0\n";
					}
				}
				else {
					foreach my $l (0..$#l_samples){
						print SAMPLE "$l_samples[$l] NA 0\n";
					}
				}
				close(SAMPLE);
			}
		}
		elsif ($sf == 0) {
			print OUT "$line\n";
		}
        next;
	}
	else {
		if ($ck_aid == 0){
			die "Cannot find any ancestral ID in the vcf.\n";
		}
		my @eles = split(/\t/, $line);
		#get nucleotide from ref. sequence
		my $nucl = &get_ancestral_allele($line, \@aids);
		$nucl = uc($nucl);
        $eles[2] = "$eles[0]\_$eles[1]";
		#check if ancestral allele exists in REF or ALT column
		my $check_existance = $eles[3];
		$check_existance .= ",$eles[4]";
		if ($nucl =~ /\*/){
			if ($check_existance !~ /\*/){
            	next;
			}
			$nucl =~ s/\*/B/g;		
		}
		if ($check_existance =~ /\*/){
			$check_existance =~ s/\*/B/g;
		}
		if ($check_existance !~ /$nucl/){
            next;
		}
		$nucl =~ s/B/\*/g;
		$check_existance =~ s/B/\*/g;
		my $anc; my $deriv; my @derivs;
		my @sorted_nucls = split(/\,/, $check_existance); #all variants (REF + ALT)
		foreach my $k (0..$#sorted_nucls){
			if ($nucl eq $sorted_nucls[$k]){
				$anc = $k;
			}
			else {
				push(@derivs, $k);
			}
		}
		unshift(@derivs, $anc); # re-sort variants
		@sorted_nucls = @sorted_nucls[@derivs];
		if ($bi == 1 && $hap == 0){ #keep only bi-allele SNPs
			if (scalar(@derivs) > 2){
				next;
			}
		}
        my @out_eles; my $exist = 0;
		for my $j (9..$#eles){
			my @out_alleles;
			if ($kaid == 0){ #if don't keep ancestral IDs, skip the IDs
				if (join(' ', @aids) =~ /\b$j\b/){
					next;
				}
			}
			if ($list){ #if only keep listed IDs, skip IDs that are not matched
				if ($indexes !~ /\b$j\b/){
					next;
				}
			}
            my @gt = split(/\:/, $eles[$j]);
            my @alleles = split(/\/|\|/, $gt[0]);
            #print "debug: @alleles\n";
            if (scalar(@alleles) != 2 && $hap == 0){
            	if ($nm == 0){
                	@alleles = ('.', '.');
                	push(@out_eles, join("\|", @alleles));
                	#print "debug2: @alleles\n";
                }
                if ($nm == 1){
                	if (scalar(@alleles) == 1){
                		@alleles = ($alleles[0], $alleles[0]);
                	}
                	else {
                		@alleles = (0, 0);
                	}
                	push(@out_eles, join("\|", @alleles));
                	#print "debug2: @alleles\n";
                }
            }
            elsif (scalar(@alleles) == 1 && $hap == 1){
                if ($alleles[0] !~ /[^0-9]/){
                	foreach my $m (0..$#derivs){
                		if ($alleles[0] == $derivs[$m]){
                			$out_alleles[0] = $m;
                		}
                	}
                    #print "debug2: $out_alleles[0]\n";
                    push(@out_eles, $out_alleles[0]);
                }
                else {
                	if ($nm == 0){
                    	push(@out_eles, '.');
                    	#print "debug2: .\n";
                    }
                    else {
                    	push(@out_eles, 0);
                    	#print "debug2: 0\n";
                    }
                }
            }
            else {
            	foreach my $l (0..1){
            		if ($alleles[$l] !~ /[^0-9]/){
                		foreach my $m (0..$#derivs){
                			if ($alleles[$l] == $derivs[$m]){
                				push(@out_alleles, $m);
                			}              			
                		}              			
                	}
                	else {
                		if ($nm == 0){
                			push(@out_alleles, '.');
                		}
                		else {
                			push(@out_alleles, 0);
                		}
                	}
                }
                if ($hap == 0){
                	push(@out_eles, join("\|", @out_alleles));      	
                }
                else {
                	push(@out_eles, $out_alleles[0]);
                }
            }
		}
		foreach my $q (0..8){
			if ($q == 0 && @rm_lists){
				$eles[$q] = &rename_chr($eles[$q], \@rm_lists);
			}
			if ($q == 3){ #NUC
				$eles[$q] = $sorted_nucls[0];
				if ($thap == 1){
					$eles[$q] = 0;
				}
				shift(@sorted_nucls);
			}
			if ($q == 4){ #ALT
				if ($bi == 0){
					if ($thap == 1){
						foreach my $x (1..scalar(@sorted_nucls)){
							$sorted_nucls[$x-1] = $x;
						}
					}
					$eles[$q] = join(/\,/, @sorted_nucls);
				}
				else {
					if ($thap == 1){
						$eles[$q] = 1;
					}
				}
			}
			if ($q == 8){
				$eles[$q] = "GT";
			}
		}
		#check if the position is a bi-allele type
		if ($bi == 1){
			my @unique = do {my %seen; grep {!$seen{$_}++} @out_eles};
			if (scalar(@unique) != 2){
				if (scalar(@unique) == 3 && $nm == 0){
					my $check_missing = -1;
					foreach my $a (0..2){
						if ($unique[$a] eq '.'){
							$check_missing = $a;
						}
					}
					if ($check_missing == -1){
						next;
					}
					else {
						@unique = splice(@unique, $check_missing, 1);
						@unique = sort {$a<=>$b} @unique;
					}
				}
				else {
					next;
				}
			}
			else {
				if ($nm == 0){
					my $check_missing = -1;
					foreach my $a (0..1){
						if ($unique[$a] eq '.'){
							$check_missing = 1;
						}
					}
					if ($check_missing == 1){
						next;
					}
				}			
			}
			$eles[4] = $sorted_nucls[$unique[1]-1];
		}
		my $joint_out;
		if ($sf == 1 || $thap == 1){ #haps format for Relate, thap format for rehh
			if ($sf == 1){
				print OUT "$eles[0] $eles[2] $eles[1] $eles[3] $eles[4] ";
			}
			elsif ($thap == 1){
				print SAMPLE "$eles[2] $eles[0] $eles[1] $eles[3] $eles[4]\n";
				
			}
			$joint_out = join(" ", @out_eles);
			$joint_out =~ s/\/|\|/ /g;
			if ($bi == 1){
				$joint_out =~ s/[2-9]/1/;
			}
		}
		else { #vcf format
			foreach my $q (0..8){
				print OUT "$eles[$q]\t";
				#print "debug\n";
			}
			$joint_out = join("\t", @out_eles);
			if ($bi == 1){
				$joint_out =~ s/[2-9]/1/;
			}
		}
		print OUT "$joint_out\n";
		$cnt++;
	}
}
if ($thap == 1){
	close(SAMPLE);
}
close(INPUT);
close(OUT);
print "Done.\n";

sub get_ancestral_allele {
	my $line = shift; my @aids = @{$_[-1]};
	my @values;
	my @eles = split(/\t/, $line);
	my @info = split(/\:/, $eles[8]);
	my @nucls = split(/\,/, $eles[4]); #ALT
	unshift(@nucls, $eles[3]); #REF
	my $ad;
	#my $anc_allele = $eles[$aid];
	my @anc_alleles;
	foreach my $n (0..$#aids){
		if ($aids[$n] !~ /[^0-9]/){
			push(@anc_alleles, $eles[$aids[$n]]);
		}
	}
	foreach my $k (0..$#info){
		if ($info[$k] eq "AD"){
			$ad = $k;
		}
	}
	my @anc_nucl_array;
	foreach my $o (0..$#anc_alleles){
		@anc_info = split(/\:/, $anc_alleles[$o]);
		my $anc_number;
		my @alleles = split(/\/|\|/, $anc_info[0]);
		my $anc_nucl;
		if ($ad =~ /[0-9]/){
			my @depths = split(/\,/, $anc_info[$ad]);
			foreach my $j (0..$#depths){
				if ($depths[$j] eq "."){
					$depths[$j] = 0;
				}
			}
			if (scalar(@alleles) == 2){
				my @sorted_depths = sort {$depths[$b] <=> $depths[$a]} 0..$#depths;
				my @sorted_alleles = @alleles[@sorted_depths];
				if ($sorted_depths[0] == 0){
					$anc_number = "N";
				}
				$anc_number = @sorted_alleles[0];
			}
			else {
				$anc_number = $alleles[0];
			}
			if ($anc_number eq "N"){
				$anc_nucl = $anc_number;
			}
			else {
				$anc_nucl = $nucls[$anc_number];
			}
		}
		else {
			$anc_number = $alleles[0];
			$anc_nucl = $nucls[$anc_number];
		}
		push(@anc_nucl_array, $anc_nucl);
	}
	#count numbers of the same elements in the array
	my %counts = ();
	foreach (@anc_nucl_array){
		$counts{$_}++;
	}
	my @keys = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
	my @vals = @counts{@keys};
	return $keys[0];
}

sub rename_chr {
	my $name = shift; my @lists = @{$_[-1]};
	my $ori_name = $name;
	foreach (@lists){
		my @eles = split(/\t+|\s+/, $_);
		if ($name eq $eles[0]){
			$name = $eles[1];
			last;
		}
	}
	if ($name eq $ori_name){
		print RED "$ori_name is unchanged.\n", RESET;
	}
	return $name;
}

sub modify_chr {
	my $name = shift;
	my $ori_name = $name;
	if ($name !~ /\D/){
		if ($name =~ /^0+/){
			$name =~ s/^0+//;
		}
		return $name;
	}
	else {
		$name =~ s/\D+/ /g;
		if ($name =~ /\s+/){
			my @tmps = split(/\s+/, $name);
			$name = $tmps[-1];
			if ($name =~ /^0+/){
				$name =~ s/^0+//;
			}
			if ($name !~ /\d/){
				$name = $tmps[-2];
					if ($name =~ /^0+/){
						$name =~ s/^0+//;
					}
				if ($name !~ /\d/){
					$name = $ori_name;
				}
			}
		}
		#print "debug: Name: $name\n";
		return $name;
	}
}

