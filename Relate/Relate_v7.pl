#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/softwares\/qsub_subroutine.pl"){
	require "$home\/softwares\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
else {
	die "Cannot find required subroutine file: qsub_subroutine.pl\n";
}

#detect server
my @server = `ip route get 1.2.3.4 \| awk \'\{print \$7\}\'`;
chomp(@server);
my $conda = '-cj_conda R-4.1 ';
my $clues_plot = 0;
foreach (@server){
	if ($_ =~ /140.112.2.73/){
		$conda = '-cj_env /home/hpc/crlee/miniconda3/envs/R-4.1/bin ';
	}
	if ($_ =~ /140.112.2.71/){
		$clues_plot = 0;
	}
	if ($_ =~ /140.112.2.90/){
		$clues_plot = 0;
		$conda = "";
	}
}

my $Relate = "\$Relate";

chomp(@ARGV);
print "The script is written by Ben Chien. Jan. 2025.\n";
print "Input command line:\n";
print "perl Relate_v7\.pl @ARGV\n";

if ($#ARGV == -1){
	&usage;
	exit;
}

my $exc; my $ran; my $sn; my $ow = 0; my $o_path; my $al;
my $inprefix; my $path; my @vcfs;
my $mask; my $rm; my $pop; my @chr_names; my $rr = 0; my $pop_check;
my $mrate = '1e-8'; my $es = 3000; my $hap = 0;
my $coal; my @popis; my $year = "--years_per_gen 1 "; my $td = "--threshold 0 ";
my $repeat = 1; my $spl = 0; my @ancs; my $bin = -1;
my $eps = 1; my $dps = 0; my $clues = 0; my @cbins; my $epsrp = 100;
my @clues_chr; my @clues_fbp; my @clues_lbp;
my $ns = "--num_samples 5 "; my $tvs = 0; my $mem = 12;
my $rc = 1; my $cf = "b"; my $tco; my $dom;
my $qout; my $replot; my @cps; my $df; my $nat;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-exc"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "-sn"){
		$ran = $ARGV[$i+1];
		$sn = 1;
	}
	if ($ARGV[$i] eq "\-h"){
		&usage;
		exit;
	}
	if ($ARGV[$i] eq "\-ow"){
        $ow = 1;
	}
	if ($ARGV[$i] eq "\-mem"){
        $mem = $ARGV[$i+1];
        unless ($mem =~ /\d/ && $mem <= 128 && $mem >= 1){
			die "-mem: only digitals can be accepted. (1-128)\n";
		}
	}
	if ($ARGV[$i] eq "\-o"){
		$o_path = $ARGV[$i+1];
		if ($o_path =~ /\/$/){
			$o_path =~ s/\/$//;
		}
		unless (-d $o_path){
			system("mkdir $o_path");
		}
		$o_path = &check_path_Relate($o_path);
	}
	if ($ARGV[$i] eq "\-mask"){
        $mask = $ARGV[$i+1];
        unless (-e $mask){
        	die "Cannot find the mask file.\n";
        }
        $mask = &check_path_Relate($mask);
	}
	if ($ARGV[$i] eq "\-rm"){
		unless (-e $ARGV[$i+1]){
			die "Cannot find the poplabels file.\n";
		}
		$ARGV[$i+1] = &check_path_Relate($ARGV[$i+1]);	
        $rm = "--remove_ids $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-pop"){
		unless (-e $ARGV[$i+1]){
			die "Cannot find the poplabels file.\n";
		}
		$ARGV[$i+1] = &check_path_Relate($ARGV[$i+1]);
		$pop_check = $ARGV[$i+1];
        $pop = "--poplabels $ARGV[$i+1] ";
	}
    if ($ARGV[$i] eq "\-vcf"){
        if (-d $ARGV[$i+1]){
            if ($ARGV[$i+1] =~ /\/$/){
                $ARGV[$i+1] =~ s/\/$//;
            }
            $path = $ARGV[$i+1];
            unless ($path) {
                $path = ".";
            }
            $path = &check_path($path);
            @vcfs = <$path\/*.vcf>;
            unless (@vcfs){
                @vcfs = <$path\/*.vcf.gz>;
                unless (@vcfs){
                    die "Cannot find the vcf\(s\).\n";
                }
            }
            @vcfs = grep {!/\.relate\.|\.NumericChr\./} @vcfs;
        }
        elsif (-e $ARGV[$i+1] && ($ARGV[$i+1] =~ /txt$/ || $ARGV[$i+1] =~ /list$/)){
            open(LIST, "<$ARGV[$i+1]") || die "Cannot open $ARGV[$i+1]: $!\n";
            @vcfs = <LIST>;
            chomp(@vcfs);
            foreach my $i (0..$#vcfs){
                $vcfs[$i] =~ s/[\x0A\x0D]//g;
            }
            if ($vcfs[0] =~ /\//){
                my @tmp = split(/\//, $vcfs[0]);
                pop(@tmp);
                $path = join("\/", @tmp);
            }
            else {
                $path = ".";
            }
			$path = &check_path_Relate($path);
			@vcfs = grep {!/\.relate\.|\.NumericChr\./} @vcfs;
        }
        else {
            die "Cannot find the vcf\(s\). The input must be a folder contains seperated vcfs or a list file.\n";
        }
        $inprefix = "Relate_output";
    }
	if ($ARGV[$i] eq "\-al"){ #ancestor id list file
        unless (-e $ARGV[$i+1]){
        	die "Cannot find the ancestor ID list file: $ARGV[$i+1]\n";
        }
        $ARGV[$i+1] = &check_path_Relate($ARGV[$i+1]);
        $al = $ARGV[$i+1];
	}
    if ($ARGV[$i] eq "\-pre"){
    	$pre = $ARGV[$i+1];
    }
	if ($ARGV[$i] eq "\-m"){
        $mrate = $ARGV[$i+1];
		if ($mrate =~ /[^0-9e\-\.]/ || $mrate =~ /^\-/){
            die "Value of -m is not right.\n";
        }
	}
    if ($ARGV[$i] eq "\-n"){ #effective size
        $es = $ARGV[$i+1];
        if ($es =~ /[^0-9]/){
            die "Value of -n is not right.\n";
        }
    }
    if ($ARGV[$i] eq "\-epsrp"){ #repeat number for estimate population numer (this can draw 95% CI)
        $epsrp = $ARGV[$i+1];
        if ($epsrp =~ /[^0-9]/){
            die "Value of -epsrp is not right.\n";
        }
    }
	if ($ARGV[$i] eq "\-map"){
        $gmap = $ARGV[$i+1];
        unless (-e $gmap){
        	die "Cannot find the genetic map.\n";
        }
        $gmap = &check_path_Relate($gmap);
	}
	if ($ARGV[$i] eq "\-rr"){ #rerun Relate main program
        $rr = 1;
	}
	if ($ARGV[$i] eq "\-coal"){ #coal file
        unless (-e $ARGV[$i+1]){
        	die "Cannot find the coal file: $ARGV[$i+1]\n";
        }
        $ARGV[$i+1] = &check_path_Relate($ARGV[$i+1]);
        $coal = "--coal $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-popi"){ #interest pop(s)
		if ($ARGV[$i+1] =~ /\,/){
			my @tmp_popis = split(/\,/, $ARGV[$i+1]);
			my $left = 0;
			my $multi_pops;
			foreach (@tmp_popis){
				if ($_ !~ /\[|\]/ && $left == 0){
					push(@popis, "--pop_of_interest $_ ");
				}
				elsif ($_ =~ /\[/){
					$_ =~ s/\[//;
					$multi_pops = $_;
					$left = 1;
				}
				elsif ($_ !~ /\[|\]/ && $left == 1){
					$multi_pops .= "\,$_";
				}
				elsif ($_ =~ /\]/){
					$_ =~ s/\]//;
					$multi_pops .= "\,$_";
					push(@popis, "--pop_of_interest $multi_pops ");
					$left = 0;
					$multi_pops = "";
				}
			}
		}
		else {
			$ARGV[$i+1] =~ s/\[|\]//g;
			@popis = "--pop_of_interest $ARGV[$i+1] ";
		}
	}
	if ($ARGV[$i] eq "\-year"){ #year of generation
		if ($ARGV[$i+1] =~ /[^0-9\.]/ || $ARGV[$i+1] <= 0){
			die "-year value is not right.\n";
		}
        $year = "--years_per_gen $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-rc"){ #repeat run for clues
        $rc = $ARGV[$i+1];
        if ($rc =~ /[^0-9]/ || $rc <= 0){
            die "Value of -n is not right.\n";
        }
	}
	if ($ARGV[$i] eq "\--force"){ #force repeat run > 1000
        $force = 1;
	}	
	if ($ARGV[$i] eq "\-spl"){ #for spline function of recombinant map
		if ($ARGV[$i+1] =~ /^\-/){
			$spl = 25;
		}
		elsif ($ARGV[$i+1] =~ /[^0-9]/){
			die "-spl value should be an integer.\n";
		}
		else {
			$spl = int($ARGV[$i+1]);
		}
	}
	if ($ARGV[$i] eq "\-hap"){ #if haploid
        $hap = 1;
	}
	if ($ARGV[$i] eq "-am"){
		if (-d $ARGV[$i+1]){
            if ($ARGV[$i+1] =~ /\/$/){
                $ARGV[$i+1] =~ s/\/$//;
            }			
		}
		@ancs = <$ARGV[$i+1]\/*.anc>;
		unless (@ancs){
			die "Cannot find anc files.\n";
		}
	}
	if ($ARGV[$i] eq "\-bins"){ #set boundary of the year for -eps
        $bin = $ARGV[$i+1];
        if ($bin !~ /\,/){
        	die "-bins format is wrong.\n";
        }
        my @bin_check = split(/\,/, $bin);
        if (scalar(@bin_check) == 3){
        	foreach (@bin_check){
        		unless ($_ > 0 && $_ <= 10){
        			die "-bin format is wrong.\n";
        		}
        	}
        }
        else {
        	die "-bin format is wrong.\n";
        }
	}
	if ($ARGV[$i] eq "\-dps"){ #detect positive selection
        $dps = 1;
	}
	if ($ARGV[$i] eq "\-clues"){ #CLUES input
        $clues = 1;
	}
	if ($ARGV[$i] eq "\-bp"){
		if ($ARGV[$i+1] !~ /\:/){
			die "-bp: Chromosome\/contig name is required.\n";
		}
		#if ($ARGV[$i+1] !~ /\-/){
		#	die "-bp: Must define a region with start and end positions.\n";
		#}
		my @bps;
		if ($ARGV[$i+1] =~ /\,/){
			@bps = split(/\,/, $ARGV[$i+1]);
		}
		else {
			@bps = $ARGV[$i+1];
		}
		foreach (@bps){
			my @eles = split(/\:|\-/, $_);
			if (scalar(@eles) != 3){
				push(@eles, $eles[1]);
				#die "-bp: Format is wrong. Please check.\n";
			}
			push(@clues_chr, $eles[0]);
			if ($eles[1] !~ /[^0-9]/){
				push(@clues_fbp, $eles[1]);
			}
			else {
				die "-bp: Format is wrong. Please check.\n";
			}
			if ($eles[2] !~ /[^0-9]/){
				push(@clues_lbp, $eles[2]);
			}
			else {
				die "-bp: Format is wrong. Please check.\n";
			}
			if ($eles[1] > $eles[2]){
				die "-bp: the start position must be equal or smaller than the end position.\n";
			}
		}
		undef(@bps);
	}
	if ($ARGV[$i] eq "\-ns" || $ARGV[$i] eq "\--num_samples"){ #for CLUES
		if ($ARGV[$i+1] =~ /[^0-9]/){
			die "-ns should be an integer.\n";
		}
		elsif ($ARGV[$i+1] <= 0){
			die "-ns should be >= 1.\n";
		}
		$ns = "--num_samples $ARGV[$i+1] ";
	}
    if ($ARGV[$i] eq "\-df"){ #effective size
        $df = $ARGV[$i+1];
        if ($df =~ /[^0-9]/){
            die "Value of -df is not right.\n";
        }
        $df = "--df $df ";
    }
	if ($ARGV[$i] eq "\-tco" || $ARGV[$i] eq "\--tCutoff"){ #CLUES plot: tCutoff
		if ($ARGV[$i+1] =~ /[^0-9]/){
			die "-tco should be an integer.\n";
		}
		$tco = "--tCutoff $ARGV[$i+1] ";
	}
=cut
	if ($ARGV[$i] eq "\-pf" || $ARGV[$i] eq "\--popFreq"){ #CLUES plot: popFreq
		if ($ARGV[$i+1] =~ /[^0-9\.]/ || $ARGV[$i+1] > 1 || $ARGV[$i+1] <= 0){
			die "-pf should be a float and between 0 and 1.\n";
		}
		$pf = "--popFreq $ARGV[$i+1] ";
	}
=cut
	if ($ARGV[$i] eq "\-d" || $ARGV[$i] eq "\--dom"){ #CLUES plot: dom
		if ($ARGV[$i+1] =~ /[^0-9\.]/ || $ARGV[$i+1] > 1 || $ARGV[$i+1] < 0){
			die "-dom should be a float and between 0 and 1.\n";
		}
		$dom = "--h $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-rp" || $ARGV[$i] eq "\--re-plot"){ #re-plot output
		if ($ARGV[$i+1] =~ /clues|all/i){
			$replot = $ARGV[$i+1];
		}
		else {
			die "-rp: Cannot recognize the input. Possible value: clues, all\n";
		}
	}
	if ($ARGV[$i] eq "\-rpp" || $ARGV[$i] eq "\--re-plot-py"){ #re-plot output *.png
		$clues_plot = 1;
	}
	if ($ARGV[$i] eq "\-nat"){ #CLUES2: --noAlleleTraj
        $nat = "--noAlleleTraj ";
	}
	if ($ARGV[$i] eq "\-cp" || $ARGV[$i] eq "\--color-palette"){ #plot color series
		my $cp = $ARGV[$i+1];
		if ($cp =~ /\,/){
			@cps = split(/\,/, $cp);
			foreach my $j (0..$#cps){
				$cps[$j] = "-c $cps[$j] ";
			}
		}
	}
	if ($ARGV[$i] eq "\-cb" || $ARGV[$i] eq "\--clue-bins"){ #set boundary of the year for CLUES plotting
        my $cbin = $ARGV[$i+1];
        @cbins = split(/\,/, $cbin);
        my $discard = 0;
        foreach my $j (0..$#cbins){
			unless (-e $cbins[$j]){
				if ($cbins[0] !~ /[^0-9\.]/){
					my $tmp = join(" ", @cbins);
					#print "debug: $tmp\n";
					if ($tmp =~ /[a-z]/i){
						die "-cb: it seems the format is worng.\n";
					}
					$cbins[0] = "--timeBins $tmp ";
					$discard = 1;
				}
			}
			else {
				open(BINS, "<$cbins[$j]") || die "Cannot open $cbins[$j]: $!\n";
				my @timeframe = <BINS>;
				close(BINS);
				chomp(@timeframe);
				foreach my $x (0..$#timeframe){
					my @tmp = split(/[\s\t\,\;]/, $timeframe[$x]);
					$timeframe[$x] = join(" ", @tmp);
				}
				$cbins[$j] = join(" ", @timeframe);
				$cbins[$j] =~ s/^\s+|\s+$//;
				$cbins[$j] =~ s/\s+/ /;
				$cbins[$j] = "--timeBins $cbins[$j] ";
			}
        }
        if ($discard == 1){
			@cbins = $cbins[0];
        }
	}
    if ($ARGV[$i] eq "\-cf"){
    	$cf = $ARGV[$i+1];
    	unless ($cf == "a" || $cf == "n"){ #CLUES: output format for Sample branch lengths function
    		die "-cf: only support format a and n now.\n";
    	}
    }
	if ($ARGV[$i] eq "\-tvs"){ #TreeViewSamples
		$tvs = 1;
		$clues = 1;
	}
}

if ($rc > 1000){
	if ($force != 1){
		die "-rc bigger than 1000. Use --force to force run the repeats.\n";
	}
}
if ($clues == 1){
	$dps = 1;
}

if ($tvs == 1){
	$cf = "a";
}
else {
	$cf = "n";
}

unless ($gmap){
	die "-map is not supplied.\n";
}

unless ($pop){
	die "-pop is not supplied.\n";
} 
else { #check if the population of interest exists
	open(POP, "<$pop_check") || die "Cannot open $pop_check: $!\n";
	my @pop_content = <POP>;
	close(POP);
	chomp(@pop_content);
	my @uni_pop;
	if ($pop){
		shift(@pop_content);
		foreach (@pop_content){
			my @tmps = split(/\t|\s+/, $_);
			push(@uni_pop, $tmps[1]);
		}
	}
	else {
		foreach (@pop_content){
			my @tmps = split(/\t|\s+/, $_);
			push(@uni_pop, $tmps[3]);
		}	
	}
	@uni_pop = do { my %seen; grep { !$seen{$_}++ } @uni_pop };
	my $uni_p = join(" ", @uni_pop);
	if (@popis){
		foreach (@popis){
			my @tmp_ck;
			if ($_ =~ /\,/){
				@tmp_ck = split(/\,/, $_);
			}
			else {
				@tmp_ck = $_;
			}
			foreach my $k (0..$#tmp_ck){
				$tmp_ck[$k] =~ s/--pop_of_interest//;
				$tmp_ck[$k] =~ s/\s+//g;
				if ($uni_p !~ /\b$tmp_ck[$k]\b/){
					die BOLD "-popi: Cannot find the population of interest: $tmp_ck[$k]\nPlease check.", RESET, "\n";
				}
			}
		}
	}
}

unless (@vcfs && $gmap && $pop){
	unless (@ancs && $pop){
		die "-am or -pop is not supplied.\n";
	}
}

if (@ancs){
    if ($ran){
        @ancs = grep(/$ran/, @ancs);
    }
}
else {
	unless ($al){
		die "Ancestral ID list \(-al\) is required.\n";
	}
	unless (@vcfs){
		die "-vcf is not supplied.\n";
	}
}

if (@clues_chr){
	if (@cps){
		if (scalar(@cps) == 1){
			foreach my $i (1..$#cps){
				push(@cps, $cps[0]);
			}
		}
		elsif (scalar(@clues_chr) != scalar(@cps)){
			die "-cp: The number of defined color should be 1 or the same with -bp.\n";
		}
	}
}

if (@cbins){
	if (scalar(@cbins) == 1){
		unless (@clues_chr){
			die "-cb: -bp must be set.\n";
		}
		foreach my $i (1..$#clues_chr){
			push(@cbins, $cbins[0]);
		}
	}
	elsif (scalar(@bps) != scalar(@cbins)){
		die "-cb: The number of defined timebins should be 1 or the same with -bp.\n";
	}
}

unless ($o_path){
	$o_path = getcwd;
}

unless (-d "qsub_files"){
	system ("mkdir qsub_files");
}
unless (-d "qsub_files\/out"){
	system ("mkdir qsub_files\/out");
}

RE:
unless ($ran){
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "my_bash_relate_input_$ran\.sh" && $sn != 1){
	$ran = undef;
	goto RE;
}
print "The qsub SN is: $ran\n";

my @r_inputs;
if (@ancs){
	foreach my $i (1..scalar(@ancs)){
		my @pre_name = split(/_chr/, $ancs[$i-1]);
		if (-e "$pre_name[0]\_chr$i.anc" && -e "$pre_name[0]\_chr$i.mut"){
			push(@r_inputs, "$pre_name[0]\_chr$i");
		}
		elsif (-e "$pre_name[0]\_chr$i.anc"){
			die "Cannot find corresponding $pre_name[0]\_chr$i.mut of $pre_name[0]\_chr$i.anc.\n";
		}
		elsif (-e "$pre_name[0]\_chr$i.mut"){
			die "Cannot find corresponding $pre_name[0]\_chr$i.anc of $pre_name[0]\_chr$i.mut.\n";
		}
	}
	if (@r_inputs){
		goto ADDON;
	}
}


@chr_names = &chr_name($vcfs[0], $pre);
if (scalar(@vcfs) == 1){
	my @vcfs_tmp;
	open (BASH, ">my_bash_seperate_vcf_$ran\.sh") || die BOLD "Cannot write my_bash_seperate_vcf_$ran\.sh: $!", RESET, "\n";
	foreach (@chr_names){
		print BASH "qsub \.\/qsub_files\/$ran\_seperate_vcf_$_\.q\n";
		$qout = "vcftools --gzvcf $vcfs[0] --recode --recode-INFO-all --chr $_ --stdout \| bgzip -c \> $path\/$ran\_$_.input.vcf.gz\\n";
		push(@vcfs_tmp, "$path\/$ran\_$_.input.vcf.gz");
		&pbs_setting("$exc\-cj_quiet -cj_qname seperate_vcf_$_ -cj_ppn 4 -cj_mem $mem -cj_sn $ran -cj_qout . $qout");
	}
	close(BASH);
	if ($exc){
		&status($ran);
	}
	@vcfs = @vcfs_tmp;
	undef(@vcfs_tmp);
}
@vcfs = &get_1_eles(\@vcfs, \@chr_names);
open(RNAME, ">$o_path\/rename_chr.list") || die "Cannot write $o_path\/rename_chr.list: $!\n";
foreach my $i (0..$#chr_names){
	my $new_idx = $i+1;
    print RNAME "$chr_names[$i] $new_idx\n"; #rename contig name list, start from 1
}
close(RNAME);
my @masks;
$qout = "";
if ($mask){
	my $tmp_mask = 1;
	unless (-d "$o_path\/masked_ref"){
		system("mkdir $o_path\/masked_ref");
	}
	foreach my $i (1..$#chr_names+1){
		unless (-e "$o_path\/masked_ref\/$i.masked.fasta" && $ow == 0){
			$tmp_mask = 0;
			last;
		}
	}	
	if ($tmp_mask == 0){
		$qout .= "perl convert_masked_fasta.pl $mask $o_path\/masked_ref $pre\\n"; #seperate contigs, convert [ATCG] to P, rename contig name to numeric
	}
	foreach my $i (1..$#chr_names+1){
		push(@masks, "$o_path\/masked_ref\/$i.masked.fasta");
	}
}
if ($gmap){
	my $tmp_gmap = 1;
	unless (-d "$o_path\/genetic_maps"){
		system("mkdir $o_path\/genetic_maps");
	}
	foreach my $i (1..$#chr_names+1){
		unless (-e "$o_path\/genetic_maps\/$i.gmap.map" && $ow == 0){
			$tmp_gmap = 0;
			last;
		}
	}
	if ($tmp_gmap == 0){
		if ($spl == 0){
			$qout .= "Rscript recomb_spline_Relate.R $gmap $o_path\/genetic_maps 1 0 $pre\\n";
		}
		elsif ($spl > 0){
			$qout .= "Rscript recomb_spline_Relate.R $gmap $o_path\/genetic_maps 0 $spl $pre\\n";
		}	
	}
}
&pbs_setting("$exc$conda\-cj_quiet -cj_qname relate_seperate -cj_sn $ran -cj_qout . $qout");
if ($exc){
    &status($ran);
}

my $check = 1;
foreach my $i (1..$#chr_names+1){
	unless (-e "$o_path\/inputs\/$ran\_input_chr$i.haps.gz" && -e "$o_path\/inputs\/$ran\_input_chr$i.sample.gz" && $ow == 0){
		$check = 0;
		last;
	}
}
if ($check == 1){
	foreach my $j (1..$#chr_names+1){
		push(@r_inputs, "$o_path\/inputs\/$ran\_input_chr$j");
	}
	goto RELATE;
}

open (BASH, ">my_bash_relate_input_$ran\.sh") || die BOLD "Cannot write my_bash_relate_input_$ran\.sh: $!", RESET, "\n";
my $cnt = 0;
print "Preparing input files...\n";
foreach (@vcfs){
	$qout = "";
	if ($pre){
		my $re = &get_first_ele_Relate($_);
		if ($re !~ /^$pre/){
			next;
		}
	}
	$cnt++;
    my $num;
    foreach my $i (0..$#chr_names){ #chr_names is still the original names, $_ is $vcf that contains the original names
    	$first_ele = &get_first_ele_Relate($_);
        if ($first_ele eq $chr_names[$i]){
            $num = $i+1;
        }
    }
    unless (-d "$o_path\/inputs"){
    	system("mkdir $o_path\/inputs");
    }
    my $hap_check;
    if (-e "$o_path\/inputs\/$ran\_input_chr$num.haps.gz" && -e "$o_path\/inputs\/$ran\_input_chr$num.sample.gz" && $ow == 0){
		$hap_check = `zcat $o_path\/inputs\/$ran\_input_chr$num.haps.gz \| head -n 5`;
		if ($hap_check =~ /\w/){
			push(@r_inputs, "$o_path\/inputs\/$ran\_input_chr$num");
			next;
        }
        else {
			system("rm $o_path\/inputs\/$ran\_input_chr$num.\*");
        }
    }
	if (-e "$o_path\/inputs\/$ran\_input_chr$num.haps" && -e "$o_path\/inputs\/$ran\_input_chr$num.sample" && $ow == 0){
		$hap_check = `head -n 5 $o_path\/inputs\/$ran\_input_chr$num.haps`;
		if ($hap_check =~ /\w/){
			push(@r_inputs, "$o_path\/inputs\/$ran\_input_chr$num");
			next;
        }
        else {
			system("rm $o_path\/inputs\/$ran\_input_chr$num.\*");
        }
    }
    #convert vcf files into haps and sample files
    unless (-e "$o_path\/inputs\/$ran\_input_chr$num.haps" &&  -e "$o_path\/inputs\/$ran\_input_chr$num.sample" && $ow == 0){
        	my $anc_hap;
        	if ($hap == 1){
        		$anc_hap = "-hap ";
        	}
        	$qout .= "perl vcf2anc_vcf.pl -vcf $_ -aid $al -keep -bi -nm -list $pop_check -sf $anc_hap\-rchr $o_path\/rename_chr.list -o $o_path\/inputs\/$ran\_input_chr$num\\n";
    }
    #prepare input files
    if (@masks){
    	unless (-e "$o_path\/inputs\/$ran\_input_chr$num.dist" && $ow == 0){
            $qout .= "$Relate\/bin\/RelateFileFormats --mode FilterHapsUsingMask --haps $o_path\/inputs\/$ran\_input_chr$num.haps --sample $o_path\/inputs\/$ran\_input_chr$num.sample --mask $masks[$num-1] -o $o_path\/inputs\/$ran\_input_mask_chr$num\\n";
        	$qout .= "mv $o_path\/inputs\/$ran\_input_mask_chr$num.haps $o_path\/inputs\/$ran\_input_chr$num.haps\\n";
        	$qout .= "mv $o_path\/inputs\/$ran\_input_mask_chr$num.dist $o_path\/inputs\/$ran\_input_chr$num.dist\\n";
        }
    }
    if ($pop){
    	unless (-e "$o_path\/inputs\/$ran\_input_chr$num.annot" && $ow == 0){
    		$qout .= "$Relate\/bin\/RelateFileFormats --mode GenerateSNPAnnotations --haps $o_path\/inputs\/$ran\_input_chr$num.haps --sample $o_path\/inputs\/$ran\_input_chr$num.sample $pop -o $o_path\/inputs\/$ran\_input_chr$num\\n";
    	}
    }
    push(@r_inputs, "$o_path\/inputs\/$ran\_input_chr$num");
    &pbs_setting("$exc\-cj_quiet -cj_qname relate_input_$cnt -cj_ppn 4 -cj_mem $mem -cj_sn $ran -cj_qout . $qout");
}
close(BASH);
if ($exc){
    &status($ran);
    print "Input files prepared.\n";
}

RELATE:
#make relate output data
print "Preparing anc and mut files...\n";
open (BASH, ">my_bash_relate_main_$ran\.sh") || die BOLD "Cannot write my_bash_relate_main_$ran\.sh: $!", RESET, "\n";
$cnt = 1;
foreach (@r_inputs){ #chromosomes
	$qout = "";
    print BASH "qsub \.\/qsub_files\/$ran\_relate_main_$cnt\.q\n";
    my $dist; my $out; my $pre_path;
    if ($_ =~ /\//){
        my @tmp = split(/\//, $_);
        $out = pop(@tmp);
        $pre_path = join("\/", @tmp);
    }
    else {
        $out = $_;
        $pre_path = ".";
    }
    unless (-e "$_.anc" && -e "$_.mut" && $ow == 0 && $rr == 0){
    	if ($mask){
    		if (-e "$_.dist"){
    			$dist = "--dist $_.dist ";
    		}
    	}
    	if (-d $out){
			system("rm -r $out");
    	}
    	my $relate_map;
    	if ($gmap){
    		my $check_map = `cat $o_path\/genetic_maps\/$cnt.gmap.map`;
    		if ($check_map !~ /[0-9]/){
    			$relate_map = "";
    		}
    		else {
    			$relate_map = "--map $o_path\/genetic_maps\/$cnt.gmap.map ";
    		}
    	}
        $qout .= "$Relate\/bin\/Relate --mode All --memory $mem -m $mrate -N $es --haps $_.haps --sample $_.sample $relate_map\--annot $_.annot $dist$coal\-o $out\\n";
        if ($pre_path ne "."){
            $qout .= "mv -f  $out.anc $pre_path\/$out.anc\\n";
            $qout .= "mv -f $out.mut $pre_path\/$out.mut\\n";
        }
    }
    &pbs_setting("$exc\-cj_quiet -cj_mem $mem -cj_qname relate_main_$cnt -cj_sn $ran -cj_qout . $qout");
    $cnt++;
}
close(BASH);
if ($exc){
    &status($ran);
    print "anc and mut files prepared.\n";
}

ADDON:
foreach my $i (0..$#r_inputs){
    $r_inputs[$i] = &check_path_Relate($r_inputs[$i]);
}

#Estimate population size (basic)
my @r_inputs_ori = @r_inputs;
if ($eps == 1){
	my $popi;
	@r_inputs = &cal_eps($epsrp, $replot, $ran, $conda, $o_path, $mrate, $pop, $popi, $year, $td, $exc, $ow, $bin, \@r_inputs);
	#update mutation rate
}
#update mutation rate
if (-e "$o_path\/popsize\/$ran\_popsize_avg.rate"){
	open(MRATE, "<$o_path\/popsize\/$ran\_popsize_avg.rate") || die "Cannot open $ran\_popsize_avg.rate file: $!\n";
	my @tmp_mrates = <MRATE>;
	chomp(@tmp_mrates);
	close(MRATE);
	my $min_yr; my $max_yr; 
	foreach my $i (1..$#tmp_mrates-1){
		my @curr_mrate = split(/\t+|\s+/, $tmp_mrates[$i]);
		my @next_mrate = split(/\t+|\s+/, $tmp_mrates[$i+1]);
		if ($i == 1){
			$min_yr = $curr_mrate[0];
		}
		if ($next_mrate[1] =~ /-nan/){
			$max_yr = $next_mrate[0];
			last;
		}
	}
	my $yr_range = $max_yr - $min_yr;
	#print "debug 1: $max_yr $min_yr $yr_range\n";
	my $mrate_tmp = 0;
	foreach my $i (1..$#tmp_mrates-1){
		my @curr_mrate = split(/\t+|\s+/, $tmp_mrates[$i]);
		my @next_mrate = split(/\t+|\s+/, $tmp_mrates[$i+1]);
		$mrate_tmp = $mrate_tmp + $curr_mrate[1] * (($next_mrate[0] - $curr_mrate[0]) / $yr_range);
		if ($next_mrate[1] =~ /-nan/){
			last;
		}
		#print "debug 3: $mrate_tmp\n";
	}
	#print "debug 4: $mrate_tmp\n";
	$mrate = $mrate_tmp;
}

#detect positive selection
my $dps_folder = "DPS";
if ($dps == 1){
	my $popi_d_name; my @popi_d_names; my @o_path_dps;
	print "Generating commend lines of detect positive selection...\n";
	if ($clues == 1){
		$dps_folder = "CLUES";
	}
	unless (-d "$o_path\/$dps_folder"){
        system("mkdir $o_path\/$dps_folder");
    }
	if (@popis){
		foreach (@popis){
			$popi_d_name = $_;
			$popi_d_name =~ s/\,/\-/g;
			$popi_d_name =~ s/--pop_of_interest//g;
			$popi_d_name =~ s/\s+//g;
			unless (-d "$o_path\/$dps_folder\/$popi_d_name"){
        		system("mkdir $o_path\/$dps_folder\/$popi_d_name");
    		}
    		push(@popi_d_names, $popi_d_name);
    		push(@o_path_dps, "$o_path\/$dps_folder\/$popi_d_name");
		}
	}
	else {
		unless (-d "$o_path\/$dps_folder\/ALL"){
        	system("mkdir $o_path\/$dps_folder\/ALL");
    	}
    	push(@popi_d_names, "ALL");
    	push(@o_path_dps, "$o_path\/$dps_folder\/ALL");
	}
	my %popis_dps_out;
	my $plot_year = $year;
	$plot_year =~ s/--years_per_gen//;
	$plot_year =~ s/\s+//g;
	foreach my $l (0..$#popi_d_names){
		my @r_inputs_dps;
		if ($popi_d_names[$l] ne "ALL"){
			@r_inputs_dps = &subset_sample($ran, $conda, $o_path_dps[$l], $pop, $popis[$l], $exc, $ow, \@r_inputs);
		}
		else {
			@r_inputs_dps = @r_inputs;
		}
		my @popi_dps_outputs;
		open (BASH, ">my_bash_relate_DPS_$popi_d_names[$l]\_$ran\.sh") || die BOLD "Cannot write my_bash_relate_DPS_$popi_d_names[$l]\_$ran\.sh: $!", RESET, "\n";
		if (@clues_chr){
			if (-e "$o_path\/rename_chr.list"){
				open (NAMES, "<$o_path\/rename_chr.list") || die "Cannot open the rename list file, $o_path\/rename_chr.list: $!\n";
				my @tmps = <NAMES>;
				chomp(@tmps);
				close(NAMES);
				foreach (@tmps){
					my @tmps2 = split(/\t+|\s+/, $_);
					push(@ori_chr_names, $tmps2[0]);
					push(@mod_chr_names, $tmps2[1]);
				}
			}
			else {
				die "-dps: Cannot find the rename list file, $o_path\/rename_chr.list.\n";
			}
			foreach my $i (0..$#clues_chr){
				foreach my $j (0..$#ori_chr_names){
					if ($clues_chr[$i] eq $ori_chr_names[$j]){
						$clues_chr[$i] = $mod_chr_names[$j];
						last;
					}
				}
				unless (-e "$o_path_dps[$l]\/$ran\_$popi_d_names[$l]\_selection_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].sele" && $clues == 0 && $ow == 0){
					print "Detect selection region \(chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\) start...\n";
					$qout = "$Relate\/scripts/DetectSelection/DetectSelection.sh -i $o_path_dps[$l]\/popsize\/$ran\_popsize_chr$clues_chr[$i] -m $mrate $year\--first_bp $clues_fbp[$i] --last_bp $clues_lbp[$i] -o $o_path_dps[$l]\/$ran\_$popi_d_names[$l]\_selection_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\\n";
				}
				print BASH "qsub \.\/qsub_files\/$ran\_relate_DPS_$popi_d_names[$l]\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\.q\n";
				&pbs_setting("$exc\-cj_quiet -cj_qname relate_DPS_$popi_d_names[$l]\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] -cj_sn $ran -cj_qout . $qout");
				if (@{$popis_dps_out{$popi_d_names[$l]}}){
					push(@{$popis_dps_out{$popi_d_names[$l]}}, "$o_path_dps[$l]\/$ran\_selection_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].sele");
				}
				else {
					@{$popis_dps_out{$popi_d_names[$l]}} = @popi_dps_outputs;
				}
			}
		}
		else {
			foreach my $i (1..scalar(@r_inputs_dps)){
				$qout = "";
				print BASH "qsub \.\/qsub_files\/$ran\_relate_DPS_$popi_d_names[$l]\_$i\.q\n";
				unless (-e "$o_path\/DPS\/$popi_d_names[$l]\/$ran\_selection_chr$i.sele" && $clues == 0 && $ow == 0){
					$qout .= "$Relate\/scripts/DetectSelection/DetectSelection.sh -i $r_inputs_dps[$i-1] -m $mrate $year\-o $o_path\/DPS\/$popi_d_names[$l]\/$ran\_selection_chr$i\\n";
				}
				push(@popi_dps_outputs, "$o_path\/DPS\/$popi_d_names[$l]\/$ran\_selection_chr$i.sele");
				&pbs_setting("$exc\-cj_quiet -cj_qname relate_DPS_$popi_d_names[$l]\_$i -cj_sn $ran -cj_qout . $qout");
			}
			@{$popis_dps_out{$popi_d_names[$l]}} = @popi_dps_outputs;
		}
		close(BASH);
	}
	if ($exc){
    	&status($ran);
	}
	foreach (@popi_d_names){
		open (BASH, ">my_bash_relate_DPS_plot_$_\_$ran\.sh") || die BOLD "Cannot write my_bash_relate_DPS_plot_$_\_$ran\.sh: $!", RESET, "\n";
		my @curr_popi_dps_outputs = @{$popis_dps_out{$_}};
		if (@clues_chr){
			foreach my $j (0..$#curr_popi_dps_outputs){
				my $dps_cnt = $j + 1;
				my $out_name = $curr_popi_dps_outputs[$j];
				$out_name =~ s/$o_path_dps[$l]\/$ran\_selection_/$o_path\/DPS\/$_./;
				$out_name =~ s/\.sele$//;
				unless (-e "$out_name.selection.history.pdf" && ($ow == 0 || $replot eq "all")){
					$qout = "Rscript DPS_plot.R $plot_year $out_name @curr_popi_dps_outputs[$i]\\n";
					print BASH "qsub \.\/qsub_files\/$ran\_relate_DPS_plot_$_\_$dps_cnt.q\n";
					&pbs_setting("$exc$conda\-cj_quiet -cj_qname relate_DPS_plot_$_\_$dps_cnt -cj_sn $ran -cj_qout . $qout");
				}
			}
		}
		else {
			unless (-e "$o_path\/DPS\/$_.selection.history.pdf" && ($ow == 0 || $replot eq "all")){
    			$qout = "Rscript DPS_plot.R $plot_year $o_path\/DPS\/$_ @curr_popi_dps_outputs\\n";
    			print BASH "qsub \.\/qsub_files\/$ran\_relate_DPS_plot_$_.q\n";
    			&pbs_setting("$exc$conda\-cj_quiet -cj_qname relate_DPS_plot_$_ -cj_sn $ran -cj_qout . $qout");
			}
		}
		close(BASH);
	}
	if ($exc){
    	&status($ran);
	}
	print "Detect positive selection finished.\n";
}

#sample branch lengths for clues
my @plot_py; my @plot_r;
if ($clues == 1 || $tvs == 1){
#print "debug: $clues $tvs\n";
    foreach my $z (1..$rc){
        print "Generating commend lines of sample branch lengths for CLUES...\n";
        my $dist;
        unless (@clues_chr){
            die "-clues: -bp must be defined.\n";
        }
        if ($clues == 1){
            unless (-d "$o_path\/CLUES"){
                system("mkdir $o_path\/CLUES");
            }
        }
        if ($tvs == 1){
            unless (-d "$o_path\/TreeViewSamples"){
                system("mkdir $o_path\/TreeViewSamples");
            }
        }
        my @ori_chr_names; my @mod_chr_names;
        if (-e "$o_path\/rename_chr.list"){
            open (NAMES, "<$o_path\/rename_chr.list") || die "Cannot open the rename list file, $o_path\/rename_chr.list: $!\n";
            my @tmps = <NAMES>;
            chomp(@tmps);
            close(NAMES);
            foreach (@tmps){
                my @tmps2 = split(/\t+|\s+/, $_);
                push(@ori_chr_names, $tmps2[0]);
                push(@mod_chr_names, $tmps2[1]);
            }
        }
        else {
            die "-clues: Cannot find the rename list file, $o_path\/rename_chr.list.\n";
        }
        unless (@popis){ #without population of interests
            open (BASH, ">my_bash_relate_CLUES_$popi_d_names[$l]\_$ran\.sh") || die BOLD "Cannot write my_bash_relate_CLUES_$popi_d_names[$l]\_$ran\.sh: $!", RESET, "\n";
            if ($clues == 1){
                unless (-d "$o_path\/CLUES\/ALL"){
                    system("mkdir $o_path\/CLUES\/ALL");
                }
            }
            if ($tvs == 1){
                unless (-d "$o_path\/TreeViewSamples\/ALL"){
                    system("mkdir $o_path\/TreeViewSamples\/ALL");
                }			
            }
            unless ($coal){
                if (-e "$o_path\/popsize\/$ran\_popsize.coal"){
                    $coal = "--coal $o_path\/popsize\/$ran\_popsize.coal ";
                }
                else {
                    die "-coal is not defined.\n";
                }
            }
            foreach my $i (0..$#clues_chr){
                foreach my $j (0..$#ori_chr_names){
                    if ($clues_chr[$i] eq $ori_chr_names[$j]){
                        $clues_chr[$i] = $mod_chr_names[$j];
                        last;
                    }
                }
                #print "debug: $clues_chr[$i]\n";
                if ($tvs == 1){
                    unless (-e "$o_path\/TreeViewSamples\/ALL\/$ran\_$z\_TreeViewSamples_chr$clues_chr[$i].anc" && $ow == 0){
                        print "Sample branch lengths for TreeViewSamples \(chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\) start...\n";
                        $qout = "$Relate\/scripts/SampleBranchLengths/SampleBranchLengths.sh -i $o_path\/popsize\/$ran\_popsize_chr$clues_chr[$i] -m $mrate $ns$coal\--format $cf -o $o_path\/TreeViewSamples\/ALL\/$ran\_$z\_TreeViewSamples_chr$clues_chr[$i]\\n";
                    }
                    if ($clues_plot == 1){
                        my $tvs_input;
                        foreach my $z (0..$#r_inputs_ori){
                            if ($r_inputs_ori[$z] =~ /chr$clues_chr[$i]\b/){
                                $tvs_input = $r_inputs_ori[$z];
                            }
                        }
                        unless (-e "$o_path\/TreeViewSamples\/ALL\/$ran\_$z\_TreeViewSamples_chr$clues_chr[$i]\_$clues_fbp[$i].pdf" && $ow == 0){
                            $qout .= "$Relate\/scripts\/TreeView/TreeViewSample.sh --haps $tvs_input.haps --sample $tvs_input.sample --poplabels $o_path\/popsize\/$ran\_popsize_chr$clues_chr[$i].poplabels --anc $o_path\/TreeViewSamples\/ALL\/$ran\_$z\_TreeViewSamples_chr$clues_chr[$i].anc --mut $o_path\/TreeViewSamples\/ALL\/$ran\_$z\_TreeViewSamples_chr$clues_chr[$i].mut --dist $o_path\/popsize\/$ran\_popsize_chr$clues_chr[$i].dist --bp_of_interest $clues_fbp[$i] $year\-o $o_path\/TreeViewSamples\/ALL\/$ran\_$z\_TreeViewSamples_chr$clues_chr[$i]\_$clues_fbp[$i]\\n";
                        }
                    }
                    print BASH "qsub \.\/qsub_files\/$ran\_relate_TreeViewSamples_chr$clues_chr[$i]\_$clues_fbp[$i].q\n";
                    &pbs_setting("$exc\-cj_quiet $conda\-cj_qname $z\_relate_TreeViewSamples_chr$clues_chr[$i]\_$clues_fbp[$i] -cj_sn $ran -cj_qout . $qout");
                    $qout = "";
                }
                if ($clues == 1){
                    unless (-e "$o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].newick" && $ow == 0){
                        print "Sample branch lengths for CLUES \(chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\) start...\n";
                        $qout = "$Relate\/scripts/SampleBranchLengths/SampleBranchLengths.sh -i $o_path\/popsize\/$ran\_popsize_chr$clues_chr[$i] -m $mrate $ns\--first_bp $clues_fbp[$i] --last_bp $clues_lbp[$i] $coal\--format $cf -o $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\\n";
                    }
                    unless (-e "$o_path\/CLUES\/ALL\/chr$clues_chr[$i]\_$clues_fbp[$i].txt" && -e "$o_path\/CLUES\/ALL\/chr$clues_chr[$i]\_$clues_fbp[$i]\_freq.txt" && $ow == 0){
                        my $haplo;
                        if ($hap == 1){
                            $haplo = "-hap";
                        }
                        else {
                            $haplo = "";
                        }
                        $qout .= "perl extract_derived_file.pl -i $o_path\/inputs\/$ran\_input_chr$clues_chr[$i] -o $o_path\/CLUES\/ALL\/ -pos $clues_fbp[$i] $haplo\\n";
                    }
                    unless (-e "$o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\_times.txt"){
                        $qout .= "python $Relate\/clues\/RelateToCLUES.py --RelateSamples $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].newick --DerivedFile $o_path\/CLUES\/ALL\/chr$clues_chr[$i]\_$clues_fbp[$i].txt --out $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\\n";
                    }
                    print BASH "qsub \.\/qsub_files\/$ran\_relate_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\.q\n";
                    &pbs_setting("$exc\-cj_quiet $conda\-cj_qname $z\_relate_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] -cj_sn $ran -cj_qout . $qout");
                    if ($z == 1 && $exc){
                        &status($ran);
                    }
                    $qout = "";
                    my $pf;
                    if ($exc){
						open(FREQ, "<$o_path\/CLUES\/ALL\/chr$clues_chr[$i]\_$clues_fbp[$i]\_freq.txt") || die "Cannot open $o_path\/CLUES\/ALL\/chr$clues_chr[$i]\_$clues_fbp[$i]\_freq.txt: $!\n";
						$pf = <FREQ>;
						chomp($pf);
						close(FREQ);
						$pf = "--popFreq $pf ";
					}
					else {
						$pf = "--popFreq $pf ";
					}
                    unless (-e "$o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_post.txt" && $ow == 0){
                        #$qout .= "python $Relate\/clues\/inference.py --times $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --out $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $tco$pf$dom$coal$cbins[$i]\\n";
                        push(@plot_py, "python $Relate\/clues\/inference.py --times $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\_times.txt --out $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $nat$tco$df$pf$dom$coal$cbins[$i]\\n");
                    }
                    elsif ($replot =~ /clues|all/i){
                        #$qout .= "python $Relate\/clues\/inference.py --times $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --out $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $tco$pf$dom$coal$cbins[$i]\\n";
                        push(@plot_py, "python $Relate\/clues\/inference.py --times $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\_times.txt --out $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $nat$tco$df$pf$dom$coal$cbins[$i]\\n");
                    }
                    if ($clues_plot == 1){
                        unless (-e "$o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].pdf" && $ow == 0){
                            #$qout .= "Rscript clues_plot.R -f $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $cps[$i]\\n";
                            push(@plot_r, "python $Relate\/clues\/plot_traj_clues2.py --figure $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --post $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_post.txt --freqs $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_freqs.txt\\n");
                        }
                        elsif ($replot =~ /clues|all/i){
                            #$qout .= "Rscript clues_plot.R -f $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $cps[$i]\\n";
                            push(@plot_r, "python $Relate\/clues\/plot_traj_clues2.py --figure $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --post $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_post.txt --freqs $o_path\/CLUES\/ALL\/$ran\_$z\_CLUES_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_freqs.txt\\n");
                        }
                    }
                }
            }
            close(BASH);
        }
        else { #with population of interests
            my $popi_d_name; my @popi_d_names; my @o_path_clues; my @o_path_tvs;
            foreach (@popis){
                $popi_d_name = $_;
                $popi_d_name =~ s/\,/\-/g;
                $popi_d_name =~ s/--pop_of_interest//g;
                $popi_d_name =~ s/\s+//g;
                if ($clues == 1){
                    unless (-d "$o_path\/CLUES\/$popi_d_name"){
                        system("mkdir $o_path\/CLUES\/$popi_d_name");
                    }
                    push(@popi_d_names, $popi_d_name);
                    push(@o_path_clues, "$o_path\/CLUES\/$popi_d_name");
                }
                if ($tvs == 1){
                    unless (-d "$o_path\/TreeViewSamples\/$popi_d_name"){
                        system("mkdir $o_path\/TreeViewSamples\/$popi_d_name");
                    }
                    push(@popi_d_names, $popi_d_name);
                    push(@o_path_tvs, "$o_path\/TreeViewSamples\/$popi_d_name");				
                }
            }
            foreach my $l (0..$#popis){
                my @r_inputs_tvs;
                my $popi_out_name;
                $popi_out_name = $popis[$l];
                $popi_out_name =~ s/--pop_of_interest /_/;
                $popi_out_name =~ s/[\s\,]//g;
                if ($tvs == 1){
                    @r_inputs_tvs = &subset_sample($ran, $conda, $o_path_tvs[$l], $pop, $popis[$l], $exc, $ow, \@r_inputs);
                    &subset_inputs($ran, $conda, $o_path_tvs[$l], $pop, $popis[$l], $exc, $ow, \@r_inputs_ori);
                }
                if ($clues == 1 && $dps == 0){
                    &subset_sample($ran, $conda, $o_path_clues[$l], $pop, $popis[$l], $exc, $ow, \@r_inputs);
                }
                open (BASH, ">my_bash_relate_CLUES_$popi_d_names[$l]\_$ran\.sh") || die BOLD "Cannot write my_bash_relate_CLUES_$popi_d_names[$l]\_$ran\.sh: $!", RESET, "\n";
                my @popi_clues_outputs;
                foreach my $i (0..$#clues_chr){
                    foreach my $j (0..$#ori_chr_names){
                        if ($clues_chr[$i] eq $ori_chr_names[$j]){
                            $clues_chr[$i] = $mod_chr_names[$j];
                            last;
                        }
                    }
                    $coal = "--coal $o_path\/popsize\/$ran\_popsize.coal ";
                    if ($tvs == 1){
                        if (-e "$o_path_tvs[$l]\/popsize\/$ran\_popsize_avg.rate"){
                            my $tmp_mrate = `head -n 1 $o_path_tvs[$l]\/popsize\/$ran\_popsize_avg.rate`;
                            chomp($tmp_mrate);
                            my @tmps = split(/\t|\s+/, $tmp_mrate);
                            $mrate = $tmps[1];
                        }
                        unless (-e "$o_path_tvs[$l]\/$ran\_$z\_TreeViewSamples$popi_out_name\_chr$clues_chr[$i].anc" && $ow == 0){
                            print "Sample branch lengths for TreeViewSamples \($popi_d_names[$l]\: chr$clues_chr[$i]\) start...\n";
                            $qout = "$Relate\/scripts/SampleBranchLengths/SampleBranchLengths.sh -i $o_path_tvs[$l]\/popsize\/$ran\_popsize_chr$clues_chr[$i] -m $mrate $ns$coal\--format $cf -o $o_path_tvs[$l]\/$ran\_$z\_TreeViewSamples$popi_out_name\_chr$clues_chr[$i]\\n"; #This will do SampleBranchLengths for all positions in a chromosome/contig
                        }
                        if ($clues_plot == 1){
                            my $tvs_input;
                            foreach my $z (0..$#r_inputs_tvs){
                                $r_inputs_tvs[$z] =~ s/_popsize_/_input_/;
                                $r_inputs_tvs[$z] =~ s/popsize/inputs/;
                                if ($r_inputs_tvs[$z] =~ /chr$clues_chr[$i]\b/){
                                    $tvs_input = $r_inputs_tvs[$z];
                                }
                            }
                            unless (-e "$o_path_tvs[$l]\/$ran\_$z\_TreeViewSamples$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i].pdf" && $ow == 0){
                                $qout .= "$Relate\/scripts\/TreeView/TreeViewSample.sh --haps $tvs_input.haps --sample $tvs_input.sample --poplabels $o_path_tvs[$l]\/popsize\/$ran\_popsize_chr$clues_chr[$i].poplabels --anc $o_path_tvs[$l]\/$ran\_$z\_TreeViewSamples$popi_out_name\_chr$clues_chr[$i].anc --mut $o_path_tvs[$l]\/$ran\_$z\_TreeViewSamples$popi_out_name\_chr$clues_chr[$i].mut --dist $o_path_tvs[$l]\/popsize\/$ran\_popsize_chr$clues_chr[$i].dist --bp_of_interest $clues_fbp[$i] $year\-o $o_path_tvs[$l]\/$ran\_$z\_TreeViewSamples$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\\n";
                            }
                        }
                        #print "debug: $clues_chr[$i]\n";
                        print BASH "qsub \.\/qsub_files\/$ran\_relate_TreeViewSamples_$popi_d_names[$l]\_chr$clues_chr[$i]\_$clues_fbp[$i].q\n";
                        &pbs_setting("$exc\-cj_quiet $conda\-cj_qname $z\_relate_TreeViewSamples_$popi_d_names[$l]\_chr$clues_chr[$i]\_$clues_fbp[$i] -cj_sn $ran -cj_qout . $qout");
                        $qout = "";
                    }
                    if ($clues == 1){
                        if (-e "$o_path_clues[$l]\/popsize\/$ran\_popsize_avg.rate"){
                            my $tmp_mrate = `head -n 1 $o_path_clues[$l]\/popsize\/$ran\_popsize_avg.rate`;
                            chomp($tmp_mrate);
                            my @tmps = split(/\t|\s+/, $tmp_mrate);
                            $mrate = $tmps[1];
                        }
                        unless (-e "$o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].newick" && $ow == 0){
                            print "Sample branch lengths for CLUES \($popi_d_names[$l]\: chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\) start...\n";
                            $qout = "$Relate\/scripts\/SampleBranchLengths\/SampleBranchLengths.sh -i $o_path_clues[$l]\/popsize\/$ran\_popsize_chr$clues_chr[$i] -m $mrate $ns\--first_bp $clues_fbp[$i] --last_bp $clues_lbp[$i] $coal\--format $cf -o $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\\n";
                        }
                        #print "debug: $o_path_clues[$l]\/chr$clues_chr[$i]\_$clues_fbp[$i].txt\n";
                        unless (-e "$o_path_clues[$l]\/chr$clues_chr[$i]\_$clues_fbp[$i].txt" && -e "$o_path_clues[$l]\/chr$clues_chr[$i]\_$clues_fbp[$i]\_freq.txt" && $ow == 0){
                            my $haplo;
                            if ($hap == 1){
                                $haplo = "-hap";
                            }
                            else {
                                $haplo = "";
                            }
                            $qout .= "perl extract_derived_file.pl -i $o_path\/inputs\/$ran\_input_chr$clues_chr[$i] -list $o_path_clues[$l]\/sample_list.txt -o $o_path_clues[$l] -pos $clues_fbp[$i] $haplo\\n";
                        }
                        unless (-e "$o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\_times.txt" && $ow == 0){
                            $qout .= "python $Relate\/clues\/RelateToCLUES.py --RelateSamples $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].newick --DerivedFile $o_path_clues[$l]\/chr$clues_chr[$i]\_$clues_fbp[$i].txt --out $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\\n";
                        }
                        print BASH "qsub \.\/qsub_files\/$ran\_relate_CLUES_$popi_d_names[$l]\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\.q\n";
                        &pbs_setting("$exc\-cj_quiet $conda\-cj_qname $z\_relate_CLUES_$popi_d_names[$l]\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] -cj_sn $ran -cj_qout . $qout");
                        if ($z == 1 && $exc){
                            &status($ran);
                        }
                        $qout = "";
                        my $pf;
                        if ($exc){
							open(FREQ, "<$o_path_clues[$l]\/chr$clues_chr[$i]\_$clues_fbp[$i]\_freq.txt") || die "Cannot open $o_path_clues[$l]\/chr$clues_chr[$i]\_$clues_fbp[$i]\_freq.txt: $!\n";
							$pf = <FREQ>;
							#print "debug1: $pf\n";
							chomp($pf);
							close(FREQ);
							$pf = "--popFreq $pf ";
							#print "debug2: $pf\n";
                        }
                        else {
							$pf = "--popFreq 1 ";
                        }
                        unless (-e "$o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_post.txt" && $ow == 0){
                            #$qout .= "python $Relate\/clues\/inference.py --times $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --out $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $tco$pf$dom$coal$cbins[$i]\\n";
                            push(@plot_py, "python $Relate\/clues\/inference.py --times $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\_times.txt --out $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $nat$tco$df$pf$dom$coal$cbins[$i]\\n");
                        }
                        elsif ($replot =~ /clues|all/i){
                            #$qout .= "python $Relate\/clues\/inference.py --times $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --out $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $tco$pf$dom$coal$cbins[$i]\\n";
                            push(@plot_py, "python $Relate\/clues\/inference.py --times $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\_times.txt --out $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $nat$tco$df$pf$dom$coal$cbins[$i]\\n");
                        }
                        #print "debug: plot\n";
                        if ($clues_plot == 1){
                            unless (-e "$o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i].png" && $ow == 0){
                                #$qout .= "Rscript clues_plot.R -f $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $cps[$i]\\n";
                                #print "debug: plot2\n";
                                push(@plot_r, "python $Relate\/clues\/plot_traj_clues2.py --figure $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --post $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_post.txt --freqs $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_freqs.txt\\n");
                            }
                            elsif ($replot =~ /clues|all/i){
                                #$qout .= "Rscript clues_plot.R -f $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] $cps[$i]\\n";
                                push(@plot_r, "python $Relate\/clues\/plot_traj_clues2.py --figure $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i] --post $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_post.txt --freqs $o_path_clues[$l]\/$ran\_$z\_CLUES$popi_out_name\_chr$clues_chr[$i]\_$clues_fbp[$i]\-$clues_lbp[$i]\_freqs.txt\\n");
                            }
                        }
                    }
                }
                close(BASH);
            }
        }
        if ($exc){
            &status($ran);
        }
        print "Sample branch lengths for CLUES\/TreeViewSamples finished. Repeat times: $z\n";
    }
    for my $k (0..$#plot_py){
        &pbs_setting("$exc\-cj_quiet $conda\-cj_qname py_plot_file_$k -cj_sn $ran -cj_qout . $plot_py[$k]");
    }
    if ($exc){
        &status($ran);
    }
	if ($clues_plot == 1){
		for my $k (0..$#plot_r){
			&pbs_setting("$exc\-cj_quiet $conda\-cj_qname plot_R_$k -cj_sn $ran -cj_qout . $plot_r[$k]");
		}
		if ($exc){
			&status($ran);
		}
	}
}

sub subset_sample {
	print "Extracting samples of interests...\n";
	my $ran = shift; my $conda = shift; my $o_path = shift; my $pop = shift, my $popi = shift; my $exc = shift; my $ow = shift; 
	my @r_inputs = @{$_[-1]};
	my @new_r_inputs;
	my $popi_out_name;
	$popi_out_name = $popi;
	$popi_out_name =~ s/--pop_of_interest /_/;
	$popi_out_name =~ s/[\s\,]//g;
	my $cnt = 1;
	#print "debug: $o_path\n";
	unless (-d "$o_path\/popsize"){
		system("mkdir $o_path\/popsize");
	}
    foreach my $i (0..$#r_inputs){
    	my $gz; my $qout;
    	if (-e "$r_inputs[$i].anc.gz" && -e "$r_inputs[$i].mut.gz"){
    		$gz = ".gz";
    	}
    	else {
    		$gz = "";
    	}
    	my @tmp = split(/\//, $r_inputs[$i]);
    	my $file_name = $tmp[-1];
    	$file_name =~ s/input/popsize/;
		unless ((-e "$o_path\/popsize\/$file_name.anc" || -e "$o_path\/popsize\/$file_name.anc.gz") && $ow == 0){
    		$qout = "$Relate\/bin\/RelateExtract --mode SubTreesForSubpopulation --anc $r_inputs[$i].anc$gz --mut $r_inputs[$i].mut$gz $pop$popi\-o $o_path\/popsize\/$file_name\\n";
    	}
    	unless (-e "$o_path\/popsize\/$file_name.dist" && $ow == 0){
    		if (-e "$o_path\/popsize\/$file_name.mut.gz"){
    			$gz = ".gz";
    		}
    		else {
    			$gz = "";
    		}
    		$qout .= "$Relate\/bin\/RelateExtract --mode ExtractDistFromMut --mut $o_path\/popsize\/$file_name.mut$gz -o $o_path\/popsize\/$file_name\\n";
    	}
    	push(@new_r_inputs, "$o_path\/popsize\/$file_name");
    	&pbs_setting("$exc\-cj_quiet -cj_qname relate_subset$popi_out_name\_$cnt -cj_sn $ran -cj_qout . $qout");
    	$cnt++;
    }
	if ($exc){
		&status($ran);
	}
	my @poplabs = <$o_path\/popsize\/*.poplabels>;
	if (@poplabs){
	    my $tmp = $poplabs[0];
	    open(LABS, "<$poplabs[0]") || die "CLUES2: Cannot open $poplabs[0]: $!\n";
	    my @plist = <LABS>;
	    close(LABS);
	    chomp(@plist);
	    shift(@plist);
	    open(OUTP, ">$o_path\/sample_list.txt") || die "CLUES2: Cannot write sample_list.txt\n";
	    foreach my $i (0..$#plist){
	        my @tmp_eles = split(/\s+|\t/, $plist[$i]);
	        print OUTP "$tmp_eles[0]\n";
	    }
	    close(OUTP);
	}
    return @new_r_inputs;
}

sub subset_inputs {
	print "Extracting input samples of interests...\n";
	my $ran = shift; my $conda = shift; my $o_path = shift; my $pop = shift, my $popi = shift; my $exc = shift; my $ow = shift; 
	my @r_inputs = @{$_[-1]}; #input .sample and .haps
	my @new_r_inputs;
	my $popi_out_name;
	$popi_out_name = $popi;
	$popi_out_name =~ s/--pop_of_interest /_/;
	$popi_out_name =~ s/[\s\,]//g;
	$pop =~ s/--poplabels //;
	$popi =~ s/--pop_of_interest //;
	my @popis = split(/\,/, $popi);
	open(S_LIST, "<$pop") || die "Cannot open $pop: $!\n";
	my @list = <S_LIST>;
	close(S_LIST);
	chomp(@list);
	shift(@list);
	my @new_list;
	foreach my $k (0..$#list){
		my @tmp = split(/\s+|\t/, $list[$k]);
		my $remove = 1;
		foreach my $l (0..$#popis){
			$tmp[1] =~ s/\s+//g;
			$popis[$l] =~ s/\s+//g;
			if ($tmp[1] eq $popis[$l]){
				$remove = 0;
			}
		}
		if ($remove == 1){
			push(@new_list, $tmp[0]);
		}
	}
	my $cnt = 1;
	unless (-d "$o_path\/inputs"){
		system("mkdir $o_path\/inputs");
	}
	open(OUT, ">$o_path\/inputs\/removed_list.txt") || die "Cannot write $o_path\/inputs\/removed_list.txt: $!\n";
	print OUT join("\n", @new_list), "\n";
	close(OUT);
    foreach my $i (0..$#r_inputs){
    	my $gz; my $qout;
    	if (-e "$r_inputs[$i].haps.gz" && -e "$r_inputs[$i].sample.gz"){
    		$gz = ".gz";
    	}
    	else {
    		$gz = "";
    	}
    	my @tmp = split(/\//, $r_inputs[$i]);
    	my $file_name = $tmp[-1];
		unless ((-e "$o_path\/inputs\/$file_name.haps" || -e "$o_path\/inputs\/$file_name.haps.gz") && $ow == 0){
			$qout = "$Relate\/bin\/RelateFileFormats --mode RemoveSamples --haps $r_inputs[$i].haps$gz --sample $r_inputs[$i].sample$gz -i $o_path\/inputs\/removed_list.txt -o $o_path\/inputs\/$file_name\\n";
    	}
    	&pbs_setting("$exc\-cj_quiet -cj_qname relate_subset_inputs$popi_out_name\_$cnt -cj_sn $ran -cj_qout . $qout");
    	$cnt++;
    }
	if ($exc){
		&status($ran);
	}
	#if (-e "$o_path\/inputs\/removed_list.txt"){
		#system("rm $o_path\/inputs\/removed_list.txt");
	#}
}

sub cal_eps {
	print "Estimate population size start...\n";
	my $rp = shift; my $replot = shift; my $ran = shift; my $conda = shift; my $o_path = shift; my $mrate = shift; my $pop = shift; my $popi = shift; my $year = shift; my $td = shift; my $exc = shift; my $ow = shift; my $bin = shift;
	my @r_inputs = @{$_[-1]};
	my $input = $r_inputs[0];
	$input =~ s/_chr1$//;
	$input =~ s/_1$//;
	my $popi_out_name;
	if ($popi){
		$popi_out_name = $popi;
		$popi_out_name =~ s/--pop_of_interest /_/;
		$popi_out_name =~ s/[\s\,]//g;
	}
	my $bin_out; my $qout;
    my $last_chr = scalar(@r_inputs);
    unless ($rp) {
    	$rp = 0;
    }
    foreach my $i (0..$rp){
    	my $mcmc; my $mcmc_n;
    	$qout = "";
    	if ($i > 0){
    		$mcmc = "\/mcmc$i";
    		$mcmc_n = "_mcmc$i";
    	}
    	unless (-d "$o_path\/popsize$mcmc"){
    		system("mkdir $o_path\/popsize$mcmc");
    	}
    	if ($bin != -1){
    		$bin_out = "--bins $bin ";
    	}
    	unless (-e "$o_path\/popsize\/$ran$mcmc_n\_popsize.coal" && $ow == 0){
    		$qout .= "$Relate\/scripts\/EstimatePopulationSize\/EstimatePopulationSize.sh --threads 8 -i $input -m $mrate --first_chr 1 --last_chr $last_chr $pop$popi$year$td$bin_out\-o $o_path\/popsize$mcmc\/$ran\_popsize\\n";
    	}
    	unless (-e "$o_path\/popsize\/$ran\_popsize.pairwise.pdf" && ($ow == 0 || $replot eq "all")){
			my $plot_year = $year;
			$plot_year =~ s/--years_per_gen|\s+//g;
			if ($i == 0){
				$qout .= "Rscript plot_population_size_new.R $o_path\/popsize $plot_year $ran\\n";
			}
    	}
    	if ($i > 0){
    		unless (-e "$o_path\/popsize\/$ran\_mcmc$i\_popsize.pairwise.coal" && -e "$o_path\/popsize\/$ran\_mcmc$i\_popsize.pairwise.bin" && $ow == 0){
    			$qout .= "mv $o_path\/popsize$mcmc\/$ran\_popsize.coal $o_path\/popsize\/$ran\_mcmc$i\_popsize.coal\\n";
    			$qout .= "mv $o_path\/popsize$mcmc\/$ran\_popsize.pairwise.coal $o_path\/popsize\/$ran\_mcmc$i\_popsize.pairwise.coal\\n";
    			$qout .= "mv $o_path\/popsize$mcmc\/$ran\_popsize.pairwise.bin $o_path\/popsize\/$ran\_mcmc$i\_popsize.pairwise.bin\\n";
    			$qout .= "mv $o_path\/popsize$mcmc\/$ran\_popsize_avg.rate $o_path\/popsize\/$ran\_mcmc$i\_popsize_avg.rate\\n";
    			$qout .= "rm -r $o_path\/popsize$mcmc\\n";
    		}
    	}
    	&pbs_setting("$exc$conda\-cj_quiet -cj_ppn 8 -cj_qname relate_EPS$popi_out_name\_$i -cj_sn $ran -cj_qout . $qout");
    }
	if ($exc){
    	&status($ran);
	}
	foreach my $i (0..$#r_inputs){
		my $num = $i+1;
		$r_inputs[$i] = "$o_path\/popsize\/$ran\_popsize$popi_out_name\_chr$num";
	}
	return @r_inputs;
} #calculate estimate population size

sub get_first_ele_Relate {
	my $vcf = shift;
	my @content; my $first_ele;
	if (-e $vcf){}
	else {
		return 2;
	}
	if ($vcf =~ /\.vcf\.gz/){
		@content = `gzip \-cd $vcf \| head \-n 1000`;
	}
	elsif ($vcf =~ /\.vcf$/){
		@content = `head -n 1000 $vcf`;
	}
	else {
		return 2;
	}
	foreach (@content){
		if ($_ =~ /^\#/){
			next;	
		}
		my @line = split(/\t/, $_);
		$first_ele = $line[0];
		last;
	}
	return $first_ele;
}
sub check_path_Relate {
	my $path = shift;
	my $dir = getcwd;
	if ($path =~ /\//){
		my @path_eles = split(/\//, $path);
		my @dir_eles = split(/\//, $dir);
		if ($path_eles[0] eq "."){
			$path =~ s/^.//;
			$path = "$dir$path";
		}
		else {
			my $dot_cnt = -1;
			foreach (@path_eles){
				if ($_ eq ".."){
					$dot_cnt++;
				}
			}
			if ($dot_cnt == -1){
                if ($path =~ /^$dir/){}
				elsif (-d "$dir\/$path" || -e "$dir\/$path"){
					$path = "$dir\/$path";
				}
				else {
                    my @tmp = split(/\//, $path);
                    pop(@tmp);
                    my $tmp_path = join("\/", @tmp);
                    if (-d "$dir\/$tmp_path"){
                        $path = "$dir\/$path";
                    }
				}
			}
			else {
				for (my $i=0; $i<=$dot_cnt; $i++){
					shift(@path_eles);
					pop(@dir_eles);
					$path = join("\/", @dir_eles)."\/".join("\/", @path_eles);
				}
			}
		}
	}
	else {
		if ($path eq "."){
			$path = "$dir";
		}
		else {
			$path = "$dir\/$path";
		}
	}
	return($path);
} #relative path to absolute path
sub modify_poplabels {
	my $pop = shift;
	my @pop_eles = split(/\s/, $pop);
	$pop = $pop_eles[1];
	if (-e $pop){
		open(POP,"<$pop") || die "Cannot open $pop: $!\n";
		my @lines = <POP>;
		close(POP);
		chomp(@lines);
		my $head = shift(@lines);
		open(OUT, ">$pop") || die "Cannot write $pop: $!\n";
		print OUT "$head\n";
		foreach (@lines){
			print OUT "$_\n$_\n";
		}
		close(OUT);
		return 1;
	}
	else {
		return 2;
	}
}
sub usage {
	print BOLD "Usage: perl Relate_6_direct_5.pl -vcf VCF_FILE -pop POPULATION_LABEL_FILE -map RECOMB_MAP_FILE -al ANCESTOR_ID_LIST [-am ANC\/MUT_FOLDER_PATH] [-hap] [-o OUTPUT_PATH] [-mask MASK_FILE] [-bins LOWER,UPPER,STEPSIZE] [-rm REMOVE_SAMPLE_ID_FILE] [-pre PREFIX] [-rr] [-coal COAL_FILE] [-dps] [-clues] [-tvs] [-bp CHR\:POS-POS] [-ns INT] [-tco INT] [-pf FLOAT] [-d FLOAT] [-cp COLOR_PALETTE] [-m VALUE] [-n VALUE] [-spl VALUE] [-popi POP_NAMES] [-year VALUE] [-rp all\|clues] [-cb FILE] [-rc INT] [-epsrp INT] [-rpp] [-nat] [--force] [-ow] [-sn SERIAL_NUMBER] [-mem MEMORY_IN_GB] [-exc] [-h]\n\n", RESET;
	print "If -am is set, -vcf, -map and -al are not required.\n";
	print "For -popi, multiple populations as an group could be indicated by \[population_1,population_2\]. Multiple independent runs can be indicated by comma as population_1,population2.\nYou can combine these two functions as [population_1,population_2],population_3\n";
	print "The first run will be population_1\+population_2, and the second run will be population_3.\n";
	return;
}
