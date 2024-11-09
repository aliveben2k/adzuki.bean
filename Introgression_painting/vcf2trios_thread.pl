#!/usr/bin/perl

use Term::ANSIColor qw(:constants);
use threads;
use threads::shared;
#use Time::HiRes qw(time);
my $home = (getpwuid $>)[7];
if (-e "$home\/softwares\/qsub_subroutine.pl"){
	require "$home\/softwares\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;
}

print "Input command line:\n";
print "perl vcf2trios_thread\.pl @ARGV\n\n";

my $path; my @vcfs; my $list; my $out; my $chr; my $ran; my $skm = 0;
my $thread_in = 1; my $reg; my $sep = 0; my $path_l;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-vcf"){
		if (-d $ARGV[$i+1]){
            $path = $ARGV[$i+1];
            if ($path =~ /\/$/){
                $path =~ s/\/$//;
            }
            @vcfs = <$path\/*.vcf.gz>;
            unless (@vcfs){
                @vcfs = <$path\/*.vcf>;
                unless(@vcfs){
                    die "Cannot find vcfs.\n";
                }
            }
		}
        elsif (-e $ARGV[$i+1]){
            @vcfs = $ARGV[$i+1];
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
	if ($ARGV[$i] eq "\-list"){ #trio_list
        if (-e $ARGV[$i+1]){
            $list = $ARGV[$i+1];
        }
        if ($list =~ /\//){
			my @tmps = split(/\//, $list);
			pop(@tmps);
			$path_l = join("\/", @tmps);
        }
		else {
			$path_l = ".";
		}
		$path_l = &check_path($path);
	}
	if ($ARGV[$i] eq "\-sn"){
        $ran = "$ARGV[$i+1]";
	}
	if ($ARGV[$i] eq "\-t"){ #only extract the target chromosome (for multi-contig vcf)
		$sep = 1;
		if (length($ARGV[$i+1]) > 0){
			$chr = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-skm"){ #skip monosite check
		$skm = 1;
	}
	if ($ARGV[$i] eq "\-n"){ #threads
        $thread_in = "$ARGV[$i+1]";
        if ($thread_in =~ /[^0-9]/){
        	die "-n parameter should be an integer number.\n";
        }
	}
	#for test only
	if ($ARGV[$i] eq "\-reg"){
		$reg = $ARGV[$i+1];
	}
	#end of test
}

#print "debug: @vcfs\n";

unless(@vcfs){
	die "-vcf argument is required.\n";
}
my @chrs = &chr_name($vcfs[0]);
my @chr_lengths = &chr_lengths($vcfs[0]);
my $chr_names = join(" ", @chrs);
if ($chr){
	if ($chr_names !~ /\b$chr\b/){
		$chr = "";
		die "-t: Cannot find the target chromosome\/contig.\n";
	}
}
my $info_out;
if ($ran){
	$info_out = "$path_l\/$ran\_genome_info.txt";
}
else {
	$info_out = "$path_l\/genome_info.txt";
}
open (INFO, ">$info_out") || die "Cannot write the genome info to $info_out: $!\n";
print INFO "Chr\tLength\n";
foreach my $i (0..$#chrs){
	print INFO "$chrs[$i]\t$chr_lengths[$i]\n";
}
close(INFO);

#print "debug: @chrs\n";

#my $begin_time = time();
my $thread; my @check;
foreach my $cnt (0..$#vcfs){
	$thread = $thread_in;
	if ($vcfs[$cnt] !~ /gz$/){
		system("bgzip $vcfs[$cnt]");
		$vcfs[$cnt] = $vcfs[$cnt]."\.gz";
	}
    my $out = $vcfs[$cnt];
    $out =~ s/vcf\.gz$/$ran.\Qtrios\E.gz/;
    if (-e $out){
    	system("rm $out");
    }
    unless (-e "$vcfs[$cnt].tbi"){
    	system("tabix $vcfs[$cnt]");
    }
    if ($sep == 1 && length($chr) > 0){
    	$out =~ s/trios.gz$/$chr.$ran.\Qtrios\E.gz/;
    }
    if (scalar(@vcfs) == 1){
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			push(@check, $_);
    		}
    	}
    }
    else {
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			$chr = $_;
                push(@check, $_);
    		}
    	}   	
    }
    my $line = `zcat $vcfs[0] \| head -n 10000 \| grep \"\#CHROM"`;
    @vcf_samples = split(/\t/, $line);
    chomp(@vcf_samples);
    my @header_for_print; my @lists;
    if ($list){
    	open(LIST, "<$list") || die "Cannot open $list: $!\n";
    	my @list_tmp = <LIST>;
    	chomp(@list_tmp);
    	foreach (@list_tmp){
    		my @line_eles = split(/\t|\s+/, $_);
    		shift(@line_eles);
    		push(@lists, @line_eles);
    	}
    	#print "debug: @lists\n";
        foreach (@lists){
            foreach my $k (0..$#vcf_samples){
                if ($_ eq $vcf_samples[$k]){
                	push(@header_for_print, $vcf_samples[$k]);
                }
            }
        }
    }
    splice(@vcf_samples, 0, 9);
    unless ($list){
    	foreach (@vcf_samples){
    		push(@header_for_print, $_);
    	}
    }
    unless ($sep == 1 && scalar(@check) > 1){
        open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
        print OUT "Allele\t";
        print OUT join("\t", @header_for_print), "\n";
    }
    #print "debug: @check\n";
    foreach my $i (0..$#check){
        if ($sep == 1 && scalar(@check) > 1){
    	    $out =~ s/trios.gz$/$check[$i].$ran\Qtrios\E.gz/;
            open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
            print OUT "Allele\t";
        	print OUT join("\t", @header_for_print), "\n";
        }
        my $chr_length; #get chromosome length
        foreach my $j (0..$#chrs){
            if ($check[$i] eq $chrs[$j]){
                my @vcf_header;
                if ($vcfs[$cnt] =~ /gz$/){
                    @vcf_header = `zcat $vcfs[$cnt] \| head -n 10000 \| grep \"\#\#contig\=\"`;
                }
                else {
                    @vcf_header = `cat $vcfs[$cnt] \| head -n 10000 \| grep \"\#\#contig\=\"`;
                }
                chomp(@vcf_header);
                foreach my $k (0..$#vcf_header){
                    my @tmp = split(/\,/, $vcf_header[$k]);
                    if ($tmp[0] =~ /$check[$i]\b/){
                        my @tmp2 = split(/\=/, $tmp[1]);
                        chomp(@tmp2);
                        $chr_length = $tmp2[1];
                    }
                }
            }
        }
        #for test only
        if ($reg){
    	    $chr_lengths[$cnt] = $reg;
        }
        #end of test
	    my $division; my $results;
	    #print "debug: $thread\n";
	    if ($thread > 1){
		    if ($vcfs[$cnt] =~ /\.gz$/){
    		    unless (-e "$vcfs[$cnt].tbi"){
    			    system("tabix $vcfs[$cnt]");
    		    }
		    }
		    else {
    		    system("bgzip -\@ $thread $vcfs[$cnt]");
    		    system("tabix $vcfs[$cnt].gz");
    		    $vcf = "$vcfs[$cnt].gz";
		    }
    	    $division = $chr_length / $thread;
    	    if ($division - int($division) > 0.5){
    		    $thread = $thread - 1;
    	    }
    	    my @thr;
    	    foreach my $d (1..$thread){
    		    my $region_start; my $region_end;
    		    if ($d < $thread){
    			    if ($d == 1){
    				    $region_start = 1;
    			    }
    			    else {
    				    $region_start = int($chr_length/$thread*($d-1))+1;
    			    }
    			    $region_end = int($chr_length/$thread*($d));
    		    }
    		    else {
    			    $region_start = int($chr_length/$thread*($d-1))+1;
    			    $region_end = $chr_length;
    		    }
    		    #print "debug: $vcfs[$cnt] $sep $chr $region_start $region_end $list\n";
    		    $thr[$d] = threads->create('v2b', ($vcfs[$cnt],$sep,$check[$i],$region_start,$region_end,$list,$skm));
    	    }
    	    foreach my $d (1..$thread){
    		    $results .= $thr[$d]->join();
    	    }
    	    undef(@thr);
	    }
	    else {
		    my $region_start = 0; my $region_end = 0;
		    $results = &v2b($vcfs[$cnt],$sep,$chr,$region_start,$region_end,$list,$skm);
        }
        print OUT $results;
        if ($sep == 1 && scalar(@check) > 1){
            close(OUT);
        }
    }
    unless ($sep == 1 && scalar(@check) > 1){
        close(OUT);
    }
}
#my $end_time = time();
#printf("%.2f\n", $end_time - $begin_time);
print "Done.\n";

sub v2b {
    my $vcf = $_[0]; my $sep = $_[1]; my $chr = $_[2]; my $region_start = $_[3]; my $region_end = $_[4]; my $list = $_[5]; my $skm = $_[6];
    my @lists;
    if ($list){
    	open(LIST, "<$list") || die "Cannot open $list: $!\n";
    	my @list_tmp = <LIST>;
    	chomp(@list_tmp);
    	foreach (@list_tmp){
    		my @line_eles = split(/\t|\s+/, $_);
    		shift(@line_eles);
    		push(@lists, @line_eles);
    	}
    }
    #print "debug: @lists\n";
    if ($vcf =~ /\.gz$/){
    	unless (-e "$vcf.tbi"){
    		system("tabix $vcf");
    	}
    	if ($region_start != 0){
        	open(INPUT, "-|", "tabix -h $vcf $chr\:$region_start\-$region_end") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    	}
    	else {
    		open(INPUT, "-|", "bgzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    	}
    }
    else {
    	die BOLD "$vcf must be compressed first.", RESET, "\n";
    }
    my $out;
    my @samples; my @index; my $sep_check = 0;
    my @header;
    while (my $line = <INPUT>){
        chomp($line);
        $line =~ s/[\x0A\x0D]//g;
        if ($line =~ /^\#/){
        	#print "debug: $line\n";
            if ($line =~ /^\#CHROM/){
                @samples = split(/\t/, $line);
                if ($list){
                	foreach (@lists){
                		foreach my $k (0..$#samples){
                			if ($_ eq $samples[$k]){
                				push(@index, $k);
                				push(@header, $samples[$k]);
                			}
                		}
                	}
                }
                splice(@samples, 0, 9);
                unless ($list){
                	@header = @samples;
                }
            }
            next;
        }
        else { #skip lines with un-targeted chromosomes
        	if ($sep == 1 && length($chr) > 0){
        		if ($line !~ /^\b$chr\b/){
        			if ($sep_check == 1){
        				last;
        			}
        			next;
        		}
        		else {
        			$sep_check = 1;
        		}
        	}
        }
        my @eles = split(/\t/, $line);
        my @ele8 = split(/\:/, $eles[8]); #split FORMAT info
        my $ds;
        foreach my $i (1..$#ele8){ #determine the position of DS in FORMAT column
            if ($ele8[$i] eq "DS"){
                $ds = $i;
            }
        }
        unless ($ds){
            $ds = 0;
        }
        #skip only one allele type position
        if ($skm == 0){
        my @check;
        if ($list){
            @check = @index;
        }
        else {
            @check = @samples;
        }
        foreach my $j (9..$#eles){
            my $sum;
            my @format = split(/\:/, $eles[$j]);
            @allele = split(/\/|\|/, $format[0]);
            foreach (@allele){
                $sum += $_;
            }
            if ($list){
                foreach my $k (0..$#index){
                    if ($j == $index[$k]){
                        $check[$k] = $sum;
                    }
                }
            }
            else {
                $check[$j-9] = $sum;
            }
        }
        @check = do { my %seen; grep { !$seen{$_}++ } @check };
        if (scalar(@check) == 1){
            next;
        }
        }
        #end skipping process
		my $out_line; my @list_array;  my @list_array2;
		if ($list){
			foreach (@index){
				push(@list_array, 0);
				push(@list_array2, 0);
			}
		}
        foreach my $j (9..$#eles){
            my @format = split(/\:/, $eles[$j]);
            my @allele; my $allele1; my $allele2;
            @allele = split(/\/|\|/, $format[0]);
            foreach my $x (0..$#allele){
                if ($allele[$x] =~ /[1-8]/){
                    $allele[$x] = $allele[$x];
                }
                elsif ($allele[$x] == 0){
                    $allele[$x] = $allele[$x];
                }
                unless ($allele[$x] =~ /[0-9]/){
                    $allele[$x] = "NA";
                }
            }
			unless ($list){ #no list
                push(@list_array, "$allele[0]");
                push(@list_array2, "$allele[1]");
            }
            else { #with list
            	foreach my $l (0..$#index){
            		if ($j == $index[$l]){
            			$list_array[$l] = "$allele[0]";
            			$list_array2[$l] = "$allele[1]";
            		}
            	}
            }
        }
        #print "debug: @eles\n";
        $out .= "$eles[1]\_$eles[0]\_a1\t";
        $out .= join("\t", @list_array);
        $out .= "\n";
        $out .= "$eles[1]\_$eles[0]\_a2\t";
        $out .= join("\t", @list_array2);
        $out .= "\n";
    }
    close(INPUT);
	return $out;
}

sub usage {
    print "Usage: perl vcf2trios_thread.pl -vcf PATH [-list LIST_FILE] [-t CHR] [-skm] [-sn SN] [-n]\n";
}
