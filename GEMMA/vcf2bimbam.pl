#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);
my $r_env = '/home/hpc/crlee/miniconda3/envs/ben/bin';
my $dir = getcwd;
use Time::HiRes qw(time);

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;
}

print "Input command line:\n";
print "perl vcf2bimbam\.pl @ARGV\n\n";

my $path; my @vcfs; my @lists; my $list; my $out; my $sep = 0; my $chr; my $ran; my $ioff = 0;
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
	if ($ARGV[$i] eq "\-list"){
        if (-e $ARGV[$i+1]){
            open (LIST, "<$ARGV[$i+1]") || die BOLD "Cannot open $ARGV[$i+1]: $!", RESET, "\n";
            @lists = <LIST>;
            chomp(@lists);
            shift(@lists);
            close(LIST);
            foreach my $j (0..$#lists){
            	my @sp_lines = split(/\t+|\s+/, $lists[$j]);
            	$lists[$j] = $sp_lines[0];
            }
            $list = $ARGV[$i+1];
        }
	}
	if ($ARGV[$i] eq "\-sn"){
        $ran = "$ARGV[$i+1]\.";
	}
	if ($ARGV[$i] eq "\-sep"){
		$sep = 1;
		if (length($ARGV[$i+1]) > 0){
			$chr = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-ioff"){
		$ioff = 1;
	}
}

unless(@vcfs){
	die "-vcf argument is required.\n";
}
my @chrs = &chr_name($vcfs[0]);
my $chr_names = join(" ", @chrs);
if ($chr_names !~ /\b$chr\b/ || length($chr) <= 0){
	$chr = "";
}
my $begin_time = time();
foreach my $vcf (@vcfs){
    if ($vcf =~ /\.gz$/){
        open(INPUT, "-|", "gzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    }
    else {
        open (INPUT, "<$vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    }
    $out = $vcf;
    $out =~ s/\.gz$//;
    $out =~ s/vcf$/$ran\Qbimbam\E.gz/;
#    if ($list && scalar(@vcfs) == 1){
#        $out = $list;
#        $out =~ s/txt$|list$/$ran\Qbimbam\E.gz/;
#    }
    if (-e $out){
    	system("rm $out");
    }
    if ($sep == 1 && length($chr) > 0){
    	$out =~ s/bimbam.gz$/$chr.$ran\Qbimbam\E.gz/;
    }
    elsif ($sep == 1){
    	my @check;
    	if ($vcf =~ /gz$/){
    		@check = `zcat $vcf \| head -n 10000`;
    	}
    	else {
    		@check = `cat $vcf \| head -n 10000`;
    	}
        foreach (@check){
            if ($_ !~ /^\#/){
                my @line_eles = split(/\t/, $_);
                $out =~ s/bimbam.gz$/$line_eles[0].$ran\Qbimbam\E.gz/;
                last;
            }
        }    	
    }
    open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
    my $cnt = 0; my @samples; my @index; my $sep_check = 0;
    while (my $line = <INPUT>){
        chomp($line);
        $line =~ s/[\x0A\x0D]//g;
        if ($line =~ /^\#/){
            if ($line =~ /^\#CHROM/){
                @samples = split(/\t/, $line);
                if ($list){
                	foreach (@lists){
                		foreach my $k (0..$#samples){
                			if ($_ eq $samples[$k]){
                				push(@index, $k);
                			}
                		}
                	}
                }
                splice(@samples, 0, 9);
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
            if ($ds == 0 || $ioff == 1){ #if there is no "DS" field or turn off imputed data, use sum of "GT" field
                foreach (@allele){
                    if ($_ =~ /[1-9]/){
                    	$_ = 1;
                        $sum += $_;
                    }
                    elsif ($_ == 0){
                        $sum += $_;
                    }
                }
            }
            else {
                if ($format[$ds] ne "."){
                    $sum = $format[$ds];
                }
                else {
                    foreach (@allele){
                        if ($_ =~ /[1-9]/){
                        	$_ = 1;
                            $sum += $_;
                        }
                        elsif ($_ == 0){
                            $sum += $_;
                        }
                    }
                }
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
        #end skipping process
#        print OUT "$eles[2]\:$eles[0]\:$eles[1]\, $eles[3]\, $eles[4]\, ";
		my $out_line; my @list_array;
		if ($list){
			foreach (@index){
				push(@list_array, 0);
			}
		}
        foreach my $j (9..$#eles){
            my @format = split(/\:/, $eles[$j]);
            my @allele; my $sum;
            @allele = split(/\/|\|/, $format[0]);
            if ($ds == 0 || $ioff == 1){ #if there is no "DS" field or turn off imputed data, use sum of "GT" field
                foreach (@allele){
                    if ($_ =~ /[1-9]/){
                    	$_ = 1;
                        $sum += $_;
                    }
                    elsif ($_ == 0){
                        $sum += $_;
                    }
                }
                unless ($sum =~ /[0-9]/){
                    $sum = "NA";
                }
            }
            else { #if there is a "DS" field, but has no value, use sum of "GT" field
                $sum = $format[$ds];
                if ($sum eq "."){
                    $sum = "";
                    foreach (@allele){
                        if ($_ =~ /[1-9]/){
                        	$_ = 1;
                            $sum += $_;
                        }
                        elsif ($_ == 0){
                            $sum += $_;
                        }
                    }
                    unless ($sum >= 0 && $sum <= 2){
                        $sum = "NA";
                    }
                }
            }
            if ($sum ne "NA" && $sum > 2){
                print "The input vcf is not bi-allele vcf, please filter first.\n";
                close(OUT);
                system("rm $out");
                close(INPUT);
                exit;
            }
			unless ($list){ #no list
                push(@list_array, $sum);
            }
            else { #with list
            	foreach my $l (0..$#index){
            		if ($j == $index[$l]){
            			$list_array[$l] = $sum;
            		}
            	}
            }
        }
        my $n_0 = 0;
        my $n_1 = 0;
        my $n_2 = 0;
        foreach my $i (0..$#list_array){
            if ($list_array[$i] >= 0 && $list_array[$i] <= 0.5){
                $n_0++;
            }
            if ($list_array[$i] > 0.5 && $list_array[$i] < 1.5){
                $n_1++;
            }
            if ($list_array[$i] >= 1.5 && $list_array[$i] <= 2){
                $n_2++;
            }
        }
        my $alleles = "$n_0 $n_1 $n_2";
        my $count = () = $alleles =~ /\b0\b/g;
        if ($count >= 2){
            next;
        }
        print OUT "$eles[2]\:$eles[0]\:$eles[1]\, $eles[3]\, $eles[4]\, ";
        print OUT join(", ", @list_array), "\n";
        $cnt++;
        print "\rProcessing $vcf: $cnt variants done.";
    }
    print "\n";
    close(INPUT);
    close(OUT);
}
my $end_time = time();
printf("%.2f\n", $end_time - $begin_time);
print "Done.\n";

sub usage {
    print "Usage: perl vcf2bimbam.pl -vcf PATH [-list LIST_FILE] [-sep CHR] [-ioff]\n";
}
sub chr_name {
	$time = scalar localtime();
	my $file = shift;
	my @content; my @line; my @id;
#	print "chr_name\n";
	if (-e $file){}
	else {
		return "no";
	}
	if ($file =~ /\.vcf\.gz/){
		@content = `gzip \-cd $file \| head \-n 10000`;
	}
	elsif ($file =~ /\.vcf$/){
		@content = `head -n 10000 $file`;
	}
	else {
		return "no";
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
            push(@id, $line[3]);
		}
	}
	return (@id);
} #get chromosome name
