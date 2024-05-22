#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
if ($#ARGV == -1){
	@ARGV = <*.assoc>;
	if ($#ARGV == -1){
		print "No association file found.\n";
		exit;
	}
}
my $maf;
if ($ARGV[-1] !~ /[^0-9\.]/){
	$maf = $ARGV[-1];
	pop(@ARGV);
}
if (-d $ARGV[0]){
	$ARGV[0] =~ s/\/$//;
    @ARGV = <$ARGV[0]\/*.assoc>;
	if ($#ARGV == -1){
		print "No association file found.\n";
		exit;
	}
}

foreach (@ARGV){
    if (-e $_){
        open(FILE, "<$_") || die "Cannot find $_\:$!\n";
        open(OUT, ">$_\.tmp") || die "Cannot find $_\.tmp:$!\n";
        my $modified = 1; my $af;
        while (my $line = <FILE>){
            chomp($line);
            if ($line =~ /p_lrt/){
                $line =~ s/p_lrt/pval/;
                $line =~ s/ps/pos/;
                $modified = 0;
                if ($maf =~ /[0-9]/){
                	my @tmp = split(/\t/, $line);;
                	foreach my $i (0..$#tmp){
                		if ($tmp[$i] =~ /\baf\b/){
                			$af = $i;
                		}
                	}
                }
            }
            my @eles = split(/\t/, $line);
            if ($line !~ /^chr\b/){
                if ($modified == 0){
                    if ($eles[0] eq "-9"){
                        my @split_eles1 = split(/\:/, $eles[1]);
                        if (scalar(@split_eles1) != 3){
                            die "Cannot modify the $_ file.\n";
                        }
                        else {
                            $eles[0] = $split_eles1[1];
                            $eles[2] = $split_eles1[2];
                            $eles[1] = $split_eles1[0];
                        }
                    }
                    if ($eles[0] !~ /[a-z]/i){
                        if ($eles[0] < 10){
                            $eles[0] = "Chr0".$eles[0];
                        }
                        else {
                            $eles[0] = "Chr".$eles[0];
                        }
                    }
                }
                if ($modified == 1 && $maf !~ /[0-9]/){
                    close(OUT);
                    close(FILE);
                    system("rm $_.tmp");
                    exit;
                }
                if ($maf =~ /[0-9]/){
                	if ($eles[$af] < $maf || $eles[$af] > (1-$maf)){
                		next;
                	}
                }
                $line = join("\t", @eles);
            }
            print OUT "$line\n";
        }
        close(FILE);
        close(OUT);
        my $out = $_;
        if ($maf =~ /[0-9]/){
            $out =~ s/.assoc$/.maf_$maf.assoc/;
        }
        system("mv $_.tmp $out");
    }
}
