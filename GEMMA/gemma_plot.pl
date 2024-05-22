#!/usr/bin/perl

use Cwd qw(getcwd);
use FindBin;
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/softwares\/qsub_subroutine.pl"){
	require "$home\/softwares\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}

my $r_env = '/home/hpc/crlee/miniconda3/envs/R-4.1/bin';
my $gemma = 'gemma-0.98.5';

chomp(@ARGV);
print "The script is written by Ben Chien. Mar. 2022\n";
print "Input command line:\n";
print "perl gemma_plot\.pl @ARGV\n";

if ($#ARGV == -1){
	&usage;
	exit;
}

my $p_dir = $FindBin::Bin;
my $exc; my $ran; my $sn; my $local;
my $path; my $o_path; my $mem = 0;
my @files;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-exc"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "\-local"){
		$local = "-cj_local ";
	}
	if ($ARGV[$i] eq "-sn"){
		$ran = $ARGV[$i+1];
		$sn = 1;
	}
	if ($ARGV[$i] eq "\-h"){
		&usage;
		exit;
	}
	if ($ARGV[$i] eq "\-a"){
		if (-d $ARGV[$i+1]){
            if ($ARGV[$i+1] =~ /\/$/){
                $ARGV[$i+1] =~ s/\/$//;
            }
            $path = $ARGV[$i+1];
            @files = <$path\/*.assoc>;
		}
		elsif (-e $ARGV[$i+1]){
            if ($ARGV[$i+1] =~ /\//){
                my @tmp = split(/\//, $ARGV[$i+1]);
                pop(@tmp);
                $path = join("\/", @tmp);
            }
            else {
                $path = ".";
            }
            @files = $ARGV[$i+1];
		}
		$path = &check_path($path);
	}
	if ($ARGV[$i] eq "\-o"){
		$o_path = $ARGV[$i+1];
		if ($o_path =~ /\/$/){
			$o_path =~ s/\/$//;
		}
		unless (-d $o_path){
			my $return = `mkdir $o_path`;
			if ($return){
				die "Cannot make the directory $o_path: $return\n";
			}
		}
		$o_path = &check_path($o_path);
	}
	if ($ARGV[$i] eq "\-mem"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$mem = $ARGV[$i+1];
		}
	}
}
if ($mem == 0){
	$mem = "";
}
else {
	$mem = "-cj_mem $mem ";
}

unless ($o_path){
    $o_path = $path;
}

unless (@files){
    die "No file is loaded.\n";
}

if ($ran){}
else {
	$ran = &rnd_str(4, "A".."Z", 0..9);
}

foreach my $i (0..$#files){
    $files[$i] = &check_path($files[$i]);
    my $out = $files[$i];
    $out =~ s/$path/$o_path/;
    $out =~ s/.assoc$//;
    &pbs_setting("$exc$mem$local\-cj_env $r_env -cj_sn $ran -cj_qout . -cj_qname gemma_plot_$i Rscript $p_dir\/qqman_v2.R -f $files[$i] -o $out\\n");
}


sub usage {
	print BOLD "Usage: perl gemma.pl -a AN_ASSOC_FOLDER|AN_ASSOC_FILE [-o OUTPUT_FOLDER_PATH] [-mem MEMORY_USE_IN_GB] [-sn SERIAL_NUMBER] [-exc] [-h]\n\n", RESET;
	return;
}
sub rnd_str {
	join("", @_[map{rand@_} 1..shift]);
} #generate serial number
sub check_path {
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
				if (-e "$dir\/$path" || -d "$dir\/$path"){
					$path = "$dir\/$path";
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
	return $path;
} #relative path to absolute path
