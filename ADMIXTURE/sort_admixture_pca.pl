#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
print "The script is written by Ben Chien. Mar. 2022.\n";
print "Input command line:\n";
print "perl sort_admixture_pca\.pl @ARGV\n";

if ($#ARGV == -1){
	&usage;
	exit;
}

my ($vcf, $table, $f_name, $path, $list, $reorder, $mtx);
my $mix = -1;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-mtx"){
		$mtx = $ARGV[$i+1];
		if (-e $mtx){}
		else {
			&usage;
			exit;		
		}
	}
	if ($ARGV[$i] eq "\-vcf"){
		$vcf = $ARGV[$i+1];
		if (-e $vcf){}
		else {
			&usage;
			exit;		
		}
	}
	if ($ARGV[$i] eq "\-q"){
		$table = $ARGV[$i+1];
		if (-e $table){
			if ($table =~ /\//){
				my @tmp1 = split(/\//, $table);
				$f_name = pop(@tmp1);
				$path = join("\/", @tmp1);
			}
			else{
				$f_name = $table;
				$path = "\.";
			}
		}
		else {
			&usage;
			exit;		
		}
	}
	if ($ARGV[$i] eq "\-l"){ #info of the samples
		$list = $ARGV[$i+1];
		if (-e $list){}
		else {
			&usage;
			exit;		
		}
	}
	if ($ARGV[$i] eq "\-reorder"){
		$reorder = 1;
	}
	if ($ARGV[$i] eq "\-m"){
		$mix = $ARGV[$i+1];
		if ($mix > 0 && $mix <= 1){}
		else {
			&usage;
			exit;		
		}
	}
}
if ($vcf !~ /\w/ || $table !~ /\w/ || $list !~ /\w/){
	print "Some required arguments are missing.\n";
	&usage;
	exit;
}
my @ids = &sample_name($vcf);
open(QFILE, "<$table") || die "Cannot open $table: $!\n";
my @q_lines = <QFILE>;
chomp(@q_lines);
close(QFILE);

#deal with the Q file
open(SORTQ, ">$path\/sorted_$f_name") || die "Cannot write $path\/sorted_$f_name: $!\n";
my $header;
my %group_arrays;
my @groups;
foreach my $i (0..$#q_lines){
	my @values = split(/\t+|\s+/, $q_lines[$i]);
	unless ($header){ #print header
		foreach my $j (1..scalar(@values)){
			print SORTQ "\tG$j";
			push(@groups, "G$j");
			$group_arrays{$groups[$j]} = ();
		}
		print SORTQ "\n";
		$header = 1;
	}
	#get value index from big to small
	my @sort_idx = sort {@values[$b] <=> @values[$a]} 0..$#values;
	my @sort_vals = @values[@sort_idx];
	my $group = $sort_idx[0] + 1;
	$group = "G".$group;
	push(@values, $ids[$i]);
	my $value = join("\t", @values);
	push(@{$group_arrays{$group}}, $value);
}
my @id_info; my @reorders;
open(LIST, "<$list") || die "Cannot open $list: $!\n";
@id_info = <LIST>;
chomp(@id_info);
close(LIST);
if ($reorder == 1){
	@reorders = @id_info;
}

foreach my $i (1..scalar(@groups)){ # $key == group
	my $key = "G$i";
	my @group_values = @{$group_arrays{$key}};
	my @admixture_values;
	foreach (@group_values){ # @group_values == id values in groups 
		my @id_values = split(/\t/, $_); # @id_values == values of each sample in admixture results + id
		if ($reorder == 1){
			foreach my $k (0..$#reorders){
				$reorders[$k] =~ s/[\x0A\x0D]//g;
				if ($reorders[$k] =~ /\b$id_values[-1]\b/){
					if ($mix > -1){
						if ($id_values[$i-1] < $mix){
							$reorders[$k] = "admix\t$_";
						}
						else {
							$reorders[$k] = "$key\t$_";
						}
					}
					else {
						$reorders[$k] = "$key\t$_";
					}						
				}
			}
		}
		else {
			push(@admixture_values, $id_values[$i-1]);
		}
	}
	unless ($reorder == 1){
		my @sort_idx = sort {@admixture_values[$b] <=> @admixture_values[$a]} 0..$#admixture_values;
		my @sort_ad_vals = @admixture_values[@sort_idx];
		my @sort_gp_vals = @group_values[@sort_idx];
		foreach my $j (0..$#sort_gp_vals){
			if ($mix > -1){
				if ($sort_ad_vals[$j] < $mix){
					print SORTQ "admix\t$sort_gp_vals[$j]\n";
				}
				else {
					print SORTQ "$key\t$sort_gp_vals[$j]\n";
				}
			}
			else {
				print SORTQ "$key\t$sort_gp_vals[$j]\n";
			}
		}
	}
}
if ($reorder == 1){
	foreach (@reorders){
		if ($_ =~ /^G[0-9]\b|^G[0-9][0-9]\b|^admix\b/){
			print SORTQ "$_\n";
		}
	}
}
close(SORTQ);

#using sorted_q file to do admixture plot
open(QIN, "<$path\/sorted_$f_name") || die "Cannot open $path\/sorted_$f_name: $!\n";
open(ADOUT, ">$path\/plot_ad_$f_name") || die "Cannot write $path\/plot_ad_$f_name: $!\n";
my @qin_lines = <QIN>;
chomp(@qin_lines);
close(QIN);
my $qin_header = shift(@qin_lines);
print ADOUT "$qin_header\n";
foreach my $i (1..scalar(@qin_lines)){
	my @qin_vals = split(/\t/, $qin_lines[$i-1]);
	pop(@qin_vals);
	$qin_vals[0] = $i;
	print ADOUT join("\t", @qin_vals), "\n";
	
}
close(ADOUT);
my $number = scalar(@groups);
#unless (-e "$path\/admixture_K$number.pdf"){
	system("Rscript admixture_plot.R $path\/plot_ad_$f_name $path $number");
#}

unless ($mtx){
	my $v2t_file = $vcf;
	$v2t_file =~ s/.gz$//;
	$v2t_file =~ s/\.vcf$//;
	unless (-e "$v2t_file.1.tmp.txt.gz"){
		system("perl vcf2table_missingNA_large.pl $vcf $v2t_file.txt");
	}
	unless (-e "$v2t_file\.rda"){
    	my $path = $v2t_file;
    	if ($path =~ /\//){
        	my @tmp = split(/\//, $v2t_file);
        	pop(@tmp);
        	$path = join("\/", @tmp);
    	}
    	else {
       		$path = '.';
    	}
    	my @files = <$path\/*.tmp.txt.gz>;
    	unless (@files){
       		die "Cannot find the files for processing.\n";
    	}
    	foreach my $i (0..$#files){
        	my $outfile = $files[$i];
        	$outfile =~ s/\.txt\.gz$//;
        	system("Rscript Calculate_pairwise_dist_simple_large1.R $files[$i] $outfile");
    	}
    	system("Rscript Calculate_pairwise_dist_simple_large2.R $path $v2t_file");
    	if (-e "$v2t_file.rda"){
        	foreach (@files){
            	system("rm $_");
            	$rda_file = $_;
            	$rda_file =~ s/txt\.gz$/rda/;
            	system("rm $rda_file");
        	}
    	}
	}
	system("Rscript Tree_PCoA.R $v2t_file.rda $v2t_file");
	$mtx = "$v2t_file.mtx";
}

#merge group and type info
my $info_out = "$path\/info_$f_name";
$info_out =~ s/Q$/txt/i;
open(INFO, ">$info_out") || die "Cannot write $info_out: $!\n";
open(QIN, "<$path\/sorted_$f_name") || die "Cannot open $path\/sorted_$f_name: $!\n";
my @qin_lines = <QIN>;
chomp(@qin_lines);
close(QIN);
shift(@qin_lines);
foreach (@qin_lines){
	@line_eles = split(/\t+|\s+/, $_);
	$id_check = 0;
	foreach my $info (@id_info){
		$info =~ s/[\x0A\x0D]//g;
		if ($info =~ /\b$line_eles[-1]\b/){
			my @info_eles = split(/\t+|\s+/, $info);
			print INFO "$line_eles[-1]\t$line_eles[0]\t$info_eles[1]\n";
			$id_check = 1;
		}
	}
	if ($id_check == 0){
		print INFO "$line_eles[-1]\t$line_eles[0]\tNA\n";
	}
}
close(INFO);

#plot PCA result
my $mix_exist = 0;
if ($mix > -1){
	$mix_exist = 1;
}
system("Rscript pca_plot_v3.R $mtx $info_out $path $mix_exist");

sub usage {
	print BOLD "Usage: perl sort_admixture_pca.pl -vcf A_VCF_FILE -q A_Q_FILE -l A_SAMILE_INFO_FILE [-mtx MATRIX_FILE] [-reorder] [-m RATIO]\n\n", RESET;
	return;
}
sub sample_name {
	$time = scalar localtime();
	my $file = shift;
	my @content; my @line;
	if (-e $file){}
	else {
		return 2;
	}
	if ($file =~ /\.vcf\.gz/){
		@content = `gzip \-cd $file \| head \-n 10000`;
	}
	elsif ($file =~ /\.vcf$/){
		@content = `head -n 10000 $file`;
	}
	else {
		return 2;
	}
	foreach (@content){
		if ($_ =~ /\#CHROM/){
			@line = split(/\t/, $_);
			for (my $i=9; $i<=$#line; $i++){
                push(@samples, $line[$i]);
			}
			last;
		}
	}
	chomp(@samples);
	return (@samples);
} #get sample name
