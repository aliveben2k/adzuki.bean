#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
my ($per, $R1, $R2);
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-p"){
		if ($ARGV[$i+1] !~ /[^0-9.]/){
			$per = $ARGV[$i+1];
			if ($per >= 1 || $per <= 0){
				die BOLD "\-p value is between 0 to 1.", RESET, "\n";
			}
		}
		else {
			die BOLD "\-p parameter is not correct.", RESET, "\n";
		}
	}
	if ($ARGV[$i] eq "\-i"){
		if ($ARGV[$i+1] !~ /^\-/){
			$R1 = $ARGV[$i+1];
			if (-e $R1){
				if ($R1 !~ /\.fq|\.fastq/){
					die BOLD "Input file\(s\) should be fastq file\(s\).", RESET, "\n";
				}
				if ($R1 =~ /_R1\.|_1\./i){
					$R2 = $ARGV[$i+2];
					if(-e $R2){}
					else {
						die BOLD "No paired file is found.", RESET, "\n";
					}
				}
			}
			else {
				die BOLD "No input file is found.", RESET, "\n";
			}
		}
	}
}

=cut
my $out;
if ($R1){
	if ($R1 !~ /gz$/){
		system("gzip $R1");
		$R1 = $R1.".gz";
		$out = 1;
	}
}
else {
	die BOLD "No input file is found.", RESET, "\n";
}
=cut
#print "Counting read number in $R1...\n";
#my $read_len = `echo \$\(cat $R1\|wc -l\)\/4\|bc`;
#print "Read number in $R1 is\: $read_len";
#my $k = int($read_len * $per);
#$k++ if (($read_len * $per) - int($read_len * $per) > 0.5);
#print "Downsampling read number to: $k\n";

open (INPUT, ">downsampling\.q") || die BOLD "Cannot write downsampling\.q: $!", RESET, "\n";
print INPUT "\#\!\/bin\/bash\n\#PBS \-l nodes\=1\:ppn\=1\n\#PBS \-o qsub_files\/out\/$ran\_gatk_04_$chr$k\.out \-j oe\n\n";
print INPUT "cd \$PBS_O_WORKDIR\n";
my $O1 = "sub_".$R1; my $O2;
my $out;
=cut
if ($R1){
	if ($R1 !~ /gz$/){
		print INPUT "gzip $R1\n";
		$R1 = $R1.".gz";
		$out = 1;
	}
}
else {
	die BOLD "No input file is found.", RESET, "\n";
}
$R1 =~ s/\.gz//;
=cut
if ($R2){
=cut
	if ($R2 !~ /gz$/){
		print INPUT "gzip $R2\n";
#		$R2 = $R2.".gz";		
	}
	$R2 =~ s/\.gz//;
=cut
	$O2 = "sub_".$R2;
	print "paste $R1 $R2 \| awk \'\{ printf\(\"\%s\"\,\$0\)\; n\+\+\; if\(n\%4\=\=0\) \{ printf\(\"\\n\"\)\;\} else \{ printf\(\"\\t\"\)\;\} \}\' \| awk \-v k\=$k \'BEGIN\{srand\(systime\(\) \+ PROCINFO\[\"pid\"\]\)\;\}\{s\=x\+\+\<k\?x1\:int\(rand\(\)\*x\)\;if\(s\<k\)R\[s\]\=\$0\}END\{for\(i in R\)print R\[i\]\}\' \| awk \-F\"\\t\" \'\{print \$1\"\\n\"\$3\"\\n\"\$5\"\\n\"\$7 \> \"$O1\"\;print \$2\"\\n\"\$4\"\\n\"\$6\"\\n\"\$8 \> \"$O2\"\}\'\n";
	#print INPUT "paste $R1 $R2 \| awk \'\{ printf\(\"\%s\"\,\$0\)\; n\+\+\; if\(n\%4\=\=0\) \{ printf\(\"\\n\"\)\;\} else \{ printf\(\"\\t\"\)\;\} \}\' \| awk \-v k\=$k \'BEGIN\{srand\(systime\(\) \+ PROCINFO\[\"pid\"\]\)\;\}\{s\=x\+\+\<k\?x1\:int\(rand\(\)\*x\)\;if\(s\<k\)R\[s\]\=\$0\}END\{for\(i in R\)print R\[i\]\}\' \| awk \-F\"\\t\" \'\{print \$1\"\\n\"\$3\"\\n\"\$5\"\\n\"\$7 \> \"$O1\"\;print \$2\"\\n\"\$4\"\\n\"\$6\"\\n\"\$8 \> \"$O2\"\}\'\n";
	#print INPUT "gzip $O1\n";
	#print INPUT "gzip $O2\n";
}
else{
	print "cat $R1 \| awk \'\{ printf\(\"\%s\"\,\$0\)\; n\+\+\; if\(n\%4\=\=0\) \{ printf\(\"\\n\"\)\;\} else \{ printf\(\"\\t\"\)\;\} \}\' \| awk \-v k\=$k \'BEGIN\{srand\(systime\(\) \+ PROCINFO\[\"pid\"\]\)\;\}\{s\=x\+\+\<k\?x1\:int\(rand\(\)\*x\)\;if\(s\<k\)R\[s\]\=\$0\}END\{for\(i in R\)print R\[i\]\}\' \| awk \-F\"\\t\" \'\{print \$1\"\\n\"\$2\"\\n\"\$3\"\\n\"\$4 \> \"$O1\"\}\'\n";
	#print INPUT "cat $R1 \| awk \'\{ printf\(\"\%s\"\,\$0\)\; n\+\+\; if\(n\%4\=\=0\) \{ printf\(\"\\n\"\)\;\} else \{ printf\(\"\\t\"\)\;\} \}\' \| awk \-v k\=$k \'BEGIN\{srand\(systime\(\) \+ PROCINFO\[\"pid\"\]\)\;\}\{s\=x\+\+\<k\?x1\:int\(rand\(\)\*x\)\;if\(s\<k\)R\[s\]\=\$0\}END\{for\(i in R\)print R\[i\]\}\' \| awk \-F\"\\t\" \'\{print \$1\"\\n\"\$2\"\\n\"\$3\"\\n\"\$4 \> \"$O1\"\}\'\n";
	#print INPUT "gzip $O1\n";
}
if ($out == 1){
	#print INPUT "rm $R1\n";
	if ($R2){
		#print INPUT "rm $R2\n";
	}
}
print "The output file is:\n$O1\.gz $O2\.gz\n";
close(INPUT);














