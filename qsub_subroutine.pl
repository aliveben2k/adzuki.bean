#!/usr/bin/perl

#list fo the function: pbs_setting, status, get_1_eles, chr_name, chr_lengths, sample_name, check_path, rnd_str

=functions for pbs_setting
[-cj_local] [-cj_env PATH] [-cj_conda ENV_NAME] [-cj_node INT] [-cj_ppn INT] [-cj_mem INT] [-cj_qname JOB_NAME] [-cj_proj PROJECT_ID] [-cj_module MODULE] [-cj_qout PATH] [-cj_sn SN] [-cj_exc] [-cj_quiet] [-cj_queue QUEUE_NAME] [-cj_mail EMAIL_ADDRESS]
-cj_local	Run the script locally.
-cj_env		Environment path that need to be set in $PATH. Can be used for multiple times.
-cj_conda	Conda environment name. Only works for h71 and h81 servers.
-cj_node	Node that will be used in the job. Default: 1.
-cj_ppn		Core that will be used in the job. Default: 1.
-cj_time	walltime of the job. Default: system default.
-cj_mem		Momory that will be used in the job. Default: system default.
-cj_qname	Name of the job, if defined, the job name will be {SN}_{job_name}. Default: {SN}_cj
-cj_queue	Select queue for the job. Only works for Taiwania 1.
-cj_proj	ID of the project. Only works for Taiwania 1.
-cj_module	A module that need to be loaded. Can be used for multiple times.
-cj_qout	The output path where the job execution info should be stored.
-cj_sn		Serial number {SN} of the job. Default: 4-digit random characters.
-cj_mail	E-mail address. It will send the notice when the job starts and is done. Only works for Taiwania 1.
-cj_exc		Send the job for execution.
-cj_quiet	Supress the message.
=cut

=functions for rnd_str
input: none
main function: generate a four-character serial number
=cut

=functions for status_tw3
input: user_ID, job_file_name, quiet mode, partition
main function: send and track the jobs if the job number reaches the limitation
=cut

=functions for status_other
input: Serial_number($ran), user_ID, job_file_name, quiet mode, partition
main function: send and track the jobs if the job number reaches the limitation (>200)
=cut

=functions for status
input: Serial_number($ran), user_ID
main function: tracking the job status by the serial number given
=cut

=function for get_1_eles
input: array_of_vcf_files, contig_id
function: giving a series of vcf file paths, and reorder them by contig_id order in the vcf
=cut

=function for chr_name
input: a_vcf_file, prefix_name_of_the_contigs
function: get the contig names in the vcf, if the second argument is given, it will return a list only with the prefixed name
=cut

=function for chr_lengths
input: a_vcf_file, prefix_name_of_the_contigs
function: get the contig lengths in the vcf, if the second argument is given, it will return a list of lengths only with the prefixed name
=cut

=function for sample_name
input: a_vcf_file
function: get the sample IDs in the vcf
=cut

=function for check_path
input: a_path(to a file or a folder)
function: change relative path to absolute path
=cut

use Term::ANSIColor qw(:constants);
use Cwd qw(getcwd);

sub pbs_setting {
my $arg = shift;
my @args = split(/\s/, $arg);
my $nodes = 1; #set nodes used for qsub job
my $ppn = 1; #set ppn used for qsub job
my $mem = 0; #set memory used for qsub job; 0 means no memory defined.
my $home1 = (getpwuid $>)[7]; #get the $HOME path

my @server = `ip route get 1.2.3.4 \| awk \'\{print \$7\}\'`;
chomp(@server);
my $serv = -1;
my $scr_dir = '#PBS'; 
my $par = '-q'; 
my $snode; 
my $wall = '-l walltime='; 
my $jenv = '-V';
my $smail = '-M '; 
my $mtype = '-m abe'; 
my $jn = '-N ';  
my $sppn; 
my $pname = '-P ';
foreach (@server){
    	if ($_ =~ /<PBS_server_IP>/){ #PBS system
		$serv = 1;
    	}
    	elsif ($_ =~ /<Slurm_server_IP>/){ #Slurm system
		$serv = 3;
		$scr_dir = '#SBATCH';
		$par = '-p';
		$snode = '-N ';
		$wall = '-t ';
		$jenv = '--export=ALL';
		$smail = '--mail-user=';
		$mtype = '-–mail-type=ALL';
		$jn = '-J '; #job name
		$sppn = '--ntasks-per-node=';
		$pname = '-A ';		
	}
	elsif ($_ =~ /<PBS_pro_server_IP>/) { #PBS_pro system (default)
		$serv = 2;
	}
}

if ($#args == -1){
	exit;
}

my $query;

my $exc; my $sn; my @envs; my $ran; my $proj; my $mail; my $user_queue; my $home;
my $qname = "cj"; my $conda; my @module; my $quiet; my $local; my $timel;
for (my $i=0;$i<=$#args;$i++){
	if ($args[$i] eq "\-cj_exc"){
		$exc = 1;
		$args[$i] = "";
	}
	if ($args[$i] eq "\-cj_quiet"){
		$quiet = 1;
		$args[$i] = "";
	}
	if ($args[$i] eq "\-cj_local"){
		$local = 1;
		$args[$i] = "";
	}
	if ($args[$i] eq "\-cj_sn"){
		if ($args[$i+1] && $args[$i+1] !~ /^\-/){
			$ran = $args[$i+1];
			$args[$i] = "";
			$args[$i+1] = "";
			$sn = 1;
		}
		else {
			exit;
		}
	}
	if ($args[$i] eq "\-cj_env"){
		my $env = 'export PATH='.$args[$i+1].':$PATH';
		push(@envs, $env);
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_conda"){
		$conda = $args[$i+1];
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_node"){
		if ($args[$i+1] !~ /[^0-9]/){
			$nodes = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_ppn"){
		if ($args[$i+1] !~ /[^0-9]/){
			$ppn = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_time"){
		if ($args[$i+1] !~ /[^0-9\:]/){
			$timel = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_mem"){
		if ($args[$i+1] !~ /[^0-9]/){
			$mem = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_qname"){
		if ($args[$i+1] =~ /\w/){
			$qname = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_proj"){
		if ($args[$i+1] !~ /[^a-z0-9]/i){
			$proj = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_module"){
		push(@module, "module load $args[$i+1]");
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_qout"){
		if (-d $args[$i+1]){
			if ($args[$i+1] =~ /\/$/){
				$args[$i+1] =~ s/\/$//;
			}
			$home = $args[$i+1];
			if ($home =~ /\/qsub_files$/){
				$home =~ s/\/qsub_files$//;
			}
			unless ($home){
				$home = '.';
			}
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_mail"){
		if ($args[$i+1] =~ /\@/){
			$mail = "$scr_dir $smail$args[$i+1]\n$scr_dir $mtype\n";
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_queue"){
		if ($ARGV[$i+1] =~ /^trans|^ct/){
			$user_queue = $ARGV[$i+1];
		}
		else {
			print "$user_queue is incorrect, skipped.\n";
			$user_queue = "";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
}
unless ($home){
	$home = $home1;
}
my $c_line = join(" ", @args);
$c_line =~ s/\s+/ /g;
$c_line =~ s/^\s+//;

if ($local == 1){
	goto COMMAND;
}

unless (-d "$home\/qsub_files"){
	system ("mkdir $home\/qsub_files");
}
unless (-d "$home\/qsub_files\/out"){
	system ("mkdir $home\/qsub_files\/out");
}
RE:
unless ($ran){
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "$home\/qsub_files\/$ran\_$qname\.q" && $sn != 1){
	$ran = undef;
	goto RE;
}

my $adjust; my $check_ppn; my $check_nodes;
if ($mem == 0){
	$mem = "";
}
else {
	if ($serv == 1 || $serv == 2){ #every core only has 6 gb memory, so if request large memory, adjust core number based on memory value.
        $check_ppn = $mem / 6;
        $check_ppn++ if ($check_ppn > int($check_ppn));
        if ($ppn < int($check_ppn)){
            $ppn = int($check_ppn);
            $adjust = "Adjusted ncpus\/mpiprocs to $ppn based on memory request.\n";
        }
        if ($serv == 1){
        	$mem = "$scr_dir \-l mem\=$mem\gb\n";
        }
        if ($serv == 2){
        	$mem = "\:mem\=$mem\gb";
        }
        if ($serv == 3){
        }
	}
	if ($serv == 3){
        $check_ppn = $mem / 12;
        $check_ppn++ if ($check_ppn > int($check_ppn));
        if ($ppn < int($check_ppn)){
            $ppn = int($check_ppn);
            $adjust = "Adjusted ncpus\/mpiprocs\/ntasks-per-node to $ppn based on memory request.\n";
        }
	    if ($mem == 0){
			$mem = "";
		}
		else {
			$mem = "$scr_dir --mem\=$mem\G\n";
		}
	}
}
if ($serv == 1){
	$check_nodes = $ppn / 12;
}
if ($serv == 2){
	$check_nodes = $ppn / 20;
}
if ($serv == 3){ #Taiwania 3
	my @avail_proj = `get_su_balance`; my $check_proj = 0; my $balance;
	chomp(@avail_proj);
	@avail_proj = grep {$_ ne ""} @avail_proj;
	foreach (@avail_proj){
		if ($_ =~ /$proj/){
			$check_proj = 1;
			$_ =~ s/[\{\}\"]//g;
			my @tmp = split(/\,/, $_); #balance is the last one
			chomp(@tmp);
			#print "debug: @tmp\n";
			my @proj_tmp = split(/\:/, $tmp[0]);
			my @bal_tmp = split(/\:/, $tmp[-1]);
			$balance = $proj_tmp[1];
			$proj = $proj_tmp[1];
		}
	}
	if ($check_proj == 0){
		unless ($quiet == 1){
			print "\nProject is not found\/defined, selecting the best project for the job...";
		}
		foreach (@avail_proj){
			$_ =~ s/[\{\}\"]//g;
			my @tmp = split(/\,/, $_);
			chomp(@tmp);
			my @proj_tmp = split(/\:/, $tmp[0]); #project ID is the first one
			#print "debug: @proj_tmp\n";
			my @bal_tmp = split(/\:/, $tmp[-1]); #balance is the last one
			unless ($proj){
				$proj = $proj_tmp[1];
				$balance = $bal_tmp[1];
			}
			else {
				if ($balance < $bal_tmp[1]){
					$balance = $bal_tmp[1];
					$proj = $proj_tmp[1];
				}
			}
		}
	}
	unless ($quiet == 1){
		print BOLD "\nProject: $proj is used. $balance balance is available.\n", RESET;
	}
	if ($balance <= 0){
		die "There is no balance for the project.\n";
		print "Project condition\(s\):\n";
		system("get_su_balance");
		exit;
	}
	elsif ($balance < 10){
		print "WARNING: The available balance for the project is less than 10.\n";
		print "Project condition\(s\):\n";
		system("get_su_balance");
	}
	elsif ($balance < 5){
		print "WARNING: The available balance for the project is less than 5.\n";
		print "Project condition\(s\):\n";
		system("get_su_balance");
	}
	$proj = "$scr_dir $pname$proj\n";
	if ($user_queue){
		$query = "$scr_dir $par $user_queue\n";
	}
	if ($ppn <= 56){
		$query = "$scr_dir $par ct56\n";
	}
	elsif ($ppn <= 224){
		$query = "$scr_dir $par ct224\n";
	}
	elsif ($ppn <= 560){
		$query = "$scr_dir $par ct560\n";
	}
	elsif ($ppn <= 2240){
		$query = "$scr_dir $par ct2k\n";
	}
	else {
		$query = "$scr_dir $par ct8k\n";
	}
	#temp_use
	my $time = scalar localtime();
	$time =~ s/[\[\]]//g;
	my @time_f = split(/\s+|\t/, $time);
	if ($time_f[1] eq 'Apr' && $time_f[-1] eq '2024'){
		if ($time_f[2] <= 27){
			$query = "$scr_dir $par trans\n";
		}
	}
	#temp_use
	$check_nodes = $ppn / 56; #56 threads in 1 node
}
$check_nodes++ if ($check_nodes > int($check_nodes));
if ($nodes < int($check_nodes)){
	$nodes = int($check_nodes);
}

my $m_out = $mem;
if ($m_out eq ""){
	$m_out = "system default";
} else {
    $m_out =~ s/\:mem\=|$scr_dir \-l mem\=|$scr_dir --mem\=//g;
}

my $t_out = $timel;
if ($timel){
	$t_out = "$scr_dir $wall$timel\n";
}
else {
	$timel = "system default";
}

my $q_type = $query;
$q_type =~ s/$scr_dir|$par|\s//g;
unless ($q_type){
    $q_type = "system default";
}
unless ($quiet == 1){
	print BOLD "\nInfo of the job:\n", RESET;
	print "serial number: $ran\njob name: $qname\n";
}
if ($serv == 3){
	my $p_name = $proj;
	$p_name =~ s/$scr_dir|-P|-A|\s//gi;
	unless ($quiet == 1){
		print "project name: $p_name\n";
	}
}
unless ($quiet == 1){
	print "queue name: $q_type\nnode\(select\): $nodes\nppn\(ncpus\/mpiprocs\): $ppn\nmem: $m_out\nwalltime: $timel\n";
	if ($adjust){
		print $adjust;
	}
}

open (INPUT, ">$home\/qsub_files\/$ran\_$qname\.q") || die BOLD "Cannot write $home\/qsub_files\/$ran\_$qname\.q: $!", RESET, "\n";
COMMAND:
unless ($local == 1){
	if ($serv == 1){
    	print INPUT "\#\!\/bin\/bash\n$scr_dir \-l nodes\=$nodes\:ppn\=$ppn\n$mem$t_out$mail$scr_dir \-o $home\/qsub_files\/out\/$ran\_$qname\.out \-j oe\n";
	}
	if ($serv == 2){
    	print INPUT "\#\!\/bin\/bash\n$scr_dir \-l select\=$nodes\:ncpus\=$ppn\:mpiprocs=$ppn$mem\n$t_out$mail$scr_dir \-o $home\/qsub_files\/out\/$ran\_$qname\.out \-j oe\n$scr_dir $jenv\n";
	}
	if ($serv == 3){
		unless ($local == 1){
			#if ($mem == 0){
			#	$mem = "";
			#}
			#else {
			#	$mem = "$scr_dir --mem\=$mem\G\n";
			#}
			my $j_name = "$scr_dir $jn$ran\_$qname\n";
    		print INPUT "\#\!\/bin\/bash\n$proj$query$j_name$scr_dir $snode$nodes\n$scr_dir $sppn$ppn\n$mem$t_out$mail$scr_dir \-o $home\/qsub_files\/out\/$ran\_$qname\.out\n$scr_dir $jenv\n\n";
    		#for taiwania 3
    		#print INPUT "ml load old-module\nml load biology\nml load biology\/GATK\/4.2.3.0\n\n";
    		#for taiwania 3
   	 	}
	}
}
if (@module){
    if ($local == 1){
    	foreach (@module){
    		system("$_");
    	}
    }
    else {
		print INPUT join("\n", @module), "\n";
	}
}
unless ($local == 1){
	unless ($serv == 3){
		print INPUT "\ncd \$PBS_O_WORKDIR\n";
	}
}
if ($conda && ($serv == 1 || $serv == 2)){
	unless ($quiet == 1){
		print "conda env: $conda\n";
    }
    if ($local == 1){
    	system("conda activate $conda");
    }
    else {
    	print INPUT "source activate $conda\n";
    }
}
if (@envs){
	foreach (@envs){
		unless ($quiet == 1){
			print "env setting: $_\n";
		}
		if ($local == 1){
			system("$_");
		}
		else {
			print INPUT "$_\n";
		}
	}
}

my @c_lines;
if ($c_line =~ /\\n/){
	my @c_lines = split(/\\n/, $c_line);
	foreach (@c_lines){
        $_ =~ s/^\s+|\s+$//;
        if ($local == 1){
        	system("$_");
        }
        else {
        	print INPUT "$_\n";
        }
	}
}
else {
	if ($local == 1){
		system("$c_line");
	}
	else {
		print INPUT "$c_line\n";
	}
}
if ($conda  && ($serv == 1 || $serv == 2)){
	if ($local == 1){
		system("conda deactivate");
	}
	else {
    	print INPUT "conda deactivate\n";
    }
}
unless ($local == 1){
	close(INPUT);
	unless ($quiet == 1){
		print BOLD "\nThe qsub file is: $home\/qsub_files\/$ran\_$qname\.q\n", RESET;
	}
}
if ($exc == 1 && $local != 1){
	my @tmp = split(/\//, $home1);
	my $uid = $tmp[-1];
	if ($serv == 3){
		&status_tw3($uid, "sbatch $home\/qsub_files\/$ran\_$qname\.q", $quiet, $q_type);
	}
	else {
		&status_other($ran, $uid, "qsub $home\/qsub_files\/$ran\_$qname\.q", $quiet, $q_type);
	}
}
#always delete qsub files and log files older than 90 days
system("find $home\/qsub_files \! -type d -mtime \+90 -exec rm -f \{\} \+");
system("find $home\/qsub_files\/out \! -type d -mtime \+90 -exec rm -f \{\} \+");
return "qsub $home\/qsub_files\/$ran\_$qname\.q";
}
sub rnd_str {
	join("", @_[map{rand@_} 1..shift]);
} #generate serial number
sub status_tw3 {
	my $time = scalar localtime();
	my $uid = shift; my $job_q = shift; my $quiet = shift; my $q_type = shift;
	my $job_count; my $check_sent = 0;# my $core_count; 
	unless ($quiet == 1){
		print "\[$time\]\: Press ctrl \+ c to terminate this script.\n";
		print "\[$time\]\: WARNING: If you terminate this script, the following qsub job\(s\)\/step\(s\) will not be generated and executed. If you run this script in background using \"nohup\", you can ignore this message.\n";
		print "\[$time\]\: Looking for $uid job\(s\)...\n";
	}
	do {
		$time = scalar localtime();
		my @stat;
		do {
			@stat = `squeue -u $uid 2\>\&1`;
		} until ($stat[0] !~ /Socket timed out/i);
		foreach (@stat){
			if ($_ =~ /$uid/i && $_ =~ /$q_type\b/){
				$job_count++; #running and padding jobs are counted together
				#@temp = split(/\s+|\t+/, $_);
				#$core_count += $temp[4];
			}
		}
		#$core_count += $ppn;
		if (($q_type eq "ct224" && $job_count >= 75) || ($q_type eq "ct56" && $job_count >= 80) || ($q_type eq "ct560" && $job_count >= 45) || ($q_type eq "ct2k" && $job_count >= 18) || ($q_type eq "ct8k" && $job_count >= 6) || ($q_type eq "trans" && $job_count >= 30)){
			unless ($quiet == 1){
				print "\rYour request is over the limitation in the waiting\/running list. Some jobs will be sent later.";
			}
			sleep(10);
		}
		else {
			my $repeat = 0; my $violate = 0;
			do {
				$violate = 0;
				my $tmp; my $resend = 0;
				do {
					$tmp = `$job_q 2\>\&1`;
					$resend++;
					sleep(3);
				} until ($tmp !~ /submission failed|error/i || $resend == 50);
				print "\n$job_q job is sent.\n";
				my $show_job_count = $job_count + 1;
				if ($tmp =~ /violates/){
					$violate = 1;
				}
				else {
					print "$tmp";
				}
				if ($show_job_count == 1){
					print "$show_job_count job is on the list.\n";
				}
				else {
					print "$show_job_count jobs are on the list.\n";
				}
				$repeat++;
			} until ($repeat == 3 || $violate == 0);
			if ($violate == 1){
				print "Job: $job_q violates queue and\/or server resource limits.\n";
				exit;
			}
			$check_sent = 1;
		}
	} while ($check_sent == 0);
	1;
}
sub status_other {
	my $time = scalar localtime();
	my $ran = shift; my $uid = shift; my $job_q = shift; my $quiet = shift; my $q_type = shift;
	my $job_count; my $core_count; my $check_sent = 0;
	unless ($quiet == 1){
		print "\[$time\]\: Press ctrl \+ c to terminate this script.\n";
		print "\[$time\]\: WARNING: If you terminate this script, the following qsub job\(s\)\/step\(s\) will not be generated and executed. If you run this script in background using \"nohup\", you can ignore this message.\n";
		print "\[$time\]\: Looking for $ran job\(s\)...\n";
	}
	do {
		$time = scalar localtime();
		my @stat = `qstat -u $uid`;
		my @temp;
		$job_count = 0;
		foreach (@stat){
			if ($_ =~ /$ran/i){
				$job_count++;
				@temp = split(/\s+/, $_);
				$core_count += $temp[6];
			}
		}
		if ($job_count >= 200){
			unless ($quiet == 1){
				print "\rYour jobs are too many. Some jobs will be sent later.";
			}
			sleep(30);
		}
		else {
			my $tmp = `$job_q`;
			print "\n$job_q job is sent.\n";
			my $show_job_count = $job_count + 1;
			print "$tmp";
			if ($show_job_count == 1){
				print "$show_job_count job is on the list.\n";
			}
			else {
				print "$show_job_count jobs are on the list.\n";
			}
			$check_sent = 1;
		}
	} while ($check_sent == 0);
	1;
}
sub status {
	my $time = scalar localtime();
	my $ran = shift; my $uid = shift;
	my $job_count;
	my @server = `ip route get 1.2.3.4 \| awk \'\{print \$7\}\'`;
	chomp(@server);
	my $serv;
	foreach (@server){
    	if ($_ =~ /172.28.111/){ #Taiwania 3
			$serv = 1;
			last;
		}
		elsif ($_ =~ /140.112.2/) { #NTU servers
			$serv = 2;
			last;
		}
		else {
			$serv = 1;
		}
	}
	print "\[$time\]\: Press ctrl \+ c to terminate this script.\n";
	print "\[$time\]\: WARNING: If you terminate this script, the following qsub job\(s\)\/step\(s\) will not be generated and executed. If you run this script in background using \"nohup\", you can ignore this message.\n";
	print "\[$time\]\: Looking for $ran job\(s\)...\n";
	do {
		$time = scalar localtime();
		my @stat; my $resend = 0; $err;
		if ($serv == 1){
			if ($uid){
				do {
					@stat = `squeue -u $uid 2\>\&1`;
					$err = 0;
					$resend++;
					sleep(3);
					foreach (@stat){
						if ($_ =~ /error/i){
							$err = 1;
						}
					}
				} until ($err == 0 || $resend == 50);
			}
			else {
				do {
					@stat = `squeue 2\>\&1`;
					$err = 0;
					$resend++;
					sleep(3);
					foreach (@stat){
						if ($_ =~ /error/i){
							$err = 1;
						}
					}
				} until ($err == 0 || $resend == 50);
			}
		}
		else {
			@stat = `qstat \-G`;
		}
		my @temp;
		my @jobs;
		$job_count = 0;
		foreach (@stat){
			if ($uid){
				my $uid_cut;
				if (length($uid) > 8){
					$uid_cut = substr $uid, 0, 8;
				}
				else {
					$uid_cut = $uid;
				}
				if ($_ =~ /$ran/ && $_ =~ /$uid_cut/){
					$job_count += 1;
					@temp = split(/\s+/, $_);
					push(@jobs, $temp[1]);
				}
			}
			elsif ($_ =~ /$ran/){
				$job_count += 1;
				@temp = split(/\s+/, $_);
				push(@jobs, $temp[1]);
			}
		}
		if ($job_count != 0){
			my $job = join("\t", @jobs);
			if ($job_count == 1){
				print "\r\[$time\]\: $job_count "; 
				print "$ran";
				print " job is still running!          ";
			}
			elsif ($job_count > 1){
				print "\r\[$time\]\: $job_count "; 
				print "$ran";
				print " jobs are still running!      ";
			}
			sleep(30);
		}
		else {
			print "\r\[$time\]\: no job is running!                           ";
		}
		@jobs = ();
	} while ($job_count != 0);
	print "\n";
	1;
} #present status of each step
sub get_1_eles{
     my @vcfs = @{$_[0]}; my @ids = @{$_[1]};
     my @ordered;
     foreach $id (@ids){
        foreach $vcf (@vcfs){
        	if ($vcf =~ /gz$/){
            	open(INPUT, "-|", "gzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
            }
            else {
            	open(INPUT, "<$vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
            }
            while (my $line = <INPUT>){
                chomp($line);
                if ($line =~ /^\#/){
                    next;
                }
                else {
                    my @eles = split(/\t+|\s+/, $line);
                    if ($eles[0] eq $id){
                        push(@ordered, $vcf);
                    }
                    last;
                }
            }
            close(INPUT);
        }
    }
    return @ordered;
} #reorder files as the contig order in the vcf file
sub chr_name {
	$time = scalar localtime();
	my $file = shift; my $pre = shift;
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
			if ($pre){
                if ($line[3] =~ /^$pre/i){
                	my @tmp = split(/\s+/, $line[3]);
                	$line[3] = $tmp[0];
                    push(@id, $line[3]);
                }
			}
			else {
                my @tmp = split(/\s+/, $line[3]);
                $line[3] = $tmp[0];
                push(@id, $line[3]);
			}
		}
	}
	return (@id);
} #get chromosome name
sub chr_lengths {
	$time = scalar localtime();
	my $vcf = shift; my $pre = shift;
	my @content; my @line; my @len;
	if ($vcf =~ /\.vcf\.gz$/){
		@content = `gzip \-cd $vcf \| head \-n 5000`;
	}
	elsif ($vcf =~ /\.vcf$/){
		@content = `head -n 1000 $vcf`;
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
			if ($pre){
                if ($line[3] =~ /^$pre/i){
                    push(@len, $line[5]);
                }
			}
			else {
                push(@len, $line[5]);
			}
		}
	}
	return (@len);
} #get interval length from vcf header
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
sub check_path {
	my $path = shift;
	my $dir = getcwd;
	if ($path =~ /\// || $path =~ /\.\./){
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
				}
				$path = join("\/", @dir_eles)."\/".join("\/", @path_eles);
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
	if ($path =~ /\/$/){
        $path =~ s/\/$//;
	}
	return($path);
} #relative path to absolute path

1;
