#!/usr/bin/perl
use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);

my $path = $ARGV[0];
if ($path =~ /\/$/){
    $path =~ s/\/$//;
}

my @f_list;
if (-d $path){
    @f_list = <$path\/*.bimbam.gz>;
}
elsif (-e $path){
    @f_list = $path;
}
else {
    die "Cannot find the bimbam files.\n";
}

foreach my $file (@f_list){
    open(INPUT, "-|", "gzip -dc $file") || die BOLD "Cannot open $file: $!", RESET, "\n";
    my $out = $file;
    $out =~ s/bimbam.gz$/m.bimbam.gz/;
    open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
    while (my $line = <INPUT>){
        chomp($line);
        $line =~ s/[\x0A\x0D]//g;
        my @eles = split(/, /, $line);
        foreach my $i (3..$#eles){
            if ($eles[$i] < 0.5){
                $eles[$i] = 0;
            }
            elsif ($eles[$i] >= 0.5 && $eles[$i] <= 1.5){
                $eles[$i] = 1;
            }
            else {
                $eles[$i] = 2;
            }
        }
        print OUT join(', ', @eles), "\n";
    }
    close(INPUT);
    close(OUT);
    system("mv $out $file");
}
