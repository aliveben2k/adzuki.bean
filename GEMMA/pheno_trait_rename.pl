#!/usr/bin/perl

my $path = $ARGV[0];
my @files = <$path\/pheno_*>;

foreach my $file (@files){
    my $new_name = $file;
    $new_name =~ s/pheno_/trait_/;
    system("mv $file $new_name");
}
