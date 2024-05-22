#!/usr/bin/perl
use Term::ANSIColor qw(:constants);

chomp(@ARGV);
if (-e $ARGV[0]){
	open (LIST, "<$ARGV[0]") || die "Cannot open $ARGV[0]: $!\n";
}
my @lists = <LIST>;
chomp(@lists);
close(LIST);

my $line; my @labels; my $m_nex;
my @nexs = <*.nex>;
chomp(@nexs);
foreach my $nex (@nexs){
	print "Working on $nex...\n";
	open (NEX, "<$nex") || die "Cannot open $nex: $!\n";
	open (OUT, ">m_$nex") || die "Cannot open m_$nex: $!\n";
	$m_nex = <m_*.nex>;
	while ($line = <NEX>){
		chomp($line);
		foreach (@lists){
			my @eles = split(/\t/, $_);
			if ($line =~ /$eles[0]/){
				$line =~ s/$eles[0]/$_/g;
				push(@labels, $_);
			}
		}
		$line =~ s/[\x0A\x0D]//g;
		print OUT "$line\n";
	}
	close(OUT);
}
my %seen;
my @unique = do { %seen; grep { !$seen{$_}++ } @labels };


my @colors = ("\#F8766D", "\#39B600", "\#00B0F6", "\#9590FF", "\#D89000", "\#00BF7D", "\#E76BF3", "\#A3A500", "\#00BFC4", "\#FF62BC");
my @ids; my @groups; 
my @pca_table = <PCA_*_table.txt>;
chomp(@pca_table);
if ($pca_table[0] =~ /txt$/){
#	print "@pca_table\n";
	foreach (@pca_table){
		print "Working on iTOL_styles_$_...\n";
		my $tree = $_;
		$tree =~ s/PCA/colored/;
		$tree =~ s/table\.txt/tree\.nex/;
		print "Working on $tree...\n";
		open (OUT2, ">iTOL_styles_$_") || die "Cannot write iTOL_styles_$_: $!\n";
		open (OUT3, ">$tree") || die "Cannot write $tree: $!\n";
		open (IN, "<$m_nex") || die "Cannot open $m_nex: $!\n";
		my @nex_c = <IN>;
		chomp(@nex_c);
		close(IN);
		print OUT2 "TREE_COLORS\n\nSEPARATOR SPACE\n\nDATA\n\n";
		open (TABLE, "<$_") || die "Cannot open $_: $!\n";
		my @table = <TABLE>;
		chomp (@table);
		shift(@table);
		close(TABLE);
		foreach my $element (@table){
			my @tmp = split(/\t+/, $element);
			push(@ids, $tmp[0]);
			$tmp[1] =~ s/G//;
			$tmp[1] = int($tmp[1]);
			push(@groups, $tmp[1]);
		}
		foreach my $label (@unique){
			$label =~ s/[\x0A\x0D]//g;
			for (my $i=0; $i<=$#ids; $i++){
				if ($label =~ /$ids[$i]/){
					my $color_pic = $groups[$i];
					$color_pic -= 1;
					if ($color_pic == 98){
						print OUT2 "$label branch \#585858 normal 1\n";
					}
					else {
						print OUT2 "$label branch $colors[$color_pic] normal 1\n";
					}
				}
			}
		}
		my $end = 0;
		foreach my $line (@nex_c){
			if ($line =~ /END/i){
				$end = 1;
			}
			for (my $i=0; $i<=$#ids; $i++){
				if ($end == 1){
					next;
				}
				elsif ($line =~ /$ids[$i]/){
					my $color_pic = $groups[$i];
					$color_pic -= 1;
					if ($color_pic == 98){
						$line = "$line\[\&\!color\=\#585858\]";
					}
					else {
						$line = "$line\[\&\!color\=$colors[$color_pic]\]";
					}
				}
			}
			print OUT3 "$line\n";
		}
		close(OUT2);
		close(OUT3);
		@ids = ();
		@groups = ();
	}
}

sub unique {
	my %seen;
	@unique = do { %seen; grep { !$seen{$_}++ } @Ls };
}