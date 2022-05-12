#!/usr/bin/perl -w

use strict;
use warnings;

my %options = &parse_commandline (@ARGV);

my $test_indvs = &read_first_columns ($options{'--pop'});
my $ref_haps = &read_first_columns ($options{'--ref'});

my ($vcf_sites, $vcf_test_indvs) = &read_vcf_spe ($test_indvs, $ref_haps, $options{'--vcf'}, $options{'--ref-mac'});
&check_chr ($vcf_sites);
my $valid_test_hap_counts = &check_indv ($vcf_sites, $vcf_test_indvs);

print ("\n", join ("\t", 'CHROM', 'POS', 'ID', @$vcf_test_indvs), "\n");
foreach my $site (@$vcf_sites) {
	print (join ("\t", $site->{'CHROM'}, $site->{'POS'}, $site->{'ID'}));
	foreach my $gt (@{$site->{'GT'}}) {
		print ("\t", $gt->[1], $gt->[0], $gt->[2]);
	}
	print ("\n");
}

print ("\nID\tS*\tset of sites\n");
foreach my $test_hap_count (@$valid_test_hap_counts) {
	my @vcf_sites_use = ();
	foreach my $site (@$vcf_sites) {
		if ($site->{'GT'}[$test_hap_count->{'indv'}][$test_hap_count->{'phase'}] eq '1') {
			push (@vcf_sites_use, $site);
		}
	}
	my @rand_hap_counts;
	if ($options{'--size'} >= @$valid_test_hap_counts) {@rand_hap_counts = @$valid_test_hap_counts;}
	else {
		my @valid_test_hap_counts_copy = @$valid_test_hap_counts;
		for (my $size = 0; $size < $options{'--size'}; $size ++) {
			push (@rand_hap_counts, splice (@valid_test_hap_counts_copy, int (rand (@valid_test_hap_counts_copy)), 1, ()));
		}
		if (!grep {$_ == $test_hap_count} @rand_hap_counts) {splice (@rand_hap_counts, 0, 1, $test_hap_count);}
	}
	print ($vcf_test_indvs->[$test_hap_count->{'indv'}], '-', $test_hap_count->{'phase'});
	&calculate_sstar (\@vcf_sites_use, \@rand_hap_counts);
}


exit 0;

sub parse_commandline {
	my @argvs = @_;
	my $argc = @argvs;
	my %options;
	my @valid_selects = ('--vcf', '--id', '--pop', '--ref', '--ref-mac', '--size', '--iter');
	for (my $i = 0; $i < $argc; $i ++) {
		if (grep (/^($argvs[$i])$/, @valid_selects) && ($i + 1 < $argc)) {
			$options{$argvs[$i]} = $argvs[$i + 1];
			$i ++;
		}
		else {
			&usage;
		}
	}
	if (!exists ($options{'--vcf'}) || !exists ($options{'--pop'}) || !exists ($options{'--ref'})) {
		&usage;
	}
	if (!exists ($options{'--ref-mac'})) {$options{'--ref-mac'} = 0;}
	if (!exists ($options{'--size'})) {$options{'--size'} = 40;}
#	if (!exists ($options{'--iter'})) {$options{'--iter'} = 5;}
	if ($options{'--ref-mac'} =~ /\D/) {&usage;}
	print ("Arguments:\n");
	foreach my $argv_names (@valid_selects) {
		if (exists ($options{$argv_names})) {
			print ($argv_names, ' ', $options{$argv_names}, "\n");
		}
	}
	print ("\n");
	return %options;
}
sub usage {
	print STDERR ("Usage: sstar.pl [--vcf | --pop | --ref | --ref-mac | --size] <args>\n");
	print STDERR ("Mandatory arguments:\n");
	print STDERR (" --vcf     VCF file\n");
#	print STDERR (" --id   haplotype ID\n");
	print STDERR (" --pop     Test population individual IDs file\n");
	print STDERR (" --ref     Reference population haplotype IDs file\n");
	print STDERR ("Optional arguments:\n");
	print STDERR (" --ref-mac Maximum limitation number of minor alleles in reference population (0)\n");
	print STDERR (" --size    number of haplotypes in same population to calculate S* (40)\n");
#	print STDERR (" --iter    times to calculate S* for one haplotype (5)\n");
	exit 1;
}

sub read_first_columns {
	my $in_filename = $_[0];
	open (my $fh_in, '<', $in_filename) or die "Cannot open $in_filename: $!";
	my @items = ();
	while (my $line = <$fh_in>) {
		chomp ($line);
		(my $item, undef) = split ("\t", $line, 2);
		push (@items, $item);
	}
	return \@items;
}

sub read_vcf_spe {
	my ($test_indvs, $ref_haps, $vcf_filename, $lim_ref_mac) = @_;
	open (my $fh_vcf, '<', $vcf_filename) or die "Cannot open $vcf_filename: $!";
	my $line_vcf;
	for ($line_vcf = <$fh_vcf>; defined ($line_vcf) && $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
	if (!(defined ($line_vcf) && $line_vcf =~ /^#CHROM/)) {
		print STDERR ("Error: No title line in $vcf_filename.\n");
		close ($fh_vcf);
		exit 2;
	}
	chomp ($line_vcf);
	my @columns = split ("\t", $line_vcf);
	if (@columns < 10) {
		print STDERR ("Error: No sample identifiers in $vcf_filename.\n");
		close ($fh_vcf);
		exit 2;
	}
	my %test_sample_counts;
	my %ref_sample_counts;
	my @valid_test_indvs;
	foreach my $column_count (9 .. $#columns) {
		if (grep {$_ eq $columns[$column_count]} @$test_indvs) {
			$test_sample_counts{$columns[$column_count]} = $column_count - 9;
			push (@valid_test_indvs, $columns[$column_count]);
		}
		if (grep {$_ eq $columns[$column_count] . '-1'} @$ref_haps) {
			$ref_sample_counts{$columns[$column_count] . '-1'}{'column'} = $column_count - 9;
			$ref_sample_counts{$columns[$column_count] . '-1'}{'phase'} = 1;
		}
		if (grep {$_ eq $columns[$column_count] . '-2'} @$ref_haps) {
			$ref_sample_counts{$columns[$column_count] . '-2'}{'column'} = $column_count - 9;
			$ref_sample_counts{$columns[$column_count] . '-2'}{'phase'} = 2;
		}
	}
	my @absent_ref_haps = ();
	foreach my $ref_hap (@$ref_haps) {
		if (!defined ($ref_sample_counts{$ref_hap})) {
			push (@absent_ref_haps, $ref_hap);
		}
	}
	if (@absent_ref_haps) {
		print STDERR ('Error: Reference haplotype ', join (', ', @absent_ref_haps), (@absent_ref_haps == 1) ? ' is' : ' are' , " not in $vcf_filename\n");
		close ($fh_vcf);
		exit 2;
	}
	
	my @vcf_sites;
	while ($line_vcf = <$fh_vcf>) {
		chomp ($line_vcf);
		@columns = split ("\t", $line_vcf);
		my $chr = $columns[0];
		my $pos = $columns[1];
		my $id  = $columns[2];
		(my $format1, undef) = split (':', $columns[8], 2);
		if ($format1 ne 'GT') {
			print ("Skipped: $id (chr$chr:$pos) GT sub-field not in the first of GT field.\n");
			next;
		}
		my @ref_alts = ($columns[3], split (',', $columns[4]));
		my @genotypes;
		foreach my $column_count (9 .. $#columns) {
			push (@genotypes, [&gt_phase ($columns[$column_count])]);
		}
		my %allele_counts;
		foreach my $genotype (@genotypes) {
			foreach my $gt (@$genotype) {
				if ($gt =~ /^\d+$/) {$allele_counts{$gt} ++;}
			}
		}
		my @alleles = keys (%allele_counts);
		if (@alleles != 2) {
			print ("Skipped: $id (chr$chr:$pos) number of allele is not 2 but ", scalar (@alleles) , ".\n");
			next;
		}
		my $skip = 0;
		foreach my $allele (@alleles) {
			if ($ref_alts[$allele] !~ /^[ATGCN]$/i) {
				print ("Skipped: $id (chr$chr:$pos) allele $allele is ", $ref_alts[$allele] , ".\n");
				$skip = 1;
			}
		}
		if ($skip) {next;}
		if ($allele_counts{$alleles[0]} < $allele_counts{$alleles[1]}) {@alleles = ($alleles[1], $alleles[0]);}
		$skip = 0;
		my $ref_mac = 0;
		foreach my $ref_hap (@$ref_haps) {
			my $phase = $ref_sample_counts{$ref_hap}{'phase'};
			my $gt = $genotypes[$ref_sample_counts{$ref_hap}{'column'}];
			if ($gt->[0] eq '/') {
				if ($gt->[1] eq $alleles[1] && $gt->[2] eq $alleles[1]) {
					$ref_mac ++;
					if ($ref_mac > $lim_ref_mac) {
						$skip = 1;
						last;
					}
				}
				elsif (($gt->[1] eq $alleles[1]) || ($gt->[2] eq $alleles[1])) {
					$ref_mac += 0.5;
					if ($ref_mac > $lim_ref_mac) {
						$skip = 1;
						last;
					}
				}
			}
			elsif ($gt->[$phase] eq $alleles[1]) {
				$ref_mac ++;
				if ($ref_mac > $lim_ref_mac) {
					$skip = 1;
					last;
				}
			}
		}
		if ($skip) {
			print ("Skipped: $id (chr$chr:$pos) reference population have more than ", $lim_ref_mac, ' rare allele', ($lim_ref_mac > 1) ? 's' : '', ".\n");
			next;
		}
		$skip = 1;
		my %site = ('CHROM' => $chr, 'POS' => $pos, 'ID' => $id, 'GT' => []);
		foreach my $valid_test_indv (@valid_test_indvs) {
			foreach my $gt (@{$genotypes[$test_sample_counts{$valid_test_indv}]}) {
				if    ($gt eq $alleles[0]) {$gt = 0;}
				elsif ($gt eq $alleles[1]) {
					$gt = 1;
					$skip = 0;
				}
			}
			push (@{$site{'GT'}}, $genotypes[$test_sample_counts{$valid_test_indv}]);
		}
		if ($skip) {
			print ("Skipped: $id (chr$chr:$pos) test population do not have the rare allele.\n");
			next;
		}
		push (@vcf_sites, \%site);
	}
	close ($fh_vcf);
	return (\@vcf_sites, \@valid_test_indvs);
}

sub gt_phase {
	(my $gt, undef) = split (':', $_[0], 2);
	if ($gt =~ /^(\d+|\.)([\|\/])(\d+|\.)$/) {return ($2, $1 , $3 );}
	elsif ($gt =~ /^(\d+)$/)                 {return ('', $1 , '' );}
	elsif ($gt eq '.')                       {return ('', '.', '.');}
	else                                     {return ('', '' , '' );}
}

sub check_chr {
	my $vcf_sites = $_[0];
	my %count_by_chr;
	foreach my $site (@$vcf_sites) {
		$count_by_chr{$site->{'CHROM'}} ++;
	}
	if (keys (%count_by_chr) > 1) {
		print STDERR ('Error: Multiple chromosomes appear: ', join (',', keys (%count_by_chr)), "\n");
		exit 3;
	}
	return;
}

sub check_indv {
	my ($vcf_sites, $indvs) = @_;
	my $valid_haps = [];
	foreach my $indv_count (0 .. $#$indvs) {
		my $ploidy = 0;
		foreach my $site (@$vcf_sites) {
			my $gt = $site->{'GT'}[$indv_count];
			if ($ploidy == 0) {
				if ($gt->[0] eq '|' || $gt->[0] eq '/') {$ploidy = 2;}
				elsif ($gt->[0] eq '' && $gt->[1] =~ /^(\d+|\.)$/) {$ploidy = 1;}
			}
			elsif ($ploidy == 1) {
				if ($gt->[0] eq '|' || $gt->[0] eq '/') {
					print ('Skipped: ',$indvs->[$indv_count], " has both haploid and diploid sites\n");
					$ploidy = -1;
					last;
				}
			}
			elsif ($ploidy == 2) {
				if ($gt->[0] eq '' && $gt->[1] =~ /^\d+$/) {
					print ('Skipped: ',$indvs->[$indv_count], " has both haploid and diploid sites\n");
					$ploidy = -1;
					last;
				}
			}
			if ($gt->[0] eq '/') {
				if ($gt->[1] eq $gt->[2]) {$gt->[0] = '|';}
				else {
					print ('Skipped: ',$indvs->[$indv_count], ' unphased at ', $site->{'ID'}, '(chr', $site->{'CHROM'}, ':', $site->{'POS'}, '), GT=', $gt->[1], $gt->[0], $gt->[2], "\n");
					$ploidy = -1;
					last;
				}
			}
			elsif ($gt->[1] . $gt->[2] . $gt->[0] eq '') {
				print ('Skipped: ',$indvs->[$indv_count], ' invalid GT at ', $site->{'ID'}, '(chr', $site->{'CHROM'}, ':', $site->{'POS'}, "), GT=\n");
				$ploidy = -1;
				last;
			}
		}
		if ($ploidy == 0) {
			print ('Skipped: ', $indvs->[$indv_count], " all sites are missing.\n");
		}
		for (my $phase = 1; $phase <= $ploidy; $phase ++) {
			push (@$valid_haps, {'indv' => $indv_count, 'phase' => $phase});
		}
	}
	return $valid_haps;
}

sub calculate_sstar {
	my ($vcf_sites, $test_hap_counts) = @_;
	my $number_of_sites = @$vcf_sites;
	if ($number_of_sites < 2) {
		print ("\t///\n");
		return;
	}
	my $sstar = [];
	$sstar->[0]{'value'} = -inf;
	$sstar->[0]{'prev'}[0] = -1;
	for (my $site_y_count = 1; $site_y_count < $number_of_sites; $site_y_count ++) {
		$sstar->[$site_y_count]{'value'} = -inf;
		$sstar->[$site_y_count]{'prev'}[0] = -1;
		for (my $site_x_count = 0; $site_x_count < $site_y_count; $site_x_count ++) {
			my $bps = $vcf_sites->[$site_y_count]{'POS'} - $vcf_sites->[$site_x_count]{'POS'};
			if ($bps < 10) {last;}
			my $distance = 0;
			my $missing = 0;
			my @mac = (0, 0);
			foreach my $test_hap_count (@$test_hap_counts) {
				my @gts;
				$gts[0] = $vcf_sites->[$site_x_count]{'GT'}[$test_hap_count->{'indv'}][$test_hap_count->{'phase'}];
				$gts[1] = $vcf_sites->[$site_y_count]{'GT'}[$test_hap_count->{'indv'}][$test_hap_count->{'phase'}];
				if    ($gts[0] eq '.') {if ($gts[1] eq '1') {$missing ++;}}
				elsif ($gts[1] eq '.') {if ($gts[0] eq '1') {$missing ++;}}
				elsif ($gts[0] eq '1' xor $gts[1] eq '1') {$distance ++;}
				if ($gts[0] eq '1') {$mac[0] ++;}
				if ($gts[1] eq '1') {$mac[1] ++;}
			}
			my $s = 0;
			if ($distance == 0 && $missing <= ((($mac[0] > 2) && ($mac[1] > 2)) ? 2 : 1)) {$s = $bps + 5000;}
			elsif (1 <= $distance && $distance <= 5) {$s = -10000;}
			elsif (5 < $distance) {$s = -inf;}
			my $sstar_x;
			if ($sstar->[$site_x_count]{'value'} < 0) {$sstar_x = $s;}
			else {
				$sstar_x = $sstar->[$site_x_count]{'value'} + $s;
			}
			if ($sstar->[$site_y_count]{'value'} < $sstar_x) {
				$sstar->[$site_y_count]{'value'} = $sstar_x;
				@{$sstar->[$site_y_count]{'prev'}} = ($site_x_count);
			}
			elsif ($sstar->[$site_y_count]{'value'} == $sstar_x) {
				push (@{$sstar->[$site_y_count]{'prev'}}, $site_x_count);
			}
		}
	}
	my $sstar_max = -inf;
	my @last_site_counts = ();
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		if ($sstar_max < $sstar->[$site_count]{'value'}) {
			$sstar_max = $sstar->[$site_count]{'value'};
			@last_site_counts = ($site_count);
		}
		elsif ($sstar_max == $sstar->[$site_count]{'value'}) {
			push (@last_site_counts, $site_count);
		}
	}
	print ("\t$sstar_max");
	if (-inf < $sstar_max) {
		my @introgressed_sites_done = ();
		my @diverged_sites = ();
		my $prev_site_count = pop (@last_site_counts);
		my @introgressed_sites = ();
		while ($prev_site_count > -1) {
			unshift (@introgressed_sites, $prev_site_count);
			if ($sstar->[$prev_site_count]{'value'} < 0 && @introgressed_sites > 1) {last;}
			push (@introgressed_sites_done, $prev_site_count);
			if (@{$sstar->[$prev_site_count]{'prev'}} > 1) {
				push (@diverged_sites, $prev_site_count);
			}
			$prev_site_count = $sstar->[$prev_site_count]{'prev'}[0];
		}
		print ("\t[Start] ", join (',', (map {$vcf_sites->[$_]{'ID'} . '(' . $vcf_sites->[$_]{'POS'} . ')'} @introgressed_sites)), " [End]\n");
		foreach my $last_site_count (reverse (@last_site_counts)) {
			@introgressed_sites = ();
			$prev_site_count = $last_site_count;
			while ($prev_site_count > -1) {
				unshift (@introgressed_sites, $prev_site_count);
				if ($sstar->[$prev_site_count]{'value'} < 0 && @introgressed_sites > 1) {
					$prev_site_count = -1;
					last;
				}
				if (grep {$_ == $prev_site_count} @introgressed_sites_done) {last;}
				push (@introgressed_sites_done, $prev_site_count);
				if (@{$sstar->[$prev_site_count]{'prev'}} > 1) {
					push (@diverged_sites, $prev_site_count);
				}
				$prev_site_count = $sstar->[$prev_site_count]{'prev'}[0];
			}
			my $start = ($prev_site_count == -1) ? '[Start]' : '...';
			print (".\t.\t$start ", join (',', (map {$vcf_sites->[$_]{'ID'} . '(' . $vcf_sites->[$_]{'POS'} . ')'} @introgressed_sites)), " [End]\n");
		}
		@diverged_sites = sort {$a <=> $b} @diverged_sites;
		while (@diverged_sites) {
			my $diverged_site = pop (@diverged_sites);
			my @prev_sites_count = @{$sstar->[$diverged_site]{'prev'}};
			shift (@prev_sites_count);
			foreach my $prev_site_count (@prev_sites_count) {
				@introgressed_sites = ($diverged_site);
				while ($prev_site_count > -1) {
					unshift (@introgressed_sites, $prev_site_count);
					if ($sstar->[$prev_site_count]{'value'} < 0) {
						$prev_site_count = -1;
						last;
					}
					if (grep {$_ == $prev_site_count} @introgressed_sites_done) {last;}
					push (@introgressed_sites_done, $prev_site_count);
					if (@{$sstar->[$prev_site_count]{'prev'}} > 1) {
						push (@diverged_sites, $prev_site_count);
					}
					$prev_site_count = $sstar->[$prev_site_count]{'prev'}[0];
				}
				my $start = ($prev_site_count == -1) ? '[Start]' : '...';
				print (".\t.\t$start ", join (',', (map {$vcf_sites->[$_]{'ID'} . '(' . $vcf_sites->[$_]{'POS'} . ')'} @introgressed_sites)), " ...\n");
			}
			@diverged_sites = sort {$a <=> $b} @diverged_sites;
		}
	}
	else {print ("\n");}
	return;
}

