#!/usr/bin/perl -w

use strict;
use warnings;

my $in_vcf_filename;
if (@ARGV == 1) {$in_vcf_filename = $ARGV[0];}
else {&usage;}
my $outfile_prefix = &outfilename ($in_vcf_filename);

print ("             Input VCF file: $in_vcf_filename\n");
print (" Output haplotype data file: $outfile_prefix.hap\n");
print ("Output map information file: $outfile_prefix.inp\n");

open (my $fh_in_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename: $!";
my $vcf_file = &read_vcf_snp ($fh_in_vcf);
close ($fh_in_vcf);

&allele_to_digit ($vcf_file->{'site'});

open (my $fh_out_map, '>', $outfile_prefix . '.inp') or die "Cannot open $outfile_prefix.inp: $!";
&print_map_file ($fh_out_map, $vcf_file->{'site'});
close ($fh_out_map);

open (my $fh_out_hap, '>', $outfile_prefix . '.hap') or die "Cannot open $outfile_prefix.hap: $!";
&print_hap_file ($fh_out_hap, $vcf_file);
close ($fh_out_map);

exit 0;

sub usage {
	print ("Usage: vcf2rehh.pl <VCF_file>\n");
	print ("<VCF_file> is converted to input files of rehh\n");
	exit 1;
}

sub outfilename {
	my $in_filename = $_[0];
	my $out_filename;
	if ($in_filename =~ /^(.+)\.vcf$/i) {$out_filename = $1;}
	else {$out_filename = $in_filename;}
	return $out_filename;
}

sub read_vcf_snp {
	my $fh_in_vcf = $_[0];
	my $line_vcf;
	my $vcf_file = {};
	for ($line_vcf = <$fh_in_vcf>; defined ($line_vcf) && ($line_vcf =~ /^##/); $line_vcf = <$fh_in_vcf>) {}
	if (!defined ($line_vcf)) {
		print STDERR ("No contents in the VCF file.\n");
		close ($fh_in_vcf);
		exit 1;
	}
	chomp ($line_vcf);
	$vcf_file->{'sample'} = &vcf_samples ($line_vcf);
	$vcf_file->{'site'} = [];
	while ($line_vcf = <$fh_in_vcf>) {
		chomp ($line_vcf);
		my @columns_vcf = split ("\t", $line_vcf);
		my %vcf_site;
		$vcf_site{'CHROM'} = $columns_vcf[0];
		$vcf_site{'POS'} = $columns_vcf[1];
		$vcf_site{'ID' } = ($columns_vcf[2] ne '.') ? $columns_vcf[2] : $columns_vcf[0] . ':' . $columns_vcf[1];
		$vcf_site{'REF'} = uc ($columns_vcf[3]);
		$vcf_site{'ALT'} = uc ($columns_vcf[4]);
		if ($columns_vcf[7] =~ /\bAA=([ATGCNatgcn-]+)\b/) {$vcf_site{'AA'} = uc ($1);}
		else {$vcf_site{'AA' } = '';}
		if (@columns_vcf > 9) {
			if ($columns_vcf[8] =~ /^GT:/i || $columns_vcf[8] =~ /^GT$/i) {
				foreach my $gt (@columns_vcf[9 .. $#columns_vcf]) {
					($gt, undef) = split (':', $gt, 2);
				}
			}
			$vcf_site{'GTs'} = [@columns_vcf[9 .. $#columns_vcf]];
		}
		if ($vcf_site{'REF'} =~ /^[ATGCN]$/ && $vcf_site{'ALT'} =~ /^[ATGCN]$/) {
			push (@{$vcf_file->{'site'}}, \%vcf_site);
		}
		else {
			print ('Chr', $vcf_site{'CHROM'}, ':', $vcf_site{'POS'}, ' (', $vcf_site{'ID'}, ') is not a 2 allelic SNP. ');
			print ('REF=', $vcf_site{'REF'}, ', ALT=', $vcf_site{'ALT'}, ". It is skipped.\n");
		}
	}
	return $vcf_file;
}

sub vcf_samples {
	my $title_line = $_[0];
	my @vcf_titles = split ("\t", $title_line);
	my $number_vcf_titles = (@vcf_titles > 9) ? 9 : @vcf_titles;
	if ($number_vcf_titles < 8) {
		print ("Warning: The header line of the VCF file is too short.\n");
	}
	my @true_titles = ('#CHROM', qw/POS ID REF ALT QUAL FILTER INFO FORMAT/);
	for (my $column_count = 0; $column_count < $number_vcf_titles; $column_count ++) {
		if ($vcf_titles[$column_count] ne $true_titles[$column_count]) {
			print ('Warning: The header ', $vcf_titles[$column_count], ' is wrong. It has to be ', $true_titles[$column_count], ".\n");
		}
	}
	return (@vcf_titles <= 9) ? [] : [@vcf_titles[9 .. $#vcf_titles]];
}

sub allele_to_digit {
	my $sites = $_[0];
	foreach my $site (@$sites) {
		if ($site->{'AA'} eq $site->{'ALT'}) {
			$site->{'REF'} = '2';
			$site->{'ALT'} = '1';
		}
		else {
			if ($site->{'AA'} ne $site->{'REF'}) {
				print ('Chr', $site->{'CHROM'}, ':', $site->{'POS'}, ' (', $site->{'ID'}, ')  REF=', $site->{'REF'}, ' is assumed to be "ancestral." ');
				print (' (originally AA=', $site->{'AA'}, ', REF=', $site->{'REF'}, ', ALT=', $site->{'ALT'}, ")\n");
			}
			$site->{'REF'} = '1';
			$site->{'ALT'} = '2';
		}
	}
	return;
}

sub print_map_file {
	my ($fh_out, $sites) = @_;
	foreach my $site (@$sites) {
		print $fh_out (join (' ', $site->{'ID'}, $site->{'CHROM'}, $site->{'POS'}, '1', '2'), "\n");
	}
	return;
}

sub print_hap_file {
	my ($fh_out, $vcf_file) = @_;
	my $number_sample = @{$vcf_file->{'sample'}};
	for (my $sample_count = 0; $sample_count < $number_sample; $sample_count ++) {
		my $ploidy = 0;
		my @out_haplotypes = ('', '');
		foreach my $site (@{$vcf_file->{'site'}}) {
			my @genotypes;
			my $present_ploidy = 0;
			if ($site->{'GTs'}[$sample_count] =~ /^([\.01])\|([\.01])$/) {
				@genotypes = ($1, $2);
				$present_ploidy = 2;
			}
			elsif ($site->{'GTs'}[$sample_count] =~ /^([\.01])\/([\.01])$/) {
				@genotypes = ('.', '.');
				$present_ploidy = 2;
			}
			elsif ($site->{'GTs'}[$sample_count] =~ /^([01])$/) {
				@genotypes = ($1, '.');
				$present_ploidy = 1;
			}
			elsif ($site->{'GTs'}[$sample_count] eq '.') {@genotypes = ('.', '.');}
			else {
				print ('Warning: Invalid GT "', $site->{'GTs'}[$sample_count], '" at chr', $site->{'CHROM'}, ':',  $site->{'POS'}, ' (',$site->{'ID'} ,'), ', $vcf_file->{'sample'}[$sample_count], "\n");
				@genotypes = ('.', '.');
			}
			if ($ploidy > 0 && $ploidy != $present_ploidy) {
				print ('Haploid diploid mixed: ', $vcf_file->{'sample'}[$sample_count], "\n");
				$ploidy = -1;
				last;
			}
			$ploidy = $present_ploidy;
			foreach my $phase (0, 1) {
				my $gt_here;
				if ($genotypes[$phase] eq '.') {$gt_here = '0';}
				elsif ($genotypes[$phase] == 0) {$gt_here = $site->{'REF'};}
				elsif ($genotypes[$phase] == 1) {$gt_here = $site->{'ALT'};}
				$out_haplotypes[$phase] .= ' ' . $gt_here;
			}
		}
		for (my $phase = 0; $phase < $ploidy; $phase ++) {
			print $fh_out ($vcf_file->{'sample'}[$sample_count], '-', $phase + 1, $out_haplotypes[$phase], "\n");
		}
	}
	return;
}
