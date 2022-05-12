#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV < 4 || 5 < @ARGV) {
	print ("Usage: maf2vcf.pl [in_VCF_file] <in_MAF_file> <chr_name> <reference_OTU_name> <output_OTU_names>\n");
	print ("<output_OTU_names> should be separated by comma without space.\n");
	exit (1);
}

my $in_vcf_filename = '';
if (@ARGV == 5) {
	$in_vcf_filename = shift (@ARGV);
}
my $in_maf_filename = $ARGV[0];
my $maf_chr = uc ($ARGV[1]);
my $reference_sample = $ARGV[2];
my @output_samples = split (',', $ARGV[3]);

# load VCF file
my $vcf_file;
if ($in_vcf_filename) {
	open (my $fh_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename: $!";
	$vcf_file = &read_vcf_file ($fh_vcf);
	close ($fh_vcf);
}
else {
	$vcf_file = &read_vcf_file (*STDIN);
}
if ($vcf_file->{'error'}) {
	print STDERR ($vcf_file->{'error'}, " in the VCF file\n");
	exit (2);
}
my $number_of_sites = @{$vcf_file->{'ID'}};

# load MAF alignment
open (my $fh_in, '<', $in_maf_filename) or die "Cannot open $in_maf_filename: $!";
my $line;
my @columns;
for($line = <$fh_in>; $line =~ /^#/; $line = <$fh_in>) {}
my $alignment = [];
my @align_pos;
my @align_length;
my $align_number = -1;
my $status = 'separator';
while ($line) {
	chomp ($line);
	@columns = split (/ +/, $line, 7);
	if ($status eq 'separator') {
		if ($line =~ /^a /) {
			$align_number ++;
			$status = 'alignment';
		}
	}
	elsif ($status eq 'alignment') {
		if ($line eq '') {
			$status = 'separator';
		}
		elsif ($columns[0] eq 's') {
			if (grep (/^$columns[1]$/, ($reference_sample, @output_samples))) {
				$alignment->[$align_number]{$columns[1]} = uc ($columns[6]);
			}
			if ($columns[1] eq $reference_sample) {
				$align_pos[$align_number] = $columns[2];
				$align_length[$align_number] = $columns[3];
			}
		}
	}
	$line = <$fh_in>;
}
close ($fh_in);
if ($align_number < 0) {
	print STDERR ("No alignment in $in_maf_filename.\n");
	exit (3);
}

# read VCF body before specified chromosome
my $site_count = 0;
while ($site_count < $number_of_sites) {
	if ($vcf_file->{'CHR'}[$site_count] eq $maf_chr) {
		last;
	}
	$site_count ++;
}

# read VCF body of specified chromosome and output VCF
print ("##fileformat=VCFv4.1\n");
print ("##source=maf2vcf.pl\n");
print ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
print ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", join ("\t", @{$vcf_file->{'sample_id'}}, @output_samples), "\n");
my $align_number_now = 0;
while ($site_count < $number_of_sites) {
	if ($vcf_file->{'CHR'}[$site_count] ne $maf_chr) {
		last;
	}
	my ($align_number_end, $letter_start, $letter_end) = &find_align_number ($vcf_file->{'POS'}[$site_count], $vcf_file->{'REF_ALT'}[$site_count][0]);
	if ($align_number_now > $align_number) {
		last;
	}
	my %out_maf;
	my %out_genotype;
	if ($align_number_end == -1) {
		foreach my $sample_id ($reference_sample, @output_samples) {
			$out_maf{$sample_id} = "-";
		}
	}
	elsif ($align_number_now == $align_number_end) {
		foreach my $sample_id ($reference_sample, @output_samples) {
			$out_maf{$sample_id} = exists ($alignment->[$align_number_now]{$sample_id}) ?
			                       substr ($alignment->[$align_number_now]{$sample_id}, $letter_start, $letter_end - $letter_start) : '-';
		}
	}
	else {
		foreach my $sample_id ($reference_sample, @output_samples) {
			$out_maf{$sample_id} = exists ($alignment->[$align_number_now]{$sample_id}) ?
			                       substr ($alignment->[$align_number_now]{$sample_id}, $letter_start) : '-';
		}
		for (my $i = $align_number_now + 1; $i < $align_number_end; $i ++) {
			foreach my $sample_id ($reference_sample, @output_samples) {
				if (exists ($alignment->[$i]{$sample_id})) {
					$out_maf{$sample_id} .= $alignment->[$i]{$sample_id};
				}
			}
		}
		foreach my $sample_id ($reference_sample, @output_samples) {
			if (exists ($alignment->[$align_number_end]{$sample_id})) {
				$out_maf{$sample_id} .= substr ($alignment->[$align_number_end]{$sample_id}, 0, $letter_end);
			}
		}
	}
	foreach my $sample_id ($reference_sample, @output_samples) {
		$out_maf{$sample_id} =~ s/-//g;
		$out_genotype{$sample_id} = -1;
		for (my $gt = 0; $gt < @{$vcf_file->{'REF_ALT'}[$site_count]}; $gt ++) {
			if ($out_maf{$sample_id} eq $vcf_file->{'REF_ALT'}[$site_count][$gt]) {
				$out_genotype{$sample_id} = $gt;
				last;
			}
		}
	}
	my $output_sample_allele_exists = 1;
	foreach my $sample_id (@output_samples) {
		if ($out_genotype{$sample_id} == -1) {
			$output_sample_allele_exists = 0;
			last;
		}
	}
	if ($output_sample_allele_exists) {
		if ($out_genotype{$reference_sample} != 0) {
			print STDERR ('Warning: chr', $vcf_file->{'CHR'}[$site_count], ':', $vcf_file->{'POS'}[$site_count], " REF allele in VCT and $reference_sample allele in MAF mismatch!\n");
		}
		my $ref = shift (@{$vcf_file->{'REF_ALT'}[$site_count]});
		my $alt = join (',', @{$vcf_file->{'REF_ALT'}[$site_count]});
		map {$_ =~ s/:.*$//g} @{$vcf_file->{'GT'}[$site_count]};
		print ($vcf_file->{'CHR'}[$site_count], "\t", $vcf_file->{'POS'}[$site_count], "\t", $vcf_file->{'ID'}[$site_count], "\t");
		print ($ref, "\t", $alt, "\t.\t.\t.\tGT\t", join ("\t", @{$vcf_file->{'GT'}[$site_count]}));
		print ("\t", join ("\t", map {$out_genotype{$_}} @output_samples), "\n");
	}
	$site_count ++;
}

exit 0;

sub read_vcf_file
{
	my $fh_vcf = $_[0];
	my $line_vcf;
	for ($line_vcf = <$fh_vcf>; $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
	if ($line_vcf !~ /^#CHROM/) {
		my %ret_val = ('error' => 'No title line');
		return \%ret_val;
	}
	chomp ($line_vcf);
	my @columns = split ("\t", $line_vcf);
	if (@columns < 10) {
		my %ret_val = ('error' => 'No sample identifiers');
		return \%ret_val;
	}
	my @sample_ids = @columns[9 .. $#columns];
	my @chrs;
	my @poss;
	my @rsids;
	my @alleles;
	my @genotypes;
	while ($line_vcf = <$fh_vcf>) {
		chomp ($line_vcf);
		my ($chr, $pos, $rsid, $ref, $alt, undef, undef, undef, undef, @gt) = split ("\t", $line_vcf);
		my @ref_alt = ($ref, split (',', $alt));
		push (@chrs, $chr);
		push (@poss, $pos);
		push (@rsids, $rsid);
		push (@alleles, [@ref_alt]);
		push (@genotypes, [@gt]);
	}
	my %ret_val = ('sample_id' => \@sample_ids, 'CHR' => \@chrs, 'POS' => \@poss, 'ID' => \@rsids, 'REF_ALT' => \@alleles, 'GT' => \@genotypes, 'error' => '');
	return \%ret_val;
}

sub find_align_number {
	my $pos_start = $_[0];
	my $ref = $_[1];
	my $pos_end = $pos_start + length ($ref) - 1;

	if ($align_pos[$align_number_now] >= $pos_start) {
		return (-1, 0, 0);
	}
	while ($align_pos[$align_number_now] + $align_length[$align_number_now] < $pos_start) {
		$align_number_now ++;
		if ($align_number_now > $align_number) {
			return (-1, 0, 0);
		}
	}
	my $align_number_end = $align_number_now;
	while ($align_pos[$align_number_end] + $align_length[$align_number_end] < $pos_end) {
		$align_number_end ++;
		if ($align_number_end > $align_number) {
			return (-1, 0, 0);
		}
	}

	my $alignment_ref = $alignment->[$align_number_now]{$reference_sample};
	my $letter_start = 0;
	my $start_offset = $pos_start - $align_pos[$align_number_now] - 1;
	if ($alignment_ref =~ /^((?:[^-]-*){$start_offset})/) {
		$letter_start = length ($1);
	}
	$alignment_ref = $alignment->[$align_number_end]{$reference_sample};
	my $letter_end = 0;
	my $end_offset = $pos_end - $align_pos[$align_number_end] - 1;
	if ($alignment_ref =~ /^((?:[^-]-*){$end_offset})/) {
		$letter_end = length ($1) + 1;
	}
	return ($align_number_end, $letter_start, $letter_end);
}
