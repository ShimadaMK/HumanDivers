#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;

# site filtering parameters +++++++++++++++++++++++++++++++++++++++++++++++++++
my $window_size = 10;
my $mismatch_maximum = 2;
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (@ARGV < 3) {
	print STDERR ("Usage: quality_read.pl BAM_file [BAM_file [...]] Reference_FASTA_file VCF_file\n");
	exit (1);
}
my $in_vcf_filename = pop (@ARGV);
my $in_ref_filename = pop (@ARGV);
my @in_bam_filenames = @ARGV;
if (grep (/^$in_vcf_filename$/, @in_bam_filenames)) {
	print STDERR ("BAM_file and VCF_file cannot be the same file.\n");
	exit (1);
}
if (grep (/^$in_ref_filename$/, @in_bam_filenames)) {
	print STDERR ("BAM_file and Reference_FASTA_file cannot be the same file.\n");
	exit (1);
}
if ($in_vcf_filename eq $in_ref_filename) {
	print STDERR ("VCF_file and Reference_FASTA_file cannot be the same file.\n");
	exit (1);
}
my $out_filename_body;
if ($in_vcf_filename =~ /\/?([^\/]+)\.vcf$/i) {
	$out_filename_body = $1;
}
else {
	$out_filename_body = $in_vcf_filename;
}

# integrate bam files by sample
my $in_bams_by_sample = {};
foreach my $in_bam_filename (@in_bam_filenames) {
	my $sample_id;
	if ($in_bam_filename =~ /\/?([^\/]+)\.(ILLUMINA|SOLID|LS454)\.(low_coverage|exome)\.bam$/i) {
		$sample_id = $1;
	}
	else {
		print STDERR ("The name of '$in_bam_filename' is not well-formed. It is skipped.\n");
		next;
	}
	push (@{$in_bams_by_sample->{$sample_id}}, $in_bam_filename);
}

# load input FASTA file
print ("Reading $in_ref_filename...\n");
my $fh_in = Bio::SeqIO->new (-file => $in_ref_filename, -format => 'Fasta');
my $reference_fa = $fh_in->next_seq;
my $reference_seq = uc ($reference_fa->seq());
my $reference_id = $reference_fa->display_id . ' ' . $reference_fa->desc();
my $reference_chr;
if ($reference_id =~ /chr([\d]{1,2}|X|Y)/ ) {
	$reference_chr = $1;
}
else {
	print STDERR ("No chromosome name in $in_ref_filename\n");
	exit (4);
}
my $reference_start_pos = 1;
my $reference_len = $reference_fa->length();
if ($reference_id =~ /chr$reference_chr:([\d]+)-([\d]+)/ ) {
	if ($2 - $1 + 1 != $reference_len) {
		print STDERR ("The title line in $in_ref_filename says chr$reference_chr:$1-$2, but sequence length is $reference_len\n");
		exit (4);
	}
	$reference_start_pos = $1;
}
my $reference_end_pos = $reference_start_pos + $reference_len - 1;

print ("Reference name=$reference_id\n");
print ("Reference length=$reference_len\n");
print ("chr$reference_chr:$reference_start_pos-$reference_end_pos\n");

# load VCF header
print ("Reading $in_vcf_filename...\n");
my $line_vcf;
my $offset = 9; #number of column before each sample's data in vcf.
open (my $fh_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename: $!";
for ($line_vcf = <$fh_vcf>; $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
if ($line_vcf !~ /^#CHROM/) {
	print STDERR ("No title line in $in_vcf_filename.\n");
	close ($fh_vcf);
	exit (3);
}
chomp ($line_vcf);
my @columns_title = split ("\t", $line_vcf);
my %col_num;
for(my $col = $offset; $col < @columns_title; $col ++) {
	$col_num{$columns_title[$col]} = $col;
}

# load VCF body and divide BAM to phased SAMs
my $pos_vcf_body = tell ($fh_vcf);

# load VCF title columns
my $site_info = [];
while (my $line = <$fh_vcf>) {
	my ($chr, $pos, $id, $ref, $alt) = split ("\t", $line, 6);
	my @ref_alt = (uc ($ref), split (',', uc ($alt)));
	my $valid = 1;
	if ($chr ne $reference_chr || grep (/([^ATGCN]|..|^$)/, @ref_alt)) {
		$valid = 0;
	}
	push (@$site_info, {'CHR' => $chr, 'POS' => $pos, 'ID' => $id, 'REF_ALT' => [@ref_alt], 'valid' => $valid});
}
my $number_of_sites = @$site_info;

# for each sample
my $genotypes = {};
my $number_of_allele = [];
foreach my $sample_id (sort {$a cmp $b} (keys (%$in_bams_by_sample))) {
	print ("--------------------------------\nSample id: $sample_id\n");
	if (!exists ($col_num{$sample_id})) {
		print STDERR ("$sample_id does not exist in $in_vcf_filename. It is skipped.\n");
		next;
	}

	## load VCF body
	seek ($fh_vcf, $pos_vcf_body, 0) or die "Seek error $in_vcf_filename: $!";
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		my $line = <$fh_vcf>;
		if (!$line) {
			print STDERR ("Fatal error occurred while loading $in_vcf_filename\n");
			close ($fh_vcf);
			exit (3);
		}
		if (!$site_info->[$site_count]{'valid'}) {
			next;
		}
		chomp ($line);
		my @columns = split ("\t", $line);
		my $gt = $columns[$col_num{$sample_id}];
		$gt =~ s/:.*$//;
		my @gt = split ('\|', $gt, 2);
		if (grep (/[\D]/, @gt)) {
			print STDERR ("Unacceptable genotype in $in_vcf_filename, pos=$site_info->[$site_count]{'POS'}, ID=$sample_id, GT=$gt.\n",
			              "Only phased genotypes at most diploid are accepted.\n");
			close ($fh_vcf);
			exit (3);
		}
		@{$genotypes->{$sample_id}[$site_count]} = map {$site_info->[$site_count]{'REF_ALT'}[$_]} (@gt);
	}
	## load BAM
	my $quality_by_sample = {};
	foreach my $in_bam_filename (@{$in_bams_by_sample->{$sample_id}}) {
		print ("Processing $in_bam_filename...\n");
		if (!(-e $in_bam_filename && -f $in_bam_filename && -r $in_bam_filename)) {
			print STDERR ("Cannot read $in_bam_filename.\nIt is skipped.\n");
			next;
		}
		my @bam_body = qx/samtools view -F 0x404 $in_bam_filename/;
		my $site_quals = {};
		foreach my $read (@bam_body) {
			chomp ($read);
			my @columns_bam = split ("\t", $read);
			my ($flag, $chr, $pos, $cigar, $seq, $qual_str) = @columns_bam[1, 2, 3, 5, 9, 10];
			if ($chr ne $reference_chr) {
				print STDERR ("$in_ref_filename is chr$reference_chr, but a read of chr$chr appears in BAM.\n");
				close ($fh_vcf);
				exit (2);
			}
			my $strand = $flag & hex('10') ? 'rev' : 'for';
			my @alleles = &split_read_by_cigar (uc ($seq), $cigar);
			my @qual_letters = &split_read_by_cigar ($qual_str, $cigar);
			my $rlen = @alleles;
			my @mismatch = (0) x $rlen;
			for (my $comp_pos = 0; $comp_pos < $rlen; $comp_pos ++) {
				if ($alleles[$comp_pos] ne '') {
					if ($pos + $comp_pos < $reference_start_pos || $reference_end_pos < $pos + $comp_pos) {
						print STDERR ("Position of read in BAM is out of $in_ref_filename.\n");
						close ($fh_vcf);
						exit (2);
					}
					my $reference_allele = substr ($reference_seq, $pos + $comp_pos - $reference_start_pos, 1);
					my $read_allele = substr ($alleles[$comp_pos], 0, 1);
					if ($read_allele ne 'N' && $reference_allele ne 'N' && $read_allele ne $reference_allele) {
						$mismatch[$comp_pos] = 1;
					}
				}
			}
			my $sum_mismatch = 0;
			my $chk_pos_head = 0;
			my $chk_pos_tail = 1 - $window_size;
			while ($chk_pos_tail + $mismatch_maximum < $rlen) {
				while ($chk_pos_head < $rlen && $alleles[$chk_pos_head] eq '') {
					$chk_pos_head ++;
				}
				if ($chk_pos_head < $rlen) {
					$sum_mismatch += $mismatch[$chk_pos_head];
				}
				if ($sum_mismatch > $mismatch_maximum) {
					for (my $rm_pos = $chk_pos_tail; $rm_pos <= $chk_pos_head; $rm_pos ++) {
						if (0 <= $rm_pos && $rm_pos < $rlen) {
							$alleles[$rm_pos] =~ s/./N/g;
							$qual_letters[$rm_pos] =~ s/./!/g;
						}
					}
				}
				while ($chk_pos_tail >= 0 && $alleles[$chk_pos_tail] eq '') {
					$chk_pos_tail ++;
				}
				if ($chk_pos_tail >= 0) {
					$sum_mismatch -= $mismatch[$chk_pos_tail];
				}
				$chk_pos_head ++;
				$chk_pos_tail ++;
			}

			for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
				if (!$site_info->[$site_count]{'valid'} || $site_info->[$site_count]{'CHR'} ne $chr || $site_info->[$site_count]{'POS'} < $pos) {
					next;
				}
				if ($site_info->[$site_count]{'POS'} > $pos + $#alleles) {
					last;
				}
				my $read_allele = $alleles[$site_info->[$site_count]{'POS'} - $pos];
				if (length ($read_allele) < 1) {
					next;
				}
				elsif (length ($read_allele) > 1) {
					substr ($read_allele, 1) = '';
				}
				my $read_allele_quality = ord ($qual_letters[$site_info->[$site_count]{'POS'} - $pos]) - 33;
				push (@{$site_quals->{$site_info->[$site_count]{'POS'}}{$read_allele}{$strand}}, $read_allele_quality);
			}
		}
		foreach my $position (keys (%$site_quals)) {
			foreach my $allele (keys (%{$site_quals->{$position}})) {
				my $qual = 0;
				foreach my $strand (keys (%{$site_quals->{$position}{$allele}})) {
					my @qualities = (sort {$b <=> $a} (@{$site_quals->{$position}{$allele}{$strand}}));
					for (my $i = 0; $i < @qualities; $i ++) {
						$qual += $qualities[$i] / ($i + 1);
					}
				}
				$quality_by_sample->{$position}{$allele} += $qual;
			}
		}
	}
	open (my $fh_out_qual, '>', $sample_id . '.qual') or die "Cannot open $sample_id: $!";
	print ('Outputting ', $sample_id, ".qual...\n");
	print $fh_out_qual ('# ', join (', ', @{$in_bams_by_sample->{$sample_id}}), "\n");
	print $fh_out_qual ("# CHR\tPOS\tVCF_GT\tA\tT\tG\tC\tN\tVCF_allele1\tVCF_allele2\thighest\t2nd\t3rd\t4th\tnewGT\n");
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		if (!$site_info->[$site_count]{'valid'}) {
			next;
		}
		print $fh_out_qual ($site_info->[$site_count]{'CHR'}, "\t", $site_info->[$site_count]{'POS'}, "\t", join ('|', @{$genotypes->{$sample_id}[$site_count]}), "\t");
		my %qual;
		foreach my $allele ('A', 'T', 'G', 'C', 'N') {
			if (exists ($quality_by_sample->{$site_info->[$site_count]{'POS'}}{$allele})) {
				$qual{$allele} = $quality_by_sample->{$site_info->[$site_count]{'POS'}}{$allele};
			}
			else {
				$qual{$allele} = 0;
			}
		}
		print $fh_out_qual (join ("\t", ($qual{'A'}, $qual{'T'}, $qual{'G'}, $qual{'C'}, $qual{'N'})), "\t");
		print $fh_out_qual ($qual{$genotypes->{$sample_id}[$site_count][0]}, "\t");
		if (@{$genotypes->{$sample_id}[$site_count]} > 1) {
			print $fh_out_qual ($qual{$genotypes->{$sample_id}[$site_count][1]});
		}
		print $fh_out_qual ("\t", join ("\t", sort {$b <=> $a} ($qual{'A'}, $qual{'T'}, $qual{'G'}, $qual{'C'})));
		my @new_gt;
		my $warning;
		(@new_gt[0..2], $warning) = &estimate_gt ($genotypes->{$sample_id}[$site_count], \%qual, $site_info->[$site_count]{'REF_ALT'});
		print $fh_out_qual ("\t", @new_gt[0, 2, 1]);
		if ($warning) {
			print $fh_out_qual ("\t", $warning);
		}
		$genotypes->{$sample_id}[$site_count] = [@new_gt];
		foreach my $v (0, 1) {
			if ($new_gt[$v] ne '') {
				if (exists ($number_of_allele->[$site_count]{$new_gt[$v]})) {
					$number_of_allele->[$site_count]{$new_gt[$v]} ++;
				}
				else {
					$number_of_allele->[$site_count]{$new_gt[$v]} = 1;
				}
			}
		}
		print $fh_out_qual ("\n");
	}
	close ($fh_out_qual);

}
close ($fh_vcf);

# Output allele counts
print ("--------------------------------\nOutputting ", $out_filename_body, ".counts...\n");
open (my $fh_out_counts, '>', $out_filename_body . '.counts') or die "Cannot open $out_filename_body.counts: $!";
print $fh_out_counts ("CHR\tPOS\tA\tT\tG\tC\tN\n");
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (!$site_info->[$site_count]{'valid'}) {
		next;
	}
	my %output_num;
	print $fh_out_counts ($reference_chr, "\t", $site_info->[$site_count]{'POS'});
	foreach my $allele ('A', 'T', 'G', 'C', 'N') {
		if (!exists ($number_of_allele->[$site_count]{$allele})) {
			$number_of_allele->[$site_count]{$allele} = 0;
		}
		print $fh_out_counts ("\t", $number_of_allele->[$site_count]{$allele});
	}
	print $fh_out_counts ("\n");
}

# Output common
my @output_samples;
for(my $col = $offset; $col < @columns_title; $col ++) {
	if (exists ($genotypes->{$columns_title[$col]})) {
		push (@output_samples, $columns_title[$col]);
	}
}

# Output to VCF file
print ("Outputting ", $out_filename_body, ".regen.vcf...\n");
open (my $fh_out_vcf, '>', $out_filename_body . '.regen.vcf') or die "Cannot open $out_filename_body.regen.vcf: $!";
print $fh_out_vcf ("##fileformat=VCFv4.1\n");
print $fh_out_vcf ("##source=quality_read.pl\n");
print $fh_out_vcf ("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">\n");
print $fh_out_vcf ("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Count\">\n");
print $fh_out_vcf ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
print $fh_out_vcf ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", join ("\t", @output_samples), "\n");
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (!$site_info->[$site_count]{'valid'}) {
		next;
	}
	my @alleles = @{$site_info->[$site_count]{'REF_ALT'}};
	my $total_allele_number = 0;
	foreach my $allele ('A', 'T', 'G', 'C') {
		$total_allele_number += $number_of_allele->[$site_count]{$allele};
		if ($number_of_allele->[$site_count]{$allele} > 0 && !grep (/^$allele$/, @alleles)) {
			push (@alleles, $allele);
		}
	}
	print $fh_out_vcf ($site_info->[$site_count]{'CHR'}, "\t",
	                   $site_info->[$site_count]{'POS'}, "\t",
	                   $site_info->[$site_count]{'ID'}, "\t",
	                   $alleles[0], "\t", join (',', @alleles[1 .. $#alleles]), "\t.\t.\t",
	                   "AN=", $total_allele_number, ";AC=", join (',', map {$number_of_allele->[$site_count]{$_}} @alleles[1 .. $#alleles]), "\t",
	                   "GT");
	foreach my $sample_id (@output_samples) {
		print $fh_out_vcf ("\t", &allele_to_gt ($genotypes->{$sample_id}[$site_count][0], @alleles));
		if ($genotypes->{$sample_id}[$site_count][2] ne '') {
			print $fh_out_vcf ($genotypes->{$sample_id}[$site_count][2]);
			print $fh_out_vcf (&allele_to_gt ($genotypes->{$sample_id}[$site_count][1], @alleles));
		}
	}
	print $fh_out_vcf ("\n");
}
close ($fh_out_vcf);

# Output to input file of Beagle unphased genotype
print ("Outputting ", $out_filename_body, ".beagle.unphased...\n");
open (my $fh_out_beagle_unphased, '>', $out_filename_body . '.beagle.unphased') or die "Cannot open $out_filename_body.beagle.unphased: $!";
print $fh_out_beagle_unphased ("I\tid");
foreach my $sample_id (@output_samples) {
	print $fh_out_beagle_unphased ("\t", $sample_id, "\t", $sample_id);
}
print $fh_out_beagle_unphased ("\n");
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (!$site_info->[$site_count]{'valid'}) {
		next;
	}
	print $fh_out_beagle_unphased ("M\t", $site_info->[$site_count]{'ID'});
	foreach my $sample_id (@output_samples) {
		my $gen_2nd = ($genotypes->{$sample_id}[$site_count][2] eq '') ?
		               $genotypes->{$sample_id}[$site_count][0] :
		               $genotypes->{$sample_id}[$site_count][1];
		print $fh_out_beagle_unphased ("\t", $genotypes->{$sample_id}[$site_count][0], "\t", $gen_2nd);
	}
	print $fh_out_beagle_unphased ("\n");
}
close ($fh_out_beagle_unphased);

#exit (0);#temporary

# Output to input files of impute2
print ("Outputting ", $out_filename_body, '.impute2.haps and ', $out_filename_body, ".impute2.gens...\n");
open (my $fh_out_impute_haps, '>', $out_filename_body . '.impute2.haps') or die "Cannot open $out_filename_body.impute2.haps: $!";
open (my $fh_out_impute_gens, '>', $out_filename_body . '.impute2.gens') or die "Cannot open $out_filename_body.impute2.gens: $!";
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (!$site_info->[$site_count]{'valid'}) {
		next;
	}
	my @alleles = sort {$number_of_allele->[$site_count]{$b} <=> $number_of_allele->[$site_count]{$a}} ('A', 'T', 'G', 'C');
	if ($number_of_allele->[$site_count]{$alleles[2]} > 0) {
		my $count = ($number_of_allele->[$site_count]{$alleles[3]} > 0) ? 4 : 3;
		print ($count, ' alleles appear at chr', $reference_chr, ':', $site_info->[$site_count]{'POS'}, ". It is skipped.\n");
		next;
	}
	elsif ($number_of_allele->[$site_count]{$alleles[1]} > 0) {
		if ( $alleles[1] eq $site_info->[$site_count]{'REF_ALT'}[0] ||
		    ($alleles[0] ne $site_info->[$site_count]{'REF_ALT'}[0] && $alleles[1] eq $site_info->[$site_count]{'REF_ALT'}[1])) {
			@alleles = @alleles[1, 0];
		}
	}
	elsif ($number_of_allele->[$site_count]{$alleles[0]} > 0) {
		if ($alleles[0] eq $site_info->[$site_count]{'REF_ALT'}[0] || $alleles[0] eq $site_info->[$site_count]{'REF_ALT'}[1]) {
			@alleles = @{$site_info->[$site_count]{'REF_ALT'}}[0, 1];
		}
	}
	else {
		@alleles = @{$site_info->[$site_count]{'REF_ALT'}}[0, 1];
	}
	printf $fh_out_impute_haps ('SNP%d %s %d %s %s', $site_count, $site_info->[$site_count]{'ID'}, $site_info->[$site_count]{'POS'}, @alleles[0, 1]);
	printf $fh_out_impute_gens ('SNP%d %s %d %s %s', $site_count, $site_info->[$site_count]{'ID'}, $site_info->[$site_count]{'POS'}, @alleles[0, 1]);
	foreach my $sample_id (@output_samples) {
		if ($genotypes->{$sample_id}[$site_count][0] eq 'N' &&
		    $genotypes->{$sample_id}[$site_count][1] eq 'N'   ) {
			print $fh_out_impute_haps (' ? ?');
			print $fh_out_impute_gens (' 0.333 0.333 0.333');
		}
		elsif (($genotypes->{$sample_id}[$site_count][0] eq $alleles[0] &&
		        $genotypes->{$sample_id}[$site_count][1] eq 'N'           ) ||
		       ($genotypes->{$sample_id}[$site_count][0] eq 'N'         &&
		        $genotypes->{$sample_id}[$site_count][1] eq $alleles[0]   )) {
			print $fh_out_impute_haps (' ? ?');
			print $fh_out_impute_gens (' 0.5 0.5 0');
		}
		elsif (($genotypes->{$sample_id}[$site_count][0] eq $alleles[1] &&
		        $genotypes->{$sample_id}[$site_count][1] eq 'N'           ) ||
		       ($genotypes->{$sample_id}[$site_count][0] eq 'N'         &&
		        $genotypes->{$sample_id}[$site_count][1] eq $alleles[1]   )) {
			print $fh_out_impute_haps (' ? ?');
			print $fh_out_impute_gens (' 0 0.5 0.5');
		}
		elsif ($genotypes->{$sample_id}[$site_count][0] eq $alleles[0] &&
		       $genotypes->{$sample_id}[$site_count][1] eq $alleles[0]   ) {
			print $fh_out_impute_haps (' 0 0');
			print $fh_out_impute_gens (' 1 0 0');
		}
		elsif ($genotypes->{$sample_id}[$site_count][0] eq $alleles[1] &&
		       $genotypes->{$sample_id}[$site_count][1] eq $alleles[1]   ) {
			print $fh_out_impute_haps (' 1 1');
			print $fh_out_impute_gens (' 0 0 1');
		}
		elsif (($genotypes->{$sample_id}[$site_count][0] eq $alleles[0] &&
		        $genotypes->{$sample_id}[$site_count][1] eq $alleles[1]   ) ||
		       ($genotypes->{$sample_id}[$site_count][0] eq $alleles[1] &&
		        $genotypes->{$sample_id}[$site_count][1] eq $alleles[0]   )) {
			if ($genotypes->{$sample_id}[$site_count][2] eq '|') {
				if ($genotypes->{$sample_id}[$site_count][0] eq $alleles[0]) {
					print $fh_out_impute_haps (' 0 1');
				}
				else {
					print $fh_out_impute_haps (' 1 0');
				}
			}
			else {
				print $fh_out_impute_haps (' 0* 1*');
			}
			print $fh_out_impute_gens (' 0 1 0');
		}
		elsif ($genotypes->{$sample_id}[$site_count][0] eq $alleles[0] &&
		       $genotypes->{$sample_id}[$site_count][1] eq ''          &&
		       $genotypes->{$sample_id}[$site_count][2] eq ''            ) {
			print $fh_out_impute_haps (' 0 -');
			print $fh_out_impute_gens (' 1 0 0');
		}
		elsif ($genotypes->{$sample_id}[$site_count][0] eq $alleles[1] &&
		       $genotypes->{$sample_id}[$site_count][1] eq ''          &&
		       $genotypes->{$sample_id}[$site_count][2] eq ''            ) {
			print $fh_out_impute_haps (' 1 -');
			print $fh_out_impute_gens (' 0 0 1');
		}
		elsif ($genotypes->{$sample_id}[$site_count][0] eq 'N' &&
		       $genotypes->{$sample_id}[$site_count][1] eq ''  &&
		       $genotypes->{$sample_id}[$site_count][2] eq ''    ) {
			print $fh_out_impute_haps (' ? -');
			print $fh_out_impute_gens (' 0.5 0 0.5');
		}
		else {
			print $fh_out_impute_haps (' ? ?');
			print $fh_out_impute_gens (' 0 0 0');
		}
	}
	print $fh_out_impute_haps ("\n");
	print $fh_out_impute_gens ("\n");
}
close ($fh_out_impute_haps);
close ($fh_out_impute_gens);
if ($reference_chr eq 'X') {
	print ("Outputting ", $out_filename_body, ".impute2.samples...\n");
	open ($fh_out_impute_haps, '>', $out_filename_body . '.impute2.samples') or die "Cannot open $out_filename_body.impute2.samples: $!";
	print $fh_out_impute_haps ("ID_1 ID_2 missing sex\n0 0 0 D\n");
	foreach my $sample_id (@output_samples) {
		print $fh_out_impute_haps ($sample_id, ' ', $sample_id, ' 0.0 ');
		my @sex = (1, 1);
		for (my $site_count = 0; ($sex[0] || $sex[1]) && ($site_count < $number_of_sites); $site_count ++) {
			if (!$site_info->[$site_count]{'valid'}) {
				next;
			}
			if ($genotypes->{$sample_id}[$site_count][0] =~ /^[ATGCN]$/ &&
			    $genotypes->{$sample_id}[$site_count][1] =~ /^[ATGCN]$/ &&
			    $genotypes->{$sample_id}[$site_count][2] =~ /^[\|\/]$/) {
				$sex[0] = 0;
			}
			elsif ($genotypes->{$sample_id}[$site_count][0] =~ /^[ATGCN]$/ &&
			       $genotypes->{$sample_id}[$site_count][1] eq '' &&
			       $genotypes->{$sample_id}[$site_count][2] eq '') {
				$sex[1] = 0;
			}
			else {
				@sex = (0, 0);
			}
		}
		if ($sex[0] && !$sex[1]) {
			print $fh_out_impute_haps ("1\n");
		}
		elsif ($sex[1] && !$sex[0]) {
			print $fh_out_impute_haps ("2\n");
		}
		else {
			print $fh_out_impute_haps ("0\n");
		}
	}
}

exit 0;

sub split_read_by_cigar {
	my ($seq, $cigar) = @_;
	my @alleles = ();
	if ($cigar eq '*') {
		return ();
	}
	while ($cigar ne '') {
		my $cigar_letter;
		my $cigar_length;
		if ($cigar =~ /^([\d]+)([MIDNSHPX=])/) {
			$cigar_letter = $2;
			$cigar_length = $1;
			$cigar =~ s/^$cigar_length$cigar_letter//;
		}
		else {
			return ();
		}
		if ($cigar_letter =~ /^[MX=]$/) {
			push (@alleles, split ('', substr ($seq, 0, $cigar_length, '')));
		}
		elsif ($cigar_letter eq 'I') {
			my $ins_seq = substr ($seq, 0, $cigar_length, '');
			if (@alleles > 0) {
				$alleles[$#alleles] .= $ins_seq;
			}
		}
		elsif ($cigar_letter eq 'S') {
			substr ($seq, 0, $cigar_length, '');
		}
		elsif ($cigar_letter =~ /^[DN]$/) {
			push (@alleles, ('') x $cigar_length);
		}
	}
	return @alleles;
}

sub estimate_gt {
	my ($in_gt, $qual, $ref_alt) = @_;
	my @new_gt = ('', '');
	my $phase = '|';
	my @warnings = ();
	my @alleles = sort {$qual->{$b} <=> $qual->{$a}} ('A', 'T', 'G', 'C');
	if (@$in_gt == 1) {
		if ($qual->{$alleles[0]} >= 40) {
			if ($qual->{$alleles[0]} - $qual->{$alleles[1]} > 5) {
				$new_gt[0] = $alleles[0];
			}
			else {
				$new_gt[0] = 'N';
				push (@warnings, '2nd_is_high');
			}
		}
		else {
			$new_gt[0] = 'N';
		}
		$phase = '';
		if (!grep (/$new_gt[0]/, (@$ref_alt, 'N'))) {
			push (@warnings, 'neither_REF_nor_ALT');
		}
		if ($in_gt->[0] ne $new_gt[0]) {
			unshift (@warnings, 'GT_changed');
		}
	}

	elsif (@$in_gt == 2) {
		if ($qual->{$alleles[0]} < 40) {
			@new_gt = ('N', 'N');
		}
		elsif ($qual->{$alleles[0]} >= $qual->{$alleles[1]} * 5) {
			@new_gt = ($alleles[0], $alleles[0]);
		}
		elsif ($qual->{$alleles[1]} < 40 || $qual->{$alleles[1]} - $qual->{$alleles[2]} <= 5) {
			if ($qual->{$alleles[1]} >= 40) {
				push (@warnings, '3rd_is_high');
			}
			if ($in_gt->[0] ne $in_gt->[1]) {
				if ($in_gt->[0] eq $alleles[0]) {
					@new_gt = ($alleles[0], 'N');
				}
				elsif ($in_gt->[1] eq $alleles[0]) {
					@new_gt = ('N', $alleles[0]);
				}
			}
			else {
				@new_gt = ($alleles[0], 'N');
				$phase = '/';
			}
		}
		else {
			if ($in_gt->[0] eq $in_gt->[1]) {
				@new_gt = ($alleles[0], $alleles[1]);
				$phase = '/';
			}
			elsif (($in_gt->[0] eq $alleles[0] && $in_gt->[1] eq $alleles[1]) xor ($in_gt->[1] eq $alleles[0] && $in_gt->[0] eq $alleles[1])) {
				@new_gt = ($in_gt->[0], $in_gt->[1]);
			}
			elsif (($in_gt->[0] eq $alleles[0] && $in_gt->[1] ne $alleles[1]) xor ($in_gt->[0] ne $alleles[0] && $in_gt->[1] eq $alleles[1])) {
				@new_gt = ($alleles[0], $alleles[1]);
			}
			elsif (($in_gt->[0] eq $alleles[1] && $in_gt->[1] ne $alleles[0]) xor ($in_gt->[0] ne $alleles[1] && $in_gt->[1] eq $alleles[0])) {
				@new_gt = ($alleles[1], $alleles[0]);
			}
			else {
				@new_gt = ($alleles[0], $alleles[1]);
				$phase = '/';
			}
		}
		if (!grep (/$new_gt[0]/, (@$ref_alt, 'N')) || !grep (/$new_gt[1]/, (@$ref_alt, 'N'))) {
			push (@warnings, 'neither_REF_nor_ALT');
		}
		if ($in_gt->[0] ne $new_gt[0] || $in_gt->[1] ne $new_gt[1] || $phase eq '/') {
			unshift (@warnings, 'GT_changed');
		}
	}
	else {
		$phase =  ('');
		$warnings[0] =  ('GT_Error');
	}
	return (@new_gt, $phase, join (',', @warnings));
}

sub allele_to_gt {
	my ($allele, @ref_alt) = @_;
	for (my $gt = 0; $gt < @ref_alt; $gt ++) {
		if ($ref_alt[$gt] eq $allele) {
			return $gt;
		}
	}
	return ('.');
}



__END__
output files are...
*.qual      : Quality values (Brad Chapman)

error code
1: argument error
2: BAM file error
3: VCF file error
4: Reference FASTA file error
