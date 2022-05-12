#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
my $window_size = 10;
my $mismatch_maximum = 2;
my $quality_minimum = 0; # Minimum quality value of valid site. ++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
my $quality_minimum_letter = chr ($quality_minimum + 33);

if (@ARGV < 3) {
	print STDERR ("Usage: phase_read.pl BAM_file [BAM_file [...]] Reference_FASTA_file VCF_file\n");
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

# integrate bam files by sample
my $in_bams_by_sample = {};
foreach my $in_bam_filename (@in_bam_filenames) {
	## load VCF body
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
if ($reference_id =~ /chr$reference_chr:([\d]*)-([\d]*)/ ) {
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
my $pos_title_line = 0;
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
my @suffixes = ('-1', '-2', '-H', '-U', '-N');

foreach my $sample_id (sort {$a cmp $b} (keys (%$in_bams_by_sample))) {
	print ("--------------------------------\nSample id: $sample_id\n");
	if (!exists ($col_num{$sample_id})) {
		print STDERR ("$sample_id does not exist in $in_vcf_filename. It is skipped.\n");
		next;
	}

	## load VCF body
	seek ($fh_vcf, $pos_vcf_body, 0) or die "Seek error $in_vcf_filename: $!";
	my $genotypes = [];
	while (my $line = <$fh_vcf>) {
		chomp ($line);
		my @columns = split ("\t", $line);
		my ($chr, $pos_start, $ref, $alt, $gt) = @columns[0, 1, 3, 4, $col_num{$sample_id}];
		my @ref_alt = (uc ($ref), split (',', uc ($alt)));
		if (grep (/[^ATGCN]/, @ref_alt)) {
			next;
		}
		if (grep (/../, @ref_alt)) {
			next;
		}
		$gt =~ s/:.*$//;
		my @gt = split ('\|', $gt, 2);
		if (grep (/[\D]/, @gt)) {
			print STDERR ("Unacceptable genotype in $in_vcf_filename, pos=$pos_start, ID=$sample_id, GT=$gt.\n",
			              "Only phased genotypes at most diploid are accepted.\n");
			exit (3);
		}
		@gt = map {$ref_alt[$_]} (@gt);
		push (@$genotypes, {'chr' => $chr, 'pos_start' => $pos_start, 'gt' => \@gt});
	}
	## load BAM
	my $quality_by_sample = {};
	foreach my $in_bam_filename (@{$in_bams_by_sample->{$sample_id}}) {
		print ("\nProcessing $in_bam_filename...\n");
		my $out_sam_suffix;
		if ($in_bam_filename =~ /\/?$sample_id\.(ILLUMINA|SOLID|LS454)\.(low_coverage|exome)\.bam$/i) {
			$out_sam_suffix = ".$1.$2.sam";
		}
		if (!(-e $in_bam_filename && -f $in_bam_filename && -r $in_bam_filename)) {
			print STDERR ("Cannot read $in_bam_filename.\nIt is skipped.\n");
			next;
		}
		my @bam = qx/samtools view -h -F 0x4 $in_bam_filename/;
		my @bam_header = grep (/^@/, @bam);
		my @bam_body = grep (/^[^@]/, @bam);
		my %fh_out;
		foreach my $suffix (@suffixes) {
			open ($fh_out{$suffix}, '>', $sample_id . $suffix . $out_sam_suffix) or die "Cannot open $sample_id$suffix$out_sam_suffix: $!";
			print {$fh_out{$suffix}} (@bam_header);
		}
		my $site_quals = {};
		foreach my $read (@bam_body) {
			chomp ($read);
			my @columns_bam = split ("\t", $read);
			my ($flag, $chr, $pos, $cigar, $seq, $qual_str) = @columns_bam[1, 2, 3, 5, 9, 10];
			if ($chr ne $reference_chr) {
				foreach my $suffix (@suffixes) {
					close ($fh_out{$suffix});
					print STDERR ("Output $sample_id$suffix$out_sam_suffix interrupted.\n");
				}
				print STDERR ("$in_ref_filename is chr$reference_chr, but a read of chr$chr appears in BAM.\n");
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
						foreach my $suffix (@suffixes) {
							close ($fh_out{$suffix});
							print STDERR ("Output $sample_id$suffix$out_sam_suffix interrupted.\n");
						}
						print STDERR ("Position of read in BAM is out of $in_ref_filename.\n");
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

			for (my $p = 0; $p < $rlen; $p ++) {
				my $pos_len = length ($alleles[$p]);
				for (my $q = 0; $q < $pos_len; $q ++) {
					if (substr ($qual_letters[$p], $q, 1) lt $quality_minimum_letter) {
						substr ($alleles[$p], $q, 1) = 'N';
						substr ($qual_letters[$p], $q, 1) = '!';
					}
				}
			}

			my @phase_match = (-1, -1);
			foreach my $genotype (@$genotypes) {
				if ($genotype->{'chr'} ne $chr || $genotype->{'pos_start'} < $pos) {
					next;
				}
				if ($genotype->{'pos_start'} > $pos + $#alleles) {
					last;
				}
				my $read_allele = $alleles[$genotype->{'pos_start'} - $pos];
				if (length ($read_allele) < 1) {
					next;
				}
				elsif (length ($read_allele) > 1) {
					substr ($read_allele, 1) = '';
				}
				my $read_allele_quality = ord ($qual_letters[$genotype->{'pos_start'} - $pos]) - 33;
				push (@{$site_quals->{$genotype->{'pos_start'}}{$read_allele}{$strand}}, $read_allele_quality);
				$read_allele =~ s/N/.?/g;
				for (my $phase = 0; $phase < @{$genotype->{'gt'}}; $phase ++) {
					if ($genotype->{'gt'}[$phase] =~ /^$read_allele$/) {
						$phase_match[$phase] = 1;
					}
					else {
						$phase_match[$phase] = 0;
					}
				}
			}
			my $suffix = &det_phase (@phase_match);
			$qual_str = join ('', @qual_letters);
			if ($suffix eq '-H') {
				$qual_str = pack ('C*', map {int (($_ - 32) / 2) + 33} (unpack ('C*', $qual_str)));
			}
			$seq = join ('', @alleles);
			$cigar =~ s/([\d]+)S//g;
			@columns_bam[5, 9, 10] = ($cigar, $seq, $qual_str);
			print {$fh_out{$suffix}} (join ("\t", @columns_bam[0..10]), "\n");
		}
		foreach my $suffix (@suffixes) {
			close ($fh_out{$suffix});
			print ("Output $sample_id$suffix$out_sam_suffix\n");
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
	print $fh_out_qual ('# ', join (', ', @{$in_bams_by_sample->{$sample_id}}), "\n");
	print $fh_out_qual ("# CHR\tPOS\tVCF_GT\tA\tT\tG\tC\tN\tVCF_allele1\tVCF_allele2\thighest\t2nd\t3rd\t4th\n");
	foreach my $genotype (@$genotypes) {
		print $fh_out_qual ($genotype->{'chr'}, "\t", $genotype->{'pos_start'}, "\t", join ('|', @{$genotype->{'gt'}}), "\t");
		my %qual;
		foreach my $allele ('A', 'T', 'G', 'C', 'N') {
			if (exists ($quality_by_sample->{$genotype->{'pos_start'}}{$allele})) {
				$qual{$allele} = $quality_by_sample->{$genotype->{'pos_start'}}{$allele};
			}
			else {
				$qual{$allele} = 0;
			}
		}
		print $fh_out_qual (join ("\t", ($qual{'A'}, $qual{'T'}, $qual{'G'}, $qual{'C'}, $qual{'N'})), "\t");
		print $fh_out_qual ($qual{$genotype->{'gt'}[0]}, "\t");
		if (@{$genotype->{'gt'}} > 1) {
			print $fh_out_qual ($qual{$genotype->{'gt'}[1]});
		}
		print $fh_out_qual ("\t", join ("\t", sort {$b <=> $a} ($qual{'A'}, $qual{'T'}, $qual{'G'}, $qual{'C'})), "\n");
	}
	close ($fh_out_qual);

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

sub det_phase {
	my @in = @_;
	my $ret = '';
	if ($in[0] > 0 && $in[1] <= 0) {
		$ret = '-1';
	}
	elsif ($in[0] == 0 && $in[1] > 0) {
		$ret = '-2';
	}
	elsif ($in[0] > 0 && $in[1] > 0) {
		$ret = '-H';
	}
	elsif ($in[0] < 0 && $in[1] < 0) {
		$ret = '-N';
	}
	elsif ($in[0] == 0 && $in[1] <= 0) {
		$ret = '-U';
	}
	return $ret;
}

__END__
output files are...
*-1.*.*.sam : reads matching to phase-1
*-2.*.*.sam : reads matching to phase-2
*-H.*.*.sam : reads matching to both phases, quality values are divided by 2
*-N.*.*.sam : reads including no polymorphic sites
*-U.*.*.sam : reads NOT matching to either phase
*.qual      : Quality values (Brad Chapman)

error code
1: argument error
2: BAM file error
3: VCF file error
4: Reference FASTA file error
