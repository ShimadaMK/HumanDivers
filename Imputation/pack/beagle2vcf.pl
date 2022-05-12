#!/usr/bin/perl -w

use strict;
use warnings;

my $phase_output = 0;
if (grep (/^--phase_output$/, @ARGV)) {
	$phase_output = 1;
	@ARGV = grep (!/^--phase_output$/, @ARGV);
}

if (@ARGV < 3 || 4 < @ARGV ) {
	print STDERR ("Usage: beagle2vcf.pl original_VCF_file modified_VCF_file Beagle_output_file [output_VCF_file]\n");
	exit (1);
}

my ($in_ovcf_filename, $in_mvcf_filename, $in_phased_filename) = @ARGV[0 .. 2];
my $out_vcf_filename;
if (@ARGV == 3) {
	if ($in_mvcf_filename =~ /^(.*)\.vcf$/) {
		$out_vcf_filename = $1 . '.new.vcf';
	}
	else {
		$out_vcf_filename = $in_mvcf_filename . '.new.vcf';
	}
}
else {
	$out_vcf_filename = $ARGV[3];
}

# check filenames
my @op_check = &search_same_values ($in_ovcf_filename, $in_phased_filename, $out_vcf_filename);
if ($op_check[1] > 0) {
	my @err_filename = ('original_VCF_file', 'Beagle_output_file', 'output_VCF_file');
	print STDERR ($err_filename[$op_check[0]], ' and ', $err_filename[$op_check[1]], " cannot be the same file.\n");
	exit (1);
}
@op_check = &search_same_values ($in_mvcf_filename, $in_phased_filename, $out_vcf_filename);
if ($op_check[1] > 0) {
	my @err_filename = ('modified_VCF_file', 'Beagle_output_file', 'output_VCF_file');
	print STDERR ($err_filename[$op_check[0]], ' and ', $err_filename[$op_check[1]], " cannot be the same file.\n");
	exit (1);
}

print ("Original VCF file: $in_ovcf_filename\n");
print ("Modified VCF file: $in_mvcf_filename\n");
print ("Beagle output file: $in_phased_filename\n");
print ("Output VCF file: $out_vcf_filename\n\n");

# read Beagle phased file
my $phased_file = &read_phased_file ($in_phased_filename);
if (@{$phased_file->{'sample_id'}} == 0) {
	print STDERR ("Sample identifier is missing in $in_phased_filename\n");
	exit (4);
}
my %phased_sample_count;
for (my $sample_count = 0; $sample_count < @{$phased_file->{'sample_id'}}; $sample_count += 2) {
	$phased_sample_count{$phased_file->{'sample_id'}[$sample_count]} = $sample_count;
}

# read original VCF file
my $ovcf_file = &read_vcf_file ($in_ovcf_filename);
my %vcf_sample_count;
for (my $sample_count = 0; $sample_count < @{$ovcf_file->{'sample_id'}}; $sample_count ++) {
	$vcf_sample_count{$ovcf_file->{'sample_id'}[$sample_count]} = $sample_count;
}

# compare sample identifiers
my %phase_map;
for (my $sample_count = 0; $sample_count < @{$phased_file->{'sample_id'}}; $sample_count += 2) {
	my $sample_id = $phased_file->{'sample_id'}[$sample_count];
	if ($sample_id ne $phased_file->{'sample_id'}[$sample_count + 1]) {
		print STDERR ('Sample identifiers are confusing at ', $sample_id,
		              ' and ', $phased_file->{'sample_id'}[$sample_count + 1], " in $in_phased_filename\n");
		exit (4);
	}
	if (exists ($vcf_sample_count{$sample_id})) {
		$phase_map{$sample_id} = '';
	}
	else {
		$vcf_sample_count{$sample_id} = -1;
	}
}

# compare phases (only heterozygous sites)
my $phased_site_count = 0;
my $ovcf_site_count = 0;
while ($phased_site_count < @{$phased_file->{'rsid'}}) {
	while ($phased_file->{'rsid'}[$phased_site_count] ne $ovcf_file->{'ID'}[$ovcf_site_count]) {
		$ovcf_site_count ++;
		if ($ovcf_site_count == @{$ovcf_file->{'ID'}}) {
			print STDERR ('The site ', $phased_file->{'rsid'}[$phased_site_count], ' does not exist in ', $in_ovcf_filename, ", or order error\n");
			exit (1);
		}
	}
	for (my $sample_count = 0; $sample_count < @{$phased_file->{'sample_id'}}; $sample_count += 2) {
		my $sample_id = $phased_file->{'sample_id'}[$sample_count];
		if ($vcf_sample_count{$sample_id} < 0) {
			next;
		}
		(my $vcf_genotype, undef) = split (':', $ovcf_file->{'GT'}[$ovcf_site_count][$vcf_sample_count{$sample_id}], 2);
		my $allele_a;
		my $allele_b;
		my $phase;
		if ($vcf_genotype =~ /^(\d+|\.)$/) {
			$allele_a = $1;
			$allele_b = '';
			$phase = '';
		}
		elsif ($vcf_genotype =~ /^(\d+|\.)(\||\/)(\d+|\.)$/) {
			$allele_a = $1;
			$allele_b = $3;
			$phase = $2;
		}
		else {
			print STDERR ("GT field in $in_ovcf_filename is invalid ");
			print STDERR ('at chr', $ovcf_file->{'CHR'}[$ovcf_site_count], ':', $ovcf_file->{'POS'}[$ovcf_site_count], ", $sample_id.\n");
			exit (1);
		}
		$allele_a = ($allele_a ne '.') ? $ovcf_file->{'REF_ALT'}[$ovcf_site_count][$allele_a] : 'N';
		if ($phase ne '') {
			$allele_b = ($allele_b ne '.') ? $ovcf_file->{'REF_ALT'}[$ovcf_site_count][$allele_b] : 'N';
		}
		my $phase_det;
		if ($phase ne '|' || $allele_a eq 'N' || $allele_b eq 'N' || $allele_a eq $allele_b ||
		    $phased_file->{'allele'}[$phased_site_count][$sample_count] eq $phased_file->{'allele'}[$phased_site_count][$sample_count + 1]) {
			$phase_det = '.';
		}
		elsif ($phased_file->{'allele'}[$phased_site_count][$sample_count]     eq $allele_a ||
		       $phased_file->{'allele'}[$phased_site_count][$sample_count + 1] eq $allele_b   ) {
			$phase_det = '=';
		}
		elsif ($phased_file->{'allele'}[$phased_site_count][$sample_count]     eq $allele_b ||
		       $phased_file->{'allele'}[$phased_site_count][$sample_count + 1] eq $allele_a   ) {
			$phase_det = 'X';
		}
		else {
			$phase_det = '.';
		}
		$phase_map{$sample_id} .= $phase_det;

	}
	$phased_site_count ++;
	$ovcf_site_count ++;
}
# determine phases (all sites)
for (my $sample_count = 0; $sample_count < @{$phased_file->{'sample_id'}}; $sample_count += 2) {
	my $sample_id = $phased_file->{'sample_id'}[$sample_count];
	if ($vcf_sample_count{$sample_id} < 0) {
		next;
	}
	foreach my $phase_det ('=', 'X') {
		if ($phase_map{$sample_id} =~ /^(\.+)$phase_det/) {
			my $rep_str = $phase_det x length ($1);
			$phase_map{$sample_id} =~ s/^(\.+)$phase_det/$rep_str$phase_det/;
		}
		if ($phase_map{$sample_id} =~ /$phase_det(\.+)$/) {
			my $rep_str = $phase_det x length ($1);
			$phase_map{$sample_id} =~ s/$phase_det(\.+)$/$phase_det$rep_str/;
		}
		while ($phase_map{$sample_id} =~ /$phase_det(\.+)$phase_det/) {
			my $rep_str = $phase_det x length ($1);
			$phase_map{$sample_id} =~ s/$phase_det(\.+)$phase_det/$phase_det$rep_str$phase_det/;
		}
	}
	$phase_map{$sample_id} =~ s/([^=])=/$1-/g;
	$phase_map{$sample_id} =~ s/=([^=])/-$1/g;
	$phase_map{$sample_id} =~ s/([^X])X/$1x/g;
	$phase_map{$sample_id} =~ s/X([^X])/x$1/g;
	if ($phase_map{$sample_id} =~ /^\.*$/) {
		$phase_map{$sample_id} =~ tr/\./o/;
	}
}

# output phase map
if ($phase_output) {
	open (my $fh_phase, '>', $in_phased_filename . '.phase') or die "Cannot open $in_phased_filename.phase': $!";
	my $out_flag = 1;
	foreach my $sample_id (@{$phased_file->{'sample_id'}}) {
		if ($out_flag){
			print $fh_phase ($sample_id, "\t", $phase_map{$sample_id}, "\n");
		}
		$out_flag = !$out_flag;
	}
	close ($fh_phase);
}

# read modified VCF file
my $mvcf_file = &read_vcf_file ($in_mvcf_filename);

# compare site identifiers
if (@{$phased_file->{'rsid'}} != @{$mvcf_file->{'ID'}}) {
	print STDERR ("The number of sites differ between $in_mvcf_filename and $in_phased_filename.\n");
	exit (4);
}
for (my $site_count = 0; $site_count < @{$phased_file->{'rsid'}}; $site_count ++) {
	if ($phased_file->{'rsid'}[$site_count] ne $mvcf_file->{'ID'}[$site_count]) {
		print STDERR ("Site IDs differ between $in_mvcf_filename and $in_phased_filename in the site No. $site_count.\n");
		exit (4);
	}
}
# check sites between the Beagle output file and the modified VCF file
my @replace_counts = (0, 0);
for (my $sample_count = 0; $sample_count < @{$mvcf_file->{'sample_id'}}; $sample_count ++) {
	my $sample_id = $mvcf_file->{'sample_id'}[$sample_count];
	if (!exists ($phased_sample_count{$sample_id}) || $vcf_sample_count{$sample_id} < 0) {
		next;
	}
	for (my $site_count = 0; $site_count < @{$phased_file->{'rsid'}}; $site_count ++) {
		(my $vcf_genotype, undef) = split (':', $mvcf_file->{'GT'}[$site_count][$sample_count], 2);
		$mvcf_file->{'GT'}[$site_count][$sample_count] = $vcf_genotype;
		if ($vcf_genotype !~ /^[\d+\.]([\|\/][\d+\.])?$/) {
			print ("$sample_id ", $mvcf_file->{'CHR'}[$site_count], ':', $mvcf_file->{'POS'}[$site_count]);
			print (" Invalid genotype. It is skipped.\n");
			next;
		}
		if ($vcf_genotype !~ /[\.\/]/) {
			next;
		}
		my @phased_genotype = map {
			&allele_to_gt ($_, @{$mvcf_file->{'REF_ALT'}[$site_count]})
		} @{$phased_file->{'allele'}[$site_count]}[$phased_sample_count{$sample_id}, $phased_sample_count{$sample_id} + 1];
		if (grep (/^\.$/, @phased_genotype)) {
			print ("$sample_id ", $mvcf_file->{'CHR'}[$site_count], ':', $mvcf_file->{'POS'}[$site_count]);
			print (" Allele in $in_phased_filename does not exist in REF nor ALT in $in_mvcf_filename. It is skipped.\n");
			$replace_counts[1] ++;
			next;
		}
		my $new_genotype;
		if ($vcf_genotype eq '.') {
			$new_genotype = $phased_genotype[0];
			if ($new_genotype ne $phased_genotype[1]) {
				print ("$in_mvcf_filename is hemi, but $in_phased_filename is hetero. An arbitrary phase has been chosen: ");
			}
		}
		elsif ($phased_genotype[0] == $phased_genotype[1]) {
			$new_genotype = $phased_genotype[0] . '|' . $phased_genotype[0];
		}
		elsif (substr ($phase_map{$sample_id}, $site_count, 1) =~ /^[=o]$/) {
			$new_genotype = $phased_genotype[0] . '|' . $phased_genotype[1];
		}
		elsif (substr ($phase_map{$sample_id}, $site_count, 1) eq 'X') {
			$new_genotype = $phased_genotype[1] . '|' . $phased_genotype[0];
		}
		else {
			print ("$sample_id ", $mvcf_file->{'CHR'}[$site_count], ':', $mvcf_file->{'POS'}[$site_count]);
			print (" Genotype in $in_mvcf_filename is $vcf_genotype, and genotype in $in_phased_filename is @phased_genotype, but the phase is unknown.");
			print (" It is skipped.\n");
			$replace_counts[1] ++;
			next;
		}
		$mvcf_file->{'GT'}[$site_count][$sample_count] = $new_genotype;
		print ("$sample_id ", $mvcf_file->{'CHR'}[$site_count], ':', $mvcf_file->{'POS'}[$site_count], " $vcf_genotype -> $new_genotype\n");
		$replace_counts[0] ++;
	}
}

# output VCF
open (my $fh_out_vcf, '>', $out_vcf_filename) or die "Cannot open $out_vcf_filename: $!";
print $fh_out_vcf ("##fileformat=VCFv4.1\n");
print $fh_out_vcf ("##source=beagle2vcf.pl\n");
print $fh_out_vcf ("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">\n");
print $fh_out_vcf ("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Count\">\n");
print $fh_out_vcf ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
print $fh_out_vcf ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", join ("\t", @{$mvcf_file->{'sample_id'}}), "\n");
for (my $site_count = 0; $site_count < @{$mvcf_file->{'ID'}}; $site_count ++) {
	my @alleles = @{$mvcf_file->{'REF_ALT'}[$site_count]};
	my %total_allele_count;
	foreach my $gt (0 .. $#alleles) {
		$total_allele_count{$gt} = 0;
	}
	$total_allele_count{'.'} = 0;
	foreach my $gt (@{$mvcf_file->{'GT'}[$site_count]}) {
		if ($gt =~ /^(\d+|\.)$/) {
			$total_allele_count{$1} ++;
		}
		elsif ($gt =~ /^(\d+|\.)[\|\/](\d+|\.)$/) {
			$total_allele_count{$1} ++;
			$total_allele_count{$2} ++;
		}
	}
	my $an = 0;
	foreach my $gt (0 .. $#alleles) {
		$an += $total_allele_count{$gt};
	}
	my $ac = join (',', map {$total_allele_count{$_}} (1 .. $#alleles));
	my $ref = shift (@alleles);
	my $alt = join (',', @alleles);
	print $fh_out_vcf ($mvcf_file->{'CHR'}[$site_count], "\t",
	                   $mvcf_file->{'POS'}[$site_count], "\t",
	                   $mvcf_file->{'ID'}[$site_count], "\t",
	                   $ref, "\t", $alt, "\t.\t.\tAN=$an;AC=$ac\tGT\t",
	                   join ("\t", @{$mvcf_file->{'GT'}[$site_count]}), "\n");
}
close ($fh_out_vcf);
print ("\n$replace_counts[0] missing or unphased genotypes have been replaced.\n");
print ("$replace_counts[1] missing or unphased genotypes have not been replaced because of phase uncertainties.\n");

exit ;

sub search_same_values
{
	my @in_v = @_;
	my $num_v = @in_v;
	if ($num_v <= 1) {
		return (0, 0);
	}
	my $i1;
	my $i2;
	for ($i1 = 0; $i1 < $num_v; $i1 ++) {
		for ($i2 = 0; $i2 < $i1; $i2 ++) {
			if ($in_v[$i1] eq $in_v[$i2]) {
				return ($i2, $i1);
			}
		}
	}
	return (0, 0);
}

sub read_vcf_file
{
	my $in_vcf_filename = $_[0];
	open (my $fh_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename: $!";
	my $line_vcf;
	for ($line_vcf = <$fh_vcf>; $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
	if ($line_vcf !~ /^#CHROM/) {
		print STDERR ("No title line in $in_vcf_filename.\n");
		close ($fh_vcf);
		exit (5);
	}
	chomp ($line_vcf);
	my @columns = split ("\t", $line_vcf);
	if (@columns < 10) {
		print STDERR ("No sample identifiers in $in_vcf_filename.\n");
		close ($fh_vcf);
		exit (5);
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
	close ($fh_vcf);
	my %ret_val = ('sample_id' => \@sample_ids, 'CHR' => \@chrs, 'POS' => \@poss, 'ID' => \@rsids, 'REF_ALT' => \@alleles, 'GT' => \@genotypes);
	return \%ret_val;
}

sub read_phased_file
{
	my $in_phased_filename = $_[0];
	open (my $fh_phased, '<', $in_phased_filename) or die "Cannot open $in_phased_filename: $!";
	my @phased_file = <$fh_phased>;
	close ($fh_phased);
	if (@phased_file == 0) {
		print STDERR ("$in_phased_filename is empty.\n");
		exit (2);
	}
	my @sample_ids = ();
	my @rsids;
	my $alleles = [];
	if ($phased_file[0] =~ /^I[ \t]/) {
		my $sample_id_str = shift (@phased_file);
		chomp ($sample_id_str);
		@sample_ids = split (/[ \t]/, $sample_id_str);
		splice (@sample_ids, 0, 2);
	}
	my @phased_file_body = grep (/^M[ \t]/, @phased_file);
	foreach my $line_str (@phased_file_body) {
		chomp ($line_str);
		$line_str =~ s/^M[ \t]//;
		my @columns = split (/[ \t]/, $line_str);
		push (@rsids, shift (@columns));
		push (@$alleles, [@columns]);
	}
	my %ret_val = ('sample_id' => \@sample_ids, 'rsid' => \@rsids, 'allele' => $alleles);
	return \%ret_val;
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
error code
1: argument error
2: original VCF file error
3: modified VCF file error
4: Beagle phased file error
5: VCF file error



文字列
$phase_map{$sample_id}

Beagle output ファイルの phase と original VCF ファイルの phase が同順か逆順かを1文字/1サイトで記述。
この値に基づいて、Beagleの出力がheteroであるサイトを使用するかどうかを決定する。
= は同順であり、Beagleの結果を使う。。
- は同順であるが、隣のheteroのサイトが逆順なので、Beagle の結果を使わない。
X は逆順であり、Beagleの結果を使う。
x は逆順であるが、隣のheteroのサイトが同順なので、Beagle の結果を使わない。
o は領域全体にわたってhomoかhemiなので（同順でも逆順でもいい）、Beagle の順序をそのまま使う。
. は両隣のheteroのサイトの phase が合わないので、Beagle の結果を使わない。

なおBeagleの出力がhomoであるサイトは上記にかかわらず使用する。
