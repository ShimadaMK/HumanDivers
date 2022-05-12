#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV != 2) {
	print ("Usage: rmsk_pos.pl [in_VCF_file] [in_rmsk_file]\n");
	exit;
}

my $filename_vcf = $ARGV[0];
my $filename_rmsk = $ARGV[1];

my $line;
my @columns;

# load RepeatMask
open (my $fh_in, '<', $filename_rmsk) or die "Cannot open $filename_rmsk: $!";
$line = <$fh_in>;
if (!$line) {
	print STDERR ("$filename_rmsk is empty.\n");
	close ($fh_in);
	exit;
}
if ($line !~ /^#/) {
	print STDERR ("No title line in $filename_rmsk.\n");
	close ($fh_in);
	exit;
}

$line =~ s/^#//;
chomp ($line);
@columns = split ("\t", $line);
my %col_num = ('genoname' => -1, 'genostart' => -1, 'genoend' => -1, 'repclass' => -1);
for (my $col = 0; $col < @columns; $col ++) {
	if ($columns[$col] eq 'genoName') {
		$col_num{'genoname'} = $col;
	}
	elsif ($columns[$col] eq 'genoStart') {
		$col_num{'genostart'} = $col;
	}
	elsif ($columns[$col] eq 'genoEnd') {
		$col_num{'genoend'} = $col;
	}
	elsif ($columns[$col] eq 'repClass') {
		$col_num{'repclass'} = $col;
	}

}
foreach my $col_name (keys (%col_num)) {
	if ($col_num{$col_name} < 0) {
		print STDERR ("$col_name column is missing in $filename_rmsk.\n");
		close ($fh_in);
		exit;
	}
}

my @repeat_chr;
my @repeat_start;
my @repeat_end;
while ($line = <$fh_in>) {
	chomp ($line);
	@columns = split ("\t", $line);
	if ($columns[$col_num{'repclass'}] eq 'Simple_repeat') {
		push (@repeat_chr,   $columns[$col_num{'genoname'}]);
		push (@repeat_start, $columns[$col_num{'genostart'}]);
		push (@repeat_end,   $columns[$col_num{'genoend'}]);
	}
}
close ($fh_in);

my $num_of_repeats = @repeat_start;
for (my $i = 0; $i < $num_of_repeats; $i ++) {
	$repeat_chr[$i] =~ s/chr//;
}

# load vcf header
open (my $fh_vcf, '<', $filename_vcf) or die "Cannot open $filename_vcf: $!";
for ($line = <$fh_vcf>; $line =~ /^##/; $line = <$fh_vcf>) {}
if ($line !~ /^#CHROM/) {
	print STDERR ("No title line in $filename_vcf.\n");
	close ($fh_vcf);
	exit;
}

#load vcf body and output masked chr & pos
print ("#CHROM\tPOS\n");
my $chr_vcf;
my $pos_vcf;
while ($line = <$fh_vcf>) {
	chomp ($line);
	($chr_vcf, $pos_vcf, undef) = split ("\t", $line, 3);
	for (my $i = 0; $i < $num_of_repeats; $i ++) {
		if ($chr_vcf eq $repeat_chr[$i] && $repeat_start[$i] < $pos_vcf && $pos_vcf <= $repeat_end[$i]) {
			print ($chr_vcf, "\t", $pos_vcf, "\n");
		}
	}
}
close ($fh_vcf);
