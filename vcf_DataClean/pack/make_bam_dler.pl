#!/usr/bin/perl -w

use strict;
use warnings;

my $uri_prefix = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/';

if (@ARGV != 3) {
	print STDERR ("Usage: make_bam_dler.pl alignment_index_file id_file region\n");
	exit (1);
}
my $in_idx_filename = $ARGV[0];
my $in_id_filename = $ARGV[1];
my $region = $ARGV[2];
if ($in_idx_filename eq $in_id_filename) {
	print STDERR ("alignment_index_file and id_file cannot be the same file.\n");
	exit (1);
}

my $line;
my @bams = ();
open (my $fh_idx, '<', $in_idx_filename) or die "Cannot open $in_idx_filename: $!";
my $bam;
while ($line = <$fh_idx>) {
	chomp ($line);
	($bam, undef) = split ("\t", $line, 2);
	if ($bam =~ /\.bam$/) {
		push (@bams, $bam);
	}
}
close ($fh_idx) or die "Cannot close $in_idx_filename: $!";

open (my $fh_id, '<', $in_id_filename) or die "Cannot open $in_id_filename: $!";
print ("#!/usr/bin/bsh\n");
while (my $id= <$fh_id>) {
	chomp ($id);
	my @uris = grep (/\/$id\.mapped\./, @bams);
	foreach my $uri (@uris) {
		if ($uri =~ /\.(ILLUMINA|SOLID|LS454)\..*\.(exome|low_coverage)\./i) {
			my $platform = uc ($1);
			my $group = lc ($2);
			print ("samtools view -b -o $id.$platform.$group.bam $uri_prefix$uri $region\n");
		}
	}
}
close ($fh_id) or die "Cannot close $in_id_filename: $!";

exit (0);
