#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV != 2 && @ARGV != 4) {
	print STDERR ("Usage: fix_fdi.pl in_FDI_file group_file [1000g_population_file color_file]\n");
	exit (1);
}
my ($in_fdi_filename, $group_filename, $population_filename, $color_filename);
my $color_mode;
if (@ARGV == 4) {
	($in_fdi_filename, $group_filename, $population_filename, $color_filename) = @ARGV[0..3];
	$color_mode = 1;
}
else {
	($in_fdi_filename, $group_filename) = @ARGV[0, 1];
	$color_mode = 0;
}

# check filenames
my @op_check = &search_same_values (@ARGV);
if ($op_check[1] > 0) {
	my @err_filename = ('in_FDI_file', 'group_file', '1000g_population_file', 'color_file');
	print STDERR ($err_filename[$op_check[0]], ' and ', $err_filename[$op_check[1]], " cannot be the same file.\n");
	exit (1);
}

my %color_population;
my %color_indv;
if ($color_mode) {
	# load color file
	open (my $fh_in, '<', $color_filename) or die "Cannot open $color_filename: $!";
	my @color_file = <$fh_in>;
	close ($fh_in);

	foreach my $color_line (@color_file) {
		$color_line =~ s/\r?\n$//;
		if ($color_line =~ /^(.+)\t(.+)$/) {
			$color_population{$1} = &color_for_network ($2);
		}
		else {
			print STDERR ($color_filename, " is invalid.\n");
			print STDERR ('> ', $color_line, "\n");
			exit (3);
		}
	}
	# load population file
	open ($fh_in, '<', $population_filename) or die "Cannot open $population_filename: $!";
	my @population_file = <$fh_in>;
	close ($fh_in);

	foreach my $population_line (@population_file) {
		$population_line =~ s/\r?\n$//;
		if ($population_line =~ /^(.+)\t(HG|NA)(\d{5})$/) {
			my $indv = $3;
			my $population = $1;
			if (exists ($color_population{$population})) {
				$color_indv{$indv} = $color_population{$population};
			}
			else {
				$color_indv{$indv} = 0;
			}
		}
		else {
			print STDERR ($population_filename, " is invalid.\n");
			print STDERR ('> ', $population_line, "\n");
			exit (4);
		}
	}
}

# load group file
open (my $fh_in, '<', $group_filename) or die "Cannot open $group_filename: $!";
my @group_file = <$fh_in>;
close ($fh_in);

my %frequencies_group;
my $members_group = {};
foreach my $group_line (@group_file) {
	if ($group_line =~ /^#/) {
		next;
	}
	$group_line =~ s/\r?\n$//;
	my ($group_id, $members) = split (': ', $group_line);
	my $group_name;
	if ($group_id =~ /^(.+)-(\d+)$/) {
		$group_name = $1;
		$frequencies_group{$group_name} = $2;
	}
	else {
		print STDERR ($group_filename, " is invalid.\n");
		print STDERR ('> ', $group_line, "\n");
		exit (2);
	}
	if ($color_mode) {
		@{$members_group->{$group_name}} = split (',', $members);
	}
}

# load FDI file and rewrite frequencies
open ($fh_in, '<', $in_fdi_filename) or die "Cannot open $in_fdi_filename: $!";
while (my $line = <$fh_in>) {
	if ($line =~ /^TAXON_NAME;([^;\s]+)\s*;/) {
		my $taxon_name = $1;
		if (exists ($frequencies_group{$taxon_name})) {
			$line =~ s/;TAXON_FREQUENCY;1;/;TAXON_FREQUENCY;$frequencies_group{$taxon_name};/;
		}
		if ($color_mode) {
			if (exists ($members_group->{$taxon_name})) {
				my %frequencies_color;
				foreach my $member (@{$members_group->{$taxon_name}}) {
					$member =~ s/[LR]$//;
					if (exists ($color_indv{$member})) {
						$frequencies_color{$color_indv{$member}} ++;
					}
					else {
						$frequencies_color{0} ++;
					}
				}
				my $color_count = 1;
				my $pie_str = '';
				foreach my $color (keys (%frequencies_color)) {
					my $freq = $frequencies_color{$color};
					$pie_str .= "TAXON_COLOR_PIE$color_count;$color;TAXON_PIE_FREQUENCY$color_count;$freq;TAXON_STYLE_PIE$color_count;SOLID;";
					$color_count ++;
				}
				$line =~ s/;(TAXON_COLOR_PIE\d+;\d+;TAXON_PIE_FREQUENCY\d+;\d+;TAXON_STYLE_PIE\d+;[A-Z]*;)+/;$pie_str/;
			}
			else {
				$taxon_name =~ s/[LR]$//;
				my $color = exists ($color_indv{$taxon_name}) ? $color_indv{$taxon_name} : 0;
				$line =~ s/;TAXON_COLOR_PIE1;\d+;/;TAXON_COLOR_PIE1;$color;/;
			}
		}
	}
	print ($line);
}
close ($fh_in);

exit (0);

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

sub color_for_network
{
	my $color_str = $_[0];
	my $color_val;
	if ($color_str !~ /[^\d]/) {
		$color_val = ($color_str < hex ('1000000')) ? $color_str : hex ('ffffff');
	}
	elsif ($color_str =~ /^#([0-9A-F])([0-9A-F])([0-9A-F])$/i) {
		$color_val = (hex ($1) + hex ('100') * hex ($2) + hex ('10000') * hex ($3)) * hex ('11');
	}
	elsif ($color_str =~ /^#([0-9A-F]{2})([0-9A-F]{2})([0-9A-F]{2})$/i) {
		$color_val = hex ($1) + hex ('100') * hex ($2) + hex ('10000') * hex ($3);
	}
	else {
		$color_val = 0;
	}
	return $color_val;
}
