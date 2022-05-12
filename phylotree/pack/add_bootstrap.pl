#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::TreeIO;

if (@ARGV != 3) {
	print ("Usage: add_bootstrap.pl [in_tree_file] [in_bootstrap_file] [out_tree_file]\n");
	exit;
}

my $filename_intree = $ARGV[0];
my $filename_inboot = $ARGV[1];
my $filename_outtree = $ARGV[2];

my @species;
my %species_order;
my $iterate;
my $pattern;
my $value;
my %bootstrap_values;

my $line;


# load original tree file

my $treeio = new Bio::TreeIO('-format' => 'newick', '-file' => $filename_intree);
my $tree = $treeio->next_tree;
my @nodes = $tree->get_nodes;


# load bootstrap consense file
open (my $fh_in, '<', $filename_inboot) or die "Cannot open $filename_inboot: $!";

## species names
$line = <$fh_in>;
while ($line ne "Species in order: \n") {
	if (!$line) {
		print STDERR ("No species in $filename_inboot.\n");
		close ($fh_in);
		exit;
	}
	$line = <$fh_in>;
}

my $count = 1;
my $name; 
while ($line ne "Sets included in the consensus tree\n") {
	if (!$line) {
		print STDERR ("No bootstrap values in $filename_inboot.\n");
		close ($fh_in);
		exit;
	}
	if ($line =~ /^ *(\d*)\. (.*)$/) {
		if ($1 == $count) {
			$name = $2;
			$name =~ s/ /_/g;
			$species[$count] = $name;
			$species_order{$name} = $count;
			$count ++;
		}
	}
	$line = <$fh_in>;
}
$count --;

print ("$count species!\n");
for (my $i = 1; $i < @species; $i++ ) {
	print ("$i\t$species[$i]\n");
}
foreach my $species (keys (%species_order)) {
	print ("$species\t$species_order{$species}\n");
}

## bootstrap values (included in the consensus tree)
while ($line !~ /Set \(species in order\)/) {
	if (!$line) {
		print STDERR ("No bootstrap values in $filename_inboot.\n");
		close ($fh_in);
		exit;
	}
	$line = <$fh_in>;
}
if ($line =~ / (\d+\.?\d*)$/) {
	$iterate = $1;
}
else {
	print STDERR ("Bootstrap iteration cannot be read in $filename_inboot.\n");
	close ($fh_in);
	exit;
}

print ("Bootstrap $iterate times!\n");
&read_bootstrap_value;

## bootstrap values (NOT included in the consensus tree)
while ($line ne "Sets NOT included in consensus tree:\n") {
	if (!$line) {
		print STDERR ("No bootstrap values of sets NOT included in consensus tree in $filename_inboot.\n");
		close ($fh_in);
		exit;
	}
	$line = <$fh_in>;
}

while ($line !~ /Set \(species in order\)/) {
	if (!$line) {
		print STDERR ("No bootstrap values of sets NOT included in consensus tree in $filename_inboot.\n");
		close ($fh_in);
		exit;
	}
	$line = <$fh_in>;
}
if ($line =~ / (\d+\.?\d*)$/) {
	if ($iterate != $1) {
		print STDERR ("iteration differs; $iterate times sets included in consensus tree, ",
		                             "while $1 times sets NOT included in consensus tree.\n");
		close ($fh_in);
		exit;
	}
}
else {
	print STDERR ("Bootstrap iteration cannot be read in $filename_inboot.\n");
	close ($fh_in);
	exit;
}

print ("Bootstrap $iterate times!\n");
&read_bootstrap_value;
close ($fh_in);

foreach $pattern (keys (%bootstrap_values)) {
	print ("$pattern\t$bootstrap_values{$pattern}\n");
}


# check both files
foreach my $node (@nodes) {
	if ($node->is_Leaf) {
		if (!exists ($species_order{ $node->id })) {
			print STDERR ($node->id, " in $filename_intree does not exist in $filename_inboot.\n");
			exit;
		}
	}
}


#search bootstrap values
my $pattern_init = '.' x $count;
my $boot;
foreach my $node (@nodes) {
	if ($node->is_Leaf) { next; }
	$pattern = $pattern_init;
	foreach my $node_dec ($node->get_Descendents) {
		if ($node_dec->is_Leaf) {
			substr ($pattern, $species_order{$node_dec->id} - 1, 1) = '*';
		}
	}
	$boot = &search_bootstrap_value($pattern);
	print ($node->id, "\t", $pattern, "\t", $boot, "\n");
	$node->bootstrap($boot);
}


# write tree file
$treeio = new Bio::TreeIO('-format' => 'newick', '-file' => ">$filename_outtree");
$treeio->write_tree($tree);

exit (0);


sub read_bootstrap_value
{
	while (($line !~ /^\./) && ($line !~ /\*/)) {
		if (!$line) {
			print STDERR ("No bootstrap values in $filename_inboot.\n");
			close ($fh_in);
			exit;
		}
		$line = <$fh_in>;
	}
	
	while (($line =~ /^\./) || ($line =~ /\*/)) {
		if ($line =~ /^([\.\* ]+[\.\*]) *(\d+\.?\d*)$/) {
			$pattern = $1;
			$value = $2;
			$pattern =~ s/ //g;
			$bootstrap_values{$pattern} = int ($value * 100 / $iterate);
		}
		$line = <$fh_in>;
	}
}

sub search_bootstrap_value
{
	my $pattern = $_[0];
	if (exists ($bootstrap_values{$pattern})) {
		return $bootstrap_values{$pattern};
	}
	$pattern =~ tr/\.\*/*./;
	if (exists ($bootstrap_values{$pattern})) {
		return $bootstrap_values{$pattern};
	}
	return '';
}


