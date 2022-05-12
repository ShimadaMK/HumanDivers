#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::TreeIO;

my %options = &parse_commandline (@ARGV);

my $node_info = {};
my $leaf_stat = {};

# load tree file
print ('Reading ', $options{'--tree'}, "\n");
my $treeio = new Bio::TreeIO('-format' => 'newick', '-file' => $options{'--tree'});
my $tree = $treeio->next_tree;
my $root_node = $tree->get_root_node;

# load tree groups file and calculate number of members
my $tree_groups = (exists ($options{'--groups'})) ? &read_groups ($options{'--groups'}) : {};

# count numberes of members and distances from the root
&set_node_info ($root_node, 0, $tree_groups, $node_info);
&set_leaf_stat ($root_node, $node_info, $leaf_stat, $options{'--place'});

# output
my $out_filename = $options{'--out'} . '.table';
print ("Writing $out_filename...\n");
open (my $fh_out, '>', $out_filename) or die "Cannot open $out_filename: $!";
print $fh_out ("# Branch information\n");
print $fh_out (join ("\t", 'Branch', 'Height', 'Length', 'Size', 'Observed'), "\n");
&output_branch_info ($fh_out, $root_node, $node_info);
print $fh_out ("\n# Leaf information\n");
print $fh_out (join ("\t", 'Leaf', 'Height', 'Observed'), "\n");
&output_leaf_info ($fh_out, $root_node, $node_info);
print $fh_out ("\n# Leaf statistics\n");
print $fh_out (join ("\t", 'Height', 'Total_Size', 'Total_Observed'), "\n");
&output_leaf_stat ($fh_out, $leaf_stat);
close ($fh_out);

$out_filename = $options{'--out'} . '.mega.txt';
print ("Writing $out_filename...\n");
open ($fh_out, '>', $out_filename) or die "Cannot open $out_filename: $!";
&output_mega_group ($fh_out, $root_node, $node_info, $leaf_stat, $options{'--long'});
close ($fh_out);
print ("Done.\n");

exit 0;

sub parse_commandline {
	my @argvs = @_;
	my $argc = @argvs;
	my %options;
	my @valid_selects = ('--tree', '--groups', '--out', '--place', '--long');
	for (my $i = 0; $i < $argc; $i ++) {
		if (grep (/^($argvs[$i])$/, @valid_selects) && ($i + 1 < $argc)) {
			$options{$argvs[$i]} = $argvs[$i + 1];
			$i ++;
		}
		else {
			&usage;
		}
	}
	if (!exists ($options{'--tree'}) || !exists ($options{'--out'})) {
		&usage;
	}
	if (exists ($options{'--place'})) {
		if ($options{'--place'} =~ /\D/) {
			&usage;
		}
	}
	else {
		$options{'--place'} = 5;
	}
	if (exists ($options{'--long'})) {
		if ($options{'--long'} !~ /(\d*\.)?\d+/) {
			&usage;
		}
	}
	else {
		$options{'--long'} = 2;
	}
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
	print ("Usage: check_tree.pl [--tree | --groups | --out] <file> --place <int> --long <float>\n");
	print ("Mandatory arguments:\n");
	print (" --tree   input tree file (Newick / New Hampshire format)\n");
	print (" --out    basename of output files\n");
	print ("Optional argument:\n");
	print (" --groups groups and members file\n");
	print (" --place  number of decimal places to round in statistics table (default=5)\n");
	print (" --long   fraction of distance of OTU to report / mode (default=2)\n");
	exit 1;
}

sub read_groups {
	my $group_filename = $_[0];
	print ("Reading $group_filename\n");
	open (my $fh_in, '<', $group_filename) or die "Cannot open $group_filename: $!";
	my %groups;
	while (my $line = <$fh_in>) {
		if ($line !~ /^#/) {
			chomp ($line);
			my ($group, $member_str) = split (/ *: */, $line);
			$groups{$group} = [split (/ *, */, $member_str)];
		}
	}
	return \%groups;
}

sub set_node_info {
	my ($parent_node, $parent_height, $groups, $node_info) = @_;
	my ($size, $observed) = (0, 0);
	if ($parent_node->is_Leaf) {
		$size = 1;
		$observed = (exists ($groups->{$parent_node->id})) ? @{$groups->{$parent_node->id}} : 1;
	}
	else {
		foreach my $child_node ($parent_node->each_Descendent) {
			$node_info->{$child_node}{'height'} = $parent_height + $child_node->branch_length;
			&set_node_info ($child_node, $node_info->{$child_node}{'height'}, $groups, $node_info);
			$size += $node_info->{$child_node}{'size'};
			$observed += $node_info->{$child_node}{'observed'};
		}
	}
	$node_info->{$parent_node}{'size'} = $size;
	$node_info->{$parent_node}{'observed'} = $observed;
	return;
}

sub set_leaf_stat {
	my ($root_node, $node_info, $leaf_stat, $place) = @_;
	my $height_format = '%.' . $place . 'f';
	foreach my $child_node ($root_node->get_Descendents) {
		if ($child_node->is_Leaf) {
			my $height = sprintf ($height_format, $node_info->{$child_node}{'height'});
			$leaf_stat->{$height}{'size'} += $node_info->{$child_node}{'size'};
			$leaf_stat->{$height}{'observed'} += $node_info->{$child_node}{'observed'};
		}
	}
	return;
}

sub output_branch_info {
	my ($fh_out, $parent_node, $node_info) = @_;
	foreach my $child_node ($parent_node->each_Descendent) {
		my $name = ($child_node->id) ? $child_node->id : '.';
		my $branch_length = $child_node->branch_length;
		my $height = $node_info->{$child_node}{'height'} - $branch_length;
		my $size = $node_info->{$child_node}{'size'};
		my $observed = $node_info->{$child_node}{'observed'};
		print $fh_out (join ("\t", $name, $height, $branch_length, $size, $observed), "\n");
		&output_branch_info ($fh_out, $child_node, $node_info);
	}
	return;
}

sub output_leaf_info {
	my ($fh_out, $parent_node, $node_info) = @_;
	if ($parent_node->is_Leaf) {
		print $fh_out (join ("\t", ($parent_node->id) ? $parent_node->id : '.', @{$node_info->{$parent_node}}{'height', 'observed'}), "\n");
	}
	else {
		foreach my $child_node ($parent_node->each_Descendent) {
			&output_leaf_info ($fh_out, $child_node, $node_info);
		}
	}
	return;
}

sub output_leaf_stat {
	my ($fh_out, $leaf_stat) = @_;
	foreach my $height (sort {$a <=> $b} (keys (%$leaf_stat))) {
		print $fh_out (join ("\t", $height, @{$leaf_stat->{$height}}{'size', 'observed'}), "\n");
	}
	return;
}

sub output_mega_group {
	my ($fh_out, $root_node, $node_info, $leaf_stat, $long) = @_;
	my $peak_height;
	my $max_observed = 0;
	foreach my $height (sort {$a <=> $b} (keys (%$leaf_stat))) {
		if ($max_observed <= $leaf_stat->{$height}{'observed'}) {
			$peak_height = $height;
			$max_observed = $leaf_stat->{$height}{'observed'};
		}
	}
	foreach my $child_node ($root_node->get_Descendents) {
		if ($child_node->is_Leaf && $node_info->{$child_node}{'height'} > $peak_height * $long) {
			print $fh_out (($child_node->id) ? $child_node->id : '.', "=Far\n");
		}
	}
}
