#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::TreeIO;

my %options = &parse_commandline (@ARGV);

my $node_info = {};

# load tree file
print ('Reading ', $options{'--tree'}, "\n");
my $treeio = new Bio::TreeIO('-format' => 'newick', '-file' => $options{'--tree'});
my $tree = $treeio->next_tree;
my $root_node = $tree->get_root_node;

# load tree groups file and calculate number of members
my $tree_groups = (exists ($options{'--tree_groups'})) ? &read_groups ($options{'--tree_groups'}) : {};
&count_members ($root_node, $tree_groups, $node_info);

# make clusters
print ("Making clusters...\n");
my $clusters = [];
&make_clusters ($options{'--size'}, $root_node, $clusters, $node_info);
print (scalar (@$clusters), " cluster(s) found.\n\n");

print ("Searching isolated branches...\n");
my $branches = [];
&set_interval_to_cluster (&interval_to_desc_cluster ($root_node, $clusters, $node_info), $root_node, $node_info);
&make_isolated_clusters ($options{'--interval'}, $root_node, $clusters, $branches, $node_info);
print (scalar (@$branches), " isolated branch(es) found.\n\n");

print ("Choosing representatives...\n");
push (@$clusters, @$branches);
undef ($branches);
foreach my $cluster (@$clusters) {
	$cluster->{'representative'} = &choose_representative ($cluster->{'root_node'}, $node_info);
}
print (scalar (@$clusters), " representative(s) chosen.\n\n");

print ("Reducing clusters...\n");
my $number_of_del_clusters = &reduce_clusters ($options{'--neighbor'}, $options{'--reps'}, $clusters, $node_info);
print ($number_of_del_clusters, " cluster(s) deleted,\n");
print (@$clusters - $number_of_del_clusters, " cluster(s) retained.\n\n");

# load Network groups file
my $nw_groups = (exists ($options{'--nw_groups'})) ? &read_groups ($options{'--nw_groups'}) : {};
my %indv_group;
foreach my $group (keys (%$nw_groups)) {
	my $group_nonum = $group;
	$group_nonum =~ s/-\d+$//;
	foreach my $id (@{$nw_groups->{$group}}) {
		$indv_group{$id} = $group_nonum;
	}
}

# output
my $out_nw_filename = $options{'--out'} . '.network.id';
print ("Writing clusters list to $out_nw_filename...\n");
open (my $fh_out_nw, '>', $out_nw_filename) or die "Cannot open $out_nw_filename: $!";
foreach my $cluster (grep {$_->{'valid'}} @$clusters) {
	my @output_otus = (exists ($tree_groups->{$cluster->{'representative'}->id})) ? @{$tree_groups->{$cluster->{'representative'}->id}} : ($cluster->{'representative'}->id);
	foreach my $otu (@output_otus) {
		$otu =~ s/^(NA|HG)(\d{5})\-1$/$2L/;
		$otu =~ s/^(NA|HG)(\d{5})\-2$/$2R/;
		$otu = substr ($otu, 0, 6);
		if (exists ($indv_group{$otu})) {
			$otu = $indv_group{$otu};
		}
	}
	print $fh_out_nw (join ("\n", @output_otus), "\n");
}
close ($fh_out_nw);
print ("Done.\n\n");

my $out_list_filename = $options{'--out'} . '.clusters';
print ("Writing clusters list to $out_list_filename...\n");
open (my $fh_out_list, '>', $out_list_filename) or die "Cannot open $out_list_filename: $!";
foreach my $cluster (@$clusters) {
	print $fh_out_list (($cluster->{'valid'}) ? 'Cluster: ' : 'Del_cluster: ');
	print $fh_out_list ($cluster->{'representative'}->id, '; ');
	print $fh_out_list (join (', ', map {$_->id} grep {$_->is_Leaf} ($cluster->{'root_node'}->get_Descendents, $cluster->{'root_node'})), "\n");
}
close ($fh_out_list);
print ("Done.\n\n");

my $out_mega_filename = $options{'--out'} . '.mega.txt';
print ("Writing MEGA group file to $out_mega_filename...\n");
open (my $fh_out_mega, '>', $out_mega_filename) or die "Cannot open $out_mega_filename: $!";
foreach my $cluster (@$clusters) {
	my ($rep_group, $member_group) = $cluster->{'valid'} ? ('R', 'Cluster') : ('Del_R', 'Del_cluster');
	print $fh_out_mega ($cluster->{'representative'}->id, "=$rep_group\n");
	foreach my $node (grep {$_->is_Leaf} ($cluster->{'root_node'}->get_Descendents, $cluster->{'root_node'})) {
		print $fh_out_mega ($node->id, "=$member_group\n");
	}
}
close ($fh_out_mega);
print ("Done.\n\n");

exit 0;

sub parse_commandline {
	my @argvs = @_;
	my $argc = @argvs;
	my %options;
	my @valid_selects = ('--tree', '--tree_groups', '--nw_groups', '--out', '--size', '--interval', '--neighbor', '--reps');
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
	if ((exists ($options{'--size'}    ) && ($options{'--size'}     !~ /^\d+$/        )) ||
	    (exists ($options{'--reps'}    ) && ($options{'--reps'}     !~ /^\d+$/        )) ||
	    (exists ($options{'--interval'}) && ($options{'--interval'} !~ /^(\d*\.)?\d+$/)) ||
	    (exists ($options{'--neighbor'}) && ($options{'--neighbor'} !~ /^(\d*\.)?\d+$/)) ) {
		&usage;
	}
	if (!exists ($options{'--size'})) {
		$options{'--size'} = 10;
	}
	if (!exists ($options{'--interval'})) {
		$options{'--interval'} = 1;
	}
	if (!exists ($options{'--neighbor'})) {
		$options{'--neighbor'} = $options{'--interval'};
	}
	if (!exists ($options{'--reps'})) {
		$options{'--reps'} = 1;
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
	print ("Usage: tree_cluster.pl [--tree | --tree_groups | --nw_groups | --out] <file> [--size | --reps] <int> [--interval | --neighbor] <float>\n");
	print ("Mandatory arguments:\n");
	print (" --tree        input tree file (Newick / New Hampshire format)\n");
	print (" --out         basename of output files\n");
	print ("Optional arguments:\n");
	print (" --tree_groups groups and members file for tree\n");
	print (" --nw_groups   groups and members file for network\n");
	print (" --size        number of members in a cluster (default=10)\n");
	print (" --interval    minimum interval length between clusters (default=1)\n");
	print (" --neighbor    range of neighbor (default= same as --interval)\n");
	print (" --reps        maximum number of clusters in a neighbor (default=1)\n");
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

sub count_members {
	my ($root_node, $groups, $node_info) = @_;
	foreach my $node ($root_node->get_Descendents, $root_node) {
		if ($node->is_Leaf) {
			$node_info->{$node}{'size'} = (exists ($groups->{$node->id})) ? @{$groups->{$node->id}} : 1;
		}
	}
	return;
}

sub make_clusters {
	my ($min_size, $parent_node, $clusters, $node_info) = @_;
	if ($parent_node->is_Leaf) {
		if ($node_info->{$parent_node}{'size'} >= $min_size) {
			my %cluster = ('root_node' => $parent_node, 'valid' => 1);
			push (@$clusters, \%cluster);
			return 0;
		}
		else {
			return $node_info->{$parent_node}{'size'};
		}
	}
	else {
		my $have_cluster = 0;
		my $number_of_desc = 0;
		foreach my $child_node ($parent_node->each_Descendent) {
			my $size = &make_clusters ($min_size, $child_node, $clusters, $node_info);
			if ($size) {
				$number_of_desc += $size;
			}
			else {
				$have_cluster ++;
			}
		}
		if ($have_cluster) {
			$node_info->{$parent_node}{'have_clusters'} = 1;
			return 0;
		}
		elsif ($number_of_desc >= $min_size) {
			my %cluster = ('root_node' => $parent_node, 'valid' => 1);
			push (@$clusters, \%cluster);
			return 0;
		}
		else {
			return $number_of_desc;
		}
	}
}

sub interval_to_desc_cluster {
	my ($parent_node, $clusters, $node_info) = @_;
	if (grep {$_ eq $parent_node} map {$_->{'root_node'}} @$clusters) {
		$node_info->{$parent_node}{'iso'} = 0;
	}
	elsif ($parent_node->is_Leaf) {
		$node_info->{$parent_node}{'iso'} = undef;
	}
	else {
		my @intervals = ();
		foreach my $child_node ($parent_node->each_Descendent) {
			my $interval = &interval_to_desc_cluster ($child_node, $clusters, $node_info);
			if (defined ($interval)) {
				push (@intervals, $child_node->branch_length() + $interval);
			}
		}
		if (@intervals) {
			my $lowest_interval = $intervals[0];
			foreach my $interval (@intervals) {
				if ($lowest_interval > $interval) {
					$lowest_interval = $interval;
				}
			}
			$node_info->{$parent_node}{'iso'} = $lowest_interval;
		}
		else {
			$node_info->{$parent_node}{'iso'} = undef;
		}
	}
	return $node_info->{$parent_node}{'iso'};
}

sub set_interval_to_cluster {
	my ($upper_interval, $parent_node, $node_info) = @_;
	if (!defined ($node_info->{$parent_node}{'iso'}) || $upper_interval < $node_info->{$parent_node}{'iso'}) {
		$node_info->{$parent_node}{'iso'} = $upper_interval;
	}
	foreach my $child_node ($parent_node->each_Descendent) {
		&set_interval_to_cluster ($node_info->{$parent_node}{'iso'} + $child_node->branch_length, $child_node, $node_info);
	}
	return;
}

sub make_isolated_clusters {
	my ($min_interval, $parent_node, $clusters, $branches, $node_info) = @_;
	if (!grep {$_ eq $parent_node} map {$_->{'root_node'}} @$clusters) {
		if (($node_info->{$parent_node}{'iso'} >= $min_interval) && !($node_info->{$parent_node}{'have_clusters'})) {
			my %branch = ('root_node' => $parent_node, 'valid' => 1);
			push (@$branches, \%branch);
		}
		else {
			foreach my $child_node ($parent_node->each_Descendent) {
				&make_isolated_clusters ($min_interval, $child_node, $clusters, $branches, $node_info);
			}
		}
	}
	return;
}

sub branch_distance {
	my @nodes = @_[0, 1];
	my $paths = [];
	foreach my $num (0, 1){
		my $node = $nodes[$num];
		while ($node) {
			unshift (@{$paths->[$num]}, $node);
			$node = $node->ancestor;
		}
	}
	while (@{$paths->[0]} && @{$paths->[1]} && ($paths->[0][0] == $paths->[1][0])) {
		shift (@{$paths->[0]});
		shift (@{$paths->[1]});
	}
	my $distance = 0;
	foreach my $num (0, 1){
		foreach my $node (@{$paths->[$num]}) {
			$distance += $node->branch_length();
		}
	}
	return $distance;
}

# Reduce by (1) size of the representative (2) size of the cluster
sub reduce_clusters {
	my ($interval, $density, $clusters, $node_info) = @_;
	my $cluster_info = [];
	my $number_of_clusters = @$clusters;
	my $distances = [];
	foreach my $cluster_count1 (0 .. $number_of_clusters - 1) {
		$cluster_info->[$cluster_count1]{'insides'} = 0;
		foreach my $cluster_count2 (0 .. $number_of_clusters - 1) {
			if ($cluster_count1 < $cluster_count2) {
				$distances->[$cluster_count1][$cluster_count2] = &branch_distance ($clusters->[$cluster_count1]{'root_node'}, $clusters->[$cluster_count2]{'root_node'});
			}
			elsif ($cluster_count1 == $cluster_count2) {
				$distances->[$cluster_count1][$cluster_count2] = 0;
			}
			else {
				$distances->[$cluster_count1][$cluster_count2] = $distances->[$cluster_count2][$cluster_count1];
			}
			if ($distances->[$cluster_count1][$cluster_count2] < $interval) {
				$cluster_info->[$cluster_count1]{'insides'} ++;
			}
		}
		$cluster_info->[$cluster_count1]{'size'} = 0;
		foreach my $member_node ($clusters->[$cluster_count1]{'root_node'}->get_Descendents, $clusters->[$cluster_count1]{'root_node'}) {
			if ($member_node->is_Leaf) {
				$cluster_info->[$cluster_count1]{'size'} += $node_info->{$member_node}{'size'};
			}
		}
	}
	my $number_of_del = 0;
	while (grep {$_->{'insides'} > $density} @$cluster_info) {
		my $min_representative_size = undef;
		my $min_size = undef;
		my @del_cluster_counts;
		foreach my $cluster_count1 (0 .. $number_of_clusters - 1) {
			if ($cluster_info->[$cluster_count1]{'insides'} > $density) {
				if (!defined ($min_representative_size) || ($node_info->{$clusters->[$cluster_count1]{'representative'}}{'size'} < $min_representative_size)) {
					$min_representative_size = $node_info->{$clusters->[$cluster_count1]{'representative'}}{'size'};
					$min_size = undef;
					@del_cluster_counts = ();
				}
				if ($node_info->{$clusters->[$cluster_count1]{'representative'}}{'size'} == $min_representative_size) {
					if (!defined ($min_size) || ($cluster_info->[$cluster_count1]{'size'} < $min_size)) {
						$min_size = $cluster_info->[$cluster_count1]{'size'};
						@del_cluster_counts = ($cluster_count1);
					}
					elsif ($cluster_info->[$cluster_count1]{'size'} == $min_size) {
						push (@del_cluster_counts, $cluster_count1);
					}
				}
			}
		}
		foreach my $cluster_count (@del_cluster_counts) {
			$clusters->[$cluster_count]{'valid'} = 0;
		}
		$number_of_del += @del_cluster_counts;
		foreach my $cluster_count1 (0 .. $number_of_clusters - 1) {
			$cluster_info->[$cluster_count1]{'insides'} = 0;
			if ($clusters->[$cluster_count1]{'valid'}) {
				foreach my $cluster_count2 (0 .. $number_of_clusters - 1) {
					if ($clusters->[$cluster_count2]{'valid'} && $distances->[$cluster_count1][$cluster_count2] < $interval) {
						$cluster_info->[$cluster_count1]{'insides'} ++;
					}
				}
			}
		}
	}
	return $number_of_del;
}

sub choose_representative {
	my ($root_node, $node_info) = @_;
	my @nodes = grep {$_->is_Leaf} ($root_node->get_Descendents, $root_node);
	my $distances = [];
	my @total_distances = (0) x @nodes;
	foreach my $node_count1 (0 .. $#nodes) {
		foreach my $node_count2 (0 .. $#nodes) {
			if ($node_count1 < $node_count2) {
				$distances->[$node_count1][$node_count2] = &branch_distance ($nodes[$node_count1], $nodes[$node_count2]);
			}
			elsif ($node_count1 == $node_count2) {
				$distances->[$node_count1][$node_count2] = 0;
			}
			else {
				$distances->[$node_count1][$node_count2] = $distances->[$node_count2][$node_count1];
			}
			$total_distances[$node_count1] += $distances->[$node_count1][$node_count2] * $node_info->{$nodes[$node_count2]}{'size'};
		}
	}
	my $central_node_count = 0;
	foreach my $node_count (1 .. $#nodes) {
		if ($total_distances[$node_count] < $total_distances[$central_node_count]) {
			$central_node_count = $node_count;
		}
		elsif ($total_distances[$node_count] == $total_distances[$central_node_count] &&
		       $node_info->{$nodes[$node_count]}{'size'} > $node_info->{$nodes[$central_node_count]}{'size'}) {
			$central_node_count = $node_count;
		}

	}
	return $nodes[$central_node_count];
}
