#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV < 2 || 4 < @ARGV) {
	print STDERR ("Usage: rm_zerodist.pl in_distance_matrix_file in_phylip_file [out_phylip_file [out_log_file]]\n");
	exit (1);
}
my $matrix_filename = $ARGV[0];
my $in_phy_filename = $ARGV[1];
my $out_phy_filename;
my $log_filename;
if (@ARGV == 2) {
	if ($in_phy_filename =~ /^(.*)\.phy$/) {
		$out_phy_filename = $1 . '.uniq.phy';
	}
	else {
		$out_phy_filename = $in_phy_filename . '.uniq.phy';
	}
}
else {
	$out_phy_filename = $ARGV[2];
}
if (@ARGV <= 3) {
	if ($out_phy_filename =~ /^(.*)\.uniq\.phy$/) {
		$log_filename = $1 . '.zerodist';
	}
	elsif ($out_phy_filename =~ /^(.*)\.phy$/) {
		$log_filename = $1 . '.zerodist';
	}
	else {
		$log_filename = $out_phy_filename . '.zerodist';
	}
}
else {
	$log_filename = $ARGV[3];
}

# check filenames
my @op_check = &search_same_values ($matrix_filename, $in_phy_filename, $out_phy_filename, $log_filename);
if ($op_check[1] > 0) {
	my @err_filename = ('in_distance_matrix_file', 'in_phylip_file', 'out_phylip_file', 'out_log_file');
	print STDERR ($err_filename[$op_check[0]], ' and ', $err_filename[$op_check[1]], " cannot be the same file.\n");
	exit (1);
}
print ("Input distance matrix file: $matrix_filename\n");
print ("Input alignment file: $in_phy_filename\n");
print ("Output alignment file: $out_phy_filename\n");
print ("Output report file: $log_filename\n\n");

my $line;
my @columns;

# load distance matrix file
open (my $fh_in, '<', $matrix_filename) or die "Cannot open $matrix_filename: $!";
$line = <$fh_in>;
my $num_sp;
if ($line =~ /^\s*(\d+)/) {
	$num_sp = $1;
}
else {
	print STDERR ("Cannot read the number of species in $matrix_filename.\n");
	close ($fh_in);
	exit (2);
}
my $dist_mat;
my @sp_names;
for (my $sp1 = 0; $sp1 < $num_sp; $sp1 ++) {
	$dist_mat->[$sp1] = [];
	$sp_names[$sp1] = '';
	while (@{$dist_mat->[$sp1]} < $num_sp) {
		$line = <$fh_in>;
		if (!$line) {
			print STDERR ("The size of the matrix is smaller than the number of species ($num_sp).\n");
			close ($fh_in);
			exit (2);
		}
		chomp ($line);
		@columns = split (' ', $line);
		if ($sp_names[$sp1] eq '') {
			$sp_names[$sp1] = shift (@columns);
		}
		push (@{$dist_mat->[$sp1]}, @columns);
	}
}
close ($fh_in);

# check distance matrix
for (my $sp1 = 0; $sp1 < $num_sp; $sp1 ++) {
	if ($dist_mat->[$sp1][$sp1] != 0) {
		print STDERR ('The distance ', $sp_names[$sp1], ' <-> ', $sp_names[$sp1], ' is ', $dist_mat->[$sp1][$sp1], ".\nIt must be 0.\n");
		exit (2);
	}
	for (my $sp2 = 0; $sp2 < $sp1; $sp2 ++) {
		if ($dist_mat->[$sp1][$sp2] != $dist_mat->[$sp2][$sp1]) {
			print STDERR ('The distances ', $sp_names[$sp1], ' <-> ', $sp_names[$sp2], ' and ', $sp_names[$sp2], ' <-> ', $sp_names[$sp1], " differ.\n",
			              'The former is ', $dist_mat->[$sp1][$sp2], ', while the latter is ', $dist_mat->[$sp2][$sp1], ".\n They must be the same.\n");
			exit (2);
		}
	}
}

# search zero-distance and check them
my $same_sp = [];
my $same_sp1;
my @is_same = (-1) x $num_sp;
for (my $sp1 = 0; $sp1 < $num_sp; $sp1 ++) {
	if ($is_same[$sp1] == -1) {
		$same_sp1 = [$sp1];
		for (my $sp2 = $sp1 + 1; $sp2 < $num_sp; $sp2 ++) {
			if ($dist_mat->[$sp1][$sp2] == 0) {
				print ($sp_names[$sp2], ' is the same as ', $sp_names[$sp1], ".\n");
				if ($is_same[$sp2] == -1) {
					push (@$same_sp1, $sp2);
					$is_same[$sp2] = $sp1;
				}
				else {
					&dist_mat_zero_error ($sp2, $is_same[$sp2], $sp1);
				}
			}
		}
		if (@{$same_sp1} > 1) {
			push (@$same_sp, $same_sp1);
		}
	}
	else {
		for (my $sp2 = $sp1 + 1; $sp2 < $num_sp; $sp2 ++) {
			if ($dist_mat->[$sp1][$sp2] == 0) {
				print ($sp_names[$sp2], ' is the same as ', $sp_names[$sp1], ".\n");
				if ($is_same[$sp1] != $is_same[$sp2]) {
					&dist_mat_zero_error ($sp1, $is_same[$sp1], $sp2);
				}
			}
			else {
				if ($is_same[$sp1] == $is_same[$sp2]) {
					&dist_mat_zero_error ($is_same[$sp1], $sp1, $sp2);
				}
			}
		}
	}
}

# set information of duplicate species
my @group_names;
my $member_names;
my $num_group = @$same_sp;
my $num_group_sp;
my %del_names;
my $uniq_num_sp = $num_sp;
for (my $group = 0; $group < $num_group; $group ++) {
	$num_group_sp = @{$same_sp->[$group]};
	$uniq_num_sp -= $num_group_sp - 1;
	$group_names[$group] = sprintf ('H%d-%d', $group + 1, $num_group_sp);
	for (my $sp = 0; $sp < $num_group_sp; $sp ++) {
		$member_names->[$group][$sp] = $sp_names[$same_sp->[$group][$sp]];
		$del_names{$member_names->[$group][$sp]} = 0;
	}
	$group_names[$group] = &check_special_member ($member_names->[$group]) . $group_names[$group];
	if (length ($group_names[$group]) > 10) {
		$group_names[$group] = substr ($group_names[$group], 0, 10);
	}
	$del_names{$member_names->[$group][0]} = $group_names[$group];
}

# output zero distance list
open (my $fh_out, '>', $log_filename) or die "Cannot open $log_filename: $!";
print $fh_out ("Group name: members\n");
for (my $group = 0; $group < $num_group; $group ++) {
	print $fh_out ($group_names[$group], ': ', join (',',@{$member_names->[$group]}), "\n");
}
close ($fh_out);

# read input phylip file and write output phylip file
open ($fh_in, '<', $in_phy_filename) or die "Cannot open $in_phy_filename: $!";
$line = <$fh_in>;
my $num_sp_phy;
if ($line =~ /^\s*(\d+)/) {
	$num_sp_phy = $1;
}
else {
	print STDERR ("Cannot read the number of species in $in_phy_filename.\n");
	close ($fh_in);
	exit (3);
}
if ($num_sp != $num_sp_phy) {
	print STDERR ("The number of species in $matrix_filename is $num_sp, but that of $in_phy_filename is $num_sp_phy.\n",
	              "They must be the same.\n");
	close ($fh_in);
	exit (4);
}

open ($fh_out, '>', $out_phy_filename) or die "Cannot open $out_phy_filename: $!";
$line =~ s/$num_sp/$uniq_num_sp/;
print $fh_out ($line);

my @is_uniq;
my $name;
my $seq;
for (my $num_line = 0; $num_line < $num_sp; $num_line ++) {
	$line = <$fh_in>;
	if (!$line) {
		print STDERR ("Cannot read sequences in $in_phy_filename.\n");
		close ($fh_in);
		exit (3);
	}
	($name, $seq) = split (/\s+/, $line, 2);
	if (exists ($del_names{$name})) {
		if ($del_names{$name} eq '0') {
			$is_uniq[$num_line] = 0;
			next;
		}
		else {
			$name = $del_names{$name};
		}
	}
	$is_uniq[$num_line] = 1;
	printf $fh_out ('%-10s %s', $name, $seq);
}
my $num_line = 0;
while ($line = <$fh_in>) {
	if ($line eq "\n") {
		print $fh_out ($line);
		next;
	}
	if ($is_uniq[$num_line]) {
		print $fh_out ($line);
	}
	$num_line ++;
	if ($num_line >= $num_sp) {
		$num_line = 0;
	}
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

sub dist_mat_zero_error
{
	my @err_sp_names;
	for (my $i = 0; $i < 3; $i ++) {
		$err_sp_names[$i] = $sp_names[$_[$i]];
	}
	my $nonzero = $dist_mat->[$_[1]][$_[2]];
	print STDERR ('Both the distances ', $err_sp_names[0], ' <-> ', $err_sp_names[1], ' and ', $err_sp_names[0], ' <-> ', $err_sp_names[2], " are 0,\n",
	              'but the distance ', $err_sp_names[1], ' <-> ', $err_sp_names[2], ' is ', $nonzero, ", not 0.\n");
	exit (2);
}

sub check_special_member
{
	my $members = $_[0];
	my @special_names  = ('Href_hg19', 'DenisovaMj');
	my @special_prefix = ('Hr'       , 'D'         );
	my $num_special = @special_names;
	my $prefix = '';
	for (my $num = 0; $num < $num_special; $num ++) {
		foreach my $member (@$members) {
			if ($special_names[$num] eq $member) {
				$prefix .= $special_prefix[$num];
				last;
			}
		}
	}
	return $prefix;
}

__END__

return code
1: Argument error.
2: Logical error in the distance matrix file.
3: Logical error in the input phylip file.
4: Logical error between the distance matrix file and the input phylip file.

