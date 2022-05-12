#!/usr/bin/perl -w

use strict;
use warnings;

my %options = &parse_commandline (@ARGV);
my @classes = &parse_classes ($options{'-C'});
if (@classes == 0) {
	print STDERR ("-C is invalid.\n");
	exit 1;
}

my $s_filename = $options{'-S'};
open (my $fh_s_in, '<', $s_filename) or die "Cannot open $s_filename: $!";
my $groups_table = &make_groups_for_mega ($fh_s_in, \@classes);
close ($fh_s_in);

if (defined ($options{'-G'})) {
	my $group_filename = $options{'-G'};
	open (my $fh_s_in, '<', $group_filename) or die "Cannot open $group_filename: $!";
	&group_name ($fh_s_in, $groups_table);
	close ($fh_s_in);
}

&output_groups ($groups_table);


exit;

sub parse_commandline {
	my @argvs = @_;
	my $argc = @argvs;
	my %options;
	my @valid_selects = ('-S', '-G', '-C');
	for (my $i = 0; $i < $argc; $i ++) {
		if (grep (/^($argvs[$i])$/, @valid_selects) && ($i + 1 < $argc)) {
			$options{$argvs[$i]} = $argvs[$i + 1];
			$i ++;
		}
		else {
			&usage;
		}
	}
	if (!exists ($options{'-S'})) {&usage;}
	if (!exists ($options{'-C'})) {$options{'-C'} = '0,5000';}
	print STDERR ("Arguments:\n");
	foreach my $argv_names (@valid_selects) {
		if (exists ($options{$argv_names})) {
			print STDERR ($argv_names, ' ', $options{$argv_names}, "\n");
		}
	}
	print STDERR ("\n");
	return %options;
}

sub usage {
	print STDERR ("Usage: sstar2mega_group.pl [-S | -G | -C] <args>\n");
	print STDERR ("Mandatory arguments:\n");
	print STDERR (" -S  S* output file\n");
	print STDERR ("Optional arguments:\n");
	print STDERR (" -G  Individual-group file (tab delimited)\n");
	print STDERR (" -C  comma separated classes of S* values (>0)\n");
	print STDERR ("EXAMPLE\n");
	print STDERR ("sstar2mega_group.pl -S a_s.txt -G a_group -C 5000,20000,50000,100000 \n");
	exit 1;
}

sub parse_classes {
	my $classes_str = $_[0];
	my @classes = split (',', $classes_str);
	if (grep (/\D/, @classes)) {@classes = ();}
	else {@classes = sort {$a <=> $b} grep (/\d/, @classes);}
	return @classes;
}

sub make_groups_for_mega {
	my ($fh_in, $classes) = @_;
	my $groups_table = [];
	my $line;
	while (defined ($line = <$fh_in>) && $line !~ /^ID\tS\*\tset of sites\W/) {}
	if (defined ($line)) {
		while ($line = <$fh_in>) {
			chomp ($line);
			my ($name, $s_value, undef) = split ("\t", $line, 3);
			if ($name ne '.') {
				my $class;
				if    ($s_value eq '///')    {$class = 'NA';}
				elsif ($s_value =~ /^-inf/i) {$class = '-inf';}
				elsif ($s_value eq '-10000') {$class = -10000;}
				elsif ($s_value =~ /\D/)     {$class = 'invalid';}
				elsif ($s_value == 0)        {$class = 0;}
				else {
					$class = 0;
					foreach my $th (@$classes) {
						if ($s_value < $th) {last;}
						else {$class = $th;}
					}
				}
				push (@$groups_table, {'name' => $name, 'class' => $class});
			}
		}
	}
	else {print ("No data\n");}
	return $groups_table;
}

sub group_name {
	my ($fh_in, $groups_table) = @_;
	my %groups;
	while (my $line = <$fh_in>) {
		chomp ($line);
		my ($indv, $group, undef) = split ("\t", $line, 3);
		if (defined ($indv) && $indv ne '') {$groups{$indv} = $group;}
	}
	foreach my $indv (@$groups_table) {
		if (defined ($groups{$indv->{'name'}})) {$indv->{'name'} = $groups{$indv->{'name'}};}
	}
	return;
}

sub output_groups {
	my $groups_table = $_[0];
	foreach my $indv (@$groups_table) {
		print ($indv->{'name'}, '=', $indv->{'class'}, "\n");
	}
	return;
}