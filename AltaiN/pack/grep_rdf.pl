#!/usr/bin/perl -w

use strict;
use warnings;

my $len_seq = 6; # number of letters of a sequence name in RDF file
my $len_char = 6; # number of letters of a character name in RDF file

if (@ARGV < 3 || 4 < @ARGV) {
	&arg_error;
}

my $edit_mode = shift (@ARGV);
if ($edit_mode ne 'in' && $edit_mode ne 'ex') {
	&arg_error;
}

my $id_mode = shift (@ARGV);
my @ids;
if ($id_mode eq '--id') {
	@ids = split (',', shift (@ARGV));
}
elsif ($id_mode eq '--idfile') {
	@ids = &ids_from_file (shift (@ARGV));
}
else {
	&arg_error;
}

my $in_rdf_filename = '';
if (@ARGV > 0) {
	$in_rdf_filename = shift (@ARGV);
}

my $fh_in_rdf;
if ($in_rdf_filename) {
	open ($fh_in_rdf, '<', $in_rdf_filename) or die "Cannot open $in_rdf_filename: $!";
}
else {
	$fh_in_rdf = *STDIN;
}

my $outputmode = 1;
for (my $line_count = 0; $line_count < $len_char; $line_count ++) {
	my $line = <$fh_in_rdf>;
	print ($line);
}
while (my $line = <$fh_in_rdf>) {
	if ($line =~ /^\r?\n$/) {
		print ($line);
		last;
	}
	(my $seq_name, undef) = split (/ +/, $line, 2);
	if ((grep {$_ eq $seq_name} @ids) xor ($edit_mode eq 'ex')) {
		print ($line);
	}
}
while (my $line = <$fh_in_rdf>) {
	print ($line);
}

if ($in_rdf_filename) {
	close ($fh_in_rdf);
}
exit (0);


sub ids_from_file
{
	my $id_filename = $_[0];
	open (my $fh_in, '<', $id_filename) or die "Cannot open $id_filename: $!";
	my @ids;
	while (my $line = <$fh_in>) {
		chomp ($line);
		push (@ids, $line);
	}
	return @ids;
}

sub arg_error
{
	print STDERR ("Usage1: grep_rfd.pl in --id <ids> [input_RDF_file]\n");
	print STDERR ("Usage2: grep_rfd.pl in --idfile <id_filename> [input_RDF_file]\n");
	print STDERR ("Usage3: grep_rfd.pl ex --id <ids> [input_RDF_file]\n");
	print STDERR ("Usage4: grep_rfd.pl ex --idfile <id_filename> [input_RDF_file]\n");
	exit (1);
}


__END__

概要
RDFファイルから、指定されたIDを持つOTUだけを抽出/除外したRDFファイルを出力する。

使い方
grep_rfd.pl in --id <ids> [input_RDF_file]
grep_rfd.pl in --idfile <id_filename> [input_RDF_file]
grep_rfd.pl ex --id <ids> [input_RDF_file]
grep_rfd.pl ex --idfile <id_filename> [input_RDF_file]

第1引数に“in”を指定すると抽出、“ex”を指定すると除外する。
IDの指定をコマンドライン引数で行うときは、第2引数に“--id”を指定し、第3引数に,区切りでIDを指定する。
IDの指定をファイルで行うときは、第2引数に“--idfile”を指定し、第3引数にIDファイル名を指定する。IDファイルには1行ごとに1つのIDを記述する。
第4引数に入力RDFファイルを指定する。省略されたときは標準入力から読み込む。

指定されたIDが入力RDFファイルに存在しない場合でも黙って処理を続ける。

