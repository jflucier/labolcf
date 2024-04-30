#!/usr/bin/perl -w

use strict;
use warnings;

my $fna = $ARGV[0];
my $job_id = $ARGV[1];

open(OUT,  "> $fna.fna") or die "Cannot generate $fna.fna. $!";
my $fna_header = find_common_element("$job_id\_contig_positions.txt");
print OUT ">gi|00000000|ref|NC_000000| $fna_header\n";
print OUT `grep -v '>' $job_id\_original.fna`;
close OUT;

# Subroutine to find the common element in the names of the contigs.
# For example, for the given file:
# 	SGQG01000020.1,Pectobacterium,carotovorum,subsp.,carotovorum,strain,Ec-006,NODE_20_length_10311_cov_7.781415,,whole,genome,shotgun,sequence	1	10311	10311
#	SGQG01000021.1,Pectobacterium,carotovorum,subsp.,carotovorum,strain,Ec-006,NODE_21_length_8273_cov_164.985237,,whole,genome,shotgun,sequence	10312	18584	8273
# We would return:
# Pectobacterium carotovorum subsp. carotovorum strain Ec-006 whole genome shotgun sequence
sub find_common_element {
	my ($contig_file) = @_;
	my $common_element = '';
	my %common_elem_hash = ();

	open (IN, $contig_file) || die "Could not open file $contig_file: $!\n";
	my $first_line = <IN>;
	chomp($first_line);
	my @fields = split(/\t/, $first_line);
	my $description = $fields[0];
	my @first_words = split(/,/, $description);

	# Iterate each line of the contig_positions file.
	while (my $line = <IN>) {
		chomp($line);
		my @fields = split(/\t/, $line);
		my $description = $fields[0];
		my @words = split(/,/, $description);

		# Iterate each word in the description.
		for (my $i = 0; $i < scalar(@first_words); $i++) {
			if ($first_words[$i] eq $words[$i]) {
				$common_elem_hash{$i} = 1;
			}
			else {
				$common_elem_hash{$i} = 0;
			}
		}
	}
	close IN;

	# Iterate the hash and build the name from common element.
	my $common_elem_count = 0;
	foreach my $key (sort {$a<=>$b} keys %common_elem_hash) {
		if ($common_elem_hash{$key} == 1) {

			# If the word is empty string, then it's likely a comma.
			if ($first_words[$key] eq '') {
				$common_element =~ s/\s$/, /;
				next;
			}

			$common_element .= $first_words[$key] . " ";
			$common_elem_count++;
		}
	}
	if ($common_elem_count < 5) {
		return $description;
	}
	else {
		return $common_element;
	}
}