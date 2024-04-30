#!/usr/bin/perl -w

# Script to call FragGeneScan in parallel, submitting the parts of the task to multiple
# cluster nodes.

use warnings;
use strict;
use File::Basename;

our $MAX_PIECES = 100; # control splitting files maxium to 100 pieces
my $t1= time;

# if ( !defined( $ENV{PHASTEST_CLUSTER_HOME} ) ) {
#     die 'The PHASTEST_CLUSTER_HOME enviroment variable is not set\n';
# }
my $scripts_dir = "$ENV{PHASTEST_HOME}/scripts";

my $fna_file = $ARGV[0];
# my $log_file = "$fna_file.log";
my $dir = dirname($fna_file);
my $fna_file_basename = basename($fna_file);

my $count = 0;

my ($rs) = `grep '>' $fna_file | wc -l` =~ m/^(\d+)/;
if ($rs > $MAX_PIECES){
	$count = $MAX_PIECES;
}else{
	$count = $rs;
}
my $mod_rest;
my $piece;
if ($count == 0){
   $mod_rest = 0;
   $piece = 1;
}else{
   $mod_rest = $rs % $count;
   $piece = int($rs/$count); # at least needed pieces in each file
}
# split fna file with multiple fasta records 
my $part = '';
my $need;
my $cou = 0;
my $piec = 0; # file postfix from 1-$count
my $start = 1;
open(IN, $fna_file) or die "Cannot open $fna_file";
while(my $l = <IN>){
	
	if($l =~ m/^>/){
		$part .= "\n" if ($part ne '');
		$cou++;
		if ($mod_rest >0){
			$need = $piece + 1;
		}else{
			$need = $piece;
		}
		if ($cou > $need){
			$mod_rest = $mod_rest - 1 if ($mod_rest >0);
			$piec++;
			$start = write_file($part, $piec, $fna_file, $start);
			$part = '';
			$cou = 0;
			
		}
		$part .= $l;
		$cou = 1 if ($cou==0);
	}	
	else{
		chomp($l);
		$part .= $l;			
	}
}
#write the last record
if ($part ne ''){
	$part .= "\n" if ($part ne '');
	$piec++;
	$start = write_file($part, $piec, $fna_file, $start);
	$part = '';
	$cou = 0;
	
}
close IN;

print "fraggenescan: piece $count\n";
# system("echo 'piece $count' >> $log_file");

chdir $dir;
print "change to $dir\n";
# system("echo change to $dir >> $log_file");

# Get list of all cluster nodes.
# Note: We only want to run on nodes that map_node.pl is aware of, just in case a new
# node is added to the cluster that does not have the needed database files on it.
my $node_list = `perl $scripts_dir/map_node.pl 0`;
my @lines = split("\n", $node_list);
my $node_set = join '|', @lines;
my $sub_programs_dir = "$ENV{PHASTEST_HOME}/sub_programs";
my $fraggenescan_exec = "$sub_programs_dir/FragGeneScan1.20/run_FragGeneScan.pl";
print "running $cm\n";
system("perl $fraggenescan_exec  -genome=$fna_file_basename\_fgs_out   -out=$fna_file_basename\_fgs_out -complete=0 -train=complete")==0 or print "died $!";

open(OUT, ">$fna_file.predict") or die "Cannot write $fna_file.predict";
print OUT ">gi|000000|ref|NC_000000|  Concatenated genome, 1-".($start-1)."\n";
my $orf_count =0;
for (my $i=1; $i<=$count; $i++){
	open(IN, "$fna_file_basename\_$i.out") or die "Cannot open $fna_file_basename\_$i.out";
	my $start_num;
	while (my $l = <IN>){
		if ($l=~/^>/){
			if($l=~ m/ (\d+)\-\d+$/){
				$start_num = $1;
			}else{
				die "Cannot get start number from $fna_file_basename\_$i.out";
			}
			next;	
		}else{
			my @arr=split(" ", $l);
			$arr[0] += $start_num-1;
			$arr[1] += $start_num-1;
			$orf_count++; 
			my $orf_str='';
			if ($orf_count >0 && $orf_count<10) {
				$orf_str='0000'.$orf_count;
			}elsif($orf_count >=10 && $orf_count<100) {
				$orf_str='000'.$orf_count;
			}elsif($orf_count >=100 && $orf_count<1000) {
				$orf_str='00'.$orf_count;
			}elsif($orf_count >=1000 && $orf_count<10000) {
				$orf_str='0'.$orf_count;
			}else{
				$orf_str=$orf_count;
			}
			my $line = '';
			($arr[0], $arr[1], $arr[2]) = check_length($arr[0], $arr[1], $arr[2]);
			if ($arr[2]=~/\+/) {
				$line = sprintf("%-10s %-10s %-10s %-5s %-20s\n", "orf$orf_str", $arr[0], $arr[1], "$arr[2]$arr[3]", $arr[4]);
			}elsif ($arr[2]=~/\-/){
				$line = sprintf("%-10s %-10s %-10s %-5s %-20s\n", "orf$orf_str", $arr[1], $arr[0], "$arr[2]$arr[3]", $arr[4]);
			}else{
				die "Format in $fna_file_basename\_$i.out not correct on line $l\n";
			}
			print OUT $line;	
		}
				
	}
	close IN;
}
close OUT;
print "run time= ". (time - $t1)." seconds\n";
print "Program exit\n";


sub write_file {
	my ($part, $piec, $fna_file,$start) = @_;
	my @tmp = split("\n", $part);
	
	open(OUT, "> $fna_file\_$piec") or die "Cannot write $fna_file\_$piec";
	for(my $i=0; $i<=$#tmp; $i++){
		if ($tmp[$i]=~m/^>/){
			$tmp[$i+1]=~s/[\n\s]//g;
			my $len = length($tmp[$i+1]);
			print OUT "$tmp[$i] $start\-". ($start - 1 + $len)."\n";
			my @arr = $tmp[$i+1]=~m/(\w{0,60})/g;
			foreach my $l (@arr){
				next if($i=~m/^\s*$/);
				print OUT $l."\n";
			}
			$start = $start + $len;
		}
	} 
	return $start;
}
sub check_length{
	my ($start, $end, $forward_or_backward) = @_;
	my $rest = ($end - $start + 1) % 3;
	if ($forward_or_backward =~ m/\+/){
		if ($rest != 0){
			$start += $rest;
		}
	}elsif($forward_or_backward =~ m/\-/){
		if ($rest != 0){
			$end -= $rest;
		}
	}
	else{
		die "forward_or_backward not clear in line with $start and $end";
	}
	return $start, $end, $forward_or_backward;
}
	
