#/usr/bin/perl -w

use strict;
use warnings;

use Cwd;

my $jobs_dir = $ARGV[0];
my $num = $ARGV[1];
make_region_DNA($jobs_dir,$num);

sub make_region_DNA{
	my($jobs_dir,$num)=@_;
	my $cur_dir = getcwd();
	chdir "$jobs_dir/$num";
	my $DNA_seq= `grep -v '>' $num.fna`;
	$DNA_seq=~s/[\s\n]//gs;
	my $sum_file_content=`cat summary.txt`;
	$sum_file_content=~s/.*---------------+\n//s;
	my @pos=();
	foreach my $line(split "\n", $sum_file_content){
		my @tmp=split " ", $line;
		push @pos, $tmp[4];
	}
	if (@pos !=0){
		open(OUT, ">region_DNA.txt") or die "Cannot write region_DNA.txt";
		for(my $i=0; $i <=$#pos; $i++){
			my ($start, $end)= $pos[$i]=~/(\d+)-(\d+)/;
			print OUT ">".($i+1)."\t $pos[$i]\n";
			my $seq= substr($DNA_seq, $start-1, $end-$start+1);
			my @tmp=$seq=~/(\w{0,60})/gs;
			foreach my $l(@tmp){
				print OUT $l."\n";
			}
			print OUT "\n";
		}
		close OUT;
	}
	chdir $cur_dir;
}