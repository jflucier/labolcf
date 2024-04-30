#!/usr/bin/perl -w

use warnings;
use strict;

my $filename = $ARGV[0];
open(IN, "$filename");
open(OUT, ">$filename.tmp");
my $head = '';
my $seq_tmp = '';
while(<IN>){
    if ($_=~/>/){
        if ($seq_tmp ne '') {
            my @tmp = $seq_tmp =~/(\w{0,60})/g;
            print OUT $head;
            print OUT join("\n", @tmp);
            $seq_tmp = '';
        }

        $head = $_;
    }else{
        if ($_ =~ /([A-Z]+)/i) {
            $seq_tmp .= $1;
        }
    }
}
close IN;
if ($seq_tmp ne '') {
    my @tmp = $seq_tmp =~/(\w{0,60})/g;
            $head =~s/>/-/g; $head=~s/^-/>/;
    print OUT $head;
    print OUT join("\n", @tmp);
    $seq_tmp = '';
}
close OUT;
system("mv $filename.tmp $filename");
