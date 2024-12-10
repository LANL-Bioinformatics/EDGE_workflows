#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $paired_fasta='';
my $single_fasta='';
my $outputDir='';
my $len_cutoff=350;

GetOptions(
    'p:s' => \$paired_fasta,
    'u:s' => \$single_fasta,
    'd=s' => \$outputDir,
    'len:s' => \$len_cutoff
);


my $short_paired_fasta="short_paired.fa";
my $short_single_fasta="short_single.fa";
my $long_fasta="long.fa";

open (my $o_paired, ">$short_paired_fasta") or die "Cannot write $short_paired_fasta\n";
open (my $o_single, ">$short_single_fasta") or die "Cannot write $short_single_fasta\n";
open (my $o_long, ">$long_fasta") or die "Cannot write $long_fasta\n";
$/ = ">";
if (-s $paired_fasta) {
    open (my $fh, $paired_fasta) or die "$! $paired_fasta";
    while (<$fh>)
    { 
        $_ =~ s/\>//g;
        my ($id, @seq) = split /\n/, $_;
        next if (!$id);
        my ($id2, @seq2) = split /\n/, <$fh>;
        my $seq = join "", @seq;
        my $seq2 = join "", @seq2;
        my $len = length($seq);
        my $len2 = length($seq2);
        if ($len > $len_cutoff and $len2 > $len_cutoff)
        {
                print $o_long ">$id\n$seq\n>$id2\n$seq2\n";
        }
        elsif ($len > $len_cutoff)
        {
                print $o_long ">$id\n$seq\n";
                print $o_single ">$id2\n$seq2\n";
        }
        elsif ($len2 > $len_cutoff)
        {
                print $o_long ">$id2\n$seq2\n";
                print $o_single ">$id\n$seq\n";
        }
        else
        {
                print $o_paired ">$id\n$seq\n>$id2\n$seq2\n";   
        }
    }
    close $fh;
}
if (-s $single_fasta)
{
    open (my $fh, $single_fasta) or die "$! $single_fasta";
    while (<$fh>)
    { 
        $_ =~ s/\>//g;
        my ($id, @seq) = split /\n/, $_;
        next if (!$id);
        my $seq = join "", @seq;
        my $len = length($seq);
        if ($len > $len_cutoff)
        {
                print $o_long ">$id\n$seq\n";
        }
        else
        {
                print $o_single ">$id\n$seq\n";
        }
    }
    close $fh;
}
$/="\n";
close $o_paired;
close $o_long;
close $o_single;