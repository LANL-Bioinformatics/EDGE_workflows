#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $id_file;
my $table;

GetOptions(
    't=s' => \$table,
    'i=s' => \$id_file
);

my %id_map;
print "$id_file";
print "$table";
open(my $fh,$id_file) or die "Cannot read id_map.txt\n";
while(<$fh>){chomp; my($new_id,$original_id)=split; $id_map{$new_id}=$original_id;}
close $fh;
open(my $ofh,">phageFinder_summary.txt") or die "Cannot write phageFinder_summary.txt\n";
open(my $result_fh,$table) or die "Cannot read PFPR_tab.txt\n";
while(<$result_fh>)
{
    my @fields=split /\s+/,$_;
    $fields[0]=$id_map{$fields[0]} if ($id_map{$fields[0]});
    print $ofh join("\t",@fields),"\n";
}
close $result_fh;       
close $ofh;