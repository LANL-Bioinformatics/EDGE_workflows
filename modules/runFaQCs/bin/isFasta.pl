#!/usr/bin/env perl

use strict;
use warnings;

$SIG{'PIPE'}=sub{};
my $file=shift @ARGV;
my ($fh,$pid)= open_file($file);
my $head=<$fh>;
close $fh;
kill 9, $pid; # avoid gunzip broken pipe
$SIG{'PIPE'} = 'DEFAULT';
($head =~/^>/)?
    print "Yes":
    print "No";


sub open_file
{
    my ($file) = @_;
    my $fh;
    my $pid;
    if ( $file=~/\.gz\.?\d?$/i ) { $pid=open($fh, "gunzip -c $file |") or die ("gunzip -c $file: $!"); }
    else { $pid=open($fh,'<',$file) or die("$file: $!"); }
    return ($fh,$pid);
}

