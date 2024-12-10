#!/usr/bin/env perl

use FindBin qw($Bin);
use Getopt::Long;
use strict;
use warnings;
use LWP::UserAgent;

my $url;

GetOptions(
    "url=s" => \$url
);


my $ua = LWP::UserAgent->new;
$ua->timeout(10);
$ua->env_proxy;
my $up=0;
my $response = $ua->get($url);

if ($response->is_success) {
        $up=1
}
else {
    die("$url is not up!")
}