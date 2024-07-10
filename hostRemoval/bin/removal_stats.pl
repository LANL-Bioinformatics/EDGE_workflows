#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;


my @host_files = ();
my @total_reads;
my @host_reads;
my @host_names;
my $outDir = '.';
my $hostRemovalPDF = "$outDir/HostRemoval/HostRemovalStats.pdf";
my $hostclean_stat_file = "$outDir/HostRemoval/hostclean.stats.txt";
my $total_host=0;
my $clean_stats;


GetOptions(
    'host:s{1,}' => \@host_files,
    's:s' => \$clean_stats
);

return 0 if ( -e $hostRemovalPDF);

foreach my $host_file (@host_files)
{
        my ($prefix, $dirs, $suffix) = fileparse($host_file,qr/\.[^.]*/);
        $prefix =~ s/\./_/g;
        push @host_names,  qq("$prefix");
        print "$prefix";
        open(my $fh, "$clean_stats") or die "$!";
        while(<$fh>)
        {
            my ($input_reads) = $_ =~ /Total reads:\s+(\d+)/;
            push @total_reads, $input_reads if (defined $input_reads);
            my ($each_host_reads) = $_ =~ /Total Host reads:\s+(\d+)/;
            if (defined $each_host_reads)
            {
                push @host_reads, $each_host_reads;
                $total_host += $each_host_reads;
            }
        }
        close $fh;
} 
@total_reads =  sort {$a<=>$b} @total_reads;
my $total_reads = pop @total_reads;
    
open (my $ofh, ">$hostclean_stat_file") or die "Cannot write $hostclean_stat_file\n";
print $ofh "Total reads: $total_reads\n";
printf $ofh ("Total non-host reads: %d (%.2f %%)\n", $total_reads - $total_host, ($total_reads - $total_host)/$total_reads*100 );
foreach my $i (0..$#host_names)
{
    printf $ofh ("%s reads: %d (%.2f %%)\n",$host_names[$i],$host_reads[$i],$host_reads[$i]/$total_reads*100);
}
close $ofh;

my $host_names_all = join (',',@host_names);
my $host_reads_all = join (',',@host_reads);
my $Rscript= "$outDir/HostRemoval/hostclean.R";
open(my $Rfh, ">$Rscript") or die "Cannot write $Rscript: $!\n";
print $Rfh <<Rscript;
pdf(file = "$hostRemovalPDF",width = 10, height = 8)
par(xpd=TRUE,mar=c(5,6,4,2))
total<-$total_reads/1000
host<-c($host_reads_all)/1000
hostnames<-c($host_names_all)
mp<-barplot(c(total,host),names.arg=c(\"Total Input\",hostnames),ylab=\"Number of Reads (K)\",col=c(\"gray\",\"red\"))
text(mp,y=c(0,host + 0.01*total) ,c("",sprintf(\"%.2f %%\",(host/total*100) )),pos=3 ) 
title(\"Host Removal Result\")
tmp<-dev.off()
Rscript

close $Rfh;
&executeCommand("R --vanilla --slave --silent < $Rscript 2>/dev/null");
unlink "$Rscript";
die "failed: No reads remain after Host Removal. \n" if ( ($total_reads - $total_host)==0);
return 0;