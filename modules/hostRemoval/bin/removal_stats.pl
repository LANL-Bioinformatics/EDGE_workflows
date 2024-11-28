#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;


my @stat_files = ();
my @total_reads;
my @host_reads;
my @host_names;
my $outDir = '.';
my $hostRemovalPDF = "$outDir/HostRemovalStats.pdf";
my $hostclean_stat_file = "$outDir/hostclean.stats.txt";
my $numHostReads;


GetOptions(
    'stats:s{1,}' => \@stat_files,
    'hostReads=s' => \$numHostReads
);

return 0 if ( -e $hostRemovalPDF);


foreach my $stat_file (@stat_files)
{
        my ($prefix, $dirs, $suffix) = fileparse($stat_file,qr/\.[^.]*/);
        $prefix =~ s/\./_/g;
        push @host_names,  qq("$prefix");
        print "$prefix.clean.stats.txt";
        open(my $fh, "$stat_file") or die "$!";
        while(<$fh>)
        {
            my ($input_reads) = $_ =~ /Total reads:\s+(\d+)/;
            push @total_reads, $input_reads if (defined $input_reads);
            my ($each_host_reads) = $_ =~ /Total Host reads:\s+(\d+)/;
            if (defined $each_host_reads)
            {
                push @host_reads, $each_host_reads;
            }
        }
        close $fh;
} 
@total_reads =  sort {$a<=>$b} @total_reads;
my $total_reads = pop @total_reads;
    
open (my $ofh, ">$hostclean_stat_file") or die "Cannot write $hostclean_stat_file\n";
print $ofh "Total reads: $total_reads\n";
printf $ofh ("Total non-host reads: %d (%.2f %%)\n", $total_reads - $numHostReads, ($total_reads - $numHostReads)/$total_reads*100 );
foreach my $i (0..$#host_names)
{
    printf $ofh ("%s reads: %d (%.2f %%)\n",$host_names[$i],$host_reads[$i],$host_reads[$i]/$total_reads*100);
}
close $ofh;

my $host_names_all = join (',',@host_names);
my $host_reads_all = join (',',@host_reads);
my $Rscript= "$outDir/hostclean.R";
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
die "failed: No reads remain after Host Removal. \n" if ( ($total_reads - $numHostReads)==0);


sub executeCommand 
{
    my $command = shift;
    if (system($command) != 0)
         { die ("the command $command failed\n");}
}
