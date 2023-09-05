sub runHostRemoval 
{
    my $pairFile=shift;
    my $unpairFile=shift;
    my $host=shift;
    my $time=time();
    my $min_score = $configuration->{"bwaMemOptions"} || "-T 50 ";
    my $ont_flag = ($configuration->{"fastq_source"} =~ /nanopore/)? "-x ont2d ": ($configuration->{"fastq_source"} =~ /pacbio/i)? "-x pacbio ": "";  
    $min_score = "-T $configuration->{min_L}" if $ont_flag;
    my $similarity_cutoff = $configuration->{"similarity"} || 90 ;
   # my $stats = "$QCoutDir/hostclean.stats.txt";
    my ($prefix, $dirs, $suffix) = fileparse($host,qr/\.[^.]*/);
    my $host_abs_path = Cwd::abs_path("$host");
    $prefix =~ s/\./_/g;
    
    my $outputDir = "$abs_outDir/HostRemoval/$prefix";
    &make_dir($outputDir);
    if ($noColorLog)
    {
        &lprint ("[Host Removal]\n $prefix\n");
    }
    else
    {
        &lprint (colored ("[Host Removal]\n $prefix\n",'yellow'));
    }
    
    if ( ! -e $host){
    	&lprint("No host file input. Skip analysis");
    	return($pairFile,$unpairFile);
    }
    my $host_index_dir=$outputDir;
    unless (  -e "$host.bwt" or ($host =~ /edge_ui/ and -w $host))
    {
        symlink("$host_abs_path", "$host_index_dir/$prefix.fa");   
        $host =   "$host_index_dir/$prefix.fa";  
    }
    
    $prefix = "$prefix.clean";
    if ( -s "$outputDir/$prefix.1.fastq" && -e "$outputDir/run${prefix}.finished" )
    {
        &lprint ("Host Removal Finished\n");
        symlink("$outputDir/$prefix.1.fastq", "$abs_outDir/HostRemoval/hostclean.1.fastq");
        symlink("$outputDir/$prefix.2.fastq", "$abs_outDir/HostRemoval/hostclean.2.fastq");
        symlink("$outputDir/$prefix.unpaired.fastq", "$abs_outDir/HostRemoval/hostclean.unpaired.fastq");
        return ("$outputDir/$prefix.1.fastq $outputDir/$prefix.2.fastq","$outputDir/$prefix.unpaired.fastq");
    }
    elsif ( -s "$outputDir/$prefix.unpaired.fastq" && -e "$outputDir/run${prefix}.finished" )
    {
        &lprint ("Host Removal Finished\n");
        symlink("$outputDir/$prefix.unpaired.fastq", "$abs_outDir/HostRemoval/hostclean.unpaired.fastq");
        return ("","$outputDir/$prefix.unpaired.fastq");
    }
    unlink "$outputDir/run${prefix}Removal.finished";
    unlink "$abs_outDir/HostRemoval/HostRemovalStats.pdf";
    unlink glob("$abs_outDir/HostRemoval/hostclean*");
  
    
    my $parameters;
    $parameters .= " -p $pairFile" if ($pairFile);
    $parameters .= " -u $unpairFile" if ( -s $unpairFile);
    $parameters .= " -ref ". $host;
    $parameters .= " -bwaMemOptions \"$min_score $ont_flag\"";
    $parameters .= " -s $similarity_cutoff";
    $parameters .= " -o $outputDir -cpu $numCPU -host";
    $parameters .= " -prefix $prefix ";
    my $command = "perl $RealBin/scripts/host_reads_removal_by_mapping.pl $parameters ";
    &lprint ("  Running \n  $command \n");
    &executeCommand($command);
    &printRunTime($time);
    &touchFile("$outputDir/run${prefix}.finished");
    unlink glob("$host_index_dir/$prefix.fa.*");
    if ( -s "$outputDir/$prefix.1.fastq")
    {
        symlink("$outputDir/$prefix.1.fastq", "$abs_outDir/HostRemoval/hostclean.1.fastq");
        symlink("$outputDir/$prefix.2.fastq", "$abs_outDir/HostRemoval/hostclean.2.fastq");
        symlink("$outputDir/$prefix.unpaired.fastq", "$abs_outDir/HostRemoval/hostclean.unpaired.fastq");
        return ("$outputDir/$prefix.1.fastq $outputDir/$prefix.2.fastq","$outputDir/$prefix.unpaired.fastq");
    }
    elsif( -s "$outputDir/$prefix.unpaired.fastq" ){
        symlink("$outputDir/$prefix.unpaired.fastq", "$abs_outDir/HostRemoval/hostclean.unpaired.fastq");
        return ("","$outputDir/$prefix.unpaired.fastq");
    }
   
}

sub runHostRemovalStat
{
    my @host_files = @{$configuration->{Host}};
    my @total_reads;
    my @host_reads;
    my @host_names;
    my $hostReomvalPDF = "$outDir/HostRemoval/HostRemovalStats.pdf";
    my $hostclean_stat_file = "$outDir/HostRemoval/hostclean.stats.txt";
    my $total_host=0;
    
    return 0 if ( -e $hostReomvalPDF);
    
    &lprint ("  Host Removal Stat and plot\n");
    
    foreach my $host_file (@host_files)
    {
          my ($prefix, $dirs, $suffix) = fileparse($host_file,qr/\.[^.]*/);
          $prefix =~ s/\./_/g;
          push @host_names,  qq("$prefix");
          open(my $fh, "$outDir/HostRemoval/$prefix/$prefix.clean.stats.txt") or die "$!";
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
pdf(file = "$hostReomvalPDF",width = 10, height = 8)
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
}
