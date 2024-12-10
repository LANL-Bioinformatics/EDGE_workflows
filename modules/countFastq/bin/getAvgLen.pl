#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $outputDir='';

my @unpairedList= ();
my @pairedList=  ();

GetOptions(
    'u:s{1,}' => \@unpairedList,
    'p:s{1,}' => \@pairedList,
    'd=s' => \$outputDir,
);


my $all_R1_fastq = "$outputDir/all.1.fastq";
my $all_R2_fastq = "$outputDir/all.2.fastq";
my $all_SE_fastq = "$outputDir/all.se.fastq";
my $count_file_list= "$outputDir/fastqCount.txt";
my $avg_read_len= 0 ;
my $PE_count = 0 ;
my $PE_total_len = 0;
my $SE_count= 0 ;
my $SE_total_len = 0 ;



open (my $fh, ">$count_file_list") or die "Cannot write $count_file_list\n";
while(my ($R1,$R2) = splice (@pairedList,0,2))
{
    ($PE_count, $PE_total_len) = &countPE_exe($R1,$R2,$all_R1_fastq,$all_R2_fastq,$all_SE_fastq,$fh);
}
foreach my $file (@unpairedList)
{
    ($SE_count,$SE_total_len)=&countFastq_exe($file,$all_SE_fastq);
    printf $fh ("%s\t%d\t%d\t%.2f\n",basename($file),$SE_count,$SE_total_len,$SE_total_len/$SE_count);
    printf ("%s\t%d\t%d\t%.2f\n",basename($file),$SE_count,$SE_total_len,$SE_total_len/$SE_count);
}
close $fh;
$avg_read_len = ($PE_count + $SE_count) > 0 ? ($PE_total_len + $SE_total_len) / ($PE_count + $SE_count) : 0 ;

sub countPE_exe{
    my $r1=shift;
    my $r2=shift;
    my $out_r1=shift;
    my $out_r2=shift;
    my $out_se=shift;
    my $count_fh=shift;
    my %seq_hash;
    my $pair_char;
    my $unpaired_count=0;
    my $read1_count=0;
    my $read2_count=0;
    my $se2_count=0;
    my $se1_count=0;
    my $paired_count=0;
    my $read1_total_len=0;
    my $read2_total_len=0;
    my $existed_id1=0;
    my $existed_id2=0;
    my ($fh1,$pid) = open_file($r1);
    open (my $ofh1, ">>$out_r1") or die "Cannot write $out_r1\n";
    open (my $ofh2, ">>$out_r2") or die "Cannot write $out_r2\n";
    open (my $ofhse, ">>$out_se") or die "Cannot write $out_se\n";
    while(<$fh1>){
        chomp;
        next unless $_ =~ /\S/;
        next if ($_ =~ /length=0/);
        my $id_line=$_;
        my ($id) = $id_line =~ /^\@(\S+).?\/?1?\s*/;
        my $seq = <$fh1>;
        chomp $seq;
        if ($seq_hash{$id}){
            $existed_id1++;
        }
        my $len = length $seq;
        $read1_total_len += $len;
        my $qual_id = <$fh1>;
        my $qual = <$fh1>;
        $seq = $seq."\n".$qual_id.$qual;
        $seq_hash{$id}++;
        $read1_count++;
    }
    close $fh1;
    my %seq_hash2;
    my ($fh2,$pid2) = open_file($r2);
    while(<$fh2>){
        chomp;
        next unless $_ =~ /\S/;
        next if ($_ =~ /length=0/);
        my $id_line=$_;
        my ($id2) = $id_line =~ /^\@(\S+)\.?\/?2?\s*/;
        $read2_count++;
        my $seq2 = <$fh2>;
        chomp $seq2;
        if ($seq_hash2{$id2}){
            $existed_id2++;
        }
        my $len = length $seq2;
        $read2_total_len += $len;
        my $qual_id = <$fh2>;
        my $qual = <$fh2>;
        $seq2 = $seq2."\n".$qual_id.$qual;
        $seq_hash2{$id2}++;
        if ($seq_hash{$id2}){
            $seq_hash{$id2}++;
            $paired_count++;
            print $ofh2 $id_line,"\n",$seq2;
        }else{
            print $ofhse $id_line,"\n",$seq2;
            $se2_count++;
        }
    }
    close $fh2;
    ($fh1,$pid) = open_file($r1);
    while(<$fh1>){
        chomp;
        next unless $_ =~ /\S/;
        next if ($_ =~ /length=0/);
        my $id_line=$_;
        my ($id) = $id_line =~ /^\@(\S+)\.?\/?1?\s*/;
        my $seq = <$fh1>;
        chomp $seq;
        my $qual_id = <$fh1>;
        my $qual = <$fh1>;
        $seq = $seq."\n".$qual_id.$qual;
        if ($seq_hash{$id} == 2){
                print $ofh1 $id_line,"\n",$seq;
        }
        if ($seq_hash{$id} == 1){
                print $ofhse $id_line,"\n",$seq;
                $se1_count++;
        }
    }
    close $fh1;
    close $ofh1;
    close $ofh2;
    close $ofhse;
    printf ("%s\t%d\t%d\t%.2f\n",basename($r1),$read1_count,$read1_total_len,$read1_total_len/$read1_count);
    printf ("%s\t%d\t%d\t%.2f\n",basename($r2),$read2_count,$read2_total_len,$read2_total_len/$read2_count);
    printf $count_fh ("%s\t%d\t%d\t%.2f\n",basename($r1),$read1_count,$read1_total_len,$read1_total_len/$read1_count);
    printf $count_fh ("%s\t%d\t%d\t%.2f\n",basename($r2),$read2_count,$read2_total_len,$read2_total_len/$read2_count);
    printf ("%d duplicate id from %s\n", $existed_id1, basename($r1)) if ($existed_id1 > 0);
    printf ("%d duplicate id from %s\n", $existed_id2, basename($r2)) if ($existed_id2 > 0);
    printf ("There are %d reads from %s don't have corresponding paired read.\n", $se1_count, basename($r1)) if ($se1_count >0); 
    printf ("There are %d reads from %s don't have corresponding paired read.\n", $se2_count, basename($r2)) if ($se2_count >0);
   
    unlink $out_se if (-z $out_se);
    return ($read1_count + $read2_count, $read1_total_len + $read2_total_len);
}

sub countFastq_exe 
{
    my $file=shift;
    my $output=shift;
    my $seq_count=0;
    my $total_length;
    my ($fh,$pid)= open_file($file);
    open (my $ofh, ">>$output") or die "Cannot write $output\n";
    while (<$fh>)
    { 
        next unless $_ =~ /\S/;
	next if ($_ =~ /length=0/);
        my $id=$_;
        $id = '@'."seq_$seq_count\n" if ($id =~ /No name/);
        my $seq=<$fh>;
        chomp $seq;
        my $q_id=<$fh>;
        my $q_seq=<$fh>;
        my $len = length $seq;
        $seq_count++;
        $total_length +=$len;
        print $ofh "$id$seq\n$q_id$q_seq";
    }
    close $fh;
    return ($seq_count,$total_length);
}

sub touchFile{
    my $file=shift;
    open (my $fh,">",$file) or die "$!";
    close $fh;
}

sub open_file
{
    my ($file) = @_;
    print "$file\n"; 
    my $fh;
    my $pid;
    if ( $file=~/\.gz\.?\d?$/i ) { $pid=open($fh, "gunzip -c $file |") or die ("gunzip -c $file: $!"); }
    else { $pid=open($fh,'<',$file) or die("$file: $!"); }
    return ($fh,$pid);
}