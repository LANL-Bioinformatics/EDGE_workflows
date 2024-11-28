#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;


my $fasta;
my $outputDir;
my $size_filter;
my $max_seq_number;
my $id_mapping;
my $contig_size_for_annotation;
my $project_name;
my $do_annotation;

GetOptions(
    'u=s{1,}' => \$fasta,
    'd=s' => \$outputDir,
    'filt=i' => \$size_filter,
    'maxseq=i' => \$max_seq_number,
    'id:s' => \$id_mapping,
    'ann_size=i' => \$contig_size_for_annotation,
    'n=s' => \$project_name,
    'ann:i' => \$do_annotation #default to false (0)
);

my $output= "$outputDir/${project_name}_contigs.fa";
my $id_map= "$outputDir/id_map.txt";
my $contig_for_annotation = "$outputDir/${project_name}_contigs_${contig_size_for_annotation}up.fa";
$max_seq_number ||= 100000;
my $serial_id= "0" x length($max_seq_number);
my $id_info;
my ($fh,$pid)=open_file($fasta);
open (my $ofh, "> $output") or die "Cannot write $output\n";
open (my $ofh2, "> $contig_for_annotation" ) or die "Cannot write $contig_for_annotation\n";
open (my $idmap_ofh, "> $id_map") or die "Cannot write $id_map\n";
$/ = ">";
while (my $line=<$fh>)
{
        $line =~ s/\>//g;
        my ($id, @seq) = split /\n/, $line;
        next if (!$id);
        ($id_info) = $id =~ /(length_\d+ read_count_\d+)/;
        my $seq = join "", @seq;
        $seq =~ s/-//g;
        $seq =~ s/ //g;
        $seq =~ s/\n//g;
        $seq =~ s/\r//g;
        $seq = uc($seq);
        my $len = length($seq);
        my $GC_num = $seq=~ tr/GCgc/GCgc/; 
        my $GC_content = sprintf("%.2f",$GC_num/$len);
        $id_info = "length_$len "if (!$id_info);
        next if ($len < $size_filter);
        $seq =~ s/(.{70})/$1\n/g; 
        chomp $seq;
        my $fasta_header;
    $id =~ s/\W/_/g;
        if($do_annotation){
            # genbank format limit the LOCUS name length
            if ($id_mapping){
                $fasta_header = ( length($id) > 20 ) ? "contig_$serial_id $id" : "$id contig_$serial_id";
                print $idmap_ofh "contig_$serial_id\t$id\n";
            }else{
                $fasta_header = ( length($project_name) > 20 || $project_name =~/^Assembly/i ) ? "contig_$serial_id $id_info GC_content_$GC_content": "${project_name}_$serial_id $id_info GC_content_$GC_content";
            }
        }else{
        $fasta_header = ($id_mapping)? "$id" : "${project_name}_$serial_id $id_info GC_content_$GC_content";
        }
        if ($len >= $contig_size_for_annotation)
        {
        print $ofh2 ">$fasta_header\n" . $seq."\n";
        }
        print $ofh ">$fasta_header\n" . $seq."\n";
        $serial_id++; 
}    
$/="\n";
close $fh;
close $ofh;
close $ofh2;
#close $idmap_ofh if ($id_mapping && $configuration->{DoAnnotation});
if ( -z $id_map) { unlink $id_map ; }

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