#!/usr/bin/perl
use warnings; use strict;
no warnings;

# Samuel Moijueh
# BS 845 Applied Statistical Modeling and Programming in R
# December 8, 2014
# Perl script for parsing the openSNP user files (SNP data) and
# creating a matrix with the SNPs of every individual in the study

open(GWAS_file, ">>genodata.txt") or die $!;  # append write mode
open(error_log, ">>error_log.txt") or die $!;
my $string;
my @files;
my $count = 1;

open (my $all_files, "<", "master_list.txt");
while (my $line = <$all_files>) {
    chomp $line;
    push (@files, $line);
}

for my $file (@files) {
    open(my $fh, "$file") or die $!;
    my $input = parse_file($fh, $file);
    print GWAS_file $input."\n";
    print "completed file $count\n";
    $count++;
}


sub parse_file{
    # takes one parameter: (1) the file handle of the file we would like to parse
    my $F1 = $_[0];
    my $file = $_[1];
    my $genome='';
    my @cols;
    
    while (<$F1>){
	$_=~ s/\r//;
	chomp;
	next if /^\s+#/; # skip comment lines
	@cols = split(/\t/, $_);  # split the line elements into an array
	# INSTEAD OF SKIPPING THE --, PUT IN AN '--'
	if ($cols[3] eq "genotype"){
	    next;
	} elsif (length($cols[3]) != 2){
	    print error_log "$file\n$cols[1]\t$cols[3]\n";
	} else{
	    $genome .= $cols[3];
	}
	
	
    }
    
    close($F1);

    
    return join(",",split(//,$genome));
}

close(GWAS_file)
