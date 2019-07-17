#!/usr/bin/perl
# Copyright 2019 Rohan Patil
# 
# Tempus Bioinformatics Challenge -
# Parsing VCF file for required data

use strict;
use warnings;
use File::Basename;
use LWP::Simple;

############################################################

print "\nThis script will read VCF file and write out a table in TSV file.\n";
print "\nFORMAT: perl parseVCF.pl VCF_file\n";

my $VCF_file = $ARGV[0] || die;
my($basename, $dirs, $suffix) = fileparse($VCF_file, qr/\.[^.]*/);

print "\nOutput files:";

my $New_VCF = $dirs.$basename."_annotated.vcf";
print "\nAnnotated VCF: $New_VCF";
open (NEWVCF,">$New_VCF");

my $TSV_file = $dirs.$basename."_output_table.tsv";
print "\nVariant Table: $TSV_file\n";
open (VARTAB,">$TSV_file");

############################################################
####  Parse VCF file to Annotated VCF and Output Table  ####
############################################################

my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$normal,$vaf5);
my ($Variant,@VarTab,$temp,$TYPE,$GT,$GQ,$DP,$DPR,$RO,$QR,$AO,$QA,$PAIREDR,$allele_freq);

print VARTAB "Variant\tType\tRead Depth\tNo. of Reads\tRead Percent\tAllele frequency";

open(my $handle1, '<', "$VCF_file") or die "\nCould not open $VCF_file $!";
print "\nParsing VCF ...\nAccessing ExAC ...\n";

while (my $line = <$handle1>)
{	  
	chomp $line;
	
	if($line =~ /\#/) 	
	{ 
		if($line =~ /\#CHROM/) 	{	print NEWVCF "$line\tNewAnnotations\n";	} 	# Adding NewAnnotations column
		else {	print NEWVCF "$line\n";	} 	# Skip lines beginning with "##"
		next;
	} 	
		
	else
	{
		print NEWVCF "$line";
		
		($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$normal,$vaf5) = split(/\t/,$line);	
		
		$Variant = $CHROM."-".$POS."-".$REF."-".$ALT; 		# VariantID (#CHROM-POS-REF-ALT)
		print VARTAB "\n$Variant";
		
		my ($AB,$ABP,$AC,$AF,$AN,$AO,$CIGAR,$DP,$DPB,$DPRA,$EPP,$EPPR,$GTI,$LEN,$MEANALT,$MQM,$MQMR,$NS,$NUMALT,$ODDS,$PAIRED,$PAIREDR,$PAO,$PQA,$PQR,$PRO,$QA,$QR,$RO,$RPL,$RPP,$RPPR,$RPR,$RUN,$SAF,$SAP,$SAR,$SRF,$SRP,$SRR,$TYPE) = split(/;/,$INFO);	
		
		($temp,$TYPE) = split(/\=/,$TYPE); 		# Type of variation.  
		## If there are multiple possibilities, annotate with the most deleterious possibility. ##
		$TYPE =~ s/,.*$//; 				# When comparing two annotations, the “most deleterious” one is shown first.
		print VARTAB "\t$TYPE";
		
		($GT,$GQ,$DP,$DPR,$RO,$QR,$AO,$QA) = split(/:/,$normal);
		print VARTAB "\t$DP"; 			# Depth of sequence coverage at the site of variation. 
		$DPR =~ s/,.*$//;
		print VARTAB "\t$DPR"; 			# No. of Reads supporting the variant. 
		($temp,$PAIREDR) = split(/\=/,$PAIREDR); 	# Proportion of observed reference alleles
		#print "\t$PAIREDR";
		print VARTAB "\t$PAIREDR"; 		# Percentage of reads supporting the variant versus those supporting reference reads.
		print NEWVCF "\tVariantType=$TYPE;ReadDepth=$DP;ReadCount=$DPR;ReadPercent=$PAIREDR;";
		
		## Searching Allele frequency of variant from Broad Institute ExAC Project API ##
		
		my $VariantURL = "http://exac.hms.harvard.edu/rest/variant/".$Variant;
		my $contents = get($VariantURL); 	# This step may take a couple of minutes depending on the system processor.
		my @Split_Content = split(/,/,$contents);
		foreach my $keyword (@Split_Content)
		{
			if ($keyword =~ /\"allele_freq/) {
			($temp, $allele_freq) = split(/:/,$keyword);
			$allele_freq =~ s/\s//;
			print VARTAB "\t$allele_freq"; 		# Allele frequency of variant from Broad Institute ExAC Project API
			print NEWVCF "AlleleFrequency=$allele_freq";
			}
			else 	{ 	next; 	}
		}
		
		print NEWVCF "\n";
	}
}
close(VARTAB);
close(NEWVCF);

############################################################
print "\n###########  DONE  ###########\n";
exit;
