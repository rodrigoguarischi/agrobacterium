#!/usr/bin/perl

print "Open .gbff file: ";
$gbff = <>;
chomp($gbff);

print "Write the name of the output file: ";
$out = <>;
chomp($out);

open(GBFF,$gbff) or die "Failed to open $gbff\n";
open(OUT,">$out") or die "Failed to write on $out\n";

$/ = "\n//\n";
print OUT "Locus_tag\tOld_locus_tag\tGenomic_position(nucleotides)\tCDS_length\tProduct\n";
while ($gbk = <GBFF>) {
	chomp($gbk);
	$gbk =~ m/LOCUS\s+(.*?) /;
	$locus = $1;
	@features = split("FEATURES",$gbk);
	@cds = split("CDS             ",$features[1]);
	foreach $x (@cds) {
		if ($x =~ m/\/translation/) {
			$x =~ m/\/locus_tag="(.+)"\n/;
			$locus_tag = $1;
			$x =~ m/\/product="(.+)"\n/;
			$product = $1;
			if ($x =~ m/\/old_locus_tag/) {
				$x =~ m/\/old_locus_tag="(.+)"\n/;
				$old_locus_tag = $1;
			}
			else {
				$old_locus_tag = NA;
			}
			if ($x =~ m/^complement/) {
				$x =~ m/^complement\(<?(\d+).+?(\d+)\)\n/;
				$start_cds = $1;
				$end_cds = $2;
			}
			else {
				$x =~ m/^<?(\d+).+?(\d+)\n/;
				$start_cds = $1;
				$end_cds = $2;
			}
			$cds_length = $end_cds-$start_cds + 1;
			print OUT "$locus_tag\t$old_locus_tag\t$locus:$start_cds-$end_cds\t$cds_length\t$product\n";
		}
	}
}
