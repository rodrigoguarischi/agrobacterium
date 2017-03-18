#!/usr/bin/perl
print "Open output from parse_gbff.pl: ";
$gbff = <>;
chomp($gbff);
print "Open .gff file: ";
$gff = <>;
chomp($gff);
print "Write the name of the output file: ";
$out = <>;
chomp($out);

open(GBFF,$gbff) or die "Failed to open $gbff\n";
%ids_gbff;
while ($input_gbff = <GBFF>) {
	chomp($input_gbff);
	@table = split("\t",$input_gbff);
	$ids_gbff{$table[0]} = 1;
}

open(GFF,$gff) or die "Failed to open $gff\n";
open(OUT,">$out") or die "Failed to write on $out\n";

print OUT "Locus_tag\tOld_Locus_Tag\tGenomic_position\tDNA_strand\tCDS_length(nucleotides)\tProduct\n";
while ($line = <GFF>) {
	chomp($line);
	if ($line =~ m/^\#/) {
		next;
	}
	else {
		@tab = split("\t",$line);
		if ($tab[2] eq "gene") {
			@gene_info = split(";",$tab[8]);
			%gene;
			foreach $x (@gene_info) {
				@locus = split("=",$x);
				if ($locus[0] eq "locus_tag") {
					$gene{$locus[0]} = $locus[1];
				}
				if ($locus[0] eq "old_locus_tag") {
					$gene{$locus[0]} = $locus[1];
				}
			}
		}
		if ($line =~ m/\tCDS\t/) {
			$region = $tab[0];
			$start_cds = $tab[3];
			$end_cds = $tab[4];
			$DNA_strand = $tab[6];
			$cds_length = $end_cds - $start_cds + 1;
			$product = $tab[8];
			$product =~ m/product=(.+?)\;/;
			$product = $1;
			$locus_tag = $gene{'locus_tag'};
			if (exists($ids_gbff{$locus_tag})) {
				if (exists($gene{'old_locus_tag'})) {
					print OUT "$gene{'locus_tag'}\t$gene{'old_locus_tag'}\t$region:$start_cds-$end_cds\t$DNA_strand\t$cds_length\t$product\n";
					undef %gene;
				}
				else {
					print OUT "$gene{'locus_tag'}\tNA\t$region:$start_cds-$end_cds\t$DNA_strand\t$cds_length\t$product\n";
					undef %gene;
				}
			}
			else {
				next;
			}
		}
	}
}
