#Abrindo arquivo com anotações do KEGG e vias metabólicas

print "Abrir arquivo com KEGG IDs associados com KEGG Pathways: ";
$arq1 = <>;
chomp($arq1);
open(ARQ1, $arq1) or die "Não foi possível abrir o arquivo $arq1\n";

#Abrindo arquivo de saída

print "Nome do arquivo de saída: ";
$out = <>;
chomp($out);
open(OUT, ">$out") or die "Não foi possível abrir o arquivo de saída $out\n";

#Parseando o arquivo e formatando

while ($line = <ARQ1>) {
	chomp($line);
	@tab = split("\t",$line);
	if ($tab[1] eq "") {
		$keggID = $tab[0];
	}
	else {
		print OUT "$keggID\t$line\n";
	}
}

close(OUT);
close(ARQ1);
