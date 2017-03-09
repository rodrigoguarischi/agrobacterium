#Abrir o arquivo do STDIN e parsear

while ($input = <>) {
	chomp($input);
	@line = split("\t",$input);
	$id = $line[0];
	#separando todas os GO IDs que o gene está associado
	@goID = split(/\|/,$line[13]);
	foreach $x (@goID) {
		print "$id\t$x\n";
	}
} 

#Saida para o STDOUT
