#!/usr/bin/perl -w

use strict;
# star from the gene lists (interatcion lists) from untangel_pairs_v1.pl and the ortholog lists from teh same output using ortholog.conversion hacked R script
# macth the human ensembl IDs with mouse IDs, but only take thsose orthologs that are 1:1, others are NA. 


## how to run: first 2 argumenats are the files below : 1: human_mouse_orthologs.txt from the R script ortholog_conversion_hacked and 
##2: human_secreted_to_receptor.txt from the untangel_orthologs_v2.pl script. Last argument is output file
## perl match_orthologs human_mouse_orthologs.txt complex_curated.csv interaction_curated.csv  gene_input.csv human_secreted_to_receptor.txt

my ($human_mouse_orthologs_file, $human_secreted_to_receptor_file, $output_file, )= @ARGV;

# eg  perl match_orthologs_v2.pl human_mouse_orthologs.txt  human_secreted_to_receptor.txt secreted_to_receptor_with_orthologs.txt
# read in ortholog file
open (OUT, ">$output_file") ||  die "Can't open  $!";

open (IN, "$human_mouse_orthologs_file") ||  die "Can't open  $!";

my %ens_orthologs; # keyed by human ensembl, assuming 1:1. Will remove all those that are multipl by use of other has later
my %symbol_orthologs; # as above but human gene symbol
my %ortholog_count; # keyed by human ensembl and ssimply coundt how many orthologs they have - for remaing stuff later. 

while (<IN>){
		chomp;
		my ($human_symbol, $human_ensembl, $mouse_symbol, $mouse_ensembl) = split;
	#	print "$human_symbol, $human_ensembl, $mouse_symbol, $mouse_ensembl\n";
		$ens_orthologs{$human_ensembl}=$mouse_ensembl;
		$symbol_orthologs{$human_symbol}=$mouse_symbol||"NA";
		$ortholog_count{$human_ensembl}++||1;
		
}
close IN;

#print "test $ens_orthologs{ENSG00000206384}  $symbol_orthologs{COL6A6} \n";# works

print OUT join ( "\t", (  'secreted_ensembl_gene_mouse' , 'secreted_gene_symbol_mouse',
		'secreted_ensembl_gene_human, $secreted_genesymbol_human' , 	
		 	'interacting_receptor_genesymbol_human',	
		'interacting_receptor_ensembl_gene_human', 
'interacting_receptor_genesymbol_mouse', 'interacting_receptor_ensembl_gene_mouse'

		)), "\n";

open (IN2, $human_secreted_to_receptor_file)||  die "Can't open  $!";
while(<IN2>){
	next if /^secreted/;
	my ($secreted_ensembl_gene_human, $secreted_genesymbol_human, $secreted_uniprot_human2, 	
		$interacting_receptor_uniprot_human, $interacting_receptor_genesymbol_human,	
		$interacting_receptor_ensembl_gene_human) = split;

my ($secreted_ensembl_gene_mouse, $secreted_gene_symbol_mouse, $interacting_receptor_ensembl_gene_mouse, $interacting_receptor_genesymbol_mouse);
	# check for multi-orthologs - partner 1

	#print "testing $secreted_ensembl_gene_human \n";
	if ($ortholog_count{$secreted_ensembl_gene_human}){

		if ($ortholog_count{$secreted_ensembl_gene_human} >1){
			$secreted_ensembl_gene_mouse ="NA";
			$secreted_gene_symbol_mouse ="NA";
		}
		else{
			$secreted_ensembl_gene_mouse =$ens_orthologs{$secreted_ensembl_gene_human};
			$secreted_gene_symbol_mouse =$symbol_orthologs{$secreted_genesymbol_human};
		}
	}
	else{ # for cases with no clear ortholog
		$secreted_ensembl_gene_mouse ="NA";
		$secreted_gene_symbol_mouse ="NA";

	}

	# partner 
		#print "testing $interacting_receptor_ensembl_gene_human \n";
	if($ortholog_count{$interacting_receptor_ensembl_gene_human}){
		if ($ortholog_count{$interacting_receptor_ensembl_gene_human} >1){
			$interacting_receptor_ensembl_gene_mouse ="NA";
			$interacting_receptor_genesymbol_mouse ="NA";
		}
		else{
			$interacting_receptor_ensembl_gene_mouse =$ens_orthologs{$interacting_receptor_ensembl_gene_human};
			$interacting_receptor_genesymbol_mouse =$symbol_orthologs{$interacting_receptor_genesymbol_human};
		}
	}
	else{
		$interacting_receptor_ensembl_gene_mouse ="NA";
			$interacting_receptor_genesymbol_mouse ="NA";
	}


	# printout
	print OUT join ( "\t", (  $secreted_ensembl_gene_mouse , $secreted_gene_symbol_mouse,
		$secreted_ensembl_gene_human, $secreted_genesymbol_human,  	
		 	$interacting_receptor_genesymbol_human,	
		$interacting_receptor_ensembl_gene_human, 
$interacting_receptor_genesymbol_mouse, $interacting_receptor_ensembl_gene_mouse

		)), "\n";






}

close (IN2);
close (IN);
close (OUT);

