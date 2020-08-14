#! /usr/bin/perl -w
use strict;
# makes pairs between secreted stuff (ligands?) and whatever they interact with 
# the trick: will untangle complexes on the way

# version 2: takes the file inputs from command line rather than being hardcoded (also output)


# missing: get ensembl gene IDs and then translate those to mouse, fit everything into same table. 
# what is not there, by design, is complexes, they are just unfolded

# ru this, then the ortholog.conversion.hacked (R) , then the 
# mathc_orthologs.pl

## how to run: first 4 argumenats are the files below from cellphonedb, in this order, last argument is output file
## perl untangle_pairs_v2.pl protein_curated.csv  complex_curated.csv interaction_curated.csv  gene_input.csv human_secreted_to_receptor.txt


my ($protein_curated_file, $complex_curated_file, $interaction_curated_file, $gene_input_file, $output_file, )= @ARGV;
# eg: protein_curated.csv  complex_curated.csv interaction_curated.csv  gene_input.csv human_secreted_to_receptor.txt

my %complex1; # hash keyed by proetine 1,2,3,4 {complex name}=1, becasue you can be part of many
my %complex2; # revers, in case
my%is_part_of_complex; # proetin to complex =1
my%complex_contains;
my %receptor; # kedy by complex name, 1 if receptor
my %complex_check;
my %interaction; # doubled hash - partern1, partner2 AND partner  partner 1 ==1 if true
my %protein_to_ensemblID; # keyed by protein, value is gene ID
my %protein_to_geneID;
init(); # read in complex and interaction files, and also traslations to gene IDs


# get secreted proteins


# tester
#my @interacting_proteins= interactions_from_protein('Q9Y6F9');
#print "test, " , "@interacting_proteins";
#print  "\nlength: ", scalar(@interacting_proteins), "\n";


#foreach my $interacting_protein(@interacting_proteins){
#		
#			print "and interacts with $interacting_protein\n";
#			print join("\t", ($interacting_protein), "test") , "\n";
#
#			}
#		
open OUT, ">$output_file" || die  "Can't open  $!";

print OUT join("\t", ( "secreted.ensembl.gene.human", "secreted.genesymbol.human", "secreted.uniprot.human2", 
	"interacting.receptor.uniprot.human", "interacting.receptor.genesymbol.human", "interacting.receptor.ensembl.gene.human")), "\n";
open(IN0, "<", $protein_curated_file)||  die "Can't open  $!";

my %seen; # just to avoid repetaing same interatcions

while (<IN0>){
	chomp;
	next unless $_;

	next if /^uniprot/;
	# this is csv
	my @words = split ",", $_;
	
	my ($acc,$uniprot,$secreted, $is_receptor)=@words[0,1,4,7];

	#print "$acc,$uniprot,$secreted, $is_receptor\n";
	if ($secreted eq "True"){
		# find stuff.
	#	print "$acc $uniprot is secreted"; 
	my @interacting_proteins= interactions_from_protein($acc);
	next unless @interacting_proteins;
	foreach my $interacting_protein(@interacting_proteins){
		#next if $interacting_protein eq " ";
			#print "and interacts with $interacting_protein\n";
			# extrac check: has to be a recptor 
			next unless $receptor{$interacting_protein};
			next if $seen{$acc}{$interacting_protein};
			print OUT  join("\t", ($protein_to_ensemblID{$acc}, $protein_to_geneID{$acc} , $acc , 
			 $interacting_protein, $protein_to_geneID{$interacting_protein}, 
			 $protein_to_ensemblID{$interacting_protein}, )) , "\n";
			$seen{$acc}{$interacting_protein}=1;
		}

	}



}

close OUT;
close IN0;

# for each of those, find interactions (that are receptors? later)








### SUBS

sub init{




	open(IN4, "<", $gene_input_file)||  die "Can't open  $!";
while (<IN4>){
	chomp;
	next unless $_;

	next if /^gene_name/;
	# this is csv
	my ($gene_name, $uniprot, $symbol, $ensembl) = split ",", $_;
#print "$gene_name, $uniprot, $symbol, $ensembl\n";
	$protein_to_ensemblID{$uniprot}=$ensembl;
	$protein_to_geneID{$uniprot}=$symbol;

}


close IN4;



	open(IN3, "<", $protein_curated_file)||  die "Can't open  $!";

while (<IN3>){
	chomp;
	next unless $_;

	next if /^uniprot/;
	# this is csv
	my @words = split ",", $_;
	
	my ($acc,$uniprot,$secreted, $is_receptor)=@words[0,1,4,7];
	$receptor{$acc}=1 if $is_receptor eq "True";
}
close IN3;

	open(IN1, "<", $complex_curated_file)||  die "Can't open  $!";

while (<IN1>){
	chomp;
	next unless $_;
	# this is csv
	my @words = split ",", $_;

	my ($name,$p1,$p2,$p3,$p4, $is_receptor)=@words[0,1,2,3,4, 10 ];  
	$receptor{$name}=1 if $is_receptor eq "TRUE";
	$complex_check{$name}=1;

	#print "$name,$p1,$p2,$p3,$p4, $is_receptor \n";
	foreach  my $i ($p1,$p2,$p3,$p4){
		next unless $i;
		$is_part_of_complex{$i}{$name}=1;
		$complex_contains{$name}{$i}=1;

#		$complex1{$i}{$name}=1;
#		$complex2{$name}{$i}=1; # just in case
	}

}
close (IN1);


# read in interaction list


open(IN2, "<", $interaction_curated_file)||  die "Can't open  $!";

while (<IN2>){
	chomp;
	next unless $_ ;
	# this is csv
	my @words = split ",", $_;

	my ($interaction_id, $partner1, $partner2)=@words[0,1,2];
	$interaction{$partner1}{$partner2}=$interaction_id;
 	$interaction{$partner2}{$partner1}=$interaction_id;

}

close(IN2);


}
##









#sub for chking if something IS a complex

sub is_complex{
	return 1 if $complex_check{$_[0]};
}

#print "test ",  is_complex("PlexinA1_complex3");


# sub for getting complex(es) for a given protein
sub get_complexes_from_protein{
	#returns array array of complexes
	my $protein=$_[0];
	my @complexes;
	foreach my $complex ( keys %{$is_part_of_complex{$protein}}){
		push(@complexes, $complex)
	}

return @complexes;

}
#print "testing P51172\n";
#my @test= get_complexes_from_protein("P07949");
#print "@test ";

# sub for getting proteins for a given complex
sub proteins_from_complex{
	my $complex=$_[0];
	my @proteins;
		foreach my $protein ( keys %{$complex_contains{$complex}}){
		push(@proteins, $protein)
	}

return @proteins;
}

# test
#my @res=proteins_from_complex("AT8B4CC50A complex");

#print "@res "; 


# sub for getting interactions starting from protein
sub interactions_from_protein{
	my $protein1=$_[0];
	# get interactions as protein, if any
	my @direct_interactions;
#testing both directionsbeasue has is 'doubled'
	foreach my $partner2 (keys %{$interaction{$protein1} }){

		# is partner 2 a complex?
		#print "$protein1 interacting  direct with with $partner2 which is ";
		if (is_complex($partner2)){
			# if so, get the member of that complex
	#		print "a complex containig: \n";
		 
			foreach my $protein2 (proteins_from_complex($partner2)){
	#			print "\t $protein2 \n";

				push(@direct_interactions, $protein2)

			}
		}
		else{
	 #print "a protein\n";
			push(@direct_interactions, $partner2)


		}
	}
	

#print  "tetsing length: ", scalar(@direct_interactions);



return (@direct_interactions);
	
}

sub interactions_from_complex{
# will this ever be neeeded? not implemenetd yet
}


	# test:
#	print "test P24821 ";#
#	my @test= interactions_from_protein('P24821');
#	print "@test";

#rint "test  P10145";
##	print "@test";



# sub for gettig interactins starting  from complex
#TODO








###### old
# test: works
#foreach my $protein(keys %complex1){
#	print "protein is $protein\n";
#	foreach my $name ( keys %{ $complex1{$protein}  }){
#		print "\t complex is $name\n";

#	}
#	print "\n";
#	
#}

#foreach my $protein(keys %complex2){
#	print "complex is $protein\n";
#	foreach my $name ( keys %{ $complex2{$protein}  }){
#		print "\t protein is $name\n";

#	}
#	print "\n";
#	
#}


#1: get secreted proteins





# find if these are in complex: if so, find what complex is interatcing with 

