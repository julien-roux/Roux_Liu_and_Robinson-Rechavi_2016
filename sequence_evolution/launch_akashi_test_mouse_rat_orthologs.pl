#!/usr/bin/env perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Bio::SimpleAlign;
use Bio::AlignIO;
require 'akashi.pm';
$|=1;

# Julien Roux
# This script launches Akashi's test on all pairs of mouse-rat 1-to-1 orthologs, and calculate summary statistics
# This script was developed for Ensembl release 79. I might not work for more recent releases

# Connect to the Ensembl API
my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db(
  -host=>"ensembldb.ensembl.org",
  -user=>"anonymous"
);

# which Ensembl version is used? 
use Bio::EnsEMBL::ApiVersion;
printf( "The API version used is %s\n", software_version() );

my $homology_adaptor = $reg->get_adaptor("Multi", "compara", "Homology");
my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
my $genomedb_adaptor = $reg->get_adaptor("Multi", "compara", "GenomeDB");

my $sp1_gdb = $genomedb_adaptor->fetch_by_registry_name("Mus musculus");
my $sp2_gdb = $genomedb_adaptor->fetch_by_registry_name("Rattus norvegicus");

my $mlss_orth = $mlss_adaptor->fetch_by_method_link_type_GenomeDBs('ENSEMBL_ORTHOLOGUES', [$sp1_gdb, $sp2_gdb]);
my @orthologies = @{$homology_adaptor->fetch_all_by_MethodLinkSpeciesSet($mlss_orth)};

open(OUT, ">", "mouse_akashi_test_scores.txt");
print OUT "member1\tmember2\tZ\tPsi\n";
ORTHOLOGY:
foreach my $orthology (@orthologies){
  if ($orthology->description eq "ortholog_one2one" and $orthology->taxonomy_level() eq 'Murinae'){
    $orthology->print_homology;
    foreach my $member (@{$orthology->get_all_Members}) {
      next ORTHOLOGY if ($member->gene_member->get_Gene->biotype() ne "protein_coding");
      # print gene ID
      print OUT $member->gene_member->stable_id, "\t";
    }

    # Lanch Akashi's test on this gene 
    my ($sum_a, $sum_e_a, $sum_v_a, $sum_ad_n, $sum_bc_n) = akashi_test($orthology);

    if (($sum_v_a ne 0) and ($sum_bc_n ne 0)) {
      my $gene_z = ($sum_a - $sum_e_a)/sqrt($sum_v_a);
      my $gene_psi = $sum_ad_n/$sum_bc_n;
      print OUT "$gene_z\t$gene_psi\n";
    }
    else {
      print OUT "NA\tNA\n";
    }
  }
}
close OUT;
print "Done\n";
exit;

