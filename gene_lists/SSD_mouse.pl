#!/usr/bin/env perl -w
use strict;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Registry;
$|=1; #no buffering of the output

# Julien Roux and Jialin Liu
# This script looks for mouse small-scale duplicates. It was developed for Ensembl release 79. I might not work for more recent releases

# Connect to the Ensembl API
my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db(
  -host=>"ensembldb.ensembl.org",
  -user=>"anonymous",
);

# which Ensembl version is used?
use Bio::EnsEMBL::ApiVersion;
printf( "The API version used is %s\n", software_version() );

# Get all adaptors needed
my $gene_tree_adaptor = $reg->get_adaptor("Multi", "compara", "GeneTree");
my $ncbi_taxon_adaptor = $reg->get_adaptor("Multi", "compara", "NCBITaxon");
my $homology_adaptor = $reg->get_adaptor("Multi", "compara", "Homology");
my $member_adaptor = $reg->get_adaptor("Multi", "compara", "Member");

# In a release database, there is a basal root that holds together all the trees.
my @children = @{$gene_tree_adaptor->fetch_all(-tree_type     => 'tree',
                                               -member_type   => 'protein',
                                               -clusterset_id => 'default',
                                           )}; # Protein tree only !


# To test the taxonomical level of the gene trees, we need first to define the deepest taxonomical level we want to keep
# Mammalia: tax_id = 40674
my $limit_node = $ncbi_taxon_adaptor->fetch_node_by_taxon_id('40674');
my $limit_left = $limit_node->left_index;
my $limit_right = $limit_node->right_index;

TREE:
foreach my $tree (sort {$a->stable_id cmp $b->stable_id} (@children)) { #Sort to always get the same order

  # Tip: will make the traversal of the tree faster
  $tree->preload();

  # Testing the taxonomical level of the gene tree, by calling the recursive subroutine
  select_sub_tree($tree->root, 0);
  $tree->root->release_tree;
}
print "Done\n";
exit;

# recursive subroutine ############################################
sub select_sub_tree {
  my ($this_tree, $level) = @_;
  print "Tree root ID: ", $this_tree->node_id, ", stable ID: ", $this_tree->tree->stable_id, ", level: ", $level, "\n";

  # get taxonomical level of the root
  my $tax_level = $this_tree->taxonomy_level(); 
  print "\tTaxonomy level: $tax_level\n";

  if ((defined $tax_level) and ($tax_level ne '')) {

    # check if the taxonomic level of the deepest node is lower/higher than the limit (using left and right indexes)
    my $deepest_node_tax = $ncbi_taxon_adaptor->fetch_node_by_name($tax_level);
    my $left_tax;
    my $right_tax;
    if (defined $deepest_node_tax) {
      $left_tax = $deepest_node_tax->left_index;
      $right_tax = $deepest_node_tax->right_index;

      # the tree has a higher taxonomic level -> reject!
      if (($left_tax > $limit_left) and ($right_tax < $limit_right)) {
        print "\tReject!\n";
      }

      # good taxonomic level, speciation node -> Test!
      elsif (($left_tax == $limit_left) and ($right_tax == $limit_right) and ($this_tree->get_tagvalue('node_type') eq 'speciation')) {
        print "\tTest!\n";
        test_SSD($this_tree);
      }

      # good taxonomic level + duplication node -> Find subtree!
      elsif (($left_tax == $limit_left) and ($right_tax == $limit_right) and ($this_tree->get_tagvalue('node_type') eq 'duplication')) {
        print "\tFind subtree (duplication)!\n";

        # retrieve all children for this node
        my @children = @{$this_tree->children()};
        $level++;
        foreach my $children_tree (@children) {
          select_sub_tree($children_tree, $level);
        }
      }

      # lower taxonomic level -> Find subtree!
      elsif (($left_tax < $limit_left) and ($right_tax > $limit_right)) {
        print "\tFind subtree!\n";
        
        # retrieve all children for this node
        my @children = @{$this_tree->children()};
        $level++;
        foreach my $children_tree (@children) {
          select_sub_tree($children_tree, $level);
        }
      }
      else {
        print "\tOut of the range!\n";
      }
    }
    else {
      print "\tCan't get the corresponding NCBI taxon! (ncbi_taxon_adaptor)\n";
    }
  }
  else {
    print "\tTaxonomic level is not defined!\n"
  }
  $this_tree->release_tree;
  return;
}
########################################################################################
# get all mouse paralogs that arose in mammals. Remove all homologies from mouse-specific duplications + low-confidence duplications + cases were paralogs are 100% identical or <10% identical (probably gene split)
sub test_SSD {
  my ($root) = @_;

  my %flaggedGenes;
  my @homologies; # to store all high-confidence homologies

  # Retrieve all mouse genes from mammalian subtree
  foreach my $node (@{$root->get_all_nodes()}) {
    if ($node->is_leaf() and ($node->taxonomy_level() eq 'Mus musculus')) {
      print "\tMouse gene: ", $node->gene_member->stable_id(), "\n";
      my $mouse_gene = $node->gene_member;

      # check all homologies and flag genes originating from low-confidence duplications / mouse-specific duplciations
    HOMOLOGY:
      foreach my $homology (@{$homology_adaptor->fetch_all_by_Member($mouse_gene, 'ENSEMBL_PARALOGUES')}) {
        # if both genes are already flagged
        next HOMOLOGY if ((exists $flaggedGenes{${$homology->gene_list}[0]->stable_id()}) and (exists $flaggedGenes{${$homology->gene_list}[1]->stable_id()}));

        # keep mammalian duplications
        next HOMOLOGY unless (($homology->taxonomy_level() eq 'Theria') or ($homology->taxonomy_level() eq 'Eutheria') or ($homology->taxonomy_level() eq 'Boreoeutheria') or ($homology->taxonomy_level() eq 'Euarchontoglires') or ($homology->taxonomy_level() eq 'Glires') or ($homology->taxonomy_level() eq 'Rodentia') or ($homology->taxonomy_level() eq 'Sciurognathi') or ($homology->taxonomy_level() eq 'Murinae') or ($homology->taxonomy_level() eq 'Mus musculus'));

        # flag mouse-specific duplicates
        if ($homology->taxonomy_level() eq 'Mus musculus') {
          $flaggedGenes{${$homology->gene_list}[0]->stable_id()}++;
          $flaggedGenes{${$homology->gene_list}[1]->stable_id()}++;
          next HOMOLOGY;
        }

        #flag low confidence duplciations
        if ($homology->gene_tree_node->get_tagvalue("duplication_confidence_score") lt 0.5) {
          $flaggedGenes{${$homology->gene_list}[0]->stable_id()}++;
          $flaggedGenes{${$homology->gene_list}[1]->stable_id()}++;
          next HOMOLOGY;
        }

        # check that 2 paralogs are not 100% identical, and are more than 10% identical (cases of gene split)
        if (($homology->get_SimpleAlign->percentage_identity() eq 100) or ($homology->get_SimpleAlign->percentage_identity() le 10)) {
          $flaggedGenes{${$homology->gene_list}[0]->stable_id()}++;
          $flaggedGenes{${$homology->gene_list}[1]->stable_id()}++;
          next HOMOLOGY;
        }

        # if all controls are fine
        push(@homologies, $homology);
      }
    }
  }

  my %reportedHomologies;
  print "Reporting homologies:\n";
  foreach my $homology (@homologies){
    my $gene1 = ${$homology->gene_list}[0]->stable_id();
    my $gene2 = ${$homology->gene_list}[1]->stable_id();

    # do not report a homology that was already printed out
    unless ((exists $reportedHomologies{$gene1}->{$gene2}) or (exists $reportedHomologies{$gene2}->{$gene1})){
      $homology->print_homology();

      open(OUT, '>>mouse_SSD.txt') or die "Can't open file";
      print OUT $gene1, "\t", $gene2, "\t", $homology->taxonomy_level(), "\t", $root->tree->stable_id, "\n";
      close OUT;
      $reportedHomologies{$gene1}->{$gene2}++;
    }
  }
  return;
}


