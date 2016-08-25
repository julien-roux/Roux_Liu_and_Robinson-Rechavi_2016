#!/usr/bin/env perl -w
use strict;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Registry;
$|=1; #no buffering of the output

# Julien Roux and Jialin Liu
# This script looks for Ensembl gene trees with topology compatible or not with the 2R vertebrate whole-genome duplication. It was developed for Ensembl release 79. I might not work for more recent releases

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

# In a release database, there is a basal root that holds together all the trees.
my @children = @{$gene_tree_adaptor->fetch_all(-tree_type     => 'tree',
                                               -member_type   => 'protein',
                                               -clusterset_id => 'default',
                                           )};

# To test the taxonomical level of the gene trees, we need first to define the appropriate range of  taxonomical level we want to keep
# The 2R duplication occurred after chordata, but before divergence of vertebrata (lamprey)
# Euteleostomi: 117571
# Vertebrata: 7742
# Chordata: 7711
# Bilateria: 33213

my $down_limit_node = $ncbi_taxon_adaptor->fetch_node_by_taxon_id('33213');
my $up_limit_node = $ncbi_taxon_adaptor->fetch_node_by_taxon_id('7711');
my $down_limit_left = $down_limit_node->left_index;
my $down_limit_right = $down_limit_node->right_index;
my $up_limit_left = $up_limit_node->left_index;
my $up_limit_right = $up_limit_node->right_index;

TREE:
foreach my $tree (sort {$a->stable_id cmp $b->stable_id} (@children)) { 
  # Tip: will make the traversal of the tree faster
  $tree->preload();
  print "Tree root:", $tree->root, "\n", $tree, "\n";

  # Testing the taxonomical level of the gene tree, by calling the recursive subroutine
  select_sub_tree($tree->root, 0);
  $tree->root->release_tree;
}
print "Done\n";
exit;

## recursive subroutine ############################################
sub select_sub_tree {
  my ($this_tree, $level) = @_;
  print "Tree root ID: ", $this_tree->node_id, ", level: ", $level, "\n";

  # get taxonomical level of the root
  my $tax_level = $this_tree->taxonomy_level(); 
  print "\tTaxonomy level: $tax_level\n";

  #check if the taxonomic level of the deepest node is lower/higher than the limit
  if ((defined $tax_level) and ($tax_level ne '')) {
    my $deepest_node_tax = $ncbi_taxon_adaptor->fetch_node_by_name($tax_level);
    my $left_tax;
    my $right_tax;
    if (defined $deepest_node_tax) {
      $left_tax = $deepest_node_tax->left_index;
      $right_tax = $deepest_node_tax->right_index;

      #limit taxonomic level (between bilateria and chordata), speciation node -> Test! or Find subtree!
      if (($left_tax <= $up_limit_left) and ($right_tax >= $up_limit_right) and ($left_tax >= $down_limit_left) and ($right_tax <= $down_limit_right) and ($this_tree->get_tagvalue('node_type') eq 'speciation')) {
        foreach my $child1 (@{$this_tree->children()}) {

          # if next node is vertebrata or euteleostomi -> Test!
          if (($child1->taxonomy_level() eq 'Vertebrata') or ($child1->taxonomy_level() eq 'Euteleostomi')) {
            print "\tTest! (", $child1->tree->stable_id, ")\n";
            test_2R_dup_mouse($child1);
            test_2R_sing_mouse($child1);
          }

          # if next node is chordata -> Find subtree!
          elsif ($child1->taxonomy_level() eq 'Chordata') {
            print "\tFind subtree (child chordata)!\n";
            # retrieve all children for this node
            $level++; 
            select_sub_tree($child1, $level);
          }
        }
      }

      # higher taxonomic level -> reject!
      elsif (($left_tax > $up_limit_left) and ($right_tax < $up_limit_right)) {
        print "\tToo high!\n";
      }

      #limit taxonomic level (between chordata and vertebrata), duplication node -> Find subtree !
      elsif (($left_tax <= $up_limit_left) and ($right_tax >= $up_limit_right) and ($left_tax >= $down_limit_left) and ($right_tax <= $down_limit_right) and ($this_tree->get_tagvalue('node_type') eq 'duplication')) {
        print "\tFind subtree (duplication)!\n";
        # retrieve all children for this node
        my @children = @{$this_tree->children()};
        $level++; 
        foreach my $children_tree (@children) {
          select_sub_tree( $children_tree, $level);
        }
      }

      # lower taxonomic level -> Find subtree!
      elsif (($left_tax < $down_limit_left) and ($right_tax > $down_limit_right)) {
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

####################################################################
sub test_2R_dup_mouse {
  my ($root) = @_;

  # We want to retrieve 1R and/or 2R duplicates
  # The duplication can be dated vertebrata (if lamprey present) or euteleostomi (if lamprey gene absent)
  # To be sure this is not a fake duplication due to bad clustering of some vertebrate genes (e.g., a few genes cluster out of vertebrate lineage, leading tree reconciliation to infer a duplication), we ask for confidence score of duplications >= 50%

  # test the first duplication node (1R)
  if ((($root->taxonomy_level() eq 'Vertebrata') or ($root->taxonomy_level() eq 'Euteleostomi')) and ($root->get_tagvalue('node_type') eq 'duplication') and ($root->get_tagvalue("duplication_confidence_score") ge 0.5)) {
    print "\tThis tree has 1R dup\n";
    foreach my $child1 (@{$root->children()}) {
      # direct speciation after the duplication (1R or 2R was lost)
      if ((($child1->taxonomy_level() eq 'Vertebrata') or ($child1->taxonomy_level() eq 'Euteleostomi')) and ($child1->get_tagvalue('node_type') eq 'speciation')) {
        check_euteleostomi_group_strict_mouse($child1, "mouse_2R_ohnologs.txt", $root->node_id());
      }
      # two rounds of WGD (1R+2R)
      if ((($child1->taxonomy_level() eq 'Vertebrata') or ($child1->taxonomy_level() eq 'Euteleostomi')) and ($child1->get_tagvalue('node_type') eq 'duplication') and ($child1->get_tagvalue("duplication_confidence_score") ge 0.5)) {
        print "\tThis tree has 2R dup\n";
        foreach my $child2 (@{$child1->children()}) {
          if ((($child2->taxonomy_level() eq 'Vertebrata') or ($child2->taxonomy_level() eq 'Euteleostomi')) and ($child2->get_tagvalue('node_type') eq 'speciation')) {
            check_euteleostomi_group_strict_mouse($child2, "mouse_2R_ohnologs.txt", $root->node_id());
          }
        }
      }
    }
  }
  return;
}

####################################################################
sub test_2R_sing_mouse {
  my ($root) = @_;

  # We want to retrieve 2R singletons in danio (no 2R, nor 3R WGD)
  # The speciation can be dated vertebrata (if lamprey present) or euteleostomi (if lamprey gene absent)

  # test the first node
  if ((($root->taxonomy_level() eq 'Vertebrata') or ($root->taxonomy_level() eq 'Euteleostomi')) and ($root->get_tagvalue('node_type') eq 'speciation')) {
    print "\tThis tree has no 1R/2R dup\n";
    check_euteleostomi_group_strict_mouse($root, "mouse_2R_singletons.txt", $root->node_id());
  }
  return;
}

###################################################################
sub check_euteleostomi_group_strict_mouse {
  my ($root, $filename, $node_id) = @_;

  # We select the group if:
  # - no duplication in the tree on branches leading to mouse (taxonomic levels: Vertebrata, Euteleostomi, Sarcopterygii, Tetrapoda, Amniota, Mammalia, Theria, Eutheria, Boreoeutheria, Euarchontoglires, Glires, Rodentia, Sciurognathi, Murinae and Mus musculus)
  # - dubious duplciations are allowed because with this criterion we find almost no tree matching (because a lot of mammal species are present in Ensembl, so some of them almost always group to the wrong place in teh mammal tree)
  # - only 1 mouse gene
  # - there is a fish outgroup with at least a few species (at least one node with taxonomic level: Clupeocephala, Acanthomorphata, Percomorphaceae, Tetraodontidae, Atherinomorphae, Poeciliinae or Otophysa)

  my $duplication = 0;
  my $outgroup = 0;
  my $mouse_gene;
  foreach my $child_node (@{$root->get_all_nodes()}) {
    if ($child_node->is_leaf()) {
      # retrieve mouse gene
      if ($child_node->gene_member->stable_id =~m/ENSMUSG/) {
        $mouse_gene = $child_node->gene_member->stable_id;
      }
    }
    else {
      # test that there was no duplication or gene split in branches leading to mouse
      if ((($child_node->get_tagvalue('node_type') eq 'duplication') or ($child_node->get_tagvalue('node_type') eq 'gene_split')) and (($child_node->taxonomy_level() eq 'Vertebrata') or ($child_node->taxonomy_level() eq 'Euteleostomi') or ($child_node->taxonomy_level() eq 'Sarcopterygii') or ($child_node->taxonomy_level() eq 'Tetrapoda') or ($child_node->taxonomy_level() eq 'Amniota') or ($child_node->taxonomy_level() eq 'Mammalia') or ($child_node->taxonomy_level() eq 'Theria') or ($child_node->taxonomy_level() eq 'Eutheria') or ($child_node->taxonomy_level() eq 'Boreoeutheria') or ($child_node->taxonomy_level() eq 'Euarchontoglires') or ($child_node->taxonomy_level() eq 'Glires') or ($child_node->taxonomy_level() eq 'Rodentia') or ($child_node->taxonomy_level() eq 'Sciurognathi') or ($child_node->taxonomy_level() eq 'Murinae') or ($child_node->taxonomy_level() eq 'Mus musculus'))) {
        $duplication = 1;
      }

      # test the fish subtree. We want several genes, but they may be higher than Clupeocephala (Acanthomorphata, Percomorphaceae...)
      if ((($child_node->taxonomy_level() eq 'Clupeocephala') or ($child_node->taxonomy_level() eq 'Acanthomorphata') or ($child_node->taxonomy_level() eq 'Percomorphaceae') or ($child_node->taxonomy_level() eq 'Tetraodontidae') or ($child_node->taxonomy_level() eq 'Atherinomorphae') or ($child_node->taxonomy_level() eq 'Poeciliinae') or ($child_node->taxonomy_level() eq 'Otophysa')) and ($child_node->get_tagvalue('node_type') eq 'speciation')) {
        $outgroup = 1;
      }
    }
  }

  if (($duplication eq 0) and (defined $mouse_gene) and ($outgroup eq 1)) {
    print "\tThis subtree was selected ($filename)\n";
    open(OUT, '>>', $filename) or die "Can't open file";
    print OUT $mouse_gene, "\t", $root->tree->stable_id, "\n";
    close OUT;
  }
  return;
}

#####################################################################
