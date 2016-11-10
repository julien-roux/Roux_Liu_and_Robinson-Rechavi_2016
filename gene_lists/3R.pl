#!/usr/bin/env perl -w
use strict;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Registry;
$|=1; # no buffering of the output

# Julien Roux and Jialin Liu
# This script looks for Ensembl gene trees with topology compatible or not with the 3R teleost-specific whole-genome duplication. It was developed for Ensembl release 79. I might not work for more recent releases

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

# To test the taxonomical level of the gene trees, we need first to define the deepest taxonomical level we want to keep
# Neopterygii: tax_id = 41665
my $limit_node = $ncbi_taxon_adaptor->fetch_node_by_taxon_id('41665');
my $limit_left = $limit_node->left_index;
my $limit_right = $limit_node->right_index;

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

# recursive subroutine ############################################
sub select_sub_tree {
  my ($this_tree, $level) = @_;
  #  print "Tree root:", $this_tree, "\n";
  print "Tree root ID: ", $this_tree->node_id, ", level: ", $level, "\n";

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
        $this_tree->print_tree();
        #test_3R_dup($this_tree);
        #test_3R_sing($this_tree);
        test_3R_dup_orth($this_tree);
        test_3R_sing_orth($this_tree);
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

# test_3R_dup ####################################################################
sub test_3R_dup {
  my ($root) = @_;

  # 3R topology: we want a neopterygii speciation node followed by a clupeocephala duplication node (high confidence) and a spotted gar ortholog
  my %genes;
  foreach my $child1 (@{$root->children()}) {
    # check that after the root there is a duplication node "clupeocephala", followed by 2 speciation nodes "clupeocephala" (that is including at least one danio gene and another fish gene)
    print "child1:", $child1 ,"\n";
    $child1->print_tree();
    # test for a 3R duplication just after the root
    if (($child1->taxonomy_level() eq 'Clupeocephala') and ($child1->get_tagvalue('node_type') eq 'duplication') and ($child1->get_tagvalue("duplication_confidence_score") ge 0.5)) {

      # check that there are 2 speciation nodes after (both including a danio gene and another fish gene)
      my $count_speciations = 0;
      foreach my $child2 (@{$child1->children()}) {
        if (($child2->taxonomy_level() eq 'Clupeocephala') and ($child2->get_tagvalue('node_type') eq 'speciation')) {
          $count_speciations++;
          foreach my $leaf (@{$child2->get_all_leaves}) {
            my $gene = $leaf->gene_member->stable_id;
            if ($gene =~m/ENSDARG/) {
              print "\t", $gene ,"\n";
              push(@{$genes{$count_speciations}}, $gene);
            }
          }
        }
      }
    }
  }

  # if both speciation tests were ok
  open(OUT, '>>zebrafish_3R_ohnologs.txt') or die "Can't open file";
  if (scalar(keys %genes) eq 2) {
    if ((scalar(@{$genes{'1'}}) eq 1) and (scalar(@{$genes{'2'}}) eq 1)) {
      print "\tThis tree was selected (3R duplication)\n";
      print OUT ${$genes{'1'}}[0], "\t", $root->tree->stable_id, "\n";
      print OUT ${$genes{'2'}}[0], "\t", $root->tree->stable_id, "\n";
    }
  }
  close OUT;
  return;
}

# test_3R_sing ############################################################
sub test_3R_sing {
  my ($root) = @_;

  foreach my $child1 (@{$root->children()}) {
    # We want to retrieve 3R singletons in danio
    #  - no duplication allowed in the tree after 3R on branches leading to zebrafish (taxonomic levels Otophysa or zebrafish)
    #  - only 1 zebrafish gene

    # test for speciation node just after the Neopterygii root. If only zebrafish and cave fish genes only present, the speciation node will be Otophysa
    if (($child1->taxonomy_level() eq 'Clupeocephala') and ($child1->get_tagvalue('node_type') eq 'speciation')) {
      my $duplication = 0;
      my $zf_gene;
      foreach my $child_node (@{$child1->get_all_nodes()}) {
        if ($child_node->is_leaf()){
          if ($child_node->gene_member->stable_id =~m/ENSDARG/) {
            $zf_gene = $child_node->gene_member->stable_id;
          }
        }
        elsif (($child_node->get_tagvalue('node_type') eq 'duplication') and ($child_node->taxonomy_level() eq 'Otophysa')){
          $duplication = 1;
        }
        elsif (($child_node->get_tagvalue('node_type') eq 'duplication') and ($child_node->taxonomy_level() eq 'Danio rerio')){
          $duplication = 1;
        }
      }

      if (($duplication eq 0) and (defined $zf_gene)){
        print "\tThis tree was selected (3R singleton)\n";
        open(OUT, '>>zebrafish_3R_singletons.txt') or die "Can't open file";
        print OUT $zf_gene, "\t", $root->tree->stable_id, "\n";
        close OUT;
      }
    }
  }
  return;
}

# test_3R_dup_orth ############################################################
sub test_3R_dup_orth {
  my ($root) = @_;
  # We want two zebrafish 3R duplicates (see test_3R_dup), plus a single ortholog in mouse / human

  # 3R topology: we want neopterygii speciation node and then a clupeocephala duplication node (high confidence) and a spotted gar ortholog
  my %zf_genes;
  foreach my $child1 (@{$root->children()}) {
    # check that after the root there is a duplication node "clupeocephala", followed by 2 speciation nodes "clupeocephala" (that is including at least one danio gene and another fish gene)

    # test for a 3R duplication just after the root
    if (($child1->taxonomy_level() eq 'Clupeocephala') and ($child1->get_tagvalue('node_type') eq 'duplication') and ($child1->get_tagvalue("duplication_confidence_score") ge 0.5)) {

      # check that there are 2 speciation nodes after (both including a danio gene and another fish gene)
      my $count_speciations = 0;
      foreach my $child2 (@{$child1->children()}) {
        if (($child2->taxonomy_level() eq 'Clupeocephala') and ($child2->get_tagvalue('node_type') eq 'speciation')) {
          $count_speciations++;
          foreach my $leaf (@{$child2->get_all_leaves}) {
            my $gene = $leaf->gene_member;
            if ($gene->stable_id =~m/ENSDARG/) {
              push(@{$zf_genes{$count_speciations}}, $gene);
            }
          }
        }
      }
    }
  }

  # if both speciation tests were ok, and there is only 1 zebrafish gene in each clupeocephala group
  if (scalar(keys %zf_genes) eq 2) {
    if ((scalar(@{$zf_genes{'1'}}) eq 1) and (scalar(@{$zf_genes{'2'}}) eq 1)) {
      # check mouse orthologs: should have only 1
      my @homologies_mouse;
      foreach my $homology (@{$homology_adaptor->fetch_all_by_Member_paired_species(${$zf_genes{'1'}}[0], "Mus_musculus")}) {
        if (($homology->taxonomy_level() eq 'Euteleostomi') and ($homology->description eq 'ortholog_one2many') and ($homology->is_tree_compliant() eq 1)) {
          push(@homologies_mouse, $homology);
        }
      }
      if (scalar(@homologies_mouse eq 1)){
        print "\tThis tree was selected (mouse 3R dup ortholog)\n";
        foreach my $member (@{$homologies_mouse[0]->gene_list}) {
          if ($member->stable_id =~m/ENSMUSG/) {
            open(OUT1, '>>mouse_orthologs_3R_ohnologs.txt') or die "Can't open file";
            print OUT1 $member->stable_id, "\t", ${$zf_genes{'1'}}[0]->stable_id, "\t", ${$zf_genes{'2'}}[0]->stable_id, "\t", $root->tree->stable_id, "\n";
            close OUT1;
          }
        }
      }

      # check human orthologs: should have only 1
      my @homologies_human;
      foreach my $homology (@{$homology_adaptor->fetch_all_by_Member_paired_species(${$zf_genes{'1'}}[0], "Homo_sapiens")}) {
        if (($homology->taxonomy_level() eq 'Euteleostomi') and ($homology->description eq 'ortholog_one2many') and ($homology->is_tree_compliant() eq 1)) {
          push(@homologies_human, $homology);
        }
      }
      if (scalar(@homologies_human eq 1)){
        print "\tThis tree was selected (human 3R dup ortholog)\n";
        foreach my $member (@{$homologies_human[0]->gene_list}) {
          if ($member->stable_id =~m/ENSG/) {
            open(OUT2, '>>human_orthologs_3R_ohnologs.txt') or die "Can't open file";
            print OUT2 $member->stable_id, "\t", ${$zf_genes{'1'}}[0]->stable_id, "\t", ${$zf_genes{'2'}}[0]->stable_id, "\t", $root->tree->stable_id, "\n";
            close OUT2;
          }
        }
      }
    }
  }
  return;
}

# test_3R_sing_orth ################################################
sub test_3R_sing_orth {
  my ($root) = @_;
  my $zf_gene;
  foreach my $child1 (@{$root->children()}) {
    # We want to retrieve 3R singletons in zebrafish
    #  - no duplication in the tree after 3R on branches leading to zebrafish (tax levels Otophysa or zebrafish)
    #  - only 1 zebrafish gene
    if (($child1->taxonomy_level() eq 'Clupeocephala') and ($child1->get_tagvalue('node_type') eq 'speciation')) {
      my $duplication = 0;
      foreach my $child_node (@{$child1->get_all_nodes()}) {
        if ($child_node->is_leaf()) {
          if ($child_node->gene_member->stable_id =~m/ENSDARG/) {
            $zf_gene = $child_node->gene_member;
          }
        }
        elsif (($child_node->get_tagvalue('node_type') eq 'duplication') and ($child_node->taxonomy_level() eq 'Otophysa')) {
          $duplication = 1;
        }
        elsif (($child_node->get_tagvalue('node_type') eq 'duplication') and ($child_node->taxonomy_level() eq 'Danio rerio')) {
          $duplication = 1;
        }
      }

      # if there is no 3R duplication, check that there is a single mouse / human ortholog
      if (($duplication eq 0) and (defined $zf_gene)) {
        # retrieve the homologies between the zebrafish gene and mouse genes
        my @homologies_mouse;
        foreach my $homology (@{$homology_adaptor->fetch_all_by_Member_paired_species($zf_gene, "Mus_musculus")}) {
          print $homology->description, "\n";
          if (($homology->taxonomy_level() eq 'Euteleostomi') and ($homology->description eq 'ortholog_one2one') and ($homology->is_tree_compliant() eq 1)) {
            push(@homologies_mouse, $homology);
          }
        }
        if (scalar(@homologies_mouse eq 1)) {
          print "\tThis tree was selected (mouse 3R sing ortholog)\n";
          foreach my $member (@{$homologies_mouse[0]->gene_list}) {
            if ($member->stable_id =~m/ENSMUSG/) {
              open(OUT1, '>>mouse_orthologs_3R_singletons.txt') or die "Can't open file";
              print OUT1 $member->stable_id, "\t", $zf_gene->stable_id, "\t", $root->tree->stable_id, "\n";
              close OUT1;
            }
          }
        }

        # retrieve the homologies between the zebrafish gene and human genes
        my @homologies_human;
        foreach my $homology (@{$homology_adaptor->fetch_all_by_Member_paired_species($zf_gene, "Homo_sapiens")}) {
          print $homology->description, "\n";
          if (($homology->taxonomy_level() eq 'Euteleostomi') and ($homology->description eq 'ortholog_one2one') and ($homology->is_tree_compliant() eq 1)) {
            push(@homologies_human, $homology);
          }
        }
        if (scalar(@homologies_human eq 1)) {
          print "\tThis tree was selected (human 3R sing ortholog)\n";
          foreach my $member (@{$homologies_human[0]->gene_list}) {
            if ($member->stable_id =~m/ENSG/) {
              open(OUT2, '>>human_orthologs_3R_singletons.txt') or die "Can't open file";
              print OUT2 $member->stable_id, "\t", $zf_gene->stable_id, "\t", $root->tree->stable_id, "\n";
              close OUT2;
            }
          }
        }
      }
    }
  }
  return;
}
####################################################################





