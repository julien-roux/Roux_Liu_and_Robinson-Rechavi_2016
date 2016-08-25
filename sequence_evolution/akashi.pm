# Julien Roux
# This script launches Akashi's test on one pair of mouse-rat 1-to-1 orthologs
# This script was developed for Ensembl release 79. I might not work for more recent releases

sub akashi_test {
  my ($orthology) = @_;

  # all AA except W and M (only one codon)
  my @amino_acids = ('G', 'P', 'A', 'V', 'L', 'I', 'C', 'F', 'Y', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T');

  # 1 for prefered codons in mouse
  # taken from drummond 2008, cell (waterston 2002)
  my %codons = (
    'GCC' => 1,
    'GCA' => 0,
    'GCG' => 0,
    'TGC' => 1,
    'GAC' => 1,
    'GAA' => 1,
    'GAG' => 1,
    'TTC' => 1,
    'GGC' => 1,
    'GGA' => 0,
    'GGG' => 0,
    'CAC' => 1,
    'ATC' => 1,
    'ATA' => 0,
    'AAA' => 1,
    'AAG' => 1,
    'CTC' => 0,
    'CTG' => 1,
    'TTG' => 0,
    'CTA' => 0,
    'TTA' => 0,
    'AAC' => 1,
    'CCC' => 1,
    'CCA' => 0,
    'CCG' => 0,
    'CAA' => 0,
    'CAG' => 1,
    'CGC' => 1,
    'CGA' => 0,
    'AGG' => 0,
    'AGA' => 0,
    'CGG' => 0,
    'TCC' => 1,
    'AGC' => 1,
    'TCG' => 0,
    'TCA' => 0,
    'ACC' => 1,
    'ACA' => 0,
    'ACG' => 0,
    'GTA' => 0,
    'GTG' => 1,
    'GTC' => 0,
    'TAC' => 1,
    # not appearing in the table (because of inosine)
    'TTT' => 0,
    'CTT' => 0,
    'ATT' => 0,
    'GTT' => 0,
    'TCT' => 0,
    'CCT' => 0,
    'ACT' => 0,
    'GCT' => 0,
    'TAT' => 0,
    'TAA' => 0,
    'TAG' => 0,
    'CAT' => 0,
    'AAT' => 0,
    'GAT' => 0,
    'TGT' => 0,
    'TGA' => 0,
    'CGT' => 0,
    'AGT' => 0,
    'GGT' => 0,
    #   'NNN' => 0,
    #   '---' => 0
  );

  my $simple_align_prot = $orthology->get_SimpleAlign();
  my $simple_align_dna = $orthology->get_SimpleAlign(-SEQ_TYPE => 'cds');
  # print "simple_align_prot:", $simple_align_prot, "\n";

  #initialize the score values
  my $sum_a = 0;
  my $sum_e_a = 0;
  my $sum_v_a = 0;
  my $sum_ad_n = 0; #ad/n
  my $sum_bc_n = 0; #bc/n

  # foreach amino acid (except W and M)
  foreach my $aa (@amino_acids) {
    # Laplace estimation (adding 1 to all counts in order to remove empty cells problem)
    my %table = ('a' => 1,
                 'b' => 1,
                 'c' => 1,
                 'd' => 1,
                );

    for (my $i=1; $i<=$simple_align_prot->length; $i++) {
      # for each occurence in the protein sequence
      my $seq1 = $simple_align_prot->get_seq_by_pos(1);
      my $seq2 = $simple_align_prot->get_seq_by_pos(2);

      my $seq_dna = $simple_align_dna->get_seq_by_pos(1);
      if ($seq1->subseq($i,$i) eq $aa) {
        # look at the conservation in the other species
        my $conservation = 0;
        $conservation = 1 if $seq2->subseq($i,$i) eq $aa;

        # look at mouse DNA sequence to see if it is a prefered codon
        my $preferred = $codons{$seq_dna->subseq(3*$i-2, 3*$i)};
        # if the conservation is perfect
        if ($conservation eq '1') {
          if ($preferred eq '1') {
            $table{'a'}++;
          }
          else {
            $table{'c'}++;
          }
        }
        else {
          if ($preferred eq '1') {
            $table{'b'}++;
          }
          else {
            $table{'d'}++;
          }
        }
      }
    }
    #    print "$aa\n";
    #    print $table{'a'}, "\t", $table{'b'}, "\n";
    #    print $table{'c'}, "\t", $table{'d'}, "\n";

    if (($table{'b'} ne 0) and ($table{'c'} ne 0) and ($table{'a'}+$table{'b'} ne 0) and  ($table{'a'}+$table{'c'} ne 0) and ($table{'b'}+$table{'d'} ne 0) and ($table{'c'}+$table{'d'} ne 0)) {
      my $n = $table{'a'}+$table{'b'}+$table{'c'}+$table{'d'};
      my $e_a = ($table{'a'}+$table{'b'})*($table{'a'}+$table{'c'})/$n;
      my $v_a = ($table{'a'}+$table{'b'})*($table{'a'}+$table{'c'})*($table{'b'}+$table{'d'})*($table{'c'}+$table{'d'})/($n*$n*($n-1));
      my $z = ($table{'a'} - $e_a)/sqrt($v_a);
      my $psi = ($table{'a'}*$table{'d'})/($table{'b'}*$table{'c'});
      # print "Z: $z\nPsi: $psi\n";

      $sum_a += $table{'a'};
      $sum_e_a += $e_a;
      $sum_v_a += $v_a;
      $sum_ad_n += ($table{'a'}*$table{'d'})/$n;
      $sum_bc_n += ($table{'b'}*$table{'c'})/$n;
    }
    else {
       print "Statistics cannot be calculated for this AA\n";
    }
  }
  return($sum_a, $sum_e_a, $sum_v_a, $sum_ad_n, $sum_bc_n)
}
1;
