#!/usr/bin/perl
use warnings;
use strict;
# Julien Roux and Jialin Liu
# July 2016
# This script extracts complex information from UniProt free text anotation (SUBUNIT field), maps UniProt IDs to Ensembl gene IDs, and classifies genes into members of different categories of complexes

# Usage: perl retrieve_monomer.pl mouse_EnsemblID_2_UniProtGeneName.txt mouse_complexes_UniProt.txt
my $infile1 = shift @ARGV;
my $infile2 = shift @ARGV;

# mapping UniProt to Ensembl
my %geneId;
open (IN1,"<", $infile1) || die "cannot open file: $!\n";
while (<IN1>) {
  chomp;
  my @array=split(/\t/, $_);
  # construct a hash: uniprotGeneName as key, ensemblID as value.
  if ((defined $array[1]) and ($array[1] ne "")){
    $geneId{$array[1]}->{$array[0]}=();
  }
}
# Keep uniprot genes mapped to a unique ensembl ID
foreach my $uniprotId (keys %geneId){
  if (length(scalar keys %{$geneId{$uniprotId}}) > 1){
    delete $geneId{$uniprotId};
  }
}
close IN1;

# Open output files
open(OUT2, ">", "mouse_heterodimers.txt") or die "cannot open file: $!\n";
open(OUT1, ">", "mouse_heteromultimers.txt") or die "cannot open file: $!\n";
open(OUT3, ">", "mouse_homomultimers.txt") or die "cannot open file: $!\n";
open(OUT4, ">", "mouse_monomers.txt") or die "cannot open file: $!\n";
open(OUT5, ">", "mouse_undefined_complexes.txt") or die "cannot open file: $!\n";

# Parse UniProt information
open (IN2, "<", $infile2) || die "cannot open file: $!\n";
while (<IN2>) {
  chomp;
  my @array = split(/\t/, $_);
  next unless (defined $array[2]);

  # set flag to 1 as soon as the gene is classified in one category
  my $flag = 0;

  ## Monomer ##############################
  if (($array[2] =~ /SUBUNIT:\sMonomer.*/) or
      ($array[2] =~ /\.\sMonomer.*/) or
      ($array[2] =~ /SUBUNIT:\sPredominantly\smonomer.*/) or
      ($array[2] =~ /and\smonomer\s.*/) or
      ($array[2] =~ /Exists(\seither)?(\sboth)?\sas(\sboth)?(\sa)?\smonomer/) or
      ($array[2] =~ /(Seems to also|A|but\sa)lso (exist|exists|detected|secreted|found) as monomer/) or
      ($array[2] =~ /Binds(\ssingle-stranded telomeric)?\sDNA as(\sa)? monomer/) or
      ($array[2] =~ /SUBUNIT: Exists in a dynamic equilibrium between monomeric/) or
      ($array[2] =~ /SUBUNIT: Forms(\sa)?\smonomer/) or
      ($array[2] =~ /Mostly monomeric/) or
      ($array[2] =~ /Active as monomer/) or
      ($array[2] =~ /Activation by cAMP releases the two active catalytic monomers/) or
      ($array[2] =~ /SUBUNIT: Present as a mixture of monomer/) or
      ($array[2] =~ /Can form homooligomers \(monomers/) or
      ($array[2] =~ /however it is likely to function as a monomer./) or
      ($array[2] =~ /SUBUNIT: Proposed to be monomeric/)
  ){
    $flag = 1;
    foreach my $ensemblId (keys %{$geneId{$array[0]}}) {
      print OUT4 $ensemblId, "\t", $array[0], "\t", $array[2],"\n";
    }
  }

  ## Homo-multimers ######################
  ## Matches to homo(?!typ)(?!log)(?!phil).* to ensure no false positives (homologous, homotypic, homolphilic, etc)
  if (($array[2] =~ /SUBUNIT:\s(H|h)omo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT:\s.+\.\sHomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT:\s[\w\-\,]+\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /and(\sa)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /or(\sa)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT:\s[\w-]+ or homomer/) or
      ($array[2] =~ /SUBUNIT:\sCan\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /(Forms|Principally)(\seither)?(\sa)?(\sneddylation-dependent)?(\slinear)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /Can(\salso)?\s(form|occur)(\sas)?(\sa)?(\sfunctional)?(\slinear)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /(May|Might)(\soligomerize\sand)?(\salso)?(\sform)?(\sbe)?(\sa)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /Exists\sas\sa(\smixture\sof)?(n octamer composed of four \w+)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /Exists\sboth\sas\sa(\scovalent|\smonomer\sand\sa)?(\sdisulfide-linked)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /Self-associates.*/) or
      ($array[2] =~ /(c|C)an self-associate (to|and) form homo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /(Can b|B)ind(s)?(\sto)?\sDNA(\sspecifically)? as(\sa)? homo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT: Antiparallel homodimer/) or
      ($array[2] =~ /SUBUNIT: Binds DNA as(\sa)?(\sdimer;)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT: Binds DNA as a dimer and can form stable heterodimers/) or
      ($array[2] =~ /SUBUNIT: Binds DNA as a dimer and can form a homodimer/) or
      ($array[2] =~ /SUBUNIT: Binds DNA as a dimer. Probably homodimerizes/) or
      ($array[2] =~ /SUBUNIT: Predominantly homodimeric/) or
      ($array[2] =~ /SUBUNIT: Monomer, homo- or heterodimer/) or
      ($array[2] =~ /SUBUNIT: Monomer, and domain-swapped homodimer/) or
      ($array[2] =~ /SUBUNIT: The active form is probably a homotrimer/) or
      ($array[2] =~ /leads to homotrimerization/) or
      ($array[2] =~ /SUBUNIT: Active as a homodimer or oligomer/) or
      ($array[2] =~ /SUBUNIT: Membrane region forms a homotetramer/) or
      ($array[2] =~ /SUBUNIT: Cis- and trans-homodimer/) or
      ($array[2] =~ /SUBUNIT: Head to tail homodimer/) or
      ($array[2] =~ /SUBUNIT: Forms both disulfide-linked homodimers and higher homooligomers/) or
      ($array[2] =~ /SUBUNIT:\sOligomeric\scomplex\sof\s\d+(\sset(s)?\sof)?(\sor\smore)?\shomo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT: Oligomer of(\sdisulfide-linked)?homo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT: May interact with [\w\d\-] and form homo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT: Double stacked ring-shaped homo(?!typ)(?!log)(?!phil).*/) or
      ($array[2] =~ /SUBUNIT: Hexamer formed by 3 homodimer/) or
      ($array[2] =~ /SUBUNIT: Dodecameric ring assembled from homodimer/) or
      ($array[2] =~ /SUBUNIT: Functional P2XRs are organized as homo(?!typ)(?!log)(?!phil).*/)
    ) {
    $flag = 1;
    foreach my $ensemblId (keys %{$geneId{$array[0]}}) {
      print OUT3 $ensemblId, "\t", $array[0], "\t", $array[2],"\n";
    }
  }
  ## Hetero-dimers ########################
  if (($array[2] =~ /SUBUNIT:\sHeterodimer.*/) or 
      ($array[2] =~ /\.\sHeterodimer.*/) or
      ($array[2] =~ /SUBUNIT:\s[\w\-\,]+\sheterodimer/) or
      ($array[2] =~ /SUBUNIT:\s[\w\d-]+\sis\sa\sheterodimer/) or
      ($array[2] =~ /SUBUNIT:\s[\w\d\s\(\)\-]+\s(form(s)?|is|proposed\sto\sbe) a heterodimer/) or
      ($array[2] =~ /Forms(\sa)?(\sstable)?(\shigh-affinity)?\sheterodimer.*/) or
      ($array[2] =~ /May(\salso)?\sform(\sa)?\sheterodimer.*/) or
      ($array[2] =~ /(?!CSC\s)and(\sa)?\sheterodimer.*/) or
      ($array[2] =~ /or(\sa)?\sheterodimer.*/) or
      ($array[2] =~ /Can\sheterodimerize.*/) or
      ($array[2] =~ /(?!subunit\s)(c|C)(an|ould)(\salso)?\sform(\sa)?(\sdisulfide\-linked)?\s(trans\-)?heterodimer.*/) or
      ($array[2] =~ /Exists also as heterodimer.*/) or
      ($array[2] =~ /Also found as heterodimer.*/) or
      ($array[2] =~ /Also forms(\sa)?\sheterodimer.*/) or
      ($array[2] =~ /SUBUNIT: Forms protease inhibiting heterodimer/) or 
      ($array[2] =~ /SUBUNIT: (In\sgeneral|Mostly)\sheterodimer/) or 
      ($array[2] =~ /SUBUNIT:(\s(A|a)ntiparallel)?(\sDisulfide-linked)? heterodimer/) or 
      ($array[2] =~ /\/heterodimer/) or
      ($array[2] =~ /SUBUNIT: PP2A consists of a common heterodimeric/) or 
      ($array[2] =~ /SUBUNIT: TFIIA is a heterodimer/) or
      ($array[2] =~ /SUBUNIT: Dimer of alpha and beta chains/) or
      ($array[2] =~ /SUBUNIT: Predominantly (forms|non\-covalently\slinked|detected as) heterodimer/) or
      ($array[2] =~ /Component of the heterodimeric/) or
      ($array[2] =~ /SUBUNIT: Probably heterodimerizes/) or
      ($array[2] =~ /(Can\sb|B)ind(s)?(\sto)?\sDNA(\sas\sa\shomodimer\sand)? as a heterodimer/) or
      ($array[2] =~ /SUBUNIT: Can occur as a homodimer or as a heterodimer/) or
      ($array[2] =~ /SUBUNIT: Binds DNA as a dimer and can form stable heterodimers/) or
      ($array[2] =~ /SUBUNIT: The interleukin\-\d+ receptor complex is a heterodimer/) or
      ($array[2] =~ /Efficient DNA binding requires dimerization with another/) or
      ($array[2] =~ /Forms hetero(t)?multimers with AHCYL/)
   ){
    $flag = 1;
    foreach my $ensemblId (keys %{$geneId{$array[0]}}) {
      print OUT2 $ensemblId, "\t", $array[0], "\t", $array[2],"\n";
    }
  }

  ## Hetero-multimers ########################
  ## Other instances of hetero*: heterotrimer, heterotetramer, etc
  ## Excluding: heterodimer
  ## Excluding: heteromer, heteromultimer, heterooligomer, heteropolymer -> other complexes
  ## Excluding: heterogeneous, heterophilic
  ## -> Regex: hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*
  if (($array[2] =~ /SUBUNIT:\sHetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or 
      ($array[2] =~ /SUBUNIT:\s[\w\-\,]+\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT:\s[\w\d-]+\sis\sa\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT:\s[\w\d\s\(\)\-]+\s(form(s)?|is|proposed\sto\sbe) a hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Part of(\sa)?(\sthe)? hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT: Subunit of(\sa)?(\sthe)? hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /\.\sHetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Forms(\sa)?(\sstable)?(\shigh-affinity)?(\s1:1:1)?\shetero(?!dimer)(?!mer)(?!multimer)(?!(-)?oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Forms(\sa)?(\sstable)?(\shigh-affinity)?(\s1:1:1)?\s(\w+)hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Forms\s\w+(and|or)\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /May(\salso)?\s(form|be\spart\sof)(\sa)?\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /and(\sprobably)?(\sa)?\shetero(?!(-)?dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /or(\sa)?\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Can\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /(c|C)(an|ould)(\salso)?\sform(\sa)?(\sdisulfide\-linked)?\s(trans\-)?hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT:\sCan form\s[\w\d\-]+\s(and|or)(\sa)? hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /(Preferentially\se|Exists)(\salso)? as(\sa)?(\s\w+\sor)? hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Also found as hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Also forms(\sa)?\shetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT: In general hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or 
      ($array[2] =~ /\/hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT: Predominantly (forms|non\-covalently\slinked|detected as) hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /(Component\sof|Identified\sin)(\sthe)?(\sa)? hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT: Probably hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /(Can\sb|B)ind(s)?(\sto)?\sDNA(\sas\sa\shomodimer\sand)? as a hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT: Binds DNA as a dimer and can form stable hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /Active as a hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/) or
      ($array[2] =~ /SUBUNIT: Hexadecamer of 4 hetero(?!dimer)(?!mer)(?!multimer)(?!oligomer)(?!polymer)(?!geneous)(?!philic).*/)
   ){
    $flag = 1;
    foreach my $ensemblId (keys %{$geneId{$array[0]}}) {
      print OUT1 $ensemblId, "\t", $array[0], "\t", $array[2],"\n";
    }
  }

  # Undefined complexes ########################
  if (($flag eq 0) and ## not already classified elsewhere
     (($array[2] =~ /.*(F|f)ound in\s[\w\d\s\(\)\-]+\scomplex\s.*/) or
      ($array[2] =~ /.*(P|p)art of\s[\w\d\s\(\)\-]+\scomplex\s.*/) or 
      ($array[2] =~ /.*(C|c)omponent of\s[\w\d\s\(\)\-]+\scomplex\s.*/) or
      ($array[2] =~ /.*(F|f)orms\s[\w\d\s\(\)\-]+\scomplex\s.*/) or

      ## Matches to hetero(mer|multimer)|oligomer|polymer).*
      ($array[2] =~ /SUBUNIT:\sHetero(mer|multimer|oligomer|polymer).*/) or 
      ($array[2] =~ /SUBUNIT:\s[\w\-\,]+\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT:\s[\w\d-]+\sis\sa\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT:\s[\w\d\s\(\)\-]+\s(form(s)?|is|proposed\sto\sbe) a hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Part of(\sa)?(\sthe)? hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT: Subunit of(\sa)?(\sthe)? hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /\.\sHetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Forms(\sa)?(\sstable)?(\shigh-affinity)?(\s1:1:1)?\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Forms\s\w+(and|or)\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /May(\salso)?\s(form|be\spart\sof)(\sa)?\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /and(\sprobably)?(\sa)?\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /or(\sa)?\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Can\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /(c|C)(an|ould)(\salso)?\sform(\sa)?(\sdisulfide\-linked)?\s(trans\-)?hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT:\sCan form\s[\w\d\-]+\s(and|or)(\sa)? hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /(Preferentially\se|Exists)(\salso)? as(\sa)?(\s\w+\sor)? hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Also found as hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Also forms(\sa)?\shetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT: In general hetero(mer|multimer|oligomer|polymer).*/) or 
      ($array[2] =~ /\/hetero(mer|multimer|(-)?oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT: Predominantly (forms|non\-covalently\slinked|detected as) hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /(Component\sof|Identified\sin)(\sthe)?(\sa)? hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /SUBUNIT: Probably hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /(Can\sb|B)ind(s)?(\sto)?\sDNA(\sas\sa\shomodimer\sand)? as a hetero(mer|multimer|oligomer|polymer).*/) or
      ($array[2] =~ /Active as a hetero(mer|multimer|oligomer|polymer).*/)
     )){
    foreach my $ensemblId (keys %{$geneId{$array[0]}}) {
      print OUT5 $ensemblId, "\t", $array[0], "\t", $array[2],"\n";
    }
  }
}
close IN2;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;



