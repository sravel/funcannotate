#!/usr/bin/perl  
use Getopt::Long;
use FindBin;
use Pod::Usage; 
use Bio::SeqIO;
use Bio::Tools::GFF;
  
my $help      = "";
my $file      = ""; 
my $gff3_file = "";
my $go_file   = "";
my $interpro_file = "";
my $product = "";
my $te_file = "";
my $tag    = "Name";
my $rename;
my %ipr;
my %go;
my %locus_info;
my $prefix   = "musa_acuminata_v2";
my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters

    -gff3_file      gff3 file to be annotated 
    -te_file        remove TE
    -prefix         prefix name for fasta output
    -product        locus tag file, required
    -go_file
    -interpro_file
    -tag            Default Name
    -rename 
    -help
/;
GetOptions(  
    'gff3_file=s'     => \$gff3_file,  
    'te_file=s'       => \$te_file,
    'go_file=s'       => \$go_file,
    'interpro_file=s' => \$interpro_file,
    'product=s'       => \$product,
    'prefix=s'        => \$prefix,
    'tag=s'           => \$tag,
    'rename'          => \$rename,
    'help|h|?'        => \$help
)  or pod2usage(-message => $message);
if ($help) { pod2usage(-message => $message); } 

if ($gff3_file eq "") {
    warn "\nWarn :: --file is empty. Please specify the locus tag file\n\n";
    warn $usage;
    exit 0;
}      
my $gff3_out_file = $prefix ."_mrna.gff3";
my $gff3_ncrna_out_file = $prefix ."_ncrna.gff3"; 
my $gff_in = new Bio::Tools::GFF(
    -file => $gff3_file,
    -gff_version => 3
);
my $gff_out = new Bio::Tools::GFF(
    -file => ">$gff3_out_file",
    -gff_version => 3
);
my $gff_ncrna_out = new Bio::Tools::GFF(
    -file => ">$gff3_ncrna_out_file",
    -gff_version => 3
); 
my %to_remove;
open(IN,$te_file);
while(<IN>){
    chomp;
    $to_remove{$_}=1;
}
close IN;

open(IN,$product); 
while (<IN>) {
    chomp; 
    my ($primary_id,$hit_name, $source, $length, $product,$ic_code,$completeness,$gene_name,$dbxref) = (split(/\t/,$_)) ;  
    my $gene_id = (split(/\./,$primary_id))[0];

    my $mrna_id = $primary_id;
   # $primary_id =~ s/g/t/;
   #$primary_id = $primary_id .".1";
    $completeness = "missing_functional_completeness" if $completeness eq "";
    my $first_gene_name = (split(/\,/,$gene_name))[0];
    my $note = join("~ ",$gene_id,$product,$first_gene_name,$completeness);
    
    $locus_info{$gene_id} = { 
        product    => $product
    };
    $locus_info{$mrna_id} = { 
        product    => $product
    }; 
    
    $locus_info{$primary_id} = { 
        product      => $product,
        completeness => $completeness,
        ic_code      => $ic_code,
        gene_name    => $gene_name,
        note         => $note,
        dbxref       => $dbxref
    };    
} 
close IN; 
if ($go_file) {
    open(GO,$go_file);
    while (<GO>) {
        chomp;
        my ($primary_id,$go_id) = (split(/\t/,$_)); 
        push @{$go{$primary_id}} , $go_id;
    }
    close GO;
}
if ($interpro_file) {
    open(IPR,$interpro_file);
    while (<IPR>) {
        chomp;
        my ($primary_id,$ipr_id) = (split(/\t/,$_)); 
        next if $ipr_id eq "-"; 
        push @{$ipr{$primary_id}} , $ipr_id;
    }
    close IPR;
}
#Contig85        EuGene  gene    10248   13220   .       +       .       ID=gene:Contig85g0000001;Name=Contig85g0000001;length=2973
#Contig85        EuGene  mRNA    10248   13220   .       +       .       ID=mRNA:Contig85g0000001;Name=Contig85g0000001;Parent=gene:Contig85g0000001;nb_exon=6;length=1379
#Contig85        EuGene  exon    10248   10402   .       +       .       ID=exon:Contig85g0000001.1;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000200
#Contig85        EuGene  exon    11031   11276   .       +       2       ID=exon:Contig85g0000001.2;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004
#Contig85        EuGene  exon    11423   11625   .       +       2       ID=exon:Contig85g0000001.3;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004
#Contig85        EuGene  exon    11697   11871   .       +       0       ID=exon:Contig85g0000001.4;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004
#Contig85        EuGene  exon    11952   12072   .       +       2       ID=exon:Contig85g0000001.5;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004
#Contig85        EuGene  exon    12742   13220   .       +       1       ID=exon:Contig85g0000001.6;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000202
#Contig85        EuGene  five_prime_UTR  10248   10275   .       +       .       ID=five_prime_UTR:Contig85g0000001.0;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000204;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  CDS     10276   10402   .       +       0       ID=CDS:Contig85g0000001.1;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000196;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  CDS     11031   11276   .       +       2       ID=CDS:Contig85g0000001.2;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  CDS     11423   11625   .       +       2       ID=CDS:Contig85g0000001.3;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  CDS     11697   11871   .       +       0       ID=CDS:Contig85g0000001.4;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  CDS     11952   12072   .       +       2       ID=CDS:Contig85g0000001.5;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000004;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  CDS     12742   12757   .       +       1       ID=CDS:Contig85g0000001.6;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000197;est_cons=100.0;est_incons=0.0
#Contig85        EuGene  three_prime_UTR 12758   13220   .       +       .       ID=three_prime_UTR:Contig85g0000001.12;Parent=mRNA:Contig85g0000001;Ontology_term=SO:0000205;est_cons=100.0;est_incons=0.0

#Contig85        EuGene  gene    682876  682947  .       -       .       ID=gene:Contig85g0000771;Name=Contig85g0000771;length=72
#Contig85        EuGene  ncRNA   682876  682947  .       -       .       ID=ncRNA:Contig85g0000771;Name=Contig85g0000771;Parent=gene:Contig85g0000771;nb_exon=1;length=72
my @tag = ("length","nb_exon","Ontology_term","est_cons","est_incons","ensembl_end_phase", "ensembl_phase","rank","exon_id","protein_id","constitutive");
my %type;
my %gene;
my %mrna;
my %ncrna;
my %trna;
my %cds;
my %utr;
my %exon; 
while (my $feature = $gff_in->next_feature) {
    foreach my $tag (@tag) {
        if ($feature->has_tag($tag)) {
            $feature->remove_tag($tag);
        }        
    }
    if ($feature->primary_tag() eq "gene") { 
        $feature->remove_tag("Alias1similar")  if $feature->has_tag("Alias1similar");   
        $feature->remove_tag("Alias2similar")  if $feature->has_tag("Alias2similar");
        unless ($feature->has_tag("Name")){
            if ($feature->has_tag("gene_id")){
                my ($gene_id) = $feature->get_tag_values("gene_id");
                $feature->remove_tag("gene_id");
                $feature->add_tag_value("Name",$gene_id);
            }
        }
        push @{$gene{$feature->seq_id}} , $feature;
    }
    elsif ($feature->primary_tag() eq "mRNA") { 
        my ($mrna_id) = $feature->get_tag_values("ID"); 
        my ($parent_id) = $feature->get_tag_values("Parent");
        unless ($feature->has_tag("Name")){
            if ($feature->has_tag("transcript_id")){
                my ($transcript_id) = $feature->get_tag_values("transcript_id");
                $feature->remove_tag("transcript_id");
                $feature->add_tag_value("Name",$transcript_id);
            }
        }
        $feature->remove_tag("_QI")  if $feature->has_tag("_QI");   
        $feature->remove_tag("_AED")  if $feature->has_tag("_AED");   
        $feature->remove_tag("_eAED")  if $feature->has_tag("_eAED");   
        $type{$parent_id} = "mRNA"; 
        $mrna{$parent_id}{$mrna_id} = $feature;
    } 
    elsif ($feature->primary_tag() eq "ncRNA") {
        my ($mrna_id) = $feature->get_tag_values("ID"); 
        my ($parent_id) = $feature->get_tag_values("Parent");
        $type{$parent_id} = "ncRNA";
        $ncrna{$parent_id}{$mrna_id} = $feature;
    }  
    elsif  ($feature->primary_tag() eq "CDS") { 
        my ($parent_id) = $feature->get_tag_values("Parent"); 
        $feature->remove_tag("ID")  if $feature->has_tag("ID"); 
        $feature->remove_tag("Name")  if $feature->has_tag("Name");        
        push @{$cds{$parent_id}} , $feature;
    } 
    elsif  ($feature->primary_tag() eq "exon") { 
        my ($parent_id) = $feature->get_tag_values("Parent");  
        $feature->remove_tag("ID") if $feature->has_tag("ID"); 
        $feature->remove_tag("Name")  if $feature->has_tag("Name");       
        push @{$exon{$parent_id}} , $feature;
    } 
    elsif  ($feature->primary_tag() =~/UTR/) { 
        my ($parent_id) = $feature->get_tag_values("Parent"); 
        $feature->remove_tag("ID")  if $feature->has_tag("ID");      
        $feature->remove_tag("Name")  if $feature->has_tag("Name");  
        push @{$utr{$parent_id}} , $feature;
    } 
}
$gff_in->close;
my %correspondance;
foreach my $seq_id (sort {$a cmp $b} keys %gene) {
    my $cpt_gene = 0;
    foreach my $feature_gene (sort {$a->start <=> $b->start} @{$gene{$seq_id}}) {
        my ($gene_id) = $feature_gene->get_tag_values("ID");
        my ($name_id) = $feature_gene->get_tag_values($tag);
        next if $to_remove{$name_id};
        if ($type{$gene_id} eq "mRNA") { 
            $cpt_gene++;
            my $new_gene_id = sprintf( "%sg%06d",$seq_id, $cpt_gene * 10 );
            if ($rename) {
                $feature_gene->remove_tag("ID");
                $feature_gene->remove_tag("Name");
                $feature_gene->add_tag_value("ID",$new_gene_id);
                $feature_gene->add_tag_value("Name",$new_gene_id);
            }
            $feature_gene->add_tag_value("Note",$locus_info{$name_id}{product}) if $locus_info{$name_id}{product};   
            $gff_out->write_feature($feature_gene);
            my $count_mrna = 0;
            foreach my $mrna_id (keys %{$mrna{$gene_id}}) {
                $count_mrna++;
                my $new_mrna_id = sprintf( "%st%06d",$seq_id, $cpt_gene * 10 ) .".".$count_mrna;
                my $feature_rna = $mrna{$gene_id}{$mrna_id}; 
                my ($mrna_name_id) = $feature_rna->get_tag_values($tag);
                if ($rename ) {
                    $feature_rna->remove_tag("ID");
                    $feature_rna->remove_tag("Name");
                    $feature_rna->remove_tag("Parent");
                    $feature_rna->add_tag_value("ID",$new_mrna_id);
                    $feature_rna->add_tag_value("Name",$new_mrna_id);
                    $feature_rna->add_tag_value("Parent",$new_gene_id);
                }                 
                if (defined $locus_info{$mrna_name_id}) {
                    $feature_rna->remove_tag("Note") if $feature_rna->has_tag("Note");
                    $feature_rna->add_tag_value("Note",$locus_info{$mrna_name_id}{product})  if $locus_info{$mrna_name_id}{product};
                    $feature_rna->add_tag_value("Dbxref",$locus_info{$mrna_name_id}{dbxref}) if $locus_info{$mrna_name_id}{dbxref} &&  $locus_info{$mrna_name_id}{dbxref} ne "N/A";
                    #$feature_rna->add_tag_value("Ontology_term","CC_functional_completeness:".$locus_info{$name_id}{completeness}); 
                    #$feature_rna->add_tag_value("Ontology_term","CC_evidence_code:".$locus_info{$name_id}{ic_code});
                }
    
                if (defined $ipr{$mrna_name_id}){
                    foreach my $ipr_id (@{$ipr{$mrna_name_id}}) {
                        $feature_rna->add_tag_value("Dbxref","InterPro:".$ipr_id);   
                    }
                }
                if (defined $go{$mrna_name_id}) {
                    foreach my $go_id (@{$go{$mrna_name_id}}) {
                        $feature_rna->add_tag_value("Ontology_term",$go_id);      
                    }                
                }
                $gff_out->write_feature($feature_rna);
                my @cds = @{$cds{$mrna_id}};
                my @exon = @{$exon{$mrna_id}}; 
                my @utr = @{$utr{$mrna_id}};   
                foreach my $feature_cds (sort {$a->start <=>$b->start} @cds , @exon,@utr) { 
                    if ($rename ) {
                        $feature_cds->remove_tag("Parent");
                        $feature_cds->add_tag_value("Parent",$new_mrna_id);
                    }
                    $gff_out->write_feature($feature_cds);
                } 
            }
        } 
        else { 
            $gff_ncrna_out->write_feature($feature_gene);
            foreach my $mrna_id (keys %{$ncrna{$gene_id}}) {
                my $feature_rna = $ncrna{$gene_id}{$mrna_id}; 
                $gff_ncrna_out->write_feature($feature_rna);
                if (@{$exon{$mrna_id}}){ 
                    my @exon = @{$exon{$mrna_id}}; 
                    foreach my $feature_cds (sort {$a->start <=>$b->start}  @exon ) { 
                        $gff_ncrna_out->write_feature($feature_cds);
                    }
                } 
            }
        }
    }
}
$gff_out->close;
$gff_ncrna_out->close;

 
