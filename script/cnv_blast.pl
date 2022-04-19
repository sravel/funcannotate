#!/usr/bin/perl 
use File::Basename;
use FindBin;
use Pod::Usage; 
use Getopt::Long;
use Data::Dumper;
use Cwd;
my $directory   = cwd();
my $verbose     = 0;
my $current_dir = getcwd; 

my $blast;
my $product;
my $gene; 
my $help;
my $output = "locus_tag_with_product.txt";
my @source;
my %cc_gene;
my %genedb_products;  

my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters
       --blast    Blast [required]
       --out      Output file name (Default locus_tag_with_product.txt)
       --help    show this help and exit
/;
GetOptions(
     'blast=s@'  => \$blast,      
     'out=s'     => \$output, 
     'help|h|?'  => \$help
)  or pod2usage(-message => $usage);
if ($help) { pod2usage(-message => $usage); }

if ($blast eq "") {
    warn "\nWarn :: --blast is empty. Please specify a tabular diamond or blastp file\n\n Use outfmt 6 qseqid qlen sseqid slen stitle pident qcovhsp scovhsp  sstart send length\n\n";
    warn $usage;
    exit 0;
}      

foreach my $blast_result ( @$blast) {
     my ($name,$path,$suffix) = fileparse($blast_result);
     ($product,$gene,$source) = &parse_blastp($blast_result,$product,$gene,$name);
     push @source , $source;
} 

open(OUT,">$output"); 

foreach my $locus (sort keys %$product) {
     my $function = "hypothetical protein";
     my $cds      = "missing_functional_completeness";
     my $evidence_code   = "ISS_5";
     my $origin   = "N/A";
     my @gene;  
     my @match;
     foreach my $source (@source) { 
		push @match, @{$product->{$locus}->{$source}};
		push @gene , @{$gene->{$locus}->{$source}};
     } 
     my @sorted =  sort { $a->{evidence_code} cmp $b->{evidence_code} } @match; 
     my $best_match = $sorted[0];
     $function      = $best_match->{product};
     $cds           = $best_match->{cds_type};
     $evidence_code = $best_match->{evidence_code};
     my $bank = $best_match->{bank};
     my $qlen = $best_match->{query_length};
     my $hit_name = $best_match->{hit_name};
     my $dbxref = $best_match->{dbxref};
     $origin  = $source;	 
     @gene = @{&sort_uniq(\@gene)};
     my $gene_name = @gene ? join(",",@gene) : "unknown_gene";
     if(($evidence_code eq "ISS_5") || ($evidence_code eq "ISS_6")) {
     	print OUT join("\t",$locus,$hit_name,$bank, $qlen,$function, $evidence_code ,$cds,"unknown_gene",$dbxref),"\n";
     }
     else {
     	print OUT join("\t",$locus,$hit_name, $bank,$qlen, $function, $evidence_code ,$cds,$gene_name,$dbxref),"\n"; 
     } 
}    
close OUT; 


#qseqid qlen sseqid slen pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp
sub parse_blastp {
    my ($file,$product,$gene_name,$bank) = @_;
    $bank = (split(/\_/,$bank))[0]; 
    
    my %dbxref;  
    my %product = %$product;
    my %gene_name = %$gene_name; 
    open(IN,$file);
    while(<IN>){
        chomp;
        next if $_ =~ /^#/;
        my ($qseqid ,$qlen ,$sseqid ,$slen ,$stitle ,$evalue,$pident,$qcovhsp ,$scovhsp ,$sstart ,$send ,$length) = (split(/\t/,$_)); 
        if ($sseqid eq "*") {
            push @{$product{$qseqid}{$bank}} , { 
                product       => "N/A", 
                cds_type      => "N/A",
                query_length  => $qlen,
                qcov 	      => "N/A",
                scov	      => "N/A",
                identity      => "N/A",
                evidence_code => "ISS_6",
                hit_name      => "N/A",
                bank          => "N/A",
                description   => "N/A"							
            };        
        }
        else {
            my $description = $stitle;
            $description =~s/$sseqid\s//;
            #if ($bank eq "TrEMBL"){
            #     
            #    $sseqid = (split(/\|/,$sseqid))[1];
            #    $description = (split(/\_\w{5}\s/,$description))[1];
            #    
            #}
            #if ($scovhsp == "") {
            #    $scovhsp = ($send - $sstart) / $slen;
            #} 
            my ($species,$gene_name,$dbxref,$name,$alias);	
            if ($bank =~ "SwissProt" || $bank eq "TrEMBL") { 
                if ($description =~ /(.*)\sOS=(.*)\sGN=(\S+)\s.*/) {
                    $description = $1;
                    $species = $2;
                    $gene_name = $3;  
                    $gene_name =~ s/,/./;
                    push @{$gene_name{$qseqid}{$bank}} , $gene_name ;
                }
                elsif ($description =~ /(.*)\sOS=(.*)\sPE=/) {
                    $description = $1;
                    $species = $2;
                } 
                    
                $dbxref = $bank.":".$sseqid; 
            }
            else {
                $dbxref = $bank.":".$sseqid;
            }
            my ($product,$cds_type, $evidence_code) = &eval_prediction($qcovhsp,$scovhsp,$pident,$description,$evalue); 
            push @{$product{$qseqid}{$bank}} , {
                product       => $product,
                query_length  => $qlen,
                cds_type      => $cds_type,
                qcov 	      => $qcovhsp,
                scov	      => $scovhsp,
                identity      => $pident,
                bank          => $bank,
                evidence_code => $evidence_code,
                hit_name      => $sseqid,
                description   => $description,
                dbxref        => $dbxref
            };   
        }  
    }
    return (\%product,\%gene_name,$bank);
}			

sub _max (@) {
    my $i = shift;
    foreach (@_) {
		$i = $_ if $_ > $i;
    }
    return $i;
}


sub _min (@) {
    my $i = shift;
    foreach (@_) {
		$i = $_ if $_ < $i;
    }
    return $i;
}

sub eval_prediction {
    my ( $QCov, $SCov, $identity, $description, $Evalue ) = @_;  
    my $evidence_code = "ISS_5";
    my $product = "hypothetical protein";
    my $cds_type = "missing_functional_completeness"; 
    #SSB comment: particularity of NCBI100 annotation for instance REF_ELAGV
    #>XP_010909149.1|LOC105035081|REF_EGV01_PF00870.1_ELAGV| PREDICTED: butyrate--CoA ligase AAE11, peroxisomal-like
    $description =~ s/PREDICTED: //;
    #SSB comment: nowaday (20161106) with all the sequences that are present in the protein databank, those having a low identity become doubtful. They could corespond to functionnal RNA such as long nc-RNA
    #SSB comment 2016111: finally I prefer to leave the ISS_6 for those with not hit found see above
    if( $Evalue <= 0.001 ) {
        if( $identity < 25 ) {
            $evidence_code = "ISS_5";
            $product = "hypothetical protein";
        }
        #SSB comment: all the hypothetical or uncharacterized prot are treated. The evidence_code of those having a Cov >= 0.75 (Q or S) become ISS_4 and one the three cds_type is assigned (complete, fragment, modules). The remaining (Qcov & Scov < 0.75) stay with ISS_5
        elsif (( $description =~ /hypothetical protein/i ) || ( $description =~ /uncharacterized protein/i ) || ( $description =~ /^unknown protein$/i ) || ( $description =~ /^predicted protein$/i )) {
            if ( $QCov >= 80 ) {
                $product = "conserved hypothetical protein";
                $evidence_code = "ISS_4";
                if ( $SCov >= 80 ) {
                	$cds_type = "complete";
                }
                else {
                    $cds_type = "fragment";
                }
            }
            elsif ( $SCov >= 80 ) {
                $product = "conserved hypothetical protein";
                $evidence_code = "ISS_4";
                $cds_type = "modules";
            }
        }
        #SSB comment: all the putative, probable, predicted or LOW QUALITY PROTEIN prot are treated. LOW QUALITY PROTEIN is a particularity of NCBI100 annotation (e.g. REF_ELAGV). The evidence_code of those having a Cov >= 0.75 (Q or S) become ISS_3 and one the three cds_type is assigned (complete, fragment, modules). The remaining (Qcov & Scov < 0.75) stay with ISS_5
        elsif (( $description =~ /putative/i ) || ( $description =~ /probable/i ) || ( $description =~ /probably/i ) || ( $description =~ /predicted/i ) || ( $description =~ /LOW QUALITY PROTEIN/i ) || ( $description =~ /^Unknown /i )) {
            if ( $QCov >= 80 ) {
                $evidence_code = "ISS_3";
                $product = $description;
                #SSB comment: 20191029 do not factorise this if with the one below only polypetides with enough coverage will become iss_3 the other stay ISS_4
                if ( $description =~ /LOW QUALITY PROTEIN: putative/i ) {
                    $product =~ s/LOW QUALITY PROTEIN: putative/putative/;
                }
                elsif ( $description =~ /LOW QUALITY PROTEIN: probable/i ) {
                    $product =~ s/LOW QUALITY PROTEIN: probable/putative/;
                }	       
                elsif ( $description =~ /LOW QUALITY PROTEIN:/i ) {
                    $product =~ s/LOW QUALITY PROTEIN:/putative/;
                }	  
                elsif ( $product =~ /putative/i ) {
                    $product =~ s/Putative/putative/;
                }
                elsif ( $description =~ /probable/i ) {
                    $product =~ s/[Pp]robable/putative/;
                }
                elsif ( $description =~ /probably/i ) {
                    $product =~ s/[Pp]robably/putative/;
                }	       
                elsif ( $description =~ /predicted/i ) {
                    $product =~ s/[Pp]redicted/putative/;
                }
                elsif ( $description =~ /unknown/i ) {
                    $product =~ s/[Uu]nknown/putative/;
                }	       
                else {$product = "BUG1 $description";}
                if ( $SCov >= 80 ) {
                    $cds_type = "complete";
                }
                else {
                    $cds_type = "fragment";
                }
            }
            elsif ( $SCov >= 80 ) {
                $evidence_code = "ISS_3";
                $product = $description;
                if ( $description =~ /LOW QUALITY PROTEIN: putative/i ) {
                    $product =~ s/LOW QUALITY PROTEIN: putative/putative/;
                }
                elsif ( $description =~ /LOW QUALITY PROTEIN: probable/i ) {
                    $product =~ s/LOW QUALITY PROTEIN: probable/putative/;
                }	       
                elsif ( $description =~ /LOW QUALITY PROTEIN:/i ) {
                    $product =~ s/LOW QUALITY PROTEIN:/putative/;
                }	  
                elsif ( $product =~ /putative/i ) {
                    $product =~ s/Putative/putative/;
                }
                elsif ( $description =~ /probable/i ) {
                    $product =~ s/[Pp]robable/putative/;
                }
                elsif ( $description =~ /probably/i ) {
                    $product =~ s/[Pp]robably/putative/;
                }	       
                elsif ( $description =~ /predicted/i ) {
                    $product =~ s/[Pp]redicted/putative/;
                }
                #putative Unknown protein from spot 77 of 2D-PAGE of etiolated coleoptile (Fragment)
                elsif ( $description =~ /unknown/i ) {
                    $product =~ s/[Uu]nknown/putative/;
                }	       
                else {$product = "BUG2 $description";}
                $cds_type = "modules";
            }
        }
        #SSB comment: all sequence identities >=0.25 not having the keywords: hypothetical, uncharacterized, putative, probable, predicted or LOW QUALITY PROTEIN prot are treated. By default they have the evidence_code ISS_4. The evidence_code of those having a Cov >= 0.75 (Q or S) become ISS_3 and one the three cds_type is assigned (complete, fragment, modules). The remaining (Qcov & Scov < 0.75) stay with ISS_4. Finally regarding the identity ISS_2 and ISS1 are assigned
        else {
            $evidence_code = "ISS_4";
            $product = "conserved hypothetical protein";		
            if ( $QCov >= 80 ) {
                $evidence_code = "ISS_3";
                $product = "putative $description";
                if ( $SCov >= 80 ) {
                    $cds_type = "complete"; 
                }
                else {
		            $cds_type = "fragment";
                }
            }
            elsif ( $SCov >= 80 ) {
                $evidence_code = "ISS_3";
                $product = "putative $description";
                $cds_type = "modules";
            } 
            #SSB comment: 20191022 bug saw by Gaetan Droc if ( $identity >= 0.45 ) alone was not a sufficient condition so too many ISS_2 were predicted
            if (( $identity >= 45 ) && ( $QCov >= 80 ) &&( $SCov >= 80 )) {
                $evidence_code = "ISS_2";
                $product = $description;
            }
            if (( $identity >= 90 ) && ( $QCov >= 90 ) &&( $SCov >= 90 )) {
                $product = $description;
                $evidence_code = "ISS_1";
            }
        }
 
        if ( $product =~ /, partial$/ ) {
            $product =~ s/, partial$//;
            $cds_type = "fragment";
        }         
        elsif ( $product =~ / \(Fragments\)$/i ) {
            $product =~ s/ \([Ff]ragments\)$//;
            $cds_type = "fragment";
        }  
        elsif ( $product =~ / \(Fragment\)$/i ) {
            $product =~ s/ \([Ff]ragment\)$//;
            $cds_type = "fragment";
        }
        #SSB comment 20200218: Truncated betaine aldehyde dehydrogenase 2
        elsif ( $product =~ / \(Truncated\)$/i ) {
            $product =~ s/ \([Tt]runcated\)$//;
            $cds_type = "fragment";
        } 
     
        if ( $product =~ / homolog \w+-like isoform X\d+$/ ) { 
            $product =~ s/ homolog \w+-like isoform X\d+$//;
        }         
        #COCNU_HT001RO01_01P009610	putative cell division cycle protein 27 homolog B isoform X1     
        elsif ( $product =~ / homolog \w+ isoform X\d+$/ ) {
            $product =~ s/ homolog \w+ isoform X\d+$//;
        }       
        elsif ( $product =~ / homolog isoform X\d+$/ ) { 
            $product =~ s/ homolog isoform X\d+$//;
        } 
        #putative Glycogen synthase kinase-3 homolog MsK-3
        elsif ( $product =~ / homolog \w+-\d$/ ) { 
            $product =~ s/ homolog \w+-\d$//;
        }    
        #vacuolar iron transporter homolog 2.1                
        elsif ( $product =~ / homolog \w+\.\d$/ ) { 
            $product =~ s/ homolog \w+\.\d$//;
        }      
        elsif ( $product =~ / homolog \w+-like$/ ) { 
            $product =~ s/ homolog \w+-like$//;
        }                    
        elsif ( $product =~ / homolog \d+$/ ) { 
            $product =~ s/ homolog \d+$//;
        }                     
        elsif ( $product =~ / homolog \w$/i )  {
            $product =~ s/ homolog \w$//;
            $product =~ s/ HOMOLOG \w$//;
        } 
        elsif ( $product =~ / homolog \w+$/ ) { 
            $product =~ s/ homolog \w+$//;
        } 
        elsif ( $product =~ / homolog$/ ) { 
            $product =~ s/ homolog$//;
        }                     
        elsif ( $product =~ / isoform X\d+$/ ) { 
            $product =~ s/ isoform X\d+$//;
        } 
        #else{print "None\t";}
        #protein HOMOLOG OF MAMMALIAN LYST-INTERACTING PROTEIN 5
        if ( $product =~ /protein HOMOLOG OF /i ) { 
            $product =~ s/protein HOMOLOG OF //;
            $product =~ s/protein homolog of //;
            $product =~ s/PROTEIN HOMOLOG OF //;
        }
        #COCNU_HT001RO01_02P015620	ATP-dependent Clp protease ATP-binding subunit ClpA homolog CD4B, chloroplastic-like
        elsif ( $product =~ / homolog \w+,/ ) { 
            $product =~ s/ homolog \w+,/,/;
        }
        #Cell division protein FtsZ homolog 2-2, chloroplastic
        elsif ( $product =~ / homolog \w+-\d,/ ) { 
            $product =~ s/ homolog \w+-\d,/,/;
        } 
        elsif ( $product =~ / homolog \w+\.\d,/ ) { 
            $product =~ s/ homolog \w+\.\d,/,/;
        }                               
        elsif ( $product =~ / homolog,/ ) { 
            $product =~ s/ homolog,/,/;
        } 
        elsif ( $product =~ / homolog protein/ ) { 
            $product =~ s/ homolog protein/-like protein/;
        } 
        elsif ( $product =~ / homologous protein/ ) { 
            $product =~ s/ homologous protein/-like protein/;
        } 
        elsif ( $product =~ / homolog subfamily/ ) { 
            $product =~ s/ homolog subfamily/-like subfamily/;
        }
        elsif ( $product =~ / homology domain/ ) { 
            $product =~ s/ homology domain/ domain/;
        }
        elsif ( $product =~ / homologous domain/ ) { 
            $product =~ s/ homologous domain/ domain/;
        }     
        if ( $product =~ /Pseudo histidine/i ) {
            $product =~ s/[Pp]seudo histidine/pseudohistidine/;
        }          
        if ( $product =~ /transferases/i ) {
            $product =~ s/transferases/transferase/;
        }
        #Vascular-related unknown protein 1 
        if ( $product =~ / unknown protein/ ) {
            $product =~ s/ unknown protein/protein/;
        } 
        #SSB comment 20200218: putative Protein LIFEGUARD 1
        $product =~ s/putative Protein/putative protein/;
        #putative 58. protein in lys 3'region https://www.uniprot.org/uniprot/A0A1D1XSG9 https://www.uniprot.org/uniprot/P51739
        $product =~ s/58\. protein/58 kDa protein/;
    }#SSB comment 20191119: end of if($Evalue <= 0.001) 
    return ( $product, $cds_type, $evidence_code );
}
 

sub sort_uniq {
    my $array = shift;
    my %seen = ();
    my @array = grep { !$seen{$_}++ } @{$array};
    return \@array;
}

sub encod {
	my  $encod = shift;
	$encod =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!-])/uc sprintf("%%%02x",ord($1))/eg; 
	return $encod;
}
