import subprocess,sys,re
from os.path import join,basename
import pandas as pd
from Bio import SearchIO
import sys
import re
import os
import snakemake
from snakemake import rules
from snakemake.io import expand, multiext, glob_wildcards, directory, temp, unpack
from collections import defaultdict, OrderedDict
from pprint import pprint as pp
from pathlib import Path
#from snakemake.utils import R
# import ipdb

def var_2_bool(key, tool, to_convert):
    """convert to boolean"""
    if isinstance(type(to_convert), bool):
        return to_convert
    elif f"{to_convert}".lower() in ("yes", "true", "t"):
        return True
    elif f"{to_convert}".lower() in ("no", "false", "f"):
        return False
    elif to_convert == 1:
        return False
    elif to_convert == 0:
        return True
    else:
        raise TypeError(
            f'CONFIG FILE CHECKING FAIL : in the "{key}" section, "{tool}" key: "{to_convert}" is not a valid boolean')


###Take the output of hmmsearch, a pvalue cutoff, and return a pandas df with matching genes in query and their associated evalue.
def parse_hmmsearch_results(hmm_file, pvalue_cutoff):
    list_evalue = []
    list_id = []
    print("Starting")
    print(hmm_file)
    for qresult in SearchIO.parse(hmm_file, 'hmmer3-text'):
        hits = qresult.hits
        print("Hits...")
        print(hits)
        if len(hits) > 0:
            for i in range(0,len(hits)):
                if hits[i].evalue < pvalue_cutoff:
                    hit_evalue = hits[i].evalue
                    hit_id = hits[i].id
                    list_evalue.append(hit_evalue)
                    list_id.append(hit_id)
    return(pd.DataFrame({'evalue': list_evalue}, index=list_id))



# Globals ---------------------------------------------------------------------
SAMPLE,FILES,   = glob_wildcards(f"{config['SAMPLES']}{{sample}}/{{files}}", followlinks=True)
SAMPLES = set(SAMPLE)
# print(FILES)

dico_samples = defaultdict(OrderedDict)
for smp in SAMPLES:
    for file_in in FILES:
        if smp in file_in:
            file_path = Path(config['SAMPLES']).joinpath(smp, file_in)
            if "_prot" in  file_in:
                dico_samples[smp]["protein"] = file_path.as_posix()
            if file_path.suffix in ".gff3":
                dico_samples[smp]["gff"] =  file_path.as_posix()
            if file_path.suffix  in ".fasta" and not "_prot" in  file_in:
                dico_samples[smp]["fasta"] =  file_path.as_posix()
# pp(dico_samples)

###Perform some verifications
fasta_list = [True if "fasta" in dico.keys() else False for key, dico in dico_samples.items() ]
gff_list = [True if "gff" in dico.keys() else False for key, dico in dico_samples.items() ]
prot_list = [True if "protein" in dico.keys() else False for key, dico in dico_samples.items() ]
if len(dico_samples) != fasta_list.count(True):
    raise SystemExit("Error, Not provided fasta for all Samples.")
if gff_list.count(True) != fasta_list.count(True):
    raise SystemExit("Error, Not provided fasta for all gff3 maybe genome names are different. Please make sure that they are identical.")
if prot_list.count(True) != fasta_list.count(True):
    raise SystemExit("Error, Not provided fasta for all protein maybe genome names are different. Please make sure that they are identical.")


blast_databases_annot = config["blast_databases"]
blast_databases_repeat = config["blast_repeat_database"]
hmm_databases = config["hmm_databases"]

### Merge all the blast databases (repeat or not) in one structure
blast_databases = blast_databases_annot.copy()

only_reannotate = var_2_bool(key="only_reannotate",tool="only_reannotate", to_convert=config["only_reannotate"])
# print(only_reannotate)
if only_reannotate:
    blast_databases.update(blast_databases_repeat)

blast_threads=config["params"]["global_blast_threads"]
# print(blast_databases_annot)
# print(blast_databases_repeat)
# print(blast_databases)
# print(hmm_databases)

#######################################
# Change workdir to output path
workdir:config["OUTPUT"]
outputdir = config["OUTPUT"]

###Indices for the fasta splited files, in the form of '001', '002'...
my_indices=[f"{i:03d}" for i in range(1, config["fasta_split_number"]+1)]

rule all:
    input:
        expand("pipeline_report_{sample}.html", sample=SAMPLES),
        expand("enriched_gff3/{sample}/{sample}_proteins.fasta.fai", sample=SAMPLES) if False not in fasta_list else expand("enriched_gff3/{sample}/{sample}_mrna.gff3", sample=SAMPLES)

def get_fasta(wildcards):
    return {"fasta": dico_samples[f"{wildcards.sample}"]["fasta"]}

def get_gff(wildcards):
    return {"gff3": dico_samples[f"{wildcards.sample}"]["gff"]}

def get_prot(wildcards):
    return {"protein": dico_samples[f"{wildcards.sample}"]["protein"]}

def get_db_name(wildcards):
    if wildcards.db_name in blast_databases:
        return {"db_name": blast_databases[f"{wildcards.db_name}"]}
    if wildcards.db_name in blast_databases_repeat:
        return {"db_name": blast_databases_repeat[f"{wildcards.db_name}"]}

wildcard_constraints:
    sample= "|".join(SAMPLES)


rule fasta_symlink:
    input:
        unpack(get_fasta),
    output:
       fasta = "input_fasta/{sample}.fasta"
    shell:
        "ln -s {input} {output}"

###########################################
### Create fasta index (used to get sequence length)
rule fasta_index:
    input:
        fasta = rules.fasta_symlink.output.fasta
    output:
        index = "input_fasta/{sample}.fasta.idx" ###FIXME : .fai... not idx !
    envmodules:
        config["modules"]["samtools"]
    log:
        out="logs/fasta_index/{sample}.log",
        err="logs/fasta_index/{sample}.err"
    shell:
        "samtools faidx --fai-idx {output.index} {input.fasta} > {log.out} 2> {log.err}"


####################################
### Split files in n pieces
rule split:
    input:
        fasta = rules.fasta_symlink.output.fasta
    output:
        list_fasta = expand("tmp_fasta/{{sample}}/{{sample}}.part_{n}.fasta",sample="{sample}", n=my_indices)
    params:
        number_of_output_fasta = config["fasta_split_number"],
        sample = "{sample}"
    envmodules:
        config["modules"]["seqkit"]
    log:
        out="logs/split/{sample}.log",
        err="logs/split/{sample}.err"
    shell:
        "seqkit split --out-dir tmp_fasta/{params.sample} --by-part {params.number_of_output_fasta} {input.fasta} > {log.out} 2> {log.err}"

####################################
### Blast Part

rule build_blast_databases:
    input:
        unpack(get_db_name)
    output:
        done=touch("blast/blast_databases/{db_name}.makeblastdb.done")
    envmodules:
        config["modules"]["blast"]
    log:
        out="logs/build_blast_databases/{db_name}.log",
        err="logs/build_blast_databases/{db_name}.err"
    shell:
        "makeblastdb -parse_seqids -in {input} -dbtype prot  -out blast/blast_databases/{wildcards.db_name}   > {log.out} 2> {log.err}"

rule build_diamond_database:
    input:
        blast_db=rules.build_blast_databases.output.done
    output:
        done=touch("blast/blast_databases/{db_name}.diamonddb.done")
    envmodules:
        config["modules"]["diamond"]
    log:
        out="logs/build_diamond_database/{db_name}.log",
        err="logs/build_diamond_database/{db_name}.err"
    shell:
        "diamond prepdb -d blast/blast_databases/{wildcards.db_name} > {log.out} 2> {log.err}"


rule diamond_analysis:
    input:
        fasta_files=outputdir+"tmp_fasta/{sample}/{sample}.part_{indice}.fasta",
        fasta_list = rules.split.output.list_fasta,
        blast_db= rules.build_diamond_database.output.done,
    output:
        tsv = temp("blast/raw_diamond_output/{sample}/{db_name}_file_{sample}.part_{indice}.tsv")
    envmodules:
        config["modules"]["diamond"]
    params:
        evalue=1e-10,
        max_target_seqs=2,
        mode="--very-sensitive",
        output_format="--outfmt 6 qseqid qlen sseqid slen stitle evalue pident qcovhsp scovhsp  sstart send length"
    threads:
        config["params"]["global_diamond_threads"]
    log:
        out="logs/diamond_computation/{sample}_{indice}_{db_name}.log",
        err="logs/diamond_computation/{sample}_{indice}_{db_name}.err"
    shell:
        "diamond blastp  --masking 0 --unal 1 {params.output_format}  {params.mode} --threads {threads}  --db blast/blast_databases/{wildcards.db_name} --query {input.fasta_files} --evalue {params.evalue} -k {params.max_target_seqs}  -o {output} > {log.out} 2> {log.err}"

rule blastp_computation:
    input:
        fasta_files=outputdir+"tmp_fasta/{sample}/{sample}.part_{indice}.fasta",
        fasta_list= rules.split.output.list_fasta
    output:
        tsv = temp("blast/raw_output/{sample}/{db_name}_file_{sample}.part_{indice}.tsv")
    params:
        blast_db = rules.build_blast_databases.output.done,
        evalue=1e-10,
        max_target_seqs=2,
        output_format="-outfmt '6 qseqid qlen sseqid slen stitle evalue pident qcovhsp scovhsp  sstart send length'",
        db_name = "{db_name}"
    envmodules:
        config["modules"]["blast"]
    threads:
        config["params"]["global_blast_threads"]
    log:
        out="logs/blastp_computation/{sample}_{indice}_{db_name}.log",
        err="logs/blastp_computation/{sample}_{indice}_{db_name}.err"
    shell:
        "blastp -db blast/blast_databases/{params.db_name} {params.output_format} -query {input.fasta_files} -out {output} -num_threads {threads} -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} > {log.out} 2> {log.err}"


rule merge_blast_results:
    input:
        expand(rules.blastp_computation.output.tsv, sample="{sample}", db_name="{db_name}", indice=my_indices) if config["use_diamond"] == 0 else  expand(rules.diamond_analysis.output.tsv, sample="{sample}", db_name="{db_name}", indice=my_indices)
    output:
        "blast/merged_output/{db_name}_{sample}.tsv"
    log:
        err="logs/merge_blast_results/{db_name}_{sample}.err"
    shell:
        "cat {input} > {output} 2> {log.err}"

####################################
### Compute ISS

###TODO : Need to find which ISS is the best before going to the filtering part. We need to have like one file "Best ISS prot" and one other "Best ISS TE"
rule compute_ISS_proteins:
    input:
        db_annot=expand("blast/merged_output/{db_name}_{{sample}}.tsv", db_name=list(blast_databases_annot.keys()))
    output:
        iss = "blast/ISS_informations/results_{sample}.iss"
    params:
        db_annotation=' --blast '.join(expand("blast/merged_output/{db_name}_{{sample}}.tsv",db_name=list(blast_databases_annot.keys())))
    envmodules:
        config["modules"]["func_annot"]
    log:
        out="logs/compute_ISS_proteins/{sample}.log",
        err="logs/compute_ISS_proteins/{sample}.err"
    shell:
        "cnv_blast.pl --blast {params.db_annotation} --out {output.iss}  > {log.out} 2> {log.err}"

rule compute_ISS_TE:
    input:
        db_repeat=expand("blast/merged_output/{db_name}_{{sample}}.tsv", db_name=list(blast_databases_repeat.keys()))
    output:
        iss = "blast/ISS_informations/results_{sample}_TE.iss"
    params:
        db_repeat_TE=' --blast '.join(expand("blast/merged_output/{db_name}_{{sample}}.tsv",db_name=list(blast_databases_repeat.keys())))
    envmodules:
        config["modules"]["func_annot"]
    log:
        out="logs/compute_ISS_TE/{sample}.log",
        err="logs/compute_ISS_TE/{sample}.err"
    shell:
        "cnv_blast.pl  --blast {params.db_repeat_TE} --out {output.iss}  > {log.out} 2> {log.err}"

###TODO : Add length information to this file to improve filtering.
###TODO : this can be skip, just read the files in a python-like rule and merge them with python.
rule merge_iss:
    input:
        proteins=rules.compute_ISS_proteins.output.iss,
        TE=rules.compute_ISS_TE.output.iss
    output:
        merge_iss = "blast/ISS_informations/ISS_merged_{sample}.iss"
    log:
        err="logs/merge_iss/{sample}.err"
    shell:
        "join -t $'\t' {input.proteins} {input.TE} > {output.merge_iss} 2> {log.err}"

####################################
### Interpro Part

rule interpro_scan:
    input:
        fasta_files = rules.split.output.list_fasta,
    output:
        tsv = "interproscan/raw/{sample}/{sample}.part_{indice}.tsv",
        xml = "interproscan/raw/{sample}/{sample}.part_{indice}.xml"
    threads: config["params"]["interproscan_threads"]
    envmodules:
        config["modules"]["interproscan"]
    log:
        out="logs/interpro_scan/{sample}_{indice}.log",
        err="logs/interpro_scan/{sample}_{indice}.err"
    shell:
        "interproscan.sh -cpu {threads} --iprlookup --goterms -f tsv,xml -T ./test_tmp_dir_interpro/ --pathways --disable-precalc -i {input} --output-file-base interproscan/raw/{wildcards.sample}/{wildcards.sample}.part_{wildcards.indice} > {log.out} 2> {log.err}"
        # "interproscan.sh -cpu 8  --iprlookup --goterms -f tsv,xml - T FUCKING_TMP --disable-precalc -i {input} --output-dir output_interpro/{sample}/result_OLEEU_F9"

rule merge_interpro_results:
    input:
        expand(rules.interpro_scan.output.tsv, sample="{sample}", indice=my_indices)
    output:
        merge_tsv = "interproscan/merged/interpro_{sample}_merged.tsv"
    log:
        err="logs/merge_interpro_results/{sample}.err"
    shell:"""
        cat {input} > {output.merge_tsv} 2> {log.err}
    """    

# xml = expand("interproscan/raw/{{sample}}/{{sample}}.part_{indice}.xml", indice=my_indices)
# cat {input.xml} > {output.file_xml} 2> {log.err}
#TODO : check that the directory is created.
#TODO : may be replaced by a bash command.
rule run_cnv_interpro:
    input:
        expand(rules.interpro_scan.output.tsv, sample="{sample}", indice=my_indices)
    output:
        go="interproscan/cnv_interpro/{sample}_go.txt",
        ipro="interproscan/cnv_interpro/{sample}_ipr.txt"
    params:
        output_prefix="interproscan/cnv_interpro/{sample}",
        output_ipro="interproscan/raw/{sample}"
    envmodules:
        config["modules"]["func_annot"]
    log:
        out="logs/run_cnv_interpro/{sample}.log",
        err="logs/run_cnv_interpro/{sample}.err"
    shell:
        "cnv_interpro.pl --result {params.output_ipro} --prefix {params.output_prefix} > {log.out} 2> {log.err}"
        # "cnv_interpro.pl --result {params.output_ipro} --prefix {params.output_prefix} > {log.out} 2> {log.err}"

rule add_annotation_to_gff3:
    input:
        unpack(get_gff),
        # to_remove= "gene_info_table/{sample}_to_remove.tsv" if only_reannotate else "",  ###Totaly arbitraty file. Just have to be sure that this file exists and will always be there, so input fasta was a good choice.
# ISS_full = "blast/ISS_informations/ISS_merged_{sample}.iss" if config["include_repeat_database"]!=0 else "blast/ISS_informations/results_{sample}.iss", only_reannotate
        ISS_full = rules.merge_iss.output.merge_iss,
        gene_ontology = {rules.run_cnv_interpro.output.go},
        interpro_res = {rules.run_cnv_interpro.output.ipro},
    output:
        processed_gff = "enriched_gff3/{sample}/{sample}_mrna.gff3"
    envmodules:
        config["modules"]["func_annot"]
    params:
        # prefix="{sample}", ###Why not wildcard.sample in the shell command ?
        tag=config["tag"],
        to_remove_string="-te_file gene_info_table/{sample}_to_remove.tsv -rename" if only_reannotate else  ""
    shell:
        """cnv_eugene_gff3.pl -product {input.ISS_full} {params.to_remove_string} \
        -go_file {input.gene_ontology} \
        -interpro_file {input.interpro_res} \
        -gff3_file {input.gff3} \
        -prefix enriched_gff3/{wildcards.sample}/{wildcards.sample} \
        -tag {params.tag}"""

if False not in fasta_list:
    rule extract_prot_cds_sequences:
        input:
            rules.fasta_symlink.output.fasta,
            processed_gff = rules.add_annotation_to_gff3.output.processed_gff,
        output:
            proteins="enriched_gff3/{sample}/{sample}_proteins.fasta",
            cds="enriched_gff3/{sample}/{sample}_cds.fasta"
        envmodules:
            config["modules"]["gffread"]
        shell:
            "gffread -g {input.fasta} -x {output.cds} -y {output.proteins} {input.processed_gff}"

    rule index_extracted_sequences:
        input:
            proteins=rules.extract_prot_cds_sequences.output.proteins,
            cds=rules.extract_prot_cds_sequences.output.cds
        output:
            proteins="enriched_gff3/{sample}/{sample}_proteins.fasta.fai",
            cds="enriched_gff3/{sample}/{sample}_cds.fasta.fai"
        envmodules:
            config["modules"]["samtools"]
        shell:
            "samtools faidx {input.proteins} && samtools faidx {input.cds}"



####################################
### HMM Part

# def hmm_parse(hmm_results_file)
# for qresult in SearchIO.parse(file_hmm_gypsy, 'hmmer3-text'):
#     hits = qresult.hits
#     #print(hits)
#     if len(hits) > 0:
#         for i in range(0,len(hits)):
#                 if hits[i].evalue <0.01:
#                     beste = hits[i].evalue
#                     hit = hits[i].id.replace('.hmm', '')
#                     bestes_gypsy.append(beste)
#                     hitss_gypsy.append(hit)

###FIXME : useless to split, hmmsearch is very fast in our case because very few input sequences.
rule hmm_search:
    input:
        rules.fasta_symlink.output.fasta,
        hmm_file=lambda wildcards: hmm_databases[wildcards.hmm_db_name]
    output:
        "hmm/hmm_search/{sample}/{hmm_db_name}.txt"
    envmodules:
        config["modules"]["hmmer"]
    shell:
        "hmmsearch {input.hmm_file} {input.fasta} > {output}"

rule hmm_parse:
    input:
        hmm_file="hmm/hmm_search/{sample}/{hmm_db_name}.txt"
    output:
        parsed_hmm="hmm/parsed_hmm/{sample}/{hmm_db_name}_parsed.tsv"
    run:
        dataframe = parse_hmmsearch_results(input.hmm_file, 0.01)
        dataframe = dataframe.add_suffix("_"+wildcards.hmm_db_name)
        dataframe.to_csv(output.parsed_hmm, sep="\t")



###TODO : fusion de tous les hmm.tsv pour chaque sample et création d'un truc "to_remove".
###TODO : OR - fusion de chaquun de ces fichiers dans le gros tableau. Possibilité d'avoir un beau tableau complet.

# rule merge_hmm_results:
#     input:
#         expand("output_hmm_search/{{sample}}/{{sample}}_{{hmm_db_name}}.part_{indice}.tsv", indice=my_indices)
#     output:
#         "merged_hmm/hmm_{sample}_{hmm_db_name}_merged.tsv"
#     shell:
#         "cat {input} > {output}"



rule busco_input:
    input:
        rules.fasta_symlink.output.fasta,
    output:
        touch("busco/input/{sample}/done")
    threads:
        config["params"]["busco_threads"]
    envmodules:
        config["modules"]["busco"]
    params:
         lineage=config["lineage"],
         busco_path=config["busco_path"]
    log:
        out="logs/busco_input/{sample}.log",
        err="logs/busco_input/{sample}.err"
    shell: "busco -i {input.fasta} -f -m prot -l {params.lineage} --offline -c {threads} --download_path {params.busco_path} --out_path busco/input/ -o {wildcards.sample} > {log.out} 2> {log.err}"



rule merge_hmm_results:
    input:
        parsed_hmm=expand(rules.hmm_parse.output.parsed_hmm,sample="{sample}", hmm_db_name = list(hmm_databases.keys()))
    output:
        merged_hmm="hmm/hmm_merged/hmm_{sample}_full_merged.tsv"
    # log:
    #     out="logs/merge_hmm_results/{sample}.log",
    #     err="logs/merge_hmm_results/{sample}.err"
    run:
        dfs = []
        # ipdb.set_trace()
        for i in input.parsed_hmm:
            df = pd.read_table(i, index_col=0)
            if len(df) >= 1:
                dfs = df.join(dfs, how = "outer")
        if len(dfs) > 1:
            print('Multime HMM results')
            dfs.to_csv(output.merged_hmm, sep="\t", na_rep = "NaN")
        if len(dfs) == 1:
            print('Single HMM results')
            dfs[0].to_csv(output.merged_hmm, sep="\t")
        else:
            print("Only one HMM result... Sad.")
            open(output.merged_hmm,'a').close() ###If the file is empty

rule busco_output:
    input:
        fasta_prot = rules.extract_prot_cds_sequences.output.proteins
    output:
        done = touch("busco/output/{sample}/done")
    threads:
        config["params"]["busco_threads"]
    envmodules:
        config["modules"]["busco"]
    params:
         lineage=config["lineage"],
         busco_path=config["busco_path"]
    log:
        out="logs/busco_output/{sample}.log",
        err="logs/busco_output/{sample}.err"
    shell: "busco -i {input.fasta_prot} -f -m prot -l {params.lineage} --offline -c {threads} --download_path {params.busco_path} --out_path busco/output/ -o {wildcards.sample} > {log.out} 2> {log.err}"

rule extract_busco_match:
    input:
        done = rules.busco_output.output.done
    output:
        txt = "gene_info_table/busco_match_{sample}.txt"
    params:
         lineage=config["lineage"],
    shell:
        "grep -v \"^#\" busco/input/{wildcards.sample}/run_{params.lineage}/full_table.tsv | cut -f3 | sort -u | grep -v \"^$\"  > {output.txt}"


###This is some kind of implementation of Laila script.
rule ISS_table_creation_and_purification:
    input:
        ISS_full=rules.merge_iss.output.merge_iss,
        cnv_interpro=rules.run_cnv_interpro.output.ipro,
        fasta_index=rules.fasta_index.output.index,
        merged_hmm=rules.merge_hmm_results.output.merged_hmm,
        busco_match=rules.extract_busco_match.output.txt
        # hmm_results="output_hmm_search/{sample}"
    output:
        ISS_raw="gene_info_table/gene_info_table_raw_{sample}.tsv",
        ISS_final="gene_info_table/gene_info_table_purified_{sample}.tsv",
        ISS_final_removed="gene_info_table/gene_info_table_removed_{sample}.tsv",
        filter_summary="gene_info_table/filter_summary_{sample}.tsv",
        genes_to_remove="gene_info_table/{sample}_to_remove.tsv"
    params:
        ipro_term_to_remove="|".join(config["ipro_term_to_remove"]),
        blast_keyword_to_remove="|".join(config["keywork_to_remove"])
    log:
        out="logs/ISS_table_creation_and_purification/{sample}.log",
        err="logs/ISS_table_creation_and_purification/{sample}.err"
    run:
        # ipdb.set_trace()
        raw_iss_table = pd.read_table(input.ISS_full, index_col=0,  names=["gene_name", "subject_prot", "results_prot",
                                                                           "length_prot", "annotation_prot",  "ISS_prot",
                                                                           "completness_prot", "gene_prot", "dbxref_prot",
                                                                           "subject_TE", "results_TE", "length_TE",
                                                                           "annotation_TE",  "ISS_TE", "completness_TE",
                                                                           "gene_TE", "dbxref_TE"])
        # raw_iss_table = pd.read_table(input.ISS_full, index_col=0,  names=["gene_name", "annotation_prot", "ISS_prot", "completness_prot", "gene_prot", "other_prot", "annotation_TE", "ISS_TE", "completness_TE", "gene_TE", "other_TE"])
        ###Add length informations.
        protein_length = pd.read_table(
            input.fasta_index,index_col=0,usecols=[0, 1],low_memory=True,names=["gene_name", "protein_length",
                                                                                      "usless_1", "usless_2",
                                                                                      "usless_3"])
        raw_iss_table = raw_iss_table.join(protein_length,on="gene_name")

        ###Add interpro informations.
        cnv_interpro_table = pd.read_table(input.cnv_interpro,names=["gene_name", "ipro_id"])  ###Read file
        cnv_interpro_table = cnv_interpro_table[cnv_interpro_table["ipro_id"] != "-"]  ###Remove useless empty rows.
        cnv_interpro_table = cnv_interpro_table.groupby('gene_name')['ipro_id'].agg([('merged_interpro', ';'.join)]) ###Merge lines with same gene id
        # cnv_interpro_table_merge = cnv_interpro_table.groupby('gene_name')['merged_interpro'].agg([('merged_interpro', ';'.join)])

        cnv_interpro_table_merge = pd.merge(left=raw_iss_table,right=cnv_interpro_table,on="gene_name",how='outer')

        ###Add hmm evalue informations
        if os.path.getsize(input.merged_hmm):
            hmm_table = pd.read_table(input.merged_hmm, index_col=0)
            # ipdb.set_trace()
            cnv_interpro_table_merge = cnv_interpro_table_merge.join(hmm_table, how="outer")

        ###Calculation of which proteins to remove
        #Filter 1 : remove all ISS_6
        removed_ISS6 = cnv_interpro_table_merge["ISS_prot"] == "ISS_6"
        #Filter 2 : Remove ISS_5 with a length < 150AA
        removed_ISS5_small = (cnv_interpro_table_merge["ISS_prot"] == "ISS_5") & (
                    cnv_interpro_table_merge["protein_length"] < 150)
        #Filter 3 : Remove if ISS_TE is greater than ISS not TE (TE 2 3 4) & (NTE 4 5)
        removed_ISS4_5_below_TE = (cnv_interpro_table_merge["ISS_TE"].isin(["ISS_2", "ISS_3", "ISS_4"])) & (
            cnv_interpro_table_merge["ISS_prot"].isin(["ISS_5", "ISS4"]))
        #Filter 4 : remove fragmented ISS_3 if the corresponding ISS_TE is 2 or 3.
        removed_incomplete_ISS3_below_TE = (cnv_interpro_table_merge["ISS_TE"].isin(["ISS_2", "ISS_3"])) & (
                    (cnv_interpro_table_merge["ISS_prot"] == "ISS_3") & (
                        cnv_interpro_table_merge["completness_prot"] == 'fragment'))
        #Filter 5 : Remove a list of keyword associated with unwanted stuff.
        removed_TE_keyword = cnv_interpro_table_merge[
            "annotation_prot"].str.contains(params.blast_keyword_to_remove,case=False,na=False)
        #Filter 6 : remove prots containing some interpro terms if they are not "modules".
        remove_TE_interpro = (cnv_interpro_table_merge[
                                  "merged_interpro"].str.contains(params.ipro_term_to_remove,case=False,na=False)) & (
                                 ~cnv_interpro_table_merge[
                                     'completness_prot'].str.contains('modules',case=False,na=False))

        #Fitler 7 : Remove IPR018289 when not associated with IPR031052 and module.
        remove_MULE = (cnv_interpro_table_merge["merged_interpro"].str.contains("IPR018289",case=False,na=False)) &(
                                ~cnv_interpro_table_merge["merged_interpro"].str.contains("IPR031052",case=False,na=False)) &(
                                ~cnv_interpro_table_merge['completness_prot'].str.contains('modules',case=False,na=False))

        #Filter 8 : Remove IPR000477 when not associated with IPR024937 or IPR003545 and not module.
        remove_RT = (cnv_interpro_table_merge["merged_interpro"].str.contains("IPR000477",case=False,na=False)) &(
                                ~cnv_interpro_table_merge["merged_interpro"].str.contains("IPR024937|IPR003545",case=False,na=False)) &(
                                ~cnv_interpro_table_merge['completness_prot'].str.contains('modules',case=False,na=False))



        ###Add a col with True when a busco is found
        busco_hit = pd.read_table(input.busco_match,header=None,index_col=0)
        busco_hit["is_busco"] = True
        cnv_interpro_table_merge = cnv_interpro_table_merge.join(busco_hit,how="outer")
        cnv_interpro_table_merge["is_busco"] = cnv_interpro_table_merge["is_busco"].fillna(False)
        busco_to_keep = cnv_interpro_table_merge["is_busco"]
        cnv_interpro_table_merge.to_csv(output.ISS_raw, sep="\t")

        # ipdb.set_trace()



        if os.path.getsize(input.merged_hmm):
            any_hmm_pvalue = ~(cnv_interpro_table_merge.filter(regex="evalue_.*").isnull()).any(axis=1)
            to_remove = (removed_ISS6 | removed_ISS5_small | removed_ISS4_5_below_TE | removed_incomplete_ISS3_below_TE | removed_TE_keyword | remove_TE_interpro | remove_MULE | remove_RT | any_hmm_pvalue)  & ~busco_to_keep
            filter_summary = pd.concat([removed_ISS6, removed_ISS5_small, removed_ISS4_5_below_TE, removed_incomplete_ISS3_below_TE, removed_TE_keyword, remove_TE_interpro, remove_MULE, remove_RT, busco_to_keep],axis=1)
            filter_summary.columns=["ISS 6", "Small ISS 5", "ISS 4-5 and TE", "Icomplete ISS3 TE", "TE Keyword", "TE Interpro", "MULE", "RT", "is_busco"]
            filter_summary = pd.concat([filter_summary,   ~(cnv_interpro_table_merge.filter(regex="evalue_.*").isnull())], axis=1)
            filter_summary.loc['Column_Total'] = filter_summary.sum(numeric_only=True,axis=0)
        else:
            to_remove = (removed_ISS6 | removed_ISS5_small | removed_ISS4_5_below_TE | removed_incomplete_ISS3_below_TE | removed_TE_keyword | remove_TE_interpro | remove_MULE | remove_RT)  & ~busco_to_keep
            filter_summary = pd.concat([removed_ISS6, removed_ISS5_small, removed_ISS4_5_below_TE, removed_incomplete_ISS3_below_TE, removed_TE_keyword, remove_TE_interpro, remove_MULE, remove_RT, busco_to_keep], axis=1)
            filter_summary.columns = ["ISS 6", "Small ISS 5", "ISS 4-5 and TE", "Icomplete ISS3 TE", "TE Keyword", "TE Interpro", "MULE", "RT", "is_busco"]
            filter_summary.loc['Column_Total'] = filter_summary.sum(numeric_only=True,axis=0)

        # ipdb.set_trace()


        # ipdb.set_trace()
        # busco_hit =  pd.read_table("./busco_generic", header=None, index_col=0)
        # busco_hit["is_busco"] = True

        # filter_summary = pd.concat([removed_ISS6, removed_ISS5_small, removed_ISS4_5_below_TE, removed_incomplete_ISS3_below_TE | removed_TE_keyword | remove_TE_interpro],axis=1)

        # to_keep = removed_ISS6 | ~removed_ISS5_small | removed_ISS4_5_below_TE | removed_incomplete_ISS3_below_TE | removed_TE_keyword | remove_dunno_interpro | remove_TE_interpro
        cnv_interpro_table_merge[to_remove].to_csv(output.ISS_final_removed, sep="\t")
        cnv_interpro_table_merge = cnv_interpro_table_merge[~to_remove]
        cnv_interpro_table_merge.to_csv(output.ISS_final, sep="\t")
        filter_summary.to_csv(output.filter_summary, sep="\t")

        to_remove[to_remove].to_csv(output.genes_to_remove,columns=[],header=False)

        # cnv_interpro_table_merge = cnv_interpro_table_merge[~(
        #             removed_ISS6 | ~removed_ISS5_small | removed_ISS4_5_below_TE | removed_incomplete_ISS3_below_TE | removed_TE_keyword)]

        # ###Now merge with all the hmm output stuff
        # for i in hmm_databases:
        #     cnv_interpro_table_merge ###FUSION




if only_reannotate:
    rule report_only_annotate:
        input:
            fasta_index=rules.fasta_index.output.index,
            ISS_protein="blast/ISS_informations/results_{sample}.iss",
        output:
            "pipeline_report_{sample}.html"
        envmodules:
            config["modules"]["R"]
        script:
            "pipeline_report.Rmd"
else:
    rule report_reannotation:
        input:
            fasta_index= rules.fasta_index.output.index,
            ISS_protein="blast/ISS_informations/results_{sample}.iss",
            ISS_TE="blast/ISS_informations/results_{sample}_TE.iss",
            processed_gff= "enriched_gff3/{sample}/{sample}_mrna.gff3",
            filter_summary="gene_info_table/filter_summary_{sample}.tsv",
            proteins_index="enriched_gff3/{sample}/{sample}_proteins.fasta.fai",
            busco_in="busco/input/{sample}/done",
            busco_out="busco/output/{sample}/done"
        output:
            "pipeline_report_{sample}.html"
        envmodules:
            config["modules"]["R"]
        params:
            to_remove_string="gene_info_table/{sample}_to_remove.tsv" if only_reannotate else "",
            filter_summary="gene_info_table/filter_summary_{sample}.tsv" if only_reannotate else ""
        script:
            "pipeline_report.Rmd"

