#!/bin/bash -l

#SBATCH -J phyluce_Aptost_all
#SBATCH -D /home/lgnewton/phyluce_Aptost_all
#SBATCH -o /home/lgnewton/phyluce_Aptost_all/out-%j.txt
#SBATCH -e /home/lgnewton/phyluce_Aptost_all/err-%j.txt
#SBATCH -t 10:00:00
#SBATCH -N 16
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1
#SBATCH -p high2
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lgnewton@ucdavis.edu # Email to which notifications will be sent

set -e
set -u

#### load modules ####
module load phyluce

#### activate phyluce environment ####
conda activate phyluce

#### spades and match contigs to probes has already been run ####

# get_match_counts
#phyluce_assembly_get_match_counts \
#    --locus-db /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/probe_match_UCEs_min65/probe.matches.sqlite \
#    --taxon-list-config /home/lgnewton/phyluce_Aptost_all/datasets.conf \
#    --taxon-group 'Aptost_Icenoglei' \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/Aptost_Icenoglei.conf \
#    --incomplete-matrix \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files

#get some snps

#get_fastas_from_match_counts
#phyluce_assembly_get_fastas_from_match_counts \
#     --contigs /home/lgnewton/Aptost_scaffolds_all \
#     --locus-db /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/probe_match_UCEs_min65/probe.matches.sqlite \
#     --match-count-output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/Aptost_Icenoglei.conf \
#     --incomplete-matrix /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/Aptost_Icenoglei.incomplete \
#     --output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/Aptost_Icenoglei.fasta \
#     --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log
     
#align with ends trimmed
#phyluce_align_seqcap_align \
#    --fasta /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/Aptost_Icenoglei.fasta \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/mafft-nexus-edge-trimmed/ \
#    --taxa 39 \
#    --aligner mafft \
#    --cores 2 \
#    --incomplete-matrix \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files

#phyluce_align_explode_alignments \
#   --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/mafft-nexus-edge-trimmed/ \
#	--input-format nexus \
#	--output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/mafft-nexus-edge-trimmed-exploded \
#	--by-taxon

##map reads against contigs
#phyluce_snp_bwa_multiple_align \
#    --config /home/lgnewton/AptoIce_phasing.conf \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/multialign-bams \
#    --cores 8 \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files \
#    --mem

#phase
#phyluce_snp_phase_uces \
#	--config /home/lgnewton/AptoIce_phasing.conf \
#	--bams /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/multialign-bams \
#	--output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/multialign-bams-phased-reads \
#	--log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files

#align phased data; taxa is now doubled
#phyluce_align_seqcap_align \
#	--fasta /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/multialign-bams-phased-reads/fastas/joined_allele_sequences_all_samples.fasta \
#	--output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed \
#	--taxa 78 \
#	--aligner mafft \
#	--cores 4 \
#	--incomplete-matrix \
#	--log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files

#phyluce_align_get_align_summary_data \
#     --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed \
#     --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files \

#phyluce_align_remove_locus_name_from_nexus_lines \
#    --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files \
#    --cores 4

#phyluce_align_get_only_loci_with_min_taxa \
#    --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean \
#    --taxa 78 \
#    --percent 0.95 \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean-min-95p-taxa \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files \
#    --cores 2

#phyluce_snp_screen_phased_alignments \
#	--alignments /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean-min-95p-taxa \
#	--output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean-min-95p-taxa-randSNP1 \
#	--input-format nexus \
#	--output-format nexus \
#	--random \
#	--cores 4 \
#	--log-path /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/log_files

#phyluce_align_format_nexus_files_for_raxml \
#	--alignments /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean-min-95p-taxa-randSNP1 \
#	--output /home/lgnewton/phyluce_Aptost_all/Aptost_Icenoglei/PHASED-DATA-mafft-nexus-edge-trimmed-clean-min-95p-taxa-randSNP1-concat-char \
#	--charsets


