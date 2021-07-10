#
# Copyright 2021 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

args = commandArgs(trailingOnly=TRUE)

if (args[1] == "-h" | args[1] == "--help") {
  cat("", sep = "\n")
  cat(paste0("Usage: Rscript Bionano_fixdups.R <scaffolds> <agp> <contigs>"), sep = "\n")
  cat(paste0("<scaffolds>: fasta file with scaffolds produced by Bionano hybrid scaffolding pipeline"), sep = "\n")
  cat(paste0("<agp>: agp file describing which contig has been included in each scaffold"), sep = "\n")
  cat(paste0("<contigs>: fasta file with contigs cut by Bionano hybrid scaffolding pipeline"), sep = "\n")
  cat(paste0("Note: before running, edit the script to set variables SAMEXTRACTCLIPS, MINIMAP2, SAMTOOLS and SEQTK to the corresponding executables"), sep = "\n")
  stop(simpleError(sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))))
}

scaffolds_file_name <- args[1]
agp_file_name <- args[2]
contigs_file_name <- args[3]

###################################################################################
SAMEXTRACTCLIPS="java -jar samextractclip.jar"
MINIMAP2="minimap2"
SAMTOOLS="samtools"
SEQTK="seqtk"
#length of gaps placed by Bionano software to identify negative gaps (contigs overlaps)
negative_gap_placeholder_length <- 13
fasta_width <- 80 #should be <20000
size_empty_bam <- 500
num_threads <- 25
###################################################################################

suppressMessages(library("Biostrings"))
working_directory <- dirname(contigs_file_name)
fixed_scaffolds_file_name <- paste0(working_directory, "/", sub(pattern = "\\.fasta", replacement = "_neg_gaps_fixed.fasta", x = basename(scaffolds_file_name)))
temp_dir <- paste0(working_directory, "/fix_scaffolding_temp")
contigs_unoriented <- readDNAStringSet(contigs_file_name, "fasta")
scaffolds <- readDNAStringSet(scaffolds_file_name, "fasta")
agp <- read.table(file = agp_file_name, stringsAsFactors = FALSE, sep = "\t")
scaffolds_names <- names(scaffolds)

#note: contigs are defined as first or second based on their relative position in the scaffold

if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}
logfile <- paste0(working_directory, "/logfile_fix_scaffolds.txt")

#cycle over the scaffolds
for (i in 1:length(scaffolds_names)) {
  merged_contigs_curr <- c()
  scaffold_name_curr <- scaffolds_names[i]
  contigs_ids_scaffold_curr_file_name <- paste0(working_directory, "/contigs_ids_scaffold_", scaffold_name_curr)
  contigs_scaffold_curr_file_name <- paste0(working_directory, "/contigs_scaffold_", scaffold_name_curr, ".fasta")
  scaffold_id_curr_file_name <- paste0(working_directory, "/scaffold_id_", scaffold_name_curr)
  scaffold_curr_file_name <- paste0(working_directory, "/scaffold_", scaffold_name_curr, ".fasta")
  alignment_contigs_scaffold_sam_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_contigs_mapped.sam")
  alignment_contigs_scaffold_bam_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_contigs_mapped.bam")
  contigs_oriented_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_contigs_oriented.fasta")
  agp_scaffold_curr <- agp[which(agp[, 1]==scaffold_name_curr), ]
  contigs_names_curr <- agp_scaffold_curr[ , 6]
  ind_negative_gaps_curr <- intersect(which(agp_scaffold_curr[, 5]=="N"), which(agp_scaffold_curr[, 6]==negative_gap_placeholder_length))
  ind_positive_gaps_curr <- intersect(which(agp_scaffold_curr[, 5]=="N"), which(agp_scaffold_curr[, 6]!=negative_gap_placeholder_length))
  ind_gaps_curr <- sort(c(ind_negative_gaps_curr, ind_positive_gaps_curr))
  positive_gaps_size <- as.double(agp_scaffold_curr[ind_positive_gaps_curr, 6])
  #if there are not negative gaps, copy the scaffold as it is
  if (length(ind_negative_gaps_curr) == 0) {
    cat(sprintf("No negative overlaps for scaffold %s\n", scaffold_name_curr), file = logfile, append = TRUE, sep = "")
    writeXStringSet(x = scaffolds[i], filepath = fixed_scaffolds_file_name, append = TRUE, width = fasta_width)
    #else, fix the overlaps
  } else {
    cat(sprintf("Solving overlaps for scaffold %s\n", scaffold_name_curr), file = logfile, append = TRUE, sep = "")
    #retrieve the correct orientation of the contigs in the scaffold
    contigs_names_curr_tmp <- paste(contigs_names_curr[-ind_gaps_curr], sep = "\n")
    cat(contigs_names_curr_tmp, file = contigs_ids_scaffold_curr_file_name, sep = "\n")
    system(paste0(SEQTK, " subseq ", contigs_file_name, " ", contigs_ids_scaffold_curr_file_name, " > ", contigs_scaffold_curr_file_name))
    system(paste0("echo ", scaffold_name_curr, " > ", scaffold_id_curr_file_name))
    system(paste0(SEQTK, " subseq ", scaffolds_file_name, " ", scaffold_id_curr_file_name, " > ", scaffold_curr_file_name))
    system(paste0(MINIMAP2, " -t ", num_threads, " -a ", scaffolds_file_name, " ", contigs_scaffold_curr_file_name, " -o ", alignment_contigs_scaffold_sam_file_name))
    system(paste0(SAMTOOLS, " view -hSb -q20 -F2304 ", alignment_contigs_scaffold_sam_file_name, " | ", SAMTOOLS, " sort > ", alignment_contigs_scaffold_bam_file_name))
    system(paste0(SAMTOOLS, " view ", alignment_contigs_scaffold_bam_file_name, " | awk '{print \">\"$1\"\\n\"$10}' > ", contigs_oriented_file_name))
    contigs <- readDNAStringSet(contigs_oriented_file_name, "fasta")
    system(paste("mv", contigs_ids_scaffold_curr_file_name, contigs_scaffold_curr_file_name, scaffold_id_curr_file_name, scaffold_curr_file_name, alignment_contigs_scaffold_sam_file_name, alignment_contigs_scaffold_bam_file_name, contigs_oriented_file_name, temp_dir))
    #cycle over the gaps in the current scaffold
    for (k in 1:length(ind_gaps_curr)) {
      #if the current gap is a positive (standard) one
      if (ind_gaps_curr[k] %in% ind_positive_gaps_curr) {
        #if it is the gap between the first two contigs, add the gap to the first contig
        if (k == 1) {
          first_contig_name <- contigs_names_curr[(ind_gaps_curr[k] - 1)]
          first_contig <- contigs[which(names(contigs) == first_contig_name)]
          second_contig_name <- contigs_names_curr[(ind_gaps_curr[k] + 1)]
          second_contig <- contigs[which(names(contigs) == second_contig_name)]
          gap_size_curr <- as.double(agp_scaffold_curr[ind_gaps_curr[k], 6])
          gap_seq_curr <- DNAStringSet(paste0(rep(x = "N", times = gap_size_curr), collapse = ""))
          merged_contigs_curr <- DNAStringSet(Map(c, first_contig, gap_seq_curr, second_contig))
          #otherwise, add the gap to the previously merged contig
        } else {
          first_contig <- merged_contigs_curr
          second_contig_name <- contigs_names_curr[(ind_gaps_curr[k] + 1)]
          second_contig <- contigs[which(names(contigs) == second_contig_name)]
          gap_size_curr <- as.double(agp_scaffold_curr[ind_gaps_curr[k], 6])
          gap_seq_curr <- DNAStringSet(paste0(rep(x = "N", times = gap_size_curr), collapse = ""))
          merged_contigs_curr <- DNAStringSet(Map(c, first_contig, gap_seq_curr, second_contig))
        }
        #else, if the current gap is a negative one ( -> 2 contigs are overlapping)  
      } else {
        reference_contig_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_reference_solving_gap_", k, "_tmp.fasta")
        query_contig_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_query_solving_gap_", k, "_tmp.fasta")
        #if the negative gap is between the first 2 contigs, then pick the first contig as the first_overlapping_contig
        if (k == 1) {
          first_overlapping_contig_name <- contigs_names_curr[(ind_gaps_curr[k] - 1)]
          first_overlapping_contig <- contigs[which(names(contigs) == first_overlapping_contig_name)]
          first_overlapping_contig_length <- width(first_overlapping_contig)
          #if the negative gap is not between the first two contigs, then pick the previously merged contig as the first_overlapping_contig
        } else {
          first_overlapping_contig_name <- names(merged_contigs_curr)
          first_overlapping_contig <- merged_contigs_curr
          first_overlapping_contig_length <- width(first_overlapping_contig)
        }
        #pick the second overlapping contig and evaluate which is the longest of the two
        second_overlapping_contig_name <- contigs_names_curr[(ind_gaps_curr[k] + 1)]
        second_overlapping_contig <- contigs[which(names(contigs) == second_overlapping_contig_name)]
        second_overlapping_contig_length <- width(second_overlapping_contig)
        ind_longest_contig <- which(c(first_overlapping_contig_length, second_overlapping_contig_length)==max(c(first_overlapping_contig_length, second_overlapping_contig_length)))
        #if the first contig is the longest, consider it as the reference; map to it the query (second contig);
        #extract the merged contig as a concatenation of the reference followed by the 3' end soft-clipped portion of the query (if there is any)
        if (ind_longest_contig == 1) {
          reference_contig <- first_overlapping_contig
          writeXStringSet(x = reference_contig, filepath = reference_contig_file_name, append = FALSE, width = fasta_width)
          query_contig <- second_overlapping_contig
          writeXStringSet(x = query_contig, filepath = query_contig_file_name, append = FALSE, width = fasta_width)
          alignment_sam_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_solving_gap_num", k, ".sam")
          alignment_bam_file_name <-  paste0(working_directory, "/", scaffold_name_curr, "_solving_gap_num", k, ".bam")
          query_sc_file_name <-  paste0(working_directory, "/", scaffold_name_curr, "_softclipped_solving_gap_", k, ".fasta")
          query_sc_file_name_tmp <-  paste0(working_directory, "/", scaffold_name_curr, "_softclipped_solving_gap_", k, "_tmp.fasta")
          system(paste0(MINIMAP2, " -t ", num_threads, " -a ", reference_contig_file_name, " ", query_contig_file_name, " -o ", alignment_sam_file_name))
          system(paste0(SAMTOOLS, " view -hSb -q20 -F2304 ", alignment_sam_file_name, " | ", SAMTOOLS, " sort > ", alignment_bam_file_name))
          system(paste0(SAMEXTRACTCLIPS, " ", alignment_bam_file_name, " | ", SEQTK, " seq -A > ", query_sc_file_name_tmp))
          soft_clipped_file_size <- file.info(query_sc_file_name_tmp)$size
          alignment_bam_file_size <- file.info(alignment_bam_file_name)$size
          #if the query is not soft-clipped then it should be fully aligned to the reference (or not aligned at all but Bionano maps say they are)
          if (soft_clipped_file_size == 0 || alignment_bam_file_size < size_empty_bam) {
            #if it is fully aligned to the reference, do nothing
            if (soft_clipped_file_size == 0 && alignment_bam_file_size > size_empty_bam) {
              cat(sprintf("WARNING: Alignment without soft-clipping between contigs %s and %s from scaffold %s (gap %s): not concatenating them\n", names(reference_contig), names(query_contig), scaffold_name_curr, k), file = logfile, append = TRUE, sep = "")
              merged_contigs_curr <- reference_contig
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
              #if it is not aligned at all, concatenate them
            } else if (alignment_bam_file_size < size_empty_bam) {
              cat(sprintf("WARNING: No alignment between contigs %s and %s from scaffold %s: concatenating them\n", names(reference_contig), names(query_contig), scaffold_name_curr), file = logfile, append = TRUE, sep = "")
              merged_contigs_curr <- DNAStringSet(Map(c, reference_contig, query_contig))
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
            }
            #if the length of the soft-clipped file is non zero
          } else {
            sc_contig_full <- readDNAStringSet(query_sc_file_name_tmp, "fasta")
            #if there is a soft-clipped region both at the 5' and at the 3' pick the one at the 3'
            if (length(sc_contig_full) > 1) {
              sc_contig <- sc_contig_full[2]
              names(sc_contig) <- substr(x = names(sc_contig), start = 1, stop = 80)
              writeXStringSet(x = sc_contig, filepath = query_sc_file_name, append = FALSE, width = fasta_width)
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
              merged_contigs_curr <- DNAStringSet(Map(c, reference_contig, sc_contig))
              #otherwise pick the only one available (if it is at the right end)
            } else {
              #check if the 3' end of the query contig is the soft-clipped portion
              end_sc <- as.double(substr(x = names(sc_contig_full), start = (nchar(names(sc_contig_full)) - 1), stop = (nchar(names(sc_contig_full)) - 1)))
              if (end_sc == 5) {
                cat(sprintf("WARNING: Alignment without soft-clipping between contigs %s and %s from scaffold %s (gap %s): not concatenating them\n", names(query_contig), names(reference_contig), scaffold_name_curr, k), file = logfile, append = TRUE, sep = "")
                merged_contigs_curr <- reference_contig
              }
              sc_contig <- sc_contig_full
              names(sc_contig) <- substr(x = names(sc_contig), start = 1, stop = 80)
              writeXStringSet(x = sc_contig, filepath = query_sc_file_name, append = FALSE, width = fasta_width)
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
              merged_contigs_curr <- DNAStringSet(Map(c, reference_contig, sc_contig))
            }
          }
          #if the second contig is the longest, consider it as the reference; map to it the query (first contig);
          #extract the merged contig as a concatenation of the the 5' end soft-clipped portion of the query (if there is any) followed by the reference
        } else {
          reference_contig <- second_overlapping_contig
          writeXStringSet(x = reference_contig, filepath = reference_contig_file_name, append = FALSE, width = fasta_width)
          query_contig <- first_overlapping_contig
          writeXStringSet(x = query_contig, filepath = query_contig_file_name, append = FALSE, width = fasta_width)
          alignment_sam_file_name <- paste0(working_directory, "/", scaffold_name_curr, "_solving_gap_num", k, ".sam")
          alignment_bam_file_name <-  paste0(working_directory, "/", scaffold_name_curr, "_solving_gap_num", k, ".bam")
          query_sc_file_name <-  paste0(working_directory, "/", scaffold_name_curr, "_softclipped_solving_gap_", k, ".fasta")
          query_sc_file_name_tmp <-  paste0(working_directory, "/", scaffold_name_curr, "_softclipped_solving_gap_", k, "_tmp.fasta")
          system(paste0(MINIMAP2, " -t ", num_threads," -a ", reference_contig_file_name, " ", query_contig_file_name, " -o ", alignment_sam_file_name))
          system(paste0(SAMTOOLS, " view -hSb -q20 -F2304 ", alignment_sam_file_name, " | ", SAMTOOLS, " sort > ", alignment_bam_file_name))
          system(paste0(SAMEXTRACTCLIPS, " ", alignment_bam_file_name, " | ", SEQTK, " seq -A > ", query_sc_file_name_tmp))
          soft_clipped_file_size <- file.info(query_sc_file_name_tmp)$size
          alignment_bam_file_size <- file.info(alignment_bam_file_name)$size
          #if the query is not soft-clipped then it should be fully aligned to the reference (or not aligned at all but Bionano maps say they are)
          if (soft_clipped_file_size == 0 || alignment_bam_file_size < size_empty_bam) {
            #if it is fully aligned to the reference, do nothing
            if (soft_clipped_file_size == 0 && alignment_bam_file_size > size_empty_bam) {
              cat(sprintf("WARNING: Alignment without soft-clipping between contigs %s and %s from scaffold %s (gap %s): not concatenating them\n", names(reference_contig), names(query_contig), scaffold_name_curr, k), file = logfile, append = TRUE, sep = "")
              merged_contigs_curr <- reference_contig
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
              #if it is not aligned at all, concatenate them
            } else if (alignment_bam_file_size < size_empty_bam) {
              cat(sprintf("WARNING: No alignment between contigs %s and %s from scaffold %s (gap %s): concatenating them\n", names(reference_contig), names(query_contig), scaffold_name_curr, k), file = logfile, append = TRUE, sep = "")
              merged_contigs_curr <- DNAStringSet(Map(c, query_contig, reference_contig))
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
            }
            #if the length of the soft-clipped file is non zero
          } else {
            sc_contig_full <- readDNAStringSet(query_sc_file_name_tmp, "fasta")
            #if there is a soft-clipped region both at the 5' and at the 3' pick the one at the 5'
            if (length(sc_contig_full) > 1) {
              sc_contig <- sc_contig_full[1]
              names(sc_contig) <- substr(x = names(sc_contig), start = 1, stop = 80)
              writeXStringSet(x = sc_contig, filepath = query_sc_file_name, append = FALSE, width = fasta_width)
              system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
              merged_contigs_curr <- DNAStringSet(Map(c, sc_contig, reference_contig))
              #otherwise pick the only one available (if it is at the right end)
            } else {
              #check if the 5' end of the query contig is the soft-clipped portion
              end_sc <- as.double(substr(x = names(sc_contig_full), start = (nchar(names(sc_contig_full)) - 1), stop = (nchar(names(sc_contig_full)) - 1)))
              if (end_sc == 3) {
                cat(sprintf("WARNING: Alignment without soft-clipping between contigs %s and %s contig from scaffold %s (gap %s): not concatenating them\n", names(query_contig), names(reference_contig), scaffold_name_curr, k), file = logfile, append = TRUE, sep = "")
                merged_contigs_curr <- reference_contig
              } else {
                sc_contig <- sc_contig_full
                names(sc_contig) <- substr(x = names(sc_contig), start = 1, stop = 80)
                writeXStringSet(x = sc_contig, filepath = query_sc_file_name, append = FALSE, width = fasta_width)
                system(paste("mv", reference_contig_file_name, query_contig_file_name, alignment_bam_file_name, alignment_sam_file_name, query_sc_file_name, query_sc_file_name_tmp, temp_dir))
                merged_contigs_curr <- DNAStringSet(Map(c, sc_contig, reference_contig))
              }
            }
          }
        }
      }
    }
    sam_files <- list.files(path = temp_dir, pattern = "\\.sam", full = TRUE)
    system(paste0("rm ", paste(sam_files, collapse = " ")))
    new_scaffold_name <- paste0(names(scaffolds[i]), "_fixed", collapse = "")
    names(merged_contigs_curr) <- new_scaffold_name
    writeXStringSet(x = merged_contigs_curr, filepath = fixed_scaffolds_file_name, append = TRUE, width = fasta_width)
  }
}
