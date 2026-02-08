// -----------------------------------------------------------------------------
// File Name:             flashGenomeAligner.cpp
// Author:                Diego Arcos Sapena
// Version:               1.0
// Copyright:             Copyright (c) 2025 Diego Arcos Sapena
//
// Description:
//   This C++ program implements an ultra-efficient, parallel algorithm for
//   rapid gene-by-gene alignment and comprehensive analysis of mutations and
//   genomic variability in sequences, meticulously designed and optimized
//   for maximum speed using SIMD instructions from Parasail library.
//
// Project:
//   Final Degree Project (TFG) in Biotechnology - ANÁLISIS DE MUTACIONES Y VARIABILIDAD GENÓMICA EN LOS VIRUS DEL ZIKA Y DENGUE MEDIANTE EL DESARROLLO DE UN ALGORITMO BIOINFORMÁTICO PARALELO PARA EL ALINEAMIENTO Y ANÁLISIS EFICIENTES DE SECUENCIAS
//   Author: Diego Arcos Sapena
//   Carried out for: Cheminformatics and Nutrition Research Group
//   Supervisor: Santiago Garcia Vallvé
//   University: Universitat Rovira i Virgili (URV)
//   Location: Tarragona, Spain
//   Date: 2025
//
// Developer contact:
//   diegoarcos33@gmail.com
//   Diegoelpercu on Github
//
// Research group contact:
//   URV-cheminformatics on Github
//
//
// -----------------------------------------------------------------------------
// License:
//   This project is licensed under the Apache License 2.0.

#include <parasail.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <omp.h>
#include <atomic>
#include <iomanip>
#include <map>
#include <cstring>
#include <algorithm>
#include <set>

const std::string valid_alphabet_global = "ACGTURYKMSWBDHVN";

// Definition of the standard genetic code (DNA to Amino Acid)
const std::map<std::string, char> genetic_code = {
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'}, // STOP
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'}, // STOP, Tryptophan
    {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'}, // Methionine
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
    {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
};

struct NtChangeDetail {
    char ref_base;
    int ref_nt_pos_gene;    // Nucleotide position relative to the gene
    char query_base;
    int ref_nt_pos_genome;  // Genomic position of ref_base
};

// --- Helper functions for TSV formatting ---
/**
* @brief Joins a vector of strings into a single string, separating elements with a delimiter.
* Performs basic escaping of problematic characters for TSV (tab, newline, double quotes).
* @param vec Vector of strings to join.
* @param delimiter Delimiter to use between elements.
* @return String resulting from the union.
*/
std::string join_strings(const std::vector<std::string>& vec, const std::string& delimiter) {
    if (vec.empty()) {
        return "";
    }
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) {
            oss << delimiter;
        }
        std::string item = vec[i];

        std::replace(item.begin(), item.end(), '\t', ' '); // Replace tabs with spaces
        std::replace(item.begin(), item.end(), '\n', ' '); // Replace newlines with spaces
        std::replace(item.begin(), item.end(), '"', '\''); // Replace double quotes with single ones (or you could double them if the TSV software prefers)
        oss << item;
    }
    return oss.str();
}

/**
* @brief Formats a base count map (char -> int) into a string for TSV.
* Example: "N:1;Y:2"
* @param base_counts Base count map.
* @return Formatted string.
*/
std::string format_base_counts_map(const std::map<char, int>& base_counts) {
    std::vector<std::string> parts;
    const std::string ambiguous_order = "NURYKMSWBDHV";
    for (char c : ambiguous_order) {
        auto it = base_counts.find(c);
        if (it != base_counts.end() && it->second > 0) {
            parts.push_back(std::string(1, c) + ":" + std::to_string(it->second));
        }
    }
    // Add other characters if for some reason they appear and are not in ambiguous_order
    for (const auto& pair : base_counts) {
        if (pair.second > 0 && ambiguous_order.find(pair.first) == std::string::npos) {
            parts.push_back(std::string(1, pair.first) + ":" + std::to_string(pair.second));
        }
    }
    return join_strings(parts, ";"); // Use semicolon as internal separator
}


// Stream processing functions
char translate_codon(const std::string& codon_str,
    const std::map<std::string, char>& code,
    bool& is_valid_for_translation) {
    is_valid_for_translation = false;

    if (codon_str.length() != 3) {
        return '?';
    }

    std::string processed_codon = codon_str;
    for (char& base : processed_codon) {
        base = std::toupper(static_cast<unsigned char>(base));
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            return '?';
        }
    }

    auto it = code.find(processed_codon);
    if (it != code.end()) {
        is_valid_for_translation = true;
        return it->second;
    }
    else {
        return '?';
    }
}

std::string read_fasta_sequence_only(const std::string& filename) {
    std::ifstream file(filename);
    std::string line_content, sequence_data;
    if (!file.is_open()) {
        std::cerr << "Error: The file could not be opened: " << filename << std::endl;
        exit(1);
    }
    while (std::getline(file, line_content)) {
        if (line_content.empty() || line_content[0] == '>') continue;
        // Remove whitespace at the end of the line (if any)
        line_content.erase(line_content.find_last_not_of(" \n\r\t") + 1);
        sequence_data += line_content;
    }
    return sequence_data;
}

std::vector<std::pair<std::string, std::string>> read_multi_fasta(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<std::string, std::string>> sequences_vector;
    std::string line_content, current_id, current_sequence;
    if (!file.is_open()) {
        std::cerr << "Error: The multi-fasta file could not be opened: " << filename << std::endl;
        exit(1);
    }
    while (std::getline(file, line_content)) {
        if (line_content.empty()) continue;
        if (line_content[0] == '>') {
            if (!current_id.empty()) { // Save the previous sequence
                sequences_vector.push_back({ current_id, current_sequence });
            }
            current_id = line_content.substr(1);
            // Remove whitespace at the beginning/end of the ID if necessary
            current_id.erase(0, current_id.find_first_not_of(" \n\r\t"));
            current_id.erase(current_id.find_last_not_of(" \n\r\t") + 1);
            current_sequence.clear();
        }
        else {
            // Remove whitespace at the end of the sequence line
            line_content.erase(line_content.find_last_not_of(" \n\r\t") + 1);
            current_sequence += line_content;
        }
    }
    if (!current_id.empty()) { // Save the last sequence
        sequences_vector.push_back({ current_id, current_sequence });
    }
    return sequences_vector;
}

/**
 * @brief Analyzes a pairwise alignment to identify mutations and their consequences.
 * * This function processes aligned query and reference (target/gene) sequences to detect
 * nucleotide mismatches, insertions, and deletions. It classifies these mutations
 * (canonical, ambiguous, invalid) and determines their impact on amino acid translation,
 * identifying missense, nonsense, and synonymous mutations, as well as frameshift events.
 * The function also calculates alignment coverage and identity, counts ambiguous/invalid
 * characters in the query, and identifies intergenic insertions occurring before the
 * main gene alignment.
 * * For UTR regions (identified by `is_utr_region` flag):
 * - Nucleotide mutations (mismatches, insertions, deletions) are reported.
 * - A counter `utr_mutation_count_out` is incremented for each such mutation.
 * - Amino acid translation and reporting are SKIPPED for UTR regions.
 * - Frameshift detection and reporting remain active for indels in UTR regions.
 *
 * @param aln_query_raw The aligned query sequence string (may include gaps '-').
 * @param aln_target_raw The aligned reference (gene/target) sequence string (may include gaps '-').
 * @param gene_genomic_start_pos_1based The 1-based genomic start position of the reference gene.
 * @param gene_len_ref_ungapped The length of the reference gene sequence without any gaps.
 * @param is_utr_region True if the current gene is a UTR region, false otherwise.
 * @param[out] coverage_out Reference to a float where the calculated coverage percentage of the reference gene by the query will be stored.
 * @param[out] identity_out Reference to a float where the calculated identity percentage of the alignment will be stored.
 * @param[out] intergenic_insertion_bases_count Reference to an integer that will store the count of bases inserted in the query before the alignment to the current gene begins.
 * @param[out] intergenic_insertion_details Reference to a vector of strings that will store details of any intergenic insertions found before the current gene.
 * @param current_gene_name The name of the current gene being analyzed.
 * @param previous_gene_name The name of the gene that was processed previously (used for context in intergenic insertion reporting, e.g., "Between PREV_GENE and CURR_GENE").
 * @param[out] ambiguous_base_counts_query_out Reference to a map that will store counts of each ambiguous IUPAC base (N, R, Y, etc.) found in the aligned query sequence.
 * @param[out] invalid_char_count_query_out Reference to an integer that will store the count of non-IUPAC characters found in the aligned query sequence.
 * @param[out] has_any_mutation_this_gene_out Reference to a boolean flag that will be set to true if any mutation (mismatch, insertion, deletion, or resulting amino acid change in coding regions) is detected within the current gene; false otherwise.
 * @param[out] distinct_insertion_events_out Reference to an integer that will store the count of distinct insertion events within the gene.
 * @param[out] distinct_deletion_events_out Reference to an integer that will store the count of distinct deletion events within the gene.
 * @param[out] distinct_mismatch_events_out Reference to an integer that will store the count of distinct mismatch events (single nucleotide substitutions) within the gene.
 * @param[out] out_mismatch_details_canonical Reference to a vector of strings that will store details of canonical mismatches (e.g., "A.123.G").
 * @param[out] out_mismatch_details_ambiguous Reference to a vector of strings that will store details of mismatches involving ambiguous IUPAC characters.
 * @param[out] out_mismatch_details_invalid Reference to a vector of strings that will store details of mismatches involving invalid (non-IUPAC) characters.
 * @param[out] out_insertion_details_canonical Reference to a vector of strings that will store details of insertions composed of canonical bases (e.g., "Ins.123.ACGT").
 * @param[out] out_insertion_details_ambiguous Reference to a vector of strings that will store details of insertions involving ambiguous IUPAC characters.
 * @param[out] out_insertion_details_invalid Reference to a vector of strings that will store details of insertions involving invalid (non-IUPAC) characters.
 * @param[out] out_deletion_details_canonical Reference to a vector of strings that will store details of deletions of canonical reference bases (e.g., "Del123", "Del123-125").
 * @param[out] out_deletion_details_ambiguous Reference to a vector of strings that will store details of deletions involving ambiguous IUPAC characters from the reference (if applicable, though typically deletions are classified based on reference which is assumed canonical).
 * @param[out] out_deletion_details_invalid Reference to a vector of strings that will store details of deletions involving invalid characters from the reference (if applicable).
 * @param[out] out_missense_aa_muts Reference to a vector of strings that will store details of missense amino acid mutations (e.g., "A123G (P41S)").
 * @param[out] out_synonymous_aa_muts Reference to a vector of strings that will store details of synonymous amino acid mutations.
 * @param[out] out_nonsense_aa_muts Reference to a vector of strings that will store details of nonsense amino acid mutations (introduction of a stop codon).
 * @param[out] out_frameshift_causing_indel_reports Reference to a vector of strings that will store concise reports of indels that cause a frameshift.
 * @param[out] out_frameshift_lifecycle_reports Reference to a vector of strings that will store detailed reports about the frameshift status (e.g., FS_started, FS_restored, FS_phase_changed).
 * @param[out] utr_mutation_count_out Reference to an integer that will store the count of nucleotide mutations within UTR regions.
 * @return int The net offset in aligned length caused by indels within the gene (total length of insertions in query minus total length of deletions in query relative to reference gene).
 * A positive value means the aligned query portion for the gene is longer than the reference gene portion due to net insertions.
 * A negative value means it's shorter due to net deletions.
 */
int analyze_alignment(
    const std::string& aln_query_raw, // Query sequence from Parasail traceback
    const std::string& aln_target_raw, // Reference (target) sequence from Parasail traceback
    int gene_genomic_start_pos_1based,
    int gene_len_ref_ungapped,
    bool is_utr_region, // True if the current gene is a UTR region
    float& coverage_out,
    float& identity_out,
    int& intergenic_insertion_bases_count,
    std::vector<std::string>& intergenic_insertion_details,
    const std::string& current_gene_name,
    const std::string& previous_gene_name,
    std::map<char, int>& ambiguous_base_counts_query_out,
    int& invalid_char_count_query_out,
    bool& has_any_mutation_this_gene_out,
    int& distinct_insertion_events_out,
    int& distinct_deletion_events_out,
    int& distinct_mismatch_events_out,
    std::vector<std::string>& out_mismatch_details_canonical,
    std::vector<std::string>& out_mismatch_details_ambiguous,
    std::vector<std::string>& out_mismatch_details_invalid,
    std::vector<std::string>& out_insertion_details_canonical,
    std::vector<std::string>& out_insertion_details_ambiguous,
    std::vector<std::string>& out_insertion_details_invalid,
    std::vector<std::string>& out_deletion_details_canonical,
    std::vector<std::string>& out_deletion_details_ambiguous,
    std::vector<std::string>& out_deletion_details_invalid,
    std::vector<std::string>& out_missense_aa_muts,
    std::vector<std::string>& out_synonymous_aa_muts,
    std::vector<std::string>& out_nonsense_aa_muts,
    std::vector<std::string>& out_frameshift_causing_indel_reports,
    std::vector<std::string>& out_frameshift_lifecycle_reports,
    int& utr_mutation_count_out // Output: count of substitutions in UTR
) {
    has_any_mutation_this_gene_out = false;
    utr_mutation_count_out = 0;
    const std::string valid_iupac_for_query = "ACGTURYKMSWBDHVN";
    const std::string canonical_bases = "ACGT";

    invalid_char_count_query_out = 0;
    ambiguous_base_counts_query_out.clear();

    int current_reading_frame_offset = 0; // This offset is now only affected by indels in non-UTR regions
    bool waiting_for_in_frame_aa_report_after_fs_restore = false;

    int aln_len = std::min(aln_query_raw.size(), aln_target_raw.size());
    int ref_pos_in_gene_1based = 1;
    int query_pos_processed_bases = 0;
    int ref_bases_aligned_count = 0;

    bool in_insertion_event = false;
    bool in_deletion_event = false;
    std::string current_insertion_sequence;
    std::string current_deletion_sequence_ref_bases;
    int indel_start_ref_pos_in_gene_1based = 0;

    int total_insertion_length_bp = 0;
    int total_deletion_length_bp = 0;

    out_mismatch_details_canonical.clear(); out_mismatch_details_ambiguous.clear(); out_mismatch_details_invalid.clear();
    out_insertion_details_canonical.clear(); out_insertion_details_ambiguous.clear(); out_insertion_details_invalid.clear();
    out_deletion_details_canonical.clear(); out_deletion_details_ambiguous.clear(); out_deletion_details_invalid.clear();
    out_missense_aa_muts.clear(); out_synonymous_aa_muts.clear(); out_nonsense_aa_muts.clear();
    out_frameshift_causing_indel_reports.clear(); out_frameshift_lifecycle_reports.clear();

    std::string current_ref_codon_bases = "";
    std::string current_query_codon_bases = "";
    int current_ref_codon_start_gene_pos_nt_1based = -1;
    std::vector<NtChangeDetail> nt_changes_in_current_codon;

    int matches_count = 0;
    distinct_mismatch_events_out = 0;
    distinct_insertion_events_out = 0;
    distinct_deletion_events_out = 0;

    int alignment_start_ref_pos_in_gene_1based = -1;
    int alignment_end_ref_pos_in_gene_1based = -1;

    auto classify_mutation_chars = [&](const std::string& chars_to_classify) -> int {
        bool has_ambiguous = false;
        for (char c : chars_to_classify) {
            if (c == '-') continue;
            if (canonical_bases.find(c) == std::string::npos) {
                if (valid_iupac_for_query.find(c) != std::string::npos) {
                    has_ambiguous = true;
                }
                else {
                    return 2;
                }
            }
        }
        return has_ambiguous ? 1 : 0;
        };

    for (int i = 0; i < aln_len && ref_bases_aligned_count < gene_len_ref_ungapped; ++i) {
        char q_char = aln_query_raw[i];
        char t_char = aln_target_raw[i];

        if (q_char != '-') {
            if (valid_iupac_for_query.find(q_char) == std::string::npos) {
                invalid_char_count_query_out++;
            }
            else if (canonical_bases.find(q_char) == std::string::npos) {
                ambiguous_base_counts_query_out[q_char]++;
            }
        }

        if (t_char == '-' && q_char != '-' && ref_bases_aligned_count == 0) {
            std::string report_key = "Between " + previous_gene_name + " and " + current_gene_name + ": ";
            bool found_report = false;
            for (std::string& report : intergenic_insertion_details) {
                if (report.rfind(report_key, 0) == 0) {
                    report += q_char;
                    found_report = true;
                    break;
                }
            }
            if (!found_report) {
                intergenic_insertion_details.push_back(report_key + std::string(1, q_char));
            }
            intergenic_insertion_bases_count++;
            query_pos_processed_bases++;
            current_ref_codon_bases.clear(); current_query_codon_bases.clear(); nt_changes_in_current_codon.clear(); current_ref_codon_start_gene_pos_nt_1based = -1;
            waiting_for_in_frame_aa_report_after_fs_restore = false;
            continue;
        }

        if (t_char != '-') {
            if (alignment_start_ref_pos_in_gene_1based == -1) alignment_start_ref_pos_in_gene_1based = ref_pos_in_gene_1based;
            alignment_end_ref_pos_in_gene_1based = ref_pos_in_gene_1based;
        }

        if (in_insertion_event && t_char != '-') {
            if (!current_insertion_sequence.empty()) {
                distinct_insertion_events_out++;
                has_any_mutation_this_gene_out = true;

                int ins_report_genomic_pos = gene_genomic_start_pos_1based + indel_start_ref_pos_in_gene_1based - 1 - 1;
                if (indel_start_ref_pos_in_gene_1based == 1 && ref_bases_aligned_count > 0) {
                    ins_report_genomic_pos = gene_genomic_start_pos_1based - 1;
                }
                else if (ref_bases_aligned_count == 0 && indel_start_ref_pos_in_gene_1based == 1) {
                    ins_report_genomic_pos = gene_genomic_start_pos_1based - 1;
                }
                std::string mut_str = "Ins." + std::to_string(ins_report_genomic_pos) + "." + current_insertion_sequence;
                int cls = classify_mutation_chars(current_insertion_sequence);
                if (cls == 0) out_insertion_details_canonical.push_back(mut_str);
                else if (cls == 1) out_insertion_details_ambiguous.push_back(mut_str);
                else out_insertion_details_invalid.push_back(mut_str);

                // Frameshift logic only if NOT in UTR
                if (!is_utr_region) {
                    int len_ins = current_insertion_sequence.length();
                    int prev_offset = current_reading_frame_offset;
                    current_reading_frame_offset = (current_reading_frame_offset + (len_ins % 3)) % 3;
                    int fs_effective_pos = gene_genomic_start_pos_1based + indel_start_ref_pos_in_gene_1based - 1;
                    if (prev_offset == 0 && current_reading_frame_offset != 0)
                        out_frameshift_lifecycle_reports.push_back(std::string("FS_started(offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_ins(" + std::to_string(len_ins) + "bp)");
                    else if (prev_offset != 0 && current_reading_frame_offset == 0) {
                        out_frameshift_lifecycle_reports.push_back(std::string("FS_restored_from_pos_") + std::to_string(fs_effective_pos) + "_by_ins(" + std::to_string(len_ins) + "bp)");
                        if (len_ins % 3 != 0) waiting_for_in_frame_aa_report_after_fs_restore = true;
                    }
                    else if (prev_offset != 0 && current_reading_frame_offset != 0 && prev_offset != current_reading_frame_offset)
                        out_frameshift_lifecycle_reports.push_back(std::string("FS_phase_changed(new_offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_ins(" + std::to_string(len_ins) + "bp)");
                    if (len_ins % 3 != 0)
                        out_frameshift_causing_indel_reports.push_back("FS.Ins." + std::to_string(ins_report_genomic_pos) + "." + current_insertion_sequence.substr(0, std::min(10, (int)current_insertion_sequence.length())) + (current_insertion_sequence.length() > 10 ? "..." : ""));
                }
            }
            in_insertion_event = false;
            current_ref_codon_bases.clear(); current_query_codon_bases.clear(); nt_changes_in_current_codon.clear(); current_ref_codon_start_gene_pos_nt_1based = -1;
        }

        if (in_deletion_event && q_char != '-') {
            if (!current_deletion_sequence_ref_bases.empty()) {
                distinct_deletion_events_out++;
                has_any_mutation_this_gene_out = true;

                int del_start_genomic_pos = gene_genomic_start_pos_1based + indel_start_ref_pos_in_gene_1based - 1;
                int del_end_genomic_pos = gene_genomic_start_pos_1based + (ref_pos_in_gene_1based - 1) - 1;
                std::string mut_str; // Defined here to be available for FS report
                if (current_deletion_sequence_ref_bases.length() == 1) {
                    mut_str = "Del" + std::to_string(del_start_genomic_pos);
                }
                else {
                    mut_str = "Del" + std::to_string(del_start_genomic_pos) + "-" + std::to_string(del_end_genomic_pos);
                }
                out_deletion_details_canonical.push_back(mut_str);

                // Frameshift logic only if NOT in UTR
                if (!is_utr_region) {
                    int len_del = current_deletion_sequence_ref_bases.length();
                    int prev_offset = current_reading_frame_offset;
                    current_reading_frame_offset = (current_reading_frame_offset - (len_del % 3) + 3) % 3;
                    int fs_effective_pos = gene_genomic_start_pos_1based + ref_pos_in_gene_1based - 1;
                    if (prev_offset == 0 && current_reading_frame_offset != 0)
                        out_frameshift_lifecycle_reports.push_back(std::string("FS_started(offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_del(" + std::to_string(len_del) + "bp)");
                    else if (prev_offset != 0 && current_reading_frame_offset == 0) {
                        out_frameshift_lifecycle_reports.push_back(std::string("FS_restored_from_pos_") + std::to_string(fs_effective_pos) + "_by_del(" + std::to_string(len_del) + "bp)");
                        if (len_del % 3 != 0) waiting_for_in_frame_aa_report_after_fs_restore = true;
                    }
                    else if (prev_offset != 0 && current_reading_frame_offset != 0 && prev_offset != current_reading_frame_offset)
                        out_frameshift_lifecycle_reports.push_back(std::string("FS_phase_changed(new_offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_del(" + std::to_string(len_del) + "bp)");
                    if (len_del % 3 != 0)
                        out_frameshift_causing_indel_reports.push_back(mut_str.replace(0, 3, "FS.Del")); // Use mut_str defined above
                }
            }
            in_deletion_event = false;
            current_ref_codon_bases.clear(); current_query_codon_bases.clear(); nt_changes_in_current_codon.clear(); current_ref_codon_start_gene_pos_nt_1based = -1;
        }

        if (t_char == '-') {
            if (!in_insertion_event) {
                in_insertion_event = true;
                current_insertion_sequence.clear();
                indel_start_ref_pos_in_gene_1based = ref_pos_in_gene_1based;
            }
            if (q_char != '-') current_insertion_sequence += q_char;
            total_insertion_length_bp++;
            if (q_char != '-') query_pos_processed_bases++;
        }
        else if (q_char == '-') {
            if (!in_deletion_event) {
                in_deletion_event = true;
                current_deletion_sequence_ref_bases.clear();
                indel_start_ref_pos_in_gene_1based = ref_pos_in_gene_1based;
            }
            current_deletion_sequence_ref_bases += t_char;
            total_deletion_length_bp++;
            ref_pos_in_gene_1based++;
            ref_bases_aligned_count++;
        }
        else { // Match or Mismatch
            bool local_allow_aa_proc = (current_reading_frame_offset == 0);
            if (local_allow_aa_proc && waiting_for_in_frame_aa_report_after_fs_restore) {
                if (current_ref_codon_bases.empty()) {
                    if ((ref_pos_in_gene_1based - 1) % 3 == 0) {
                        waiting_for_in_frame_aa_report_after_fs_restore = false;
                    }
                    else {
                        local_allow_aa_proc = false;
                    }
                }
                else {
                    local_allow_aa_proc = false;
                }
            }

            if (local_allow_aa_proc && !is_utr_region) {
                if (current_ref_codon_bases.empty()) {
                    current_ref_codon_start_gene_pos_nt_1based = ref_pos_in_gene_1based;
                }
                current_ref_codon_bases += t_char;
                current_query_codon_bases += q_char;

                if (t_char != q_char) {
                    nt_changes_in_current_codon.push_back({ t_char, ref_pos_in_gene_1based, q_char, gene_genomic_start_pos_1based + ref_pos_in_gene_1based - 1 });
                }

                if (current_ref_codon_bases.length() == 3) {
                    bool ref_valid, query_valid;
                    char ref_aa = translate_codon(current_ref_codon_bases, genetic_code, ref_valid);
                    char query_aa = translate_codon(current_query_codon_bases, genetic_code, query_valid);
                    if (ref_valid && query_valid) {
                        int aa_pos = (current_ref_codon_start_gene_pos_nt_1based - 1) / 3 + 1;
                        bool codon_nucs_changed = !nt_changes_in_current_codon.empty();
                        std::string nt_detail_str_for_aa = "";
                        if (codon_nucs_changed) {
                            for (size_t k = 0; k < nt_changes_in_current_codon.size(); ++k) {
                                if (k > 0) nt_detail_str_for_aa += "_";
                                nt_detail_str_for_aa += std::string(1, nt_changes_in_current_codon[k].ref_base) +
                                    std::to_string(nt_changes_in_current_codon[k].ref_nt_pos_genome) +
                                    std::string(1, nt_changes_in_current_codon[k].query_base);
                            }
                        }
                        std::string aa_report_str = (nt_detail_str_for_aa.empty() ? "" : nt_detail_str_for_aa + " ") +
                            "(" + std::string(1, ref_aa) + std::to_string(aa_pos) + std::string(1, query_aa) + ")";
                        if (ref_aa != query_aa) {
                            has_any_mutation_this_gene_out = true;
                            if (query_aa == '*') out_nonsense_aa_muts.push_back(aa_report_str);
                            else out_missense_aa_muts.push_back(aa_report_str);
                        }
                        else if (codon_nucs_changed) {
                            has_any_mutation_this_gene_out = true;
                            out_synonymous_aa_muts.push_back(aa_report_str);
                        }
                    }
                    current_ref_codon_bases.clear(); current_query_codon_bases.clear(); nt_changes_in_current_codon.clear(); current_ref_codon_start_gene_pos_nt_1based = -1;
                }
            }
            else {
                if (!current_ref_codon_bases.empty()) {
                    current_ref_codon_bases.clear(); current_query_codon_bases.clear(); nt_changes_in_current_codon.clear(); current_ref_codon_start_gene_pos_nt_1based = -1;
                }
            }

            if (t_char != q_char) {
                distinct_mismatch_events_out++;
                has_any_mutation_this_gene_out = true;
                if (is_utr_region) {
                    utr_mutation_count_out++;
                }
                int genomic_pos_mismatch = gene_genomic_start_pos_1based + ref_pos_in_gene_1based - 1;
                std::string mut_str = std::string(1, t_char) + "." + std::to_string(genomic_pos_mismatch) + "." + std::string(1, q_char);
                int cls = classify_mutation_chars({ q_char });
                if (t_char == 'N') cls = std::max(cls, classify_mutation_chars({ t_char }));
                if (cls == 0) out_mismatch_details_canonical.push_back(mut_str);
                else if (cls == 1) out_mismatch_details_ambiguous.push_back(mut_str);
                else out_mismatch_details_invalid.push_back(mut_str);
            }
            else {
                matches_count++;
            }
            query_pos_processed_bases++;
            ref_pos_in_gene_1based++;
            ref_bases_aligned_count++;
        }
    }

    if (in_insertion_event) {
        if (!current_insertion_sequence.empty()) {
            distinct_insertion_events_out++;
            has_any_mutation_this_gene_out = true;

            int ins_report_genomic_pos = gene_genomic_start_pos_1based + indel_start_ref_pos_in_gene_1based - 1 - 1;
            if (indel_start_ref_pos_in_gene_1based == 1 && ref_bases_aligned_count > 0) {
                ins_report_genomic_pos = gene_genomic_start_pos_1based - 1;
            }
            else if (ref_bases_aligned_count == 0 && indel_start_ref_pos_in_gene_1based == 1) {
                ins_report_genomic_pos = gene_genomic_start_pos_1based - 1;
            }
            std::string mut_str_report = "Ins." + std::to_string(ins_report_genomic_pos) + "." + current_insertion_sequence;
            int cls = classify_mutation_chars(current_insertion_sequence);
            if (cls == 0) out_insertion_details_canonical.push_back(mut_str_report);
            else if (cls == 1) out_insertion_details_ambiguous.push_back(mut_str_report);
            else out_insertion_details_invalid.push_back(mut_str_report);

            // Frameshift logic only if NOT in UTR
            if (!is_utr_region) {
                int len_ins = current_insertion_sequence.length();
                int prev_offset_val = current_reading_frame_offset;
                current_reading_frame_offset = (current_reading_frame_offset + (len_ins % 3)) % 3;

                int fs_effective_pos = gene_genomic_start_pos_1based + indel_start_ref_pos_in_gene_1based - 1;
                if (prev_offset_val == 0 && current_reading_frame_offset != 0)
                    out_frameshift_lifecycle_reports.push_back(std::string("FS_started(offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_ins(" + std::to_string(len_ins) + "bp)_at_aln_end");
                else if (prev_offset_val != 0 && current_reading_frame_offset == 0) {
                    out_frameshift_lifecycle_reports.push_back(std::string("FS_restored_from_pos_") + std::to_string(fs_effective_pos) + "_by_ins(" + std::to_string(len_ins) + "bp)_at_aln_end");
                    // waiting_for_in_frame_aa_report_after_fs_restore is not typically set at alignment end as processing stops.
                }
                else if (prev_offset_val != 0 && current_reading_frame_offset != 0 && prev_offset_val != current_reading_frame_offset)
                    out_frameshift_lifecycle_reports.push_back(std::string("FS_phase_changed(new_offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_ins(" + std::to_string(len_ins) + "bp)_at_aln_end");

                if (len_ins % 3 != 0)
                    out_frameshift_causing_indel_reports.push_back("FS.Ins." + std::to_string(ins_report_genomic_pos) + "." + current_insertion_sequence.substr(0, std::min(10, (int)current_insertion_sequence.length())) + (current_insertion_sequence.length() > 10 ? "..." : ""));
            }
        }
    }
    if (in_deletion_event) {
        if (!current_deletion_sequence_ref_bases.empty()) {
            distinct_deletion_events_out++;
            has_any_mutation_this_gene_out = true;

            int del_start_genomic_pos = gene_genomic_start_pos_1based + indel_start_ref_pos_in_gene_1based - 1;
            int del_end_genomic_pos = gene_genomic_start_pos_1based + (ref_pos_in_gene_1based - 1) - 1;
            std::string mut_str_report; 
            if (current_deletion_sequence_ref_bases.length() == 1) {
                mut_str_report = "Del" + std::to_string(del_start_genomic_pos);
            }
            else {
                mut_str_report = "Del" + std::to_string(del_start_genomic_pos) + "-" + std::to_string(del_end_genomic_pos);
            }
            out_deletion_details_canonical.push_back(mut_str_report);

            // Frameshift logic only if NOT in UTR
            if (!is_utr_region) {
                int len_del = current_deletion_sequence_ref_bases.length();
                int prev_offset_val = current_reading_frame_offset;
                current_reading_frame_offset = (current_reading_frame_offset - (len_del % 3) + 3) % 3;

                int fs_effective_pos = gene_genomic_start_pos_1based + (ref_pos_in_gene_1based - 1);
                if (prev_offset_val == 0 && current_reading_frame_offset != 0)
                    out_frameshift_lifecycle_reports.push_back(std::string("FS_started(offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_del(" + std::to_string(len_del) + "bp)_at_aln_end");
                else if (prev_offset_val != 0 && current_reading_frame_offset == 0) {
                    out_frameshift_lifecycle_reports.push_back(std::string("FS_restored_from_pos_") + std::to_string(fs_effective_pos) + "_by_del(" + std::to_string(len_del) + "bp)_at_aln_end");
                }
                else if (prev_offset_val != 0 && current_reading_frame_offset != 0 && prev_offset_val != current_reading_frame_offset)
                    out_frameshift_lifecycle_reports.push_back(std::string("FS_phase_changed(new_offset ") + (current_reading_frame_offset == 1 ? "+1" : "+2") + ")_from_pos_" + std::to_string(fs_effective_pos) + "_by_del(" + std::to_string(len_del) + "bp)_at_aln_end");

                if (len_del % 3 != 0)
                    out_frameshift_causing_indel_reports.push_back(mut_str_report.replace(0, 3, "FS.Del"));
            }
        }
    }

    if (alignment_start_ref_pos_in_gene_1based != -1 && gene_len_ref_ungapped > 0) {
        int covered_len = alignment_end_ref_pos_in_gene_1based - alignment_start_ref_pos_in_gene_1based + 1;
        coverage_out = 100.0f * static_cast<float>(covered_len) / gene_len_ref_ungapped;
        coverage_out = std::min(coverage_out, 100.0f);
    }
    else {
        coverage_out = 0.0f;
    }

    int total_aligned_positions_for_identity_calc = matches_count + distinct_mismatch_events_out + total_insertion_length_bp + total_deletion_length_bp;
    if (total_aligned_positions_for_identity_calc > 0) {
        identity_out = 100.0f * static_cast<float>(matches_count) / total_aligned_positions_for_identity_calc;
    }
    else if (coverage_out > 0 && ref_bases_aligned_count > 0 && matches_count == ref_bases_aligned_count && distinct_mismatch_events_out == 0 && total_insertion_length_bp == 0 && total_deletion_length_bp == 0) {
        identity_out = 100.0f;
    }
    else if (ref_bases_aligned_count == 0 && coverage_out == 0.0f) {
        identity_out = 0.0f;
    }
    else {
        identity_out = 0.0f;
    }
    return total_insertion_length_bp - total_deletion_length_bp;
}



int main(int argc, char* argv[]) {
    omp_set_dynamic(0);
    omp_set_num_threads(omp_get_num_procs());
    std::cout << "Using " << omp_get_max_threads() << " threads." << std::endl;
    namespace fs = std::filesystem;
    double total_estimated_sequential_time_ns = 0.0;

    // Configuring file paths
    fs::path base_path = fs::current_path();

    fs::path genes_layout_file;
    std::string ref_path_str;
    std::string query_path_str;
    std::string tsv_output_path_str;

    if (argc >= 4) {
        // Mode: Setup via arguments
        genes_layout_file = argv[1];
        ref_path_str = argv[2];
        query_path_str = argv[3];
        
        if (argc >= 5) {
             tsv_output_path_str = argv[4];
        } else {
             // Automatic output name based on query file
             fs::path qp(query_path_str);
             fs::path results_folder = base_path / "results";
             if (!fs::exists(results_folder)) try { fs::create_directory("results"); } catch(...) {}
             tsv_output_path_str = (results_folder / ("analysis_" + qp.stem().string() + ".tsv")).string();
        }
        
        std::cout << "--- Configuration from Arguments ---" << std::endl;
        std::cout << "Genes Layout: " << genes_layout_file << std::endl;
        std::cout << "Reference:    " << ref_path_str << std::endl;
        std::cout << "Query:        " << query_path_str << std::endl;
        std::cout << "Output:       " << tsv_output_path_str << std::endl;
        std::cout << "------------------------------------" << std::endl;

    } else {
        // Mode: Default (Fallback to SARS-CoV-2 Test)
        std::cout << "No arguments provided. Using default configuration (SARS-CoV-2 Test Example)." << std::endl;
        std::cout << "Usage: " << argv[0] << " <genes_layout_file> <reference_file> <query_file> [output_file]" << std::endl;

        fs::path data_folder = base_path / "data";
        if (!fs::exists(data_folder)) {
            data_folder = base_path / "../data";
        }

        genes_layout_file = data_folder / "examples/tests/SARS_CoV2_test/pilot_genes.tsv"; 
        ref_path_str = (data_folder / "examples/tests/SARS_CoV2_test/SARS_CoV2_reference_genome.fasta").string();
        query_path_str = (data_folder / "examples/tests/SARS_CoV2_test/gisaid_hcov-19_2025_05_26_16.fasta").string();
        
        fs::path results_folder = base_path / "results";
        if (!fs::exists(results_folder)) {
            results_folder = base_path / "../results";
            // Attempt to create results directory if it doesn't exist
            if (!fs::exists(results_folder)) {
                 try { fs::create_directory("results"); results_folder = base_path / "results"; } catch(...) {}
            }
        }
        tsv_output_path_str = (results_folder / "analysis_results_SARS_CoV2_test.tsv").string();
    }
    
    fs::path gene_file_path = genes_layout_file; // Alias for compatibility with rest of main


    std::vector<int> gene_starts_1based;
    std::vector<int> gene_ends_1based;
    std::vector<std::string> gene_names_list;

    std::ifstream gene_info_file(gene_file_path);
    if (!gene_info_file.is_open()) {
        std::cerr << "Fatal error: Could not open gene information file: " << gene_file_path << std::endl;
        return 1;
    }

    std::string file_line;

    if (std::getline(gene_info_file, file_line)) {
        // Header line read and discarded
    }
    else {
        std::cerr << "Error: Could not read header line or file is empty: " << gene_file_path << std::endl;
        gene_info_file.close();
        return 1;
    }

    while (std::getline(gene_info_file, file_line)) {
        if (file_line.empty() || file_line[0] == '#') {
            continue;
        }

        std::istringstream iss(file_line);
        std::string name_field, start_field, end_field, description_field;

        if (!std::getline(iss, name_field, '\t') ||
            !std::getline(iss, start_field, '\t') ||
            !std::getline(iss, end_field, '\t')) {
            std::cerr << "Error parsing first 3 fields in gene file line: " << file_line << std::endl;
            std::cerr << "Expected format: Gene_Name<TAB>Start_Pos<TAB>End_Pos<TAB>Description" << std::endl;
            continue;
        }

        if (!std::getline(iss, description_field)) {
            description_field = "";
        }

        try {
            gene_names_list.push_back(name_field);
            gene_starts_1based.push_back(std::stoi(start_field));
            gene_ends_1based.push_back(std::stoi(end_field));
        }
        catch (const std::exception& e) {
            std::cerr << "Error converting positions or processing data for gene file line: " << file_line
                << "\n  -> " << e.what() << std::endl;
            continue;
        }
    }
    gene_info_file.close();
    // Reading FASTA sequences
    std::string ref_sequence = read_fasta_sequence_only(ref_path_str);
    auto query_sequences_list = read_multi_fasta(query_path_str);

    if (ref_sequence.empty()) {
        std::cerr << "Fatal error: The reference sequence is empty." << std::endl;
        return 1;
    }
    if (query_sequences_list.empty()) {
        std::cerr << "Warning: No query sequences loaded." << std::endl;
        exit(1);
    }

    std::atomic<int> processed_queries_count(0);
    const int total_queries_to_process = static_cast<int>(query_sequences_list.size());
    auto overall_start_time = std::chrono::high_resolution_clock::now();

    // Setting up the Parasail scoring matrix
    const char* alphabet = "ACGTURYKMSWBDHVN";
    std::map<char, std::string> ambiguity_map = {
    {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"}, {'U', "T"},
    {'R', "AG"}, {'Y', "CT"}, {'K', "GT"}, {'M', "AC"},
    {'S', "GC"}, {'W', "AT"}, {'B', "GTC"}, {'D', "GAT"},
    {'H', "ACT"}, {'V', "GCA"}, {'N', "ACGT"}
    };

    parasail_matrix_t* matrix = parasail_matrix_create(alphabet, 3, -3);
    int size = matrix->size;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            char a = alphabet[i];
            char b = alphabet[j];
            int score;

            if (a == b) {
                score = 3;
            }
            else {
                const std::string& group_a = ambiguity_map[a];
                const std::string& group_b = ambiguity_map[b];

                bool compatible = false;
                for (char x : group_a) {
                    if (group_b.find(x) != std::string::npos) {
                        compatible = true;
                        break;
                    }
                }

                if (a == 'N' && b == 'N') score = 0;
                else score = compatible ? 1 : -3;
            }

            matrix->user_matrix[i * size + j] = score;
        }
    }


    // Preparing the TSV output file
    std::ofstream tsv_output_file(tsv_output_path_str);
    if (!tsv_output_file.is_open()) {
        std::cerr << "Fatal error: Could not open TSV output file: " << tsv_output_path_str << std::endl;
        parasail_matrix_free(matrix);
        return 1;
    }
    // Write TSV header
    tsv_output_file << "Query_ID\tGene_Name\tGene_Ref_Start\tGene_Ref_End\tCoverage_Percent\tIdentity_Percent\tAlignment_Score\t"
        << "Query_Alignment_Starts_At\tQuery_Consumed_Up_To\t"
        << "Distinct_Mismatch_Events\tCanonical_Mismatch_Details\tAmbiguous_Mismatch_Details\tInvalid_Mismatch_Details\t"
        << "Distinct_Insertion_Events\tCanonical_Insertion_Details\tAmbiguous_Insertion_Details\tInvalid_Insertion_Details\t"
        << "Distinct_Deletion_Events\tCanonical_Deletion_Details\tAmbiguous_Deletion_Details\tInvalid_Deletion_Details\t"
        << "Missense_Mutation_Count\tMissense_Mutation_Details\t"
        << "Nonsense_Mutation_Count\tNonsense_Mutation_Details\t"
        << "Synonymous_Mutation_Count\tSynonymous_Mutation_Details\t"
        << "Frameshift_Causing_Indel_Count\tFrameshift_Causing_Indel_Details\t"
        << "Frameshift_Lifecycle_Summary\t"
        << "UTR_Mutation_Count\t"
        << "Intergenic_Insertions_Before_Gene\tNon_ACGT_Query_Base_Counts\tInvalid_Query_Base_Count\n";

    // Pre-calculation of the initial offset for each query using SSW
    int ref_seq_len = ref_sequence.length();
    int ssw_ref_window_len = std::max(100, ref_seq_len / 60);
    int ssw_query_window_len = std::max(50, ref_seq_len / 300);
    std::vector<int> query_initial_offsets(total_queries_to_process);

#pragma omp parallel for reduction(+:total_estimated_sequential_time_ns)
    for (int i = 0; i < total_queries_to_process; ++i) {
        const std::string& current_query_s = query_sequences_list[i].second;
        int actual_ssw_ref_win = std::min(ssw_ref_window_len, ref_seq_len);
        int actual_ssw_query_win = std::min(ssw_query_window_len, (int)current_query_s.length());

        std::string query_sub_ssw = actual_ssw_query_win > 0 ? current_query_s.substr(0, actual_ssw_query_win) : "";
        std::string ref_sub_ssw = actual_ssw_ref_win > 0 ? ref_sequence.substr(0, actual_ssw_ref_win) : "";

        int ssw_offset = 0;
        if (!query_sub_ssw.empty() && !ref_sub_ssw.empty()) {
            auto ssw_task_start_time = std::chrono::high_resolution_clock::now();
            parasail_result_ssw_t* ssw_res = parasail_ssw(
                query_sub_ssw.c_str(), query_sub_ssw.length(),
                ref_sub_ssw.c_str(), ref_sub_ssw.length(),
                7, 2, matrix 
            );
            ssw_offset = ssw_res->read_begin1 - ssw_res->ref_begin1;
            parasail_result_ssw_free(ssw_res);
            auto ssw_task_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> ssw_task_duration_ns =
                ssw_task_end_time - ssw_task_start_time;
            total_estimated_sequential_time_ns += ssw_task_duration_ns.count();

        }
        query_initial_offsets[i] = ssw_offset;
    }

    // Global counters for final summary
    const size_t num_total_genes = gene_starts_1based.size();
    std::atomic<long long> total_genes_analyzed_count(0);
    std::map<char, std::atomic<long long>> overall_ambiguous_base_counts;
    for (char c : std::string("NURYKMSWBDHV")) overall_ambiguous_base_counts[c] = 0;
    std::atomic<long long> overall_invalid_chars_count(0);
    std::atomic<long long> overall_intergenic_bases_count(0);
    std::atomic<long long> total_genes_with_any_mutation(0);


    // Main processing loop (parallelized by query)
#pragma omp parallel for schedule(dynamic) reduction(+:total_estimated_sequential_time_ns)
    for (int q_idx = 0; q_idx < total_queries_to_process; ++q_idx) {
        std::ostringstream tsv_lines_for_this_query_buffer;
        int query_current_consumed_end_coord = query_initial_offsets[q_idx];

        const std::string& current_query_id = query_sequences_list[q_idx].first;
        const std::string& current_query_seq_data = query_sequences_list[q_idx].second;

        for (size_t g_idx = 0; g_idx < num_total_genes; ++g_idx) {
            auto gene_processing_start_time = std::chrono::high_resolution_clock::now();
            bool can_process_this_gene_iteration = true;

            int gene_start_1based = gene_starts_1based[g_idx];
            int gene_len = gene_ends_1based[g_idx] - gene_start_1based + 1;

            if (gene_start_1based - 1 < 0 || (gene_start_1based - 1 + gene_len) > static_cast<int>(ref_sequence.size())) {
                can_process_this_gene_iteration = false;
            }

            std::string gene_ref_segment;
            if (can_process_this_gene_iteration) {
                gene_ref_segment = ref_sequence.substr(gene_start_1based - 1, gene_len);
            }

            int align_query_from_pos = query_current_consumed_end_coord;
            if (can_process_this_gene_iteration) {
                if (align_query_from_pos < 0) align_query_from_pos = 0;
                if (align_query_from_pos >= static_cast<int>(current_query_seq_data.size())) {
                    can_process_this_gene_iteration = false;
                }
            }

            std::string query_segment_for_alignment;
            if (can_process_this_gene_iteration) {
                int padding_min = 50;
                int padding_dynamic = gene_len / 2;
                int padding_extra_right = std::min(200, std::max(padding_min, padding_dynamic));
                if (gene_len < 50) padding_extra_right = std::min(200, std::max(padding_extra_right, gene_len));

                int query_segment_len = gene_len + padding_extra_right;
                if (align_query_from_pos + query_segment_len > static_cast<int>(current_query_seq_data.size())) {
                    query_segment_len = current_query_seq_data.size() - align_query_from_pos;
                }
                if (query_segment_len <= 0) {
                    can_process_this_gene_iteration = false;
                }
                else {
                    query_segment_for_alignment = current_query_seq_data.substr(align_query_from_pos, query_segment_len);
                }
            }

            if (can_process_this_gene_iteration) {
                parasail_result_t* alignment_result = parasail_sg_trace_scan_16(
                    gene_ref_segment.c_str(), gene_ref_segment.length(),
                    query_segment_for_alignment.c_str(), query_segment_for_alignment.length(),
                    7, 3, matrix);

                if (alignment_result && alignment_result->score > -20000) {
                    parasail_traceback_t* traceback_data = parasail_result_get_traceback(
                        alignment_result,
                        gene_ref_segment.c_str(), gene_ref_segment.length(),
                        query_segment_for_alignment.c_str(), query_segment_for_alignment.length(),
                        matrix, '|', '.', ' ');

                    if (traceback_data) {
                        float coverage_val = 0.0f, identity_val = 0.0f;
                        std::vector<std::string> intergenic_details_list;
                        int intergenic_bases = 0;
                        std::map<char, int> query_amb_bases_map;
                        int query_invalid_chars = 0;
                        bool gene_had_mutation = false;
                        int d_mismatches = 0, d_insertions = 0, d_deletions = 0;

                        std::vector<std::string> mis_canon, mis_amb, mis_inv;
                        std::vector<std::string> ins_canon, ins_amb, ins_inv;
                        std::vector<std::string> del_canon, del_amb, del_inv;
                        std::vector<std::string> aa_miss, aa_syn, aa_non;
                        std::vector<std::string> fs_indel_reps, fs_life_reps;

                        // Determine if current gene is UTR
                        std::string current_gene_name_str = gene_names_list[g_idx];
                        bool is_current_gene_utr = (current_gene_name_str == "5_UTR" || current_gene_name_str == "3_UTR");
                        int utr_mut_count_for_gene = 0;


                        int gene_net_offset = analyze_alignment(
                            traceback_data->ref, traceback_data->query, // query is aln_query, ref is aln_target
                            gene_start_1based, gene_len,
                            is_current_gene_utr, 
                            coverage_val, identity_val,
                            intergenic_bases, intergenic_details_list,
                            current_gene_name_str, (g_idx > 0 ? gene_names_list[g_idx - 1] : "START_OF_GENOME"),
                            query_amb_bases_map, query_invalid_chars, gene_had_mutation,
                            d_insertions, d_deletions, d_mismatches,
                            mis_canon, mis_amb, mis_inv,
                            ins_canon, ins_amb, ins_inv,
                            del_canon, del_amb, del_inv,
                            aa_miss, aa_syn, aa_non,
                            fs_indel_reps, fs_life_reps,
                            utr_mut_count_for_gene 
                        );

                        total_genes_analyzed_count++;
                        if (gene_had_mutation) total_genes_with_any_mutation++;
                        for (const auto& pair : query_amb_bases_map) overall_ambiguous_base_counts[pair.first] += pair.second;
                        overall_invalid_chars_count += query_invalid_chars;
                        overall_intergenic_bases_count += intergenic_bases;

                        int next_query_start_coordinate = align_query_from_pos + intergenic_bases + gene_len + gene_net_offset;

                        tsv_lines_for_this_query_buffer
                            << current_query_id << "\t" << current_gene_name_str << "\t"
                            << gene_start_1based << "\t" << gene_ends_1based[g_idx] << "\t"
                            << std::fixed << std::setprecision(2) << coverage_val << "\t"
                            << std::fixed << std::setprecision(2) << identity_val << "\t"
                            << alignment_result->score << "\t"
                            << align_query_from_pos + 1 << "\t"      // Query_Alignment_Starts_At
                            << next_query_start_coordinate << "\t"   // Query_Consumed_Up_To
                            << d_mismatches << "\t" << join_strings(mis_canon, ";") << "\t" << join_strings(mis_amb, ";") << "\t" << join_strings(mis_inv, ";") << "\t"
                            << d_insertions << "\t" << join_strings(ins_canon, ";") << "\t" << join_strings(ins_amb, ";") << "\t" << join_strings(ins_inv, ";") << "\t"
                            << d_deletions << "\t" << join_strings(del_canon, ";") << "\t" << join_strings(del_amb, ";") << "\t" << join_strings(del_inv, ";") << "\t"
                            << aa_miss.size() << "\t" << join_strings(aa_miss, ";") << "\t"
                            << aa_non.size() << "\t" << join_strings(aa_non, ";") << "\t"
                            << aa_syn.size() << "\t" << join_strings(aa_syn, ";") << "\t"
                            << fs_indel_reps.size() << "\t" << join_strings(fs_indel_reps, ";") << "\t"
                            << join_strings(fs_life_reps, "; ") << "\t"
                            << utr_mut_count_for_gene << "\t" 
                            << join_strings(intergenic_details_list, " ") << "\t"
                            << format_base_counts_map(query_amb_bases_map) << "\t"
                            << query_invalid_chars << "\n";

                        query_current_consumed_end_coord = next_query_start_coordinate;
                        parasail_traceback_free(traceback_data);
                    }
                    else { // Traceback Failure
                        query_current_consumed_end_coord = align_query_from_pos + gene_len;
                    }
                    parasail_result_free(alignment_result);
                }
                else { // Alignment failure or very low score
                    if (alignment_result) parasail_result_free(alignment_result); // Release if created
                    query_current_consumed_end_coord = align_query_from_pos + gene_len;
                }
            }
            auto gene_processing_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::nano> gene_processing_duration_ns =
                gene_processing_end_time - gene_processing_start_time;
            total_estimated_sequential_time_ns += gene_processing_duration_ns.count();
        }
#pragma omp critical
        {
            tsv_output_file << tsv_lines_for_this_query_buffer.str();
            processed_queries_count++;
            std::cout << "\rQueries Processed: " << processed_queries_count << "/" << total_queries_to_process
                << " (" << std::fixed << std::setprecision(1)
                << (100.0 * processed_queries_count / total_queries_to_process) << "%)" << std::flush;
        }
    }

    tsv_output_file.close();
    parasail_matrix_free(matrix);

    auto overall_end_time = std::chrono::high_resolution_clock::now();

    auto total_execution_time = std::chrono::duration<double>(
        overall_end_time - overall_start_time
    );

    std::cout << "\n\nProcess completed." << std::endl;
    std::cout << "Results saved in: " << tsv_output_path_str << std::endl;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Total execution time (parallel): "
        << total_execution_time.count() << " seconds" << std::endl;
    double final_sequential_time_s = total_estimated_sequential_time_ns / 1e9; 
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Estimated total sequential execution time (SSW + Gene Alignments): "
        << final_sequential_time_s << " seconds" << std::endl;
    // Optional global summary in console
    std::cout << "\n--- Global Summary (Console) ---" << std::endl;
    std::cout << "Total queries processed: " << processed_queries_count.load() << std::endl;
    std::cout << "Total gene analysis attempted: " << total_genes_analyzed_count.load() << std::endl;
    std::cout << "Total genes with any mutation detected: " << total_genes_with_any_mutation.load() << std::endl;
    std::cout << "Global counts of ambiguous bases in queries:" << std::endl;
    for (char c_rep : std::string("NURYKMSWBDHV")) {
        if (overall_ambiguous_base_counts[c_rep].load() > 0) {
            std::cout << "  " << c_rep << ": " << overall_ambiguous_base_counts[c_rep].load() << std::endl;
        }
    }
    std::cout << "Total bases in intergenic insertions (global): " << overall_intergenic_bases_count.load() << std::endl;
    std::cout << "Total invalid characters in queries (global): " << overall_invalid_chars_count.load() << std::endl;

    return 0;
}
