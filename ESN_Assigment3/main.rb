#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment3: MAIN SCRIPT
# author: Enrique Solera Navarro
#-----------------------------------------------------------------
# The use of BioRuby is focused on specific tasks related to gene sequence analysis, 
# without necessitating a complex object-oriented architecture that would require multiple classes
# BioRuby itself provides a set of classes and modules for biological computation. In the script, 
# we are leveraging these existing classes (like Bio::Sequence) without needing to encapsulate the functionality further. 
# This means we are already using object-oriented principles where they are most beneficial, provided by the BioRuby library

#Import the different modules needed
require 'rest-client'
require 'bio'
require 'uri'

# Reads gene identifiers from a file, removing duplicates.
# @param filename [String] The file containing gene identifiers.
# @return [Array<String>] An array of unique gene identifiers.
def load_from_file(filename)
  genes = []
  File.open(filename, "r") do |file|
    file.each_line { |line| genes << line.chomp }
  end
  genes.uniq
end

# Fetches data from a specified URI.
# @param uri_str [String] The URI to fetch data from.
# @return [String] The response body if successful.
# @raise [RuntimeError] Raises an error if the fetch fails.
def uri_fetch(uri_str)
  begin
    response = RestClient.get(uri_str)
    response.body
  rescue RestClient::ExceptionWithResponse => e
    raise "Failed to fetch URI #{uri_str}: #{e.response}"
  end
end

# Fetches and converts gene information to a Bio::Sequence object.
# @param gene_id [String] The gene's identifier.
# @return [Bio::Sequence, nil] The gene information, or nil if not found.
def obtain_gene_info(gene_id)
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id}"
  begin
    response_body = uri_fetch(address)
    return nil if response_body.nil? || response_body.empty?

    entry = Bio::EMBL.new(response_body)
    entry.to_biosequence
  rescue => e
    puts "Error obtaining gene info for #{gene_id}: #{e.message}"
    nil
  end
end


# Determines if any of the target's matches are included in a given exon
# @param [String] exon_id Identifier of the exon
# @param [Array<Integer>] target_sequence_matches Array of match start positions
# @param [Integer] len_seq Length of the sequence
# @param [Array<Integer>] exon_position Start and end positions of the exon
# @param [String] strand Strand of the exon ('+' or '-')
# @return [Hash] Match positions and corresponding exon details
def find_target_in_exon(exon_id, target_sequence_matches, len_seq, exon_position, strand)
  target_in_exon = {}

  target_sequence_matches.each do |match_init|
    match_end = match_init + $len_target - 1

    if match_within_exon?(match_init, match_end, exon_position)
      match_init, match_end = convert_positions(match_init, match_end, len_seq) if strand == '-'
      target_in_exon[[match_init, match_end]] = [exon_id, strand]
    end
  end

  target_in_exon unless target_in_exon.empty?
end

private

# Checks if the match is within the exon boundaries
def match_within_exon?(match_init, match_end, exon_position)
  (match_init >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
end

# Converts positions for the reverse strand
def convert_positions(match_init, match_end, len_seq)
  [len_seq - match_end, len_seq - match_init]
end

# Finds target sequences within the exons of a given Bio::EMBL object
# @param [Bio::EMBL] bio_seq_object The Bio::EMBL object to process
# @return [Hash] Positions of targets within exons and their strand information
def get_exons_targets(bio_seq_object)
  len_bio_seq = bio_seq_object.length
  target_positions_in_exon = {}

  target_matches_in_seq_forward = find_target_matches(bio_seq_object, $target)
  target_matches_in_seq_reverse = find_target_matches(bio_seq_object.reverse_complement, $target)

  bio_seq_object.features.each do |feature|
    next unless feature.feature == 'exon' && feature.position !~ /[A-Z]/

    exon_id, position, strand = parse_exon_feature(feature, len_bio_seq)
    target_pos_in_exon = find_target_in_exon(exon_id, strand == '+' ? target_matches_in_seq_forward : target_matches_in_seq_reverse, len_bio_seq, position, strand)

    target_positions_in_exon.merge!(target_pos_in_exon) if target_pos_in_exon
  end

  target_positions_in_exon
end

private

# Finds target sequence matches in a given sequence
# Finds target sequence matches in a given nucleotide sequence
# @param [Bio::Sequence::NA] sequence Nucleotide sequence to search
# @param [String] target Target sequence to find
# @return [Array<Integer>] Array of match start positions
def find_target_matches(sequence, target)
  match_positions = []
  regex = Regexp.new(target)

  sequence.to_s.scan(regex) do |match|
    match_positions << Regexp.last_match.begin(0)
  end

  match_positions
end


# Parses an exon feature to extract relevant information
def parse_exon_feature(feature, len_bio_seq)
  exon_id = feature.qualifiers[0].value.gsub('exon_id=', '')
  position = feature.position
  strand = position.include?('complement') ? '-' : '+'

  position = parse_position(position, len_bio_seq, strand)
  [exon_id, position, strand]
end

# Parses the position string to get start and end positions
def parse_position(position, len_bio_seq, strand)
  position = position.gsub('complement', '').delete('()').split('..').map(&:to_i)
  strand == '-' ? position.map { |pos| len_bio_seq - pos }.reverse : position
end

# Adds target match features to the Bio::EMBL object and writes to a GFF3 file
# @param [String] gene_id The gene identifier
# @param [Hash] targets The hash containing target positions and exon information
# @param [Bio::EMBL] bioseq The Bio::EMBL object to which features are added
def add_features(gene_id, targets, bioseq)
  targets.each do |target, exonid_strand|
    add_feature_to_bioseq(target, exonid_strand, bioseq)
    write_feature_to_gff(gene_id, target, exonid_strand)
  end
end

private

# Adds a new feature to the Bio::EMBL object
def add_feature_to_bioseq(target, exonid_strand, bioseq)
  feat = Bio::Feature.new("#{feature_name(target)}_in_exon", "#{target[0]}..#{target[1]}")
  feat.append(Bio::Feature::Qualifier.new('nucleotide_motif', "#{feature_name(target)}_in_#{exonid_strand[0]}"))
  feat.append(Bio::Feature::Qualifier.new('strand', exonid_strand[1]))
  bioseq.features << feat
end

# Writes the feature to the GFF3 file
def write_feature_to_gff(gene_id, target, exonid_strand)
  $gff_genes.puts format_gff3_entry(gene_id, target, exonid_strand)
end

# Formats the GFF3 entry
def format_gff3_entry(gene_id, target, exonid_strand)
  "#{gene_id}\t.\t#{feature_name(target)}\t#{target[0]}\t#{target[1]}\t.\t#{exonid_strand[1]}\t.\tID=#{exonid_strand[0]}"
end

# Returns the feature name
def feature_name(target)
  "#{$target.upcase}"
end

# Extracts chromosome information from a Bio::Sequence object and writes to a GFF3 file
# @param [String] gene_id The gene identifier
# @param [Bio::Sequence] bio_seq_object The Bio::Sequence object to process
# @return [Array<String, String, String>, false] Chromosome number, start and end positions or false if not found
def get_chromosome(gene_id, bio_seq_object)
  primary_accession = bio_seq_object.primary_accession
  return false unless primary_accession

  chrom_array = primary_accession.split(":")
  write_chromosome_info_to_gff(gene_id, chrom_array)

  [chrom_array[2], chrom_array[3], chrom_array[4]]
end

private

# Writes chromosome information to the GFF3 file
def write_chromosome_info_to_gff(gene_id, chrom_array)
  $gff_chr.puts format_gff3_chromosome_entry(gene_id, chrom_array)
end

# Formats the GFF3 chromosome entry
def format_gff3_chromosome_entry(gene_id, chrom_array)
  "#{chrom_array[2]}\t.\tgene\t#{chrom_array[3]}\t#{chrom_array[4]}\t.\t+\t.\tID=#{gene_id}"
end

# Opens a file for writing, creating a new file or overwriting if it already exists
# @param [String] filename The name of the file to be opened
# @return [File] The File object
def create_open_file(filename)
  File.delete(filename) if File.exists?(filename)
  File.open(filename, "a+")
end

# Converts target positions to chromosome coordinates and writes to a GFF3 file
# @param [String] gene The gene identifier
# @param [Hash] targets Hash containing target positions and exon information
# @param [Array] chr Array with chromosome information
def convert_to_chr(gene, targets, chr)
  targets.each do |positions, exon_strand|
    pos_ini_chr = chr[1].to_i + positions[0].to_i
    pos_end_chr = chr[1].to_i + positions[1].to_i
    
    write_target_to_gff(chr[0], pos_ini_chr, pos_end_chr, exon_strand, gene)
  end
end

private

# Writes target information to the GFF3 file
def write_target_to_gff(chromosome, start_pos, end_pos, exon_strand, gene)
  $gff_chr.puts "#{chromosome}\t.\tnucleotide_motif\t#{start_pos}\t#{end_pos}\t.\t#{exon_strand[1]}\t.\tID=#{exon_strand[0]};parent=#{gene}"
end

# Ensures a gene list file is provided as an argument
abort "USAGE: main.rb geneList.txt" unless ARGV[0]

# Checks if the provided file exists
abort "Error: File #{ARGV[0]} does not exist" unless File.exist?(ARGV[0])

# Sets up global variables
$target = "cttctt"
$len_target = $target.length

# Define the process_gene method
def process_gene(gene)
seq_obj = obtain_gene_info(gene)
return if seq_obj.nil?

target_hash = get_exons_targets(seq_obj)
if target_hash.empty?
  $no_targets.puts gene
else
  add_features(gene, target_hash, seq_obj)
  chr = get_chromosome(gene, seq_obj)
  convert_to_chr(gene, target_hash, chr)
end
end

# Main execution with updated progress control
puts "Working on the tasks..."

$gff_genes = create_open_file("genes.gff3: Enrique Solera Navarro")
$gff_chr = create_open_file("chromosomes.gff3: Enrique Solera Navarro")
$no_targets = create_open_file("genes_without_target.txt: Enrique Solera Navarro")

# Printing headers
$gff_genes.puts "##gff-version 3"
$gff_chr.puts "##gff_chr-version 3"
$no_targets.puts "Genes without #{$target.upcase} in exons\n\n"

genes_ids = load_from_file(ARGV[0])
total_genes = genes_ids.length

genes_ids.each_with_index do |gene, index|
  process_gene(gene)

  # Update progress every 20 genes
  if (index + 1) % 20 == 0 || index + 1 == total_genes
    percentage_complete = ((index + 1) / total_genes.to_f * 100).round(2)
    puts "Processed #{index + 1}/#{total_genes} genes (#{percentage_complete}%)"
  end
end

puts "DONE!"
puts "Check the output in genes.gff3, chromosomes.gff3, and genes_without_target.txt"