#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment3: MAIN SCRIPT
# author: Enrique Solera Navarro
#-----------------------------------------------------------------

#Import the different modules needed
require 'net/http'
require 'bio'

def load_from_file(filename)
    # Method to read the genes from the file
    File.readlines(filename).map(&:chomp).uniq
end

def fetch_gene_sequences(gene_ids)
    Bio::NCBI.default_email = "enrique.solera.navarro@alumnos.upm.es"  # Replace with your email
    ncbi = Bio::NCBI::REST.new
  
    gene_ids.map do |id|
      ncbi.efetch(id, {"db" => "nucleotide", "rettype" => "fasta", "retmode" => "text"})
    end
end

def scan_exons_for_target_sequence(gene_sequence, target_sequence = "CTTCTT")
    Bio::FlatFile.new(Bio::GenBank, StringIO.new(gene_sequence)).each_entry do |entry|
      entry.features.each do |feature|
        if feature.feature == 'exon'
          exon_seq = entry.naseq.splicing(feature.position)
          if exon_seq.include?(target_sequence)
            puts "Target sequence found in exon of gene #{entry.accession}"
          end
        end
      end
    end
end

gene_ids = load_from_file(ARGV[0])
annotated_gene_sequences = fetch_gene_sequences(gene_ids)

annotated_gene_sequences.each do |gene_seq|
  scan_exons_for_target_sequence(gene_seq)
end