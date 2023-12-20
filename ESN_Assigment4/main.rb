#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment4: MAIN SCRIPT
# author: Enrique Solera Navarro
#-----------------------------------------------------------------
# Searching for Orthologues
#-----------------------------------------------------------------
# Since the A. thaliana database consists of coding sequences (CDS), it will be converted into a protein database. 
# This conversion facilitates the use of standard blastp queries rather than relying on computationally intensive 
# translated tblastn or blastx methods. It's important to note that this approach is only viable when the 
# nucleotide database is specifically limited to CDS. Nevertheless, it is conserved in the program the option to perform
# a blastx / tblastn if the user want to do it.
#-----------------------------------------------------------------
# PARAMETERS CHOSEN
#-----------------------------------------------------------------
# The E-value, or Expectation value, is a key metric in BLAST searches. It estimates the number of times an alignment 
# with a given score or better could occur by chance in a database search. A lower E-value indicates a more significant match. 
# An E-value of 10e-6 suggests that the likelihood of the alignment occurring by random chance is very low, thus implying 
# a more significant and biologically relevant match
# The bit score (S’) is another important measure in BLAST results. It represents the quality of the alignment between the 
# query sequences and the target sequence. A higher bit score indicates a higher similarity and, consequently, a more reliable alignment. 
# The bit score is advantageous because it is normalized (adjusted for scoring parameters such as the substitution matrix and gap penalty), 
# making it useful for comparing alignments obtained using different scoring parameters. A minimum bit score of 50 is a reasonable threshold 
# to ensure that the alignments are of a decent quality,
# Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, 
# Bioinformatics, Volume 24, Issue 3, February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585
#------------------------------------------------------------------
# CONTINUED ANALYSIS
#------------------------------------------------------------------
# 1. Phylogenetic Analysis: Construct a phylogenetic tree using the sequences of the identified orthologues along with other related sequences. 
# This helps in understanding the evolutionary relationships and confirming orthology. The code would involve aligning the sequences, selecting 
# an appropriate model of evolution, and then building and visualizing the phylogenetic tree.
# 2. Conservation of Protein Domains and Motifs: Examine the protein domains and motifs in the orthologue candidates to check for conservation. 
# This step often involves using tools like Pfam or InterProScan to identify and compare functional domains. The analysis would focus on whether 
# the critical domains and motifs are conserved across the species, which supports the hypothesis of orthology.
# 3. Functional Studies: Investigate the function of the genes in their respective organisms. This might involve looking at gene expression data, 
# knockout studies, or literature research to see if the genes have similar roles in their respective organisms. The code for such an analysis might 
# involve querying gene expression databases, performing statistical analysis on expression data, or text mining scientific literature for functional information.
# 4. Synteny Analysis: Examine the genomic context of these genes in different species. Orthologous genes often exhibit conserved synteny, meaning 
# they are located in similar genomic contexts in different species. This analysis would involve comparing the genomic regions surrounding the orthologue 
# candidates across different species.
# 5. Evolutionary Rate Analysis: Compare the rate of evolution of these genes across species. Orthologues are expected to have similar rates of evolution. 
# This would involve calculating and comparing the non-synonymous (dN) and synonymous (dS) substitution rates for these genes.
#-----------------------------------------------------------------

require 'bio'

# Constants and Default Parameters
E_VAL = 10e-6  # E-value threshold for BLAST
MIN_BIT_SCORE = 50  # Minimum bit score to consider a hit significant
DB_DIR = "Databases"  # Directory to store BLAST databases

# Function to ask user for input
def ask_input(prompt)
  puts prompt
  STDOUT.flush
  gets.chomp
end

# Function to check if file exists
def check_file(filename)
  # Aborts the script if the file does not exist
  abort "Error: File #{filename} does not exist" unless File.exist?(filename)
  filename
end

# Function to guess the sequence type (protein or nucleotide)
def guess_sequence_type(file)
  # Reads the first sequence from the file and estimates if it is nucleotide or protein based on ATCG content
  first_seq = Bio::FastaFormat.open(file).first.seq
  if first_seq.count('ATCGatcg') > 0.85 * first_seq.size  # More than 85% ATCG content
    'nucl'  # It's a nucleotide sequence
  else
    'prot'  # It's a protein sequence
  end
end

# Function to create BLAST database
def create_blast_db(filename, dbtype, dbname)
  # Creates a BLAST database from the given file, silencing the output
  puts "Creating BLAST database for #{filename}..."
  system("makeblastdb -in '#{filename}' -dbtype #{dbtype} -out #{DB_DIR}/#{dbname} > /dev/null 2>&1")
  puts "Database created successfully."
end

# Function to translate nucleotide sequence to protein
def translate_to_protein(nucleotide_file, output_file)
  # Translates a nucleotide sequence file to protein sequence file
  puts "Translating nucleotide sequence from #{nucleotide_file}..."
  nucleotide_db = Bio::FastaFormat.open(nucleotide_file)
  File.open(output_file, "w") do |file|
    nucleotide_db.each do |entry|
      file.puts(">#{entry.entry_id}")
      file.puts(Bio::Sequence.auto(entry.seq).translate)
    end
  end
  puts "Translation completed. Protein sequences saved to #{output_file}."
end

# Main Script
puts "Reciprocal Best BLAST Orthologue Search using BioRuby"

# Get input file names from the user
search_file = check_file(ask_input("Enter the name of the search file (either protein or nucleotide sequence):"))
target_file = check_file(ask_input("Enter the path of the target file (either protein or nucleotide sequence):"))

# Determine the sequence types of the input files
search_file_type = guess_sequence_type(search_file)
target_file_type = guess_sequence_type(target_file)

# Check if files are compatible (one protein, one nucleotide)
if [search_file_type, target_file_type].sort == ['nucl', 'prot']
  puts "You have a protein and nucleotide sequence. How do you want to perform the analysis?"
  method_choice = ask_input("Choose BLAST method: 1 for blastx/tblastn, 2 for translating nucleotide to protein and blastp:")

  # Create BLAST database directory
  Dir.mkdir(DB_DIR) unless Dir.exist?(DB_DIR)

  # Choose the BLAST method based on user input
  if method_choice == "1"
    # User has chosen the blastx/tblastn method
    puts "You have chosen the blastx/tblastn method"

    # Create unique database names based on the file names
    db_search = "#{File.basename(search_file)}_db"
    db_target = "#{File.basename(target_file)}_db"

    # Create BLAST databases for the input files
    # For the search file (protein sequence), create a protein database
    # For the target file (nucleotide sequence), create a nucleotide database
    if search_file_type == 'prot' and target_file_type == 'nucl'
      create_blast_db(search_file, 'prot', db_search)
      create_blast_db(target_file, 'nucl', db_target)

      # Initialize BLAST factories for searching
      # blastx: Compare protein query against nucleotide database
      # tblastn: Compare nucleotide database against protein query
      factory_search = Bio::Blast.local('blastx', "#{DB_DIR}/#{db_search}", "-e 10e-6 -sorthits 1")
      factory_target = Bio::Blast.local('tblastn', "#{DB_DIR}/#{db_target}", "-e 10e-6 -sorthits 1")
    elsif search_file_type == 'nucl' and target_file_type == 'prot'
      create_blast_db(search_file, 'nucl', db_search)
      create_blast_db(target_file, 'prot', db_target)
      factory_search = Bio::Blast.local('tblastn', "#{DB_DIR}/#{db_search}", "-e 10e-6 -sorthits 1")
      factory_target = Bio::Blast.local('blastx', "#{DB_DIR}/#{db_target}", "-e 10e-6 -sorthits 1")
    end
  elsif method_choice == "2"
    # User has chosen the translate and blastp method
    puts "You have chosen the translate and blastp method"

    # Translate the nucleotide sequence to protein and save to a new file
    translated_target_file = "#{target_file}_translated.fa"
    translate_to_protein(target_file, translated_target_file)

    # Create unique database names based on the file names
    db_search = "#{File.basename(search_file)}_db"
    db_target = "#{File.basename(translated_target_file)}_db"

    # Create BLAST databases for the input files (both protein sequences)
    create_blast_db(search_file, 'prot', db_search)
    create_blast_db(translated_target_file, 'prot', db_target)

    # Initialize BLAST factories for searching (blastp: protein-protein comparison)
    factory_search = Bio::Blast.local('blastp', "#{DB_DIR}/#{db_search}", "-e 10e-6 -sorthits 1")
    factory_target = Bio::Blast.local('blastp', "#{DB_DIR}/#{db_target}", "-e 10e-6 -sorthits 1")
  else
    # Exit the script if an invalid method choice is made
    abort "Invalid method choice. Please choose 1 or 2."
  end

  # Exit the script if both files are not of expected types (one protein and one nucleotide)
  else
    abort "Both files need to be either a protein or a nucleotide sequence."
  end

  # Open an output file to write the results of the BLAST search
  File.open('output_orthologues.txt', 'w') do |out_file|
    # Write a header line indicating the files being compared
    out_file.puts "Orthologues between files: #{search_file} and #{target_file}"

    # Open the search file and (depending on the method) the target file or its translated version
    ff_search = Bio::FastaFormat.open(search_file)
    ff_target = Bio::FastaFormat.open(method_choice == "1" ? target_file : translated_target_file)

    # Create a hash to store the target sequences, indexed by their IDs
    target_hash = ff_target.each_with_object({}) { |seq, hash| hash[seq.entry_id] = seq.seq }

    # Iterate through each sequence in the search file to find orthologues
    ff_search.each do |seq_search|
      # Perform BLAST search and proceed only if hits are found that meet the criteria (e-value and bit score)
      report_target = factory_target.query(seq_search)
      next if report_target.hits.empty? || report_target.hits.first.evalue > E_VAL || report_target.hits.first.bit_score < MIN_BIT_SCORE

      # Extract the ID of the best hit from the target
      target_id = report_target.hits.first.definition.split(/\s|\|/).first

      # Perform a reciprocal BLAST search using the best hit's sequence
      report_search = factory_search.query(">#{target_id}\n#{target_hash[target_id]}")
      next if report_search.hits.empty? || report_search.hits.first.evalue > E_VAL || report_search.hits.first.bit_score < MIN_BIT_SCORE

      # Check if the best hit from the reciprocal search matches the original search sequence
      match_id = report_search.hits.first.definition.split(/\s|\|/).first
      next unless seq_search.entry_id == match_id

      # Write the pair of matching IDs (orthologues) to the output file
      out_file.puts "#{seq_search.entry_id} <-> #{target_id}"

      # Print a message indicating a reciprocal hit was found
      puts "Reciprocal hit found between #{seq_search.entry_id} and #{target_id}."
      puts "..."
    end
end

# Indicate completion of the BLAST search and the location of the results file
puts "Reciprocal BLAST search complete. Results in 'output_orthologues.txt'"




