require 'bio'

# Constants and Default Parameters
E_VAL = 10e-6
MIN_BIT_SCORE = 50
DB_DIR = "Databases"

# Function to ask user for input
def ask_input(prompt)
  puts prompt
  STDOUT.flush
  gets.chomp
end

# Function to check if file exists
def check_file(filename)
  abort "Error: File #{filename} does not exist" unless File.exist?(filename)
  filename
end

# Function to guess the sequence type (protein or nucleotide)
def guess_sequence_type(file)
    first_seq = Bio::FastaFormat.open(file).first.seq
    if first_seq.count('ATCGatcg') > 0.85 * first_seq.size # More than 85% ATCG content
      'nucl'
    else
      'prot'
    end
end

# Function to create BLAST database
def create_blast_db(filename, dbtype, dbname)
  system("makeblastdb -in '#{filename}' -dbtype #{dbtype} -out #{DB_DIR}/#{dbname} > /dev/null 2>&1")
end

# Function to translate nucleotide sequence to protein
def translate_to_protein(nucleotide_file, output_file)
  nucleotide_db = Bio::FastaFormat.open(nucleotide_file)
  File.open(output_file, "w") do |file|
    nucleotide_db.each do |entry|
      file.puts(">#{entry.entry_id}")
      file.puts(Bio::Sequence.auto(entry.seq).translate)
    end
  end
end

# Main Script
puts "Reciprocal Best BLAST Orthologue Search using BioRuby"

# Get input file names
search_file = check_file(ask_input("Enter the path to the first file (either protein or nucleotide sequence):"))
target_file = check_file(ask_input("Enter the path to the second file (either protein or nucleotide sequence):"))

# Guess the sequence type
search_file_type = guess_sequence_type(search_file)
target_file_type = guess_sequence_type(target_file)

# Check if one is protein and the other nucleotide
if [search_file_type, target_file_type].sort == ['nucl', 'prot']
  puts "You have a protein and nucleotide sequence. How do you want to perform the analysis?"
  method_choice = ask_input("Choose BLAST method: 1 for blastx/tblastn, 2 for translating nucleotide to protein and blastp:")

  # Create BLAST database directory
  Dir.mkdir(DB_DIR) unless Dir.exist?(DB_DIR)

  if method_choice == "1"
    puts " You have chosen the blastx/tblastn method"
    db_search = "#{File.basename(search_file)}_db"
    db_target = "#{File.basename(target_file)}_db"
    puts " Generating the databases for the analysis"
    create_blast_db(search_file, 'prot', db_search)
    create_blast_db(target_file, 'nucl', db_target)
    puts "DONE ! Now performing the analysis"
    factory_search = Bio::Blast.local('tblastn', "#{DB_DIR}/#{db_search}", "-e 10e-6 -sorthits 1")
    factory_target = Bio::Blast.local('blastx', "#{DB_DIR}/#{db_target}"), "-e 10e-6 -sorthits 1"
  elsif method_choice == "2"
    puts "You have chosen the translate and blastp method"
    translated_target_file = "#{target_file}_translated.fa"
    puts "Translating nucleotide sequence to protein sequence"
    translate_to_protein(target_file, translated_target_file)
    puts "DONE! Now generating the databases to perform the analysis"
    db_search = "#{File.basename(search_file)}_db"
    db_target = "#{File.basename(translated_target_file)}_db"
    create_blast_db(search_file, 'prot', db_search)
    create_blast_db(translated_target_file, 'prot', db_target)
    puts "DONE! Now performing the analysis"
    factory_search = Bio::Blast.local('blastp', "#{DB_DIR}/#{db_search}", "-e 10e-6 -sorthits 1")
    factory_target = Bio::Blast.local('blastp', "#{DB_DIR}/#{db_target}", "-e 10e-6 -sorthits 1")
  else
    abort "Invalid method choice. Please choose 1 or 2."
  end
else
  abort "Both files need to be either a protein or a nucleotide sequence."
end

# Open output file and conduct BLAST search
File.open('output_orthologues.txt', 'w') do |out_file|
  out_file.puts "Orthologues between files: #{search_file} and #{target_file}"

  ff_search = Bio::FastaFormat.open(search_file)
  ff_target = Bio::FastaFormat.open(method_choice == "1" ? target_file : translated_target_file)
  target_hash = ff_target.each_with_object({}) { |seq, hash| hash[seq.entry_id] = seq.seq }

  ff_search.each do |seq_search|
    report_target = factory_target.query(seq_search)
    next if report_target.hits.empty? || report_target.hits.first.evalue > E_VAL || report_target.hits.first.bit_score < MIN_BIT_SCORE

    target_id = report_target.hits.first.definition.split(/\s|\|/).first
    report_search = factory_search.query(">#{target_id}\n#{target_hash[target_id]}")
    next if report_search.hits.empty? || report_search.hits.first.evalue > E_VAL || report_search.hits.first.bit_score < MIN_BIT_SCORE

    match_id = report_search.hits.first.definition.split(/\s|\|/).first
    next unless seq_search.entry_id == match_id

    out_file.puts "#{seq_search.entry_id} <-> #{target_id}"
  end
end

puts "Reciprocal BLAST search complete. Results in 'output_orthologues.txt'"



