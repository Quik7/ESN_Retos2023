#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: GENE Object
# author: Enrique Solera Navarro
#-----------------------------------------------------
require 'net/http'
require 'json'
require './protein.rb'
require './annotation.rb'

class Gene
  attr_accessor :gene_id, :prot_id, :kegg, :go, :uso_general

  # Class variable to store all gene instances
  @@total_gene_objects = {}

  # Class variable to cache network responses
  @@cached_responses = {}

  # Initialize a new Gene object with given parameters
  # Initialize a new Gene object with given parameters.
  # @param params [Hash] Parameters for initializing a Gene object with keys :gene_id, :prot_id, :kegg, :go, :uso_general
  def initialize(params = {})
    @gene_id = validate_gene_id(params.fetch(:gene_id, "AT0G00000"))
    @prot_id = params.fetch(:prot_id, "XXXXXX") # Default UniProt ID
    @kegg = params.fetch(:kegg, {}) # KEGG annotations, if any
    @go = params.fetch(:go, {}) # GO annotations, if any
    # Store the new object in the class variable
    @@total_gene_objects[@prot_id] = self
    @uso_general = UsoGeneral.new
    
  end

  # Class method to return all gene objects.
  # @return [Hash] All stored Gene objects.
  def self.all_genes
    @@total_gene_objects
  end

  # Class method to fetch UniProt ID for a given gene ID.
  # @param gene_id [String] The gene ID to fetch the UniProt ID for.
  # @return [String] The fetched UniProt ID.
  def self.get_prot_id(gene_id)
    uri = URI("http://togows.org/entry/ebi-uniprot/#{gene_id}/entry_id.json")
    response = get_response(uri)
    #puts "Response from UniProt ID fetch: #{response}"  # Debugging
    response[0]
  end

  # Class method to create gene objects from a file.
  # @param filename [String] Name of the file containing gene IDs.
  def self.load_from_file(filename)
    File.open(filename, 'r') do |file|
      file.each_line do |line|
        gene_id = line.chomp
        #puts "Read gene ID from file: #{gene_id}"  # Debugging
        prot_id = get_prot_id(gene_id)
        #puts "Fetched UniProt ID for gene: #{prot_id}"  # Debugging
        Gene.new(gene_id: gene_id, prot_id: prot_id)
        #puts "Created Gene object with ID: #{@gene_id}, associated with Protein ID: #{@prot_id}"
        Protein.create_prot(prot_id, 0, gene_id)
      end
    end
  end

  # Instance method to annotate gene with KEGG and GO data
  def annotate
    annotate_kegg
    annotate_go
    @kegg.each { |id, name| @uso_general.add_annotation('KEGG', id, name) }
    @go.each { |id, term| @uso_general.add_annotation('GO', id, term) }
  end


  
  # Validates the gene ID format.
  # @param gene_id [String] The gene ID to validate.
  # @return [String] The validated gene ID.
  # @raise [RuntimeError] If the gene ID is not in the correct format.
  def validate_gene_id(gene_id)
    if gene_id.match?(/A[Tt]\d[Gg]\d\d\d\d\d/)
      gene_id
    else
      puts "Error: the Gene ID #{gene_id} is not in the correct format. It should be ATXGXXXXX"
      exit
    end
  end

  # # Fetches response from a URI and caches it.
  # @param uri [URI] The URI to fetch the response from.
  # @return [Hash, Array] Parsed JSON response.
  # @raise [JSON::ParserError] If there's an error parsing the response.
  def self.get_response(uri)
    @@cached_responses[uri] ||= Net::HTTP.get_response(uri)
    JSON.parse(@@cached_responses[uri].body)
  rescue JSON::ParserError
    raise "Error parsing response from #{uri}"
  end

  # Annotates the gene with KEGG pathways
  def annotate_kegg
    uri = URI("http://togows.org/entry/kegg-genes/ath:#{@gene_id}/pathways.json")
    response = self.class.get_response(uri)
    response[0]&.each do |path_id, path_name|
      # We insert a new key (KEGG Pathway ID) with its corresponding value (pathway name)
      @kegg[path_id] = path_name
    end
  end

  # Annotates the gene with GO terms
  def annotate_go
    uri = URI("http://togows.org/entry/ebi-uniprot/#{@gene_id}/dr.json")
    response = self.class.get_response(uri)
  
    # Check if the response has the expected data structure
    if response.is_a?(Array) && response[0].is_a?(Hash) && response[0]['GO']
      response[0]['GO'].each do |num|
        # We must check the GO refers to a biological process (it will start with a 'P')
        # We insert a new key (GO ID) with its corresponding value (GO term name)
        @go[num[0]] = num[1].sub(/^P:/, "") if num[1].start_with?("P:")
      end
    end
  end
end
