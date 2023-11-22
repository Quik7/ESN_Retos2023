#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: Protein Object
# author: Enrique Solera Navarro
#-----------------------------------------------------

require 'net/http'
require 'json'

# The Protein class represents a protein with its unique identifiers and network associations
class Protein
  attr_accessor :prot_id, :intact_id, :network

  # Class variables to store all protein instances
  @@total_protein_objects = {}
  @@total_protwithintact_objects = {}

  # Initialize a new Protein object with given parameters
  # @param params [Hash] optional parameters for initializing a Protein object.
  def initialize(params = {})
    @prot_id = params.fetch(:prot_id, "XXXXXX") # Default UniProt ID
    @intact_id = params.fetch(:intact_id, nil) # IntAct ID
    @network = params.fetch(:network, nil) # Associated network ID

    # Store the new object in class variables
    @@total_protein_objects[@prot_id] = self
    @@total_protwithintact_objects[@intact_id] = self if @intact_id
  end

  # Returns all protein objects
  # @return [Hash] all stored Protein objects.
  def self.all_prots
    @@total_protein_objects
  end

  # Returns all proteins with an IntAct ID
  # @return [Hash] all Protein objects with an IntAct ID.
  def self.all_prots_withintact
    @@total_protwithintact_objects
  end

  # Creates a Protein object with provided details and checks for interactions
  # @param prot_id [String] the UniProt ID of the protein.
  # @param level [Integer] the interaction depth level.
  # @param gene_id [String] optional gene ID.
  # @param intact [String] optional IntAct ID.
  def self.create_prot(prot_id, level, gene_id = nil, intact = nil)
    #puts "Creating Protein object with UniProt ID: #{prot_id}, IntAct ID: #{intact}, at level: #{level}"
    intact ||= get_prot_intactcode(gene_id) if level.zero?
    if intact && (level < $MAX_LEVEL)
      PPI.create_ppis(intact, level)
    end
    new(prot_id: prot_id, intact_id: intact)
    level += 1
  end

  # Checks if a protein with the given IntAct ID exists
  # @param intact_id [String] the IntAct ID to check.
  # @return [Boolean] whether the Protein exists.
  def self.exists(intact_id)
    @@total_protwithintact_objects.has_key?(intact_id)
  end

  # Retrieves the IntAct code for a given gene ID
  # @param gene_id [String] the gene ID.
  # @return [String, nil] the IntAct code, or nil if not found.
  def self.get_prot_intactcode(gene_id)
    uri = URI("http://togows.org/entry/ebi-uniprot/#{gene_id}/dr.json")
    response = Net::HTTP.get_response(uri)
    data = JSON.parse(response.body)
    data[0]['IntAct']&.first&.first
  end

  # Converts an IntAct accession code to a UniProt ID
  # @param intact_accession_code [String] the IntAct accession code.
  # @return [String] the corresponding UniProt ID.
  def self.intactcode2protid(intact_accession_code)
    intact_accession_code.chomp!
    if valid_uniprot_accession?(intact_accession_code)
      address = URI("http://togows.org/entry/ebi-uniprot/#{intact_accession_code}/entry_id.json")
      response = Net::HTTP.get_response(address)
      data = JSON.parse(response.body)
      data[0] || "UniProt ID not found"
    else
      "UniProt ID not found"
    end
  end

  private
  # https://www.rubyguides.com/2018/10/method-visibility/
  # private keyword is used to denote that methods like validate_gene_id, get_response, annotate_kegg, and annotate_go are meant 
  # for internal use within the class. They are part of the class's internal logic and are not intended to be accessed directly 
  # from outside the class.




  # Validates a UniProt accession code
  # @param code [String] the UniProt accession code to validate.
  # @return [Boolean] whether the code is valid.
  def self.valid_uniprot_accession?(code)
    code.match?(/[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/)
  end
end

# The PPI class deals with the protein-protein interactions
class PPI
  
  @@ppis = []

  # Returns all protein-protein interactions
  # @return [Array] all stored protein-protein interactions.
  def self.all_ppis
    @@ppis
  end

  # Fetches and filters protein-protein interactions for a given protein IntAct code
  # @param intact_code [String] the IntAct code of the protein.
  # @return [Array, nil] a list of protein-protein interactions, or nil if none found
  def self.get_ppis(intact_code)
    uri = URI("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{intact_code}?format=tab25")
    res = Net::HTTP.get_response(uri)
    return nil unless res

    parse_ppi_response(res.body)
  end

  # Creates protein objects for each interacting protein pair
  # @param intact_id [String] the IntAct ID of the protein.
  # @param level [Integer] the interaction depth level.
  def self.create_ppis(intact_id, level)
    level += 1
    return if level >= $MAX_LEVEL

    ppis = get_ppis(intact_id)
    ppis&.each do |id1, id2|
      #Preventing Infinite Loops in Interactions
      next if Protein.exists(id1) || Protein.exists(id2)
      Protein.create_prot(Protein.intactcode2protid(id1), level, nil, id1) if id2 == intact_id
      Protein.create_prot(Protein.intactcode2protid(id2), level, nil, id2) if id1 == intact_id
    end
  end

  private

  # Parses the response from PPI API and filters interactions based on criteria
  # @param body [String] the response body from PPI API.
  # @return [Array] parsed and filtered list of interactions.
  def self.parse_ppi_response(body)
    body.split("\n").each_with_object([]) do |line, ppis|
      fields = line.split("\t")
      p1, p2, intact_miscore = fields[0].sub(/uniprotkb:/, ""), fields[1].sub(/uniprotkb:/, ""), fields[14].sub(/intact-miscore:/, "").to_f

      next unless valid_interaction?(p1, p2, intact_miscore, fields)

      interaction = [p1, p2].sort
      ppis << interaction unless @@ppis.include?(interaction)
      @@ppis << interaction
    end
  end

  # Parses a line from the PPI response into protein IDs and interaction score
  # @param line [String] a line of text from the PPI response.
  # @return [Array] parsed protein IDs and interaction score.
  def self.parse_line(line)
    fields = line.split("\t")
    [fields[0].sub(/uniprotkb:/, ""), fields[1].sub(/uniprotkb:/, ""), fields[14].sub(/intact-miscore:/, "").to_f]
  end

  # Checks if the interaction is valid based on specific criteria
  # @param p1 [String] the first protein ID.
  # @param p2 [String] the second protein ID.
  # @param intact_miscore [Float] the interaction score.
  # @param fields [Array] additional fields from the PPI response.
  # @return [Boolean] whether the interaction is valid.
  def self.valid_interaction?(p1, p2, intact_miscore, fields)
    # Filter 1: we remove self interactions
    # Filter 2: we remove iteractions that have a score below the value we want
    # Filter 3: we remove interactions detected by "two hybrib" and "two hybrid array" (too many false positives)
    # Filter 4: we remove those interactions from proteins which are not from Arabidopsis Thaliana
    return false if p1 == p2 || intact_miscore < 0.4
    fields[6] !~ /psi-mi:"MI:0018"/ && fields[6] !~ /psi-mi:"MI:0397"/ && fields[9] =~ /taxid:3702/ && fields[10] =~ /taxid:3702/
  end
end

