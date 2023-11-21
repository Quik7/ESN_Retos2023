require './gene.rb' # Import of Gene object
require './protein.rb' # Import of Protein object

class Network
  attr_accessor :network_id, :num_nodes, :members

  # Class variable to store all InteractionNetwork instances
  @@total_network_objects = {}
  @@number_of_networks = 0

  # Initialize a new InteractionNetwork object with given parameters
  def initialize(params = {})
    @network_id = params.fetch(:network_id, @@number_of_networks + 1)
    @num_nodes = params.fetch(:num_nodes, 0)
    @members = params.fetch(:members, Hash.new)

    # Store the new object and increment network count
    @@total_network_objects[@network_id] = self
    @@number_of_networks += 1
  end

  # Returns all InteractionNetwork objects
  def self.all_networks
    @@total_network_objects
  end

  # Creates a new InteractionNetwork
  def self.create_network
    network_id = @@number_of_networks + 1
    new(network_id: network_id, num_nodes: 2, members: Hash.new)
    network_id
  end

  # Adds nodes to a specified network
  def self.add_nodes2network(network_id)
    network = @@total_network_objects[network_id]
    network.num_nodes += 1 if network
  end

  # Adds a gene to a specified network
  def self.add_gene2network(network_id, gene_object)
    #puts "Adding Gene to Network: Gene ID - #{gene_object.gene_id}, Network ID - #{network_id}"  # Debugging
    network = @@total_network_objects[network_id]
    if network
      network.members[gene_object.gene_id] = gene_object
      gene_object.annotate
    end
  end

  # Assigns a protein to a network and updates the network recursively
  def self.assign(protein_object, network_id)
    protein_object.network = network_id
    update_network_for_ppis(protein_object, network_id)

    # If the protein has an associated gene, add it to the network
    gene = Gene.all_genes[protein_object.prot_id]
    add_gene2network(network_id, gene) if gene
  end

  

  # Helper method to update network based on protein-protein interactions (PPIs)
  def self.update_network_for_ppis(protein_object, network_id)
    # Return immediately if there are no protein-protein interactions stored in $PPIS.
    return unless $PPIS && $PPIS.any?
    $PPIS.each do |ppi|
      # Determine the interactant ID by checking if either of the PPI pair matches the current protein's IntAct ID.
      # If neither matches, interactant_id is set to nil.
      interactant_id = ppi[0] == protein_object.intact_id ? ppi[1] : ppi[1] == protein_object.intact_id ? ppi[0] : nil
      next unless interactant_id

      interactant = Protein.all_prots_withintact[interactant_id]
      next unless interactant && interactant.network.nil?

      add_nodes2network(network_id)
      # Recursively assign the interactant protein to the same network.
      # This will also trigger updates for any proteins interacting with this interactant.
      assign(interactant, network_id)
    end
  end
end

  
  