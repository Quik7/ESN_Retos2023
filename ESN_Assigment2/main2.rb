puts ""

# Import the different modules needed
require './gene.rb'
require './protein.rb'
require './network.rb'
require 'net/http'

# Check the input file
unless ARGV[0] && ARGV[1]
    abort "USAGE: main.rb geneList.txt output.txt"
end

unless File.exists?(ARGV[0])
    abort "Error: File #{ARGV[0]} does not exist"
end

# Function for data retrieval from URI
def uri_fetch(uri_str)
    address = URI(uri_str)  
    response = Net::HTTP.get_response(address)
    
    case response
      when Net::HTTPSuccess then
        return response
      else
        abort "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
    end 
end

# Function to write results to an output file
def write_record(networks, filename)
    if File.exists?(filename) 
        File.delete(filename) 
    end
  
    File.open(filename, "a+") do |fnr|
      fnr.puts "--------------------------------"
      fnr.puts "Assigment2 RECORD"
      fnr.puts "author: Enrique Solera Navarro"
      fnr.puts "--------------------------------\n\n"
      fnr.puts "For every network created we display:"
      fnr.puts "  1) Network ID"
      fnr.puts "  2) Number of nodes in the network"
      fnr.puts "  3) Genes form the list that are included in the network"
      fnr.puts "      (Depth of interactions selected: #{$MAX_LEVEL})"
      fnr.puts "\n\n\n"
  
      networks.each do |id_net, network_obj|  
          fnr.puts "NETWORK NUMBER: #{id_net}"
          fnr.puts "\tNodes_number: #{network_obj.num_nodes}"
          fnr.puts "\tGenes associated and annotations:"
          network_obj.members.each do |id, gene|
              fnr.puts "\t\t-- #{id}"
              gene.kegg.each do |kegg_id, kegg_name|
                  fnr.puts "\t\t\t\t\t KEGG Pathways ID: #{kegg_id};\tKEGG Pathways Name: #{kegg_name}"
              end
              gene.go.each do |go_id, go_term_name|
                  fnr.puts "\t\t\t\t\t GO ID: #{go_id};\tGO Term Name: #{go_term_name}"
              end
          end
          fnr.puts "\n-----------------------------------------------\n"
      end
  end
end

$MAX_LEVEL = 2
puts "Obtaining all interacting proteins deeping #{$MAX_LEVEL.to_s} levels ..."
Gene.load_from_file(ARGV[0])

$PPIS = PPI.all_ppis
puts "DONE!\n\n"

puts "Obtaining networks ..."
Protein.all_prots_withintact.each do |id, prot_object|
    if not prot_object.network
        new_network = Network.create_network
        Network.assign(prot_object, new_network)
    end
end
puts "DONE!\n\n"

puts "Recording report..."
write_record(Network.all_networks, ARGV[1])
puts "TASKS COMPLETED, output #{ARGV[1]} recorded!"


if ARGV[0] == "-help" || ARGV[0] == "-h"
  ayuda
  exit
end