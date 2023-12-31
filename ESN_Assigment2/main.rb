#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: Main script
# author: Enrique Solera Navarro
#-----------------------------------------------------
# Import the required modules
require './gene.rb'
require './protein.rb'
require './network.rb'
require 'net/http'
require './annotation.rb'

# Define a constant for maximum interaction depth
$MAX_LEVEL = 2

# Function to fetch data from a URI with basic error handling.
# @param uri_str [String] the URI from which data is to be fetched.
# @return [Net::HTTPResponse] HTTP response object if successful.
# @raise [RuntimeError] if unable to fetch data from the URI.
def uri_fetch(uri_str)
  address = URI(uri_str)  
  response = Net::HTTP.get_response(address)

  case response
  when Net::HTTPSuccess then
    return response
  else
    abort "Error: Unable to fetch data from #{uri_str}. Error type: #{response.class}"
  end 
end

# Function to write network information to a file.
# @param networks [Hash] hash containing network ID and network objects.
# @param filename [String] the name of the file to write to.
def write_record(networks, filename)
  File.open(filename, "w") do |file|
    file.puts "Network Analysis Report"
    file.puts "Depth of Interaction Analysis: #{$MAX_LEVEL}\n\n"
    file.puts "--------------------------------------------\n\n"

    networks.each do |id_net, network_obj|
      file.puts "Network ID: #{id_net}, Number of Nodes: #{network_obj.num_nodes}"
      file.puts "Genes in Network:"
      network_obj.members.each do |id, gene|
        file.puts "\tGene ID: #{id}"
        # Access annotations from uso_general
        gene.uso_general.get_annotations('KEGG').each do |kegg_id, name|
          file.puts "\t\tKEGG Pathway: #{kegg_id} - #{name}"
        end

        gene.uso_general.get_annotations('GO').each do |go_id, term|
          file.puts "\t\tGO Term: #{go_id} - #{term}"
        end
        #gene.kegg.each { |kegg_id, name| file.puts "\t\tKEGG Pathway: #{kegg_id} - #{name}" }
        #gene.go.each { |go_id, term| file.puts "\t\tGO Term: #{go_id} - #{term}" }
      end
      file.puts "\n--------------------------------------------\n"
    end
  end
end


# Main execution block
# Prints usage information to the console.
def ayuda
  puts 'Assigment #2 - Interaction Network'
  puts 'Enrique Solera Navarro 2023'
  puts "\nUsage: ruby main.rb geneList.txt output.txt\n"
  puts 'geneList.txt : co-expressed gene list'
  puts 'output.txt : generated file with the report of the information linking these predicted sub-sets into known regulatory networks'
end

# Main function to process input and output files for network analysis.
# @param input_file [String] the name of the input file.
# @param output_file [String] the name of the output file.
def main(input_file, output_file)
  if ARGV[0] == "-help" || ARGV[0] == "-h"
    ayuda
    exit
  end
  
  abort "USAGE: ruby main.rb geneList.txt output.txt" unless input_file && output_file
  abort "Error: File #{input_file} does not exist" unless File.exist?(input_file)
  if File.exist?(output_file)
    print "The file #{output_file} already exist. Do you want to overwrite? (y/n): "
    response = STDIN.gets.chomp.downcase
    unless response == 'y'
      puts "Change the name of the output file."
      abort
      return
    end
    puts "The file #{output_file} has been updated."
  end
  puts "Processing gene list and building interaction networks..."
  Gene.load_from_file(input_file)

  $PPIS = PPI.all_ppis
  # puts "Total PPIs loaded: #{$PPIS.size}" DEBUGGING
  Protein.all_prots_withintact.each do |id, prot_object|
    unless prot_object.network
      new_network = Network.create_network
      Network.assign(prot_object, new_network)
    end
  end
  puts "Finished"
  puts "Writing network analysis to #{output_file}..."
  write_record(Network.all_networks, output_file)
  puts "Process completed. Report generated in #{output_file}."
end

# Check for command-line arguments and execute the main function
main(ARGV[0], ARGV[1])

