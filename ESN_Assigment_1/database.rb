#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: database, main script
# author: Enrique Solera Navarro
#-----------------------------------------------------------------
#Import the diferent modules needed
require './gene.rb'
require './seedStock.rb'
require './cross.rb'

# Input arguments
gene_file, stock_file, cross_file, update_file = ARGV

unless ARGV.length == 4 # We check user inputs the filenames correctly
    abort " USAGE: ruby database.rb gene_information.tsv seed_stock_data.tsv cross_data.tsv output_file.tsv"
end

