require './gene.rb'
require './seed_stock.rb'
require './hybrid_cross.rb'
require './functions.rb'

def main(gene_file, seed_stock_file, cross_data_file, new_stock_file)
    genes = load_genes(gene_file)
    seed_stocks = load_seed_stocks(seed_stock_file, genes)
    crosses = load_cross_data(cross_data_file, genes, seed_stocks)
    plant_seeds_and_update_genebank(seed_stocks, new_stock_file)
    calculate_chi_square_and_update_genes(crosses)
    print_final_report(genes)
end
  
  # Call to the main function with the program arguments
unless ARGV.length == 4
    puts "ERROR: Usage: ruby #{$0} gene_information.tsv seed_stock_data.tsv cross_data.tsv new_stock_file.tsv"
    exit
end
  
# The splat (*) operator converts the ARGV array into an argument list, allowing the ARGV elements to be passed as arguments to the main function.
main(*ARGV)
