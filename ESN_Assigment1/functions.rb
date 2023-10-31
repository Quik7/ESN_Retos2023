#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: Functions 
# author: Enrique Solera Navarro
#-----------------------------------------------------
# The purpose of this script is that it contains various functions that handle the loading of data 
# from files, processing of data, and updating of genebank information. These functions are crucial 
# for the main program to execute its tasks.
#-----------------------------------------------------

def load_genes(file)
    genes = {}
    File.open(file).each_with_index do |line, index|
      next if index == 0  # Skip header
      fields = line.chomp.split("\t")
      gene_id, gene_name, mutant_phenotype = fields
      genes[gene_id] = Gene.new(
        gene_id: gene_id,
        gene_name: gene_name,
        mutant_phenotype: mutant_phenotype
      )
    end
    genes
  end
  
  def load_seed_stocks(file, genes)
    seed_stocks = {}
    File.open(file).each_with_index do |line, index|
      next if index == 0  # Skip header
      fields = line.chomp.split("\t")
      seed_stock_id, mutant_gene_id, last_planted, storage, grams_remaining = fields
      seed_stocks[seed_stock_id] = SeedStock.new(
        seed_stock_id: seed_stock_id,
        mutant_gene: genes[mutant_gene_id],
        last_planted: last_planted,
        storage: storage,
        grams_remaining: grams_remaining.to_i
      )
    end
    seed_stocks
  end
  
  def load_cross_data(file, genes, seed_stocks)
    crosses = []
    File.open(file).each_with_index do |line, index|
      next if index == 0  # Skip header
      fields = line.chomp.split("\t")
      parent1_id, parent2_id, f2_wild, f2_p1, f2_p2, f2_p1p2 = fields
      parent1 = seed_stocks[parent1_id].mutant_gene
      parent2 = seed_stocks[parent2_id].mutant_gene
      crosses << HybridCross.new(
        parent1: parent1,
        parent2: parent2,
        f2_wild: f2_wild,
        f2_p1: f2_p1,
        f2_p2: f2_p2,
        f2_p1p2: f2_p1p2
      )
    end
    crosses
  end
  
  def plant_seeds_and_update_genebank(seed_stocks, file)
    #checks if the file already exist and ask the user if he wants to delete it or not
    if File.exist?(file)
      print "The file #{file} already exist. Do you want to overwrite? (y/n): "
      response = STDIN.gets.chomp.downcase
      unless response == 'y'
        puts "Change the name of the output file."
        abort
        return
      end
      puts "The file #{file} has been updated."
    end
  
    File.open(file, 'w') do |f|
      f.puts "Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining"
      seed_stocks.each do |seed_stock_id, seed_stock|
        seed_stock.plant_seeds(7)
        f.puts [
          seed_stock.seed_stock_id,
          seed_stock.mutant_gene.gene_id,
          seed_stock.last_planted,
          seed_stock.storage,
          seed_stock.grams_remaining
        ].join("\t")
      end
    end
  
  end
  
  # This function processes each cross in the list of crosses, calculates the chi-square and determines whether the genes of the parents are genetically linked.
  # If the chi-square value is equal to or greater than 3.841 (critical value for a significance level of 0.05), then the genes are considered to be linked.
  # In that case, the genes are added as linked genes and a message is printed on the console indicating that the genes are linked.
  
  def calculate_chi_square_and_update_genes(crosses)
    crosses.each do |cross|
      chi_square = cross.calculate_chi_square
      if chi_square >= 3.841  # Critical value for a significance level of 0.05
        cross.parent1.add_linked_gene(cross.parent2)
        cross.parent2.add_linked_gene(cross.parent1)
        puts "Recording: #{cross.parent1.gene_name} is genetically linked to #{cross.parent2.gene_name} with chi-square score #{chi_square}"
      end
    end
  end
  

  # This function prints the final report of the genetically linked genes.
  # It iterates through each gene, and for each gene, it checks if it has any linked genes.
  # If it has linked genes, it prints the gene name and its linked genes.
  
  def print_final_report(genes)
    puts "Final Report:"
    genes.each do |gene_id, gene|
      next if gene.linked_genes.empty?
      linked_genes = gene.linked_genes.uniq.map(&:gene_name).join(', ')
      puts "#{gene.gene_name} is linked to #{linked_genes}"
    end
  end
  