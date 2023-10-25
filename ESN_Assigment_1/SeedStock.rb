require 'date'

class SeedStock
    #Seed_Stock	Mutant_Gene_ID	Last_Planted Storage Grams_Remainin
    attr_accessor :seed_stock, :mutant_gene_id, :last_planted, :storage, :grams_remaining
    
    filename = ARGV[0]
    outputname = ARGV[1]

    unless filename and outputname
        puts "ERROR: proper execution of this script is: ruby main.rb gene_information.tsv  seed_stock_data.tsv  cross_data.tsv outputfile.tsv"
        abort
    end

    def initialize(seed_stock, mutant_gene_id, last_planted, storage, grams_remaining)
        @seed_stock = seed_stock
        @mutant_gene_id = mutant_gene_id
        #@gene = Gene.new(mutant_gene_id)
        #abort "gene must be a gene object" unless @mutant_gene_id.is_a? Gene
        @last_planted = last_planted
        @storage = storage
        @grams_remaining = grams_remaining.to_i #converts to integer
    end

    # Plant seeds and update the grams remaining
    def plant(grams)
        @grams_remaining -= grams
        # https://ruby-doc.org/stdlib-2.5.1/libdoc/date/rdoc/Date.html#method-c-today
        @last_planted = Date.today.strftime('%d/%m/%Y')
        if @grams_remaining <= 0 #sanity check
        puts "WARNING: Seed stock #{@seed_stock} has been completely used up!"
        @grams_remaining = 0
        end
    end

    # Convert the object to its string representation
    def to_s
        "#{@seed_stock}\t#{@mutant_gene_id}\t#{@last_planted}\t#{@storage}\t#{@grams_remaining}"
    end

    # Read the seed stock data from file
    seed_stocks = []
    File.open(filename, 'r').each_with_index do |line, index|
    # Skip the header
        next if index == 0
    
        data = line.chomp.split("\t")
        seed_stocks << SeedStock.new(*data)
    end

    # Plant 7 grams of seeds from each record
    seed_stocks.each do |stock|
        stock.plant(7)
    end

    # Write the updated data to a new file
    File.open(outputname, 'w') do |file|
        file.puts "Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining"
        seed_stocks.each do |stock|
            file.puts stock.to_s
        end
    end
end