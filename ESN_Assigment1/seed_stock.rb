#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: Seed_Stock OBJECT
# author: Enrique Solera Navarro
#-----------------------------------------------------

require 'date'
class SeedStock
    attr_accessor :seed_stock_id, :mutant_gene, :last_planted, :storage, :grams_remaining
  
    def initialize(params = {})
      @seed_stock_id = params.fetch(:seed_stock_id, '')
      @mutant_gene = params.fetch(:mutant_gene, nil)
      @last_planted = params.fetch(:last_planted, '')
      @storage = params.fetch(:storage, '')
      @grams_remaining = params.fetch(:grams_remaining, 0)
    end

    def validate_params
        unless @mutant_gene.is_a? Gene
          raise "Error: the provided gene is invalid."
          abort
        end
        unless @grams_remaining.is_a? Integer
          raise "Error: the number of remaining seeds must be a whole number."
          abort
        end
    end
  
    def plant_seeds(grams)
      @grams_remaining -= grams
      #https://ruby-doc.org/stdlib-2.5.1/libdoc/date/rdoc/Date.html#method-c-today
      @last_planted = Date.today.strftime('%d/%m/%Y')
      if @grams_remaining <= 0
        @grams_remaining = 0
        puts "WARNING: we have run out of Seed Stock #{@seed_stock_id}"
      end
    end
  end
  