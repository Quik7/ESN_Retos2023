class Gene
    @@genes = []
    attr_reader :gene_id, :gene_name, :mutant_phenotype

    filename = ARGV[0]

    unless filename
        puts "ERROR: proper execution of this script is: ruby main.rb gene_information.tsv  seed_stock_data.tsv  cross_data.tsv outputfile.tsv"
        abort
    end


    def initialize(gene_id)
      if gene_id =~ /A[Tt]\d[Gg]\d\d\d\d\d/   # Check if the gene_id matches the desired format
        @gene_id = gene_id
      else
        puts "ERROR: Invalid Gene Identifier format '#{gene_id}'. Expected format: A[Tt]\\d[Gg]\\d\\d\\d\\d\\d"
        abort
      end

      @@genes << self

    end

    def self.all_genes
      return @@genes
    end
end