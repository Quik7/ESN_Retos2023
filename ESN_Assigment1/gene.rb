class Gene
  attr_accessor :gene_id, :gene_name, :mutant_phenotype, :linked_genes

  def initialize(params = {})
    @gene_id = params.fetch(:gene_id, '')
    @gene_name = params.fetch(:gene_name, '')
    @mutant_phenotype = params.fetch(:mutant_phenotype, '')
    @linked_genes = []

    validate_gene_id
  end

  def validate_gene_id
    unless @gene_id.match?(/A[Tt]\d[Gg]\d\d\d\d\d/)
      puts "Error: the gene identifier #{@gene_id} is not in the correct format. It has to be ATXGXXXXX"
      exit
    end
  end

  def add_linked_gene(gene)
    unless gene.is_a? Gene
        puts "Error: the object provided is not a gene."
        abort
    end
    @linked_genes << gene
  end
end
