class InteractionNetwork
    attr_accessor :genes, :kegg_pathways, :go_terms
  
    def initialize(genes)
      @genes = genes
      @kegg_pathways = {} # { kegg_id => kegg_name }
      @go_terms = {} # { go_id => go_name }
    end
  
    # Method to add KEGG pathway
    def add_kegg_pathway(kegg_id, kegg_name)
      @kegg_pathways[kegg_id] = kegg_name
    end
  
    # Method to add GO Term
    def add_go_term(go_id, go_name)
      @go_terms[go_id] = go_name
    end
end  