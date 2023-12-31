#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: Hybrid Cross OBJECT 
# author: Enrique Solera Navarro
#-----------------------------------------------------

class HybridCross
    attr_accessor :parent1, :parent2, :f2_wild, :f2_p1, :f2_p2, :f2_p1p2
  
    def initialize(params = {})
      @parent1 = params.fetch(:parent1, nil)
      @parent2 = params.fetch(:parent2, nil)
      @f2wild = params.fetch(:f2_wild, 0)

      @f2_wild = @f2wild.to_f #convert data to float
      @f2p1 = params.fetch(:f2_p1, 0)
      @f2_p1 = @f2p1.to_f
      @f2p2 = params.fetch(:f2_p2, 0)
      @f2_p2 = @f2p2.to_f
      @f2p1p2 = params.fetch(:f2_p1p2, 0)
      @f2_p1p2 = @f2p1p2.to_f

    
    end
  
    def calculate_chi_square
      total = @f2_wild + @f2_p1 + @f2_p2 + @f2_p1p2
      expected_wild = total * 9 / 16
      expected_p1 = total * 3 / 16
      expected_p2 = total * 3 / 16
      expected_p1p2 = total * 1 / 16
  
      chi_square = ((@f2_wild - expected_wild)**2 / expected_wild) +
                   ((@f2_p1 - expected_p1)**2 / expected_p1) +
                   ((@f2_p2 - expected_p2)**2 / expected_p2) +
                   ((@f2_p1p2 - expected_p1p2)**2 / expected_p1p2)
  
      chi_square
    end
  end
  