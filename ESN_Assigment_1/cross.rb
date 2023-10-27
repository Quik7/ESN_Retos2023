require 'math'
class cross
  attr_accessor :parent1, :parent2, :f2_wild, :f2_p1, :f2_p2, :f2_p1p2

  @@all_cross_object = Array.new # Class variable that will save in an array all the instances of HybridCross created
  
  def initialize(parent1, parent2, f2_wild, f2_p1, f2_p2, f2_p1p2)
    @parent1 = parent1
    @parent2 = parent2
    @f2_wild = f2_wild.to_i
    @f2_p1 = f2_p1.to_i
    @f2_p2 = f2_p2.to_i
    @f2_p1p2 = f2_p1p2.to_i
    # Everytime an HybridCross object is initialized we add it to the array that contains all the instances of this object
    @@all_cross_object << self 

end

end