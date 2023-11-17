class Annotation
    attr_accessor :annotations
  
    def initialize
      @annotations = {} # { type => { id => name } }
    end
  
    def add_annotation(type, id, name)
      @annotations[type] ||= {}
      @annotations[type][id] = name
    end
  end
  