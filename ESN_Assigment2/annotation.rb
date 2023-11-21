class UsoGeneral
  def initialize
    @annotations = {}
  end

  def add_annotation(type, identifier, value)
    @annotations[type] ||= {}
    @annotations[type][identifier] = value
  end

  def get_annotations(type)
    @annotations[type] || {}
  end

  def all_annotations
    @annotations
  end
end
  