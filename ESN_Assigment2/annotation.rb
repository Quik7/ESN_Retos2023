# Class UsoGeneral is used to manage a collection of annotations.
# Annotations are stored based on their type and identifier.
class UsoGeneral
  # Initializes a new UsoGeneral instance.
  def initialize
    @annotations = {}
  end

  # Adds an annotation to the collection.
  # @param type [Symbol, String] the type of the annotation.
  # @param identifier [Symbol, String] the identifier for the annotation.
  # @param value [Object] the value of the annotation.
  def add_annotation(type, identifier, value)
    @annotations[type] ||= {}
    @annotations[type][identifier] = value
  end

  # Retrieves annotations of a specific type.
  # @param type [Symbol, String] the type of annotations to retrieve.
  # @return [Hash] a hash of annotations for the given type.
  def get_annotations(type)
    @annotations[type] || {}
  end

  # Returns all annotations.
  # @return [Hash] a hash containing all annotations.
  def all_annotations
    @annotations
  end
end

  