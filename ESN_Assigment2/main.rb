# Check the input file

unless ARGV[0] && ARGV[1] # We check user provides the gene list file
    abort "USAGE: main.rb genelist.txt output.txt"
end

unless File.exists?(ARGV[0]) # We check the given file exists
    abort "Error: File #{ARGV[0]} does not exist"
end

#checks if the file already exist and ask the user if he wants to delete it or not
if File.exist?(ARGV[1])
    print "The file #{ARGV[1]} already exist. Do you want to overwrite? (y/n): "
    response = STDIN.gets.chomp.downcase
    unless response == 'y'
      puts "Change the name of the output file."
      abort
      return
    end
    puts "The file #{ARGV[1]} has been updated."
  end