require 'csv'
require 'date'
class SeedStock
    # Read the seed_stock_data.tsv file
    seed_stock_data = CSV.read("seed_stock_data.tsv", col_sep: "\t")
    header = seed_stock_data.shift

    # Get today's date
    today = Date.today.strftime("%-m/%-d/%Y")

    # Update the seed stock data
    warnings = []
    seed_stock_data.each do |row|
        grams_remaining = row[4].to_i
        grams_remaining -= 7
        
        if grams_remaining <= 0
            warnings << "Warning: Seed stock #{row[0]} has run out of seeds!"
            grams_remaining = 0
        end
    
        row[4] = grams_remaining.to_s
        row[2] = today
    end

    # Display the warnings
    warnings.each { |warning| puts warning }

    # Write the updated data to a new file
    CSV.open("updated_seed_stock_data.tsv", "wb", col_sep: "\t") do |csv|
    csv << header
    seed_stock_data.each { |row| csv << row }
    end

end