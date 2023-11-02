#!/usr/bin/env ruby
# == NAME
# fastngsadmix2ethnicity
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A skeleton script for Ruby
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'
require 'rubyXL'
require 'rubyXL/convenience_methods/cell'
require 'rubyXL/convenience_methods/color'
require 'rubyXL/convenience_methods/font'
require 'rubyXL/convenience_methods/workbook'
require 'rubyXL/convenience_methods/worksheet'

### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
    puts opts
    exit
}

opts.parse! 

abort "Must provide output name (--outfile)" unless options.outfile

files = Dir["*.txt"].sort

color = {
	"even" => "EBEBEB",
	"uneven" => "D6D6D6"
}

headers = []
data = []

files.each do |file|

    sample = file.split(".")[0]
    calls = {}
    IO.readlines(file).each do |line|
        line.strip!
        ethnic,perc = line.split(":").collect{|l| l.strip }
        perc = perc.to_f
        calls[ethnic] = perc
        headers << ethnic unless headers.include?(ethnic)
    end
    
    data << { "sample" => sample, "calls" => calls   }

end

headers.sort!

workbook = RubyXL::Workbook.new
sheet = workbook.worksheets[0]
sheet.sheet_name = "Ethnicity predictions"

row = 0
col = 0

# Construct the header row
sheet.add_cell(row,col,"Sample")
col += 1
headers.each do |h|
    sheet.add_cell(row,col,h)
    sheet.sheet_data[0][col].change_font_bold(true)
    col += 1
end

# iterate over all data objects and fill the rows
data.each do |d|
    row += 1
    col = 0
    sheet.add_cell(row,col,d["sample"])
    dcalls = d["calls"]
    headers.each_with_index do |h,i|
        e = h
        p = 0
        if dcalls.has_key?(h)
            p = dcalls[h] 
        end
        sheet.add_cell(row,i+1,p)
        
    end
    row.even? ? bg = color["even"] : bg = color["uneven"]
    sheet.change_row_fill(row, bg)
end

workbook.write(options.outfile)