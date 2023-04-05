#!/usr/bin/env ruby
# == NAME
# script_skeleton.rb
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


### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-i","--infile", "=Infile","fastNGSadmix output file") {|argument| options.infile = argument }
opts.on("-t","--threshold", "=THRESHOLD","Reporting threshold (float)") {|argument| options.threshold = argument }
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse!

options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

options.threshold ? threshold = options.threshold : threshold = 0.1

population = {
    "GBR" => "British" ,
    "FIN" => "Finnish",
    "CHS" => "Han Chinese (South)",
    "PUR" => "Puerto Rican",
    "CDX" => "Chinese (Dai)",
    "CLM" => "Colombian",
    "IBS" => "Iberian/Spanish",
    "PEL" => "Peruvian (Lima)",
    "PJL" => "Punjabi (Lahore)",
    "KHV" => "Kinh (Vietnam)",
    "ACB" => "African Caribbean",
    "GWD" => "Gambian (Mandinka)",
    "ESN" => "Esan (Nigeria)",
    "BEB" => "Bengali (Bangladesh)",
    "MSL" => "Mende (Sierra Leone)",
    "STU" => "Sri Lankan Tamil",
    "ITU" => "Indian Telugu",
    "CEU" => "Western European",
    "YRI" => "Yoruba (Nigeria)",
    "CHB" => "Han Chinese",
    "JPT" => "Japanese",
    "LWK" => "Luhya (Kenya)",
    "ASW" => "African",
    "MXL" => "Mexican",
    "TSI" => "Toscany (Italy)",
    "GIH" => "Gujarati (India)"
}

lines = IO.readlines(options.infile)

header = lines.shift.split(" ")
values = lines.shift.split(" ")

values.each_with_index do |v,i|

	v = v.to_f
	if v > threshold
		pop = header[i]
		p_name = population[pop]
		output_stream.puts "#{p_name}: #{v}"
	end
end

output_stream.close
