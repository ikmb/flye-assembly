#!/bin/env ruby

require 'optparse'
require 'ostruct'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
options.same = false
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-s","--[no-]same") {|argument| options.same = true } 
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

date = Time.now.strftime("%Y-%m-%d")
options.centre ? center = options.centre : center = "IKMB"

options.folder.gsub!(/\/$/, '')
movies = Dir["#{options.folder}/*.bam"].collect{|m| m.gsub(/\/$/, '') }

sample_name = "AssemblyProject1"

puts "ProjectID;Movie;MovieIndex"
movies.each do |movie|

	options.same ? project = sample_name : project = movie.split("/")[-1].split(".")[0]

	puts "#{project};#{movie};#{movie}.pbi"

end
	
