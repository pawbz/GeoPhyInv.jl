
using Literate

freadme=joinpath(@__DIR__, "page1.jl")
#output_folder=joinpath(@__DIR__, "..","..","..","notebooks","modeling") 
output_folder=joinpath(@__DIR__, ".") 


#Literate.notebook(freadme, output_folder)
Literate.markdown(freadme, output_folder)

freadme=joinpath(@__DIR__, "page2.jl")
Literate.markdown(freadme, output_folder)
#Literate.notebook(freadme, output_folder)

