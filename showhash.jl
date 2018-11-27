using JLD

@load "kmers.jld"

#@show kmers

for (key, counts) in pairs(kmers)
    if length(counts) < sum(counts) # there's at least 1 kmer in every file
        print("$key\t")
        println(join(counts, "\t"))
    end
end
