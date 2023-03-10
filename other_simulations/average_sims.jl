using Serialization
using CSV
using DataFrames

thread_count = Threads.nthreads()

for j in 1:11

    file_extend = "/result_r$(j).jls"
    
    filename = "data/thread1/result_r1.jls"
    data = deserialize(filename)
    shockdensity = zeros(size(data))

    for id in 1:thread_count
        filename = "data/thread$(id)" * file_extend
        data = deserialize(filename)

        shockdensity .+= data
    
    end

    shockdensity ./= thread_count

    filename_output = "data/result_r$(j).csv"
    data_out = DataFrame(shockdensity, :auto)

    CSV.write(filename_output, data_out)
end
