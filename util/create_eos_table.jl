@info "Create EOS table"

@info "Creating file \"eos.quantities\"..."

# regular_quantities = [1, 2, 12, 20, 21, 22, 23, 24]
# additional_quantities = []
# derivative_quantities = []
# pair_quantities = []
# quadruple_quantities = []
# microscopic_quantities = []
# error_quantities = []

# These are all the quantities (why not?)
regular_quantities = collect(1:24)
additional_quantities = collect(1:1)
derivative_quantities = collect(1:10)
pair_quantities = [10, 11, 4002, 0]
quadruple_quantities = [1]
microscopic_quantities = [10040, 11040, 10050, 11050]
error_quantities = collect(1:8)
output_format = 0               # 1=ASCII, 0=HDF5

open("eos.quantities", "w") do fh
    println(fh,
            "# number of regular, additional and derivative quantities (see table 7.1)")
    println(fh,
            join([length(regular_quantities), length(additional_quantities),
                  length(derivative_quantities)], " "))
    println(fh, "# indices of regular, additional and derivative quantities")
    println(fh,
            join([regular_quantities; additional_quantities;
                  derivative_quantities], " "))
    println(fh,
            "# number of pairs and quadruples for composition data (see table 3.2/3.3)")
    println(fh,
            join([length(pair_quantities), length(quadruple_quantities)], " "))
    println(fh, "# indices for pairs and quadruples for composition data")
    println(fh, join([pair_quantities; quadruple_quantities], " "))
    println(fh, "# number of microscopic quantities (see section 4.2.4)")
    println(fh, length(microscopic_quantities))
    println(fh, "# indices of microscopic quantities")
    println(fh, join(microscopic_quantities, " "))
    println(fh, "# number of error quantities")
    println(fh, length(error_quantities))
    println(fh, "# indices of error quantities")
    println(fh, join(error_quantities, " "))
    println(fh, "# format of output file, ASCII (1) or HDF5 (else)")
    println(fh, output_format)
    return nothing
end

@info "Creating file \"eos.parameters\"..."

# This selects all the points in the original table, presumably without interpolation (why not?)

file_t = split(read("eos.t", String))
log_t = parse(Int, file_t[1])
num_t = parse(Int, file_t[2])
min_t = parse(Float64, file_t[3])
max_t = parse(Float64, file_t[end])
@assert length(file_t) == num_t + 2

file_nb = split(read("eos.nb", String))
log_nb = parse(Int, file_nb[1])
num_nb = parse(Int, file_nb[2])
min_nb = parse(Float64, file_nb[3])
max_nb = parse(Float64, file_nb[end])
@assert length(file_nb) == num_nb + 2

file_yq = split(read("eos.yq", String))
log_yq = parse(Int, file_yq[1])
num_yq = parse(Int, file_yq[2])
min_yq = parse(Float64, file_yq[3])
max_yq = parse(Float64, file_yq[end])
@assert length(file_yq) == num_yq + 2

open("eos.parameters", "w") do fh
    println(fh, "# order of interpolation in first, second and third index")
    println(fh, "1 1 1")
    println(fh,
            "# calculation of beta-equilibrium (1: yes, else: no) and for given entropy (1: yes, else: no)")
    println(fh, "0 0")
    println(fh,
            "# tabulation scheme (0 = explicit listing, 1 = loops, see manual)")
    println(fh, "1")
    println(fh,
            "# parameter values (first, second and third index) depending on tabulation scheme")
    println(fh, join([min_t, min_nb, min_yq], " "))
    println(fh, join([max_t, max_nb, max_yq], " "))
    println(fh, join([num_t, num_nb, num_yq], " "))
    println(fh, join([log_t, log_nb, log_yq], " "))
    return nothing
end

@info "Generate EOS file..."
run(pipeline(`echo 3`, `/Users/eschnett/src/codehdf5/compose`))

# h5repack -v -f SHUF -f GZIP=9 eoscompose.h5 eoscompose_compressed.h5

@info "Done."
