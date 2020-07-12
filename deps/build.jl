using Conda
pip = joinpath(Conda.BINDIR, "pip")
run(`$pip install ecmwf-api-client`)
