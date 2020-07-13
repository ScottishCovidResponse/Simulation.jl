using Conda
using Pkg
pip = joinpath(Conda.BINDIR, "pip")
run(`$pip install ecmwf-api-client`)
Pkg.build("ORCA")
