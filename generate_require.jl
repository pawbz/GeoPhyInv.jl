
import Pkg
const PT = Pkg.Types

Pkg.activate(pwd())             # current directory as the project
ctx = PT.Context()
pkg = ctx.env.pkg
if pkg ≡ nothing
    @error "Not in a package, I won't generate REQUIRE."
    exit(1)
else
    @info "found package" pkg = pkg
end

deps = PT.get_deps(ctx)
non_std_deps = sort(collect(setdiff(keys(deps), values(ctx.stdlibs))))

open("REQUIRE", "w") do io
    for d in non_std_deps
        println(io, d)
        @info "listing $d"
    end
end
