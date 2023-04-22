using Revise
using GeoPhyInv
using Statistics
using Test
using LossFunctions

@init_parallel_stencil(2, false, Float32, 2)


pa = SeisForwExpt(:acou_homo2D)

for field in [:p, :vx, :vz]
    println("############ Testing Backprop for source type ", field)

    srcwav = pa[:srcwav][1]
    GeoPhyInv.update!(srcwav[1], [field])

    for src_types in [[1, -1], [2, -2]]
        pa.c.backprop_flag = :save # do backpropagation

        GeoPhyInv.update!(pa, [srcwav], [src_types[1]])

        update!(pa)
        rec1 = deepcopy(pa.c.data[1])
        rec1 = rec1[1].d[1]
        rec1 ./= std(rec1)

        # change source flag and update wavelets in pa
        GeoPhyInv.update!(pa, [srcwav], [src_types[2]])
        pa.c.backprop_flag = :force # do backpropagation

        update!(pa)
        rec2 = deepcopy(pa.c.data[1])

        # time reverse
        reverse!(rec2)
        rec2 = rec2[1].d[1]
        rec2 ./= std(rec2)

        # compare results
        # compute L2dist
        @show err = value(L2DistLoss(), rec1, rec2, AggMode.Mean())

        # desired accuracy?
        @test err < 1e-10
    end
end
