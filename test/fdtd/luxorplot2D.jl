function mlaxor(attrib_mod)
    b=GeoPhyInv.Fields(attrib_mod)
    cols = ["red","blue","black","green","blue2","orange","brown","tan1","blue4","grey","green4","purple4","orange4","brown4"]
    i=1
    @pdf begin
    xlim=(0,100)
    zlim=(0,100)
        for c in b
        mgrid=GeoPhyInv.get_mgrid(eval(:(GeoPhyInv.$c())),attrib_mod, range(xlim[1], stop=xlim[2], step=10), range(zlim[1], stop=zlim[2], step=10))
        println(mgrid)
         setcolor(cols[i])
         println(cols[i])
        setmode("overlap")
            for z in mgrid[2], x in mgrid[1]
            circle(Point(z, x),1,:fill) 
            # fontsize(2)
            # fontface("Georgia")
            # [Luxor.text(string(z),Point(z,x))]
        if (z==mgrid[2][1] && x==mgrid[1][1])
            fontsize(5)
            fontface("Georgia")
            [Luxor.text(string(c),Point(xlim[2]+10,zlim[1]+(i-1)*5))]
            end
        end
        i+=1
    end 
    end 600 600
    finish()
    preview()
end 