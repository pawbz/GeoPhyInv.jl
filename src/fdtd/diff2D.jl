# borrowed 2-point stencil macros from ParallelStencil.FiniteDifferences; extended to higher-order differences

import ..ParallelStencil: INDICES, WITHIN_DOC
iz, ix = INDICES[1], INDICES[2]
@static if (_fd.order == 2)
    izi, ixi = :($iz + 1), :($ix + 1)
elseif (_fd.order == 4)
    izi, ixi = :($iz + 3), :($ix + 3)
elseif (_fd.order == 6)
    izi, ixi = :($iz + 5), :($ix + 5)
end

@static if (_fd.order == 2)

    macro within(macroname::String, A::Symbol)
        if macroname == "@all"
            esc(:($iz <= size($A, 1) && $ix <= size($A, 2)))
        elseif macroname == "@inn"
            esc(
                :(
                    $iz <= size($A, 1) - 2 &&
                    $ix <= size($A, 2) - 2
                ),
            )
        else
            error(
                "unkown macroname: $macroname. If you want to add your own assignement macros, overwrite the macro 'within(macroname::String, A::Symbol)'; to still use the exising macro within as well call ParallelStencil.FiniteDifferences{1|2|3}D.@within(macroname, A) at the end.",
            )
        end
    end

    macro d_za(A::Symbol)
        esc(:($A[$iz+1, $ix] - $A[$iz, $ix]))
    end
    macro d_xa(A::Symbol)
        esc(:($A[$iz, $ix+1] - $A[$iz, $ix]))
    end
    macro d_zi(A::Symbol)
        esc(:($A[$iz+1, $ixi] - $A[$iz, $ixi]))
    end
    macro d_xi(A::Symbol)
        esc(:($A[$izi, $ix+1] - $A[$izi, $ix]))
    end
    # dummy y
elseif (_fd.order == 4)

    macro within(macroname::String, A::Symbol)
        if macroname == "@all"
            esc(:($iz <= size($A, 1) && $ix <= size($A, 2)))
        elseif macroname == "@inn"
            esc(
                :(
                    $iz <= size($A, 1) - 6 &&
                    $ix <= size($A, 2) - 6
                ),
            )
        else
            error(
                "unkown macroname: $macroname. If you want to add your own assignement macros, overwrite the macro 'within(macroname::String, A::Symbol)'; to still use the exising macro within as well call ParallelStencil.FiniteDifferences{1|2|3}D.@within(macroname, A) at the end.",
            )
        end
    end

    macro d_za(A::Symbol)
        esc(
            :(
                $A[$iz+2,  $ix] * 27.0 - $A[$iz+1,  $ix] * 27.0 +
                $A[$iz+3,  $ix] - $A[$iz,  $ix]
            ),
        )
    end
    macro d_xa(A::Symbol)
        esc(
            :(
                $A[$iz,  $ix+2] * 27.0 - $A[$iz,  $ix+1] * 27.0 +
                $A[$iz,  $ix+3] - $A[$iz,  $ix]
            ),
        )
    end
    macro d_zi(A::Symbol)
        esc(
            :(
                $A[$iz+2, $ixi] * 27.0 - $A[$iz+1, $ixi] * 27.0 +
                $A[$iz+3, $ixi] - $A[$iz, $ixi]
            ),
        )
    end
    macro d_xi(A::Symbol)
        esc(
            :(
                $A[$izi, $ix+2] * 27.0 - $A[$izi, $ix+1] * 27.0 +
                $A[$izi, $ix+3] - $A[$izi, $ix]
            ),
        )
    end
elseif (_fd.order == 6)
    macro within(macroname::String, A::Symbol)
        if macroname == "@all"
            esc(:($iz <= size($A, 1) && $ix <= size($A, 2)))
        elseif macroname == "@inn"
            esc(
                :(
                    $iz <= size($A, 1) - 10 &&
                    $ix <= size($A, 2) - 10
                ),
            )
        else
            error(
                "unkown macroname: $macroname. If you want to add your own assignement macros, overwrite the macro 'within(macroname::String, A::Symbol)'; to still use the exising macro within as well call ParallelStencil.FiniteDifferences{1|2|3}D.@within(macroname, A) at the end.",
            )
        end
    end

    macro d_za(A::Symbol)
        esc(
            :(
                $A[$iz+3,  $ix] * 2250.0 - $A[$iz+2,  $ix] * 2250.0 +
                $A[$iz+1,  $ix] * 125.0 - $A[$iz+4,  $ix] * 125.0 +
                $A[$iz+5,  $ix] * 9.0 - $A[$iz,  $ix] * 9.0
            ),
        )
    end
    macro d_xa(A::Symbol)
        esc(
            :(
                $A[$iz,  $ix+3] * 2250.0 - $A[$iz,  $ix+2] * 2250.0 +
                $A[$iz,  $ix+1] * 125.0 - $A[$iz,  $ix+4] * 125.0 +
                $A[$iz,  $ix+5] * 9.0 - $A[$iz,  $ix] * 9.0
            ),
        )
    end
    macro d_zi(A::Symbol)
        esc(
            :(
                $A[$iz+3,  $ixi] * 2250.0 - $A[$iz+2,  $ixi] * 2250.0 +
                $A[$iz+1,  $ixi] * 125.0 - $A[$iz+4,  $ixi] * 125.0 +
                $A[$iz+5,  $ixi] * 9.0 - $A[$iz,  $ixi] * 9.0
            ),
        )
    end
    macro d_xi(A::Symbol)
        esc(
            :(
                $A[$izi,  $ix+3] * 2250.0 - $A[$izi,  $ix+2] * 2250.0 +
                $A[$izi,  $ix+1] * 125.0 - $A[$izi,  $ix+4] * 125.0 +
                $A[$izi,  $ix+5] * 9.0 - $A[$izi,  $ix] * 9.0
            ),
        )
    end
end

macro all(A::Symbol)
    esc(:($A[$iz,  $ix]))
end
macro inn(A::Symbol)
    esc(:($A[$izi, $ixi]))
end
macro av(A::Symbol)
    esc(
        :(
            (
                $A[$iz,  $ix] +
                $A[$iz+1,  $ix] +
                $A[$iz,  $ix+1] +
                $A[$iz+1,  $ix+1] 
            ) * 0.25
        ),
    )
end
macro av_zi(A::Symbol)
    esc(:(($A[$iz, $ixi] + $A[$iz+1, $ixi]) * 0.5))
end
macro av_xi(A::Symbol)
    esc(:(($A[$izi, $ix] + $A[$izi, $ix+1]) * 0.5))
end


macro d_yi(args...) end
macro av_yi(args...) end
macro d_ya(args...) end
macro av_xzi(args...) end
macro av_yzi(args...) end
macro av_xyi(args...) end