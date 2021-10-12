# borrowed 2-point stencil macros from ParallelStencil.FiniteDifferences; extended to higher-order differences

import ..ParallelStencil: INDICES, WITHIN_DOC
iz, iy, ix = INDICES[1], INDICES[2], INDICES[3]
@static if (FD_ORDER == 2)
    izi, iyi, ixi = :($iz + 1), :($iy + 1), :($ix + 1)
elseif (FD_ORDER == 4)
    izi, iyi, ixi = :($iz + 3), :($iy + 3), :($ix + 3)
end

@static if (FD_ORDER == 2)

    macro within(macroname::String, A::Symbol)
        if macroname == "@all"
            esc(:($iz <= size($A, 1) && $iy <= size($A, 2) && $ix <= size($A, 3)))
        elseif macroname == "@inn"
            esc(
                :(
                    $iz <= size($A, 1) - 2 &&
                    $iy <= size($A, 2) - 2 &&
                    $ix <= size($A, 3) - 2
                ),
            )
        else
            error(
                "unkown macroname: $macroname. If you want to add your own assignement macros, overwrite the macro 'within(macroname::String, A::Symbol)'; to still use the exising macro within as well call ParallelStencil.FiniteDifferences{1|2|3}D.@within(macroname, A) at the end.",
            )
        end
    end

    macro d_za(A::Symbol)
        esc(:($A[$iz+1, $iy, $ix] - $A[$iz, $iy, $ix]))
    end
    macro d_ya(A::Symbol)
        esc(:($A[$iz, $iy+1, $ix] - $A[$iz, $iy, $ix]))
    end
    macro d_xa(A::Symbol)
        esc(:($A[$iz, $iy, $ix+1] - $A[$iz, $iy, $ix]))
    end
    macro d_zi(A::Symbol)
        esc(:($A[$iz+1, $iyi, $ixi] - $A[$iz, $iyi, $ixi]))
    end
    macro d_yi(A::Symbol)
        esc(:($A[$izi, $iy+1, $ixi] - $A[$izi, $iy, $ixi]))
    end
    macro d_xi(A::Symbol)
        esc(:($A[$izi, $iyi, $ix+1] - $A[$izi, $iyi, $ix]))
    end
elseif (FD_ORDER == 4)

    macro within(macroname::String, A::Symbol)
        if macroname == "@all"
            esc(:($iz <= size($A, 1) && $iy <= size($A, 2) && $ix <= size($A, 3)))
        elseif macroname == "@inn"
            esc(
                :(
                    $iz <= size($A, 1) - 6 &&
                    $iy <= size($A, 2) - 6 &&
                    $ix <= size($A, 3) - 6
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
                $A[$iz+2, $iy, $ix] * 27.0 - $A[$iz+1, $iy, $ix] * 27.0 +
                $A[$iz+3, $iy, $ix] - $A[$iz, $iy, $ix]
            ),
        )
    end
    macro d_ya(A::Symbol)
        esc(
            :(
                $A[$iz, $iy+2, $ix] * 27.0 - $A[$iz, $iy+1, $ix] * 27.0 +
                $A[$iz, $iy+3, $ix] - $A[$iz, $iy, $ix]
            ),
        )
    end
    macro d_xa(A::Symbol)
        esc(
            :(
                $A[$iz, $iy, $ix+2] * 27.0 - $A[$iz, $iy, $ix+1] * 27.0 +
                $A[$iz, $iy, $ix+3] - $A[$iz, $iy, $ix]
            ),
        )
    end
    macro d_zi(A::Symbol)
        esc(
            :(
                $A[$iz+2, $iyi, $ixi] * 27.0 - $A[$iz+1, $iyi, $ixi] * 27.0 +
                $A[$iz+3, $iyi, $ixi] - $A[$iz, $iyi, $ixi]
            ),
        )
    end
    macro d_yi(A::Symbol)
        esc(
            :(
                $A[$izi, $iy+2, $ixi] * 27.0 - $A[$izi, $iy+1, $ixi] * 27.0 +
                $A[$izi, $iy+3, $ixi] - $A[$izi, $iy, $ixi]
            ),
        )
    end
    macro d_xi(A::Symbol)
        esc(
            :(
                $A[$izi, $iyi, $ix+2] * 27.0 - $A[$izi, $iyi, $ix+1] * 27.0 +
                $A[$izi, $iyi, $ix+3] - $A[$izi, $iyi, $ix]
            ),
        )
    end
end
macro all(A::Symbol)
    esc(:($A[$iz, $iy, $ix]))
end
macro inn(A::Symbol)
    esc(:($A[$izi, $iyi, $ixi]))
end
macro av(A::Symbol)
    esc(
        :(
            (
                $A[$iz, $iy, $ix] +
                $A[$iz+1, $iy, $ix] +
                $A[$iz+1, $iy+1, $ix] +
                $A[$iz+1, $iy+1, $ix+1] +
                $A[$iz, $iy+1, $ix+1] +
                $A[$iz, $iy, $ix+1] +
                $A[$iz+1, $iy, $ix+1] +
                $A[$iz, $iy+1, $ix]
            ) * 0.125
        ),
    )
end
macro av_zi(A::Symbol)
    esc(:(($A[$iz, $iyi, $ixi] + $A[$iz+1, $iyi, $ixi]) * 0.5))
end
macro av_yi(A::Symbol)
    esc(:(($A[$izi, $iy, $ixi] + $A[$izi, $iy+1, $ixi]) * 0.5))
end
macro av_xi(A::Symbol)
    esc(:(($A[$izi, $iyi, $ix] + $A[$izi, $iyi, $ix+1]) * 0.5))
end
macro av_yzi(A::Symbol)
    esc(
        :(
            (
                $A[$iz, $iy, $ixi] +
                $A[$iz+1, $iy, $ixi] +
                $A[$iz, $iy+1, $ixi] +
                $A[$iz+1, $iy+1, $ixi]
            ) * 0.25
        ),
    )
end
macro av_xzi(A::Symbol)
    esc(
        :(
            (
                $A[$iz, $iyi, $ix] +
                $A[$iz+1, $iyi, $ix] +
                $A[$iz, $iyi, $ix+1] +
                $A[$iz+1, $iyi, $ix+1]
            ) * 0.25
        ),
    )
end
macro av_xyi(A::Symbol)
    esc(
        :(
            (
                $A[$izi, $iy, $ix] +
                $A[$izi, $iy+1, $ix] +
                $A[$izi, $iy, $ix+1] +
                $A[$izi, $iy+1, $ix+1]
            ) * 0.25
        ),
    )
end
