#ifndef MACROUTILS_HPP
#define MACROUTILS_HPP

/// Iterate over a row of data,
/// (optional) compute the physical position on the given direction
#define ITERATE_ROW(ix, min, max, compute_physics, phy)                                 \
    for(index_t ix = min; ix < max; ++ix) {                                             \
        Real phy = 0;                                                                   \
        if constexpr(compute_physics)                                                   \
            phy = real(ix + params.st_n##phy) * params.d##phy + params.origin##phy;

/// Iterate over a face of data, just a wrapper of nested call of ITERATE_ROW
#define ITERATE_FACE(ixA, minA, maxA, ixB, minB, maxB, compute_physics, phyA, phyB)     \
    ITERATE_ROW(ixA, minA, maxA, compute_physics, phyA)                                 \
        _Pragma("omp simd")                                                             \
        ITERATE_ROW(ixB, minB, maxB, compute_physics, phyB)

/// Iterate over a volume of data, just a wrapper of nested call of ITERATE_ROW
#define ITERATE_DOMAIN(ixA, minA, maxA, ixB, minB, maxB, ixC, minC, maxC, compute_physics, phyA, phyB, phyC)    \
    ITERATE_ROW(ixA, minA, maxA, compute_physics, phyA)                                                         \
        ITERATE_FACE(ixB, minB, maxB, ixC, minC, maxC, compute_physics, phyB, phyC)


/// Whether to skip the last index of ...
/// ... no one
#define NO_SKIP 0
/// ... x
#define SKIP_U 1
/// ... y
#define SKIP_V 2
/// ... z
#define SKIP_W 4

/// Iterate over velocity domain, allowing not to iterate on direction-last nodes
/// (skip is useful for optimization, since direction-last nodes are on physical domain, so their value is set by boundaries)
#define ITERATE_DOMAIN_VELOCITY(i, j, k, compute_physics, skip_C)                                           \
    ITERATE_DOMAIN(k, 0, ((params.isOnRight && (skip_C & SKIP_W)) ? params.loc_nZ - 1 : params.loc_nZ),     \
                   j, 0, ((params.isOnTop && (skip_C & SKIP_V)) ? params.loc_nY - 1 : params.loc_nY),       \
                   i, 0, ((skip_C & SKIP_U) ? params.loc_nX - 1 : params.loc_nX),                           \
                   compute_physics, Z, Y, X)


/// Iterate on pressure domain, since all pressure point are into domain skip is not allowed
#define ITERATE_DOMAIN_PRESSURE(i, j, k, compute_physics)                                   \
    ITERATE_DOMAIN(k, 0, params.loc_nZ,                                                     \
                   j, 0, params.loc_nY,                                                     \
                   i, 0, params.loc_nX,                                                     \
                   compute_physics, Z, Y, X)


/// Iterate over a face on plane XZ
#define ITERATE_XZ_FACE(i, k, compute_physics)                                              \
    ITERATE_FACE(k, 0, params.loc_nZ, i, 0, params.loc_nX, compute_physics, Z, X)

/// Iterate over a face on plane XY
#define ITERATE_XY_FACE(i, j, compute_physics)                                              \
    ITERATE_FACE(j, 0, params.loc_nY, i, 0, params.loc_nX, compute_physics, Y, X)

/// Iterate over a face on plane YZ
#define ITERATE_YZ_FACE(j, k, compute_physics)                                              \
    ITERATE_FACE(k, 0, params.loc_nZ, j, 0, params.loc_nY, compute_physics, Z, Y)


/// Iterate on "physical" domain, that is, it excludes all boundary points, since they are all on real boundary,
/// and their value is computed with exact functions
#define ITERATE_DOMAIN_PHYSICAL(i, j, k, compute_physics)                                   \
    ITERATE_DOMAIN(k, 1, params.phy_nZ - 1,                                                 \
                   j, 1, params.phy_nY - 1,                                                 \
                   i, 1, params.phy_nX - 1,                                                 \
                   compute_physics, Z, Y, X)

/// Iterate on "physical" face on XY plane, that is, it excludes all boundary points, since they are all on real boundary,
#define ITERATE_XY_PHYSICAL_FACE(i, j, compute_physics)                                     \
    ITERATE_FACE(j, 1, params.phy_nY - 1, i, 1, params.phy_nX - 1, compute_physics, Y, X)

/// Iterate on "physical" domain, that is, it excludes all boundary points, since they are all on real boundary,
#define ITERATE_XZ_PHYSICAL_FACE(i, k, compute_physics)                                     \
    ITERATE_FACE(k, 1, params.phy_nZ - 1, i, 1, params.phy_nX - 1, compute_physics, Z, X)

/// Close row iteration
#define ITERATE_ROW_END() }

/// Close face iteration
#define ITERATE_FACE_END()  \
    }                       \
}

/// Close domain iterations
#define ITERATE_DOMAIN_END()    \
        }                       \
    }                           \
}

#endif //MACROUTILS_HPP
