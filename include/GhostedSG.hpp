#ifndef AEROHPC_A_GHOSTEDSG_H
#define AEROHPC_A_GHOSTEDSG_H

#include "StaggeredGrid.hpp"

template<Addressing_T A>
class GhostedSG : public StaggeredGrid<A> {
public:

    GhostedSG() = delete;

    /**
     * Construct a staggered grid with ghost nodes,
     * @param nodes number of non-ghost nodes
     * @param ghosts_p_axis_p_direction number of ghost nodes,
     *                                  3 couples (one per axis),
     *                                  each couple define the number of ghost nodes (before, after) non-ghost nodes
     */
    GhostedSG(const std::array<index_t, 3> &nodes,
              const std::array<index_t , 2> &ghosts_x,
              const std::array<index_t , 2> &ghosts_y,
              const std::array<index_t , 2> &ghosts_z)
            : StaggeredGrid<A>({nodes[0] + ghosts_x[0] + ghosts_x[1],
                                nodes[1] + ghosts_y[0] + ghosts_y[1],
                                nodes[2] + ghosts_z[0] + ghosts_z[1]}),
              _px(ghosts_x[0]),
              _py(ghosts_y[0]),
              _pz(ghosts_z[0]) {
        this->_nodes = nodes;
    };

    /**
     * Construct a staggered grid with ghost nodes, number of ghost is the same in axis-independent
     */
    GhostedSG(const std::array<index_t , 3> &nodes, index_t ghosts)
            : GhostedSG(nodes,
                        {ghosts, ghosts}, {ghosts, ghosts}, {ghosts, ghosts}) {};

#pragma inline

    /**
     * Operator that accesses the memory using a 3D view of the object
     */
    Real &operator()(Component c, index_t i, index_t j, index_t k) override;

#pragma inline

    /**
     * Operator that accesses 3D view for read-only operations
     */
    const Real &operator()(Component c, index_t i, index_t j, index_t k) const override;

private:
    /**
     * Padding induced by ghost point on axis x
     */
    index_t _px;

    /**
     * Padding induced by ghost point on axis y
     */
    index_t _py;

    /**
     * Padding induced by ghost point on axis z
     */
    index_t _pz;
};

#endif //AEROHPC_A_GHOSTEDSG_H
