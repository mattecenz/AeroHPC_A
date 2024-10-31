#ifndef AEROHPC_A_OPERATORS_H
#define AEROHPC_A_OPERATORS_H

// Sorry for the quite unreadable code :))
// Many parts are commented 'cause I'm still working on the generalization
// of the boundary conditions

// In case of redundant classes pls tell me

class VectorialFunction
{
public:
    virtual Real u(Real x, Real y, Real z){return 0;}
    virtual Real v(Real x, Real y, Real z){return 0;}
    virtual Real w(Real x, Real y, Real z){return 0;}
};



template<Addressing_T A>
class Operator
{
protected:

    enum BoundaryCase
    {
        // Face cases
        North, South, East, West, Up, Down,

        // Edge cases
        NE, NW, SE, SW, UN, US, UE, UW, DN, NS, DE, DW,

        // Vertex cases
        UNE, UNW, USE, USW, DNE, DNW, DSE, DSW
    };

    // each subarray contains the 3 components of the velocity
    // [ijk, i+1jk, i-1jk, ij+1k, ij-1k, ijk+1, ijk-1]

    // These are the methods to override in order to create a specific operator
    virtual void apply_locally(
        std::array< std::array<Real, 3>, 13 > &dst_array, 
        std::array< std::array<Real, 3>, 13 > &src_array){}
    
    virtual void apply_locally_face_boundary(
        std::array< std::array<Real, 3>, 14 > &dst_array, 
        std::array< std::array<Real, 3>, 14 > &src_array,
        BoundaryCase bcase){}
    
    virtual void apply_locally_edge_boundary(
        std::array< std::array<Real, 3>, 9 > &dst_array, 
        std::array< std::array<Real, 3>, 9 > &src_array,
        BoundaryCase bcase){}
    
    virtual void apply_locally_vertex_boundary(
        std::array< std::array<Real, 3>, 10 > &dst_array, 
        std::array< std::array<Real, 3>, 10 > &src_array,
        BoundaryCase bcase){}
    
    


    void apply_locally(
        Grid<A> &dst_grid, const Grid<A> &src_grid,
        size_t i, size_t j, size_t k)
    {
        // Create a local copy of the input data
        std::array< std::array<Real, 3>, 13 > input;

        input[0] = {src_grid(Component::U, i, j, k), src_grid(Component::V, i, j, k), src_grid(Component::W, i, j, k)};
        input[1] = {src_grid(Component::U, i + 1, j, k), src_grid(Component::V, i + 1, j, k), src_grid(Component::W, i + 1, j, k)};
        input[2] = {src_grid(Component::U, i - 1, j, k), src_grid(Component::V, i - 1, j, k), src_grid(Component::W, i - 1, j, k)};
        input[3] = {src_grid(Component::U, i, j + 1, k), src_grid(Component::V, i, j + 1, k), src_grid(Component::W, i, j + 1, k)};
        input[4] = {src_grid(Component::U, i, j - 1, k), src_grid(Component::V, i, j - 1, k), src_grid(Component::W, i, j - 1, k)};
        input[5] = {src_grid(Component::U, i, j, k + 1), src_grid(Component::V, i, j, k + 1), src_grid(Component::W, i, j, k + 1)};
        input[6] = {src_grid(Component::U, i, j, k - 1), src_grid(Component::V, i, j, k - 1), src_grid(Component::W, i, j, k - 1)};

        // We need this other data in order to achieve 2nd order taylor prolongation of the components :((

        input[7] = {src_grid(Component::U, i + 1, j - 1, k), src_grid(Component::V, i + 1, j - 1, k), 0.0};
        input[8] = {src_grid(Component::U, i - 1, j + 1, k), src_grid(Component::V, i - 1, j + 1, k), 0.0};

        input[9] = {src_grid(Component::U, i + 1, j, k - 1), 0.0, src_grid(Component::W, i + 1, j, k - 1)};
        input[10] = {src_grid(Component::U, i - 1, j, k + 1), 0.0, src_grid(Component::W, i - 1, j, k + 1)};

        input[11] = {0.0, src_grid(Component::V, i, j + 1, k - 1), src_grid(Component::W, i, j + 1, k - 1)};
        input[12] = {0.0, src_grid(Component::V, i, j - 1, k + 1), src_grid(Component::W, i, j - 1, k + 1)};


        std::array< std::array<Real, 3>, 13 > output;

        apply_locally(output, input);

        dst_grid(Component::U, i, j, k) = output[0][0];
        dst_grid(Component::V, i, j, k) = output[0][1];
        dst_grid(Component::W, i, j, k) = output[0][2];
    }

    void apply_locally_face_boundary(
        Grid<A> &dst_grid, const Grid<A> &src_grid,
        size_t i, size_t j, size_t k, BoundaryCase bcase)
    {
        // Create a local copy of the input data
        std::array< std::array<Real, 3>, 14 > input;
        Component c;
        size_t comp;

        input[0] = {src_grid(Component::U, i, j, k), src_grid(Component::V, i, j, k), src_grid(Component::W, i, j, k)};
        input[1] = {src_grid(Component::U, i + 1, j, k), src_grid(Component::V, i + 1, j, k), src_grid(Component::W, i + 1, j, k)};
        input[2] = {src_grid(Component::U, i - 1, j, k), src_grid(Component::V, i - 1, j, k), src_grid(Component::W, i - 1, j, k)};
        input[3] = {src_grid(Component::U, i, j + 1, k), src_grid(Component::V, i, j + 1, k), src_grid(Component::W, i, j + 1, k)};
        input[4] = {src_grid(Component::U, i, j - 1, k), src_grid(Component::V, i, j - 1, k), src_grid(Component::W, i, j - 1, k)};
        input[5] = {src_grid(Component::U, i, j, k + 1), src_grid(Component::V, i, j, k + 1), src_grid(Component::W, i, j, k + 1)};
        input[6] = {src_grid(Component::U, i, j, k - 1), src_grid(Component::V, i, j, k - 1), src_grid(Component::W, i, j, k - 1)};


        // We need this other data in order to achieve 2nd order taylor prolongation of the components :((

        input[7] = {src_grid(Component::U, i + 1, j - 1, k), src_grid(Component::V, i + 1, j - 1, k), 0.0};
        input[8] = {src_grid(Component::U, i - 1, j + 1, k), src_grid(Component::V, i - 1, j + 1, k), 0.0};

        input[9] = {src_grid(Component::U, i + 1, j, k - 1), 0.0, src_grid(Component::W, i + 1, j, k - 1)};
        input[10] = {src_grid(Component::U, i - 1, j, k + 1), 0.0, src_grid(Component::W, i - 1, j, k + 1)};

        input[11] = {0.0, src_grid(Component::V, i, j + 1, k - 1), src_grid(Component::W, i, j + 1, k - 1)};
        input[12] = {0.0, src_grid(Component::V, i, j - 1, k + 1), src_grid(Component::W, i, j - 1, k + 1)};

        switch (bcase)
        {
        case BoundaryCase::West:
            input[13] = {src_grid(Component::U, i + 2, j, k), src_grid(Component::V, i + 2, j, k), src_grid(Component::W, i + 2, j, k)};
            c = Component::U;
            comp = 0;
            break;
        
        case BoundaryCase::South:
            input[13] = {src_grid(Component::U, i, j + 2, k), src_grid(Component::V, i, j + 2, k), src_grid(Component::W, i, j + 2, k)};
            c = Component::V;
            comp = 1;
            break;
        
        case BoundaryCase::Down:
            input[13] = {src_grid(Component::U, i, j, k + 2), src_grid(Component::V, i, j, k + 2), src_grid(Component::W, i, j, k + 2)};
            c = Component::W;
            comp = 2;
            break;
        
        case BoundaryCase::East:
            input[13] = {src_grid(Component::U, i - 2, j, k), src_grid(Component::V, i - 2, j, k), src_grid(Component::W, i - 2, j, k)};
            c = Component::U;
            comp = 0;
            break;

        case BoundaryCase::North:
            input[13] = {src_grid(Component::U, i, j - 2, k), src_grid(Component::V, i, j - 2, k), src_grid(Component::W, i, j - 2, k)};
            c = Component::V;
            comp = 1;
            break;

        case BoundaryCase::Up:
            input[13] = {src_grid(Component::U, i, j, k - 2), src_grid(Component::V, i, j, k - 2), src_grid(Component::W, i, j, k - 2)};
            c = Component::W;
            comp = 2;
            break;
        
        default:
            throw std::invalid_argument("This boundary type is not supported by this method;");


        }

        std::array< std::array<Real, 3>, 14 > output;

        apply_locally_face_boundary(output, input, bcase);

        dst_grid(c, i, j, k) = output[0][comp];
    }


    void apply_locally_vertex_boundary(
        Grid<A> &dst_grid, const Grid<A> &src_grid,
        size_t i, size_t j, size_t k, BoundaryCase bcase)
    {
        // Create a local copy of the input data
        std::array< std::array<Real, 3>, 10 > input;
        size_t i_, j_, k_;

        input[0] = {src_grid(Component::U, i, j, k), src_grid(Component::V, i, j, k), src_grid(Component::W, i, j, k)};
        input[1] = {src_grid(Component::U, i + 1, j, k), src_grid(Component::V, i + 1, j, k), src_grid(Component::W, i + 1, j, k)};
        input[2] = {src_grid(Component::U, i - 1, j, k), src_grid(Component::V, i - 1, j, k), src_grid(Component::W, i - 1, j, k)};
        input[3] = {src_grid(Component::U, i, j + 1, k), src_grid(Component::V, i, j + 1, k), src_grid(Component::W, i, j + 1, k)};
        input[4] = {src_grid(Component::U, i, j - 1, k), src_grid(Component::V, i, j - 1, k), src_grid(Component::W, i, j - 1, k)};
        input[5] = {src_grid(Component::U, i, j, k + 1), src_grid(Component::V, i, j, k + 1), src_grid(Component::W, i, j, k + 1)};
        input[6] = {src_grid(Component::U, i, j, k - 1), src_grid(Component::V, i, j, k - 1), src_grid(Component::W, i, j, k - 1)};


        switch (bcase)
        {
        case BoundaryCase::UNE:
            i_ = i - 2;
            j_ = j - 2;
            k_ = k - 2;
            break;

        default:
            throw std::invalid_argument("This boundary type is not supported by this method;");
        }

        input[7] = {src_grid(Component::U, i_, j_, k_), src_grid(Component::V, i_, j_, k_), src_grid(Component::W, i_, j_, k_)};

        std::array< std::array<Real, 3>, 10 > output;

        apply_locally_vertex_boundary(output, input, bcase);

        dst_grid(Component::U, i, j, k) = output[0][0];
        dst_grid(Component::V, i, j, k) = output[0][1];
        dst_grid(Component::W, i, j, k) = output[0][2];

    }


public:

    void apply(Grid<A> &dst_grid, Grid<A> &src_grid)
    {
        // Since the computation of each element of the dst_grid
        // is independent, intra-node parallelization could be
        // implemented, i.e. using openMP
        
        for (unsigned int i = 1; i < src_grid.nx - 2; ++i)
        {
            for (unsigned int j = 1; j < src_grid.ny - 2; ++j)
            {
                for (unsigned int k = 1; k < src_grid.nz - 2; ++k)
                {
                    apply_locally(dst_grid, src_grid, i, j, k);
                }
            }
        }
        // Boundary conditions:
        // west face : i = 0
        for (unsigned int j = 1; j < src_grid.ny - 2; ++j)
        {
            for (unsigned int k = 1; k < src_grid.nz - 2; ++k)
            {
                apply_locally_face_boundary(dst_grid, src_grid, 0, j, k, BoundaryCase::West);
            }
        }

        // South face : j = 0
        for (unsigned int i = 1; i < src_grid.nx - 2; ++i)
        {
            for (unsigned int k = 1; k < src_grid.nz - 2; ++k)
            {
                apply_locally_face_boundary(dst_grid, src_grid, i, 0, k, BoundaryCase::South);
            }
        }

        // Down face : k = 0;
        
        for (unsigned int i = 1; i < src_grid.nx - 2; ++i)
        {
            for (unsigned int j = 1; j < src_grid.ny - 2; ++j)
            {
                apply_locally_face_boundary(dst_grid, src_grid, i, j, 0, BoundaryCase::Down);
            }
        }
        

        // east face : i = nx - 2
        
        for (unsigned int j = 1; j < src_grid.ny - 2; ++j)
        {
            for (unsigned int k = 1; k < src_grid.nz - 2; ++k)
            {
                apply_locally_face_boundary(dst_grid, src_grid, src_grid.nx - 2, j, k, BoundaryCase::East);
            }
        }

        // North face : j = ny - 2
        for (unsigned int i = 1; i < src_grid.nx - 2; ++i)
        {
            for (unsigned int k = 1; k < src_grid.nz - 2; ++k)
            {
                apply_locally_face_boundary(dst_grid, src_grid, i, src_grid.ny - 2, k, BoundaryCase::North);
            }
        }
        
        // Up face : k = 0;
        
        for (unsigned int i = 1; i < src_grid.nx - 2; ++i)
        {
            for (unsigned int j = 1; j < src_grid.ny - 2; ++j)
            {
                apply_locally_face_boundary(dst_grid, src_grid, i, j, src_grid.nz - 2, BoundaryCase::Up);
            }
        }

        // Up north east : i = nx - 2, j = ny - 2, k = nz - 2
        // NOT WORKING YET
        //apply_locally_vertex_boundary(dst_grid, src_grid, src_grid.nx - 2, src_grid.ny - 2, src_grid.nz - 2, BoundaryCase::UNE);
    }

    virtual void attach_to_model(Grid<A> &model){}
    
};




template<Addressing_T A>
class OneStepOperator : public Operator<A>
{
protected:


    virtual void apply_locally(
        std::array< std::array<Real, 3>, 13 > &dst_array, 
        std::array< std::array<Real, 3>, 13 > &src_array)
    {
        // First let's interpolate the value of the velocity 
        // (using second order taylor expansion)
        Real v_expansion_i = 0.75 * src_array[0][1] + 
            0.375 * src_array[7][1] - 0.125 * src_array[8][1];
        
        Real u_expansion_j = 0.75 * src_array[0][0] + 
            0.375 * src_array[8][0] - 0.125 * src_array[7][0];

        Real w_expansion_i = 0.75 * src_array[0][2] + 
            0.375 * src_array[9][2] - 0.125 * src_array[10][2];

        Real u_expansion_k = 0.75 * src_array[0][0] + 
            0.375 * src_array[10][0] - 0.125 * src_array[9][0];

        Real v_expansion_k = 0.75 * src_array[0][1] + 
            0.375 * src_array[12][1] - 0.125 * src_array[11][1];
        
        Real w_expansion_j = 0.75 * src_array[0][2] + 
            0.375 * src_array[11][2] - 0.125 * src_array[12][2];

        // Then apply the values
        dst_array[0][0] = coefficients[2] * (
            src_array[0][0] * (src_array[1][0] - src_array[2][0]) +
            v_expansion_i * (src_array[3][0] - src_array[4][0]) +
            w_expansion_i * (src_array[5][0] - src_array[6][0])) +
            coefficients[1] * (
            src_array[1][0] + src_array[2][0] + src_array[3][0] + 
            src_array[4][0] + src_array[5][0] + src_array[6][0]) +
            coefficients[0] * src_array[0][0];
        
        dst_array[0][1] = coefficients[2] * (
            u_expansion_j * (src_array[1][1] - src_array[2][1]) +
            src_array[0][1] * (src_array[3][1] - src_array[4][1]) +
            w_expansion_j * (src_array[5][1] - src_array[6][1])) +
            coefficients[1] * (
            src_array[1][1] + src_array[2][1] + src_array[3][1] + 
            src_array[4][1] + src_array[5][1] + src_array[6][1]) +
            coefficients[0] * src_array[0][1];
        
        dst_array[0][2] = coefficients[2] * (
            u_expansion_k * (src_array[1][2] - src_array[2][2]) +
            v_expansion_k * (src_array[3][2] - src_array[4][2]) +
            src_array[0][2] * (src_array[5][2] - src_array[6][2])) +
            coefficients[1] * (
            src_array[1][2] + src_array[2][2] + src_array[3][2] + 
            src_array[4][2] + src_array[5][2] + src_array[6][2]) +
            coefficients[0] * src_array[0][2];
        
    }

    virtual void apply_locally_face_boundary(                           // I'm still working on this method
        std::array< std::array<Real, 3>, 14 > &dst_array, 
        std::array< std::array<Real, 3>, 14 > &src_array,
        Operator<A>::BoundaryCase bcase)
    {
        std::array<size_t, 12> indices;
        size_t component;

        switch (bcase)
        {
            case Operator<A>::BoundaryCase::West:
            indices[5] = 1;     // i + 1
            indices[4] = 2;     // i - 1
            indices[3] = 3;     // j + 1
            indices[2] = 4;     // j - 1
            indices[1] = 5;     // k + 1
            indices[0] = 6;     // k - 1
            indices[6] = 7;     // i + 1 j - 1
            indices[7] = 8;     // i - 1 j + 1
            indices[8] = 9;     // i + 1 k - 1
            indices[9] = 10;    // i - 1 k + 1
            indices[10] = 11;   // j + 1 k - 1
            indices[11] = 12;   // k - 1 k + 1
            component = 0;      // U
            break;
        }

        // First let's interpolate the value of the velocity 
        // (using second order taylor expansion)
        Real v_expansion_i = 0.75 * src_array[0][1] + 
            0.375 * src_array[indices[6]][1] - 0.125 * src_array[indices[7]][1];
        
        Real u_expansion_j = 0.75 * src_array[0][0] + 
            0.375 * src_array[indices[7]][0] - 0.125 * src_array[indices[6]][0];

        Real w_expansion_i = 0.75 * src_array[0][2] + 
            0.375 * src_array[indices[8]][2] - 0.125 * src_array[10][2];

        Real u_expansion_k = 0.75 * src_array[0][0] + 
            0.375 * src_array[10][0] - 0.125 * src_array[9][0];

        Real v_expansion_k = 0.75 * src_array[0][1] + 
            0.375 * src_array[12][1] - 0.125 * src_array[11][1];
        
        Real w_expansion_j = 0.75 * src_array[0][2] + 
            0.375 * src_array[11][2] - 0.125 * src_array[12][2];




    }


    const std::array<bool, N_COMPONENTS> component_mask;
    Real coefficients[3];
    Real boundary_coefficients[4];
    Real face_boundary_diagonal_value;
    Real vertex_boundary_diagonal_value;

public:

    // Constructor
    explicit OneStepOperator(Real diagonal_value, Real direct_neighbor_linear, Real direct_neighbor_nonlinear,
        Real a = 0, Real b = 0, Real c = 0, Real d = 0)
        : component_mask({true, true, true, false})
    {
        coefficients[0] = diagonal_value;
        coefficients[1] = direct_neighbor_linear;
        coefficients[2] = direct_neighbor_nonlinear;

        boundary_coefficients[0] = a;
        boundary_coefficients[1] = b;
        boundary_coefficients[2] = c;
        boundary_coefficients[3] = d;

        face_boundary_diagonal_value = (- 2 * 2 * direct_neighbor_linear) + c;
        vertex_boundary_diagonal_value = 3 * c;
    }

    // Some useful operators

    OneStepOperator<A> operator+(const OneStepOperator &obj)
    {
        OneStepOperator<A> res(
            coefficients[0] + obj.coefficients[0], 
            coefficients[1] + obj.coefficients[1], 
            coefficients[2] + obj.coefficients[2],
            boundary_coefficients[0] + obj.boundary_coefficients[0],
            boundary_coefficients[1] + obj.boundary_coefficients[1],
            boundary_coefficients[2] + obj.boundary_coefficients[2],
            boundary_coefficients[3] + obj.boundary_coefficients[3]);
        
        return res;
    }

    OneStepOperator<A> operator*(const Real &val)
    {
        OneStepOperator<A> res(
            coefficients[0] * val,
            coefficients[1] * val, 
            coefficients[2] * val,
            boundary_coefficients[0] * val,
            boundary_coefficients[1] * val,
            boundary_coefficients[2] * val,
            boundary_coefficients[3] * val);
        return res;
    }

    friend OneStepOperator<A> operator*(const Real &val, const OneStepOperator<A> &obj)
    {
        OneStepOperator<A> res(
            obj.coefficients[0] * val,
            obj.coefficients[1] * val, 
            obj.coefficients[2] * val,
            obj.boundary_coefficients[0] * val,
            obj.boundary_coefficients[1] * val,
            obj.boundary_coefficients[2] * val,
            obj.boundary_coefficients[3] * val);
        return res;
    }

};






// This class is only useful for debugging
template<Addressing_T A>
class IdentityOperator : public OneStepOperator<A>
{
public:
    IdentityOperator(Grid<A> &src)
        : OneStepOperator<A>::OneStepOperator(1.0, 0.0, 0.0)
    {}

protected:
    void apply_locally(
        std::array< std::array<Real, 3>, 13 > &dst_array, 
        std::array< std::array<Real, 3>, 13 > &src_array) override
    {
        dst_array[0][0] = src_array[0][0];
        dst_array[0][1] = src_array[0][1];
        dst_array[0][2] = src_array[0][2];
    }

    void apply_locally_face_boundary(
        std::array< std::array<Real, 3>, 14 > &dst_array, 
        std::array< std::array<Real, 3>, 14 > &src_array,
        OneStepOperator<A>::BoundaryCase bcase) override
    {
        dst_array[0][0] = src_array[0][0];
        dst_array[0][1] = src_array[0][1];
        dst_array[0][2] = src_array[0][2];
    }

    void apply_locally_vertex_boundary(
        std::array< std::array<Real, 3>, 10 > &dst_array, 
        std::array< std::array<Real, 3>, 10 > &src_array,
        OneStepOperator<A>::BoundaryCase bcase) override
    {
        dst_array[0][0] = src_array[0][0];
        dst_array[0][1] = src_array[0][1];
        dst_array[0][2] = src_array[0][2];
    }
};









template<Addressing_T A>
class LaplacianOperator : public OneStepOperator<A>
{
public:
    LaplacianOperator(Grid<A> &src)
        : OneStepOperator<A>::OneStepOperator(- 6 / (src.dx * src.dx), 1 / (src.dx * src.dx), 0.,
            - 1 / (5 * src.dx * src.dx), 2 / (src.dx * src.dx), - 5 / (src.dx * src.dx), 16 / (5 * src.dx * src.dx))
    {}

protected:
    void apply_locally(
        std::array< std::array<Real, 3>, 13 > &dst_array, 
        std::array< std::array<Real, 3>, 13 > &src_array) override
    {
        dst_array[0][0] = this->coefficients[0] * src_array[0][0] +
            this->coefficients[1] * (
            src_array[1][0] + src_array[2][0] + src_array[3][0] + 
            src_array[4][0] + src_array[5][0] + src_array[6][0]);
        
        dst_array[0][1] = this->coefficients[0] * src_array[0][1] +
            this->coefficients[1] * (
            src_array[1][1] + src_array[2][1] + src_array[3][1] + 
            src_array[4][1] + src_array[5][1] + src_array[6][1]);
        
        dst_array[0][2] = this->coefficients[0] * src_array[0][2] +
            this->coefficients[1] * (
            src_array[1][2] + src_array[2][2] + src_array[3][2] + 
            src_array[4][2] + src_array[5][2] + src_array[6][2]);
    }


    void apply_locally_face_boundary(
        std::array< std::array<Real, 3>, 14 > &dst_array, 
        std::array< std::array<Real, 3>, 14 > &src_array,
        OneStepOperator<A>::BoundaryCase bcase) override
    {
        // Remember that boundary_coefficients[0] = a, [1] = b, [2] = c, [3] = d
        // Also in the vector we have the data ordered in this way:
        // [ijk, i+1jk, i-1jk, ij+1k, ij-1k, ijk+1, ijk-1], and the last one is +2
        size_t indices[6];
        size_t component;
        
        switch (bcase)
        {
        case OneStepOperator<A>::BoundaryCase::West:
            indices[5] = 1;     // i + 1
            indices[4] = 2;     // i - 1
            indices[3] = 3;     // j + 1
            indices[2] = 4;     // j - 1
            indices[1] = 5;     // k + 1
            indices[0] = 6;     // k - 1
            component = 0;      // U
            break;
        
        case OneStepOperator<A>::BoundaryCase::South:
            indices[3] = 1;     // i + 1
            indices[2] = 2;     // i - 1
            indices[5] = 3;     // j + 1
            indices[4] = 4;     // j - 1
            indices[1] = 5;     // k + 1
            indices[0] = 6;     // k - 1
            component = 1;      // V
            break;
        
        case OneStepOperator<A>::BoundaryCase::Down:
            indices[1] = 1;     // i + 1
            indices[0] = 2;     // i - 1
            indices[3] = 3;     // j + 1
            indices[2] = 4;     // j - 1
            indices[5] = 5;     // k + 1
            indices[4] = 6;     // k - 1
            component = 2;      // W
            break;

        case OneStepOperator<A>::BoundaryCase::East:
            indices[4] = 1;     // i + 1
            indices[5] = 2;     // i - 1
            indices[3] = 3;     // j + 1
            indices[2] = 4;     // j - 1
            indices[1] = 5;     // k + 1
            indices[0] = 6;     // k - 1
            component = 0;      // U
            break;
        
        case OneStepOperator<A>::BoundaryCase::North:
            indices[3] = 1;     // i + 1
            indices[2] = 2;     // i - 1
            indices[4] = 3;     // j + 1
            indices[5] = 4;     // j - 1
            indices[1] = 5;     // k + 1
            indices[0] = 6;     // k - 1
            component = 1;      // V
            break;
        
        case OneStepOperator<A>::BoundaryCase::Up:
            indices[1] = 1;     // i + 1
            indices[0] = 2;     // i - 1
            indices[3] = 3;     // j + 1
            indices[2] = 4;     // j - 1
            indices[4] = 5;     // k + 1
            indices[5] = 6;     // k - 1
            component = 2;      // W
            break;
        
        }
        
        dst_array[0][component] = this->face_boundary_diagonal_value * src_array[0][component] +
            this->coefficients[1] * (
            src_array[indices[0]][component] + src_array[indices[1]][component] +
            src_array[indices[2]][component] + src_array[indices[3]][component]) +
            this->boundary_coefficients[3] * src_array[indices[4]][component] +
            this->boundary_coefficients[1] * src_array[indices[5]][component] +
            this->boundary_coefficients[0] * src_array[13][component];
        
        
        //dst_array[0][component] = - 3 * src_array[0][component];
    }

    void apply_locally_vertex_boundary(                         // Doesn't work yet
        std::array< std::array<Real, 3>, 10 > &dst_array, 
        std::array< std::array<Real, 3>, 10 > &src_array,
        OneStepOperator<A>::BoundaryCase bcase) override
    {
        // Remember that boundary_coefficients[0] = a, [1] = b, [2] = c, [3] = d
        // Also in the vector we have the data ordered in this way:
        // [ijk, i+1jk, i-1jk, ij+1k, ij-1k, ijk+1, ijk-1], and the last one is +2
        size_t indices[6];

        switch (bcase)
        {
        case OneStepOperator<A>::BoundaryCase::UNE:
            indices[4] = 1;     // i + 1
            indices[5] = 2;     // i - 1
            indices[2] = 3;     // j + 1
            indices[3] = 4;     // j - 1
            indices[0] = 5;     // k + 1
            indices[1] = 6;     // k - 1
            break;
        }
        
        dst_array[0][0] = this->vertex_boundary_diagonal_value * src_array[0][0] +
            this->boundary_coefficients[3] * (src_array[indices[4]][0] + src_array[indices[2]][0] + src_array[indices[0]][0]) +
            this->boundary_coefficients[1] * (src_array[indices[5]][0] + src_array[indices[3]][0] + src_array[indices[1]][0]) +
            this->boundary_coefficients[0] * src_array[13][0];
        /*
        dst_array[0][1] = this->vertex_boundary_diagonal_value * src_array[0][1] +
            this->boundary_coefficients[3] * (src_array[indices[5]][1] + src_array[indices[3]][1] + src_array[indices[1]][1]) +
            this->boundary_coefficients[1] * (src_array[indices[4]][1] + src_array[indices[2]][1] + src_array[indices[0]][1]) +
            this->boundary_coefficients[0] * src_array[7][1];
        
        dst_array[0][2] = this->vertex_boundary_diagonal_value * src_array[0][2] +
            this->boundary_coefficients[3] * (src_array[indices[5]][2] + src_array[indices[3]][2] + src_array[indices[1]][2]) +
            this->boundary_coefficients[1] * (src_array[indices[4]][2] + src_array[indices[2]][2] + src_array[indices[0]][2]) +
            this->boundary_coefficients[0] * src_array[7][2];
        */
        //dst_array[0][0] = - 3 * src_array[0][0];
        //dst_array[0][1] = - 3 * src_array[0][1];
        //dst_array[0][2] = - 3 * src_array[0][2];

    }
    
};







template<Addressing_T A>
class ConvectiveOperator : public OneStepOperator<A>
{
public:
    ConvectiveOperator(Grid<A> &src)
        : OneStepOperator<A>::OneStepOperator(0.0f, 0.0f, 0.5f / src.dx)
    {}

protected:
    void apply_locally(
        std::array< std::array<Real, 3>, 13 > &dst_array, 
        std::array< std::array<Real, 3>, 13 > &src_array) override
    {
        // First let's interpolate the value of the velocity 
        // (using second order taylor expansion)
        Real v_expansion_i = 0.75 * src_array[0][1] + 
            0.375 * src_array[7][1] - 0.125 * src_array[8][1];
        
        Real u_expansion_j = 0.75 * src_array[0][0] + 
            0.375 * src_array[8][0] - 0.125 * src_array[7][0];

        Real w_expansion_i = 0.75 * src_array[0][2] + 
            0.375 * src_array[9][2] - 0.125 * src_array[10][2];

        Real u_expansion_k = 0.75 * src_array[0][0] + 
            0.375 * src_array[10][0] - 0.125 * src_array[9][0];

        Real v_expansion_k = 0.75 * src_array[0][1] + 
            0.375 * src_array[12][1] - 0.125 * src_array[11][1];
        
        Real w_expansion_j = 0.75 * src_array[0][2] + 
            0.375 * src_array[11][2] - 0.125 * src_array[12][2];

        // Then apply the values
        dst_array[0][0] = this->coefficients[2] * (
            src_array[0][0] * (src_array[1][0] - src_array[2][0]) +
            v_expansion_i * (src_array[3][0] - src_array[4][0]) +
            w_expansion_i * (src_array[5][0] - src_array[6][0]));
        
        dst_array[0][1] = this->coefficients[2] * (
            u_expansion_j * (src_array[1][1] - src_array[2][1]) +
            src_array[0][1] * (src_array[3][1] - src_array[4][1]) +
            w_expansion_j * (src_array[5][1] - src_array[6][1]));
        
        dst_array[0][2] = this->coefficients[2] * (
            u_expansion_k * (src_array[1][2] - src_array[2][2]) +
            v_expansion_k * (src_array[3][2] - src_array[4][2]) +
            src_array[0][2] * (src_array[5][2] - src_array[6][2]));

    }

    
};






template<Addressing_T A>
Real computeError(const Grid<A> &grid, VectorialFunction &fun)
{
    Real sum = 0.0;
    Real sum_exact = 0.0;
    
    Real sdx = grid.dx/2;
    Real sdy = grid.dy/2;
    Real sdz = grid.dz/2;

    //const auto &grid = model.grid;

    // Loop through the entire grid (no boundary)
    for (index_t i = 1; i < grid.nx - 2; ++i) {
        for (index_t j = 1; j < grid.ny - 2; ++j) {
            for (index_t k = 1; k < grid.nz - 2; ++k) {

                // Convert grid indices to real space coordinates
                Real x = static_cast<Real>(i) * grid.dx;
                Real y = static_cast<Real>(j) * grid.dy;
                Real z = static_cast<Real>(k) * grid.dz;

                // Calculate the exact solution for each component
                //Real exactU = fun.u(x + sdx, y, z);
                //Real exactV = fun.v(x, y + sdy, z);
                //Real exactW = fun.w(x, y, z + sdz);

                // Calculate the exact solution for each component
                Real exactU = fun.u(x + grid.dx, y + sdy, z + sdz);
                Real exactV = fun.v(x + sdx, y + grid.dy, z + sdz);
                Real exactW = fun.w(x + sdx, y + sdy, z + grid.dz);

                // Access the computed grid components
                Real gridU = grid(Component::U, i, j, k);
                Real gridV = grid(Component::V, i, j, k);
                Real gridW = grid(Component::W, i, j, k);

                // Calculate the differences
                Real diffU = gridU - exactU;
                Real diffV = gridV - exactV;
                Real diffW = gridW - exactW;

                // Add the squares of the differences to sum
                sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);
                sum_exact += (exactU * exactU) + (exactV * exactV) + (exactW * exactW);
                //std::cout << i << " " << j << " " << k << std::endl;
                //std::cout << "\tExpected U = " << exactU << "\tObtained U = " << gridU << std::endl;
                //std::cout << "\tExpected V = " << exactV << "\tObtained V = " << gridV << std::endl;
                //std::cout << "\tExpected W = " << exactW << "\tObtained W = " << gridW << std::endl;
            }
        }
    }
    // Commented 'cause I'm currently trying to generalize the bc's
    // Loop on the boundary faces
    // west face: i = 0
    /*
    for (index_t j = 1; j < grid.ny - 2; ++j)
    {
        for (index_t k = 0; k < grid.nz - 2; ++k)
        {
            // Convert grid indices to real space coordinates
            Real x = static_cast<Real>(0) * model.dx;
            Real y = static_cast<Real>(j) * model.dy;
            Real z = static_cast<Real>(k) * model.dz;

            // Calculate the exact solution for each component
            Real exactU = fun.u(x + sdx, y, z);

            // Access the computed grid components
            Real gridU = grid(Component::U, 0, j, k);

            // Calculate the differences
            Real diffU = gridU - exactU;

            // Add the squares of the differences to sum
            sum += (diffU * diffU);

            //std::cout << 0 << " " << j << " " << k << std::endl;
            //std::cout << "\tExpected U = " << exactU << "\tObtained U = " << gridU << std::endl;
        }
    }
    // South face : j = 0
    
    for (index_t i = 1; i < grid.nx - 2; ++i)
    {
        for (index_t k = 1; k < grid.nz - 2; ++k)
        {
            // Convert grid indices to real space coordinates
            Real x = static_cast<Real>(i) * model.dx;
            Real y = static_cast<Real>(0) * model.dy;
            Real z = static_cast<Real>(k) * model.dz;

            // Calculate the exact solution for each component
            Real exactV = fun.v(x, y + sdy, z);

            // Access the computed grid components
            Real gridV = grid(Component::V, i, 0, k);

            // Calculate the differences
            Real diffV = gridV - exactV;

            // Add the squares of the differences to sum
            sum += (diffV * diffV);

            //std::cout << i << " " << 0 << " " << k << std::endl;
            //std::cout << "\tExpected V = " << exactV << "\tObtained V = " << gridV << std::endl;
        }
    }
    
    // Down face : k = 0
    for (index_t i = 1; i < grid.nx - 2; ++i)
    {
        for (index_t j = 1; j < grid.ny - 2; ++j)
        {
            // Convert grid indices to real space coordinates
            Real x = static_cast<Real>(i) * model.dx;
            Real y = static_cast<Real>(j) * model.dy;
            Real z = static_cast<Real>(0) * model.dz;

            // Calculate the exact solution for each component
            Real exactW = fun.w(x, y, z + sdz);

            // Access the computed grid components
            Real gridW = grid(Component::W, i, j, 0);

            // Calculate the differences
            Real diffW = gridW - exactW;

            // Add the squares of the differences to sum
            sum += (diffW * diffW);

            //std::cout << i << " " << j << " " << 0 << std::endl;
            //std::cout << "\tExpected W = " << exactW << "\tObtained W = " << gridW << std::endl;
        }
    }

    // East face : i = nx - 2
    for (index_t j = 1; j < grid.ny - 2; ++j)
    {
        for (index_t k = 1; k < grid.nz - 2; ++k)
        {
            // Convert grid indices to real space coordinates
            Real x = static_cast<Real>(grid.nz - 2) * model.dx;
            Real y = static_cast<Real>(j) * model.dy;
            Real z = static_cast<Real>(k) * model.dz;

            // Calculate the exact solution for each component
            Real exactU = fun.u(x + sdx, y, z);

            // Access the computed grid components
            Real gridU = grid(Component::U, grid.nz - 2, j, k);

            // Calculate the differences
            Real diffU = gridU - exactU;

            // Add the squares of the differences to sum
            sum += (diffU * diffU);

            //std::cout << 0 << " " << j << " " << k << std::endl;
            //std::cout << "\tExpected U = " << exactU << "\tObtained U = " << gridU << std::endl;
        }
    }

    // North face : j = ny - 2
    
    for (index_t i = 1; i < grid.nx - 2; ++i)
    {
        for (index_t k = 1; k < grid.nz - 2; ++k)
        {
            // Convert grid indices to real space coordinates
            Real x = static_cast<Real>(i) * model.dx;
            Real y = static_cast<Real>(grid.ny - 2) * model.dy;
            Real z = static_cast<Real>(k) * model.dz;

            // Calculate the exact solution for each component
            Real exactV = fun.v(x, y + sdy, z);

            // Access the computed grid components
            Real gridV = grid(Component::V, i, grid.ny - 2, k);

            // Calculate the differences
            Real diffV = gridV - exactV;

            // Add the squares of the differences to sum
            sum += (diffV * diffV);

            //std::cout << i << " " << 0 << " " << k << std::endl;
            //std::cout << "\tExpected V = " << exactV << "\tObtained V = " << gridV << std::endl;
        }
    }

    // Up face : k = nz - 2
    for (index_t i = 1; i < grid.nx - 2; ++i)
    {
        for (index_t j = 1; j < grid.ny - 2; ++j)
        {
            // Convert grid indices to real space coordinates
            Real x = static_cast<Real>(i) * model.dx;
            Real y = static_cast<Real>(j) * model.dy;
            Real z = static_cast<Real>(grid.nz - 2) * model.dz;

            // Calculate the exact solution for each component
            Real exactW = fun.w(x, y, z + sdz);

            // Access the computed grid components
            Real gridW = grid(Component::W, i, j, grid.nz - 2);

            // Calculate the differences
            Real diffW = gridW - exactW;

            // Add the squares of the differences to sum
            sum += (diffW * diffW);

            //std::cout << i << " " << j << " " << 0 << std::endl;
            //std::cout << "\tExpected W = " << exactW << "\tObtained W = " << gridW << std::endl;
        }
    }
    */
    
    // Up north east
    /*
    {
        // Convert grid indices to real space coordinates
        Real x = static_cast<Real>(grid.nx - 2) * model.dx;
        Real y = static_cast<Real>(grid.ny - 2) * model.dy;
        Real z = static_cast<Real>(grid.nz - 2) * model.dz;

        // Calculate the exact solution for each component
        Real exactU = fun.u(x + sdx, y, z);
        Real exactV = fun.v(x, y + sdy, z);
        Real exactW = fun.w(x, y, z + sdz);

        // Access the computed grid components
        Real gridU = grid(Component::U, grid.nx - 2, grid.ny - 2, grid.nz - 2);
        Real gridV = grid(Component::V, grid.nx - 2, grid.ny - 2, grid.nz - 2);
        Real gridW = grid(Component::W, grid.nx - 2, grid.ny - 2, grid.nz - 2);

        // Calculate the differences
        Real diffU = gridU - exactU;
        Real diffV = gridV - exactV;
        Real diffW = gridW - exactW;

        // Add the squares of the differences to sum
        sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);
        std::cout << grid.nx - 2 << " " << grid.ny - 2 << " " << grid.nz - 2 << std::endl;
        std::cout << "\tExpected U = " << exactU << "\tObtained U = " << gridU << std::endl;
        std::cout << "\tExpected V = " << exactV << "\tObtained V = " << gridV << std::endl;
        std::cout << "\tExpected W = " << exactW << "\tObtained W = " << gridW << std::endl;

    }
    */


    return std::sqrt(sum * (grid.dx * grid.dy * grid.dz))/std::sqrt(sum_exact * (grid.dx * grid.dy * grid.dz));
}




// Just for the moment for testing
void set_boundary_values(Grid<STANDARD> &grid, const VectorFunction &dirichlet)
{
    //auto &grid = model.grid;
    // For now i try with the west face
    for (unsigned int j = 0; j < grid.ny; ++j)
    {
        for (unsigned int k = 0; k < grid.nz; ++k)
        {
            Real x = 0;
            Real y = static_cast<Real>(j) * grid.dy;
            Real z = static_cast<Real>(k) * grid.dz;

            grid(Component::U, -1, j, k) = dirichlet(x, y, z)[0];
            //grid(Component::V, -1, j, k) = dirichlet(x, y, z)[1];
            //grid(Component::W, -1, j, k) = dirichlet(x, y, z)[2];
        }
    }

    // South face 
    for (unsigned int i = 0; i < grid.nx; ++i)
    {
        for (unsigned int k = 0; k < grid.nz; ++k)
        {
            Real x = static_cast<Real>(i) * grid.dx;
            Real y = 0;
            Real z = static_cast<Real>(k) * grid.dz;

            //grid(Component::U, i, -1, k) = dirichlet(x, y, z)[0];
            grid(Component::V, i, -1, k) = dirichlet(x, y, z)[1];
            //grid(Component::W, i, -1, k) = dirichlet(x, y, z)[2];
        }
    }

    // Down face
    for (unsigned int i = 0; i < grid.nx; ++i)
    {
        for (unsigned int j = 0; j < grid.ny; ++j)
        {
            Real x = static_cast<Real>(i) * grid.dx;
            Real y = static_cast<Real>(j) * grid.dy;
            Real z = 0;

            //grid(Component::U, i, j, -1) = dirichlet(x, y, z)[0];
            //grid(Component::V, i, j, -1) = dirichlet(x, y, z)[1];
            grid(Component::W, i, j, -1) = dirichlet(x, y, z)[2];
        }
    }
    
    // East face
    for (unsigned int j = 0; j < grid.ny; ++j)
    {
        for (unsigned int k = 0; k < grid.nz; ++k)
        {
            Real x = 1;
            Real y = static_cast<Real>(j) * grid.dy;
            Real z = static_cast<Real>(k) * grid.dz;

            grid(Component::U, grid.nx - 1, j, k) = dirichlet(x, y, z)[0];
            //grid(Component::V, -1, j, k) = dirichlet(x, y, z)[1];
            //grid(Component::W, -1, j, k) = dirichlet(x, y, z)[2];
        }
    }

    // North face 
    for (unsigned int i = 0; i < grid.nx; ++i)
    {
        for (unsigned int k = 0; k < grid.nz; ++k)
        {
            Real x = static_cast<Real>(i) * grid.dx;
            Real y = 1;
            Real z = static_cast<Real>(k) * grid.dz;

            //grid(Component::U, i, -1, k) = dirichlet(x, y, z)[0];
            grid(Component::V, i, grid.ny - 1, k) = dirichlet(x, y, z)[1];
            //grid(Component::W, i, -1, k) = dirichlet(x, y, z)[2];
        }
    }

    // Up face
    for (unsigned int i = 0; i < grid.nx; ++i)
    {
        for (unsigned int j = 0; j < grid.ny; ++j)
        {
            Real x = static_cast<Real>(i) * grid.dx;
            Real y = static_cast<Real>(j) * grid.dy;
            Real z = 1;

            //grid(Component::U, i, j, -1) = dirichlet(x, y, z)[0];
            //grid(Component::V, i, j, -1) = dirichlet(x, y, z)[1];
            grid(Component::W, i, j, grid.nz - 1) = dirichlet(x, y, z)[2];
        }
    }
    
    
}




#endif