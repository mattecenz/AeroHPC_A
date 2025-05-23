#include "C2Decomp.hpp"

void C2Decomp::allocX(C2D_DTYPE *&var) const {

    int xsize = decompMain.xsz[0];
    int ysize = decompMain.xsz[1];
    int zsize = decompMain.xsz[2];

    var = new Real[xsize * ysize * zsize];

}

void C2Decomp::allocY(C2D_DTYPE *&var) const {

    int xsize = decompMain.ysz[0];
    int ysize = decompMain.ysz[1];
    int zsize = decompMain.ysz[2];

    var = new Real[xsize * ysize * zsize];

}

void C2Decomp::allocZ(C2D_DTYPE *&var) const {

    int xsize = decompMain.zsz[0];
    int ysize = decompMain.zsz[1];
    int zsize = decompMain.zsz[2];

    var = new Real[xsize * ysize * zsize];

}

void C2Decomp::deallocXYZ(C2D_DTYPE *&var) {
    delete[] var;
    var = NULL;
}


