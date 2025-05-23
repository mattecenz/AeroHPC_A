#ifndef _2DECOMPCH_
#define _2DECOMPCH_

#include "Traits.hpp"

#include "mpi.h"
#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <cmath>
#include <string>
#include <memory.h>

// Allow to match library datatypes
typedef Real C2D_DTYPE;
#define C2D_MPI_DTYPE Real_MPI

class C2Decomp {

public:

    int myTypeBytes;

    //Global Size
    int nxGlobal, nyGlobal, nzGlobal;

    //MPI rank info
    int nRank, nProc;


public:
    //parameters for 2D Cartesian Topology
    int dims[2], coord[2];
    int periodic[2];


public:
    MPI_Comm DECOMP_2D_COMM_CART_X, DECOMP_2D_COMM_CART_Y, DECOMP_2D_COMM_CART_Z;
    MPI_Comm DECOMP_2D_COMM_ROW, DECOMP_2D_COMM_COL;


	//Defining neighboring blocks
	int neighbor[3][6];

	private:
	//Flags for periodic condition in 3D
	bool periodicX, periodicY, periodicZ;


public:
    //Struct used to store decomposition info for a given global data size
    typedef struct decompinfo {
        int xst[3], xen[3], xsz[3];
        int yst[3], yen[3], ysz[3];
        int zst[3], zen[3], zsz[3];

        int *x1dist, *y1dist, *y2dist, *z2dist;
        int *x1cnts, *y1cnts, *y2cnts, *z2cnts;
        int *x1disp, *y1disp, *y2disp, *z2disp;

        int x1count, y1count, y2count, z2count;

        bool even;
    } DecompInfo;


public:
    //main default decomposition information for global size nx*ny*nz
    DecompInfo decompMain;
    int decompBufSize;


public:
    //Starting/ending index and size of data held by the current processor
    //duplicate 'decompMain', needed by apps to define data structure
    int xStart[3], xEnd[3], xSize[3]; //x-pencil
    int yStart[3], yEnd[3], ySize[3]; //y-pencil
    int zStart[3], zEnd[3], zSize[3]; //z-pencil


private:
    //These are the buffers used by MPI_ALLTOALL(V) calls
    C2D_DTYPE *work1_r, *work2_r; //Only implementing real for now...


public:

    C2Decomp(const int nx, const int ny, const int nz,
            const int pRow, const int pCol, const bool periodicBC[3]) {

        nxGlobal = nx;
        nyGlobal = ny;
        nzGlobal = nz;

        periodicX = periodicBC[0];
        periodicY = periodicBC[1];
        periodicZ = periodicBC[2];

        decompBufSize = 0;
        work1_r = NULL;
        work2_r = NULL;

        decomp2DInit(pRow, pCol);
    }

    void decomp2DInit(int pRow, int pCol);

    void best2DGrid(int nProc, int &pRow, int &pCol);

    void FindFactor(int num, int *factors, int &nfact);

    void decomp2DFinalize();

    //Just get it running without the optional decomp for now...
    void transposeX2Y(C2D_DTYPE *src, C2D_DTYPE *dst);

    void transposeY2Z(C2D_DTYPE *src, C2D_DTYPE *dst);

    void transposeZ2Y(C2D_DTYPE *src, C2D_DTYPE *dst);

    void transposeY2X(C2D_DTYPE *src, C2D_DTYPE *dst);

    //Get Transposes but with array indexing with the major index of the pencil...
    void transposeX2Y_MajorIndex(C2D_DTYPE *src, C2D_DTYPE *dst);

    void transposeY2Z_MajorIndex(C2D_DTYPE *src, C2D_DTYPE *dst);

    void transposeZ2Y_MajorIndex(C2D_DTYPE *src, C2D_DTYPE *dst);

    void transposeY2X_MajorIndex(C2D_DTYPE *src, C2D_DTYPE *dst);


    //calls for overlapping communication and computation...
    void transposeX2Y_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeX2Y_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2Z_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2Z_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeZ2Y_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeZ2Y_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2X_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2X_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    //calls for overlapping communication and computation...
    void transposeX2Y_MajorIndex_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeX2Y_MajorIndex_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2Z_MajorIndex_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2Z_MajorIndex_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeZ2Y_MajorIndex_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeZ2Y_MajorIndex_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2X_MajorIndex_Start(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);

    void transposeY2X_MajorIndex_Wait(MPI_Request &handle, C2D_DTYPE *src, C2D_DTYPE *dst, C2D_DTYPE *sbuf, C2D_DTYPE *rbuf);


    void decompInfoInit();

    void decompInfoFinalize();


    //only doing real
    void allocX(C2D_DTYPE *&var) const;

    void allocY(C2D_DTYPE *&var) const;

    void allocZ(C2D_DTYPE *&var) const;

    static void deallocXYZ(C2D_DTYPE *&var);

    void updateHalo(C2D_DTYPE *in, C2D_DTYPE *&out, int level, int ipencil);

    void decomp2DAbort(int errorCode, std::string msg);

    void initNeighbor();

    void getDist();

    void distribute(int data1, int proc, int *st, int *en, int *sz);

    void partition(int nx, int ny, int nz, int *pdim, int *lstart, int *lend, int *lsize);

    void prepareBuffer(DecompInfo *dii);

    void getDecompInfo(DecompInfo dcompinfo_in);

    void memSplitXY(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeXY(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeXY_YMajor(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memSplitYZ(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memSplitYZ_YMajor(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeYZ(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeYZ_ZMajor(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memSplitZY(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memSplitZY_ZMajor(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeZY(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeZY_YMajor(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memSplitYX(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memSplitYX_YMajor(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    void memMergeYX(C2D_DTYPE *in, int n1, int n2, int n3, C2D_DTYPE *out, int iproc, int *dist);

    //IO
    void writeOne(int ipencil, C2D_DTYPE *var, std::string &filename);

    void writeVar(MPI_File &fh, MPI_Offset &disp, int ipencil, C2D_DTYPE *var);

    void writeScalar(MPI_File &fh, MPI_Offset &disp, int n, C2D_DTYPE *var);

    void writePlane(int ipencil, C2D_DTYPE *var, int iplane, int n, const std::string& filename);

    void writeEvery(int ipencil, C2D_DTYPE *var, int iskip, int jskip, int kskip, std::string &filename, bool from1);

    void readOne(int ipencil, C2D_DTYPE *var, const std::string& filename);

    void readVar(MPI_File &fh, MPI_Offset &disp, int ipencil, C2D_DTYPE *var);

    void readScalar(MPI_File &fh, MPI_Offset &disp, int n, C2D_DTYPE *var);


};


#endif
