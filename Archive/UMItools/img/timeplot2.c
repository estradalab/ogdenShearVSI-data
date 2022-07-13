/*
 * MATLAB Compiler: 3.0
 * Date: Tue Feb  4 10:01:39 2003
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "timeplot2" 
 */
#include "timeplot2.h"
#include "errormesg.h"
#include "libmatlbm.h"
#include "pwd.h"
#include "read_hdr.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;

static mxChar _array3_[5] = { '*', '.', 'i', 'm', 'g' };
static mxArray * _mxarray2_;

static double _array5_[2] = { 0.0, 1.0 };
static mxArray * _mxarray4_;
static mxArray * _mxarray6_;

static mxChar _array8_[23] = { '%', 's', '-', '-', '-', '-', '-', 'i',
                               'm', 'a', 'g', 'e', 's', ' ', 'n', 'o',
                               't', ' ', 'f', 'o', 'u', 'n', 'd' };
static mxArray * _mxarray7_;

static mxChar _array10_[5] = { '*', '.', 'h', 'd', 'r' };
static mxArray * _mxarray9_;
static mxArray * _mxarray11_;
static mxArray * _mxarray12_;

static mxChar _array14_[4] = { 'i', 'n', 't', '8' };
static mxArray * _mxarray13_;
static mxArray * _mxarray15_;

static mxChar _array17_[5] = { 'u', 'i', 'n', 't', '8' };
static mxArray * _mxarray16_;
static mxArray * _mxarray18_;

static mxChar _array20_[5] = { 's', 'h', 'o', 'r', 't' };
static mxArray * _mxarray19_;

static mxChar _array22_[3] = { 'i', 'n', 't' };
static mxArray * _mxarray21_;
static mxArray * _mxarray23_;

static mxChar _array25_[5] = { 'f', 'l', 'o', 'a', 't' };
static mxArray * _mxarray24_;
static mxArray * _mxarray26_;

static mxChar _array28_[34] = { 'D', 'a', 't', 'a', ' ', 'T', 'y', 'p', 'e',
                                ' ', '%', 'd', ' ', 'U', 'n', 's', 'u', 'p',
                                'p', 'o', 'r', 't', 'e', 'd', '.', ' ', 'A',
                                'b', 'o', 'r', 't', 'i', 'n', 'g' };
static mxArray * _mxarray27_;

static mxChar _array30_[3] = { 'b', 'o', 'f' };
static mxArray * _mxarray29_;

void InitializeModule_timeplot2(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeDouble(8.0);
    _mxarray2_ = mclInitializeString(5, _array3_);
    _mxarray4_ = mclInitializeDoubleVector(1, 2, _array5_);
    _mxarray6_ = mclInitializeDouble(0.0);
    _mxarray7_ = mclInitializeString(23, _array8_);
    _mxarray9_ = mclInitializeString(5, _array10_);
    _mxarray11_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray12_ = mclInitializeDouble(-1.0);
    _mxarray13_ = mclInitializeString(4, _array14_);
    _mxarray15_ = mclInitializeDouble(2.0);
    _mxarray16_ = mclInitializeString(5, _array17_);
    _mxarray18_ = mclInitializeDouble(4.0);
    _mxarray19_ = mclInitializeString(5, _array20_);
    _mxarray21_ = mclInitializeString(3, _array22_);
    _mxarray23_ = mclInitializeDouble(16.0);
    _mxarray24_ = mclInitializeString(5, _array25_);
    _mxarray26_ = mclInitializeDouble(32.0);
    _mxarray27_ = mclInitializeString(34, _array28_);
    _mxarray29_ = mclInitializeString(3, _array30_);
}

void TerminateModule_timeplot2(void) {
    mxDestroyArray(_mxarray29_);
    mxDestroyArray(_mxarray27_);
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray23_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mtimeplot2(int nargout_,
                            mxArray * path,
                            mxArray * file,
                            mxArray * x,
                            mxArray * y,
                            mxArray * z);

_mexLocalFunctionTable _local_function_table_timeplot2
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfTimeplot2" contains the normal interface for the
 * "timeplot2" M-function from file
 * "/net/quentin/home/hernan/matlab/img/timeplot2.m" (lines 1-122). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
mxArray * mlfTimeplot2(mxArray * path,
                       mxArray * file,
                       mxArray * x,
                       mxArray * y,
                       mxArray * z) {
    int nargout = 1;
    mxArray * tdata = NULL;
    mlfEnterNewContext(0, 5, path, file, x, y, z);
    tdata = Mtimeplot2(nargout, path, file, x, y, z);
    mlfRestorePreviousContext(0, 5, path, file, x, y, z);
    return mlfReturnValue(tdata);
}

/*
 * The function "mlxTimeplot2" contains the feval interface for the "timeplot2"
 * M-function from file "/net/quentin/home/hernan/matlab/img/timeplot2.m"
 * (lines 1-122). The feval function calls the implementation version of
 * timeplot2 through this function. This function processes any input arguments
 * and passes them to the implementation version of the function, appearing
 * above.
 */
void mlxTimeplot2(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[5];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: timeplot2 Line: 1 Column:"
            " 1 The function \"timeplot2\" was called with m"
            "ore than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: timeplot2 Line: 1 Column:"
            " 1 The function \"timeplot2\" was called with m"
            "ore than the declared number of inputs (5)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 5 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 5; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    mplhs[0]
      = Mtimeplot2(nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    mlfRestorePreviousContext(
      0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mtimeplot2" is the implementation version of the "timeplot2"
 * M-function from file "/net/quentin/home/hernan/matlab/img/timeplot2.m"
 * (lines 1-122). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function tdata = timeplot2(path, file,x,y,z)
 */
static mxArray * Mtimeplot2(int nargout_,
                            mxArray * path,
                            mxArray * file,
                            mxArray * x,
                            mxArray * y,
                            mxArray * z) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_timeplot2);
    mxArray * tdata = NULL;
    mxArray * data = NULL;
    mxArray * n = NULL;
    mxArray * ydim = NULL;
    mxArray * xdim = NULL;
    mxArray * bytes = NULL;
    mxArray * fmt = NULL;
    mxArray * mesg = NULL;
    mxArray * fp = NULL;
    mxArray * position = NULL;
    mxArray * k = NULL;
    mxArray * j = NULL;
    mxArray * i = NULL;
    mxArray * endz = NULL;
    mxArray * endy = NULL;
    mxArray * endx = NULL;
    mxArray * startz = NULL;
    mxArray * starty = NULL;
    mxArray * startx = NULL;
    mxArray * num = NULL;
    mxArray * hdr = NULL;
    mxArray * hfiles = NULL;
    mxArray * files = NULL;
    mxArray * root = NULL;
    mxArray * sz = NULL;
    mxArray * ans = NULL;
    mxArray * oldpath = NULL;
    mclCopyArray(&path);
    mclCopyArray(&file);
    mclCopyArray(&x);
    mclCopyArray(&y);
    mclCopyArray(&z);
    /*
     * 
     * % function result = timeplot2(path,file,x,y,z)
     * %
     * % Returns the intensity of the voxels in the cube determined
     * % by the two points (x(1), y(1)), z(1)) and (y(1), y(2)) z(2)) in a time
     * % series of fmri???.img.  
     * % the data comes back as a row vector
     * %
     * oldpath=pwd;
     */
    mlfAssign(&oldpath, mlfPwd());
    /*
     * cd (path);
     */
    mclAssignAns(&ans, mlfNCd(0, mclVa(path, "path")));
    /*
     * 
     * sz = size(file);
     */
    mlfAssign(&sz, mlfSize(mclValueVarargout(), mclVa(file, "file"), NULL));
    /*
     * root = file(1,1:sz(2)-8);
     */
    mlfAssign(
      &root,
      mclArrayRef2(
        mclVa(file, "file"),
        _mxarray0_,
        mlfColon(
          _mxarray0_,
          mclMinus(mclIntArrayRef1(mclVv(sz, "sz"), 2), _mxarray1_),
          NULL)));
    /*
     * 
     * files = dir(strcat(root,'*.img'));
     */
    mlfAssign(
      &files, mlfNDir(1, mlfStrcat(mclVv(root, "root"), _mxarray2_, NULL)));
    /*
     * if (size(files)==[0 1])
     */
    if (mclEqBool(
          mlfSize(mclValueVarargout(), mclVv(files, "files"), NULL),
          _mxarray4_)) {
        /*
         * tdata=0;
         */
        mlfAssign(&tdata, _mxarray6_);
        /*
         * fprintf('%s-----images not found',file);
         */
        mclAssignAns(
          &ans, mlfNFprintf(0, _mxarray7_, mclVa(file, "file"), NULL));
        /*
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * 
     * hfiles = dir(strcat(root,'*.hdr'));
     */
    mlfAssign(
      &hfiles, mlfNDir(1, mlfStrcat(mclVv(root, "root"), _mxarray9_, NULL)));
    /*
     * sz = size(files);
     */
    mlfAssign(&sz, mlfSize(mclValueVarargout(), mclVv(files, "files"), NULL));
    /*
     * hfiles(1).name;
     */
    mclAssignAns(
      &ans, mlfIndexRef(mclVv(hfiles, "hfiles"), "(?).name", _mxarray0_));
    /*
     * hdr = read_hdr(hfiles(1).name);
     */
    mlfAssign(
      &hdr,
      mclFeval(
        mclValueVarargout(),
        mlxRead_hdr,
        mlfIndexRef(mclVv(hfiles, "hfiles"), "(?).name", _mxarray0_),
        NULL));
    /*
     * 
     * 
     * 
     * % determine which positions are included into the ROI
     * num = 1;
     */
    mlfAssign(&num, _mxarray0_);
    /*
     * 
     * startx = min(x);
     */
    mlfAssign(&startx, mlfMin(NULL, mclVa(x, "x"), NULL, NULL));
    /*
     * starty = min(y);
     */
    mlfAssign(&starty, mlfMin(NULL, mclVa(y, "y"), NULL, NULL));
    /*
     * startz = min(z);
     */
    mlfAssign(&startz, mlfMin(NULL, mclVa(z, "z"), NULL, NULL));
    /*
     * endx = max(x);
     */
    mlfAssign(&endx, mlfMax(NULL, mclVa(x, "x"), NULL, NULL));
    /*
     * endy = max(y);
     */
    mlfAssign(&endy, mlfMax(NULL, mclVa(y, "y"), NULL, NULL));
    /*
     * endz = max(z);
     */
    mlfAssign(&endz, mlfMax(NULL, mclVa(z, "z"), NULL, NULL));
    /*
     * 
     * 
     * for i=startx:endx
     */
    {
        mclForLoopIterator viter__;
        for (mclForStart(
               &viter__, mclVv(startx, "startx"), mclVv(endx, "endx"), NULL);
             mclForNext(&viter__, &i);
             ) {
            mclForLoopIterator viter__0;
            /*
             * for j=starty:endy
             */
            for (mclForStart(
                   &viter__0,
                   mclVv(starty, "starty"),
                   mclVv(endy, "endy"),
                   NULL);
                 mclForNext(&viter__0, &j);
                 ) {
                mclForLoopIterator viter__1;
                /*
                 * for k=startz:endz
                 */
                for (mclForStart(
                       &viter__1,
                       mclVv(startz, "startz"),
                       mclVv(endz, "endz"),
                       NULL);
                     mclForNext(&viter__1, &k);
                     ) {
                    /*
                     * position(num) = (k-1)*(hdr.xdim*hdr.ydim) + (j-1)*hdr.xdim + i-1;
                     */
                    mclArrayAssign1(
                      &position,
                      mclMinus(
                        mclPlus(
                          mclPlus(
                            mclMtimes(
                              mclMinus(mclVv(k, "k"), _mxarray0_),
                              mclFeval(
                                mclValueVarargout(),
                                mlxMtimes,
                                mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
                                mlfIndexRef(mclVv(hdr, "hdr"), ".ydim"),
                                NULL)),
                            mclFeval(
                              mclValueVarargout(),
                              mlxMtimes,
                              mclMinus(mclVv(j, "j"), _mxarray0_),
                              mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
                              NULL)),
                          mclVv(i, "i")),
                        _mxarray0_),
                      mclVv(num, "num"));
                    /*
                     * num = num +1;
                     */
                    mlfAssign(&num, mclPlus(mclVv(num, "num"), _mxarray0_));
                /*
                 * end
                 */
                }
                mclDestroyForLoopIterator(viter__1);
            /*
             * end
             */
            }
            mclDestroyForLoopIterator(viter__0);
        /*
         * end
         */
        }
        mclDestroyForLoopIterator(viter__);
    }
    /*
     * 
     * num = num -1;
     */
    mlfAssign(&num, mclMinus(mclVv(num, "num"), _mxarray0_));
    /*
     * 
     * % extract data from files
     * tdata = zeros(1,sz(1));
     */
    mlfAssign(
      &tdata, mlfZeros(_mxarray0_, mclIntArrayRef1(mclVv(sz, "sz"), 1), NULL));
    /*
     * 
     * for i=1:sz(1)
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclIntArrayRef1(mclVv(sz, "sz"), 1));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray11_);
        } else {
            /*
             * 
             * [fp mesg]= fopen(files(i).name);
             * %	disp(i)
             * %disp  ( files(i).name)
             * 
             * if fp == -1
             * disp(mesg);
             * return
             * end
             * 
             * 
             * switch hdr.datatype     
             * case 0
             * fmt = 'int8';
             * bytes = 1;
             * 
             * case 2
             * fmt = 'uint8';
             * bytes = 1;
             * case 4
             * fmt = 'short';
             * bytes = 2;
             * case 8
             * fmt = 'int';
             * bytes = 2;
             * case 16
             * fmt = 'float';
             * bytes = 4;
             * case 32
             * fmt = 'float';
             * xdim = hdr.xdim * 2;
             * ydim = hdr.ydim * 2;
             * bytes = 8;
             * 
             * otherwise
             * errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.bits));
             * return
             * 
             * end
             * 
             * 
             * % average the ROI positions
             * for n=1:num      
             * fseek(fp,(position(n))*bytes,'bof');     
             * data = fread(fp,1,fmt);
             * tdata(i) = tdata(i) + data;
             * end
             * 
             * tdata(i) = tdata(i) / num;
             * 
             * fclose(fp);
             * end
             */
            for (; ; ) {
                mclFeval(
                  mlfVarargout(&fp, &mesg, NULL),
                  mlxFopen,
                  mlfIndexRef(mclVv(files, "files"), "(?).name", mlfScalar(v_)),
                  NULL);
                if (mclEqBool(mclVv(fp, "fp"), _mxarray12_)) {
                    mlfDisp(mclVv(mesg, "mesg"));
                    goto return_;
                }
                {
                    mxArray * v_0
                      = mclInitialize(
                          mlfIndexRef(mclVv(hdr, "hdr"), ".datatype"));
                    if (mclSwitchCompare(v_0, _mxarray6_)) {
                        mlfAssign(&fmt, _mxarray13_);
                        mlfAssign(&bytes, _mxarray0_);
                    } else if (mclSwitchCompare(v_0, _mxarray15_)) {
                        mlfAssign(&fmt, _mxarray16_);
                        mlfAssign(&bytes, _mxarray0_);
                    } else if (mclSwitchCompare(v_0, _mxarray18_)) {
                        mlfAssign(&fmt, _mxarray19_);
                        mlfAssign(&bytes, _mxarray15_);
                    } else if (mclSwitchCompare(v_0, _mxarray1_)) {
                        mlfAssign(&fmt, _mxarray21_);
                        mlfAssign(&bytes, _mxarray15_);
                    } else if (mclSwitchCompare(v_0, _mxarray23_)) {
                        mlfAssign(&fmt, _mxarray24_);
                        mlfAssign(&bytes, _mxarray18_);
                    } else if (mclSwitchCompare(v_0, _mxarray26_)) {
                        mlfAssign(&fmt, _mxarray24_);
                        mlfAssign(
                          &xdim,
                          mclFeval(
                            mclValueVarargout(),
                            mlxMtimes,
                            mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
                            _mxarray15_,
                            NULL));
                        mlfAssign(
                          &ydim,
                          mclFeval(
                            mclValueVarargout(),
                            mlxMtimes,
                            mlfIndexRef(mclVv(hdr, "hdr"), ".ydim"),
                            _mxarray15_,
                            NULL));
                        mlfAssign(&bytes, _mxarray1_);
                    } else {
                        mlfErrormesg(
                          mlfSprintf(
                            NULL,
                            _mxarray27_,
                            mlfIndexRef(mclVv(hdr, "hdr"), ".bits"),
                            NULL));
                        mxDestroyArray(v_0);
                        goto return_;
                    }
                    mxDestroyArray(v_0);
                }
                {
                    int v_1 = mclForIntStart(1);
                    int e_0 = mclForIntEnd(mclVv(num, "num"));
                    if (v_1 > e_0) {
                        mlfAssign(&n, _mxarray11_);
                    } else {
                        for (; ; ) {
                            mclAssignAns(
                              &ans,
                              mlfFseek(
                                mclVv(fp, "fp"),
                                mclMtimes(
                                  mclIntArrayRef1(
                                    mclVv(position, "position"), v_1),
                                  mclVv(bytes, "bytes")),
                                _mxarray29_));
                            mlfAssign(
                              &data,
                              mlfFread(
                                NULL,
                                mclVv(fp, "fp"),
                                _mxarray0_,
                                mclVv(fmt, "fmt"),
                                NULL));
                            mclIntArrayAssign1(
                              &tdata,
                              mclPlus(
                                mclIntArrayRef1(mclVv(tdata, "tdata"), v_),
                                mclVv(data, "data")),
                              v_);
                            if (v_1 == e_0) {
                                break;
                            }
                            ++v_1;
                        }
                        mlfAssign(&n, mlfScalar(v_1));
                    }
                }
                mclIntArrayAssign1(
                  &tdata,
                  mclMrdivide(
                    mclIntArrayRef1(mclVv(tdata, "tdata"), v_),
                    mclVv(num, "num")),
                  v_);
                mclAssignAns(&ans, mlfFclose(mclVv(fp, "fp")));
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * 
     * cd(oldpath);
     */
    mclAssignAns(&ans, mlfNCd(0, mclVv(oldpath, "oldpath")));
    /*
     * %close
     * 
     * 
     * %plot(tdata/max(tdata))  
     * %save timeseries tdata  
     * 
     * return
     * 
     * 
     * 
     * 
     * 
     */
    return_:
    mclValidateOutput(tdata, 1, nargout_, "tdata", "timeplot2");
    mxDestroyArray(oldpath);
    mxDestroyArray(ans);
    mxDestroyArray(sz);
    mxDestroyArray(root);
    mxDestroyArray(files);
    mxDestroyArray(hfiles);
    mxDestroyArray(hdr);
    mxDestroyArray(num);
    mxDestroyArray(startx);
    mxDestroyArray(starty);
    mxDestroyArray(startz);
    mxDestroyArray(endx);
    mxDestroyArray(endy);
    mxDestroyArray(endz);
    mxDestroyArray(i);
    mxDestroyArray(j);
    mxDestroyArray(k);
    mxDestroyArray(position);
    mxDestroyArray(fp);
    mxDestroyArray(mesg);
    mxDestroyArray(fmt);
    mxDestroyArray(bytes);
    mxDestroyArray(xdim);
    mxDestroyArray(ydim);
    mxDestroyArray(n);
    mxDestroyArray(data);
    mxDestroyArray(z);
    mxDestroyArray(y);
    mxDestroyArray(x);
    mxDestroyArray(file);
    mxDestroyArray(path);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return tdata;
}
