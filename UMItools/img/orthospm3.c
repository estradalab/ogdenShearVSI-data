/*
 * MATLAB Compiler: 3.0
 * Date: Tue Feb  4 10:03:45 2003
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "orthospm3" 
 */
#include "orthospm3.h"
#include "mwservices.h"
#include "axis.h"
#include "close.h"
#include "colormap.h"
#include "errorbar.h"
#include "event_avg.h"
#include "gca.h"
#include "ginput.h"
#include "hold.h"
#include "interp3.h"
#include "libmatlbm.h"
#include "mean.h"
#include "meshgrid.h"
#include "ov.h"
#include "pwd.h"
#include "read_hdr.h"
#include "read_img2.h"
#include "subplot.h"
#include "timeplot2.h"
#include "title.h"

extern mxArray * SPM_scale_factor;

static mxChar _array1_[3] = { 'a', 'l', 'l' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;

static mxChar _array4_[5] = { '*', '.', 'i', 'm', 'g' };
static mxArray * _mxarray3_;

static mxChar _array6_[21] = { 'S', 'e', 'l', 'e', 'c', 't', ' ',
                               'S', 'P', 'M', ' ', '*', '.', 'i',
                               'm', 'g', ' ', 'f', 'i', 'l', 'e' };
static mxArray * _mxarray5_;
static mxArray * _mxarray7_;
static mxArray * _mxarray8_;

static mxChar _array10_[4] = { '.', 'i', 'm', 'g' };
static mxArray * _mxarray9_;

static mxChar _array12_[4] = { '.', 'h', 'd', 'r' };
static mxArray * _mxarray11_;
static double _ieee_nan_;
static mxArray * _mxarray13_;

static mxChar _array15_[28] = { 'S', 'e', 'l', 'e', 'c', 't', ' ',
                                'a', 'n', 'a', 't', 'o', 'm', 'i',
                                'c', 'a', 'l', ' ', '*', '.', 'i',
                                'm', 'g', ' ', 'f', 'i', 'l', 'e' };
static mxArray * _mxarray14_;

static mxChar _array17_[7] = { 'n', 'e', 'a', 'r', 'e', 's', 't' };
static mxArray * _mxarray16_;
static mxArray * _mxarray18_;
static mxArray * _mxarray19_;

static double _array21_[384] = { 0.0, .007874015748031496, .015748031496062992,
                                 .023622047244094488, .031496062992125984,
                                 .03937007874015748, .047244094488188976,
                                 .05511811023622047, .06299212598425197,
                                 .07086614173228346, .07874015748031496,
                                 .08661417322834646, .09448818897637795,
                                 .10236220472440945, .11023622047244094,
                                 .11811023622047244, .12598425196850394,
                                 .13385826771653542, .14173228346456693,
                                 .14960629921259844, .15748031496062992,
                                 .1653543307086614, .1732283464566929,
                                 .18110236220472442, .1889763779527559,
                                 .19685039370078738, .2047244094488189,
                                 .2125984251968504, .2204724409448819,
                                 .22834645669291337, .23622047244094488,
                                 .2440944881889764, .25196850393700787,
                                 .25984251968503935, .26771653543307083,
                                 .2755905511811024, .28346456692913385,
                                 .29133858267716534, .2992125984251969,
                                 .30708661417322836, .31496062992125984,
                                 .3228346456692913, .3307086614173228,
                                 .33858267716535434, .3464566929133858,
                                 .3543307086614173, .36220472440944884,
                                 .3700787401574803, .3779527559055118,
                                 .3858267716535433, .39370078740157477,
                                 .4015748031496063, .4094488188976378,
                                 .41732283464566927, .4251968503937008,
                                 .4330708661417323, .4409448818897638,
                                 .44881889763779526, .45669291338582674,
                                 .4645669291338583, .47244094488188976,
                                 .48031496062992124, .4881889763779528,
                                 .49606299212598426, .5039370078740157,
                                 .5118110236220472, .5196850393700787,
                                 .5275590551181102, .5354330708661417,
                                 .5433070866141733, .5511811023622047,
                                 .5590551181102362, .5669291338582677,
                                 .5748031496062992, .5826771653543308,
                                 .5905511811023623, .5984251968503937,
                                 .6062992125984252, .6141732283464567,
                                 .6220472440944882, .6299212598425197,
                                 .6377952755905512, .6456692913385826,
                                 .6535433070866141, .6614173228346456,
                                 .6692913385826772, .6771653543307087,
                                 .6850393700787402, .6929133858267716,
                                 .7007874015748031, .7086614173228347,
                                 .7165354330708662, .7244094488188977,
                                 .7322834645669292, .7401574803149606,
                                 .7480314960629921, .7559055118110236,
                                 .7637795275590551, .7716535433070866,
                                 .7795275590551181, .7874015748031495,
                                 .7952755905511811, .8031496062992126,
                                 .8110236220472441, .8188976377952756,
                                 .8267716535433071, .8346456692913387,
                                 .8425196850393701, .8503937007874016,
                                 .8582677165354331, .8661417322834646,
                                 .8740157480314961, .8818897637795275,
                                 .889763779527559, .8976377952755905,
                                 .905511811023622, .9133858267716536,
                                 .9212598425196851, .9291338582677166,
                                 .937007874015748, .9448818897637795,
                                 .952755905511811, .9606299212598425,
                                 .9685039370078741, .9763779527559056,
                                 .984251968503937, .9921259842519685, 1.0, 0.0,
                                 .007874015748031496, .015748031496062992,
                                 .023622047244094488, .031496062992125984,
                                 .03937007874015748, .047244094488188976,
                                 .05511811023622047, .06299212598425197,
                                 .07086614173228346, .07874015748031496,
                                 .08661417322834646, .09448818897637795,
                                 .10236220472440945, .11023622047244094,
                                 .11811023622047244, .12598425196850394,
                                 .13385826771653542, .14173228346456693,
                                 .14960629921259844, .15748031496062992,
                                 .1653543307086614, .1732283464566929,
                                 .18110236220472442, .1889763779527559,
                                 .19685039370078738, .2047244094488189,
                                 .2125984251968504, .2204724409448819,
                                 .22834645669291337, .23622047244094488,
                                 .2440944881889764, .25196850393700787,
                                 .25984251968503935, .26771653543307083,
                                 .2755905511811024, .28346456692913385,
                                 .29133858267716534, .2992125984251969,
                                 .30708661417322836, .31496062992125984,
                                 .3228346456692913, .3307086614173228,
                                 .33858267716535434, .3464566929133858,
                                 .3543307086614173, .36220472440944884,
                                 .3700787401574803, .3779527559055118,
                                 .3858267716535433, .39370078740157477,
                                 .4015748031496063, .4094488188976378,
                                 .41732283464566927, .4251968503937008,
                                 .4330708661417323, .4409448818897638,
                                 .44881889763779526, .45669291338582674,
                                 .4645669291338583, .47244094488188976,
                                 .48031496062992124, .4881889763779528,
                                 .49606299212598426, .5039370078740157,
                                 .5118110236220472, .5196850393700787,
                                 .5275590551181102, .5354330708661417,
                                 .5433070866141733, .5511811023622047,
                                 .5590551181102362, .5669291338582677,
                                 .5748031496062992, .5826771653543308,
                                 .5905511811023623, .5984251968503937,
                                 .6062992125984252, .6141732283464567,
                                 .6220472440944882, .6299212598425197,
                                 .6377952755905512, .6456692913385826,
                                 .6535433070866141, .6614173228346456,
                                 .6692913385826772, .6771653543307087,
                                 .6850393700787402, .6929133858267716,
                                 .7007874015748031, .7086614173228347,
                                 .7165354330708662, .7244094488188977,
                                 .7322834645669292, .7401574803149606,
                                 .7480314960629921, .7559055118110236,
                                 .7637795275590551, .7716535433070866,
                                 .7795275590551181, .7874015748031495,
                                 .7952755905511811, .8031496062992126,
                                 .8110236220472441, .8188976377952756,
                                 .8267716535433071, .8346456692913387,
                                 .8425196850393701, .8503937007874016,
                                 .8582677165354331, .8661417322834646,
                                 .8740157480314961, .8818897637795275,
                                 .889763779527559, .8976377952755905,
                                 .905511811023622, .9133858267716536,
                                 .9212598425196851, .9291338582677166,
                                 .937007874015748, .9448818897637795,
                                 .952755905511811, .9606299212598425,
                                 .9685039370078741, .9763779527559056,
                                 .984251968503937, .9921259842519685, 1.0, 0.0,
                                 .007874015748031496, .015748031496062992,
                                 .023622047244094488, .031496062992125984,
                                 .03937007874015748, .047244094488188976,
                                 .05511811023622047, .06299212598425197,
                                 .07086614173228346, .07874015748031496,
                                 .08661417322834646, .09448818897637795,
                                 .10236220472440945, .11023622047244094,
                                 .11811023622047244, .12598425196850394,
                                 .13385826771653542, .14173228346456693,
                                 .14960629921259844, .15748031496062992,
                                 .1653543307086614, .1732283464566929,
                                 .18110236220472442, .1889763779527559,
                                 .19685039370078738, .2047244094488189,
                                 .2125984251968504, .2204724409448819,
                                 .22834645669291337, .23622047244094488,
                                 .2440944881889764, .25196850393700787,
                                 .25984251968503935, .26771653543307083,
                                 .2755905511811024, .28346456692913385,
                                 .29133858267716534, .2992125984251969,
                                 .30708661417322836, .31496062992125984,
                                 .3228346456692913, .3307086614173228,
                                 .33858267716535434, .3464566929133858,
                                 .3543307086614173, .36220472440944884,
                                 .3700787401574803, .3779527559055118,
                                 .3858267716535433, .39370078740157477,
                                 .4015748031496063, .4094488188976378,
                                 .41732283464566927, .4251968503937008,
                                 .4330708661417323, .4409448818897638,
                                 .44881889763779526, .45669291338582674,
                                 .4645669291338583, .47244094488188976,
                                 .48031496062992124, .4881889763779528,
                                 .49606299212598426, .5039370078740157,
                                 .5118110236220472, .5196850393700787,
                                 .5275590551181102, .5354330708661417,
                                 .5433070866141733, .5511811023622047,
                                 .5590551181102362, .5669291338582677,
                                 .5748031496062992, .5826771653543308,
                                 .5905511811023623, .5984251968503937,
                                 .6062992125984252, .6141732283464567,
                                 .6220472440944882, .6299212598425197,
                                 .6377952755905512, .6456692913385826,
                                 .6535433070866141, .6614173228346456,
                                 .6692913385826772, .6771653543307087,
                                 .6850393700787402, .6929133858267716,
                                 .7007874015748031, .7086614173228347,
                                 .7165354330708662, .7244094488188977,
                                 .7322834645669292, .7401574803149606,
                                 .7480314960629921, .7559055118110236,
                                 .7637795275590551, .7716535433070866,
                                 .7795275590551181, .7874015748031495,
                                 .7952755905511811, .8031496062992126,
                                 .8110236220472441, .8188976377952756,
                                 .8267716535433071, .8346456692913387,
                                 .8425196850393701, .8503937007874016,
                                 .8582677165354331, .8661417322834646,
                                 .8740157480314961, .8818897637795275,
                                 .889763779527559, .8976377952755905,
                                 .905511811023622, .9133858267716536,
                                 .9212598425196851, .9291338582677166,
                                 .937007874015748, .9448818897637795,
                                 .952755905511811, .9606299212598425,
                                 .9685039370078741, .9763779527559056,
                                 .984251968503937, .9921259842519685, 1.0 };
static mxArray * _mxarray20_;

static double _array23_[126] = { 0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0,
                                 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0,
                                 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0,
                                 72.0, 75.0, 78.0, 81.0, 84.0, 87.0, 90.0, 93.0,
                                 96.0, 99.0, 102.0, 105.0, 108.0, 111.0, 114.0,
                                 117.0, 120.0, 123.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static mxArray * _mxarray22_;

static double _array25_[126] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 6.0,
                                 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0,
                                 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0,
                                 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0, 78.0,
                                 81.0, 84.0, 87.0, 90.0, 93.0, 96.0, 99.0,
                                 102.0, 105.0, 108.0, 111.0, 114.0, 117.0,
                                 120.0, 123.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static mxArray * _mxarray24_;

static double _array27_[126] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 3.0, 6.0, 9.0, 12.0, 15.0,
                                 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0,
                                 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0,
                                 66.0, 69.0, 72.0, 75.0, 78.0, 81.0, 84.0, 87.0,
                                 90.0, 93.0, 96.0, 99.0, 102.0, 105.0, 108.0,
                                 111.0, 114.0, 117.0, 120.0, 123.0 };
static mxArray * _mxarray26_;
static mxArray * _mxarray28_;

static mxChar _array30_[3] = { 'o', 'f', 'f' };
static mxArray * _mxarray29_;

static mxChar _array32_[39] = { 0x005c, 'n', 'f', 'i', 'r', 's', 't',
                                ' ', 'd', 'i', 's', 'p', 'l', 'a', 'y',
                                ' ', 't', 'h', 'e', ' ', 'o', 'r', 't',
                                'h', 'o', 'g', 'o', 'n', 'a', 'l', ' ',
                                's', 'e', 'c', 't', 'i', 'o', 'n', 's' };
static mxArray * _mxarray31_;

static mxChar _array34_[10] = { 'd', 'u', 'm', 'm', 'm',
                                'y', 'n', 'a', 'm', 'e' };
static mxArray * _mxarray33_;
static mxArray * _mxarray35_;

static mxChar _array37_[55] = { 'S', 'e', 'l', 'e', 'c', 't', ' ', 'a',
                                ' ', 'f', 'i', 'l', 'e', ' ', 'f', 'r',
                                'o', 'm', ' ', 't', 'h', 'e', ' ', 'n',
                                'e', 'x', 't', ' ', 'r', 'u', 'n', '.',
                                ' ', ' ', 'c', 'a', 'n', 'c', 'e', 'l',
                                ' ', 'w', 'h', 'e', 'n', ' ', 'f', 'i',
                                'n', 'i', 's', 'h', 'e', 'd', ' ' };
static mxArray * _mxarray36_;
static mxArray * _mxarray38_;

static mxChar _array40_[66] = { '(', 'x', ',', 'y', ',', 'z', ')', '=', ' ',
                                ' ', '(', '%', '3', '.', '2', 'f', ' ', '%',
                                '3', '.', '2', 'f', ' ', '%', '3', '.', '2',
                                'f', ')', ' ', 'm', 'm', ',', ' ', ' ',
                                0x005c, 'n', ' ', '(', '%', '3', 'd', ' ',
                                '%', '3', 'd', ' ', '%', '3', 'd', ')', 'v',
                                'x', ' ', ',', ' ', 0x005c, 'n', ' ', 'v',
                                'a', 'l', '=', ' ', '%', 'd' };
static mxArray * _mxarray39_;
static mxArray * _mxarray41_;

static mxChar _array43_[6] = { 'e', 'v', 'e', 'n', 't', 's' };
static mxArray * _mxarray42_;

static mxChar _array45_[6] = { 'e', 'v', '_', 'a', 'v', 'g' };
static mxArray * _mxarray44_;

static mxChar _array47_[5] = { 't', 'd', 'a', 't', 'a' };
static mxArray * _mxarray46_;

static mxChar _array49_[3] = { 'e', 'v', '*' };
static mxArray * _mxarray48_;

static mxChar _array51_[3] = { '2', '2', '4' };
static mxArray * _mxarray50_;

static mxChar _array53_[5] = { 't', 'i', 'g', 'h', 't' };
static mxArray * _mxarray52_;

static mxChar _array55_[5] = { 'X', 't', 'i', 'c', 'k' };
static mxArray * _mxarray54_;

static mxChar _array57_[15] = { 'a', 'v', 'g', 'R', 'e', 's', 'p', 'o',
                                'n', 's', 'e', '.', 'd', 'a', 't' };
static mxArray * _mxarray56_;

static mxChar _array59_[9] = { 't', 'd', 'a', 't', 'a', '.', 'd', 'a', 't' };
static mxArray * _mxarray58_;

void InitializeModule_orthospm3(void) {
    _mxarray0_ = mclInitializeString(3, _array1_);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeString(5, _array4_);
    _mxarray5_ = mclInitializeString(21, _array6_);
    _mxarray7_ = mclInitializeDouble(1.0);
    _mxarray8_ = mclInitializeDouble(4.0);
    _mxarray9_ = mclInitializeString(4, _array10_);
    _mxarray11_ = mclInitializeString(4, _array12_);
    _ieee_nan_ = mclGetNaN();
    _mxarray13_ = mclInitializeDouble(_ieee_nan_);
    _mxarray14_ = mclInitializeString(28, _array15_);
    _mxarray16_ = mclInitializeString(7, _array17_);
    _mxarray18_ = mclInitializeDouble(3.0);
    _mxarray19_ = mclInitializeDouble(128.0);
    _mxarray20_ = mclInitializeDoubleVector(128, 3, _array21_);
    _mxarray22_ = mclInitializeDoubleVector(42, 3, _array23_);
    _mxarray24_ = mclInitializeDoubleVector(42, 3, _array25_);
    _mxarray26_ = mclInitializeDoubleVector(42, 3, _array27_);
    _mxarray28_ = mclInitializeDouble(2.0);
    _mxarray29_ = mclInitializeString(3, _array30_);
    _mxarray31_ = mclInitializeString(39, _array32_);
    _mxarray33_ = mclInitializeString(10, _array34_);
    _mxarray35_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray36_ = mclInitializeString(55, _array37_);
    _mxarray38_ = mclInitializeDouble(-10.0);
    _mxarray39_ = mclInitializeString(66, _array40_);
    _mxarray41_ = mclInitializeDouble(221.0);
    _mxarray42_ = mclInitializeString(6, _array43_);
    _mxarray44_ = mclInitializeString(6, _array45_);
    _mxarray46_ = mclInitializeString(5, _array47_);
    _mxarray48_ = mclInitializeString(3, _array49_);
    _mxarray50_ = mclInitializeString(3, _array51_);
    _mxarray52_ = mclInitializeString(5, _array53_);
    _mxarray54_ = mclInitializeString(5, _array55_);
    _mxarray56_ = mclInitializeString(15, _array57_);
    _mxarray58_ = mclInitializeString(9, _array59_);
}

void TerminateModule_orthospm3(void) {
    mxDestroyArray(_mxarray58_);
    mxDestroyArray(_mxarray56_);
    mxDestroyArray(_mxarray54_);
    mxDestroyArray(_mxarray52_);
    mxDestroyArray(_mxarray50_);
    mxDestroyArray(_mxarray48_);
    mxDestroyArray(_mxarray46_);
    mxDestroyArray(_mxarray44_);
    mxDestroyArray(_mxarray42_);
    mxDestroyArray(_mxarray41_);
    mxDestroyArray(_mxarray39_);
    mxDestroyArray(_mxarray38_);
    mxDestroyArray(_mxarray36_);
    mxDestroyArray(_mxarray35_);
    mxDestroyArray(_mxarray33_);
    mxDestroyArray(_mxarray31_);
    mxDestroyArray(_mxarray29_);
    mxDestroyArray(_mxarray28_);
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray20_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static void Morthospm3(mxArray * threshold,
                       mxArray * onsets,
                       mxArray * window,
                       mxArray * spm_file,
                       mxArray * anat_file,
                       mxArray * tseries_path,
                       mxArray * func_root);

_mexLocalFunctionTable _local_function_table_orthospm3
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfOrthospm3" contains the normal interface for the
 * "orthospm3" M-function from file
 * "/net/quentin/home/hernan/matlab/img/orthospm3.m" (lines 1-222). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlfOrthospm3(mxArray * threshold,
                  mxArray * onsets,
                  mxArray * window,
                  mxArray * spm_file,
                  mxArray * anat_file,
                  mxArray * tseries_path,
                  mxArray * func_root) {
    mlfEnterNewContext(
      0,
      7,
      threshold,
      onsets,
      window,
      spm_file,
      anat_file,
      tseries_path,
      func_root);
    Morthospm3(
      threshold, onsets, window, spm_file, anat_file, tseries_path, func_root);
    mlfRestorePreviousContext(
      0,
      7,
      threshold,
      onsets,
      window,
      spm_file,
      anat_file,
      tseries_path,
      func_root);
}

/*
 * The function "mlxOrthospm3" contains the feval interface for the "orthospm3"
 * M-function from file "/net/quentin/home/hernan/matlab/img/orthospm3.m"
 * (lines 1-222). The feval function calls the implementation version of
 * orthospm3 through this function. This function processes any input arguments
 * and passes them to the implementation version of the function, appearing
 * above.
 */
void mlxOrthospm3(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[7];
    int i;
    if (nlhs > 0) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: orthospm3 Line: 1 Column:"
            " 1 The function \"orthospm3\" was called with m"
            "ore than the declared number of outputs (0)."),
          NULL);
    }
    if (nrhs > 7) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: orthospm3 Line: 1 Column:"
            " 1 The function \"orthospm3\" was called with m"
            "ore than the declared number of inputs (7)."),
          NULL);
    }
    for (i = 0; i < 7 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 7; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0,
      7,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6]);
    Morthospm3(
      mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5], mprhs[6]);
    mlfRestorePreviousContext(
      0,
      7,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6]);
}

/*
 * The function "Morthospm3" is the implementation version of the "orthospm3"
 * M-function from file "/net/quentin/home/hernan/matlab/img/orthospm3.m"
 * (lines 1-222). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function orthospm3(threshold, onsets, window ,spm_file, anat_file, tseries_path, func_root )
 */
static void Morthospm3(mxArray * threshold,
                       mxArray * onsets,
                       mxArray * window,
                       mxArray * spm_file,
                       mxArray * anat_file,
                       mxArray * tseries_path,
                       mxArray * func_root) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_orthospm3);
    int nargin_
      = mclNargin(
          7,
          threshold,
          onsets,
          window,
          spm_file,
          anat_file,
          tseries_path,
          func_root,
          NULL);
    mxArray * events = NULL;
    mxArray * ev_std = NULL;
    mxArray * ev_avg = NULL;
    mxArray * n = NULL;
    mxArray * tdata = NULL;
    mxArray * str = NULL;
    mxArray * zs = NULL;
    mxArray * ys = NULL;
    mxArray * xs = NULL;
    mxArray * fig = NULL;
    mxArray * button = NULL;
    mxArray * j = NULL;
    mxArray * i = NULL;
    mxArray * numRuns = NULL;
    mxArray * p = NULL;
    mxArray * op = NULL;
    mxArray * file = NULL;
    mxArray * f = NULL;
    mxArray * fig3 = NULL;
    mxArray * fig2 = NULL;
    mxArray * fig1 = NULL;
    mxArray * stretch = NULL;
    mxArray * d = NULL;
    mxArray * mymap = NULL;
    mxArray * tmp = NULL;
    mxArray * myhot = NULL;
    mxArray * mygray = NULL;
    mxArray * out_data = NULL;
    mxArray * spm_data2 = NULL;
    mxArray * zi = NULL;
    mxArray * yi = NULL;
    mxArray * xi = NULL;
    mxArray * z = NULL;
    mxArray * y = NULL;
    mxArray * x = NULL;
    mxArray * anat_data = NULL;
    mxArray * hdr = NULL;
    mxArray * spm_scale = NULL;
    mxArray * spm_data = NULL;
    mxArray * spm_hdr = NULL;
    mxArray * hdrname = NULL;
    mxArray * imgname = NULL;
    mxArray * sz = NULL;
    mxArray * path = NULL;
    mxArray * name = NULL;
    mxArray * roi = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&threshold);
    mclCopyArray(&onsets);
    mclCopyArray(&window);
    mclCopyArray(&spm_file);
    mclCopyArray(&anat_file);
    mclCopyArray(&tseries_path);
    mclCopyArray(&func_root);
    /*
     * % function result = orthospm3( threshold , onsets, window [,spm_file , anat_file, tseries_path, func_root_name] )
     * %
     * % this function takes a file in analyze format and 
     * % displays orthogonal sections thru the planes interscting
     * % at x,y,z.  It then overlays a Statistical Parameter Map 
     * % thresholded at "threshold".
     * %
     * % additionally, this function allows you to extract a time series
     * % from a chose pixel, just by clicking on it.  The time series is saved
     * % in the file "tdata.dat" if you use the RIGHT MOUSE BUTTON.
     * % This file gets overwritten everytime you select a new pixel.
     * % 
     * % *** This version of the program calls event_avg.m 
     * % in order to average all the events into one.  No deconvolution!!
     * %
     * % onsets and window must be in scan units
     * 
     * global SPM_scale_factor 
     * close all
     */
    mclPrintAns(&ans, mlfNClose(0, _mxarray0_, NULL));
    /*
     * 
     * roi =0;  % this is the number of pixels around the pixel of interest
     */
    mlfAssign(&roi, _mxarray2_);
    /*
     * % that we'll include in the time series avg.
     * if nargin < 4
     */
    if (nargin_ < 4) {
        /*
         * [name path] = uigetfile('*.img','Select SPM *.img file');
         */
        mlfAssign(
          &name, mlfUigetfile(&path, _mxarray3_, _mxarray5_, NULL, NULL));
        /*
         * spm_file = strcat(path,name)
         */
        mlfAssign(
          &spm_file, mlfStrcat(mclVv(path, "path"), mclVv(name, "name"), NULL));
        mclPrintArray(mclVa(spm_file, "spm_file"), "spm_file");
    /*
     * end
     */
    }
    /*
     * 
     * % figure out name for spm file
     * sz = size(spm_file);
     */
    mlfAssign(
      &sz, mlfSize(mclValueVarargout(), mclVa(spm_file, "spm_file"), NULL));
    /*
     * imgname = strcat(spm_file(1,1:sz(2)-4) , '.img');
     */
    mlfAssign(
      &imgname,
      mlfStrcat(
        mclArrayRef2(
          mclVa(spm_file, "spm_file"),
          _mxarray7_,
          mlfColon(
            _mxarray7_,
            mclMinus(mclIntArrayRef1(mclVv(sz, "sz"), 2), _mxarray8_),
            NULL)),
        _mxarray9_,
        NULL));
    /*
     * hdrname = strcat(spm_file(1,1:sz(2)-4) , '.hdr');
     */
    mlfAssign(
      &hdrname,
      mlfStrcat(
        mclArrayRef2(
          mclVa(spm_file, "spm_file"),
          _mxarray7_,
          mlfColon(
            _mxarray7_,
            mclMinus(mclIntArrayRef1(mclVv(sz, "sz"), 2), _mxarray8_),
            NULL)),
        _mxarray11_,
        NULL));
    /*
     * %label=sprintf('FUNCTIONAL:  ....%s',imgname(end-15:end));
     * 
     * spm_hdr = read_hdr(hdrname);
     */
    mlfAssign(&spm_hdr, mlfRead_hdr(mclVv(hdrname, "hdrname")));
    /*
     * spm_data = read_img2(spm_hdr,imgname);
     */
    mlfAssign(
      &spm_data,
      mlfRead_img2(mclVv(spm_hdr, "spm_hdr"), mclVv(imgname, "imgname")));
    /*
     * spm_data(find(spm_data==NaN))=0;
     */
    mclArrayAssign1(
      &spm_data,
      _mxarray2_,
      mlfFind(NULL, NULL, mclEq(mclVv(spm_data, "spm_data"), _mxarray13_)));
    /*
     * spm_scale =  SPM_scale_factor;
     */
    mlfAssign(&spm_scale, mclVg(&SPM_scale_factor, "SPM_scale_factor"));
    /*
     * 
     * if nargin < 4
     */
    if (nargin_ < 4) {
        /*
         * [name path] = uigetfile('*.img','Select anatomical *.img file');
         */
        mlfAssign(
          &name, mlfUigetfile(&path, _mxarray3_, _mxarray14_, NULL, NULL));
        /*
         * anat_file= strcat(path,name)
         */
        mlfAssign(
          &anat_file,
          mlfStrcat(mclVv(path, "path"), mclVv(name, "name"), NULL));
        mclPrintArray(mclVa(anat_file, "anat_file"), "anat_file");
    /*
     * end
     */
    }
    /*
     * 
     * % figure out name for anatomical file
     * sz = size(anat_file);
     */
    mlfAssign(
      &sz, mlfSize(mclValueVarargout(), mclVa(anat_file, "anat_file"), NULL));
    /*
     * 
     * imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
     */
    mlfAssign(
      &imgname,
      mlfStrcat(
        mclArrayRef2(
          mclVa(anat_file, "anat_file"),
          _mxarray7_,
          mlfColon(
            _mxarray7_,
            mclMinus(mclIntArrayRef1(mclVv(sz, "sz"), 2), _mxarray8_),
            NULL)),
        _mxarray9_,
        NULL));
    /*
     * hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');
     */
    mlfAssign(
      &hdrname,
      mlfStrcat(
        mclArrayRef2(
          mclVa(anat_file, "anat_file"),
          _mxarray7_,
          mlfColon(
            _mxarray7_,
            mclMinus(mclIntArrayRef1(mclVv(sz, "sz"), 2), _mxarray8_),
            NULL)),
        _mxarray11_,
        NULL));
    /*
     * %label=sprintf('%s\nANATOMY: ....%s', label, imgname(end-15:end));
     * 
     * hdr = read_hdr(hdrname);
     */
    mlfAssign(&hdr, mlfRead_hdr(mclVv(hdrname, "hdrname")));
    /*
     * anat_data = read_img2(hdr,imgname);
     */
    mlfAssign(
      &anat_data, mlfRead_img2(mclVv(hdr, "hdr"), mclVv(imgname, "imgname")));
    /*
     * 
     * % interpolate the paramter map to fit the anatomical one
     * % note the transpose....
     * [x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
     */
    mlfAssign(
      &x,
      mlfNMeshgrid(
        3,
        &y,
        &z,
        mclFeval(
          mclValueVarargout(),
          mlxColon,
          _mxarray7_,
          mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".ydim"),
          NULL),
        mclFeval(
          mclValueVarargout(),
          mlxColon,
          _mxarray7_,
          mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".xdim"),
          NULL),
        mclFeval(
          mclValueVarargout(),
          mlxColon,
          _mxarray7_,
          mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".zdim"),
          NULL)));
    /*
     * [xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);
     */
    mlfAssign(
      &xi,
      mlfNMeshgrid(
        3,
        &yi,
        &zi,
        mclFeval(
          mclValueVarargout(),
          mlxColon,
          _mxarray7_,
          mlfIndexRef(mclVv(hdr, "hdr"), ".ydim"),
          NULL),
        mclFeval(
          mclValueVarargout(),
          mlxColon,
          _mxarray7_,
          mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
          NULL),
        mclFeval(
          mclValueVarargout(),
          mlxColon,
          _mxarray7_,
          mlfIndexRef(mclVv(hdr, "hdr"), ".zdim"),
          NULL)));
    /*
     * 
     * yi = yi * spm_hdr.ydim/hdr.ydim;
     */
    mlfAssign(
      &yi,
      mclFeval(
        mclValueVarargout(),
        mlxMrdivide,
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mclVv(yi, "yi"),
          mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".ydim"),
          NULL),
        mlfIndexRef(mclVv(hdr, "hdr"), ".ydim"),
        NULL));
    /*
     * xi = xi * spm_hdr.xdim/hdr.xdim;
     */
    mlfAssign(
      &xi,
      mclFeval(
        mclValueVarargout(),
        mlxMrdivide,
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mclVv(xi, "xi"),
          mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".xdim"),
          NULL),
        mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
        NULL));
    /*
     * zi = zi * spm_hdr.zdim/hdr.zdim;
     */
    mlfAssign(
      &zi,
      mclFeval(
        mclValueVarargout(),
        mlxMrdivide,
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mclVv(zi, "zi"),
          mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".zdim"),
          NULL),
        mlfIndexRef(mclVv(hdr, "hdr"), ".zdim"),
        NULL));
    /*
     * 
     * spm_data2 = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
     */
    mlfAssign(
      &spm_data2,
      mlfInterp3(
        mclVv(x, "x"),
        mclVv(y, "y"),
        mclVv(z, "z"),
        mclVv(spm_data, "spm_data"),
        mclVv(xi, "xi"),
        mclVv(yi, "yi"),
        mclVv(zi, "zi"),
        _mxarray16_,
        NULL));
    /*
     * spm_data2(find(isnan(spm_data2)))=0;
     */
    mclArrayAssign1(
      &spm_data2,
      _mxarray2_,
      mlfFind(NULL, NULL, mlfIsnan(mclVv(spm_data2, "spm_data2"))));
    /*
     * %spm_data = spm_data2;
     * 
     * %whos
     * 
     * % threshold the map at the 50%
     * if nargin ==0
     */
    if (nargin_ == 0) {
        /*
         * threshold = mean(mean(mean(spm_data2))) * 3
         */
        mlfAssign(
          &threshold,
          mclMtimes(
            mlfMean(
              mlfMean(mlfMean(mclVv(spm_data2, "spm_data2"), NULL), NULL),
              NULL),
            _mxarray18_));
        mclPrintArray(mclVa(threshold, "threshold"), "threshold");
    /*
     * end
     */
    }
    /*
     * spm_data2(find(spm_data2 <= threshold )) = 0;
     */
    mclArrayAssign1(
      &spm_data2,
      _mxarray2_,
      mlfFind(
        NULL,
        NULL,
        mclLe(mclVv(spm_data2, "spm_data2"), mclVa(threshold, "threshold"))));
    /*
     * 
     * % scale the maps to use the whole colormap
     * anat_data = anat_data * 2^7 / max(max(max(anat_data)));
     */
    mlfAssign(
      &anat_data,
      mclMrdivide(
        mclMtimes(mclVv(anat_data, "anat_data"), _mxarray19_),
        mlfMax(
          NULL,
          mlfMax(
            NULL,
            mlfMax(NULL, mclVv(anat_data, "anat_data"), NULL, NULL),
            NULL,
            NULL),
          NULL,
          NULL)));
    /*
     * spm_data2 = spm_data2 * 2^7 / max(max(max(spm_data2)));
     */
    mlfAssign(
      &spm_data2,
      mclMrdivide(
        mclMtimes(mclVv(spm_data2, "spm_data2"), _mxarray19_),
        mlfMax(
          NULL,
          mlfMax(
            NULL,
            mlfMax(NULL, mclVv(spm_data2, "spm_data2"), NULL, NULL),
            NULL,
            NULL),
          NULL,
          NULL)));
    /*
     * 
     * out_data = anat_data;
     */
    mlfAssign(&out_data, mclVv(anat_data, "anat_data"));
    /*
     * out_data(find(spm_data2)) =  2^7 + spm_data2(find(spm_data2));
     */
    mclArrayAssign1(
      &out_data,
      mclPlus(
        _mxarray19_,
        mclArrayRef1(
          mclVv(spm_data2, "spm_data2"),
          mlfFind(NULL, NULL, mclVv(spm_data2, "spm_data2")))),
      mlfFind(NULL, NULL, mclVv(spm_data2, "spm_data2")));
    /*
     * 
     * %configure the colormap:
     * mygray = [0:1/127:1]' * [1 1 1]; 
     */
    mlfAssign(&mygray, _mxarray20_);
    /*
     * 
     * myhot = [0:3:125]' * [1 0 0] ;
     */
    mlfAssign(&myhot, _mxarray22_);
    /*
     * tmp =   [0:3:125]' * [0 1 0] ;
     */
    mlfAssign(&tmp, _mxarray24_);
    /*
     * myhot = [myhot; tmp];
     */
    mlfAssign(
      &myhot, mlfVertcat(mclVv(myhot, "myhot"), mclVv(tmp, "tmp"), NULL));
    /*
     * tmp =   [0:3:125]' * [0 0 1];
     */
    mlfAssign(&tmp, _mxarray26_);
    /*
     * myhot =  [myhot;  tmp]/128;
     */
    mlfAssign(
      &myhot,
      mclMrdivide(
        mlfVertcat(mclVv(myhot, "myhot"), mclVv(tmp, "tmp"), NULL),
        _mxarray19_));
    /*
     * 
     * myhot(round(128/3): 128, 1) = 1;
     */
    mclArrayAssign2(
      &myhot,
      _mxarray7_,
      mlfColon(
        mlfScalar(svDoubleScalarRound(42.666666666666664)), _mxarray19_, NULL),
      _mxarray7_);
    /*
     * myhot(round(128*2/3):128,2) = 1;
     */
    mclArrayAssign2(
      &myhot,
      _mxarray7_,
      mlfColon(
        mlfScalar(svDoubleScalarRound(85.33333333333333)), _mxarray19_, NULL),
      _mxarray28_);
    /*
     * 
     * mymap = [mygray; myhot];
     */
    mlfAssign(
      &mymap, mlfVertcat(mclVv(mygray, "mygray"), mclVv(myhot, "myhot"), NULL));
    /*
     * 
     * colormap(mymap)
     */
    mclPrintAns(&ans, mlfNColormap(0, mclVv(mymap, "mymap"), NULL));
    /*
     * axis off
     */
    mclPrintAns(&ans, mlfNAxis(0, NULL, NULL, _mxarray29_, NULL));
    /*
     * 
     * 
     * d=out_data;
     */
    mlfAssign(&d, mclVv(out_data, "out_data"));
    /*
     * 
     * fprintf('\nfirst display the orthogonal sections');
     */
    mclAssignAns(&ans, mlfNFprintf(0, _mxarray31_, NULL));
    /*
     * x=ceil(hdr.xdim/2);
     */
    mlfAssign(
      &x,
      mlfCeil(
        mclFeval(
          mclValueVarargout(),
          mlxMrdivide,
          mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
          _mxarray28_,
          NULL)));
    /*
     * y=ceil(hdr.ydim/2);
     */
    mlfAssign(
      &y,
      mlfCeil(
        mclFeval(
          mclValueVarargout(),
          mlxMrdivide,
          mlfIndexRef(mclVv(hdr, "hdr"), ".ydim"),
          _mxarray28_,
          NULL)));
    /*
     * z=ceil(hdr.zdim/2);
     */
    mlfAssign(
      &z,
      mlfCeil(
        mclFeval(
          mclValueVarargout(),
          mlxMrdivide,
          mlfIndexRef(mclVv(hdr, "hdr"), ".zdim"),
          _mxarray28_,
          NULL)));
    /*
     * 
     * %colordef black 	
     * stretch = hdr.zsize/hdr.xsize;
     */
    mlfAssign(
      &stretch,
      mclFeval(
        mclValueVarargout(),
        mlxMrdivide,
        mlfIndexRef(mclVv(hdr, "hdr"), ".zsize"),
        mlfIndexRef(mclVv(hdr, "hdr"), ".xsize"),
        NULL));
    /*
     * 
     * %colordef black
     * [fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);
     */
    mlfAssign(
      &fig1,
      mlfOv(
        &fig2,
        &fig3,
        mclVv(hdr, "hdr"),
        mclVv(d, "d"),
        mclVv(x, "x"),
        mclVv(y, "y"),
        mclVv(z, "z"),
        _mxarray2_));
    /*
     * 
     * 
     * if nargin < 4
     */
    if (nargin_ < 4) {
        /*
         * % determine which files contain the time series..:
         * f='dummmyname';
         */
        mlfAssign(&f, _mxarray33_);
        /*
         * file=[];
         */
        mlfAssign(&file, _mxarray35_);
        /*
         * path = [];
         */
        mlfAssign(&path, _mxarray35_);
        /*
         * op=pwd;
         */
        mlfAssign(&op, mlfPwd());
        /*
         * while f ~= 0
         */
        while (mclNeBool(mclVv(f, "f"), _mxarray2_)) {
            /*
             * [f p] = uigetfile('*.img','Select a file from the next run.  cancel when finished ');
             */
            mlfAssign(
              &f, mlfUigetfile(&p, _mxarray3_, _mxarray36_, NULL, NULL));
            /*
             * if f~=0
             */
            if (mclNeBool(mclVv(f, "f"), _mxarray2_)) {
                /*
                 * file = [file ; f];
                 */
                mlfAssign(
                  &file, mlfVertcat(mclVv(file, "file"), mclVv(f, "f"), NULL));
                /*
                 * path = [path ; p];
                 */
                mlfAssign(
                  &path, mlfVertcat(mclVv(path, "path"), mclVv(p, "p"), NULL));
                /*
                 * cd(p);
                 */
                mclAssignAns(&ans, mlfNCd(0, mclVv(p, "p")));
            /*
             * end
             */
            }
        /*
         * end
         */
        }
        /*
         * cd(op);
         */
        mclAssignAns(&ans, mlfNCd(0, mclVv(op, "op")));
        /*
         * numRuns = size(file,1);
         */
        mlfAssign(
          &numRuns,
          mlfSize(mclValueVarargout(), mclVv(file, "file"), _mxarray7_));
    /*
     * end
     */
    }
    /*
     * 
     * [i j button] = ginput(1);
     */
    mlfAssign(&i, mlfNGinput(3, &j, &button, _mxarray7_));
    /*
     * i=round(i);j=round(j);
     */
    mlfAssign(&i, mlfRound(mclVv(i, "i")));
    mlfAssign(&j, mlfRound(mclVv(j, "j")));
    /*
     * 
     * 
     * %i=1;
     * while i >= -10
     */
    while (mclGeBool(mclVv(i, "i"), _mxarray38_)) {
        /*
         * 
         * fig = floor(gca);
         */
        mlfAssign(&fig, mlfFloor(mlfGca(NULL)));
        /*
         * switch(fig)
         */
        {
            mxArray * v_ = mclInitialize(mclVv(fig, "fig"));
            if (mclSwitchCompare(v_, mlfFloor(mclVv(fig1, "fig1")))) {
                /*
                 * case floor(fig1)
                 * x=j;
                 */
                mlfAssign(&x, mclVv(j, "j"));
                /*
                 * y=i;
                 */
                mlfAssign(&y, mclVv(i, "i"));
            /*
             * case floor(fig2)
             */
            } else if (mclSwitchCompare(v_, mlfFloor(mclVv(fig2, "fig2")))) {
                /*
                 * z=j;
                 */
                mlfAssign(&z, mclVv(j, "j"));
                /*
                 * x=i;
                 */
                mlfAssign(&x, mclVv(i, "i"));
            /*
             * case floor(fig3)
             */
            } else if (mclSwitchCompare(v_, mlfFloor(mclVv(fig3, "fig3")))) {
                /*
                 * y=i;
                 */
                mlfAssign(&y, mclVv(i, "i"));
                /*
                 * z=j;
                 */
                mlfAssign(&z, mclVv(j, "j"));
            /*
             * end
             */
            }
            mxDestroyArray(v_);
        }
        /*
         * 
         * 
         * xs = ceil(x*spm_hdr.xdim/hdr.xdim);
         */
        mlfAssign(
          &xs,
          mlfCeil(
            mclFeval(
              mclValueVarargout(),
              mlxMrdivide,
              mclFeval(
                mclValueVarargout(),
                mlxMtimes,
                mclVv(x, "x"),
                mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".xdim"),
                NULL),
              mlfIndexRef(mclVv(hdr, "hdr"), ".xdim"),
              NULL)));
        /*
         * ys = ceil(y*spm_hdr.ydim/hdr.ydim);
         */
        mlfAssign(
          &ys,
          mlfCeil(
            mclFeval(
              mclValueVarargout(),
              mlxMrdivide,
              mclFeval(
                mclValueVarargout(),
                mlxMtimes,
                mclVv(y, "y"),
                mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".ydim"),
                NULL),
              mlfIndexRef(mclVv(hdr, "hdr"), ".ydim"),
              NULL)));
        /*
         * zs = ceil(z*spm_hdr.zdim/hdr.zdim);
         */
        mlfAssign(
          &zs,
          mlfCeil(
            mclFeval(
              mclValueVarargout(),
              mlxMrdivide,
              mclFeval(
                mclValueVarargout(),
                mlxMtimes,
                mclVv(z, "z"),
                mlfIndexRef(mclVv(spm_hdr, "spm_hdr"), ".zdim"),
                NULL),
              mlfIndexRef(mclVv(hdr, "hdr"), ".zdim"),
              NULL)));
        /*
         * 
         * str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm,  \n (%3d %3d %3d)vx , \n val= %d',...
         */
        mlfAssign(
          &str,
          mlfSprintf(
            NULL,
            _mxarray39_,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfIndexRef(mclVv(hdr, "hdr"), ".xsize"),
              mclFeval(
                mclValueVarargout(),
                mlxMinus,
                mclVv(xs, "xs"),
                mlfIndexRef(
                  mclVv(spm_hdr, "spm_hdr"), ".origin(?)", _mxarray7_),
                NULL),
              NULL),
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfIndexRef(mclVv(hdr, "hdr"), ".ysize"),
              mclFeval(
                mclValueVarargout(),
                mlxMinus,
                mclVv(ys, "ys"),
                mlfIndexRef(
                  mclVv(spm_hdr, "spm_hdr"), ".origin(?)", _mxarray28_),
                NULL),
              NULL),
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfIndexRef(mclVv(hdr, "hdr"), ".zsize"),
              mclFeval(
                mclValueVarargout(),
                mlxMinus,
                mclVv(zs, "zs"),
                mlfIndexRef(
                  mclVv(spm_hdr, "spm_hdr"), ".origin(?)", _mxarray18_),
                NULL),
              NULL),
            mclVv(xs, "xs"),
            mclVv(ys, "ys"),
            mclVv(zs, "zs"),
            mclMtimes(
              mlfIndexRef(
                mclVv(spm_data, "spm_data"),
                "(?,?,?)",
                mclVv(xs, "xs"),
                mclVv(ys, "ys"),
                mclVv(zs, "zs")),
              mclVv(spm_scale, "spm_scale")),
            NULL));
        /*
         * hdr.xsize*(xs-spm_hdr.origin(1)), ...
         * hdr.ysize*(ys-spm_hdr.origin(2)), ...
         * hdr.zsize*(zs-spm_hdr.origin(3)), ...
         * xs,ys,zs, ...
         * spm_data(xs,ys,zs)*spm_scale );
         * 
         * %colordef black
         * [ x y z] 
         */
        mclPrintAns(
          &ans, mlfHorzcat(mclVv(x, "x"), mclVv(y, "y"), mclVv(z, "z"), NULL));
        /*
         * [fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);
         */
        mlfAssign(
          &fig1,
          mlfOv(
            &fig2,
            &fig3,
            mclVv(hdr, "hdr"),
            mclVv(d, "d"),
            mclVv(x, "x"),
            mclVv(y, "y"),
            mclVv(z, "z"),
            _mxarray2_));
        /*
         * subplot(221), title(str) %, xlabel(label)
         */
        mclPrintAns(&ans, mlfNSubplot(0, _mxarray41_, NULL, NULL, NULL));
        mclPrintAns(&ans, mlfNTitle(0, mclVv(str, "str"), NULL));
        /*
         * 
         * %        pos = length(path)-8;
         * %         tdata = [];
         * %         for k=1:3
         * %             path(pos) = num2str(k)
         * %             tmp = timeplot2(path, file,[x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);
         * %             tdata=[tdata; tmp];
         * %             
         * %         end
         * 
         * tmp=mean(mean(spm_data(xs-roi:xs+roi, ys-roi:ys+roi)));
         */
        mlfAssign(
          &tmp,
          mlfMean(
            mlfMean(
              mclArrayRef2(
                mclVv(spm_data, "spm_data"),
                mlfColon(
                  mclMinus(mclVv(xs, "xs"), mclVv(roi, "roi")),
                  mclPlus(mclVv(xs, "xs"), mclVv(roi, "roi")),
                  NULL),
                mlfColon(
                  mclMinus(mclVv(ys, "ys"), mclVv(roi, "roi")),
                  mclPlus(mclVv(ys, "ys"), mclVv(roi, "roi")),
                  NULL)),
              NULL),
            NULL));
        /*
         * 
         * if nargin < 4
         */
        if (nargin_ < 4) {
            /*
             * tdata=[];
             */
            mlfAssign(&tdata, _mxarray35_);
            /*
             * for n=1:numRuns
             */
            {
                int v_ = mclForIntStart(1);
                int e_ = mclForIntEnd(mclVv(numRuns, "numRuns"));
                if (v_ > e_) {
                    mlfAssign(&n, _mxarray35_);
                } else {
                    /*
                     * tmp = timeplot2(path(n,:), file(n,:),[xs-roi xs+roi],[ys-roi ys+roi],[zs-roi zs+roi]);
                     * tmp = (tmp -tmp(1))./ mean(tmp);
                     * tdata = [tdata tmp];
                     * end
                     */
                    for (; ; ) {
                        mlfAssign(
                          &tmp,
                          mlfTimeplot2(
                            mclArrayRef2(
                              mclVv(path, "path"),
                              mlfScalar(v_),
                              mlfCreateColonIndex()),
                            mclArrayRef2(
                              mclVv(file, "file"),
                              mlfScalar(v_),
                              mlfCreateColonIndex()),
                            mlfHorzcat(
                              mclMinus(mclVv(xs, "xs"), mclVv(roi, "roi")),
                              mclPlus(mclVv(xs, "xs"), mclVv(roi, "roi")),
                              NULL),
                            mlfHorzcat(
                              mclMinus(mclVv(ys, "ys"), mclVv(roi, "roi")),
                              mclPlus(mclVv(ys, "ys"), mclVv(roi, "roi")),
                              NULL),
                            mlfHorzcat(
                              mclMinus(mclVv(zs, "zs"), mclVv(roi, "roi")),
                              mclPlus(mclVv(zs, "zs"), mclVv(roi, "roi")),
                              NULL)));
                        mlfAssign(
                          &tmp,
                          mclRdivide(
                            mclMinus(
                              mclVv(tmp, "tmp"),
                              mclIntArrayRef1(mclVv(tmp, "tmp"), 1)),
                            mlfMean(mclVv(tmp, "tmp"), NULL)));
                        mlfAssign(
                          &tdata,
                          mlfHorzcat(
                            mclVv(tdata, "tdata"), mclVv(tmp, "tmp"), NULL));
                        if (v_ == e_) {
                            break;
                        }
                        ++v_;
                    }
                    mlfAssign(&n, mlfScalar(v_));
                }
            }
        /*
         * else 
         */
        } else {
            /*
             * path=tseries_path;
             */
            mlfAssign(&path, mclVa(tseries_path, "tseries_path"));
            /*
             * file=func_root;
             */
            mlfAssign(&file, mclVa(func_root, "func_root"));
            /*
             * %file = 'sravol_exxxxxx';
             * tdata=[];
             */
            mlfAssign(&tdata, _mxarray35_);
            /*
             * 
             * tmp = timeplot2(path, file,[xs-roi xs+roi],[ys-roi ys+roi],[zs-roi zs+roi]);
             */
            mlfAssign(
              &tmp,
              mlfTimeplot2(
                mclVv(path, "path"),
                mclVv(file, "file"),
                mlfHorzcat(
                  mclMinus(mclVv(xs, "xs"), mclVv(roi, "roi")),
                  mclPlus(mclVv(xs, "xs"), mclVv(roi, "roi")),
                  NULL),
                mlfHorzcat(
                  mclMinus(mclVv(ys, "ys"), mclVv(roi, "roi")),
                  mclPlus(mclVv(ys, "ys"), mclVv(roi, "roi")),
                  NULL),
                mlfHorzcat(
                  mclMinus(mclVv(zs, "zs"), mclVv(roi, "roi")),
                  mclPlus(mclVv(zs, "zs"), mclVv(roi, "roi")),
                  NULL)));
            /*
             * tmp = (tmp -tmp(1))./ mean(tmp);
             */
            mlfAssign(
              &tmp,
              mclRdivide(
                mclMinus(
                  mclVv(tmp, "tmp"), mclIntArrayRef1(mclVv(tmp, "tmp"), 1)),
                mlfMean(mclVv(tmp, "tmp"), NULL)));
            /*
             * tdata = [tdata tmp];
             */
            mlfAssign(
              &tdata,
              mlfHorzcat(mclVv(tdata, "tdata"), mclVv(tmp, "tmp"), NULL));
        /*
         * end
         */
        }
        /*
         * 
         * whos events ev_avg tdata
         */
        mclPrintAns(&ans, mlfWhos(_mxarray42_, _mxarray44_, _mxarray46_, NULL));
        /*
         * [ev_avg ev_std] = event_avg(tdata,onsets,window,1);
         */
        mlfAssign(
          &ev_avg,
          mlfEvent_avg(
            &ev_std,
            mclVv(tdata, "tdata"),
            mclVa(onsets, "onsets"),
            mclVa(window, "window"),
            _mxarray7_));
        /*
         * whos ev*
         */
        mclPrintAns(&ans, mlfWhos(_mxarray48_, NULL));
        /*
         * 
         * 
         * %subplot 224, plot(events','y.');hold on
         * %subplot 224, plot(10*spm_hrf(2),'y'),hold on
         * %keyboard
         * % subplot 224, plot(ev_avg,'r');axis tight ; hold off; 
         * subplot 224, errorbar(ev_avg,ev_std);axis tight ; hold off; 
         */
        mclPrintAns(&ans, mlfNSubplot(0, _mxarray50_, NULL, NULL, NULL));
        mclAssignAns(
          &ans,
          mlfNErrorbar(
            0,
            mclVv(ev_avg, "ev_avg"),
            mclVv(ev_std, "ev_std"),
            NULL,
            NULL,
            NULL));
        mclAssignAns(&ans, mlfNAxis(0, NULL, NULL, _mxarray52_, NULL));
        mlfHold(_mxarray29_);
        /*
         * set(gca, 'Xtick',[0:2:window])
         */
        mclPrintAns(
          &ans,
          mlfNSet(
            0,
            mlfGca(NULL),
            _mxarray54_,
            mlfColon(_mxarray2_, _mxarray28_, mclVa(window, "window")),
            NULL));
        /*
         * 
         * if button==3
         */
        if (mclEqBool(mclVv(button, "button"), _mxarray18_)) {
            /*
             * save avgResponse.dat events -ASCII
             */
            mlfSave(_mxarray56_, "w", "events", events, NULL);
            /*
             * save tdata.dat tdata -ASCII
             */
            mlfSave(_mxarray58_, "w", "tdata", tdata, NULL);
        /*
         * end
         */
        }
        /*
         * 
         * [i j button] = ginput(1);
         */
        mlfAssign(&i, mlfNGinput(3, &j, &button, _mxarray7_));
        /*
         * i=round(i);j=round(j);
         */
        mlfAssign(&i, mlfRound(mclVv(i, "i")));
        mlfAssign(&j, mlfRound(mclVv(j, "j")));
    /*
     * 
     * end 
     */
    }
    mxDestroyArray(ans);
    mxDestroyArray(roi);
    mxDestroyArray(name);
    mxDestroyArray(path);
    mxDestroyArray(sz);
    mxDestroyArray(imgname);
    mxDestroyArray(hdrname);
    mxDestroyArray(spm_hdr);
    mxDestroyArray(spm_data);
    mxDestroyArray(spm_scale);
    mxDestroyArray(hdr);
    mxDestroyArray(anat_data);
    mxDestroyArray(x);
    mxDestroyArray(y);
    mxDestroyArray(z);
    mxDestroyArray(xi);
    mxDestroyArray(yi);
    mxDestroyArray(zi);
    mxDestroyArray(spm_data2);
    mxDestroyArray(out_data);
    mxDestroyArray(mygray);
    mxDestroyArray(myhot);
    mxDestroyArray(tmp);
    mxDestroyArray(mymap);
    mxDestroyArray(d);
    mxDestroyArray(stretch);
    mxDestroyArray(fig1);
    mxDestroyArray(fig2);
    mxDestroyArray(fig3);
    mxDestroyArray(f);
    mxDestroyArray(file);
    mxDestroyArray(op);
    mxDestroyArray(p);
    mxDestroyArray(numRuns);
    mxDestroyArray(i);
    mxDestroyArray(j);
    mxDestroyArray(button);
    mxDestroyArray(fig);
    mxDestroyArray(xs);
    mxDestroyArray(ys);
    mxDestroyArray(zs);
    mxDestroyArray(str);
    mxDestroyArray(tdata);
    mxDestroyArray(n);
    mxDestroyArray(ev_avg);
    mxDestroyArray(ev_std);
    mxDestroyArray(events);
    mxDestroyArray(func_root);
    mxDestroyArray(tseries_path);
    mxDestroyArray(anat_file);
    mxDestroyArray(spm_file);
    mxDestroyArray(window);
    mxDestroyArray(onsets);
    mxDestroyArray(threshold);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    /*
     * %colordef white
     * return
     */
}
