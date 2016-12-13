/*extern "C" {
    #include "ds_ssort.h"
    #include "bwt_aux.h"
    #include "lcp_aux.h"
};
*/
#include<time.h>
#include<ctime>
#include<math.h>
#include<stdio.h>
#include<iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;

inline bool leq(int a1, int a2, int b1, int b2); 
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) ;
//static void radixPass(int* a, int* b, int* r, int n, int K);
void sa_for_ints(int* s, int* SA, int n, int K);
void SuffixArrayKS(unsigned char* text, int* SA, int n);
void SuffixArraySelPos(unsigned char* text, int* SA, unsigned char* LCP, int nPos);
void LongestCommonPrefixDirect(unsigned char* text, int* SA, unsigned char* LCP, int l, int witLength);
void LongestCommonPrefixKasai (unsigned char *text, int* SA, unsigned char* LCP, int n);
int ACTGto0123(unsigned char a);

