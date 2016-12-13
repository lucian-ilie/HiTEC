#include<math.h>
using namespace std;

inline bool leq(int a1, int a2, int b1, int b2) // lexicographic order 
{ return(a1 < b1 || a1 == b1 && a2 <= b2); }		// for pairs
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) 
{ return(a1 < b1 || a1 == b1 && leq(a2, a3, b2, b3)); }	// and triples

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int* a, int* b, int* r, int n, int K)
{ // count occurrences
	int* c = new int[K+1], i=0, sum=0, t=0;		// counter array
	for (i=0; i<=K; i++) c[i] = 0;				// reset counters
	for (i=0; i<n; i++) c[r[a[i]]]++;			// count occurrences
	for (i=0, sum=0; i<=K; i++)						// exclusive prefix sums
	{	t = c[i]; c[i] = sum; sum += t; }
	for (i=0; i<n; i++) b[c[r[a[i]]]++] = a[i];	// sort
	delete [] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}^n,    // require s[n] = s[n+1] = s[n+2] = 0, n>=2
void sa_for_ints(int* s, int* SA, int n, int K)
{	int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2, i=0, j=0, k=0, p=0, t=0;
	int* s12  = new int[n02+3];  s12[n02] =  s12[n02+1] =  s12[n02+2] = 0;
	int* SA12 = new int[n02+3]; SA12[n02] = SA12[n02+1] = SA12[n02+2] = 0;
	int* s0   = new int[n0];
	int* SA0  = new int[n0];

	// generate positions of mod 1 and mod 2 suffixes
	// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
	for (i=0, j=0; i<n+(n0-n1); i++) if (i%3 != 0) s12[j++] = i;

	// lsb radix sort the mod 1 and mod 2 triples
	radixPass(s12, SA12, s+2, n02, K);
	radixPass(SA12, s12, s+1, n02, K);
	radixPass(s12, SA12, s,   n02, K);

	// find lexicographic names of triples
	int name = 0, c0 = -1, c1 = -1, c2 = -1;
	for (i=0; i<n02; i++)
	{	if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2)
		{ name++; c0 = s[SA12[i]]; c1 = s[SA12[i]+1]; c2 = s[SA12[i]+2]; }
		if (SA12[i] % 3 == 1) { s12[SA12[i]/3]     = name; }  // left half
		else 				  { s12[SA12[i]/3 + n0] = name; }  // right half
	}

	// recurse if names are not yet unique
	if (name < n02) {
		sa_for_ints(s12, SA12, n02, name);
		// store unique names in s12 using the suffix array
		for (i=0; i<n02; i++) s12[SA12[i]] = i+1;
	} else // generate the suffix array of s12 directly
		for (i=0; i<n02; i++) SA12[s12[i]-1] = i;

	// stably sort the mod 0 suffixes from SA12 by their first character
	for (i=0, j=0; i<n02; i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
	radixPass(s0, SA0, s, n0, K);

	// merge sorted SA0 suffixes and sorted SA12 suffixes
	for (p=0, t=n0-n1, k=0; k<n; k++) {
		#define GetI() (SA12[t] < n0 ? SA12[t]*3+1 : (SA12[t]-n0)*3+2)
		i = GetI(); // pos of current offset 12 suffix
		j = SA0[p]; // pos of current offset 0  suffix
		if (SA12[t] < n0 ?  // different compares for mod 1 and mod 2 suffixes
			leq(s[i],         s12[SA12[t]+n0],   s[j],         s12[j/3]) :
			leq(s[i], s[i+1], s12[SA12[t]-n0+1], s[j], s[j+1], s12[j/3+n0]))
		{  // suffix from SA12 is smaller
			SA[k] = i; t++;
			if (t == n02)  // done --- only SA0 suffixes left
				for (k++; p<n0; p++, k++) SA[k] = SA0[p];
		} else { // suffix from SA0 is smaller
			SA[k] = j; p++;
			if (p == n0)  // done --- only SA12 suffixes left
				for (k++; t<n02; t++, k++) SA[k] = GetI();
		}
	}
	delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}

void SuffixArrayKS(unsigned char* text, int* SA, int n) // computes SA
{
	// convert to integer string
	int K = 0, i=0;
	int* s = new int[n+3]; s[n] = s[n+1] = s[n+2] = 0; 
	for (i=0; i<n; i++) 
	{	s[i] = (int)text[i]; if (s[i] > K) K = s[i]; }
	// find suffix array SA 
	// ASCII order: # $ .. 0 1 .. 9 .. A B .. Z .. a b .. z
	sa_for_ints(s, SA, n, K);
}


void SuffixArraySelPos(unsigned char* text, int* SA, unsigned char* LCP, int nPos) {  
		// compute the suffix array and LCP for nPos positions in text
		// SA must have the same size nPos and MUST contain the positions in text to be sorted
		// SA[i] = i for all 0<=i<=nPos-1 = strlen(text)-1 to compute the entire SA
		// quicksort the positions; LCP stops at z
	int low = 0, high = nPos-1, stackIndex = -1, i=0, tmp=0, less=0, curPos=0, otherPos=0, maxlcp1=0, maxlcp2=0;			
	int SAcurPos, SAotherPos, maxStackSize;			// stack holds pairs of positions (low,high) for quicksort
	maxStackSize = 10*(int)log(nPos);
	int *stack = new int[maxStackSize];				// stackIndex = -1  <=>  stack is empty
	int *stack2, *tmpStack;
	stack[0] = 0; stack[1] = nPos-1; stackIndex = 1;	// initially stack contains (0,nPos-1)
	LCP[0] = (unsigned char)0;
	while (stackIndex > -1) {
		low = stack[stackIndex-1];  high = stack[stackIndex]; stackIndex-=2;		// pop stack
		curPos = low; otherPos = high; //cout << "pop (" << low << ","  << high << ")\n";
//for (i=0; i<nPos; i++) cout << SA[i] << " "; cout << endl;
		maxlcp1 = maxlcp2 = 0;		// maxlcp1(2) stores lcps between seq at curPos and smaller (larger) ones
		while (curPos != otherPos) {
			i = 0;			// compare text[SA[curPos] ...] ? text[SA[otherPos] ...]; less = -1,0,1 for <,=,>
			SAcurPos = SA[curPos]; SAotherPos = SA[otherPos];
			while ((text[SAcurPos+i] != 'z') && (text[SAotherPos+i] != 'z') && (text[SAcurPos+i] == text[SAotherPos+i]))
				i++;
			if (text[SAcurPos+i] == 'z') 
				if (text[SAotherPos+i] == 'z') less = 0;
				else less = 1;
			else if (text[SAotherPos+i] == 'z') less = -1; 
			else if (text[SAcurPos+i] < text[SAotherPos+i]) less = -1; 
			else less = 1; 
//cout << "less = " << less << "  curPos = " << curPos << "  otherPos = " << otherPos << endl;
			if (less == -1) maxlcp2 = (maxlcp2 < i) ? i : maxlcp2;
			else if (less == 1) maxlcp1 = (maxlcp1 < i) ? i : maxlcp1;
			else if (curPos < otherPos) maxlcp2 = (maxlcp2 < i) ? i : maxlcp2;
			else maxlcp1 = (maxlcp1 < i) ? i : maxlcp1;
			if (curPos < otherPos) {
				if (less == 1) {
					tmp = SA[curPos]; SA[curPos] = SA[otherPos]; SA[otherPos] = tmp;
					tmp = curPos; curPos = otherPos; otherPos = tmp+1; 
//for (i=0; i<nPos; i++) cout << SA[i] << " "; cout << endl;
				}
				else  otherPos--; //cout << "otherPos " << otherPos << " --\n"; 
			}
			else 
				if (less == -1) {
					tmp = SA[curPos]; SA[curPos] = SA[otherPos]; SA[otherPos] = tmp;
					tmp = curPos; curPos = otherPos; otherPos = tmp-1; 
//for (i=0; i<nPos; i++) cout << SA[i] << " "; cout << endl;
				}
				else otherPos++; //cout << "otherPos " << otherPos << " ++\n"; 
		}		// here curPos = otherPos
		if (curPos > low) LCP[curPos] = maxlcp1; //cout << "set LCP[" << curPos << "] = " << maxlcp1 << endl;
		if (curPos < high) LCP[curPos+1] = maxlcp2; //cout << "set LCP[" << curPos+1 << "] = " << maxlcp2 << endl;
		if (curPos > low+1) { //cout << "push (" << low << "," << curPos-1 << ")\n";
			if (stackIndex > maxStackSize - 3) {			// increase stack
				maxStackSize = 2*maxStackSize;
				stack2 = new int[maxStackSize];
				for (i=0; i<stackIndex+1; i++) stack2[i] = stack[i];
				tmpStack = stack;
				stack = stack2;
				delete[] tmpStack;
			}
			stack[stackIndex+1] = low; stack[stackIndex+2] = curPos-1;  stackIndex+=2;
		}
		if (curPos < high-1) { //cout << "push (" << curPos+1 << "," << high << ")\n";
			if (stackIndex > maxStackSize - 3) {			// increase stack
				maxStackSize = 2*maxStackSize;
				stack2 = new int[maxStackSize];
				for (i=0; i<stackIndex+1; i++) stack2[i] = stack[i];
				tmpStack = stack;
				stack = stack2;
				delete[] tmpStack;
			}
			stack[stackIndex+1] = curPos+1; stack[stackIndex+2] = high;  stackIndex+=2;
		}
	}
	delete[] stack;
}
// LCP using direct comparison; stops at 'z', that is, it assumes 'z' != 'z'
void LongestCommonPrefixDirect(unsigned char* text, int* SA, unsigned char* LCP, int l, int witLength) { // LCP up to the first z character
	int i=0, i1=0, j1=0, k=0;
	unsigned char i1c, j1c;
	LCP[0] = 0;
	for (i=1; i<l; i++) {
		k = 0; i1 = SA[i]; j1 = SA[i-1];
		while ((k <= witLength+1) && (i1+k<l) && (j1+k<l) && ((i1c = text[i1+k]) != 'z') && ((j1c = text[j1+k]) != 'z') && (text[i1+k] == text[j1+k])) 
			k++;
		LCP[i] = k;
	}
}
// LCP using the alg of Kasai et al. -- linear time but requires SAinv
void LongestCommonPrefixKasai (unsigned char *text, int* SA, unsigned char* LCP, int n)
{										
	int i=0, j=0;
		// compute ISA (inverse SA)
	int* ISA = new int[n]; 
	for (i=0; i<n; i++) ISA[SA[i]] = i;
	// find LCP
	LCP[0] = 0;
	for (i=0; i<n; i++) // compute LCP[ISA[i]]
		if (ISA[i] != 0) {
			if (i==0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while (text[i+j] == text[SA[ISA[i]-1]+j])
				j++;
			LCP[ISA[i]] = j;
		}
}


int ACTGto0123(unsigned char a) { // it is the inverse of characters (from main): ACTGfrom0123[ACTGto0123(c)] == c
	switch ( a ) {				
		case 'A': return(0);
        case 'C': return(1);
        case 'T': return(2);
        case 'G': return(3);
		default: return(4);		// for 'z'
	}
}

