#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef STATS
#define STATS 0
#endif

#include <cstring>
#include <fstream>

#include <sys/time.h>
#include <sys/resource.h>

#include "work.hh"
#include "divsufsort.h"
using namespace std;

template<class T>
inline void zaparr(T & x) {

//	{
//		assert(x != NULL);
//	}
	delete[] x;
	x = NULL;
}

double dbinom(int64_t k, int64_t n, double p){
	return gsl_cdf_binomial_P(k, p, n)- gsl_cdf_binomial_P(k-1, p, n);
}

//########################################
//### f(k,l,w) = no. of ways k errors can be placed in an l-read s.t. MDE < w
//### MDE = Maximum Distance Between Errors
//### exp.no.reads with k err.and MDE < w = f(k,l,w) * p^k * (1-p)^(l-k) * n
//########################################

double f(int64_t k, int64_t l, int64_t w) {
	double **fw = new double*[k + 2];
	for (int i = 0; i <= k + 1; i++)
		fw[i] = new double[l + 2];
	int64_t j = 1;
	while (j <= w) { fw[1][j] = 1; j = j + 1; }
	j = w + 1;
	while (j <= l + 1) { fw[1][j] = 0; j = j + 1; }
	int i = 2;
	while (i <= k + 1) {
		j = 1;
		while (j <= w + 1) {
			//fw[i][j] = (gsl_cdf_binomial_P(i - 1, 0.5, j - 1)-gsl_cdf_binomial_P(i - 2, 0.5, j - 1)) * 2 ^ (j - 1);
			fw[i][j] = dbinom(i-1,j-1,.5)* pow(2, (j-1));
			j = j + 1;
		}
		j = w + 2;
		while (j <= l + 1) {
			double sum = 0;
			for (int j1 = 1; j1 <= w; j1++) {
				sum = sum + fw[i - 1][j - j1];
			}
			fw[i][j] = sum;
			j = j + 1;
		}
		i = i + 1;
	}
	return fw[k + 1][l + 1];
}

//#######################################
//### U(w,l,n,p) gives the total number of reads uncorrectable with witLength = w
//#######################################
double U(int64_t w, int64_t l, int64_t n, double p) {
	double total = 0;
	int j = 1;
	int finish = 0; //### finish becomes 1 when f(j,l,w) is v.small
	while (finish == 0) {
		double tmp = f(j, l, w) * pow(p, j) * pow((1 - p), (l - j)) * n;
		total = total + tmp;
		if ((tmp > 0) && (tmp < 1))
			finish = 1;
		j = j + 1;
	}
	return (total);
}

//########################################
//### expErrReads(l,n,p) = expected erroneous reads
//########################################
double expErrReads(int64_t l, int64_t n, double p) {
	return ((1 - pow((1 - p), l)) * n);
}

//########################################
//### expCorrReads(l, n, p) = expected correct reads
//########################################
double expCorrReads(int64_t l, int64_t n, double p) {
	return (pow((1 - p), l) * n);
}

//########################################
//### D(w,l,n,L,p) = expected number of destructible reads, that is, correct reads
//### that are turned wrong because they contain va with v containing errors,
//### a correct but v appearing elsewhere as correct vb
//########################################
double D(int64_t w, int64_t l, int64_t n, int64_t L, double p) {
	// ## v has some errros: (1-(1-p)^w)
	// ## a is correct: (1-p)
	//  ## v appears elsewhere: (1-(1-1/4^w)^L)
	//  ## followed by b<>a: 3/4

	double qw = (1 - pow(1-p, w)) * (1 - p) * (1 - pow(1 - pow(0.25, w), L)) * 3 / 4; //## prob 1 wit gone wrong
	
	//  ## 1 - (1 - qw)^(l-w) = prob 1 read gone wrong

	return ((1 - pow(1 - qw, l - w)) * expCorrReads(l, n, p));
}


int main(int argc, char* argv[]) {

	// Check Arguments
	if (argc != 5) {
		cout << "ERROR: Wrong number of arguments!\n";
		cout << " Usage: ./hitec <inputReads> <correctedReads> <genomeLength> <perBaseErrorRate>" << endl;
		cout << "Please consult the readme file for more information how to run program" << endl;
		exit( EXIT_FAILURE);
	}

	int64_t minCorr, //any letter with minCorr or more occurences in a cluster is considered correct
			maxErr, //any letter with maxErr or less occurences in a cluster is considered erroneous
			minMatch2 = 2, //minimum matched letters required to correct an error
			minClusterSize; // only clusters this size or larger are considered (lower size correspond to erroneous witnesses)

	int64_t minMatch, witLength;  // witLength = length of "witness"; a string preceding a position is a "witness" supporting the character in that position

	double error = atof(argv[4]); // the error is error/100

	ofstream fout;
	fout.open(argv[2]);


//	int64_t maxIteration; 
	int64_t genomeLength = atoi(argv[3]), 
			numberOfReads = 0, // will be set from file 
			readLength,		// will be set from file 
			coverage = 0; // set later as numberOfReads * readLength / genomeLength,

	int64_t size; 

	unsigned char ACTGfrom0123[4] = { 'A', 'C', 'T', 'G' }; // +2%4 gives the reverse complement; works also with ACTGto0123
	int64_t i = 0, j = 0, k = 0, curRead = 0, curPos = 0, curSAPos = 0, 
		allReadsSize = 0, allReadsTotalSize = 0, endCluster = 0;
//	int64_t numberOfErrReads = 0, //number of reads containing one or more errors; computed when reads are created
//			numberOfCorrectedReads = 0, //number of reads which had error before our correction procedure but became error free afterwards
//			numberOfNewErrReads = 0, //number of reads which were initially correct but had errors introduced during the correction procedure
//			rightErrReadCorrected = 0; //number of reads which had error before our correction procedure but became error free or changed to another read that is error free afterwards
	unsigned char curChar;
	clock_t t1, t2, tStart, tEnd;

	
	//******************************************
	//Reading reads from the file   ************
	//******************************************

	unsigned long memoryUsed = 0;

	int64_t numberOfReadsPerSplit,
		numSplit,
		reminder, 
		curSplit,
		boolFlag,
		indCorrReadsFile = 0;

	
	FILE *pfile2 = fopen(argv[1], "r");
	if (!pfile2) {
		cout << "ERROR: Unable to open text file " << endl;
		exit( EXIT_FAILURE);
	}
	fseek(pfile2, 0, SEEK_END); // read file
	int64_t readsSize = ftell(pfile2);
	rewind(pfile2);
	unsigned char *buf = new unsigned char[readsSize];
	memoryUsed += readsSize * sizeof(unsigned char);

	fread((void *) &buf[0], sizeof(unsigned char), readsSize, pfile2); 
	fclose(pfile2);

	i = 0;
	j = 0;

	string fileFormat = "";

	//Determine file format and read Length

	readLength = 0;
	while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
		i++; //blank lines
	
	if (buf[i] == '>')
		fileFormat = "fasta";
	else if (buf[i] == '@')
		fileFormat = "fastq";
	else{
		cerr << "Wrong Input File Format (must be fasta or fastq)" << endl;
		exit(1);
	}
	while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n'))
		i++; //read Header
	while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
		i++; //blank lines
	while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) {
			i++;
			readLength++;
	}

	int64_t readInd = 0;

	unsigned char *** reads;
	int64_t* numReadsPerSplit;

	i =0;

	bool flagDiscRead;
	if(fileFormat == "fasta"){
	// Go through buf array to find number of reads
	while (i < readsSize - 1) {
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		if (buf[i] == '>'){
			while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n'))
				i++; //read Header
		}
		else{
			cerr << "Wrong Header in the Input File Format" << endl; 
			exit(1);
		}
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines

		int len = 0;
		flagDiscRead = false;
		while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) {
			if((buf[i] == 'A')	|| (buf[i] == 'C') || (buf[i] == 'G') || (buf[i] == 'T') ||(buf[i] == 'a') || (buf[i] == 'c') || (buf[i] == 'g')	|| (buf[i] == 't')) {
				i++;
				len++;
			}
			else {
				flagDiscRead = true;
				i++;
			}
		}
		if(!flagDiscRead){
			if (len !=readLength){ cerr << "Error: a read with different length exist." << endl; exit(1);}
			j += len;
		}
	}

	// process reads to keep only reads information *******************
	size = j / readLength; // number of reads containing no N's

	//if readLength > 100, Reads are splited into subset of reads and processed separately. 
	numberOfReadsPerSplit = size; //Each subset of reads has a coverage >=70. Minimum number of reads in each subset is numberOfReadsPerSplit=(70 * genomeLength) / readLength;
	numSplit = 1; // number of split subsets of reads
	//The first numSplit reads are stored in the first array, reads[0], and so on.
	boolFlag = 0; 	//In order to have array with same size, only the reminder first array have an extra element. 
	               // Size of each array is numberOfReadsPerSplit+boolFlag
	reminder = 0;
	
	if(readLength >= 100){
		numberOfReadsPerSplit = (70 * genomeLength) / readLength; 
		numSplit = size / numberOfReadsPerSplit;
		if (numSplit == 0) numSplit = 1;
		numberOfReadsPerSplit = size / numSplit;
		reminder = size % numSplit;
	}


	reads = new unsigned char**[numSplit];
	numReadsPerSplit = new int64_t[numSplit]; // array to keep number of reads in each subset of reads.

	curSplit = 0;
	
	for (int64_t ind = 0; ind < numSplit; ind++){ 
		numReadsPerSplit[ind] = 0;
		reads[ind] = new unsigned char* [numberOfReadsPerSplit+1];
		memoryUsed += (numberOfReadsPerSplit+1) * sizeof(unsigned char);
	}

	for (int64_t ind = 0; ind < numSplit; ind++) {
		for (i = 0; i < numberOfReadsPerSplit+1; i++) {
			reads[ind][i] = new unsigned char[readLength+1];
			reads[ind][i][0] = '\0';
			memoryUsed += readLength * sizeof(unsigned char);
		}
	}

	i = 0;
	j = 0;

	readInd = 0;

	while (i < readsSize - 1) {
		while ((i < readsSize - 1) &&(buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		if (buf[i] == '>'){
			j = 0;
			while ((i < readsSize - 1) &&(buf[i] != '\r') && (buf[i] != '\n'))
				i++; //read Header
		}
		else{
			cerr << "Wrong Header in the Input Fasta File Format." << endl; 
			exit(1);
		}
		while ((i < readsSize - 1) &&(buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		while ((i < readsSize - 1) &&(buf[i] != '\r') && (buf[i] != '\n')) {
			if((buf[i] == 'A')	|| (buf[i] == 'C') || (buf[i] == 'G') || (buf[i] == 'T') ||(buf[i] == 'a') || (buf[i] == 'c') || (buf[i] == 'g')	|| (buf[i] == 't')) {
			reads[curSplit][readInd][j++] = toupper(buf[i++]);
			}
			else
				i++;
		}
	
		//Discard reads containing character not in {A,C,G,T}
		if (j == readLength)
		{
			reads[curSplit][readInd][j] = '\0';
			readInd++;
			numReadsPerSplit[curSplit]++;
			if(curSplit < reminder) boolFlag = 1;
			else  boolFlag = 0;
			if( curSplit < numSplit-1)
				if ( readInd >= numberOfReadsPerSplit + boolFlag) {readInd = 0; curSplit++; }				
		}
	}
	}
	

	else if(fileFormat == "fastq"){

	// Go through buf array to find number of reads
	while (i < readsSize - 1) {
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		if (buf[i] == '@'){
			while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n'))
				i++; //read Header
		}
		else{
			cerr << "Wrong Input File Format (must be fasta or fastq)" << endl; 
			exit(1);
		}
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines

		int len = 0;
		flagDiscRead = false;
		while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) {
			if((buf[i] == 'A')	|| (buf[i] == 'C') || (buf[i] == 'G') || (buf[i] == 'T') ||(buf[i] == 'a') || (buf[i] == 'c') || (buf[i] == 'g')	|| (buf[i] == 't')) {
				i++;
				len++;
			}
			else {
				flagDiscRead = true;
				i++;
			}
		}
		if(!flagDiscRead){
			if (len !=readLength){ cerr << "Error: a read with different length exist." << endl; exit(1);}
			j += len;
		}
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) {
				i++; //line start with +
		}
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) 
				i++; // quality line

	}

	// process reads to keep only reads information *******************
	size = j / readLength; // number of reads

	//if readLength > 100, Reads are splited into subset of reads and processed separately. 
	numberOfReadsPerSplit = size; //Each subset of reads has a coverage 70. Minimum number of reads in each subset is numberOfReadsPerSplit=(70 * genomeLength) / readLength;
	numSplit = 1; // number of splited subset of reads
	//The first numSplit reads are stored in the first array, reads[0], and so on.
	boolFlag = 0; 	//In order to have array with same size, only the reminder first array have an extra element. 
	               // Size of each array is numberOfReadsPerSplit+boolFlag
	reminder = 0;
	
	if(readLength >= 100){
		numberOfReadsPerSplit = (70 * genomeLength) / readLength; 
		numSplit = size / numberOfReadsPerSplit;
		
		if (numSplit == 0) numSplit = 1;
		
		numberOfReadsPerSplit = size / numSplit;
		reminder = size % numSplit;
	}


	reads = new unsigned char**[numSplit];
	numReadsPerSplit = new int64_t[numSplit]; // array to keep number of reads in each subset of reads.

	curSplit = 0;
	

	for (int64_t ind = 0; ind < numSplit; ind++){ 
		numReadsPerSplit[ind] = 0;
		reads[ind] = new unsigned char* [numberOfReadsPerSplit+1];
		memoryUsed += (numberOfReadsPerSplit+2) * sizeof(unsigned char);
	}



	for (int64_t ind = 0; ind < numSplit; ind++) {
		if(ind < reminder) boolFlag = 1;
		else boolFlag = 0;
		for (i = 0; i < numberOfReadsPerSplit+1; i++) {
			reads[ind][i] = new unsigned char[readLength+1];
			reads[ind][i][0] = '\0';
			memoryUsed += readLength * sizeof(unsigned char);
		}
	}


	i = 0;
	j = 0;

	readInd = 0;

	while (i < readsSize - 1) {
		while ((i < readsSize - 1) &&(buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		if (buf[i] == '@'){
			j = 0;
			while ((i < readsSize - 1) &&(buf[i] != '\r') && (buf[i] != '\n'))
				i++; //read Header
		}
		else{
			cerr << "Wrong Input File Format (must be fasta or fastq)" << endl; 
			exit(1);
		}

		while ((i < readsSize - 1) &&(buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines


		while ((i < readsSize - 1) &&(buf[i] != '\r') && (buf[i] != '\n')) {
			if((buf[i] == 'A')	|| (buf[i] == 'C') || (buf[i] == 'G') || (buf[i] == 'T') ||(buf[i] == 'a') || (buf[i] == 'c') || (buf[i] == 'g')	|| (buf[i] == 't')) {
				reads[curSplit][readInd][j++] = toupper(buf[i++]);
			}
			else
				i++;
		}
	
		if (j == readLength)
		{
			reads[curSplit][readInd][j] = '\0';
			readInd++;
			numReadsPerSplit[curSplit]++;
			if(curSplit < reminder) boolFlag = 1;
			else  boolFlag = 0;
			if( curSplit < numSplit-1)
				if ( readInd >= numberOfReadsPerSplit + boolFlag) {readInd = 0; curSplit++; }	
		}

		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) 
				i++; // line start with +
		while ((i < readsSize - 1) && (buf[i] == '\r') || (buf[i] == '\n'))
			i++; //blank lines
		while ((i < readsSize - 1) && (buf[i] != '\r') && (buf[i] != '\n')) 
				i++; //quality line
	}
	}
	else{
		cerr << "Wrong Input File Format, either fasta or fastq is supported" << endl; 
		exit(0);
	}


	delete[] buf;
	memoryUsed -= readsSize * sizeof(unsigned char);



	numberOfReads = numberOfReadsPerSplit;
	coverage = numberOfReads * readLength / genomeLength;

	//########################################
	//### compute wM, wm = witness lengths
	//### wM = min (w | D(w) < 0.0001 * expErrReads)
	//### wm = arg min (U(w) + D(w))
	//########################################
	double p = error/double(100); //assume p is given --- if not then work with p=0.03

	int64_t wM = 1;
	while (D(wM, readLength, numberOfReads, genomeLength, p) > 0.0001 * expErrReads(readLength, numberOfReads, p))
		wM = wM + 1;


	int64_t wm = wM;
	double temp = U(wm, readLength, numberOfReads, p) + D(wm, readLength, numberOfReads, genomeLength, p);
	double temp1 = U(wm - 1, readLength, numberOfReads, p) + D(wm - 1, readLength, numberOfReads, genomeLength, p);
	while (temp1 < temp) {
		wm = wm - 1;
		temp = U(wm, readLength, numberOfReads, p) + D(wm, readLength, numberOfReads, genomeLength, p);
		temp1 = U(wm - 1, readLength, numberOfReads, p) + D(wm - 1, readLength, numberOfReads, genomeLength, p);
	}

	//########################################
	//### compute T (done for each iteration)
	//### T = min (k | Wc(k) > We(k) + 1)
	//########################################

	int64_t w = wM;
	//### here w = wM or ww = m depending on the iteration number

	double qc = (readLength - w) * pow((1 - p), (w + 1)) / (double) genomeLength;
	double qe = (readLength - w) * pow((1 - p), w) * p / (double) (3 * genomeLength);
	int64_t T = 1;
	while (dbinom(T,numberOfReads,qc)*genomeLength <= dbinom(T,numberOfReads,qe)*genomeLength ) {
		T = T + 1;
	}
	T = T + 2;

	minCorr = T;
	maxErr = T-1;
	minClusterSize = T+1;

	
	
	//********************************************************************
	// build string of all reads (and reverse complements), SA, SAinv, LCP
	//********************************************************************

	if(numSplit > 1)
		allReadsTotalSize = (readLength + 1) * (numberOfReads+1) * 2; // total size including z's
	else
		allReadsTotalSize = (readLength + 1) * (numberOfReads) * 2; // total size including z's

	int64_t* SA = new int64_t[allReadsTotalSize];
	memoryUsed += allReadsTotalSize * sizeof(int64_t);


	unsigned char* allReads = new unsigned char[allReadsTotalSize+1]; // all reads and their reverse complements

	for(curSplit = 0; curSplit < numSplit; curSplit++){
		numberOfReads = numReadsPerSplit[curSplit];

	allReadsTotalSize = (readLength + 1) * numberOfReads * 2; // total size including z's

	allReadsSize = readLength * numberOfReads * 2; // size without z's; in SA we stop at allReadsSize-1

	string revRead;
	

	memoryUsed += (allReadsTotalSize + 1) * sizeof(unsigned char);

	curPos = 0; // current position in allReads			// separated by z

	for (i = 0; i < numberOfReads; i++) {
		for (j = 0; j < readLength; j++) 
			allReads[curPos++] = reads[curSplit][i][j]; // the read reads[i]
		
		allReads[curPos++] = 'z';
		for (j = 0; j < readLength; j++)
			allReads[curPos++] = ACTGfrom0123[(ACTGto0123(reads[curSplit][i][readLength
					- 1 - j]) + 2) % 4];
		// the reverse complement of reads[i]
		allReads[curPos++] = 'z';
	}
	allReads[curPos] = '\0';


	//the array "reads[curSplit]" could be deallocated here to save space
	for (curRead = 0; curRead < numberOfReadsPerSplit+1; curRead++) {
		delete[] reads[curSplit][curRead];
		memoryUsed -= readLength * sizeof(unsigned char);
	}

	delete[] reads[curSplit];

	memoryUsed -= numberOfReads * sizeof(unsigned char);

	double used = (double) memoryUsed / 1024;
	used = used / 1024;

	tStart = clock();

	int64_t changedNuc, changedNucTotal = 0 ;
	double changedRatio = 1;
	int64_t iteration = 0;
	int64_t numCluster;
	
	printf("%-20s %20s", "Reads File: ", argv[1]);cout << endl;
	printf("%-20s %20ld", "numSplit: ", numSplit);cout << endl;
	printf("%-20s %17ld", "Total Number Of Reads: ", size);cout << endl;
	printf("%-20s %20ld", "readLength: ", readLength);cout << endl;
	printf("%-20s %20ld", "genomeLength: ", genomeLength); cout << endl;
	printf("%-20s %20ld", "coverage: ", coverage);cout << endl;
	printf("%-19s %20.2f%s", "per-base error: ", error, "%");cout << endl;


	cout << "StartCorrecting...." << endl;

	printf("%-20s %19ld", "correcting curSplit: ", curSplit+1);cout << endl;
	printf("%-20s %16ld", "numberOfReadsPerSplit : ", numReadsPerSplit[curSplit] );cout << endl;


//do procedure till we are able to correct 
	while(changedRatio >= 0.0001 && iteration < 9){
		iteration++;
		changedNuc = 0;

	if(iteration == 1)		 witLength = wm+1;
	else if (iteration == 2) witLength = wM+1;
	else if (iteration == 3) witLength = wM+1;
	else if (iteration == 4) witLength = wm;
	else if (iteration == 5) witLength = wM;
	else if (iteration == 6) witLength = wM;
	else if (iteration == 7) witLength = wm-1;
	else if (iteration == 8) witLength = wM-1;
	else witLength = wM-1;

	cout << "iteration: " << iteration << "..." << endl;
	
	for (i = 0; i < allReadsTotalSize; i++) {
			SA[i] = i;
		}

		if(divsufsort(allReads, SA, allReadsTotalSize) != 0) {
		    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
		    exit(EXIT_FAILURE);
		}
	
		int64_t clusterInd = 0;

		int i1=0, j1=0;
		unsigned char i1c, j1c;

		curSAPos = 0; // curSAPos = current positions in SA
		endCluster = -1; // to start searching the initial cluster from i=0 below

		//*************************************
		//*** find start and end of cluster ***
		//*************************************

		
		double clusterSizeEstimated = 2*genomeLength + 2*numberOfReads*(readLength-witLength+1)*(1-pow((1-p),witLength));
		int64_t currentSize = (int64_t) clusterSizeEstimated /10; //an estimation for size of Cluster	
		

		int64_t *startClusterArr = new int64_t[currentSize] ;
		int64_t *endClusterArr =  new int64_t[currentSize] ;		

		
		i = 0;

		while (i < allReadsSize) {	
			while ((i < allReadsSize) && (SA[i] % (readLength + 1) > readLength	- witLength - 1))
				i++;
			if (i > allReadsSize - 1)
				break; // end of SA -- exit do-while (and terminate program)

			// Double size of Cluster if number of clusters is bigger 
			if(clusterInd >= currentSize){
				
				currentSize = currentSize*2;
				
				// Reallocating startClusterArr and endClusterArr 
			
				int64_t *startClusterArr2 = new int64_t[currentSize] ;
				int64_t *endClusterArr2 =  new int64_t[currentSize] ;
		
				
				for(int64_t ind = 0; ind <= clusterInd-1; ind++){
					startClusterArr2[ind] = startClusterArr[ind];
					endClusterArr2[ind] = endClusterArr[ind];
				}
				
				int64_t *startClusterArrTmp;
				startClusterArrTmp = startClusterArr;
				startClusterArr = startClusterArr2;
				delete[] startClusterArrTmp;
				
				int64_t *endClusterArrTmp;
				endClusterArrTmp = endClusterArr;
				endClusterArr = endClusterArr2;
				delete[] endClusterArrTmp;
				
			}
			
			startClusterArr[clusterInd] = i;
			i++;

			k = 0; i1 = SA[i]; j1 = SA[i-1];
			while ((k <= witLength+1) && (i1+k<allReadsTotalSize) && (j1+k<allReadsTotalSize) && ((i1c = allReads[i1+k]) != 'z') && ((j1c = allReads[j1+k]) != 'z') && (allReads[i1+k] == allReads[j1+k])) 
				k++;

			if ((k < witLength) || (SA[i] % (readLength + 1) == readLength - witLength)){
			 // single position cluster
				endClusterArr[clusterInd] = i - 1;
				if(endClusterArr[clusterInd]-startClusterArr[clusterInd] +1 >= minClusterSize)
					clusterInd++;	
			}
			else {

				do{
				k = 0; i1 = SA[i]; j1 = SA[i-1];
				while ((k <= witLength+1) && (i1+k<allReadsTotalSize) && (j1+k<allReadsTotalSize) && ((i1c = allReads[i1+k]) != 'z') && ((j1c = allReads[j1+k]) != 'z') && (allReads[i1+k] == allReads[j1+k])) 
					k++;

				if ((k > witLength - 1) && (SA[i] % (readLength + 1) != readLength - witLength))
					i++; // continue as long as LCP >= witLength and there is no 'z' right after
				else
					break;
				}while(true);

				endClusterArr[clusterInd] = i - 1;
				if(endClusterArr[clusterInd]-startClusterArr[clusterInd] +1 >= minClusterSize)
					clusterInd++;
			}
		}

		
		numCluster = clusterInd ;
		
#if DEBUG	
		unsigned char* LCPKasai = new unsigned char[allReadsTotalSize];
		t1 = clock();
		LongestCommonPrefixKasai(allReads, SA, LCPKasai, allReadsTotalSize);
		t2 = clock();
		cout << "LCP (Kasai) time: " << (double)(t2-t1)/(double)CLOCKS_PER_SEC << endl;
#endif

		//*******************************
		//print allreads, SA, SAinv, LCP
		//*******************************
#if DEBUG
		cout << "allReads: " << endl;
		cout << " "; for(i=1; i<allReadsTotalSize/10+1; i++) printf("%*d", 10, i); cout << endl;
		for(i=0; i<allReadsTotalSize; i++) cout << i%10; cout << endl;
		cout << allReads << endl;

		cout << "allReadsTotalSize = " << allReadsTotalSize << endl;
		cout << "allReadsSize = " << allReadsSize << endl;

		// compting SA inverse
		int64_t* SAinv = new int64_t[allReadsTotalSize];
		for(i=0; i<allReadsTotalSize; i++)
		SAinv[SA[i]] = i;

		cout << "   i  SA SAi LCP  " << endl;
		for (i=0; i<allReadsTotalSize; i++) {
			printf ("%*d%*d%*d%*d ", 4, i, 4, SA[i], 4, SAinv[i], 4, LCP[i]);
			for(j=0; j<LCP[i]; j++) cout << allReads[SA[i]+j];
			if ((j=SA[i]+LCP[i]) < allReadsTotalSize) cout << (char)tolower(allReads[j]);
			cout << endl;
		}
#endif

		//*************************************************
		// Correcting
		//*************************************************

//		int pos;

		int64_t freqChar[4];
		int64_t SASubCluster[4];

		for (i = 0; i < 4; i++)
			SASubCluster[i] = -1;

		// a cluster in SA is either a single position i = startCluster = endCluster with
		// text[SA[i]..SA[i]+witLength] contains no 'z' and LCP[i] < witLength
		// or an interval [startCluster..endCluster] such that
		// LCP[i]>= witLength for all i in startsCluster+1 .. endCluster
		// text[SA[i]+witLength] != 'z', for all i = startCluster..endCluster
		// LCP[startsCluster] < witLength
		// LCP[endCluster+1] < witLength or = witLength but there is a 'z' right after

//		int corrInt;

		//************************************
		//*** start processing in SA order ***
		//************************************

		curSAPos = 0; // curSAPos = current positions in SA
		endCluster = -1; // to start searching the initial cluster from i=0 below

//		int changed = 0;
//		int LCP_i;

		clusterInd = 0;

		// *************************
		// **** process cluster ****
		//************** reset findings for this cluster; ACTG=0123 -- if not reset, then they remember everything

		for(clusterInd = 0; clusterInd < numCluster ; clusterInd++){
{			
			for (i = 0; i < 4; i++)
				freqChar[i] = 0;

			for (curSAPos = startClusterArr[clusterInd]; curSAPos < endClusterArr[clusterInd] + 1; curSAPos++) {
				curChar = allReads[SA[curSAPos] + witLength];

				if(curChar != 'A' && curChar != 'C' && curChar != 'G' && curChar != 'T' )
					cerr << "wrong char: " << curChar << endl;

				if (freqChar[ACTGto0123(curChar)] == 0){
					SASubCluster[ACTGto0123(curChar)] = curSAPos; //save start SA position for subCluster of curChar
				}

				//save frequency of each nucleotides in a cluster, size of each sub-cluster
				freqChar[ACTGto0123(curChar)]++;
			}

			int i1, j1, maxLCP, minLCP, corrChar = 0; //, readsErr;
			int64_t pos1, pos2;

//			bool flag = false;

			//set the parameter based on the clusters

//			int corrList[4], indCorrList = 0, noLetters = 0;
			int  indCorrList = 0, noLetters = 0;

			for (i = 0; i < 4; i++) {
				if (freqChar[i] > 0)
					noLetters++;
				if (freqChar[i] >= minCorr)
					indCorrList++;
			}

			if(noLetters > 2) minMatch = minMatch2;
			else      minMatch = 0;

			if (indCorrList > 0) {//changing procedure is possible
				for (i = 0; i < 4; i++) {
					if (freqChar[i] <= maxErr) { // it is believed to be an error
						for (i1 = 0; i1 < freqChar[i]; i1++) { //for all positions in sub cluster i which is belived to be an error

							if(SASubCluster[i] + i1 >= allReadsTotalSize)
							{
								cerr << "exceed allReadsTotalSize!! " << clusterInd << endl;
								exit(1);
							}
							pos1 = SA[SASubCluster[i] + i1] + witLength + 1; //position of first nucleotide after current char which is suspected to be error

							//checking matching string after correct possition and error position
							maxLCP = 0;
							for (j = 0; j < 4; j++) { //for each correct subCluster j
								if (freqChar[j] >= minCorr) { // j is believed to be correct
									for (j1 = 0; j1 < freqChar[j]; j1++) { //for all reads in subcluster j
										pos2 = SA[SASubCluster[j] + j1] + witLength + 1; //position of first nucleotide after current char which is believed to be correct
										//find match length of wrong nucleotide and correct one
										minLCP = 0;
										unsigned char i1c, j1c;
										int errorTolerance = 1;

										if (allReads[pos1] == 'z') {
											minLCP = 1;
											corrChar = j;
											break;
										}

										//let one misMatch, if errorTolerance is negative means there are more than one misMatch
										while ((errorTolerance >= 0) && (pos1 + minLCP < allReadsTotalSize)	&& (pos2 + minLCP < allReadsTotalSize) 
											&& ((i1c = allReads[pos1 + minLCP]) != 'z')	&& ((j1c = allReads[pos2 + minLCP]) != 'z')) {
												//	while ((pos1+minLCP< allReadsTotalSize) && (pos1+minLCP< allReadsTotalSize) && ((i1c = allReads[pos1+minLCP]) != 'z') && ((j1c = allReads[pos1+minLCP]) != 'z') && (allReads[pos1+minLCP] == allReads[pos2+minLCP]))
												if (allReads[pos1 + minLCP] != allReads[pos2 + minLCP]){
													errorTolerance--;
												}
												if (errorTolerance >= 0)
													minLCP++;

										}
										if (maxLCP < minLCP) {
											maxLCP = minLCP;
											corrChar = j;
										}
									}
								}
							}
							if (SA[SASubCluster[i] + i1] + witLength < allReadsTotalSize) {
								int64_t posInSA = SA[SASubCluster[i] + i1] + witLength;
								//change if there is a match
								if (maxLCP >= min(minMatch, max((int64_t) 1, readLength - (posInSA % (readLength + 1)) - 1))) {
									//curPos = SA[SASubCluster[i] + i1] % ((readLength + 1) * 2) + witLength;									
									int64_t revComPos = 0;
									curPos = SA[SASubCluster[i] + i1] % ((readLength + 1) * 2) + witLength; //position in read
									if (curPos < readLength){ //valid position in forward read
										revComPos = SA[SASubCluster[i] + i1] + (readLength - 1 - curPos) * 2 + witLength + 2;
									}
									else if (curPos < 2 * readLength + 1){
										revComPos = SA[SASubCluster[i] + i1] - ((curPos	- witLength	- readLength - 1) * 2 + witLength + 2);
									}
									if(posInSA < allReadsTotalSize ){
										allReads[posInSA] = ACTGfrom0123[corrChar]; //correct error
									}
									changedNuc++;
									changedNucTotal++;
									if(revComPos < allReadsTotalSize){
										allReads[revComPos]	= ACTGfrom0123[(corrChar + 2) % 4]; //correct error in reverse complement
									}
									changedNuc++;
									changedNucTotal++;
								}
							}
						}
					}
				}
			}
			}
		} // done with all positions in SA (up to allReadsSize)

		t2 = clock();
		
//		double changedRatioByNow = (double) changedNucTotal / (double) allReadsSize;
		changedRatio = (double) changedNuc / (double) allReadsSize;

		zaparr(startClusterArr);
		zaparr(endClusterArr);

	}


	cout << "correction is done:" << endl;

	//writing corrected reads int64_to a file
	unsigned char* corrReads = new unsigned char[allReadsSize];
	int64_t indLast = 0;

	char strTmp[1000];

	for (i = 0; i < numReadsPerSplit[curSplit]; i++) {
		curPos = i * ((readLength + 1) * 2);
		//write the header
		corrReads[indLast++] = (unsigned char) '>';
		sprintf(strTmp,"%ld",indCorrReadsFile);
		indCorrReadsFile++;
		for(size_t j = 0; j < strlen(strTmp); j++) 
			corrReads[indLast++] = strTmp[j];
		corrReads[indLast++] = (unsigned char) '\n';	

		for (int64_t j = 0; j < readLength; j++){				
			unsigned char ch = allReads[j + curPos];
			corrReads[indLast++] = ch;
		}
		corrReads[indLast++] = (unsigned char) '\n';
	}
	corrReads[indLast-1] = '\0';
	fout << corrReads << endl;
	}

	tEnd = clock();

	double time = (double) (tEnd - tStart) / (double) CLOCKS_PER_SEC;

	t1 = clock();


	fout.close();

	t2 = clock();
	time = (double) (t2 - t1) / (double) CLOCKS_PER_SEC;


	struct timeval tim;
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	tim = ru.ru_utime;
	double t = (double) tim.tv_sec + (double) tim.tv_usec / 1000000.0;
	tim = ru.ru_stime;
	t += (double) tim.tv_sec + (double) tim.tv_usec / 1000000.0;

	printf("%-20s %20s", "Corrected Reads has been copied in ", argv[2]);cout << endl;


	return 0;
}

