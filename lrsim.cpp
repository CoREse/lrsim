#include "optutils/OptHelper.h"
#include <string>
#include <vector>
#include "htslib/htslib/faidx.h"
#include <assert.h>
#include <random>
#include <stdio.h>
#include <sstream>
#include "htslib/htslib/thread_pool.h"
using namespace std;

const char Bases[4]={'A', 'T', 'G', 'C'};
int RegionSize=300;
double RegionFloatVariance=0.05;
int BlockSize=100;
double BlockFloatVariance=0.05;
double ReadFloatVariance=0.05;
double MaximumError=0.4;
double MinError=0.01;
string RunHash="";
int ThreadN;
pthread_mutex_t m_Output;

htsThreadPool p = {NULL, 0};

string simRead(const string & RefSeq, double ErrorRate, int Start, int Length, int Min, int Max, const vector<double> & RefRegionErrorFloat, mt19937 & Generator)
{
    int Ins=22, Del=49, Sub=0;
    if (ReadFloatVariance!=0.0)
    {
        normal_distribution<double> ReadErrorFloatDist(0,ReadFloatVariance);
        ErrorRate+=ReadErrorFloatDist(Generator);
    }
    int All=Ins+Del+Sub;
    int InsT=Ins;
    int DelT=Ins+Del;
    string Read="";
    int i=Start;
    int End=Start+Length;
    if (End>RefSeq.length()) End=RefSeq.length();
    uniform_int_distribution<int> BaseDist(0,3);
    uniform_int_distribution<int> TypeDist(0,All-1);
    uniform_real_distribution<double> ErrorDist;
    int BlockI=-1;
    normal_distribution<double> BlockFloatDist(0,BlockFloatVariance);
    double BlockFloat=0;
    while (i < End)
    {
        if (BlockI!=(i-Start)/BlockSize)
        {
            BlockI=(i-Start)/BlockSize;
            BlockFloat=BlockFloatDist(Generator);
        }
        if (ErrorDist(Generator)<max(min(MaximumError,ErrorRate+RefRegionErrorFloat[i/RegionSize]+BlockFloat),MinError))
        {
            int TypeI=TypeDist(Generator);
            if (TypeI<InsT)
            {
                Read+=Bases[BaseDist(Generator)];
                continue;
            }
            else if (TypeI<DelT) ;
            else
            {
                Read+=Bases[BaseDist(Generator)];
            }
        }
        else
        {
            Read+=RefSeq[i];
        }
        ++i;
    }
    if (Read.length()<Min) Read="";
    if (Read.length()>Max) Read=Read.substr(0,Max);
    return Read;
}

void outputReads(int SN, vector<string> & Reads, string & Quals)
{
    int i=SN-Reads.size();
    for (int j=0;j<Reads.size();++j)
    {
        fprintf(stdout, "@Simulated_%s_%d\n%s\n+\n%s\n",RunHash.c_str(),i,Reads[j].c_str(),Quals.substr(0,Reads[j].length()).c_str());
        ++i;
    }
}

int getRefInd(const vector<unsigned long long> & RefAccL, unsigned long long &Position)
{
    for (int i=0;i<RefAccL.size();++i)
    {
        if (Position<RefAccL[i])
        {
            Position-=(i==0?0:RefAccL[i-1]);
            return i;
        }
    }
}

struct RunSimArgs
{
    vector<string> * pRef;
    vector<string> * pRefNames;
    unsigned long long TotalLength;
    vector<unsigned long long> * pRefAccL;
    double ErrorRate;
    unsigned long long NumberOfBases;
    vector<vector<double>> *pRegionErrorFloat;
    int ThisSeed;
    int DSize;
    unsigned * DistLengths;
    unsigned * DistFormers;
};

void runSim(const vector<string> & Ref, vector<string> & RefNames, unsigned long long TotalLength, const vector<unsigned long long> & RefAccL, double ErrorRate, unsigned long long NumberOfBases, const vector<vector<double>> &RegionErrorFloat, int ThisSeed, int DSize=0, unsigned * DistLengths=nullptr, unsigned * DistFormers=nullptr)
{
    double Mean=8000;
    double Variance=7000;
    int Min=1;
    int Max=1000000;
    string Quals="";
    for (int i=0;i<Max;++i) Quals+='!'-int(10.0*log10(ErrorRate));
    unsigned long long BN=0;
    int SN=0;
    vector<string> Reads;
    mt19937 Generator;
    Generator.seed(ThisSeed);
    normal_distribution<double> LenDist(Mean,Variance);
    uniform_int_distribution<unsigned> LenDDist(0, DSize==0?0:DSize-1);
    uniform_int_distribution<unsigned long long> PosDist(0,TotalLength-1);
    while (BN < NumberOfBases)
    {
        unsigned long long Pos=PosDist(Generator);
        int RefIndex=getRefInd(RefAccL,Pos);
        int Length=0;
        int Former=0;
        if (DSize==0) Length=LenDist(Generator);
        else
        {
            unsigned LengthI=LenDDist(Generator);
            Length=DistLengths[LengthI];
            Former=DistFormers[LengthI];
            uniform_int_distribution<unsigned> LDist(Former, Length);
            Length=LDist(Generator);
            normal_distribution<double> LVariance(1,0.015);
            double LV=LVariance(Generator);
            Length=int (((double)Length)*LV);
        }
        string Read=simRead(Ref[RefIndex], ErrorRate, Pos, Length, Min, Max, RegionErrorFloat[RefIndex], Generator);
        if (Read.length()==0) continue;
        BN+=Read.length();
        SN+=1;
        Reads.push_back(Read);
        if (ThreadN==1)
        {
            if (Reads.size()>1000)
            {
                outputReads(SN,Reads,Quals);
                Reads.clear();
            }
        }
        else
        {
            if (Reads.size()>1000)
            {
                // if (pthread_mutex_trylock(&m_Output)==0)
                // {
                pthread_mutex_lock(&m_Output);
                outputReads(SN,Reads,Quals);
                pthread_mutex_unlock(&m_Output);
                Reads.clear();
                // }
            }
        }
    }
    pthread_mutex_lock(&m_Output);
    outputReads(SN,Reads,Quals);
    pthread_mutex_unlock(&m_Output);
}

void * runSimHandler(void * Args)
{
    RunSimArgs * pArgs=(RunSimArgs *) Args;
    runSim(*pArgs->pRef, * pArgs->pRefNames, pArgs->TotalLength, * pArgs->pRefAccL, pArgs->ErrorRate, pArgs->NumberOfBases, *pArgs->pRegionErrorFloat, pArgs->ThisSeed, pArgs->DSize, pArgs-> DistLengths, pArgs->DistFormers);
    delete pArgs;
    return NULL;
}

void sim(vector<const char *> &RefFileNames, double ErrorRate, double Depth, string BasesString, string DistFileName, int Seed, int ThreadN)
{
    vector<string> RefNames;
    vector<string> Ref;
    vector<unsigned long long> RefAccL;
    unsigned long long TotalLength=0;
    for (int i=0;i<RefFileNames.size();++i)
    {
        faidx_t * RefFile=fai_load(RefFileNames[i]);
        int NSeq=faidx_nseq(RefFile);
        char * Seq;
        for (int i=0;i<(NSeq);++i)
        {
            const char * Name=faidx_iseq(RefFile,i);
            RefNames.push_back(Name);
            int Length;
            Seq =fai_fetch(RefFile,Name,&Length);
            assert(Length>=0);
            Ref.push_back(Seq);
            free(Seq);
            TotalLength+=Length;
            RefAccL.push_back(TotalLength);
        }
        fai_destroy(RefFile);
    }
    vector<vector<double>> RegionErrorFloat;
    normal_distribution<double> RegionFloatDist(0,RegionFloatVariance);
    mt19937 Generator;
    Generator.seed(Seed);
    for (int i=0;i<Ref.size();++i)
    {
        RegionErrorFloat.push_back(vector<double>());
        for (int j=0;j<Ref[i].size()/RegionSize+1;++j)
        {
            RegionErrorFloat[i].push_back(RegionFloatDist(Generator));
        }
    }
    unsigned * DistLengths=nullptr;
    unsigned * DistFormers=nullptr;
    unsigned DSize=0;
    if (DistFileName!="")
    {
        fprintf(stderr, "Reading dist file %s...\n",DistFileName.c_str());
        FILE* df=fopen(DistFileName.c_str(),"rb");
        fread(&DSize,4,1,df);
        DistLengths=(unsigned *)malloc(sizeof(unsigned)*DSize);
        DistFormers=(unsigned *)malloc(sizeof(unsigned)*DSize);
        fread(DistLengths,4,DSize,df);
        fread(DistFormers,4,DSize,df);
        fprintf(stderr,"Dist file read, size:%u\n", DSize, DistLengths[0]);
    }
    unsigned long long NBases=0;
    if (BasesString=="")
    {
        NBases=(unsigned long long) (double(TotalLength)*Depth);
    }
    else
    {
        if ((BasesString[BasesString.length()-1]&char(0b11011111))=='M') NBases=(unsigned long long)(((double)atof(BasesString.substr(0,BasesString.length()-1).c_str()))*1000000.0);
        else if ((BasesString[BasesString.length()-1]&char(0b11011111))=='K') NBases=(unsigned long long)(((double)atof(BasesString.substr(0,BasesString.length()-1).c_str()))*1000.0);
        else if ((BasesString[BasesString.length()-1]&char(0b11011111))=='G') NBases=(unsigned long long)(((double)atof(BasesString.substr(0,BasesString.length()-1).c_str()))*1000000000.0);
        else NBases=(unsigned long long)((double)atof(BasesString.substr(0,BasesString.length()).c_str()));
    }
    fprintf(stderr,"Generating reads...Ref length: %llu, TotalBases=%llu\n",TotalLength, NBases);
    if (ThreadN==1)
    {
        runSim(Ref, RefNames, TotalLength, RefAccL, ErrorRate, NBases, RegionErrorFloat, Seed, DSize, DistLengths, DistFormers);
    }
    else
    {
        if (!(p.pool = hts_tpool_init(ThreadN))) {
			fprintf(stderr,"Error creating thread pool\n");
            exit(1);
		}
		p.qsize=ThreadN*2;
        hts_tpool_process *RunProcess=hts_tpool_process_init(p.pool,p.qsize,1);
		pthread_mutex_init(&m_Output,NULL);
        unsigned long long EachBases=NBases/p.qsize;
        vector<mt19937> Generators;
        for (int i=0;i<p.qsize;++i)
        {
            // Generators.push_back(mt19937());
            // Generators[Generators.size()-1].seed(Seed+i);
            // RunSimArgs *A=new RunSimArgs{&Ref, &RefNames, TotalLength, &RefAccL, ErrorRate, EachBases, &RegionErrorFloat, &(Generators[Generators.size()-1]), DSize, DistLengths, DistFormers};
            RunSimArgs *A=new RunSimArgs{&Ref, &RefNames, TotalLength, &RefAccL, ErrorRate, EachBases, &RegionErrorFloat, Seed+i, DSize, DistLengths, DistFormers};
            hts_tpool_dispatch(p.pool,RunProcess,runSimHandler,(void *)A);
        }
        hts_tpool_process_flush(RunProcess);
        hts_tpool_process_destroy(RunProcess);
    }
    if (DistLengths!=nullptr) free(DistLengths);
    if (DistFormers!=nullptr) free(DistFormers);
}

int main(int argc, const char* argv[])
{
    string RunString="";
    const char * const Version="v0.3.1";
	for (int i=1;i<argc;++i) RunString+=string(" ")+argv[i];
	size_t Hash=hash<string>()(RunString);
	stringstream ss;
	ss<<std::hex<<Hash;
	ss>>RunHash;
    ThreadN=1;
    double Depth=1;
    int Seed=0;
    string DistFileName="";
    string BasesString="";
    double ErrorRate=0.15;
	OptHelper OH=OptHelper("lrsim [Options] fa [fa2] [fa3] ...");
    OH.addOpt('t', "threads", 1, "Number", "Number of threads (different threadn will result in different results).",'i',&(ThreadN));
    OH.addOpt('d', "depth", 1, "Number", "Sequencing depth.",'F',&(Depth));
    OH.addOpt('s', "seed", 1, "Number", "Random seed.",'i',&(Seed));
    OH.addOpt('f', "distfile", 1, "FileName", "Distribution binary file generated by extractDistData.py.", 'S', &DistFileName);
    OH.addOpt('b', "bases", 1, "Number", "Number of bases to be simulated, will omit -d. Format: 500, 500K, 500M, 5.1G (1K=1000).",'S',&(BasesString));
    OH.addOpt('e', "error", 1, "Number", "Error rate.",'F',&(ErrorRate));
    OH.addOpt('r', "rfvariance", 1, "Number", "Region error rate float variance.",'F',&(RegionFloatVariance));
    OH.addOpt(0, "regionsize", 1, "Number", "Region size for error rate floating.",'i',&(RegionSize));
    OH.getOpts(argc,argv);
    if (OH.Args.size()<=0) {OH.showhelp(); return 1;}
    fprintf(stderr, "lrsim %s\n",Version);
    fprintf(stderr, "Generating simulated long reads for");
    for (int i=0;i<OH.Args.size();++i) fprintf(stderr, " %s", OH.Args[i]);
    if (BasesString!="") fprintf(stderr, " of total %s of bases.", BasesString.c_str());
    else fprintf(stderr, " at depth %.2lf.", Depth);
    fprintf(stderr," Error rate: %lf. Seed: %d.\n", ErrorRate, Seed);
    sim(OH.Args, ErrorRate, Depth, BasesString, DistFileName, Seed, ThreadN);
    fprintf(stderr,"Done reads simulation for");
    for (int i=0;i<OH.Args.size();++i) fprintf(stderr, " %s", OH.Args[i]);
    if (BasesString!="") fprintf(stderr, " of total %s of bases.", BasesString.c_str());
    else fprintf(stderr, " at depth %.2lf.", Depth);
    return 0;
}