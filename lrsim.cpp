#include "optutils/OptHelper.h"
#include <string>
#include <vector>
#include "htslib/htslib/faidx.h"
#include <assert.h>
#include <random>
#include <stdio.h>
using namespace std;

const char Bases[4]={'A', 'T', 'G', 'C'};

string simRead(string & RefSeq, double ErrorRate, int Start, int Length, int Min, int Max, mt19937 & Generator)
{
    int Ins=22, Del=49, Sub=0;
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
    while (i< End)
    {
        if (ErrorDist(Generator)<ErrorRate)
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
        fprintf(stdout, "@Simulated_%d\n%s\n+\n%s\n",i,Reads[j].c_str(),Quals.substr(0,Reads[j].length()).c_str());
        ++i;
    }
}

int getRefInd(vector<unsigned long long> & RefAccL, unsigned long long &Position)
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

void runSim(vector<string> & Ref, vector<string> & RefNames, unsigned long long TotalLength, vector<unsigned long long> & RefAccL, unsigned long long NumberOfBases, mt19937 & Generator, int DSize=0, unsigned * DistLengths=nullptr, unsigned * DistFormers=nullptr)
{
    fprintf(stderr,"Generating reads...Ref length: %d, TotalBases=%llu\n",TotalLength, NumberOfBases);
    double Mean=8000;
    double Variance=7000;
    int Min=1;
    int Max=100000;
    string Quals="";
    for (int i=0;i<Max;++i) Quals+='!';
    unsigned long long BN=0;
    int SN=0;
    vector<string> Reads;
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
        string Read=simRead(Ref[RefIndex], 0.15, Pos, Length, Min, Max, Generator);
        if (Read.length()==0) continue;
        BN+=Read.length();
        SN+=1;
        Reads.push_back(Read);
        if (Reads.size()>1000)
        {
            outputReads(SN,Reads,Quals);
            Reads.clear();
        }
    }
    outputReads(SN,Reads,Quals);
}

void sim(vector<const char *> &RefFileNames, double Depth, string DistFileName, mt19937 & Generator, int ThreadN)
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
    runSim(Ref, RefNames, TotalLength, RefAccL, (unsigned long long) (double(TotalLength)*Depth), Generator, DSize, DistLengths, DistFormers);
    if (DistLengths!=nullptr) free(DistLengths);
    if (DistFormers!=nullptr) free(DistFormers);
}

int main(int argc, const char* argv[])
{
    int ThreadN=1;
    double Depth=1;
    int Seed=0;
    string DistFileName="";
	OptHelper OH=OptHelper("lrsim [Options] fa [fa2] [fa3] ...");
    OH.addOpt('t', "threads", 1, "Number", "Number of threads.",'i',&(ThreadN));
    OH.addOpt('d', "depth", 1, "Number", "Sequencing depth.",'F',&(Depth));
    OH.addOpt('s', "seed", 1, "Number", "Random seed.",'i',&(Seed));
    OH.addOpt('f', "distfile", 1, "FileName", "Distribution binary file generated by extractDistData.py.", 'S', &DistFileName);
    OH.getOpts(argc,argv);
    mt19937 Generator;
    Generator.seed(Seed);
    if (OH.Args.size()<=0) {OH.showhelp(); return 1;}
    fprintf(stderr, "Generating simulated long reads for %s... at depth %.2lf. Seed: %d.\n", OH.Args[0], Depth, Seed);
    sim(OH.Args, Depth, DistFileName, Generator, ThreadN);
    return 0;
}