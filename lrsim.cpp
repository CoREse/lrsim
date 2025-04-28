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
double RegionFloatVariance=0.0;
double RegionFloatVarianceRatio=0.0;
int BlockSize=100;
double BlockFloatVariance=0.05;
double BlockFloatVarianceRatio=1.0;
double ReadFloatVariance=0.05;
double AccurateSequencingThreshold=0.03;
double ReadFloatVarianceRatio=1.0;
double MaxError=0.4;
double MaxErrorRatio=4.0;
double MinError=0.01;
double MinErrorRatio=0.1;
double LengthFloatStd=0.015;
bool NoLengthFloat=false;
int FixedReadLength=0;
bool TailingN=false;
double Linkage=1.0;
string RunHash="";
int ThreadN;
int Ins=22, Del=49, Sub=0;
pthread_mutex_t m_Output;

htsThreadPool p = {NULL, 0};

string simRead(const string & RefSeq, double ErrorRate, int Start, int Length, int Min, int Max, const vector<double> & RefRegionErrorFloat, int QueueN, mt19937 & Generator)
{
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
    int End=RefSeq.length();
    uniform_int_distribution<int> BaseDist(0,3);
    uniform_int_distribution<int> TypeDist(0,All-1);
    uniform_real_distribution<double> ErrorDist;
    int BlockI=-1;
    normal_distribution<double> BlockFloatDist(0,BlockFloatVariance);
    double BlockFloat=0;
    int LastError=0;//0: no error, 1: insertion, 2: deletion, 3: sub
    int LinkageInsT, LinkageDelT;
    if (Linkage!=0)
    {
        LinkageInsT=100;
        LinkageDelT=(int)(((double)LinkageInsT)/Linkage);
    }
    int CurrentInsT,CurrentDelT;
    uniform_int_distribution<int> LinkageTypeDist(0,LinkageInsT+LinkageDelT+Sub-1);
    int TypeI;
    while (i < End)
    {
        if (BlockI!=(i-Start)/BlockSize)
        {
            BlockI=(i-Start)/BlockSize;
            BlockFloat=BlockFloatDist(Generator);
        }
        if (ErrorDist(Generator)<max(min(MaxError,ErrorRate+RefRegionErrorFloat[i/RegionSize]+BlockFloat),MinError))
        {
            if ((LastError==1 || LastError==2) &&Linkage!=0)
            {
                CurrentDelT=LinkageDelT;
                CurrentInsT=LinkageInsT;
                TypeI=LinkageTypeDist(Generator);
            }
            else
            {
                CurrentDelT=DelT;
                CurrentInsT=InsT;
                TypeI=TypeDist(Generator);
            }
            if (TypeI<CurrentInsT)
            {
                Read+=Bases[BaseDist(Generator)];
                LastError=1;
                continue;
            }
            else if (TypeI<CurrentDelT) LastError=2;
            else
            {
                Read+=Bases[BaseDist(Generator)];
                LastError=3;
            }
        }
        else
        {
            Read+=RefSeq[i];
            LastError=0;
        }
        if (Read.length()>=Length) break;
        ++i;
    }
    if (TailingN && Read.length()<Length)
    {
        for (int i=Read.length();i<Length;++i) Read+='N';
    }
    if (Read.length()<Min) Read="";
    if (Read.length()>Max) Read=Read.substr(0,Max);
    return Read;
}

void outputReads(int QueueN, int SN, vector<string> & Reads, string & Quals)
{
    int i=SN-Reads.size();
    for (int j=0;j<Reads.size();++j)
    {
        fprintf(stdout, "@Simulated_%s_%d_%d\n%s\n+\n%s\n",RunHash.c_str(), QueueN,i,Reads[j].c_str(),Quals.substr(0,Reads[j].length()).c_str());
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
    int QueueN;
    int DSize;
    unsigned * DistLengths;
    unsigned * DistFormers;
};

void runSim(const vector<string> & Ref, vector<string> & RefNames, unsigned long long TotalLength, const vector<unsigned long long> & RefAccL, double ErrorRate, unsigned long long NumberOfBases, const vector<vector<double>> &RegionErrorFloat, int Seed, int QueueN, int DSize=0, unsigned * DistLengths=nullptr, unsigned * DistFormers=nullptr)
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
    Generator.seed(Seed+QueueN);
    normal_distribution<double> LenDist(Mean,Variance);
    uniform_int_distribution<unsigned> LenDDist(0, DSize==0?0:DSize-1);
    uniform_int_distribution<unsigned long long> PosDist(0,TotalLength-1);
    while (BN < NumberOfBases)
    {
        unsigned long long Pos=PosDist(Generator);
        int RefIndex=getRefInd(RefAccL,Pos);
        int Length=0;
        int Former=0;
        if (FixedReadLength!=0)
        {
            Length=FixedReadLength;
            if (!NoLengthFloat)
            {
                normal_distribution<double> LVariance(1,LengthFloatStd);
                double LV=LVariance(Generator);
                Length=int (((double)Length)*LV);
            }
        }
        else
        {
            if (DSize==0) Length=LenDist(Generator);
            else
            {
                unsigned LengthI=LenDDist(Generator);
                Length=DistLengths[LengthI];
                Former=DistFormers[LengthI];
                uniform_int_distribution<unsigned> LDist(Former, Length);
                Length=LDist(Generator);
                if (!NoLengthFloat)
                {
                    normal_distribution<double> LVariance(1,LengthFloatStd);
                    double LV=LVariance(Generator);
                    Length=int (((double)Length)*LV);
                }
            }
        }
        string Read=simRead(Ref[RefIndex], ErrorRate, Pos, Length, Min, Max, RegionErrorFloat[RefIndex], QueueN, Generator);
        if (Read.length()==0) continue;
        BN+=Read.length();
        SN+=1;
        Reads.push_back(Read);
        if (ThreadN==1)
        {
            if (Reads.size()>1000)
            {
                outputReads(QueueN,SN,Reads,Quals);
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
                outputReads(QueueN,SN,Reads,Quals);
                pthread_mutex_unlock(&m_Output);
                Reads.clear();
                // }
            }
        }
    }
    pthread_mutex_lock(&m_Output);
    outputReads(QueueN,SN,Reads,Quals);
    pthread_mutex_unlock(&m_Output);
}

void * runSimHandler(void * Args)
{
    RunSimArgs * pArgs=(RunSimArgs *) Args;
    runSim(*pArgs->pRef, * pArgs->pRefNames, pArgs->TotalLength, * pArgs->pRefAccL, pArgs->ErrorRate, pArgs->NumberOfBases, *pArgs->pRegionErrorFloat, pArgs->ThisSeed, pArgs->QueueN, pArgs->DSize, pArgs-> DistLengths, pArgs->DistFormers);
    delete pArgs;
    return NULL;
}

void sim(vector<const char *> &RefFileNames, double ErrorRate, double Depth, string BasesString, string ModelFileName, string DistFileName, int Seed, int ThreadN)
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
    if (ModelFileName!="")
    {
        fprintf(stderr, "Reading model file %s...\n",ModelFileName.c_str());
        FILE* df=fopen(ModelFileName.c_str(),"rb");
        int LRSMVersion=0;
        fread(&LRSMVersion,4,1,df);
        if (LRSMVersion==0)
        {
            fread(&Ins,4,1,df);
            fread(&Del,4,1,df);
            fread(&DSize,4,1,df);
            DistLengths=(unsigned *)malloc(sizeof(unsigned)*DSize);
            DistFormers=(unsigned *)malloc(sizeof(unsigned)*DSize);
            fread(DistLengths,4,DSize,df);
            fread(DistFormers,4,DSize,df);
            fprintf(stderr,"Model %s (v%d) read, error ratio: %d:%d:%d, length dist size: %u.\n", ModelFileName.c_str(), LRSMVersion, Ins, Del, Sub, DSize, DistLengths[0]);
        }
        else if (LRSMVersion==1)
        {
            fread(&Ins,4,1,df);
            fread(&Del,4,1,df);
            fread(&Linkage,8,1,df);
            fread(&DSize,4,1,df);
            DistLengths=(unsigned *)malloc(sizeof(unsigned)*DSize);
            DistFormers=(unsigned *)malloc(sizeof(unsigned)*DSize);
            fread(DistLengths,4,DSize,df);
            fread(DistFormers,4,DSize,df);
            fprintf(stderr,"Model %s (v%d) read, error ratio: %d:%d:%d, linkage: %lf, length dist size: %u.\n", ModelFileName.c_str(), LRSMVersion, Ins, Del, Sub, Linkage, DSize, DistLengths[0]);
        }
    }
    else if (DistFileName!="")
    {
        fprintf(stderr, "Reading dist file %s...\n",DistFileName.c_str());
        FILE* df=fopen(DistFileName.c_str(),"rb");
        fread(&DSize,4,1,df);
        DistLengths=(unsigned *)malloc(sizeof(unsigned)*DSize);
        DistFormers=(unsigned *)malloc(sizeof(unsigned)*DSize);
        fread(DistLengths,4,DSize,df);
        fread(DistFormers,4,DSize,df);
        fprintf(stderr,"Dist file read, size: %u.\n", DSize, DistLengths[0]);
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
        runSim(Ref, RefNames, TotalLength, RefAccL, ErrorRate, NBases, RegionErrorFloat, Seed, 0, DSize, DistLengths, DistFormers);
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
            RunSimArgs *A=new RunSimArgs{&Ref, &RefNames, TotalLength, &RefAccL, ErrorRate, EachBases, &RegionErrorFloat, Seed, i, DSize, DistLengths, DistFormers};
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
    const char * const Version="v0.7";
	for (int i=1;i<argc;++i) RunString+=string(" ")+argv[i];
	size_t Hash=hash<string>()(RunString);
	stringstream ss;
	ss<<std::hex<<Hash;
	ss>>RunHash;
    ThreadN=1;
    double Depth=1;
    int Seed=0;
    string ModelFileName="";
    string DistFileName="";
    string BasesString="";
    double ErrorRate=0.15;
    string IERatio="22:49:0";
    bool ShowVersion=false;
	OptHelper OH=OptHelper("lrsim [Options] fa [fa2] [fa3] ...");
    OH.addOpt('b', "bases", 1, "Number", "Number of bases to be simulated, will omit -d. Format: 500, 500K, 500M, 5.1G (1K=1000).",'S',&(BasesString));
    OH.addOpt('d', "depth", 1, "Number", "Sequencing depth.",'F',&(Depth));
    OH.addOpt('s', "seed", 1, "Number", "Random seed.",'i',&(Seed));
    OH.addOpt('t', "threads", 1, "Number", "Number of threads (different threadn will result in different results).",'i',&(ThreadN));
    OH.addOpt('f', "distfile", 1, "FileName", "Distribution binary file generated by extractDistData.py. (Deprecated)", 'S', &DistFileName);
    OH.addOpt('m', "modelfile", 1, "FileName", "Model binary file generated by samtools stats real.bam | extractModel.py.", 'S', &ModelFileName);
    OH.addOpt('e', "error", 1, "Number", "Error rate.",'F',&(ErrorRate));
    OH.addOpt(0, "rvarianceratio", 1, "Number", "Multiplies errorrate to form read error rate float variance.",'F',&(ReadFloatVarianceRatio));
    OH.addOpt(0, "bfvarianceratio", 1, "Number", "Multiplies errorrate to form read block error rate float variance.",'F',&(BlockFloatVarianceRatio));
    OH.addOpt(0, "blocksize", 1, "Number", "Region size for error rate floating.",'i',&(BlockSize));
    OH.addOpt(0, "rfvarianceratio", 1, "Number", "Multiplies errorrate to form region error rate float variance.",'F',&(RegionFloatVarianceRatio));
    OH.addOpt(0, "regionsize", 1, "Number", "Region size for error rate floating.",'i',&(RegionSize));
    OH.addOpt(0, "eratio", 1, "INS:DEL:SUB", "Indel error ratio. Substitutions could be seen as Ins+Del. -m will override the INS (with 100) and DEL but not the SUB",'S',&(IERatio));
    OH.addOpt(0, "lengthfloatratio", 1, "Number", "Ratio std for read length float.",'F',&(LengthFloatStd));
    OH.addOpt(0, "fixedreadlength", 1, "Length", "Specify a fixed read length, if not 0, will omit length distribution but will still be influenced by length float, along with --nolengthfloat --tailingn to get absolutely fixed read length.",'i',&(FixedReadLength));
    OH.addOpt(0, "nolengthfloat", 0, "", "Do not randomly float the read length.",'b',&(NoLengthFloat));
    OH.addOpt(0, "tailingn", 0, "", "Padding the read with N at the tail if the read reaches the end of the chromosome (useful if you want absolute fixed length reads)",'b',&(TailingN));
    OH.addOpt(0, "version", 0, "", "Show version and exit.",'b',&(ShowVersion));
    OH.getOpts(argc,argv);
    if (ShowVersion) {fprintf(stderr, "lrsim %s\nBy CRE\n", Version);return 0;}
    if (OH.Args.size()<=0) {OH.showhelp(); return 1;}
    fprintf(stderr, "lrsim %s\n",Version);
    fprintf(stderr, "Generating simulated long reads for");

    //get IERatio
    for (int i=0;i<3;++i)
    {
        int End=IERatio.find(":");
        if (End==string::npos)
        {
            if (i==0) Ins=atoi(IERatio.c_str());
            else if (i==1) Del=atoi(IERatio.c_str());
            else Sub=atoi(IERatio.c_str());
            break;
        }
        if (i==0) Ins=atoi(IERatio.substr(0, End).c_str());
        else if (i==1) Del=atoi(IERatio.substr(0, End).c_str());
        else Sub=atoi(IERatio.substr(0, End).c_str());
        IERatio.erase(0, End + 1);
    }

    for (int i=0;i<OH.Args.size();++i) fprintf(stderr, " %s", OH.Args[i]);
    if (BasesString!="") fprintf(stderr, " of total %s of bases.", BasesString.c_str());
    else fprintf(stderr, " at depth %.2lf.", Depth);
    MinError=ErrorRate*MinErrorRatio;
    MaxError=ErrorRate*MaxErrorRatio;
    BlockFloatVariance=BlockFloatVarianceRatio*ErrorRate;
    if (ErrorRate<AccurateSequencingThreshold) ReadFloatVarianceRatio=2.0;//Accurate sequencing correction
    ReadFloatVariance=ReadFloatVarianceRatio*ErrorRate;
    RegionFloatVariance=RegionFloatVarianceRatio*ErrorRate;
    fprintf(stderr," Error rate: %lf. Seed: %d.\n", ErrorRate, Seed);
    sim(OH.Args, ErrorRate, Depth, BasesString, ModelFileName, DistFileName, Seed, ThreadN);
    fprintf(stderr,"Done reads simulation for");
    for (int i=0;i<OH.Args.size();++i) fprintf(stderr, " %s", OH.Args[i]);
    if (BasesString!="") fprintf(stderr, " of total %s of bases.", BasesString.c_str());
    else fprintf(stderr, " at depth %.2lf.", Depth);
    return 0;
}
