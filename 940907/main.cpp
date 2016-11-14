#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <set>
#include <math.h>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <limits.h>
#include <fstream>
#include <vector>
#include <queue>
#include <utility> 
#include <time.h>
#include <string>

#define cerdl cerr << endl
#define CE cerr << endl << 

#define _FILE_OFFSET_BITS 64
#include <stdlib.h>
#include <unistd.h>
#include <sys/resource.h>
#define errExit(msg)do { perror(msg); exit(EXIT_FAILURE);} while (0) 
#define MAX_CHILD 15
using namespace std;

int NUM_Data;
int NUM_Feature;
int NUM_TEST;
double ADD;
double SEMILAR;
char* name_store;
clock_t difft;
double **train;
double **test;
double *output_label;
double *C_test;
double *C;
set<double>*ux;
set<double>uc;
double *maxValF;
int **adj; // this is adjacent matrix
pair<double , double> **IPAXCF;
vector<int>*pa; // parent of Nodes
bool swVStructure;
bool rmv_;
double lambda;
int mxLevel;

struct Informations
{
	double IFjPafj;
	double IFjGPafj;
	double IFjGPafjSPa;
	double IFjCSGPa;
};

struct ent
{
    ent(int n , double info)
    {
        xj = n;
        IXiXj = info;
    }
    int xj;
    double IXiXj;
};

struct Node
{
    vector<Node*>childs;
    bool vis;
    int szBrother;
    int szChild;
    int fnum;
    double IXC;
    double H;
    int level ;
    double IFiCSpa;
    double IFiPa;
    Node(int fum , double IC,double Hf)
    {
		H = Hf;
        vis = true;
        fnum = fum;
        IXC = IC;   
    }
};

double* intersect(int numf , double valf , double* isect , int szIsect)
{
	if(!szIsect)
	{
		for( int i = 0 ; i < NUM_Data ; i++)
		{
			if(train[i][numf] == valf)
				isect[i] = 1;
			else
				isect[i] = 0;
		}
	}
	else
	{
		for( int i = 0 ; i < NUM_Data ; i++)
		{
			if(train[i][numf] == valf && isect[i] == 1) /// dddddddddd add && isect[i] 
				isect[i] = 1;
			else
				isect[i] = 0;
		}
	}
	return isect;
}

double summ(double *isect , int sz)
{
	double num = 0;
	for (int i  =0 ; i < sz ; i++)
	{
		if(isect[i])
			num++;
	}
	return num;
}

double learning_()
{
	double  numTrue = 0;
	FILE *of = fopen(name_store,"a+");
	FILE *ftime;
	FILE *facc; 
	ofstream outAcc("res.txt");
	ftime = fopen("Time.txt","a");
	facc = fopen("acc.txt","a");
	for( int i = 0 ; i < NUM_TEST ; i++ )
	{
		long double cmax = -1;
		long double pcmax = -1;
		for(set<double>::iterator it = uc.begin(); it != uc.end() ; it++)
	    {			
			long double pc = 1;
			double *isect = new double[NUM_Data];
			//intersect
			memset(isect,0,sizeof(double)*NUM_Data);
			for (int j = 0 ; j < NUM_Feature ; j++)
			{
				long double pcc = 1;
			
				int szIsect = 0;
				bool sw = false;
				for (vector<int>::iterator iter =  pa[j].begin() ; iter != pa[j].end() ; iter++)
				{
					sw = true;
					if((*iter) == NUM_Feature)
					{
						isect = intersect(NUM_Feature , (*it) , isect , szIsect);
						szIsect = NUM_Data;
					}
					else
					{
						isect = intersect((*iter) , test[i][(*iter)] , isect , szIsect);
						szIsect = NUM_Data;
					}
				}
				double pCondition = summ(isect,szIsect);
				isect = intersect(j , test[i][j] , isect , szIsect);
				szIsect = NUM_Data;
				double pjoint = summ(isect,szIsect);
				if(sw)
				{
					pcc = (pjoint / (pCondition+1));
				}
				else
					pcc = pjoint / NUM_Data;
				pc = pc*pcc;
			}
			if( pc > pcmax)
			{
				pcmax = pc;
				cmax = (*it);
			}
		}
		output_label[i] = cmax;
		outAcc << cmax << endl;
		if(cmax == C_test[i])
			numTrue++;
	}
	cout << "NUM true " << numTrue << endl;
	fprintf(of,"accuracy is %f \n",numTrue  / NUM_TEST);
	double acc__ = (double)numTrue  / (double)NUM_TEST;
	fprintf(facc, "%s %lf\n" , name_store , numTrue  / NUM_TEST);
	fprintf(of,"total Time is  %f",(((float)difft)/CLOCKS_PER_SEC));
	fprintf(ftime, "%s %lf\n" ,name_store, (((float)difft)/CLOCKS_PER_SEC));
	fclose(ftime);
	fclose(facc);
	fclose(of);
	return acc__;
}


double IxixjC(int i , int j)
{
    double sm = 0;
    for(set<double>::iterator itc = uc.begin(); itc != uc.end() ; itc++)
        for(set<double>::iterator it = ux[i].begin(); it != ux[i].end() ; it++)
            for(set<double>::iterator it2 = ux[j].begin(); it2 != ux[j].end() ; it2++)
            {
                double cnt = 0;
                double cnti = 0;
                double cntj = 0;
                double numC = 0;
                for( int f = 0 ; f < NUM_Data ; f++)
                {
                    if(C[f] == (*itc))numC++;
                    if(train[f][i] == (*it) && train[f][j] == (*it2) && C[f] == (*itc))
                        cnt++;
                    if(train[f][i] == (*it) && C[f] == (*itc))
                        cnti++;
                    if(train[f][j] == (*it2) && C[f] == (*itc))
                        cntj++;
                }
                cnt /= numC;
                cnti /= numC;
                cntj /= numC;
                if(cnt == 0);
                else
                {
                    sm += cnt*log((cnt)/(cnti*cntj));
                }
            }
    return sm;
}

Informations IXjXiPa(int pa , int j,int grand_father)
{
	Informations result;
    double sm = 0;
    double smIFiPafi = 0;
    double smIFjGpa = 0;
	for(set<double>::iterator it = ux[pa].begin(); it != ux[pa].end() ; it++)
		for(set<double>::iterator itGPa = ux[grand_father].begin(); itGPa != ux[grand_father].end() ; itGPa++)
            for(set<double>::iterator it2 = ux[j].begin(); it2 != ux[j].end() ; it2++)
            {
				double rnd =  rand() / (double)INT_MAX;
				if(rnd < 0.4) continue;
                double cnt = 0;
                double cntGPa = 0;
                double cntj = 0;
                double numPa = 0;
                double numFi = 0;
                double cntjGpa = 0;
                double numGpa =  0;
                
                for( int f = 0 ; f < NUM_Data ; f++)
                {
                    if(train[f][pa] == (*it))numPa++;
                    if(train[f][j] == (*it2))numFi++;
                    if(train[f][grand_father] == (*itGPa))numGpa++;
                    if(train[f][pa] == (*it) && train[f][j] == (*it2) && train[f][grand_father] == (*itGPa))
                        cnt++;
                    if(train[f][pa] == (*it) && train[f][j] == (*it2))
                        cntj++;
                    if(train[f][grand_father] == (*itGPa) && train[f][j] == (*it2))
                        cntjGpa++;
                    if(train[f][pa] == (*it) && train[f][grand_father] == (*itGPa))
                        cntGPa++;
                }
                double numConPFiFi = cntj;
                cnt /= numPa;
                cntGPa /= numPa;
                cntj /= numPa;
                if(cnt == 0);
                else
                {
                    sm += cnt*log((cnt)/(cntGPa*cntj));
                    smIFiPafi += (numConPFiFi / NUM_Data)*log((numConPFiFi / NUM_Data)/((numFi/NUM_Data)*(numPa/NUM_Data)));
                    smIFjGpa += (cntjGpa / NUM_Data)*log((cntjGpa / NUM_Data)/((numFi/NUM_Data)*(numGpa/NUM_Data)));
                }
            }
            result.IFjGPafj = smIFjGpa;
            result.IFjGPafjSPa = sm;
            result.IFjPafj = smIFiPafi;
    return result;
}

Informations IXjCPa(int pa , int j)
{
	Informations res;
    double sm = 0;
    double smIFiPafi = 0;
    double smFjC = 0;
	for(set<double>::iterator it = ux[pa].begin(); it != ux[pa].end() ; it++)
	{
		double rnd = rand() / (float)INT_MAX;
		if(rnd < 0.4) continue;
		for(set<double>::iterator itc = uc.begin(); itc != uc.end() ; itc++)
            for(set<double>::iterator it2 = ux[j].begin(); it2 != ux[j].end() ; it2++)
            {
                double cnt = 0;
                double cntC = 0;
                double cntj = 0;
                double numPa = 0;
                double numFi = 0;
                double numFjC = 0;
                double numC = 0;
                for( int f = 0 ; f < NUM_Data ; f++)
                {
					if(C[f] == (*itc))numC++;
                    if(train[f][pa] == (*it))numPa++;
                    if(train[f][j] == (*it2))numFi++;
                    if(train[f][pa] == (*it) && train[f][j] == (*it2) && C[f] == (*itc))
                        cnt++;
                    if(train[f][pa] == (*it) && train[f][j] == (*it2))
                        cntj++;
                    if(train[f][pa] == (*it) && C[f] == (*itc))
                        cntC++;
                    if(train[f][j] == (*it2) && C[f] == (*itc))
						numFjC++;
                }
                
                if(cnt == 0);
                else
                {
					smIFiPafi += (cntj / NUM_Data)*log((cntj/NUM_Data)/((numFi/NUM_Data)*(numPa/NUM_Data)));
					smFjC += (numFjC / NUM_Data)*log((numFjC/NUM_Data)/((numFi/NUM_Data)*(numC/NUM_Data)));
                    sm += (cnt/(numPa))*log((cnt/(numPa))/((cntC/numPa)*(cntj/numPa)));
                }
            }
	}
    res.IFjGPafjSPa = sm;
    res.IFjPafj = smIFiPafi;
    res.IFjGPafj = smFjC;
    
    return res;
}

double cxic(int j)
{
    double sm = 0;
    for(set<double>::iterator it = uc.begin(); it != uc.end() ; it++)
    {
        for(set<double>::iterator it2 = ux[j].begin(); it2 != ux[j].end() ; it2++)
        {
            double cnt  = 0;
            double cnti = 0;
            double cntj = 0;
            for( int f = 0 ; f < NUM_Data ; f++)
            {
                if((C[f] == (*it)) && ( train[f][j] == (*it2)))
                    cnt++;
                if((C[f] == (*it)))
                    cnti++;
                if((train[f][j] == (*it2)))
                    cntj++;
            }
            cnt  = cnt / NUM_Data;
            cnti /= NUM_Data;
            cntj /= NUM_Data;
            if(cnt == 0);
            else
            {
                sm += cnt*log((cnt)/(cnti*cntj));
            }
        }
    }
    return sm;
}

double HF(int f)
{
	double out = 0;
	for(set<double>::iterator it = ux[f].begin() ; it != ux[f].end() ; it++)
	{
		double num =0;
		for(int i = 0 ; i < NUM_Data ;i++)
		{
			if((*it) == train[i][f])
				num++;
		}
		out += -1 * (num/NUM_Data)*(log((num/NUM_Data))/log(2));
	}
	return out / ux[f].size();
}

pair<vector<Node*>::iterator,Informations> search(Node* &parent,Node * nd , int level , Node *grand_father_node)
{
	
	vector<Node*> &vc = parent->childs;
	rmv_ = false;
    double mx = -10000;
    vector<Node*>::iterator itMx;
    
    Informations tem;
    Informations mxInfo;
    tem.IFjGPafj = 0;
    tem.IFjPafj = 0;
    tem.IFjGPafjSPa = 0;
	mxInfo.IFjGPafj = 0;
    mxInfo.IFjPafj = 0;
    mxInfo.IFjGPafjSPa = 0;

    bool szLimetedB = false;
    bool swHaveChild = false;
    if (vc.size() > MAX_CHILD)
		szLimetedB = true;
    for(vector<Node*>::iterator it = vc.begin() ; it != vc.end() ; it++)
    {
		double rnd =  rand() / (double)INT_MAX;
		if(rnd < 0.2) continue;
		swHaveChild = true;
		tem = IXjCPa((*it)->fnum , nd->fnum);
		(*it)->IFiCSpa = tem.IFjGPafjSPa;
		(*it)->IFiPa = tem.IFjPafj;	

        if(tem.IFjPafj - tem.IFjGPafjSPa > mx)
        {
            mx = tem.IFjGPafjSPa;
            itMx = it;
            mxInfo = tem;
        }
    }
    double hfi = nd->H;
	if( mxInfo.IFjPafj - mxInfo.IFjGPafjSPa < 0 && mxInfo.IFjGPafj - hfi < 0)
	{
		rmv_ = true;
		cerr << "errr rmv_ " << endl;
	}
    if(level < 2 && swHaveChild && (szLimetedB || mxInfo.IFjPafj - mxInfo.IFjGPafjSPa > lambda * mxInfo.IFjGPafj - hfi ))
        return make_pair(itMx,mxInfo);
    else if(level > 1 && swHaveChild && (szLimetedB || parent-> IFiPa - parent->IFiCSpa < (*itMx)->IFiPa - (*itMx)->IFiCSpa ))   
		return make_pair(itMx,mxInfo);
    else
        return make_pair(vc.end(),mxInfo);
}

struct Graph
{
    Node *first;
    int sz;
    Graph()
    {
        first = new Node(NUM_Feature,999999,0);
        sz = 10000;
    }

    bool insert_(int fnum , double IC , double H)
    {
        int level = 0;
        Node *newNode = new Node(fnum,IC,H);
        Node *cur = first;
        cur-> level = level;
        vector<Node*>::iterator itCur;
        pair<vector<Node*>::iterator,Informations> outSearch;
        if(IC > ADD)
        {
            Node* prev = first;
            outSearch = search(cur,newNode,1,first);
            itCur = outSearch.first;
            if(rmv_)
				return false;
            if(itCur != cur->childs.end() && !swVStructure)
            {
                cur = (*itCur);
                cur-> level = level;
            }
            if(itCur != cur->childs.end() && swVStructure)
            {
                cur = (*itCur);
                cur-> level = level;
                newNode->childs.push_back(cur);
                first->childs.push_back(newNode);
                swVStructure = false;
                return true;
            }
            if(itCur == cur->childs.end())
            {
                first->childs.push_back(newNode);
                return true;
            }
            else
            {
                level++;
                cur = (*itCur);
                cur-> level = level;
                while(cur != NULL && (level < mxLevel || prev->IFiPa - prev->IFiCSpa < cur->IFiPa - cur->IFiCSpa))
                {
                    if(level > mxLevel)
                    {
                        cerr << "Rejected" << endl;
                        cur->childs.clear();
                        return false;
                    }
                    outSearch = search(cur,newNode,level,prev);
					itCur = outSearch.first;
					if(rmv_)
						return false;
                    
                    if(itCur != cur->childs.end() && !swVStructure)
                    {
                        level++;
                        prev = cur;
                        cur = (*itCur);
                    }
                    else if(itCur != cur->childs.end() && swVStructure)
                    {
                        cur = (*itCur);
						newNode->childs.push_back(cur);
						first->childs.push_back(newNode);
						swVStructure = false;
						return true;
                    }
                    else
                        break;
                }
                if(itCur == cur->childs.end())
                {
                    cur->childs.push_back(newNode);
                    return true;
                }
                else if(!swVStructure && prev->IFiPa - prev->IFiCSpa > cur->IFiPa - cur->IFiCSpa && cur->IXC < newNode->IXC)
                {
                    newNode->childs.push_back(cur);
                    prev->childs.erase(itCur);
                    prev->childs.push_back(newNode);
                    return true;
                }
                else if (swVStructure)
                {
					cur = (*itCur);
					newNode->childs. push_back(cur);
					first->childs.push_back(newNode);
					swVStructure = false;
					return true;
				}
            }
        }
        return true;
    }
};

int operator <(ent x1 , ent x2)
{
    return x1.IXiXj < x2.IXiXj;
}

double Ixixj(int i , int j)
{
    double sm = 0;
    for(set<double>::iterator it = ux[i].begin(); it != ux[i].end() ; it++)
        for(set<double>::iterator it2 = ux[j].begin(); it2 != ux[j].end() ; it2++)
        {
            double cnt = 0;
            double cnti = 0;
            double cntj = 0;
            for( int f = 0 ; f < NUM_Data ; f++)
            {
                if(train[f][i] == (*it) && train[f][j] == (*it2))
                    cnt++;
                if(train[f][i] == (*it))
                    cnti++;
                if(train[f][j] == (*it2))
                    cntj++;
            }
            cnt /= NUM_Data;
            cnti /= NUM_Data;
            cntj /= NUM_Data;
            if(cnt == 0);
            else
            {
                sm += cnt*log((cnt)/(cnti*cntj));
            }
        }
    return sm;
}

void findParent(Graph g_)
{
	Node *cur = g_.first;
	queue<Node*>qu;
	qu.push(cur);
	cur->level = 0;
	while(!qu.empty())
	{
		cur = qu.front();
		qu.pop();
		if(cur -> level > mxLevel )
			continue;
		for(vector<Node*>::iterator iter = cur->childs.begin() ; iter != cur->childs.end() ; iter++)
		{
			Node* ch = (*iter);
			pa[ch->fnum].push_back(cur->fnum);
			ch->level = cur->level + 1;
			qu.push(ch);
		}
	}
}


void bfs_(Graph g_)
{
	FILE *of;
	of = fopen(name_store,"a+");
    queue<Node*>qug;
    Node* cur = g_.first;
    qug.push(cur);
    int numSelected = 0 ;
    for(int i = 0 ; i < NUM_Feature+1 ; i++)
        for(int j = 0 ; j < NUM_Feature +1 ; j++)
            adj[i][j] = 0;
    cur->level = 0;
    fprintf(of,"\n");
    vector<double>out_Info;
    while(!qug.empty())
    {
        cur = qug.front();
        qug.pop();
        fprintf(of,"%d ",cur->fnum);
        out_Info.push_back(cur->IXC);
        numSelected++;
        if(cur->level > mxLevel )continue;
        for(int i = 0 ; i < (int)cur->childs.size();i++)
        {
			cur->childs[i]->level = cur->level +1;
            qug.push(cur->childs[i]);
            adj[cur->fnum][cur->childs[i]->fnum]=1;
        }
    }
    fprintf(of,"\n");
    for(int i = 0 ; i < (int)out_Info.size();i++)
    {
		fprintf(of,"%f ",out_Info[i]);
	}
    fprintf(of,"\n%d\n",numSelected);
    fclose(of);
    freopen("adj.txt","w",stdout);
    //~ for ( int i = 0 ; i < NUM_Feature+1 ; i++)
    //~ {
        //~ for ( int j = 0 ; j < NUM_Feature+1 ; j++)
            //~ cout << adj[i][j]<<',';
        //~ cout << endl;
    //~ }
    cerr << "bfs end " << endl;
}
int swapdata(int firstInx,int lastInx)
{
	for(int rt = firstInx,tinx=0 ; tinx < NUM_TEST ; tinx++,rt++)
	{
		for(int c = 0 ; c < NUM_Feature+1 ; c++)
		{
			if(c == NUM_Feature)
			{
				double t = C[rt];
				C[rt] = C_test[tinx];
				C_test[tinx] = t;
			}
			double t = train[rt][c];
			train[rt][c]  = test[tinx][c];
			test[tinx][c] = t;
		}
	}
	return 0;
}
int main(int argc , char **argv)
{
	struct rlimit old, new_;
	struct rlimit *newp;
	pid_t pid;
	cerr << argc << endl;

	pid = 0;

	ADD = 0.02;
	SEMILAR = 0.02;
	new_.rlim_cur = (long long)1024*1024;
	new_.rlim_max = (long long)1024*1024*1024*1024;
	newp = &new_;

	if (prlimit(pid, RLIMIT_CPU, newp, &old) == -1)
		errExit("prlimit-1");
	printf("Previous limits: soft=%lld; hard=%lld\n",

	(long long) old.rlim_cur, (long long) old.rlim_max);

	/* Retrieve and display new CPU time limit */

	if (prlimit(pid, RLIMIT_CPU, NULL, &old) == -1)
		errExit("prlimit-2");
	printf("New limits: soft=%lld; hard=%lld\n",
	(long long) old.rlim_cur, (long long) old.rlim_max);

	//~ exit(EXIT_FAILURE);
	
	swVStructure = false;
	
	lambda = 1;
	mxLevel = 5;
	name_store = argv[2];
	ofstream of(name_store);
	of.close();
	ifstream input(argv[1]);
	string s;
	int numLine = 0;
	int numCol = -1;
	bool sw = true;
	while(getline(input,s))
	{
		if(sw)
		{
			stringstream ss;
			ss << s;
			string t;
			while(ss >> t)numCol++;
			sw = false;
		}
		numLine++;
	}
	input.close();
	CE 1;
	NUM_Data = (numLine*90) / 100+1;
	NUM_TEST = (numLine*10) / 100;
	int step = NUM_TEST;
	NUM_Feature = numCol;
	
	adj = new int*[NUM_Feature+1];
	for (int i = 0 ; i  < NUM_Feature+1 ; i++)
		adj[i] = new int[NUM_Feature+1];
	
	output_label = new double [NUM_TEST];
	CE 2;
	C_test = new double[NUM_TEST];
	CE 3;
	C = new double [NUM_Data];
	CE 4;
	ux = new set<double>[NUM_Feature + 1];
	CE 5;
	maxValF = new double [NUM_Feature];
	CE 6;
	pa = new vector<int>[NUM_Feature+1];
	CE 7;
    IPAXCF  = new pair<double , double>*[NUM_Feature+1];
    train = new double*[NUM_Data];
    for(int i = 0 ; i < NUM_Data ; i++)
        train[i] = new double[NUM_Feature+1];
    CE 8;
    test = new double*[NUM_TEST];
	for( int i = 0 ; i < NUM_TEST ; i++)
	{
		test[i] = new double[NUM_Feature+1];
	}
    CE 9 ;
    memset(maxValF,INT_MIN,NUM_Feature);
	int fold = 0 ;
	input.open(argv[1]);
	string line;
	int tInsert=0;
	int vInsert=0;
	cerr << "NTrain " << NUM_Data << " " << "Ntest " << NUM_TEST << " NF " << NUM_Feature << endl;
	while(getline(input,line))
	{
		stringstream ss;
		ss << line;
		bool insertInTrain = true;
		for(int f=0; f < NUM_Feature + 1 ; f++)
		{
			double fd;
			ss >> fd;
			ux[f].insert(fd);
			double rnd = (double)rand() / (double)INT_MAX;
			cerr << rnd << endl;
			//~ if((rnd < 0.8 || vInsert >= NUM_TEST) && tInsert < NUM_Data)
			if(tInsert < NUM_Data)
			{
				train[tInsert][f] = fd;
				if(f == NUM_Feature){
					C[tInsert] = fd;
					uc.insert(fd);
				}
			}
			else if(vInsert < NUM_TEST)
			{
				test[vInsert][f] =fd;
				if(f == NUM_Feature)
					C_test[tInsert]= fd;
				insertInTrain = false;
			}
			else
				cerr << "error read file " << endl;
		}
		if(insertInTrain)
			tInsert++;
		else
			vInsert++;
		ss.clear();
	}
	input.close();
	bool firstFold = false;
	double acc = 0;
	int firstInx = 0;
	int lastInx = firstInx+step;
    do
    {
		Graph g_;
		cerr << fold ;cerdl;
		cerr << firstInx << endl;
		cerr << lastInx << endl;
		cerr << NUM_Data << endl;
		cerr << tInsert << endl;
		cerr << NUM_TEST << endl;
		cerr << vInsert << endl;
		
		if(!firstFold)
		{
			swapdata(firstInx,lastInx);
			firstInx = lastInx;
			lastInx += step;
		}
		firstFold = false;
		vector<ent> ve;
		clock_t stTime = clock();
		
		for (int f = 0 ; f < NUM_Feature ; f++)
		{
			double icx = cxic(f);
			double H = HF(f);
			g_.insert_(f,icx,H);
		}
		clock_t edTime = clock();
		difft = edTime -  stTime;
		
		cerr << "find Parent" << endl;
		findParent(g_);
		cerr << "learning_ " << endl;
		acc += learning_();
		cerr << "bfs_ " << endl;
		bfs_(g_);
		fold++;
	}while(fold < 9);
	delete []train;
	delete []test;
	delete [] IPAXCF;
	delete []ux;
	delete []output_label;
	delete []pa;
	cerr << acc/9 << endl;
    return 0;
}
