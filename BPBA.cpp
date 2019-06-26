// BPBA.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ekk_c_api.h"
#include "afxtempl.h"
#include <malloc.h>                //包含分配内存的函数或操作符号
#include <time.h>
#include <set>                     //不知是否有用
#include "Ob.h"
#include "MakeStack.h"
#include "Furnace.h"
#include <stdlib.h>
#include <assert.h>
#include <ilcplex/ilocplex.h>
#include "CLibSVM.h"
#include<windows.h>
//#include <vld.h>
//#include <vldapi.h>
using namespace std;               //使用名字空间
#define C_MAX 300
#define F_MAX 100
#define MAX 100000
#define EPS 0.0001
#define Permit_Gap   0.003                  //permit integer gap
#define XN 4
typedef IloArray<IloBoolVarArray>     IloBoolVarArray22;
typedef IloArray<IloNumVarArray>     IloNumVarArray22;

int Col_numb,Fur_numb;       //板卷和炉子的数量
int NHB,NHS,HHB,HHS;
float C1[C_MAX+1][F_MAX+1];  //板卷与炉子间的关联系数
float C2[C_MAX+1][C_MAX+1];  //板卷与板卷间的关联系数
int   arc_matrix[C_MAX+1][C_MAX+1]; //板卷与板卷之间的分支
int   fur_coil[F_MAX+1][C_MAX+1]; //路子与板卷之间的分支
int median[C_MAX+1];  //是否在方案中的矩阵
int Ykk[C_MAX][F_MAX]; //是否是用心板卷的矩阵
double TestData[10000][60];//机器学习训练数据
int RowsNumOfTestDataSVM;//测试数据txt的数据行数
int Sub_problem_Num;
float totalcost;

COb coil_inf[C_MAX+1];
struct coils
{
	CString number;               //板卷序号
	float h;                  //板卷高度
	float w;                  //板卷重量
	float f;                  //板卷未使用的惩罚费用
	float ww;                 //板卷的综合权重
	int flag;                 //板卷是否被使用
	float PRI;
} coil[C_MAX+1];

CFurnace batch[F_MAX+1];
struct Knap
{
	float cost;
	float weight;
	int ingot[301];
	int ingot_num;
	int MIni_temperature;
	CString ingot_type;
} pp;//[200000];

//struct furnace *fur_pool;      //负费用列池
//struct furnace *furlist;       //列表用于存储当前解得列
CList<CFurnace,CFurnace> pit_pool;  //负费用列池
CList<CFurnace,CFurnace> pitlist;   //初始解

int con_m,var_n;              //限制主问题的函数及列数
double *a;                    //限制主问题的矩阵
double dual_var[C_MAX+F_MAX+1+10000];     //对偶变量
int       NEW_COL_NUM;             //all column number be generated
int old_col_num;

double   run_time;           //运行时间
clock_t   begin_time;         //开始时间
clock_t   finish_time;        //结束时间


double   sub_time;           //运行时间
clock_t   sub_begin_time;         //开始时间
clock_t   sub_finish_time;        //结束时间

double   final_gap;          //对偶间隙
double   lowerb;             //下界
double   uperb;              //上界即初始解
bool change_flag=false;
bool is_optimal;
double last_lower_bound=0;

struct AddCut
{
	int coil[11][2];
	int right;
	int number;
};

struct n 
{
	double lower_bound;
	int arc_maximal[C_MAX][C_MAX];
	int fur_coil[F_MAX][C_MAX];
//	CFurnace pit_knap[F_MAX];
	COb ingot_knap[C_MAX];
	int level;
	int node_type;
	int tree;
	int *Col;
	int column_num;
	int solu_numb;
	double solution[C_MAX+2];//每个节点上
	int soluindex[C_MAX+2];
	int model_modify[C_MAX];
	int s;
	int r;
	int fur_select;
	int coil_select;
	int solu_type;
	bool is_change;
	int ingot_num_knap;      //求解背包问题时的钢锭
	int median[C_MAX+1];   //是否是中心板卷
	int is_median;          //分支所选出来的中心
	int NodeNum;
	int ParentNum;
	int Ykk[C_MAX][F_MAX];

	//IloEnv env;
	//IloModel mod;
	//IloNumVarArray Lamada;//线性规划变量
	//IloRangeArray range;
	//IloNumArray C;
	//IloObjective cost;

	CList<AddCut,AddCut> CutList;
	struct n *next;
};

CList<AddCut,AddCut> CurrentCutList;//当前加入的割

//EKKVector  *opt_col ;         //columns correspond to  optimal solution结构体，后面有定义
//EKKModel   * Linear_Mode;     //linear mode 
//EKKContext * env;             //osl enviorment
CList<CFurnace,CFurnace> fur_list;

struct n *nodelist;//?定义了两个struct n类型的结构体
struct n *result;
struct n *nownode;
double LBinBoot;
int nodenum;
int       *ColInParent;
//bool change_flag=false;
COb ingot_knap[C_MAX];   //求解背包时用到的钢锭数组
COb branch_coil[C_MAX];   //分支时用于辨别相似的板卷
int ingot_num_knap;             //背包问题中的钢锭数
CFurnace pit_knap[F_MAX]; 
CList<double,double> soln_of_each_node;     //每个节点中对应线性规划的最优解
int col_num;                       //每个节点中对应线性规划的列数
double opt_soln_of_eachnode;       //每个节点中对应线性规划的最优目标函数
int numincol=0;
float feasible_value=10000;
int model_modify[C_MAX];
int median_from;     //记录中心枚举从哪开始
int median_to;     //记录中心枚举从哪结束
double LowerBound[1000];
int i_LowerBound = 0;

int valid_number;  //在主问题中加入的优先级有效不等式个数
COb valid_maste[C_MAX]; //在主问题加入的优先级有效不等式
int valid_total_num;  //优先级约束总数
int number_upperbound;  //装入板卷个数上界
int NumOfValidInequality;

CList<int,int> FNlist;

float **dkk;   //d[k][k][j];
float **eij;   //e[i][j];


FILE  *A;
int total_column_number=0;
/*----------------------------------------函数声明--------------------------------------*/

void init_data(int dream);                                      //读数据
void coil_data(FILE *read_file,FILE *fp,int dream);              //数据初始化
void coil_coil_coefficient();                          //计算C2参数
void coil_furnace_coefficient();                       //计算C1参数

int Heuristic(int dream);
void GenFirstRMP(IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost);
int DynamicPro(double dual, double *dual_vars);
//struct Knap  Knapsack(struct Knap *p,int head[],double *price,int center,int start,int furnace,int dream);
float  Knapsack(float p[][2],int head[],float *price,int center,int start,int furnace,int dream,COb *Tempcoil1);
struct Knap FindBiggest(double *price,int furnace);
float traceback(float **p,int head[],double *dual_vars,int center,float *price,int fur);
void FindReducost(int num,struct furnace *tempfur,double reducost);
void up_linear_mode(IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost);
int CacuDert(struct furnace *nodepath, double *dual_vars);

void ekk_init_linear();
int  ekk_solve_linear( );
void WriteMps( double *a, int m ,int n);
int GenerateCuts(int NumbInCol,int old_col_num,int NEW_COL_NUM,int possible,IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost, int b,int c);
float Imp_Upbound(CList<CFurnace,CFurnace> &Flist,float cur_value,int *UCoil);


float absf(float a);                                    //绝对值函数
struct furnace* get_fur(int No);
double Round( double x,int len);
void pop_path(struct furnace **nodepath);
void   push_path( struct furnace **nodepath);
void Set_sequence(int *sequence);
void Set_wsequence(int *sequence);
void add_fur( struct furnace**nodepath);
void add_fur(int *coilin, float w,float h,double cost,int len);
void release_workspace(void);
void ekk_release( );
double GetValue(CFurnace tempfur);
void branch_and_price();
void   GenerateNode(struct n *nownode);
void   pop(struct  n **node);
void   push(struct n **node);
int linear_relax(IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost);
void write_node(struct n *node,int NumbInCol,int old_col_num,int NEW_COL_NUM,struct n *nownode,int possible);
void    Network_arc( void);
void update_arc(int node_num, struct n *node);
void sort(struct n **list);
void GetBranchArc( struct n *node);
int GenNewRmp( int tree ,struct n *parentnode,IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost);
int CheckBranch(CFurnace p,int r,int s,int be_median, int fur, int coi);
int solve_linear(IloEnv *env,IloModel *mod,IloRangeArray range,IloNumVarArray Lamada,IloNumArray C,IloObjective *cost);
void output(float gap, float time,int dream);
CFurnace Getpit(int No);
float give_cost_in_knap(int *a,int b,int n,int center);
void get_collect();
float CheckConstraint();//在子问题求解过程中，用于判断是否经其方案中加入割约束所对应的对偶变量
int Get_Number_Upper();
int GenerateCutss(int NumbInCol,int old_col_num,int NEW_COL_NUM,int possible,IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost,int b,int c);
void Variable_Fixing();
FILE *ffp;
int sett;
/***************************************主程序*******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>
#define  CRITICAL_SECTION   pthread_mutex_t
#define  _vsnprintf         vsnprintf
#endif
//Log{
#define MAXLOGSIZE 20000000
#define MAXLINSIZE 16000
#include <time.h>
#include <sys/timeb.h>
#include <stdarg.h>
char logfilename1[]="MyLog1.log";
char logfilename2[]="MyLog2.log";
static char logstr[MAXLINSIZE+1];
char datestr[16];
char timestr[16];
char mss[4];
CRITICAL_SECTION cs_log;
FILE *flog;
#ifdef WIN32
void Lock(CRITICAL_SECTION *l) {
	EnterCriticalSection(l);
}
void Unlock(CRITICAL_SECTION *l) {
	LeaveCriticalSection(l);
}
#else
void Lock(CRITICAL_SECTION *l) {
	pthread_mutex_lock(l);
}
void Unlock(CRITICAL_SECTION *l) {
	pthread_mutex_unlock(l);
}
#endif
void LogV(const char *pszFmt,va_list argp) {
	struct tm *now;
	struct timeb tb;

	if (NULL==pszFmt||0==pszFmt[0]) return;
	_vsnprintf(logstr,MAXLINSIZE,pszFmt,argp);
	ftime(&tb);
	now=localtime(&tb.time);
	sprintf(datestr,"%04d-%02d-%02d",now->tm_year+1900,now->tm_mon+1,now->tm_mday);
	sprintf(timestr,"%02d:%02d:%02d",now->tm_hour     ,now->tm_min  ,now->tm_sec );
	sprintf(mss,"%03d",tb.millitm);
	printf("%s %s.%s %s",datestr,timestr,mss,logstr);
	flog=fopen(logfilename1,"a");
	if (NULL!=flog) {
		fprintf(flog,"%s %s.%s %s",datestr,timestr,mss,logstr);
		if (ftell(flog)>MAXLOGSIZE) {
			fclose(flog);
			if (rename(logfilename1,logfilename2)) {
				remove(logfilename2);
				rename(logfilename1,logfilename2);
			}
		} else {
			fclose(flog);
		}
	}
}
void Log(const char *pszFmt,...) {
	va_list argp;

	Lock(&cs_log);
	va_start(argp,pszFmt);
	LogV(pszFmt,argp);
	va_end(argp);
	Unlock(&cs_log);
}
//Log}
int _tmain(int argc, _TCHAR* argv[])
{
    //sett<=5是有效不等式+fixing，6-10，+子问题有效不等式，11-15+动态上界
	
	//int nFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
	//nFlag |= _CRTDBG_LEAK_CHECK_DF;
	//_CrtSetDbgFlag(nFlag);

#ifdef WIN32
	InitializeCriticalSection(&cs_log);
#else
	pthread_mutex_init(&cs_log,NULL);
#endif

	//	int a[2][31]={{0,100,235,100,246,114,175,145,107,162,201,239,216,190,132,88,145,135,87,120,290},{0,45,49,10,40,16,27,20,22,24,35,38,39,34,16,15,16,10,22,16,68}};
	int a[2][31]={{0,54,52,40,56,65,77,45,47,57,79,73,72,58,68,56,48,65,51,56,75},{0,45,49,10,40,16,27,20,22,24,35,38,39,34,16,15,16,10,22,16,68}};
	for(sett=7;sett<=7;sett++)
	{
		    char* filename;
			char dd[100]="Result_set";
			filename=dd;
			char C1[2],C2[3];
			itoa(sett,C1,10);
			strcpy(filename,"Result_set");//把后面指向的字符串拷贝到前面的指针指向的地方
			strcat(filename,"_");
			strcat(filename,C1);
			strcat(filename,".txt");
			//if(sett==5||sett==10||sett==15)
			//	continue;


			if((ffp=fopen(filename,"wt+"))==NULL)
			{
				printf("Cannot open file strike any key exit 1!");
				exit(1);
			}
			//if(sett==5||sett==10||sett==15)
			//	continue;
		for(int dream=1;dream<=1;dream++)
		{
			//if(dream==5||dream==9||dream==12||dream==13||dream==14)

			/*if(dream!=5&&dream!=10&&dream!=15)
				continue;*/

			//if(dream==3||dream==6||dream==8||dream==12||dream==13)
			//	continue;
			is_optimal=true;
			Sub_problem_Num=0;

			//Fur_numb=5;
			//Col_numb=a[0][dream];

			if(sett==1||sett==6|sett==11)
			{
				Fur_numb=4;
				Col_numb=40;
			}
			else if(sett==2||sett==7||sett==12)
			{
				Fur_numb=6;
				Col_numb=60;
			}
			else if(sett==3||sett==8||sett==13)
			{
				Fur_numb=8;
				Col_numb=80;
			}
			else if(sett==4||sett==9||sett==14)
			{
				Fur_numb=10;
				Col_numb=100;
			}
			else if(sett==5||sett==10||sett==15)
			{
				Fur_numb=12;
				Col_numb=120;
			}

			NHB=0;
			NHS=0;
			HHB=0;
			HHS=0;
			sub_time=0;
			valid_number=0;
			valid_total_num=0;
			total_column_number=0;
			int rtn;
			int rt;
			clock_t tick1,tick2;      //时间变量用于记录算法运行时间

			dkk=(float **)malloc((Col_numb+1)*sizeof(float*));
			eij=(float **)malloc((Col_numb+1)*sizeof(float*));
			for (int i=0;i<=Col_numb;i++)
			{
				dkk[i]=(float *)malloc((Fur_numb+1)*sizeof(float));
				eij[i]=(float *)malloc((Fur_numb+1)*sizeof(float));
			}

			init_data(dream);              //读入数据

			nodenum=0;
			NEW_COL_NUM=0;
			result=NULL;
			nodelist=NULL;
			//		Get_Number_Upper();
			Heuristic(dream);
			begin_time=clock();
			branch_and_price();

			finish_time=clock();

			float Gap=(Round((LBinBoot - result->lower_bound )/LBinBoot ,4));
			printf("\nresult is %f,lowerbound is: %f",result->lower_bound,LBinBoot);
			printf("\nGAP is: %f",Gap);
			printf("\n run time is: %f",(float)(finish_time-begin_time)/1000);
			output(Gap,(float)(finish_time-begin_time)/1000,dream);

			fprintf(ffp,"group-%d coil-%d NHB-%d NHS-%d HHB-%d HHS-%d:\t",dream,Col_numb,NHB,NHS,HHB,HHS);
			fprintf(ffp,"对偶间隙为 %f	%d\t",Gap,is_optimal);
			fprintf(ffp,"运行时间为 %f\t",(float)(finish_time-begin_time)/1000);
			fprintf(ffp,"目标 %f,上界 %f\t",result->lower_bound,LBinBoot);
			fprintf(ffp,"节点数为 %d\t",nodenum);
			fprintf(ffp,"子问题迭代次数 %d\t",Sub_problem_Num);
			fprintf(ffp,"\n");
			release_workspace( );
		}
		fclose(ffp);
	}
	
#ifdef WIN32
	DeleteCriticalSection(&cs_log);
#else
	pthread_mutex_destroy(&cs_log);
#endif
}

/*************************************读文件函数*****************************************/
void init_data(int dream)
{

	int empty_num;
	FILE *fp1=NULL;
	FILE *fp2=NULL;
	FILE *fp3=NULL;

	char file_name[600]="coil";
	char C[2];
	itoa(dream,C,10);
	strcat(file_name,C);
	strcat(file_name,".txt");
	if((fp1=fopen(file_name,"r"))==NULL)
	{
		printf("板卷数据未打开!\n");
		exit(0);
	}
	if((fp2=fopen("板卷数据.txt","w+"))==NULL)
		return;

	coil_data(fp1,fp2,dream);                        //数据读入数组中

	coil_coil_coefficient();     //C2
	/*	for(int m=1;m<=Col_numb;m++)
	for(int n=1;n<=Col_numb;n++)
	printf("C2[%d][%d]=%8f",m,n,C2[m][n]);*/

	//     printf("请输入此刻空闲炉子的数目:\n");
	//     scanf("%3d",&empty_num);

	//empty_num=Fur_numb;
	empty_num=5;

	char fur_name[600]="炉子数据";
	strcat(fur_name,C);
	strcat(fur_name,".txt");
	//if((fp3=fopen(fur_name,"r"))==NULL) 
	//{
	//	printf("炉子数据未打开\n");	
	//	exit(0);
	//}

	for(int i=1;i<=empty_num;i++)
	{
		//fscanf(fp3,"%d   %d \n",&batch[i].No,&batch[i].type_num);// 输入炉子各值
		batch[i].type_num=i;
		//if(batch[i].type_num==1)
		//	NHB++;
		//if(batch[i].type_num==2)
		//	HHB++;
		//if(batch[i].type_num==3)
		//	NHS++;
		//if(batch[i].type_num==4)
		//	HHS++;
	}
	//fclose(fp3);
	coil_furnace_coefficient();  //C1
	Network_arc();

	/*	T[1]=160;Height[1]=4685;pj[1]=50;
	T[2]=160;Height[2]=4685;pj[2]=50;
	T[3]=105;Height[3]=4685;pj[3]=50;
	T[4]=105;Height[4]=4685;pj[4]=50;*/
}  

/***************************************读入参数******************************************/

void coil_data(FILE *read_file,FILE *fp,int dream)
{
	/*	for(int i=1;i<=Col_numb;i++)                                   //读取coilsi.text文件
	{
	fscanf(read_file,"%d %f %f %f %d %s  %f  %d  %d",&(coil_inf[i].number),&(coil_inf[i].thickness),&(coil_inf[i].weight),&(coil_inf[i].width),&(coil_inf[i].curve),&(coil_inf[i].contract),&(coil_inf[i].outerdia),&(coil_inf[i].duedate),&coil_inf[i].degree);
	coil[i].number=coil_inf[i].number;
	coil[i].w=coil_inf[i].weight;
	coil[i].h=coil_inf[i].width;
	coil[i].f=50;  // 惩罚费用，可变化的
	fprintf(fp,"%10d  %7.3f  %7.3f  %3f\n",coil[i].number,coil[i].w,coil[i].h,coil[i].f);//生成需要的文本文件：板卷数据.txt
	coil[i].flag=0;
	}

	fclose(read_file);
	fclose(fp);*/
	char buffer1[100];
	char buffer2[100];
	char buffer3[100];
	char buffer4[100];
	char buffer5[100];
	char buffer6[100];
	char buffer7[100];
	char buffer8[100];
	char buffer9[100];
	char buffer10[100];
	char buffer11[100];
	float thick,weight,Outer,m_nProcessLess,m_nPermitPass;
	for(int i=1;i<=Col_numb;i++)                                   //读取coilsi.text文件&(pc[i].m_nThick)
	{
		//		coil_inf[i]=new COb;
		//		fscanf(read_file,"%s %f %f %d %s   %s  %s  %s  %d  %s  %s  %s  %s	%s",buffer1,&thick,&weight,&(coil_inf[i].m_nWidth),buffer2,buffer3,buffer4,buffer5,&(coil_inf[i].m_nOuter),buffer6,buffer7,buffer8,buffer9,buffer10);
		fscanf(read_file,"%s	%f	%d	%f	%s	%s	%s	%s	%f	%s	%s	%s	%s	%s	%f	%f	%s",buffer1,&(thick),&(coil_inf[i].m_nWidth),&(weight),buffer4,buffer3,
			buffer6,buffer7,&(Outer),buffer2,buffer8,buffer5,buffer9,buffer11,&(m_nPermitPass),&(m_nProcessLess),buffer10);
		coil_inf[i].m_strCoilNum=buffer1;
		coil_inf[i].m_strCurve=buffer2;
		coil_inf[i].m_strTrademark=buffer3;
		coil_inf[i].m_strOutmark=buffer4;
		coil_inf[i].m_strOrderNum=buffer5;
		coil_inf[i].m_nStuffGroup=buffer6;
		coil_inf[i].m_nSupplNum=buffer7;
		coil_inf[i].m_strDate=buffer8;
		coil_inf[i].m_strDegree=buffer9;
		coil_inf[i].m_strDegrease=buffer10;
		coil_inf[i].m_strProductDate=buffer11;
		coil_inf[i].index=i;
		coil_inf[i].number=i;
		coil_inf[i].SetWeight(weight);
		coil_inf[i].SetThick(thick);
		coil_inf[i].SetOuter(Outer);
		coil_inf[i].SetProcessLess(m_nProcessLess);
		coil_inf[i].SetPermitPass(m_nPermitPass);
		coil_inf[i].index=i;
		coil_inf[i].SetPRI(1,dream);

		coil[i].number=coil_inf[i].m_strCoilNum;
		coil[i].w=coil_inf[i].m_nWeight;
		coil[i].h=coil_inf[i].m_nWidth;
		coil[i].f=50;  // 惩罚费用，可变化的
		coil[i].flag=0;
		coil[i].PRI=coil_inf[i].PRI;
	}
	fclose(read_file);
}

/*************************************计算C2参数******************************************/
void coil_coil_coefficient()
{
	/*	int i,j;
	float w1,w2,w3,w4,w5; //退火曲线  合同号  外径差  宽度差 厚度差 各个因素的影响系数(待修改)
	w1=1;
	w2=1;
	w3=1;
	w5=1;
	int t1=1,t2=1,t3,t4=10;

	for(i=1;i<=Col_numb;i++)  // 退火曲线
	for(j=1;j<=Col_numb;j++)
	{
	if(coil_inf[i].CurveMatch(coil_inf[j]))
	C2[i][j]+=w1*0;
	else
	C2[i][j]+=w1*Big_number;
	}
	for(i=1;i<=Col_numb;i++)      
	for(j=1;j<=Col_numb;j++)
	{
	if(((coil_inf[i].m_nThick>=1.1)&&(coil_inf[j].m_nThick>=1.1))||((coil_inf[i].m_nThick<1.1)&&(coil_inf[j].m_nThick<1.1)) )    // 厚度差
	C2[i][j]+=w5*t4*absf(coil_inf[i].m_nThick-coil_inf[j].m_nThick);
	else if(((coil_inf[i].m_nThick>=1.1)&&(coil_inf[j].m_nThick<1.1))||((coil_inf[j].m_nThick>=1.1)&&(coil_inf[i].m_nThick<1.1)))
	C2[i][j]+=w5*Big_number;
	}*/
	FILE *fp;
	if((fp=fopen("coil_coil_cost.txt","wt+"))==NULL)
	{
		printf("Cannot open file strike any key exit 1!");
		exit(1);
	}
	for(int i=1;i<=Col_numb;i++)
	{
		for(int j=i;j<=Col_numb;j++)
		{
			if(coil_inf[j].MixMatch(coil_inf[i]))
			{
				if(coil_inf[j].CurveMatch(coil_inf[i]))
					C2[i][j]+=0;
				else
					C2[i][j]+=0.5*abs(atoi(coil_inf[i].m_strCurve)-atoi(coil_inf[j].m_strCurve));
				C2[i][j]+=2*abs(coil_inf[i].m_nThick-coil_inf[j].m_nThick);
				//if(abs(coil_inf[i].m_nOuter-coil_inf[j].m_nOuter)>50)
				    C2[i][j]+=0.017*abs(coil_inf[i].m_nOuter-coil_inf[j].m_nOuter);
			}
			else
			{
				if(coil_inf[j].CurveMatch(coil_inf[i]))
					C2[i][j]+=0;
				else
					C2[i][j]+=MAX;
				if(coil_inf[j].ThickMatch(coil_inf[i]))
					C2[i][j]+=2*abs(coil_inf[i].m_nThick-coil_inf[j].m_nThick);
				else
					C2[i][j]+=MAX;
				//if(abs(coil_inf[i].m_nOuter-coil_inf[j].m_nOuter)>50)
				    C2[i][j]+=0.017*abs(coil_inf[i].m_nOuter-coil_inf[j].m_nOuter);
			}
			if(i==j)
				C2[i][j]=0;
			C2[j][i]=C2[i][j];
			C2[j][i]=C2[j][i]*1000;
			C2[i][j]=C2[i][j]*1000;
			if(C2[j][i]>=10000000)
			{
				arc_matrix[i][j]=4;
				arc_matrix[j][i]=4;
			}
		}
	}
	for(int i=1;i<=Col_numb;i++)
	{
		for(int j=1;j<=Col_numb;j++)
		{
			if(C2[i][j]<MAX)
			    fprintf(fp,"%d	%d	%d	%lf	%lf\n",i,j,coil_inf[j].m_nWidth,C2[i][j],coil_inf[j].PRI+coil_inf[j].m_nWeight);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


/****************************************计算C1参数************************************/
void coil_furnace_coefficient()  //C1系数的产生
{
	/*	int t5=5,t6=10,t7=20,t8=3,t9=7,t10=17,t11=8,t12=13,t13=23,t14=1;
	int i,j,temp;
	int beginning_date=20040801;

	for(i=1;i<=Col_numb;i++)      
	for(j=1;j<=Fur_numb;j++)
	{
	temp=coil_inf[i].duedate-beginning_date;
	if(coil_inf[i].degree==4)
	{
	if(temp<=9)
	{
	if(coil_inf[i].curve==11 ||coil_inf[i].curve==12 ||coil_inf[i].curve==13||coil_inf[i].curve==33||coil_inf[i].curve==34||coil_inf[i].curve==66||coil_inf[i].curve==67||coil_inf[i].curve==68)
	C1[i][j]+=t8*0.7*temp;
	else
	C1[i][j]+=t8*temp;
	}
	else if((9<temp)&&(temp<=16))
	{
	if(coil_inf[i].curve==11 ||coil_inf[i].curve==12 ||coil_inf[i].curve==13||coil_inf[i].curve==33||coil_inf[i].curve==34||coil_inf[i].curve==66||coil_inf[i].curve==67||coil_inf[i].curve==68)
	C1[i][j]+=t9*0.7*temp;
	else
	C1[i][j]+=t9*temp;
	}
	else
	{
	if(coil_inf[i].curve==11 ||coil_inf[i].curve==12 ||coil_inf[i].curve==13||coil_inf[i].curve==33||coil_inf[i].curve==34||coil_inf[i].curve==66||coil_inf[i].curve==67||coil_inf[i].curve==68)
	C1[i][j]+=t10*0.7*temp;
	else
	C1[i][j]+=t10*temp;
	}
	}
	else
	{
	if(temp<=9)
	{
	if(coil_inf[i].curve==11 ||coil_inf[i].curve==12 ||coil_inf[i].curve==13||coil_inf[i].curve==33||coil_inf[i].curve==34||coil_inf[i].curve==66||coil_inf[i].curve==67||coil_inf[i].curve==68)
	C1[i][j]+=t11*0.7*temp;
	else
	C1[i][j]+=t11*temp;
	}
	else if((9<temp)&&(temp<=16))
	{
	if(coil_inf[i].curve==11 ||coil_inf[i].curve==12 ||coil_inf[i].curve==13||coil_inf[i].curve==33||coil_inf[i].curve==34||coil_inf[i].curve==66||coil_inf[i].curve==67||coil_inf[i].curve==68)
	C1[i][j]+=t12*0.7*temp;
	else
	C1[i][j]+=t12*temp;
	}
	else
	{
	if(coil_inf[i].curve==11 ||coil_inf[i].curve==12 ||coil_inf[i].curve==13||coil_inf[i].curve==33||coil_inf[i].curve==34||coil_inf[i].curve==66||coil_inf[i].curve==67||coil_inf[i].curve==68)
	C1[i][j]+=t13*0.7*temp;
	else
	C1[i][j]+=t13*temp;
	}
	}

	if((batch[j].type_num==1)||(batch[j].type_num==2))
	{
	if((1950<=coil_inf[i].outerdia)&&(coil_inf[i].outerdia<2050))
	C1[i][j]+=t14*(2050-coil_inf[i].outerdia);
	if(coil_inf[i].outerdia<1950.0)
	C1[i][j]+=Big_number;
	}
	else
	{
	if(coil_inf[i].outerdia>2050)
	C1[i][j]+=Big_number;
	}
	}*/
	int a=0;
	for(int i=1;i<=Col_numb;i++)
	{
		for(int j=1;j<=Fur_numb;j++)
		{
			if((batch[j].type_num==1)||(batch[j].type_num==2))
			{
				C1[i][j]+=0.01*(2550-coil_inf[i].m_nOuter);//0
			}
			else
			{
				if(coil_inf[i].m_nOuter>=2050)
					C1[i][j]+=MAX;
				else
					C1[i][j]+=0.01*(2050-coil_inf[i].m_nOuter);//0
			}
			if((batch[j].type_num==1)||(batch[j].type_num==3))
			{
				if(atoi(coil_inf[i].m_strCurve)>=60)
					C1[i][j]+=MAX;
			}
			else if((batch[j].type_num==2)||(batch[j].type_num==4))
			{
				if(atoi(coil_inf[i].m_strCurve)<60)
					C1[i][j]+=0.1*(68-atoi(coil_inf[i].m_strCurve));
			}
			C1[i][j]=C1[i][j]*1000;

		}
	}
}

/*************************************启发式获得初始解***********************************/
int Heuristic(int dream)
{
	double cost=0;
	int *coilin;
	float height=0;
	float weight=0;
	int X[C_MAX][F_MAX];
	int parater[8];
	parater[0]=Col_numb;
	int aver=Col_numb/3;
	aver=aver/4;

	if(sett==1||sett==6||sett==11)
	{
		parater[1]=1;//aver;//NHB
	    parater[2]=1;//aver;//NHS
	    parater[3]=1;//aver;//HHB
	    parater[4]=1;//aver;//HHS
	    parater[5]=0;
	}
	else if(sett==2||sett==7||sett==12)
	{
		parater[1]=2;//aver;//NHB
	    parater[2]=1;//aver;//NHS
	    parater[3]=2;//aver;//HHB
	    parater[4]=1;//aver;//HHS
	    parater[5]=0;
	}
	else if(sett==3||sett==8||sett==13)
	{
		parater[1]=2;//aver;//NHB
	    parater[2]=2;//aver;//NHS
	    parater[3]=2;//aver;//HHB
	    parater[4]=2;//aver;//HHS
	    parater[5]=0;
	}
	else if(sett==4||sett==9||sett==14)
	{
		parater[1]=4;//aver;//NHB
	    parater[2]=3;//aver;//NHS
	    parater[3]=2;//aver;//HHB
	    parater[4]=1;//aver;//HHS
	    parater[5]=0;
	}
	else if(sett==5||sett==10||sett==15)
	{
		parater[1]=3;//aver;//NHB
	    parater[2]=3;//aver;//NHS
	    parater[3]=3;//aver;//HHB
	    parater[4]=3;//aver;//HHS
	    parater[5]=0;
	}
	CList<COb,COb> list1,list2;
	CList<CStack,CStack> list3;
	CMakeStack make;
	totalcost=0;
	make.mainS(list1,list2,list3,coil_inf,parater,dream);
	float maxweight=0;
	int maxw=0;
	int num=0;
	Fur_numb=list3.GetCount();
	for(int i=1;i<=Col_numb;i++)
	{
		coil_inf[i].PRI=coil_inf[i].PRI*1000;
		coil_inf[i].m_nWeight=coil_inf[i].m_nWeight*1000;
	}
	for(int i=0;i<Fur_numb;i++)
	{
		num=0;
		POSITION pos=list3.GetHeadPosition();
		CStack ts=list3.GetAt(pos);
		coilin=new int[301];
		maxweight=ts.coil[0].m_nWeight;
		maxw=0;
		for(int j=0;j<XN;j++)
		{
			if(ts.coil[j].m_nOuter!=0)
			{
				coilin[j+1]=ts.coil[j].index;
				height=ts.SumWidth();
				weight=ts.SumWeight();
				if(ts.coil[j].m_nWeight>maxweight)
				{
					maxweight=ts.coil[j].m_nWeight;
					maxw=j;
				}
				num++;
				ASSERT(coilin[j+1]>0);
			}
		}

		CFurnace newpath;
		newpath.No=pitlist.GetCount()+1;

		newpath.lenth=num;
		if(ts.furnace.Fur_type=="NHB")
		{
			NHB++;
			newpath.index=1;
		}
		else if(ts.furnace.Fur_type=="HHB")
		{
			newpath.index=2;
			HHB++;
		}
		else if(ts.furnace.Fur_type=="NHS")
		{
			newpath.index=3;
			NHS++;
		}
		else if(ts.furnace.Fur_type=="HHS")
		{
			newpath.index=4;
			HHS++;
		}
		// 	newpath.coils_in=coilin;
		newpath.weight=weight;
		newpath.height=height;
		newpath.center=coilin[1];
		for(int m=1;m<=num;m++)
		{
			newpath.ingot.AddTail(coil_inf[coilin[m]]);
			newpath.coils_in[m]=coilin[m];
		}

		cost=GetValue(newpath);
		newpath.cost=cost;
		totalcost+=cost;
		total_column_number++;
		pitlist.AddTail(newpath);
		pit_pool.AddTail(newpath); 
		delete []coilin;
		list3.RemoveHead();
	}
 
	printf("%f\n",totalcost);
	Fur_numb = 4;
	//      Fur_numb=list3.GetCount();
	return 1;
}

/***********************************获得第一个限制主问题********************************/
void GenFirstRMP(IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost)
{


	*env = mod->getEnv();

	int column_set[C_MAX+1][F_MAX+1];
	int i,j;
	//	struct furnace *temp_fur;
	double temp_cost;
	double *RMP_cost;
	int RN;

	int size=pitlist.GetCount();
	NEW_COL_NUM=size;
	con_m=Col_numb+Fur_numb;  //变量数加1[对应于线性规划矩阵的列]
	var_n=size;     //[对应于线性规化的行]
	RN=var_n+1;        //用于存储矩阵

	struct   n * feasible;

	
	a=(double *)calloc((con_m+1)*(var_n+1),sizeof(double));
	RMP_cost =new double[size+1]; //费用列
	POSITION pos=pitlist.GetHeadPosition();
	i=1;
	while(pos)
	{
		CFurnace temp_fur=pitlist.GetAt(pos);
		for(int a=1;a<=temp_fur.lenth;a++)
			printf("  coil[%d]  ",temp_fur.coils_in[a]);

		RMP_cost[i]=temp_fur.cost;
		for(j=1;j<=Col_numb;j++)
			column_set[j][i]=0;
		for(j=1;j<=temp_fur.lenth;j++)
		{
			column_set[temp_fur.coils_in[j]][i]++;
		}
		FNlist.AddTail(temp_fur.No);
		pitlist.GetNext(pos);
		i++;

	}//写到这了

	temp_cost=0.0;
	for(j=1;j<=size;j++)
		temp_cost+=RMP_cost[j];

	feasible=new struct n;
	feasible->lower_bound=temp_cost;
	feasible->level=0;
	feasible->column_num=NEW_COL_NUM;
	feasible->Col=NULL;
	result=feasible;

	for(i=1;i<=size;i++)
		for(int j=1;j<=Col_numb;j++)
			if(column_set[j][i]==1)
				printf("\ncolumn_set[%d][%d]=1",j,i);

	for(j=1;j<=size;j++)
	{
		C->add(IloNum(RMP_cost[j]));
	}


	//约束的左右边项
	//IloNumArray  ConstraintMin(env), ConstraintMax(env);
	IloNumArray  ConstraintMin(*env);
	IloNumArray  ConstraintMax(*env);

	//IloObjective cost;
	//for(i=0;i<Fur_numb;i++)      //方案唯一性约束左右端项
	//{
		ConstraintMin.add(IloNum(-IloInfinity));
		ConstraintMax.add(IloNum(NHB));//1

		ConstraintMin.add(IloNum(-IloInfinity));
		ConstraintMax.add(IloNum(HHB));

		ConstraintMin.add(IloNum(-IloInfinity));
		ConstraintMax.add(IloNum(NHS));

		ConstraintMin.add(IloNum(-IloInfinity));
		ConstraintMax.add(IloNum(HHS));


	//}
	for(i=Fur_numb;i<Fur_numb+Col_numb;i++) //钢锭唯一性约束左右端项
	{
		ConstraintMin.add(IloNum(-IloInfinity));//-INFINITE-IloInfinity
		ConstraintMax.add(IloNum(1));
	}
	

	*range= IloRangeArray(*env, ConstraintMin, ConstraintMax);
	mod->add(*range);

	*cost = IloAdd(*mod, IloMaximize(*env));



	pos=pitlist.GetHeadPosition();
	int fur_num=1;
	while (pos)
	{
		CFurnace temp_fur=pitlist.GetAt(pos);
		IloNumColumn col = (*cost)((*C)[fur_num-1]);
		for (i = 0; i < Fur_numb+Col_numb; i++) 
		{
			if(i==temp_fur.index-1)
			{
				col+=(*range)[i](1);      //方案唯一性约束
			}
			else if(i<Fur_numb)            
			{
				col += (*range)[i](0);
			}
			else if(i<Fur_numb+Col_numb)                         //每个钢锭只能放在一个炉子中
			{
				col += (*range)[i](column_set[i-Fur_numb+1][fur_num]);
			}
		}

		Lamada->add(IloNumVar(col,0,1,ILOFLOAT));

		col.end();
		pitlist.GetNext(pos);
		fur_num++;
	}

	delete []RMP_cost;
	free(a);
}

struct furnace* get_fur(int No)
{    	
	struct furnace *now;
	
	return now;
}

void ekk_init_linear()
{
	
}

void WriteMps( double *a, int m ,int n)
{
	int  i,j;
	double  up,low;
	FILE  *simp;
	simp=fopen("simp.mps","w");

	fprintf(simp,"NAME          EXAMPLE\n");
	fprintf(simp,"OBJSENSE\n");
	fprintf(simp,"MAXIMIZE\n");
	fprintf(simp,"ROWS\n");
	for(i=0;i<=m;i++)
	{

		if(i==0)
			fprintf(simp," N  OBJ\n");
		else if (i<=Fur_numb)
			fprintf(simp," E  ROW%d\n",i);
		// 		    else if (nownode->model_modify[i-Pit_num]==1&&(change_flag||nownode->is_change))
		// 			    fprintf(simp," E  ROW%d\n",i);
		else
			fprintf(simp," L  ROW%d\n",i);

	}

	fprintf(simp,"COLUMNS\n");

	for(j=1;j<=n;j++){

		// 		 if(a[j]!=0)
		fprintf(simp,"    COL%-5d   OBJ        %12lf\n",j,a[j]);


		for(i=1;i<=m;i++){
			if(a[i*(n+1)+j]>=0.0000001)
				fprintf(simp,"    COL%-5d   ROW%-5d   %12lf\n",j,i,a[i*(n+1)+j]);
		}
		//fprintf(simp,"\n");

	}

	fprintf(simp,"RHS \n");

	for(i=1;i<=m;i++)
		fprintf(simp,"    RHS1       ROW%-5d   %12lf\n",i,a[i*(n+1)]);


	fprintf(simp,"RANGES \n");
	fprintf(simp,"BOUNDS \n");

	up=1.0;
	low=0.0;
	for(j=1;j<=n;j++)
	{
		fprintf(simp," LO BND1      COL%-5d   %12lf\n",j,low);
		fprintf(simp," UP BND1      COL%-5d   %12lf\n",j,up);

	}
	fprintf(simp,"ENDATA");
	fclose(simp);
	////////////////////////
	FILE  *A;
	A=fopen("a.txt","w");

	fprintf(A,"NAME          EXAMPLE\n");
	fprintf(A,"OBJSENSE\n");
	fprintf(A,"MAXIMIZE\n");
	fprintf(A,"ROWS\n");
	for(i=0;i<=m;i++){

		if(i==0)
			fprintf(A," N  OBJ\n");
		else
			fprintf(A," E  ROW%d\n",i);

	}

	fprintf(A,"COLUMNS\n");

	for(j=1;j<=n;j++){

		if(a[j]!=0)
			fprintf(A,"    COL%-5d   OBJ        %12lf\n",j,a[j]);


		for(i=1;i<=m;i++){
			if(a[i*(n+1)+j]>=0.0000001)
				fprintf(A,"    COL%-5d   ROW%-5d   %12lf\n",j,i,a[i*(n+1)+j]);
		}
		//fprintf(simp,"\n");

	}

	fprintf(A,"RHS \n");

	for(i=1;i<=m;i++)
		fprintf(A,"    RHS1       ROW%-5d   %12lf\n",i,a[i*(n+1)]);


	fprintf(A,"RANGES \n");
	fprintf(A,"BOUNDS \n");

	up=1.0;
	low=0.0;
	for(j=1;j<=n;j++)
	{
		fprintf(A," LO BND1      COL%-5d   %12lf\n",j,low);
		fprintf(A," UP BND1      COL%-5d   %12lf\n",j,up);

	}
	fprintf(A,"ENDATA");
	fclose(A);
	/*    int  i,j;
	double  up,low;
	FILE  *simp;
	simp=fopen("simp.mps","w");

	fprintf(simp,"NAME          EXAMPLE\n");
	fprintf(simp,"ROWS\n");
	for(i=0;i<=m;i++)
	{

	if(i==0)
	fprintf(simp," N  OBJ\n");
	else if (i==1)
	fprintf(simp," E  ROW%d\n",i);
	else
	fprintf(simp," L  ROW%d\n",i);

	}

	fprintf(simp,"COLUMNS\n");

	for(j=1;j<=n;j++){

	if(a[j]!=0)
	fprintf(simp,"    COL%-5d   OBJ        %12lf\n",j,a[j]);


	for(i=1;i<=m;i++){
	if(a[i*(n+1)+j]>=0.0000001)
	fprintf(simp,"    COL%-5d   ROW%-5d   %12lf\n",j,i,a[i*(n+1)+j]);
	}
	//fprintf(simp,"\n");

	}

	fprintf(simp,"RHS \n");

	for(i=1;i<=m;i++)
	fprintf(simp,"    RHS1       ROW%-5d   %12lf\n",i,a[i*(n+1)]);


	fprintf(simp,"RANGES \n");
	fprintf(simp,"BOUNDS \n");

	up=1.0;
	low=0.0;
	for(j=1;j<=n;j++)
	{
	fprintf(simp," LO BND1      COL%-5d   %12lf\n",j,low);
	//fprintf(simp," UP BND1      COL%-5d   %12lf\n",j,up);

	}
	fprintf(simp,"ENDATA");
	fclose(simp);
	////////////////////////
	FILE  *A;
	A=fopen("a.txt","w");

	fprintf(A,"NAME          EXAMPLE\n");
	fprintf(A,"ROWS\n");
	for(i=0;i<=m;i++){

	if(i==0)
	fprintf(A," N  OBJ\n");
	else
	fprintf(A," E  ROW%d\n",i);

	}

	fprintf(A,"COLUMNS\n");

	for(j=1;j<=n;j++){

	if(a[j]!=0)
	fprintf(A,"    COL%-5d   OBJ        %12lf\n",j,a[j]);


	for(i=1;i<=m;i++){
	if(a[i*(n+1)+j]>=0.0000001)
	fprintf(A,"    COL%-5d   ROW%-5d   %12lf\n",j,i,a[i*(n+1)+j]);
	}
	//fprintf(simp,"\n");

	}

	fprintf(A,"RHS \n");

	for(i=1;i<=m;i++)
	fprintf(A,"    RHS1       ROW%-5d   %12lf\n",i,a[i*(n+1)]);


	fprintf(A,"RANGES \n");
	fprintf(A,"BOUNDS \n");

	up=1.0;
	low=0.0;
	for(j=1;j<=n;j++)
	{
	fprintf(A," LO BND1      COL%-5d   %12lf\n",j,low);
	fprintf(A," UP BND1      COL%-5d   %12lf\n",j,up);

	}
	fprintf(A,"ENDATA");
	fclose(A);*/
}

int   ekk_solve_linear( )
{
		return 0;
}



int DynamicPro(double dual, double *dual_vars)
{
	double *price;
	float *price1;
	int reduce_num=0;
	price=new double[Col_numb+1];
	price1=new float[Col_numb+1];
	bool is_break=false;
	int cc[C_MAX];
	float p[20000][2];

	int mini;
	float cur_dual[C_MAX+1][C_MAX+1]; //计算每两个板卷有效不等式所对应的对偶值

	bool is_extract=true;


	for (int s=0;s<=Col_numb+1;s++)  //初始化
	{
		for (int r=0;r<=Col_numb;r++)
		{
			cur_dual[r][s]=0;
			cur_dual[s][r]=0;
		}
		cc[s]=0;
	}
	int counter=0;
	POSITION cut_pos=CurrentCutList.GetHeadPosition();//此处在求解价格子问题时，计算有效不等式所对应的对偶变量
	while (cut_pos)
	{
		counter++;
		struct AddCut st=CurrentCutList.GetAt(cut_pos);
		for(int s=0;s<st.number;s++)
		{
			cur_dual[st.coil[s][1]][st.coil[s][0]]+=dual_var[Fur_numb+Col_numb+valid_total_num+counter-1];
			cur_dual[st.coil[s][0]][st.coil[s][1]]+=dual_var[Fur_numb+Col_numb+valid_total_num+counter-1];
		}
		CurrentCutList.GetNext(cut_pos);
	}
	NumOfValidInequality = counter;

	COb *temp_coil;//记录暂时每个炉子内可装的板卷的信息
	temp_coil=new COb[Col_numb+1];

	COb *temp_coil1;

	for(int j=1;j<=Fur_numb;j++)//炉子
	{
		for(int k=median_from;k<=median_to;k++)//中心板卷
		{
			bool is_median=false;
			int total_height=coil_inf[k].m_nWidth;
			float total_value=0;
			if(C1[k][j]>C_MAX*1000||fur_coil[j][k]==4||median[k]==4||cc[k]==1||Ykk[k][j]==4)//
				continue;
			int feasible_coil_num=0;
			for (int i=1;i<=Col_numb;i++) //每个板卷设置费用
			{
				if(C1[i][j]<C_MAX*1000&&fur_coil[j][i]!=4&&C2[i][k]<C_MAX*1000&&arc_matrix[i][k]!=4&&i!=k&&median[i]!=4) //找到满足条件的板卷
				{
					feasible_coil_num++;
					total_height+=coil_inf[i].m_nWidth+70;
					temp_coil[feasible_coil_num]=coil_inf[i];
					price[feasible_coil_num]=coil_inf[i].m_nWeight+coil_inf[i].PRI-C1[i][j]-dual_var[i+Fur_numb-1]-C2[i][k]-cur_dual[k][i];//基本模型中的费用，没有考虑有效不等式和分支约束
					total_value+=price[feasible_coil_num];
					if(arc_matrix[k][1]==3)
						is_median=true;
				}
			}
			temp_coil1=new COb[feasible_coil_num+1];
			int n1=0;
			int tt=0;
			if(is_median)
				goto L1;
			for(int m=1;m<=feasible_coil_num/2;m++)
			{
				bool flag=false;
				{
					for (int y=1;y<=feasible_coil_num;y++)
					{
						if (arc_matrix[temp_coil[m].index][temp_coil[y].index]==3||arc_matrix[temp_coil[m].index][temp_coil[y].index]==4||cur_dual[temp_coil[m].index][temp_coil[y].index]!=0)
						{
							flag=true;
						}
					}
					if(flag==false)
					{
						n1++;
						price1[n1]=price[m];
						temp_coil1[n1]=temp_coil[m];
					}
				}
			}

			int head[C_MAX+2];
			float ineq=Knapsack(p,head,price1,n1,1,j,k,temp_coil1)-tt*10000;
			delete []temp_coil1;
            /*建立模型求解价格子问题*/
L1:			sub_begin_time=clock();


			IloEnv env2;
			IloModel mod2(env2);	
			IloBoolVarArray22 Zik(env2);              //辅助决策变量，判断两个板卷是否在一起
			IloExpr CC(env2);                         //目标函数，最大化，总的装包费用
			IloBoolVarArray Ai(env2,feasible_coil_num+1);//决策变量ai，决策板卷是否加入到炉子中
			IloRangeArray Constraints(env2);

			for (int i = 0; i <=feasible_coil_num; i++)
			{
				Zik.add(IloBoolVarArray(env2,feasible_coil_num+1)); //初始化决策变量

			}

			//加入约束
			
			IloExpr availExpr(env2);
			for (int s=1;s<=feasible_coil_num;s++)
			{

				CC+=price[s]*Ai[s];/*加入目标函数*/

				availExpr+=(temp_coil[s].m_nWidth+70)*Ai[s];/*能力约束*/

				for(int r=s+1;r<=feasible_coil_num;r++)/*分支约束*/
				{
					CC=CC-cur_dual[temp_coil[s].index][temp_coil[r].index]*Zik[s][r];/*加入目标函数*/

					if(arc_matrix[temp_coil[r].index][temp_coil[s].index]==4)    //两个板卷一定不再一起
						Constraints.add(Ai[r]+Ai[s]<=1);//mod2.add(Ai[r]+Ai[s]<=1);
					else if(arc_matrix[temp_coil[r].index][temp_coil[s].index]==3)  //两个板卷一定在一起
					{
						Constraints.add(Ai[r]-Ai[s]==0);//mod2.add(Ai[r]==Ai[s]);
					}

					if(cur_dual[temp_coil[s].index][temp_coil[r].index]>0)/*判断两个板卷是否在一起的约束*/
						Constraints.add(Zik[s][r]-Ai[s]-Ai[r]>=-1);//mod2.add(Zik[s][r]>=Ai[s]+Ai[r]-1);
				}
				if(arc_matrix[k][temp_coil[s].index]==4)    //两个板卷一定不再一起
					Constraints.add(Ai[s]==0);//mod2.add(Ai[s]==0);
				else if(arc_matrix[k][temp_coil[s].index]==3)  //两个板卷一定在一起
					Constraints.add(Ai[s]==1);//mod2.add(Ai[s]==1);

			}
			Constraints.add(availExpr<=4700-coil_inf[k].m_nWidth);//mod2.add(availExpr<=4700-coil_inf[k].m_nWidth);
			availExpr.end();
			mod2.add(Constraints);

			if(sett>=6)
				mod2.add(CC>=ineq-0.1);

			sub_finish_time=clock();
			sub_time+=float(sub_finish_time-sub_begin_time)/1000;

			IloObjective obj=IloMaximize(env2, CC);
			mod2.add(obj);
			
			IloCplex cplex(env2);//采用CPLEX求解
			cplex.setParam(IloCplex::AdvInd, 0);
			cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
			cplex.setParam(IloCplex::MIPDisplay,0);
			cplex.setParam(IloCplex::TiLim,3600);
//			cplex.setParam(IloCplex::SolnPoolAGap,0.5);
			cplex.extract(mod2);
			cplex.exportModel("sub_problem.lp");
			
			int dd=cplex.solve();
			

			if(dd!=0)//如果有可行解
			{
				float reduce_cost=cplex.getObjValue();
				reduce_cost+=coil_inf[k].m_nWeight+coil_inf[k].PRI-C1[k][j]-dual_var[k+Fur_numb-1]-dual_var[j-1];
				int num=0;
				for(int l=1;l<=valid_number;l++)
				{
					for(int ll=1;ll<=valid_maste[l].collect_num;ll++)
					{
						num++;
						if(k==valid_maste[l].index)
							reduce_cost+=-dual_var[Fur_numb+Col_numb+num-1];
						else if(k==valid_maste[l].ingotlist[ll])
							reduce_cost+=dual_var[Fur_numb+Col_numb+num-1];
					}

				}
				dkk[k][j]=reduce_cost+dual_var[j-1];
				if(Round(reduce_cost,6)>(double)0.1)
				{
					CFurnace reduce_pit;  //负费用列
					/*添加中心板卷*/
					reduce_pit.lenth=1;
					reduce_pit.ingot.AddTail(coil_inf[k]);  
					reduce_pit.coils_in[reduce_pit.lenth]=k;
					reduce_pit.index=j;
					reduce_pit.type_num=batch[j].type_num;
					reduce_pit.center=k;
					
					/*添加其他的板卷*/
					for(int i=feasible_coil_num;i>=1;i--)
					{
						if(cplex.getValue(Ai[i])>0.9)
						{
							//printf("%d	%lf, ",temp_coil[i].index,price[i]);
							reduce_pit.lenth++;
							reduce_pit.ingot.AddTail(temp_coil[i]);
							reduce_pit.coils_in[reduce_pit.lenth]=temp_coil[i].index;
							
						}
					}
					//printf("\n");

					float minimize=MAX;
					int centerr=-1;
					for(int i=1;i<=reduce_pit.lenth;i++)
					{
						float temp=0;
						for(int i1=1;i1<=reduce_pit.lenth;i1++)
						{
							temp+=C2[reduce_pit.coils_in[i]][reduce_pit.coils_in[i1]];
						}
						if(temp<minimize)
						{
							centerr=reduce_pit.coils_in[i];
							minimize=temp;
						}
					}
					if(centerr!=-1);
					reduce_pit.center=centerr;
					cc[centerr]=1;
					if (reduce_pit.lenth==5)
						int is=0;

					reduce_pit.cost=GetValue(reduce_pit); //获得列所对应的目标函数

					pit_pool.AddTail(reduce_pit);  //将该负费用列加入到列池中
					reduce_num++;

					/*输出所找到的负费用列*/
					if(nownode)
					{
//						fprintf(fp1,"r=%d,s=%d,tree=%d\n",nownode->r,nownode->s,nownode->tree);
					}
//					Log("This is a Fur: %d, median: %d	reduce:%f\n",j,k,reduce_cost);
					for(int i=0;i<reduce_pit.ingot.GetCount();i++)
					{
						COb ing=reduce_pit.ingot.GetAt(reduce_pit.ingot.FindIndex(i));

					}

				}
			}
			//清理空间

			cplex.end();

			mod2.end();

			env2.end();		
		}
		for(int i=median_from;i<=median_to;i++)
			cc[i]=0;
	}
	delete []temp_coil;
//	fclose(fp1);
	delete []price;
	delete []price1;



	if(reduce_num==0)
		return 0;
	else
		return 1;
	
}

float  Knapsack(float p[][2],int head[],float *price,int center,int start,int furnace,int dream,COb *Tempcoil1)
{
	int i,j,k;
	int n=center;
	int left=0,right=0,next=1;
	float y,m;
	float we=0;

	float c;

	c=4700-coil_inf[dream].m_nWidth;

	//	bool is_rs=false;
	head[n+1]=0;
	head[n]=1;
	p[0][0]=0;
	p[0][1]=0;
	for(i=n;i>=1;i--)
	{
	   /*if(i=1)
		int z=0; 
			if(num_code==1)  
				int nn;*/
		k=left;
		for(j=left;j<=right;j++)
		{
			if(p[j][0]+Tempcoil1[i].m_nWidth+70>c)
				break;              //超过限制的时候停止
			y=p[j][0]+Tempcoil1[i].m_nWidth+70; //记录
   			m=p[j][1]+price[i];
			while(k<=right&&p[k][0]<y)
			{
				p[next][0]=p[k][0];
				p[next++][1]=p[k++][1];
			}
			
			if(k<=right && p[k][0]==y)
			{
				if(m<p[k][1])
					m=p[k][1];
				k++;
			}
			if(m>p[next-1][1])
			{
				p[next][0]=y;
				p[next++][1]=m;
			}
			while(k<=right && p[k][1]<=p[next-1][1])
				k++;
		}
		while (k<=right)
		{
			p[next][0]=p[k][0];
			p[next++][1]=p[k++][1];
            
		}
		left=right+1;
		right=next-1;
		head[i-1]=next;
	}
	return p[next-1][1];
}


float traceback(float **p,int head[],double *dual_vars,int center,float *price,int fur)
{
	double reducost=0;
	return reducost;

}

void FindReducost(int num,struct furnace *tempfur,double reducost)
{
}

void up_linear_mode(IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost)
{

	double  objectiveCoefficient;
	double  *element;
	int     *TempCol;
	int     reduce_falg;
	int     *indexRow;
	int     startOfEachColumn;

	int iter=0;
	int k,rtcod,i;
	int *pp;

	int Cut_Num=CurrentCutList.GetCount();
	TempCol= new int[Col_numb+Fur_numb+Cut_Num+valid_total_num+1];
	indexRow= (int *)calloc(con_m,sizeof(int));
	element = (double *)calloc(con_m,sizeof(double));

	while(1)
	{
		CFurnace tempfur;
		if(pit_pool.GetCount()==0)
			break;
		tempfur=pit_pool.GetHead();
		pit_pool.RemoveHead();

		if(tempfur.ingot.GetCount()>0)
		{
			if(1==1)//reduce_falg
			{
				TempCol[0]=1;
				for(i=0;i<=Col_numb+Fur_numb+valid_total_num+Cut_Num-1;i++)
				{
					TempCol[i]=0;
				}
				TempCol[tempfur.index-1]=1;
				for(i=1;i<=tempfur.lenth;i++)
				{
					TempCol[tempfur.coils_in[i]+Fur_numb-1]+=1;
				}


				POSITION cut_pos=CurrentCutList.GetHeadPosition();//判断当先列在哪个有效不等式中
				i=1;
				while (cut_pos)
				{
					struct AddCut st=CurrentCutList.GetAt(cut_pos);
					for(int k=0;k<st.number;k++)
					{
						bool flag1=false;
						bool flag2=false;
						for(int l=1;l<=tempfur.lenth;l++)
						{
							if(st.coil[k][0]==tempfur.coils_in[l])
								flag1=true;
							if(st.coil[k][1]==tempfur.coils_in[l])
								flag2=true;
						}
						if(flag1&&flag2)
						{
							TempCol[Fur_numb+Col_numb+valid_total_num+i-1]=1; //如果当前方案中包括割中的板卷对，则该变量在有效不等式中的系数是1
							break;
						}
					}
					CurrentCutList.GetNext(cut_pos);
					i++;
				}

				C->add(IloNum(tempfur.cost));

				IloNumColumn col = (*cost)(tempfur.cost);

				for (i = 0; i < Fur_numb+Col_numb+valid_total_num+Cut_Num; i++) 
				{

					col+=(*range)[i](TempCol[i]);
				}


				Lamada->add(IloNumVar(col,0,1,ILOFLOAT));

				col.end();


				var_n=var_n+1;
				NEW_COL_NUM+=1;
				total_column_number++;
				tempfur.No=pitlist.GetCount()+1;
				pitlist.AddTail(tempfur); //将此时的负费用列加入到pathlist中
				iter++;
				FNlist.AddTail(tempfur.No);
			}

		}
		else
			break;
		//  		if(iter>=20)//若加的列数大于机器数量则跳出，即每次从列池中选出M_numb个列加到模型中
		//  		    break;
	}
	free(indexRow );
	free(element );
	delete []TempCol;

}


int CacuDert(struct furnace *nodepath, double *dual_vars)
{
	int     flag=0;
	return flag;
}

double Round( double x,int len)
/*-------------------------------------------------------------------------*/
{
	double result;
	double eps;
	double temp;
	double reduce;

	eps=1/(pow(double(10),len));
	result=x/eps;

	if(result>=0)
	{   temp=floor(result);
	reduce=result - temp  ;
	if(reduce<=0.5)
		temp=temp;	   
	else
		temp=temp +1;
	}
	else
	{   temp=ceil(result);
	reduce=temp -result ;
	if(reduce<=0.5)
		temp=temp;	   
	else
		temp=temp -1;
	}

	result=temp*eps;
	return result;

}

float absf(float a)
{
	float c;
	if ( a < 0.0)
		c=-a;
	else
		c=a;
	return(c);
}

void pop_path(struct furnace **nodepath)//从列池中取出列
{
	/*	*nodepath=fur_pool;
	if(fur_pool)
	fur_pool=fur_pool->before;
	if(*nodepath)
	(*nodepath)->before=NULL;*/
}


void push_path (struct furnace **nodepath)
{
}

void add_fur( struct furnace**nodepath)
{
}

void add_fur(int *coilin, float w,float h,double cost,int len)
{

}

void Set_sequence(int *sequence)
{
	float temp;
	int seq;
	float outdia[C_MAX];
	for(int k=1;k<=Col_numb;k++)
	{
		sequence[k]=0;
		outdia[k]=coil_inf[k].m_nOuter;
	}
	int i=1;
	while(1)
	{
		if(i>Col_numb)
			break;

		temp=outdia[1];
		seq=1;
		for(int j=1;j<=Col_numb;j++)
		{
			if(outdia[j]>temp)
			{
				temp=outdia[j];
				seq=j;
			}
		}
		outdia[seq]=0;
		sequence[i]=seq;
		i++;
	}
}

void Set_wsequence(int *sequence)
{
	float temp;
	int seq;
	float weight[C_MAX];
	for(int k=1;k<=Col_numb;k++)
	{
		sequence[k]=0;
		weight[k]=coil_inf[k].m_nWeight;
	}
	int i=1;
	while(1)
	{
		if(i>Col_numb)
			break;

		temp=weight[1];
		seq=1;
		for(int j=1;j<=Col_numb;j++)
		{
			if(weight[j]>temp)
			{
				temp=weight[j];
				seq=j;
			}
		}
		weight[seq]=0;
		sequence[i]=seq;
		i++;
	}
}


void release_workspace(void)
{

	//	free(solution);
	NEW_COL_NUM=0;
	con_m=0;
	var_n=0;
	if(pit_pool.GetCount()!=0)
	    pit_pool.RemoveAll();
	if(pitlist.GetCount()!=0)
	    pitlist.RemoveAll();
	if(CurrentCutList.GetCount()!=0)
		CurrentCutList.RemoveAll();

	if(result->Col!=NULL)
		delete []result->Col;
	if(result->CutList.GetCount()!=0)
		result->CutList.RemoveAll();
	delete result;

	if(soln_of_each_node.GetCount()!=0)
		soln_of_each_node.RemoveAll();

	while (1)
	{
		pop(&nownode);
		if(!nownode)
			break;

		if(nownode->Col!=NULL)
			delete []nownode->Col;
		if(nownode->CutList.GetCount()!=0)
			nownode->CutList.RemoveAll();
		free (nownode);
		nownode=NULL;
	}


	for(int i=1;i<=Col_numb;i++)
	{
		for(int j=1;j<=Col_numb;j++)
		{
			C2[i][j]=0;
		}
		for(int k=1;k<=Fur_numb;k++)
		{
			C1[i][k]=0;
		}
	}
	for(int i=0;i<=Col_numb;i++)
	{
		free(dkk[i]);
		dkk[i]=NULL;
	}
	free(dkk);
	dkk=NULL;

	for(int i=0;i<=Col_numb;i++)
	{
		free(eij[i]);
		eij[i]=NULL;
	}
	free(eij);
	eij=NULL;
}


void ekk_release( )
{

	int rtcod;
	// Delete the model               
	//rtcod = ekk_deleteModel(Linear_Mode);
	//// OSL Destructor                 
	//ekk_endContext(env);
}

double GetValue(CFurnace tempfur)
{
	double cost=0;
	POSITION pos=tempfur.ingot.GetHeadPosition();
	while (pos)
	{
		COb ing=tempfur.ingot.GetAt(pos);
		cost+=ing.m_nWeight+ing.PRI-C2[ing.index][tempfur.center]-C1[ing.index][tempfur.index];
		tempfur.ingot.GetNext(pos);
	}

	return cost;
}

void branch_and_price()
{
	int rtn;
	int nowlevel;

	struct n *first_node;  //The first node in the branch tree
	first_node=NULL;
	nownode =NULL;
	GenerateNode(first_node);


	while(1)
	{
		pop(&nownode);
		if(!nownode)
			return;
		if(nownode->solu_type==1)
		{
			if(nownode->Col!=NULL)
			    delete []nownode->Col;
			if(nownode->CutList.GetCount()!=0)
			    nownode->CutList.RemoveAll();
			free (nownode);
			nownode=NULL;
			sort(&nodelist);
			continue;
		}

		if(nownode->lower_bound-result->lower_bound>0.5) //if lowerb<Comp_min,there is the   
		{                                  //possibility to branch
			if((Round((nownode->lower_bound - result->lower_bound)/result->lower_bound,4)<(double)Permit_Gap) &&(nownode->level>Col_numb)) 
			{
				//but if the gap between best solution and this lowound at this node is 
				//very minimal,we can cut this node
				sort(&nodelist);
				if(nownode->Col!=NULL)
				    delete []nownode->Col;
				if(nownode->CutList.GetCount()!=0)
				    nownode->CutList.RemoveAll();
				free (nownode);
			}			
			else
			{
				rtn=nownode->solu_type;   //judge the type of the current solution
				if(rtn==0)    //the solution is integer, stop branching
				{             //save the current solution, backtracking
					if(result->Col!=NULL)
					delete []result->Col;
					if(result->CutList.GetCount()!=0)
					    result->CutList.RemoveAll();
					free (result);
					result=nownode;
					sort(&nodelist);
				}		
				if(rtn==-1)  //the optimal solution is fractional, continue branching
				{	
					finish_time=clock();
					if(float(finish_time-begin_time)/1000>=3600)//10800
					{
						if(nownode->Col!=NULL)
							delete []nownode->Col;
						if(nownode->CutList.GetCount()!=0)
							nownode->CutList.RemoveAll();
						free (nownode);
						is_optimal=false;
						break;
					}
					GenerateNode(nownode); //generate the next level nodes
					if(nownode->Col!=NULL)
					    delete []nownode->Col;
					if(nownode->CutList.GetCount()!=0)
					    nownode->CutList.RemoveAll();
					free (nownode);
				}
			}
		}
		else if( nownode->lower_bound==result->lower_bound)
		{			
			rtn=nownode->solu_type;   //judge the type of the current solution
			if(rtn==0)    //the solution is integer, stop branching
			{           //save the current solution, backtracking
				if(result->Col!=NULL)
				    delete []result->Col;
				if(result->CutList.GetCount()!=0)
				    result->CutList.RemoveAll();
				free (result);
				result=nownode;
			}
			if(rtn==-1)
			{
				if(nownode->Col!=NULL)
				    delete []nownode->Col;
				if(nownode->CutList.GetCount()!=0)
				    nownode->CutList.RemoveAll();
				free (nownode);
				nownode=NULL;
			}
			sort(&nodelist);
		}
		else                       //if lowerb>Comp_min,prune the same level nodes
		{                          //continue backtracking
			nowlevel=nownode->level;
			
			do
			{
				if(nownode->Col!=NULL)
				    delete []nownode->Col;
				if(nownode->CutList.GetCount()!=0)
				    nownode->CutList.RemoveAll();
				free (nownode);
				nownode=NULL;
				sort(&nodelist);
				pop(&nownode);
				if(nownode==NULL)
					return;
			}
			while(nownode->lower_bound<result->lower_bound);//nownode->level==nowlevel
			push(&nownode);
		}


	}
}

void   GenerateNode(struct n *nownode)
{
	int rtn;
	int NumbInCol;
	struct n *nextnode, *tempnode, *tempsort;
	nextnode=NULL;
	tempsort=NULL;
	tempnode=NULL;
	int possible;
	if(nownode==NULL)
	{
		remove("./DataSVM/TestDataSVM.txt");
		RowsNumOfTestDataSVM = 0;
		IloEnv env;
		IloModel mod(env);
		IloNumVarArray Lamada(env);//线性规划变量
		IloRangeArray range;
		IloNumArray C(env);
		IloObjective cost;

		if(FNlist.GetCount()!=0)
			FNlist.RemoveAll();

		GenFirstRMP(&env,&mod,&Lamada,&C,&range,&cost);
		int flag=1;
		int b=3;
		int c=2;
		possible=linear_relax(&env,&mod,&Lamada,&C,&range,&cost);
		Variable_Fixing();
		nextnode=new struct n;
		write_node(nextnode,0,0,NEW_COL_NUM,nownode,possible);
		LBinBoot=nextnode->lower_bound;               //Get the LowerBound at the Boot node
		nextnode->level=1;
		nextnode->next=NULL;
		nextnode->is_change=false;
//#pragma region
//		i_LowerBound = 0;
//		printf("\nnextnode->level=%d", nextnode->level);
//		printf("\nnownode->level=%d", nownode->level);
//		printf("\ntempnode->level=%d", tempnode->level);
//		LowerBound[i_LowerBound] = nextnode->lower_bound;
//		i_LowerBound++;
//#pragma endregion


		tempsort=nextnode;
		while(tempsort)
		{
			tempnode=tempsort;
			tempsort=tempsort->next;
			tempnode->next=NULL;
			push(&tempnode);
		}		
		LowerBound[i_LowerBound] = tempnode->lower_bound;
		i_LowerBound++;

#pragma region 选择合适的变量
		/*FILE *OutputDualValue;
		OutputDualValue = fopen("./OutputDualValue.txt", "a+");
		int TotalNumOfInequality = Fur_numb + Col_numb + valid_total_num + NumOfValidInequality;
		for (int i = 0; i < TotalNumOfInequality; i++)
		{
			fprintf(OutputDualValue, "%f		", dual_var[i]);
		}
		fprintf(OutputDualValue, "%f		", tempnode->lower_bound);
		fprintf(OutputDualValue, "%d", tempnode->level);
		fprintf(OutputDualValue, "\n");
		fclose(OutputDualValue);*/

#pragma endregion

		Lamada.clear();
		Lamada.end();
		C.clear();
		C.end();
		range.clear();
		range.end();
		mod.removeAllProperties();
		mod.end();
		env.removeAllProperties();
		env.end();

		return;
	}
	else if(nownode)
	{
		for(int i=1;i<=2;i++)
		{
			remove("./DataSVM/TestDataSVM.txt");
			RowsNumOfTestDataSVM = 0;
			IloEnv env;
			IloModel mod(env);
			IloNumVarArray Lamada(env);//线性规划变量
			IloRangeArray range;
			IloNumArray C(env);
			IloObjective cost;
			if(FNlist.GetCount()!=0)
				FNlist.RemoveAll();
			if(nownode->level==4&&nownode->tree==1)
				int a=0;
			NumbInCol=0;
			numincol=0;
			numincol=var_n;
			old_col_num=NEW_COL_NUM;
			ColInParent=new int[nownode->column_num+1];
			update_arc(i,nownode);    
			rtn=GenNewRmp(i,nownode,&env,&mod,&Lamada,&C,&range,&cost);   // update RMP's model;
			printf("step3!");
			NumbInCol=var_n;
			if(rtn==0)
			{
				printf("node level is: %d	tree is %d\n",(nownode->level)+1,i);
				if((nownode->level+1)==17&&i==1)
					int a=0;
				int flag=1;
				int cou=0;
				if((nownode->s==10&&nownode->r==33)||(nownode->s==33&&nownode->r==10))
					int a=0;
				possible=linear_relax(&env,&mod,&Lamada,&C,&range,&cost);
				Variable_Fixing();
				if((nownode->s==10&&nownode->r==33)||(nownode->s==33&&nownode->r==10))
					int a=0;
				nextnode=new n;
				nextnode->level=(nownode->level)+1;
				nextnode->tree=i;
				write_node(nextnode,NumbInCol,old_col_num,NEW_COL_NUM,nownode,possible);
				if(i==2)
					nextnode->is_change=true;
				else
					nextnode->is_change=nownode->is_change;
				
				nextnode->next=tempsort;

				LowerBound[i_LowerBound] = nextnode->lower_bound;
				i_LowerBound++;
				tempsort=nextnode;
				// add new coumn to this node
			}
			else
			{
				printf("\n******************************************************\n");
				printf("\nError occur during call GenNewRmp,return incident error \n");
				exit(0);
			}
#pragma region
			LowerBound[i_LowerBound] = nextnode->lower_bound;
			i_LowerBound++;
#pragma endregion
#pragma region 选择合适的变量
			//FILE *OutputDualValue;
			//OutputDualValue = fopen("./OutputDualValue.txt", "a+");
			////fprintf(OutputDualValue, "%d\t%d\t", Fur_numb, Col_numb);
			//int TotalNumOfInequality = Fur_numb + Col_numb + valid_total_num + NumOfValidInequality;
			//for (int i = 0; i < TotalNumOfInequality; i++)
			//{
			//	for (int j = TotalNumOfInequality - 1; j>i; j--)
			//	{
			//		if (fabs(dual_var[j]) > fabs(dual_var[j - 1]))
			//		{
			//			int temp_dual_var = dual_var[j];
			//			dual_var[j] = dual_var[j - 1];
			//			dual_var[j - 1] = temp_dual_var;
			//		}
			//	}
			//}
			//for (int i = 0; i < 44; i++)
			//{
			//	fprintf(OutputDualValue, "%f\t", dual_var[i]);
			//	printf("%f\t", dual_var[i]);
			//}
			//int rs_Binary[16] = { 0 };
			//int i_r_Binary = 0;
			//while (nownode->r!=0)
			//{
			//	rs_Binary[i_r_Binary] = nownode->r % 2;
			//	i_r_Binary = i_r_Binary + 1;
			//	nownode->r = nownode->r / 2;
			//}
			//int i_s_Binary = 8;
			//while (nownode->s != 0)
			//{
			//	rs_Binary[i_s_Binary] = nownode->s % 2;
			//	i_s_Binary = i_s_Binary + 1;
			//	nownode->s = nownode->s / 2;
			//}
			////printf("nownnode->lower_bound=%lf, result->lower_bound", nownode->lower_bound, result->lower_bound);
			//double ImproveRateOfOptimalValue = (nextnode->lower_bound - result->lower_bound) / result->lower_bound;
			//fprintf(OutputDualValue, "%lf\t", ImproveRateOfOptimalValue);
			//printf("%lf\t", ImproveRateOfOptimalValue);
			//fprintf(OutputDualValue, "%d\t", nextnode->level);
			//printf("%d\t", nextnode->level);
			//for (int i = 0; i < 16; i++)
			//{
			//	fprintf(OutputDualValue, "%d\t", rs_Binary[i]);
			//	printf("%d\t", rs_Binary[i]);
			//}
			//fprintf(OutputDualValue, "%f\t", nextnode->lower_bound);
			//printf( "%lf", nextnode->lower_bound);
			//fprintf(OutputDualValue, "\n");
			//fclose(OutputDualValue);

#pragma endregion

//			free(ColInParent);
			delete []ColInParent;
//			nodenum++;
			Lamada.clear();
			Lamada.end();
			C.clear();
			C.end();
			range.clear();
			range.end();
			mod.removeAllProperties();
			mod.end();
			env.removeAllProperties();
			env.end();
		}
		if(tempsort->solu_type==-1 &&tempsort->next->solu_type==-1)
		{
			sort(&tempsort);
			while(tempsort)
			{
				tempnode=tempsort;
				tempsort=tempsort->next;
				tempnode->next=NULL;
				push(&tempnode);
			}
			return;
		}
		if(tempsort->solu_type==-1 &&tempsort->next->solu_type==0)
		{
			while(tempsort)
			{
				tempnode=tempsort;
				tempsort=tempsort->next;
				tempnode->next=NULL;
				push(&tempnode);
			}	   
			return;
		}
		if(tempsort->solu_type==0 &&tempsort->next->solu_type==-1)
		{
			tempnode=tempsort;
			tempsort=tempsort->next;
			tempsort->next=tempnode;
			tempnode->next=NULL;
			while(tempsort)
			{
				tempnode=tempsort;
				tempsort=tempsort->next;
				tempnode->next=NULL;
				push(&tempnode);
			}
			return;
		}
		if(tempsort->solu_type==0 &&tempsort->next->solu_type==0)
		{
			sort(&tempsort);
			while(tempsort)
			{
				tempnode=tempsort;
				tempsort=tempsort->next;
				tempnode->next=NULL;
				push(&tempnode);
			}
			return;
		}
		else
		{
			while(tempsort)
			{
				tempnode=tempsort;
				tempsort=tempsort->next;
				tempnode->next=NULL;
				push(&tempnode);
			}	   
			return;

		}
	}

}

void push(struct n **node)
/*---------------------------------------------------------------------------*/
/* <    >                                                                    */
/*****************************************************************************/
{
	if(*node)
	{
		(*node)->next=nodelist;
		nodelist=*node;
	}
}

/*****************************************************************************/
/* [      ]                                                                  */
/* [      ]    void pop()                                                    */
/* [      ]                                                                  */
/*---------------------------------------------------------------------------*/
void pop(struct n **node)
/*---------------------------------------------------------------------------*/
/* <    >                                                                    */
/*****************************************************************************/
{
	*node=nodelist;
	if(nodelist)
		nodelist=nodelist->next;
	if(*node)
		(*node)->next=NULL;
}

int linear_relax(IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost)
{
	int     rtn;
	int     rt;
	int     i=0;
	if(pit_pool.GetCount()!=0)
		pit_pool.RemoveAll();

	median_from=0;
	median_to=0;
	int add_value=Col_numb/2;
	old_col_num=NEW_COL_NUM;
	int NumbInCol=0;
	numincol=0;
	numincol=var_n;
	
	NumbInCol=var_n;
	rt=solve_linear(env,mod,*range,*Lamada,*C,cost);//什么情况，返回值能为-1？//??
	if(nownode==NULL)
	{
		NumbInCol=0;
		old_col_num=0;
	}
   
	if (rt==-1)
	{   
		printf("error occour during solve linear program\n");
		return -1;
	}	    
	int non_number=0;
	while(1)
	{	
		finish_time=clock();
		if(float(finish_time-begin_time)/1000>=3600)//10800
		{
			is_optimal=false;
			break;
		}
		if (pit_pool.GetCount()==0)
		{
			median_from=median_to+1;
			if(median_from>=Col_numb)
				median_from=1;
			median_to=median_from+add_value;
			if (median_to>=Col_numb)
			{
				median_to=Col_numb;
			}
			finish_time=clock();
		    if(float(finish_time-begin_time)/1000>=3600)//10800
		    {
			    is_optimal=false;
			    break;
		    }

			rtn=DynamicPro(*dual_var,dual_var);
			

			Sub_problem_Num++;
			
			if(rtn==0)
				non_number++;
			else
				non_number=0;
			if(non_number==2)     //the current solution is optimal
				break;   
			//terminate the pricing algorithm
			if(rt==-999)
				break;

		}

		up_linear_mode(env,mod,Lamada,C,range,cost);

#pragma region
		FILE *OutputDualValue;
		OutputDualValue = fopen("./DataSVM/TestDataSVM.txt", "a+");
		//fprintf(OutputDualValue, "%d\t%d\t", Fur_numb, Col_numb);
		int TotalNumOfInequality = Fur_numb + Col_numb + valid_total_num + NumOfValidInequality;
		for (int i = 0; i < TotalNumOfInequality; i++)
		{
			for (int j = TotalNumOfInequality - 1; j>i; j--)
			{
				if (fabs(dual_var[j]) > fabs(dual_var[j - 1]))
				{
					int temp_dual_var = dual_var[j];
					dual_var[j] = dual_var[j - 1];
					dual_var[j - 1] = temp_dual_var;
				}
			}
		}
		for (int i = 0; i < 44; i++)
		{
			fprintf(OutputDualValue, "%lf\t", dual_var[i]);
			printf("%lf\t", dual_var[i]);
		}
		int rs_Binary[16] = { 0 };
		int i_r_Binary = 0;
		int r_nownode = nownode->r;
		int s_nownode = nownode->s;
		while (r_nownode!=0)
		{
			rs_Binary[i_r_Binary] = r_nownode % 2;
			i_r_Binary = i_r_Binary + 1;
			r_nownode = r_nownode / 2;
		}
		int i_s_Binary = 8;
		while (s_nownode != 0)
		{
			rs_Binary[i_s_Binary] = s_nownode % 2;
			i_s_Binary = i_s_Binary + 1;
			s_nownode = s_nownode / 2;
		}
		//printf("nownnode->lower_bound=%lf, result->lower_bound", nownode->lower_bound, result->lower_bound);
		//double ImproveRateOfOptimalValue = (result->lower_bound - last_lower_bound) / result->lower_bound;
		//fprintf(OutputDualValue, "%lf\t", ImproveRateOfOptimalValue);
		//printf("%lf\t", ImproveRateOfOptimalValue);
		//fprintf(OutputDualValue, "%d\t", nextnode->level);
		//printf("%d\t", nextnode->level);
		for (int i = 0; i < 16; i++)
		{
			fprintf(OutputDualValue, "%d\t", rs_Binary[i]);
			printf("%d\t", rs_Binary[i]);
		}
		fprintf(OutputDualValue, "%f\t", result->lower_bound);
		printf( "%lf", result->lower_bound);
		fprintf(OutputDualValue, "\n");
		fclose(OutputDualValue);
		RowsNumOfTestDataSVM++;
#pragma endregion
		/*here use the ekk api solve the linear program*/
		if(rtn!=0)
		{
			bool flag=true;
			while(flag)
			{
				rt=solve_linear(env,mod,*range,*Lamada,*C,cost);
				flag=GenerateCutss(NumbInCol,old_col_num,NEW_COL_NUM,rt,env,mod,Lamada,C,range,cost,3, 2);
				finish_time=clock();
				if(float(finish_time-begin_time)/1000>=3600)//10800
				{
					is_optimal=false;
					break;
				}
			}
		}
		else
			rt=solve_linear(env,mod,*range,*Lamada,*C,cost);
//		rt=solve_linear(env,mod,*range,*Lamada,*C);
		if (rt==-1)
		{  printf("error occour during solve linear program\n");
		return -1;
		// 			getchar();
		// 		    exit(0);
		}
	}
	
	return 0;
}

void write_node(struct n *node,int NumbInCol,int old_col_num,int NEW_COL_NUM,struct n *nownode,int possible)
{
	int       solu_numb,i;
	int       *Col;
	double   Flow[C_MAX+1];
	nodenum++;
	node->NodeNum=nodenum;

	for(int i=1;i<=Col_numb;i++)
	{
		if(nownode)
			node->model_modify[i]=model_modify[i]; //记录需要在模型中改变的行
		else
			node->model_modify[i]=0;
		for(int j=1;j<=Col_numb;j++)
			node->arc_maximal[i][j]=arc_matrix[i][j];
		node->median[i]=median[i];
		for(int j=1;j<=Fur_numb;j++)
		{
			node->fur_coil[j][i]=fur_coil[j][i];
			node->Ykk[i][j]=Ykk[i][j];
		}
	}
	if(nownode)
		node->ParentNum=nownode->NodeNum;
	//for(i=1;i<=ingot_num_knap;i++)
	//	node->ingot_knap[i]=ingot_knap[i];
	node->ingot_num_knap=ingot_num_knap;

	Col =new int[NEW_COL_NUM-old_col_num+NumbInCol+1];
	for(int i=1;i<=NumbInCol;i++)
		Col[i]=ColInParent[i];

	for(i=NumbInCol+1;i<=NEW_COL_NUM - old_col_num +NumbInCol;i++)//??
		Col[i]=old_col_num+ i - NumbInCol;

	node->Col=Col;

	POSITION cut_pos=CurrentCutList.GetHeadPosition(); //记录各个节点中加入的割平面
	while (cut_pos)
	{
		struct AddCut st=CurrentCutList.GetAt(cut_pos);
		node->CutList.AddTail(st);
		CurrentCutList.RemoveAt(cut_pos);
		cut_pos=CurrentCutList.GetHeadPosition();
	}

	/*Here We Will Get The Solution of Linear Model*/ 
	int coil_in[C_MAX];
	if(possible!=-1)
	{
		FILE *fp1;
		if((fp1=fopen("middle_result.txt","wt+"))==NULL)
		{
			printf("Cannot open file strike any key exit 7!");
			exit(1);
		}
		solu_numb=0; 
		int n=col_num;
		POSITION sol=soln_of_each_node.GetHeadPosition();
		for(int m=0;m<n;m++)
		{
			double solu=soln_of_each_node.GetAt(sol);
			printf("\n");
			printf("solu = %f		",solu);
			if(Round(solu,6)>(double)EPS)  //only record the number of current optimal solutions greater than ZERO
			{
				solu_numb++;
				node->solution [solu_numb]=solu;
				node->soluindex[solu_numb]=m+1;

				CFurnace Now_pit;
				Now_pit=Getpit(Col[node->soluindex[solu_numb]]);
				//int lenth=Now_pit.ingot.GetCount();
				fprintf(fp1,"%d	%lf\n",Now_pit.index,Now_pit.cost);
				printf("%d	%lf\n", Now_pit.index, Now_pit.cost);
				POSITION pos=Now_pit.ingot.GetHeadPosition();
				while (pos)
				{
					COb ob=Now_pit.ingot.GetAt(pos);
					fprintf(fp1,"%s	%s	%d	%lf\n",ob.m_strCoilNum,ob.m_strCurve,ob.m_nWidth,ob.m_nWeight);
					printf("%s	%s	%d	%lf   %d\n", ob.m_strCoilNum, ob.m_strCurve, ob.m_nWidth, ob.m_nWeight,ob.index);
					if(coil_in[ob.index]!=1)
						coil_in[ob.index]=1;
					Now_pit.ingot.GetNext(pos);
				}
				fprintf(fp1,"%f  ",node->solution[solu_numb]);
				printf("%f  ", node->solution[solu_numb]);
				fprintf(fp1,"\n\n");
				printf("\n");
			}
			soln_of_each_node.GetNext(sol);
		}

		fprintf(fp1,"\n");
		printf("\n");
		for(int i=1;i<=Col_numb;i++)
		{
			if (coil_in[i] == 1)
				printf("%d\n", i);
				fprintf(fp1,"%d\n",i);
		}
		fprintf(fp1,"\n");
		printf("\n");

		//soln_of_each_node.RemoveAll();
		node->lower_bound=Round(opt_soln_of_eachnode,0);
		node->column_num=n;
		node->solu_numb=solu_numb;    // record the solution number
		fclose(fp1);
		
	}

	if(possible==-1)//||flag
	{
		node->solu_type=1;
		node->r=0;
		node->s=0;
	}
	else if(solu_numb==(NHB+NHS+HHB+HHS))
	{  // record the solution type 	
		node->solu_type =0;
		node->r=0;
		node->s=0;
		//printf("opt in this node is %lf",opt);
		//getchar();
	}
	else
	{
		node->solu_type =-1;
		GetBranchArc( node);
	}
	
	// 	node->model_modify[node->r]=1;
	// 	node->model_modify[node->s]=1;
}
void    Network_arc( void)
{
	for(int i=1;i<=Col_numb;i++)
	{
		for(int j=1;j<=Col_numb;j++)
		{
			arc_matrix[i][j]=1;
		}
		ingot_knap[i]=coil_inf[i];
		branch_coil[i]=coil_inf[i];
		model_modify[i]=0;
		median[i]=1;
		for(int j=1;j<=Fur_numb;j++)
		{
			fur_coil[j][i]=1;
			Ykk[i][j]=1;
		}
	}
	ingot_num_knap=Col_numb;
//	get_collect();
	for(int k=1;k<=Fur_numb;k++)
		pit_knap[k]=batch[k];
}

void update_arc(int node_num, struct n *node)
{
	int i,j;
	int r,s;
	r=node->r;
	s=node->s;
	int pos1,pos2,pos1_in,pos2_in;

	if((s==23)||(r==23))
		int a=0;

	/*更新割平面*/
	CurrentCutList.RemoveAll();
	POSITION pos=node->CutList.GetHeadPosition();
	while (pos)
	{
		struct AddCut st=node->CutList.GetAt(pos);
		CurrentCutList.AddTail(st);
		node->CutList.GetNext(pos);
	}

	for(j=1;j<=Col_numb;j++)
	{
		for (i = 1; i <= Col_numb; i++)
		{
			arc_matrix[i][j] = node->arc_maximal[i][j];
			if (arc_matrix[i][j] != 0)
			{
				printf("		arc_martix[%d][%d]:%f", i,j,arc_matrix[i][j]);
			}
		}
		median[j]=node->median[j];
		printf("\n");
		for(int k=1;k<=Fur_numb;k++)
		{
			fur_coil[k][j]=node->fur_coil[k][j];
			printf("fur_coil[%d][%d]:%f		", k, j, fur_coil[k][j]);
			Ykk[j][k]=node->Ykk[j][k];
			printf("Ykk[%d][%d]:%f		", j, k, Ykk[j][k]);
		}
		model_modify[j]=node->model_modify[j];
		printf("\nmodel_modify[%d]:%f		", model_modify[j]);

	}
	ingot_num_knap=node->ingot_num_knap;
	if(r==20&&s==21)
		int a=0;
	if(node_num==1)
	{
		printf("		node->is_median=%f", node->is_median);
		if(node->is_median!=-1)
		{
			median[node->is_median]=4;
			model_modify[node->is_median]=-1;
			if(coil_inf[node->is_median].collect_num!=0)//如果好的都不在方案里，那么不好的肯定不在方案里
			{
				for(int k=1;k<=coil_inf[node->is_median].collect_num;k++)
				{
					median[coil_inf[node->is_median].ingotlist[k]]=4;
					model_modify[coil_inf[node->is_median].ingotlist[k]]=-1;
				}
			}
		}
		else if(node->fur_select!=-1)
		{
			printf("		%d", node->fur_select);
			fur_coil[node->fur_select][node->coil_select]=4;  //板卷一定不在该种类型炉子中
		}
		else
		{
			arc_matrix[r][s]=4;     //3：一定在一起，4：一定不在一起
			arc_matrix[s][r]=4;
			printf("\n		arc_matrix[%d][%d]:%f", r, s, arc_matrix[r][s]);
			bool is_in=false;
			printf("\n branch_coil[%d].collect_num=%d", r, branch_coil[r].collect_num);
			for(int k=1;k<branch_coil[r].collect_num;k++)
			{
				printf("\n branch_coil[%d].ingotlist[%d]=%f", r, k, branch_coil[r].ingotlist[k]);
				if(s==branch_coil[r].ingotlist[k])
				{
					is_in=true;
					break;
				}
			}
			for(int k=1;k<=Col_numb;k++)
			{
				printf("   arc_matrix[%d][%d]=%f", r, k, arc_matrix[r][k]);
				printf("   arc_matrix[%d][%d]=%f", s, k, arc_matrix[s][k]);
				if (arc_matrix[r][k] == 3)
				{
					arc_matrix[s][k]=4;
					arc_matrix[k][s]=4;
				}
				else if(arc_matrix[s][k]==3)
				{
					arc_matrix[r][k]=4;
					arc_matrix[k][r]=4;
				}
			}
		}
		
	}
	else if(node_num==2)
	{
		printf("		node->is_median=%f", node->is_median);
		if(node->is_median!=-1)
		{
			median[node->is_median]=3;
			model_modify[node->is_median]=1;

			if(ingot_knap[node->is_median].collect_num!=0)//如果不好的都在方案里，那么好的肯定在方案里
			{
				for(int k=1;k<=ingot_knap[node->is_median].collect_num;k++)
				{
					median[ingot_knap[node->is_median].ingotlist[k]]=3;
					model_modify[ingot_knap[node->is_median].ingotlist[k]]=1;
				}
			}
		}
		else if(node->fur_select!=-1)
		{
			printf("		%d", node->fur_select);
			model_modify[node->coil_select]=1;
			for(j=1;j<=Fur_numb;j++)
			{
				if(j!=node->fur_select)
					fur_coil[j][node->coil_select]=4;  //板卷一定不能放在别的类型的炉子里
				else 
					fur_coil[j][node->coil_select]=3;
			}
		}
		else
		{
			printf("		model_modify[%d]=%f", r, model_modify[r]);
			model_modify[r]=1;
			printf("		model_modify[%d]=%f", s, model_modify[s]);
			model_modify[s]=1;
			printf("		ingot_knap[%d].collect_num=%d", r, ingot_knap[r].collect_num);
			if(ingot_knap[r].collect_num!=0)//如果不好的都在方案里，那么好的肯定在方案里
			{
				printf("		ingot_knap[%d].collect_num=%d", r, ingot_knap[r].collect_num);
				for(int k=1;k<=ingot_knap[r].collect_num;k++)
				{
					median[ingot_knap[r].ingotlist[k]]=3;
					model_modify[ingot_knap[r].ingotlist[k]]=1;
				}
			}
			if(ingot_knap[s].collect_num!=0)//如果不好的都在方案里，那么好的肯定在方案里
			{
				for(int k=1;k<=ingot_knap[s].collect_num;k++)
				{
					median[ingot_knap[s].ingotlist[k]]=3;
					model_modify[ingot_knap[s].ingotlist[k]]=1;
				}
			}
			arc_matrix[r][s]=3;
			arc_matrix[s][r]=3;
			for(int k=1;k<=Col_numb;k++)
			{
				if(arc_matrix[r][k]==4)
				{
					arc_matrix[s][k]=4;
					arc_matrix[k][s]=4;
				}
				if(arc_matrix[s][k]==4)
				{
					arc_matrix[r][k]=4;
					arc_matrix[k][r]=4;
				}
			}
			if((s==21&&r==28)||(s==28&&r==21))
				int a=0;
			
		}
		
	}
	else
	{	
		printf("Para Err in up_arc\n");
		exit(0);
	}
}

void sort(struct n **list)
/*---------------------------------------------------------------------------*/
/* <    >                                                                    */
/*****************************************************************************/
{
	struct n *node1, *min, *before, *list0;
	before=NULL;
	list0=NULL;
	node1=*list;
	min=*list;
	while(*list)
	{
		while(node1)
		{
			if(node1->next)
			{
				if(node1->next->lower_bound>min->lower_bound)//<
				{
					min=node1->next;
					before=node1;
				}
				else if(node1->next->lower_bound==min->lower_bound)
				{
					if(node1->next->solu_type==0)
					{
						min=node1->next;
						before=node1;
					}
				}               
				/*Here add my idear sort by integer*/
			}
			node1=node1->next;
		}
		if(before==NULL)
		{
			*list=(*list)->next;
			min->next=list0;
			list0=min;
		}
		else
		{
			before->next=min->next;
			min->next=list0;
			list0=min;
		}
		before=NULL;
		node1=*list;
		min=*list;
	}
	*list=list0;
}

void GetBranchArc( struct n *node)
{
#pragma region


	int      i,j;
	int      l;
	int      lenth;
	int      solution_number;
//	double   Flow[C_MAX+1][F_MAX+1];
	float flow[C_MAX+1][C_MAX+1];
	float is_median[C_MAX+1];//各个板卷中心分量
	double   MostFraArc;
	double   value;
	int      r,s;
	float f_c[F_MAX+1][C_MAX+1];

	for(i=1;i<Col_numb+1;i++)
	{
		for (j=1;j<Col_numb+1;j++)
		{
			flow[i][j] = 0;
		}
		for(int k=1;k<=Fur_numb;k++)
			f_c[k][i]=0;
		is_median[i]=0;
	}

	solution_number=node->solu_numb;

	POSITION pos;
	CFurnace Now_pit;
	pos=pitlist.GetHeadPosition();
	for (l=1;l<=solution_number;l++)
	{
		//		CPit Now_pit;
		//		Now_pit=Getpit(node->Col[node->soluindex[l]]);
		Now_pit=pitlist.GetAt(pos);
		while(Now_pit.No!=node->Col[node->soluindex[l]])
		{
			pitlist.GetNext(pos);
			Now_pit=pitlist.GetAt(pos);
		}
		if(Now_pit.ingot.GetCount()==0)
		{
			printf("Error When Getcol in GetBranchArc \n");
			for(i=1;i<=node->column_num;i++)
				printf("%d  ",node->Col[i]);
			getchar( );
			exit(0);
		}
		POSITION pos=Now_pit.ingot.GetHeadPosition(); 
		
		printf("furnace:%d	",Now_pit.index);
//		fprintf(A,"\nf:%d	s:%f	c:%f	cen:%d	",Now_pit.index,node->solution[l],Now_pit.cost,Now_pit.center);
		while(pos)
		{
			COb ings=Now_pit.ingot.GetAt(pos);
			f_c[Now_pit.index][ings.index]+=node->solution[l];
			is_median[ings.index]+=node->solution[l];
			POSITION pos_next=pos;
			Now_pit.ingot.GetNext(pos_next);
//			fprintf(A,"%d	",ings.index);
			printf("%d	",ings.index);
			while (pos_next)
			{
				COb ingt=Now_pit.ingot.GetAt(pos_next);
				printf("\n%d_%d:	", ings.index, ingt.index);
				flow[ings.index][ingt.index]+=node->solution[l];
				printf("%f	", flow[ings.index][ingt.index]);
				flow[ingt.index][ings.index]+=node->solution[l];
				printf("%f	", flow[ingt.index][ings.index]);
				
				Now_pit.ingot.GetNext(pos_next);
			}
            //timet[ings.index]++;
//			times[ings.index][Now_pit.index]++;			
			Now_pit.ingot.GetNext(pos);
		}
//		printf("\n");
		Now_pit.ingot.RemoveAll();
	}

	
#pragma region
	TestData[10000][61] = { 0 };
	FILE * InputData1;
	if ((InputData1 = fopen("./DataSVM/TestDataSVM.txt", "rt+")) == NULL)
	{
		printf("rs数据未打开!\n");
		exit(0);
	}
	for (int i1 = 0; i1 < RowsNumOfTestDataSVM; i1++)
	{
		for (int j1 = 0; j1 < 61; j1++)
		{
			fscanf(InputData1, "%lf\t", &TestData[i1][j1]);
		}

	}
	fclose(InputData1);

	int counter1 = 0;
	for (int i3 = 1; i3 < 301; i3++)
	{
		for (int j3 = 1; j3 < 301; j3++)
		{
			if (flow[i3][j3]>0 && flow[i3][j3] < 1)
			{
				int i3_0 = i3;
				int j3_0 = j3;
				printf("i=%d  j=%d   counter1=%d\n", i3_0, j3_0, counter1);
				int counter2 = 0;
				while (j3_0 != 0)
				{
					TestData[counter1][59 - counter2] = j3_0 % 2;
					j3_0 = j3_0 / 2;
					counter2++;
				}
				while (counter2 < 8)
				{
					TestData[counter1][59 - counter2] = 0;
					counter2++;
				}
				while (i3_0 != 0)
				{
					TestData[counter1][59 - counter2] = i3_0 % 2;
					i3_0 = i3_0 / 2;
					counter2++;
				}
				while (counter2 < 8)
				{
					TestData[counter1][59 - counter2] = 0;
					counter2++;
				}
				counter1++;
				
			}
		}
	}
	int i_select = RowsNumOfTestDataSVM - 1;
	for (int i2 = 0; i2 < counter1; i2++)
	{
		for (int j2 = 0; j2 < 44; j2++)
		{
			TestData[i2][j2] = TestData[i_select][j2] / TestData[i_select][0];
		}
	}
	FILE *OutputTestData;
	if ((OutputTestData = fopen("./DataSVM/TestDataSVMForTest.txt", "wt+")) == NULL)
	{
		printf("打开失败\n");
		exit(0);
	}
	printf("counter1=%d\t", counter1);
	for (int i4 = 0; i4 < counter1; i4++)
	{
		for (int j4 = 0; j4 < 44; j4++)
		{
			fprintf(OutputTestData, "%lf\t", TestData[i4][j4]);
		}

		for (int j4 = 44; j4 < 60; j4++)
		{
			fprintf(OutputTestData, "%d\t", TestData[i4][j4]);
		}
		fprintf(OutputTestData, "%lf\n", TestData[i4][60]);
	}
	fclose(OutputTestData);
	int featureNum = 60;
	CLibSVM libSVR;
	libSVR.initialize("./DataSVM/TestDataSVMForTest.txt", featureNum);

	libSVR.predictTest();
	double MaxValue = 0;
	int CountOfPredictValueTest = 0;//计数，遍历了几个预测值
	int IndexInPredictValueTest = 0;
	vector<double>::iterator t;
	for (t = libSVR.predictValueTest.begin(); t != libSVR.predictValueTest.end(); t++)
	{
		if ((*t > MaxValue) && (CountOfPredictValueTest<libSVR.predictValueTest.size()-1))
		{
			MaxValue = *t;
			printf("MaxValue=%lf\t", MaxValue);
			IndexInPredictValueTest = CountOfPredictValueTest;
			printf("IndexInPredictValueTest=%d\t", IndexInPredictValueTest);
		}
		CountOfPredictValueTest++;
	}
	svm_node* RowItemOfFeature = new svm_node[featureNum];
	list<svm_node*>::iterator it_RowItem;
	r = 0;
	s = 0;
	int CountOfxListTest = 0;
	int rs_Binary[16] = { 0 };
	for (int i5 = 44; i5 < 60; i5++)
	{
		rs_Binary[i5 - 44] = TestData[IndexInPredictValueTest][i5];
		printf("rs_Binary[%d]=%d\t", i5 - 44, rs_Binary[i5 - 44]);
	}
	/*for (it_RowItem = libSVR.xListTest.begin(); it_RowItem != libSVR.xListTest.end(); it_RowItem++)
	{
		for (int i5_0 = 44; i5_0 < 60; i5_0++)
		{
			rs_Binary[i5_0 - 44] = (int)(*it_RowItem)[i5_0].value + 0.1;
			printf("rs_Binary[%d]=%d\t", i5_0 - 44, (*it_RowItem)[i5_0].value);
		}
		CountOfxListTest++;
		if (CountOfxListTest == IndexInPredictValueTest)
		{
			for (int i5 = 44; i5 < 60; i5++)
			{
				rs_Binary[i5 - 44] = (int)(*it_RowItem)[i5].value+0.1;
				printf("CountOfxListTest=%d\t", CountOfxListTest);
				printf("rs_Binary[%d]=%d\t", i5 - 44, (*it_RowItem)[i5].value);
			}
		}
	}*/
	for (int i6 = 0; i6 < 8; i6++)
	{
		if (rs_Binary[i6] == 1)
		{
			
			r += pow(2, (8 - i6));
			printf("r=%d\t", r);
		}
	}
	for (int i7 = 8; i7 < 16; i7++)
	{
		if (rs_Binary[i7] == 1)
		{
			
			s += pow(2, (16 - i7));
			printf("s=%d\t", s);
		}
	}
	node->r = r;
	node->s = s;
	printf("node->r=%d\t", node->r);
	printf("node->s=%d\t", node->s);
#pragma endregion


	/*Check if there exist such arc[i][j]=0 and Flow[i][j]>0*/    //??
	float most = 0.5;
	node->fur_select=-1;
	node->coil_select=-1;
	node->is_median=-1;
//#pragma region 获得两个板卷Xij中更靠近1，并且C2[i][j]更小的值
//	most=0;
//	int i_indexs=0;
//	int i_indext=0;
//	float via_most=MAX;
//	for(i=1;i<=Col_numb;i++)
//	{
//		for (j=i+1;j<=Col_numb;j++)
//		{
//			printf("\nflow[%d][%d]=%f", i, j, flow[i][j]);
//			value=Round(flow[i][j],5);
//			value=fabs(value);
//			if(Round(value,5)>Round(most,5)&&Round(value,5)!=0&&Round(value,5)!=1)
//			{
//				most=value;
//				i_indexs=i;
//				i_indext=j;
//				via_most=C2[i][j];
//			}
//			else if(Round(value,5)==Round(most,5)&&Round(value,5)!=0&&Round(value,5)!=1)
//			{
//				if(C2[i][j]<via_most)
//				{
//					most=value;
//					i_indexs=i;
//					i_indext=j;
//					via_most=C2[i][j];
//				}
//			}
//		}
//	}
//	r=i_indexs;
//	s=i_indext;
//	node->r=r;
//	node->s=s;
//#pragma endregion
	if(r>0&&s>0)
		int a=0;
	printf("level:%d	tree:%d	median:%d	fur:%d	coil:%d	r:%d	s%d:	%lf	%d	%d\n",node->level,node->tree,node->is_median,node->fur_select,node->coil_select,r,s,node->lower_bound,node->NodeNum,node->ParentNum);
	if(r==60&&s==62)
		int a=0;
	printf("median=%d,fur=%d,coil=%d, r=%d, s=%d\n",node->is_median,node->fur_select,node->coil_select,r,s);
	if(node->is_median!=-1)
	{
		printf("median[%d]=%f",node->is_median,is_median[node->is_median]);
		if(node->median[node->is_median]==3||node->median[node->is_median]==4)
			getchar();
	}
	else if(node->fur_select!=-1)
	{
		printf("fur_coil[%d][%d]=%d",node->fur_select,node->coil_select,node->fur_coil[node->fur_select][node->coil_select]);
		if(node->fur_coil[node->fur_select][node->coil_select]==3||node->fur_coil[node->fur_select][node->coil_select]==4)
		{
			printf("Error Occur fur=%d,coi=%d",node->fur_select,node->coil_select);
			printf("level number is : %d\n",node->level);
			getchar( );
		}
	}
	else
	{
		printf("arc_matrix[%d][%d]= %d",r,s,node->arc_maximal[r][s]);
		if(node->arc_maximal[r][s]==3||node->arc_maximal[r][s]==4)
		{
			int a=0;
			printf("Error Occur r=%d, s=%d",r,s);
			printf("level number is : %d\n",node->level);
		}
	}
	if(node->fur_select==4&&node->coil_select==86)
		int a=0;
#pragma endregion
}

int GenNewRmp( int tree ,struct n *parentnode,IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost)
{
	
	int     r, s,i;
	int     l;
	int     columnNumb;
	int     *Col;
	int     RN;
	double  objectiveCoefficient;
	int     *TempCol;
	int     rtcod;
	int     rtn;
	int be_median,fur,coi;

	r=parentnode->r;
	s=parentnode->s;
	be_median=parentnode->is_median;
	fur=parentnode->fur_select;
	coi=parentnode->coil_select;

	TempCol = new int[Col_numb+1];

	columnNumb=parentnode->column_num;
	Col=parentnode->Col;

	int **La; //记录模型中A矩阵的系数
	La=(int **)malloc((columnNumb+1)*sizeof(int*));
	for (i=0;i<=columnNumb;i++)
	{
		La[i]=(int *)malloc((Col_numb+2)*sizeof(int));
	}


	var_n=0;
	if(tree==1)
	{
		POSITION pos;
		CFurnace temp;
		pos=pitlist.GetHeadPosition();

		printf("\nparentnode->column_num=%d,", parentnode->column_num);
		for(int i=1;i<=parentnode->column_num;i++)//
		{
			temp=pitlist.GetAt(pos);
			printf("		temp.No=%d  Col[%d]=%d", temp.No, i, Col[i]);
			while(temp.No!=Col[i])
			{
				pitlist.GetNext(pos);
				temp=pitlist.GetAt(pos);
			}
			printf("		temp.No=%d  Col[%d]=%d", temp.No, i, Col[i]);
			printf("		r=%d, s=%d, be_median=%d, fur=%d, coi=%d\n ", r, s, be_median, fur, coi);
			rtn=CheckBranch(temp,r,s,be_median,fur,coi);
			if(rtn==0||rtn==2)
			{
				var_n++;
				ColInParent[var_n]=Col[i];
				POSITION pos=temp.ingot.GetHeadPosition();
				while (pos)
				{
					COb ob=temp.ingot.GetAt(pos);
					printf("		La[%d][%d]=%f", var_n,ob.index,La[var_n][ob.index]);
					La[var_n][ob.index]=1;
					temp.ingot.GetNext(pos);
				}
			}
		}
	}
	else if(tree==2)
	{
		POSITION pos;
		CFurnace temp;
		pos=pitlist.GetHeadPosition();
		printf("\nparentnode->column_num=%d,", parentnode->column_num);
		for(int i=1;i<=parentnode->column_num;i++)
		{ 
			temp=pitlist.GetAt(pos);
			printf("		temp.No=%d  Col[%d]=%d", temp.No, i, Col[i]);
			while(temp.No!=Col[i])
			{
				pitlist.GetNext(pos);
				temp=pitlist.GetAt(pos);
			}
			printf("		temp.No=%d  Col[%d]=%d", temp.No, i, Col[i]);
			rtn=CheckBranch(temp,r,s,be_median,fur,coi);
			if(rtn==1||rtn==2)
			{
				var_n++;
				ColInParent[var_n]=Col[i];
//				fprintf(fp1,"%d	%d	",temp.index,temp.center);
				POSITION pos=temp.ingot.GetHeadPosition();
				while (pos)
				{
					COb ob=temp.ingot.GetAt(pos);
					printf("		La[%d][%d]=%f", var_n, ob.index, La[var_n][ob.index]);
					La[var_n][ob.index]=1;
//					fprintf(fp1,"%d	",ob.index);
					temp.ingot.GetNext(pos);
				}
//				fprintf(fp1,"\n");
				for (int k = 1; k <= temp.lenth; k++)
				{
					La[var_n][temp.coils_in[k]] = 1;
					printf("\nLa[%d][%d]=%d", var_n, temp.coils_in[k], La[var_n][temp.coils_in[k]]);
				}
			}
			//			temp.ingot.RemoveAll();
		}
//		fclose(fp1);
	}

	RN=var_n+1;

	IloNumArray  ConstraintMin(*env);
	IloNumArray  ConstraintMax(*env);

	//IloObjective cost;
	ConstraintMin.add(IloNum(-IloInfinity));
	ConstraintMax.add(IloNum(NHB));//1

	ConstraintMin.add(IloNum(-IloInfinity));
	ConstraintMax.add(IloNum(HHB));

	ConstraintMin.add(IloNum(-IloInfinity));
	ConstraintMax.add(IloNum(NHS));

	ConstraintMin.add(IloNum(-IloInfinity));
	ConstraintMax.add(IloNum(HHS));

#pragma region 添加分支约束，此处分支针对的是板卷， 被分支的板卷置0或1，这里置0或1的方式是给定约束：1<=()<=1。原来主问题模型中第（炉子数+1）行到第（炉子数+板卷数）行约束描述的是一个板卷最多被分到一个炉子里。现在确定了需要分支的板卷也就是子问题模型中的xi时，找到xi所对应的唯一性约束，将不等式改为等式，即表示该板卷一定被选或者一定不被选。
	for(int i=Fur_numb;i<Fur_numb+Col_numb;i++) //钢锭唯一性约束左右端项
	{
		if(model_modify[i-Fur_numb+1]==1)
		{
			ConstraintMin.add(IloNum(1));
			ConstraintMax.add(IloNum(1));
		}
		else if(model_modify[i-Fur_numb+1]==-1)
		{
			ConstraintMin.add(IloNum(0));
			ConstraintMax.add(IloNum(0));
		}
		else
		{
		    ConstraintMin.add(IloNum(-IloInfinity));//-INFINITE-IloInfinity
		    ConstraintMax.add(IloNum(1));
		}
	}
#pragma endregion

	*range= IloRangeArray(*env, ConstraintMin, ConstraintMax);
	mod->add(*range);

	*cost = IloAdd(*mod, IloMaximize(*env));

	POSITION pos;
	//	    CPit temp_pit;
	pos=pitlist.GetHeadPosition();

	for(int j=1;j<=var_n;j++)
	{
		CFurnace temp_pit=pitlist.GetAt(pos);
		printf("\n");
		while(temp_pit.No!=ColInParent[j])
		{
			pitlist.GetNext(pos);
			temp_pit=pitlist.GetAt(pos);
		}
		//		temp_pit=Getpit( ColInParent[j]);
		TempCol[0]=0;
		for( i=1;i<=Col_numb;i++)
			TempCol[i]=0;
		//		fprintf(fp,"number is: %d\n",temp_pit.index);
		POSITION pos=temp_pit.ingot.GetHeadPosition();
		int loop=temp_pit.ingot.GetCount();
		for(i=1;i<=loop;i++)
		{
			//			CIngot ing=temp_pit.ingot.GetAt(temp_pit.ingot.FindIndex(i-1));
			COb ing=temp_pit.ingot.GetAt(pos);
			//			fprintf(fp,"%d	",ing.m_nIndex);
			TempCol[ing.index]=1;
			temp_pit.ingot.GetNext(pos);
		}
		objectiveCoefficient=temp_pit.cost;
		FNlist.AddTail(temp_pit.No);
		C->add(IloNum(temp_pit.cost));
		IloNumColumn col = (*cost)(temp_pit.cost);

		for (i = 0; i < Fur_numb+Col_numb; i++) 
		{
			if(i==temp_pit.index-1)
			{
				col+=(*range)[i](1);      //方案唯一性约束
			}
			else if(i<Fur_numb)            
			{
				col += (*range)[i](0);
			}
			else if(i<Fur_numb+Col_numb)                        //每个钢锭只能放在一个炉子中
			{
				col += (*range)[i](TempCol[i-Fur_numb+1]);
			}
		}

		Lamada->add(IloNumVar(col,0,1,ILOFLOAT));

		col.end();
	}

	/*加入后产生的割平面*/

	int total=var_n;
	POSITION pos_st=parentnode->CutList.GetHeadPosition();
	while(pos_st)
	{
		struct AddCut st=parentnode->CutList.GetAt(pos_st);
		IloExpr Exprr(*env);
		for(int i=0;i<total;i++)
		{
			for (int k=0;k<st.number;k++)
			{
				if(La[i+1][st.coil[k][1]]==1&&La[i+1][st.coil[k][0]]==1)
				{
					Exprr+=(*Lamada)[i];
					break;
				}
			}
		}
		IloRange rang(*env,IloNum(-IloInfinity),Exprr,st.right);
		range->add(rang);
		mod->add(rang);
		Exprr.end();
		parentnode->CutList.GetNext(pos_st);
	}


	//    fclose(fp);
	delete []TempCol;
	for(i=0;i<=columnNumb;i++)
	{
		free(La[i]);
		La[i]=NULL;
	}
	free(La);
	La=NULL;
	
	return 0;

}

int CheckBranch(CFurnace p,int r,int s,int be_median, int fur, int coi)
{
	printf("\nbe_median=%d		", be_median);
	if(be_median!=-1)
	{
		POSITION pos=p.ingot.GetHeadPosition();
		COb ing;
		bool is_in=false;
		while(pos)
		{
			ing=p.ingot.GetAt(pos);
			printf("ing.index=%d", ing.index);
			if(ing.index==be_median)
			{
				is_in=true;
				break;
			}
			printf("\nmedian[%d]=%d", ing.index, median[ing.index]);
			if(median[ing.index]==4)
				return 4;
			p.ingot.GetNext(pos);
		}
		if(is_in==false)
			return 2;
		else
			return 1;

	}
	else if(fur!=-1)
	{
		bool is_coil=false;
		bool is_furnace=false;
		printf("\np.index=%d", p.index);
		if(p.index==fur)
			is_furnace=true;
		POSITION pos=p.ingot.GetHeadPosition();
		COb ing;
		bool is_in=false;
		while(pos)
		{
			ing=p.ingot.GetAt(pos);
			if(ing.index==coi)
			{
				is_coil=true;
				break;
			}
			p.ingot.GetNext(pos);
		}
		if(is_coil==false&&is_furnace==false)
			return 2;
		else if(is_coil==true&&is_furnace==true)
			return 1;
		else if(is_coil==true&&is_furnace==false)
			return 0;
		else if(is_coil==false&&is_furnace==true)
			return 2;
	}
	else
	{
		bool r_flag=false;
		bool s_flag=false;
		int coil_is_in[C_MAX];
		for (int i=1;i<=Col_numb;i++)
		{
			coil_is_in[i]=0;
		}
		POSITION pos=p.ingot.GetHeadPosition();
		COb ing;
		while(pos)
		{
			ing=p.ingot.GetAt(pos);
			coil_is_in[ing.index]=1;
			printf("\ning.index=%d", ing.index);
			p.ingot.GetNext(pos);
		}
		bool is_find=false;
		for(int i=1;i<=Col_numb;i++)
		{
			for(int j=i+1;j<=Col_numb;j++)
			{
				if(coil_is_in[96]==1&&coil_is_in[64]==1&&i==64&&j==96&&coil_is_in[63]==1&&coil_is_in[80]==1)
					int a=0;
				if(arc_matrix[i][j]==4&&coil_is_in[i]==1&&coil_is_in[j]==1)
				{
					r_flag=false;
					is_find=true;
					break;
				}
				else if(arc_matrix[i][j]==4&&(coil_is_in[i]!=1||coil_is_in[j]!=1))
					r_flag=true;
				else if(arc_matrix[i][j]==3&&coil_is_in[i]==1&&coil_is_in[j]==1)
					r_flag=true;
				else if (arc_matrix[i][j]==3&&((coil_is_in[i]!=1&&coil_is_in[j]==1)||(coil_is_in[j]!=1&&coil_is_in[i]==1)))
				{
					r_flag=false;
					is_find=true;
					break;
				}
				else
					r_flag=true;
					
			}
			if(is_find)
				break;
		}
		if(r_flag==false)
			return 3;
		else
		{
			if(coil_is_in[r]!=1&&coil_is_in[s]!=1)
				return 2;
			else if(coil_is_in[r]==1&&coil_is_in[s]==1)
			    return 1;
			else
				return 0;
		}

	}
}

int solve_linear(IloEnv *env,IloModel *mod,IloRangeArray range,IloNumVarArray Lamada,IloNumArray C,IloObjective *cost)
{
	int i;
	IloCplex cplex(*env);//

	IloNumArray solns(*env);
	//	cplex.setParam(IloCplex::WorkMem,2000.0);
	cplex.setParam(IloCplex::AdvInd, 0);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
	cplex.setParam(IloCplex::BarDisplay,0);
	cplex.setParam(IloCplex::MIPDisplay,0);
	cplex.setParam(IloCplex::SimDisplay,0);
	cplex.setParam(IloCplex::TiLim,3600);
	cplex.extract(*mod);
	//cplex.setParam(IloCplex::TiLim,18000);
	cplex.exportModel("steel.lp");
	int a=cplex.solve();
	if(soln_of_each_node.GetCount()!=0)
		soln_of_each_node.RemoveAll();
	if(a==0)
	{
		printf("there is no feasible solution!");
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		return -1;
	}
	IloNumArray duals(*env);
	duals.clear();


	int dd=range.getSize();
	printf("%d-", dd);
	cplex.getDuals(duals,range);

	int coun=duals.getSize();
	printf("%d-",duals.getSize());
	for(i=0;i<coun;i++)
	{
		dual_var[i]=Round(duals[i],6);
		printf("dual[%d]=%f	",i,dual_var[i]);
		if(i%3==0)
		printf("\n");

	}
	i=0;
	int colnum;
	//	IloNumArray solns(env);
	solns.clear();
	cplex.getValues(Lamada,solns);


	col_num=solns.getSize();
	printf("%d-", colnum);

	//soln_of_each_node.clear();
	float value_lower=0;
	for(i=0;i<col_num;i++)
	{
		soln_of_each_node.AddTail(solns[i]);
		printf("\n%d:%f_%f		", i, C[i], solns[i]);
		value_lower+=C[i]*solns[i];
	}

	opt_soln_of_eachnode=value_lower;

	int Un_coil[C_MAX+1];//表示在某个炉子里面放了哪些板卷，被放进去的板卷的序号位置为1，如Un_coil[12]=1，表示第12号板卷放进了这个炉子里
	for(i=1;i<=Col_numb;i++)
	{
		Un_coil[i]=0;
	}
	CList<CFurnace,CFurnace> Furnacelist;

	FILE  *A;
	A=fopen("solu.txt","w");
	POSITION pos;
	CFurnace Now_pit;//当前的解的情况，即每个炉子的情况以及每个炉子内装炉情况
	pos=pitlist.GetHeadPosition();
	for(i=0;i<col_num;i++)
	{
		if(solns[i]>0)
		{
			printf("soln[%d]=%f	,value[%d]=%f	\n",i,solns[i],i,C[i]);
			Now_pit=pitlist.GetAt(pos);
			while(Now_pit.No!=FNlist.GetAt(FNlist.FindIndex(i)))
			{
				pitlist.GetNext(pos);
				Now_pit=pitlist.GetAt(pos);
			}
			if(Now_pit.ingot.GetCount()==0)
			{
				printf("Error When Getcol in GetBranchArc \n");
				getchar( );
				exit(0);
			}
			Furnacelist.AddTail(Now_pit);
			printf("%d	fu rnace:%d	%f	", i, Now_pit.index, solns[i]);
			fprintf(A,"%d	fu rnace:%d	%f	",i,Now_pit.index,solns[i]);
			printf("\n");
			for(int k=1;k<=Now_pit.lenth;k++)
			{
			    Un_coil[Now_pit.coils_in[k]]=1;
				printf("%d	", Now_pit.coils_in[k]);
				fprintf(A,"%d	",Now_pit.coils_in[k]);
			}
			fprintf(A,"	\n");
		}
	}
	fclose(A);
	printf("Total Profit=%f	feasible=%f\n",cplex.getObjValue(),result->lower_bound);


	int solu_num=0;
	double sum=0;
	for(i=0;i<col_num;i++)
	{
		if (solns[i]>=0.000001)
		{
			sum+=solns[i];
			solu_num++;
			//          printf("opt[%d]=%lf  obj[%d]=%lf\n",i,soln[i],i,obj[i]);
		}
	}
	if(nownode)//&&stage==2
	{
		if((NHB+NHS+HHB+HHS)==solu_num&&opt_soln_of_eachnode>result->lower_bound)
		{
			float im=Imp_Upbound(Furnacelist,opt_soln_of_eachnode,Un_coil);
			result->lower_bound=opt_soln_of_eachnode;
			if(im>opt_soln_of_eachnode)
			{
				result->lower_bound=im;
				up_linear_mode(env,mod,&Lamada,&C,&range,cost);
			}
		}
	}
	else
	{
		if((NHB+NHS+HHB+HHS)==solu_num&&opt_soln_of_eachnode>feasible_value)//&&stage==2
		{
			float im=Imp_Upbound(Furnacelist,opt_soln_of_eachnode,Un_coil);
			result->lower_bound=opt_soln_of_eachnode;
			feasible_value=opt_soln_of_eachnode;
			if(im>opt_soln_of_eachnode)
			{
				result->lower_bound=im;
				feasible_value=im;
				up_linear_mode(env,mod,&Lamada,&C,&range,cost);
			}
		}
	}

	if(((double)(sum-(double)(NHB+NHS+HHB+HHS)))>EPS||((double)(sum-(double)(NHB+NHS+HHB+HHS)))<-EPS)
	{
		printf("sum=%lf\n",sum);
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		return -1;
	}
	cplex.clearModel();
	cplex.clear();
	cplex.end();
	finish_time=clock();

	return 0;

}

void output(float gap, float time,int dream)
{

	char file_name[600]="result";
	char CC1[2],CC2[2];
	itoa(dream,CC2,10);
	itoa(sett,CC1,10);
	strcat(file_name,"_");
	strcat(file_name,CC1);
	strcat(file_name,"_");
	strcat(file_name,CC2);
	strcat(file_name,".txt");

	FILE *fp1;
	if((fp1=fopen(file_name,"wt+"))==NULL)
	{
		printf("Cannot open file strike any key exit 7!");
		exit(1);
	}

	for(int m=1;m<=result->solu_numb;m++)
	{
		if(result->solution[m]>0)
		{
			CFurnace Now_pit;
			Now_pit=Getpit(result->Col[result->soluindex[m]]);
			//int lenth=Now_pit.ingot.GetCount();
			fprintf(fp1,"%d\n",Now_pit.index);
			POSITION pos=Now_pit.ingot.GetHeadPosition();
			while (pos)
			{
				COb ob=Now_pit.ingot.GetAt(pos);
				fprintf(fp1,"%d	%s	%s	%d	%lf	%lf	%lf	%lf\n",ob.index,ob.m_strCoilNum,ob.m_strCurve,ob.m_nWidth,ob.m_nWeight,ob.PRI,C2[ob.index][Now_pit.center],C1[ob.index][Now_pit.index]);
				Now_pit.ingot.GetNext(pos);
			}
			fprintf(fp1,"%f  ",result->solution[m]);
			fprintf(fp1,"\n\n");
		}

	}

	fprintf(fp1,"对偶间隙为 %f\n",gap);
	fprintf(fp1,"运行时间为 %f\n",time);
	fprintf(fp1,"节点数为 %d\n",nodenum);
	fprintf(fp1,"上界界为 %lf\t",result->lower_bound);
	fprintf(fp1,"下界为 %lf\t",LBinBoot);
	fprintf(fp1,"子问题运行时间为 %lf\t",sub_time);
	fprintf(fp1,"迭代次数为 %lf\t",Sub_problem_Num);
	fprintf(fp1,"产生列数为 %d\t",total_column_number);
	fclose(fp1);
}

CFurnace Getpit(int No)
{
	POSITION pos;
	CFurnace temp;
	pos=pitlist.GetHeadPosition();
	while(pos)
	{
		temp=pitlist.GetAt(pos);
		if(temp.No==No)
			break;
		pitlist.GetNext(pos);
		// 		temp.ingot.RemoveAll();
	}
	return temp;
}

float give_cost_in_knap(int *a,int b,int n,int center)//所有a[i][k]均要求i>k
{
	// 	int flag=0;
	// 	if(arc_matrix[center][b]==3||arc_matrix[b][center]==3)
	// 		return 5000;
	// 	if(arc_matrix[center][b]==4||arc_matrix[b][center]==4)
	// 		return 1;
	// 	for(int i=0;i<=20;i++)
	// 	{
	// 		if(a[i]>=1&&a[i]<=Ingot_Num)
	// 		{
	// 			if(arc_matrix[a[i]][b]==3||arc_matrix[a[i]][b]==3)
	// 				return 5000;
	// 			if(arc_matrix[a[i]][b]==4||arc_matrix[b][a[i]]==4)
	// 				return 1;
	// //			break;
	// 		}
	// 	}
	if(arc_matrix[ingot_knap[b].index][ingot_knap[center].index]==4||arc_matrix[ingot_knap[center].index][ingot_knap[b].index]==4)
		return -1;
	else if(ingot_knap[center].collect_num>0&&ingot_knap[b].collect_num==0)
	{
		int k=0;
		while(k<ingot_knap[center].collect_num)
		{
			int ing=ingot_knap[center].ingotlist[k];
			if((arc_matrix[ing][ingot_knap[b].index]==4||arc_matrix[ingot_knap[b].index][ing]==4))
				return -1;
			k++;
		}
	}
	else if(ingot_knap[b].collect_num>0&&ingot_knap[center].collect_num==0)
	{
		int k=0;
		while(k<ingot_knap[b].collect_num)
		{
			int ing=ingot_knap[b].ingotlist[k];
			if((arc_matrix[ing][ingot_knap[center].index]==4||arc_matrix[ingot_knap[center].index][ing]==4))
				return -1;
			k++;
		}
	}
	else
	{
		int k_1=0;
		while(k_1<ingot_knap[center].collect_num)
		{
			int ing_1=ingot_knap[center].ingotlist[k_1];
			int k_2=0;
			while (k_2<ingot_knap[b].collect_num)
			{
				int ing_2=ingot_knap[b].ingotlist[k_2];
				if(arc_matrix[ing_1][ing_2]==4||arc_matrix[ing_2][ing_1]==4)
					return -1;
				k_2++;
			}
			k_1++;
		}
	}
	for(int i=0;i<=n;i++)
	{
		if(a[i]>=1&&a[i]<=Col_numb)
		{
			if(arc_matrix[ingot_knap[a[i]].index][ingot_knap[b].index]==4||arc_matrix[ingot_knap[b].index][ingot_knap[a[i]].index]==4)
				return -1;
			else if(ingot_knap[a[i]].collect_num>0&&ingot_knap[b].collect_num==0)
			{
				int k=0;
				while(k<ingot_knap[a[i]].collect_num)
				{
					int ing=ingot_knap[a[i]].ingotlist[k];
					if(arc_matrix[ing][ingot_knap[b].index]==4||arc_matrix[ingot_knap[b].index][ing]==4)
						return -1;
					k++;
				}
			}
			else if(ingot_knap[b].collect_num>0&&ingot_knap[a[i]].collect_num==0)
			{
				int k=0;
				while(k<ingot_knap[b].collect_num)
				{
					int ing=ingot_knap[b].ingotlist[k];
					if(arc_matrix[ing][ingot_knap[a[i]].index]==4||arc_matrix[ingot_knap[a[i]].index][ing]==4)
						return -1;
					k++;
				}
			}
			else
			{
				int k_1=0;
				while(k_1<ingot_knap[a[i]].collect_num)
				{
					int ing_1=ingot_knap[a[i]].ingotlist[k_1];
					int k_2=0;
					while (k_2<ingot_knap[b].collect_num)
					{
						int ing_2=ingot_knap[b].ingotlist[k_2];
						if(arc_matrix[ing_1][ing_2]==4||arc_matrix[ing_2][ing_1]==4)
							return -1;
						k_2++;
					}
					k_1++;
				}
			} 
		}
	}
	return 0;
}

void get_collect()
{
	float we=0;
	for(int i=1;i<=Col_numb;i++)
	{
		coil_inf[i].collect_num=0;
		ingot_knap[i].collect_num=0;
		for(int j=1;j<=Col_numb;j++)
		{
			if(i!=j)
			{
				//比i不好的
				if(coil_inf[i].m_strCurve==coil_inf[j].m_strCurve&&C2[i][j]<MAX&&coil_inf[i].m_nWidth<=coil_inf[j].m_nWidth&&abs(coil_inf[i].m_nOuter-coil_inf[j].m_nOuter)<=5&&coil_inf[i].m_nThick==coil_inf[j].m_nThick&&(coil_inf[j].m_nWeight+coil_inf[j].PRI<=coil_inf[i].m_nWeight+coil_inf[i].PRI||(coil_inf[j].m_nWeight+coil_inf[j].PRI==coil_inf[i].m_nWeight+coil_inf[i].PRI&&i<j)))//
				{
					coil_inf[i].collect_num++;
					coil_inf[i].ingotlist[coil_inf[i].collect_num]=j;
					valid_total_num++;
					
				}
				//比i好的
				if(coil_inf[i].m_strCurve==coil_inf[j].m_strCurve&&C2[i][j]<MAX&&coil_inf[i].m_nWidth>=coil_inf[j].m_nWidth&&abs(coil_inf[i].m_nOuter-coil_inf[j].m_nOuter)<=5&&coil_inf[i].m_nThick==coil_inf[j].m_nThick&&(coil_inf[j].m_nWeight+coil_inf[j].PRI>=coil_inf[i].m_nWeight+coil_inf[j].PRI||(coil_inf[j].m_nWeight+coil_inf[j].PRI==coil_inf[i].m_nWeight+coil_inf[j].PRI&&j>i))) //
				{
					ingot_knap[i].collect_num++;
					ingot_knap[i].ingotlist[ingot_knap[i].collect_num]=j;
					
				}
				if(coil_inf[i].m_strCurve==coil_inf[j].m_strCurve&&C2[i][j]<=500&&coil_inf[i].m_nWidth==coil_inf[j].m_nWidth)//&&(coil_inf[i].m_nWeight+coil_inf[i].PRI-coil_inf[j].m_nWeight-coil_inf[j].PRI>=0||(Round(abs(coil_inf[i].m_nWeight+coil_inf[i].PRI-coil_inf[j].m_nWeight-coil_inf[j].PRI),5)==0&&i<j))
				{
					branch_coil[i].collect_num++;
					branch_coil[i].ingotlist[branch_coil[i].collect_num]=j;
				}
			}
		}
		if(coil_inf[i].collect_num>0)
		{
			valid_number++;
			valid_maste[valid_number]=coil_inf[i];
		}
	}
}

struct Knap FindBiggest(double *price,int furnace)
{
	struct Knap pp;
	pp.ingot_num=0;
	IloEnv env2;
	IloModel mod2(env2);
	IloBoolVarArray A(env2,ingot_num_knap);
	IloExpr CC(env2);
	int *aa;
	IloExpr availExpr(env2);
	for(int i=0;i<ingot_num_knap;i++)
	{
		CC=CC+price[i+1]*A[i];
		availExpr+=price[i+1]*A[i];
	}
	mod2.add(availExpr>=dual_var[Fur_numb-1]+EPS);
	availExpr.end();
	for (int i=0;i<ingot_num_knap;i++)
	{
		for(int j=i+1;j<ingot_num_knap;j++)
		{
			int a=give_cost_in_knap(aa,i+1,-1,j+1);
			if(a==-1)
				mod2.add(A[i]+A[j]<=1);
		}
		if(price[i+1]==0)
			mod2.add(A[i]==0);
	}
	mod2.add(IloMinimize(env2, CC));

	IloCplex cplex(env2);//
	cplex.setParam(IloCplex::AdvInd, 0);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
	cplex.setParam(IloCplex::MIPDisplay,0);
	cplex.extract(mod2);
	cplex.setParam(IloCplex::EpGap,0.005);
	cplex.exportModel("steel_add.lp");
	int dd=cplex.solve();
	if(dd==0)
	{
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		mod2.removeAllProperties();
		mod2.end();
		env2.end();
		return pp;
	}
	for(int i=0;i<ingot_num_knap;i++)
	{
		if(cplex.getValue(A[i])>0.99)
		{
			pp.ingot_num++;
			pp.ingot[pp.ingot_num]=i+1;
		}
	}
	cplex.clearModel();
	cplex.clear();
	cplex.end();
	mod2.removeAllProperties();
	mod2.end();
	env2.end();
	return pp;

}

int GenerateCuts(int NumbInCol,int old_col_num,int NEW_COL_NUM,int possible,IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost,int b,int c)
{
	int flag=0;
	int       solu_numb,i;
	int       *Col;
	float   Flow[C_MAX+1][C_MAX+1];
	CList<AddCut,AddCut> TempCuts;

	/*记录模型中各个列在列池中的编号*/
	Col =new int[NEW_COL_NUM-old_col_num+NumbInCol+1];
	for(i=1;i<=NumbInCol;i++)
		Col[i]=ColInParent[i];

	for(i=NumbInCol+1;i<=NEW_COL_NUM - old_col_num +NumbInCol;i++)//??
		Col[i]=old_col_num+ i - NumbInCol;

	/*Here We Will Get The Solution of Linear Model*/     
	if(possible!=-1)
	{
/*		FILE *fp1;
		if((fp1=fopen("middle_result.txt","wt+"))==NULL)
		{
			printf("Cannot open file strike any key exit 7!");
			exit(1);
		}*/
		solu_numb=0;            //非零解个数
		int n=col_num;          //总列数
		CList<double,double> solution;   //非零解值
		CList <int,int>soluindex;     //非零解序号
		CList <CFurnace,CFurnace> temp_list;  //记录非零解所对应的组炉方案

		int **Clique;    //记录各个团约束，每个板卷所对应的行中非零变量，矩阵大小是Col_numb*col_num
		Clique=(int **)malloc((Col_numb+1)*sizeof(int*));
		for (i=0;i<=Col_numb;i++)
		{
			Clique[i]=(int *)malloc((n+1)*sizeof(int));
		}

		for (int cc=0;cc<=Col_numb;cc++)//初始化
		{
			for(int dd=0;dd<=Col_numb;dd++)
			{
				Flow[cc][dd]=0;
			}
			for(int ee=0;ee<n;ee++)
				Clique[cc][ee]=0;
		}

		int **La; //记录模型中A矩阵的系数
		La=(int **)malloc((col_num+1)*sizeof(int*));
		for (i=0;i<=col_num;i++)
		{
			La[i]=(int *)malloc((Col_numb+2)*sizeof(int));
		}

        POSITION pos=pitlist.GetHeadPosition();//
		POSITION sol=soln_of_each_node.GetHeadPosition();
		for(i=1;i<=n;i++)
		{
			CFurnace Now_pit=pitlist.GetAt(pos);
			double solu=soln_of_each_node.GetAt(sol);
			while(Now_pit.No!=Col[i])              //在列池中找到当前模型中所涉及的列
			{
//				Now_pit.ingot.RemoveAll();
				pitlist.GetNext(pos);
				Now_pit=pitlist.GetAt(pos);
			}
//			fprintf(fp1,"\nf:%d	%d	s:%f	c:%f	cen:%d	",i-1,Now_pit.index,solu,Now_pit.cost,Now_pit.center);
			if(Round(solu,6)>(double)EPS)          //记录模型中最优解值不为零的列
			{
//				fprintf(fp1,"\nf:%d	%d	s:%f	c:%f	cen:%d	",i-1,Now_pit.index,solu,Now_pit.cost,Now_pit.center);
				solu_numb++;   
				solution.AddTail(solu);
				soluindex.AddTail(i);
			}

			for(int co=1;co<=Now_pit.lenth;co++)   //记录模型A矩阵中的系数
			{
				La[i][Now_pit.coils_in[co]]=1;
//				fprintf(fp1,"%d	",Now_pit.coils_in[co]);
				if(Round(solu,6)>(double)EPS)
				{
					Clique[Now_pit.coils_in[co]][solu_numb]=1;
//					fprintf(fp1,"%d	",Now_pit.coils_in[co]);
				}
			}
			soln_of_each_node.GetNext(sol);
			Now_pit.ingot.RemoveAll();
			
		}		

		int **A;  //求割模型中系数矩阵,记录的是两个方案之间的连线
		float **Cij;//求割模型中费用系数

		A=(int **)malloc((solu_numb+1)*sizeof(int*));
		for (i=0;i<=solu_numb;i++)
		{
			A[i]=(int *)malloc((solu_numb+2)*sizeof(int));
		}

		Cij=(float **)malloc((solu_numb+1)*sizeof(float*));
		for (i=0;i<=solu_numb;i++)
		{
			Cij[i]=(float *)malloc((solu_numb+2)*sizeof(float));
		}

		for(i=0;i<=solu_numb;i++)//初始化
			for (int j=0;j<=solu_numb;j++)
			{
				Cij[i][j]=0;
				A[i][j]=0;
			}

		for(int i1=1;i1<=Col_numb;i1++)//给参数赋值
		{
			for(int k=1;k<=solu_numb;k++)
			{
				for(int k1=k+1;k1<=solu_numb;k1++)
				{
					if(Clique[i1][k]==1&&Clique[i1][k1]==1)
					{
						A[k][k1]=1;
						A[k1][k]=1;
						Cij[k][k1]=solution.GetAt(solution.FindIndex(k-1))+solution.GetAt(solution.FindIndex(k1-1));
						Cij[k1][k]=Cij[k][k1];
					}
				}
			}
		}
/*		for(i=1;i<=Col_numb;i++)
		{
			fprintf(fp1,"\nCliaue	%d:	",i);
			for(int j=1;j<=solu_numb;j++)
			{
				if(Clique[i][j]>0)
					fprintf(fp1,"%d	",j);
			}
		}
		fclose(fp1);*/

		//建立模型，求解割
		IloEnv env2;
		IloModel mod2(env2);
		IloBoolVarArray22 Make(env2);
		IloExpr CC(env2);

		for (i = 0; i <=solu_numb; i++)
		{
			Make.add(IloBoolVarArray(env2,solu_numb+1));
			
		}

		for(int i=1;i<=solu_numb;i++)
		{
			for(int j=1;j<=solu_numb;j++)
			{
				CC=CC+Cij[i][j]*Make[i][j];
			}
		}
		IloExpr availExpr2(env2);
		for (int i=1;i<=solu_numb;i++)
		{
			IloExpr availExpr(env2);
			IloExpr availExpr3(env2);
			for(int i1=1;i1<=solu_numb;i1++)
			{
				availExpr+=Make[i][i1];
				availExpr3+=Make[i1][i];
			}
			mod2.add(availExpr==availExpr3);
//			mod2.add(availExpr<=1);
			availExpr.end();
			availExpr3.end();

			for(int j=1;j<=solu_numb;j++)
			{
				mod2.add(Make[i][j]<=A[i][j]);
				availExpr2+=Make[i][j];
			}
		}
		mod2.add(availExpr2==b);
		availExpr2.end();

		mod2.add(IloMaximize(env2, CC));

		IloCplex cplex(env2);//
		cplex.setParam(IloCplex::AdvInd, 0);
		cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
		cplex.setParam(IloCplex::MIPDisplay,0);
		cplex.extract(mod2);
		cplex.exportModel("steel_add.lp");
		int dd=cplex.solve();

		int l[11][2];// 方案对
		int counter=0;
		
		if(dd!=0)
		{
//			struct AddCut st;
//			int count_number=0;
			printf("\nvalue:%f\n",cplex.getObjValue());
			for(int ii=1;ii<= solu_numb;ii++)
			{
				for (int j=1;j<= solu_numb;j++)
				{
					float a=cplex.getValue(Make[ii][j]);
					if(a>0.99)
					{
						printf("\n%d->%d	%f\n",ii,j,a);

						l[++counter][0]=ii;
						l[counter][1]=j;
					}
				}
			}
			float cos=cplex.getObjValue();
			int Coil_point[C_MAX];
			counter=0;
			if(Round(cos-c,6)>EPS)
			{
				for(int j=1;j<=b;j++)  //找出团约束所对应的板卷
				{
					for(int i=1;i<=Col_numb;i++)
					{
						if(Clique[i][l[j][0]]==1&&Clique[i][l[j][1]]==1)
						{
							Coil_point[++counter]=i;
//							printf("%d:%d\n",counter,i);
							break;
						}
					}
				}

				//for(i=1;i<=counter;i++)
				//	printf("%d	\n",Coil_point[i]);
				struct AddCut st;
				for(int j=1;j<=b;j++)
				{
					int coil1[2];
					int cont=0;
					for(int i=1;i<=b;i++)
					{
						if(Clique[Coil_point[i]][l[j][0]]==1)
						{
							coil1[cont]=Coil_point[i];
							cont++;
						}
						if(cont>=2)
							break;
					}
					st.coil[j-1][0]=coil1[0];
					st.coil[j-1][1]=coil1[1];
				}
				st.number=b;
				st.right=c/2;

				TempCuts.AddTail(st);
				CurrentCutList.AddTail(st);
			}
			
			if(dd!=0&&cos>c)
			{
				int total=soln_of_each_node.GetCount();
				POSITION pos=TempCuts.GetHeadPosition();
				while(pos)
				{
					struct AddCut st=TempCuts.GetAt(pos);
					IloExpr Exprr(*env);
					for(int i=0;i<total;i++)
					{
						for (int k=0;k<st.number;k++)
						{
							if(La[i+1][st.coil[k][1]]==1&&La[i+1][st.coil[k][0]]==1)
							{
								Exprr+=(*Lamada)[i];
								printf("%d	",i);
								break;
							}
						}
					}
					printf("\n");
					IloRange rang(*env,IloNum(-IloInfinity),Exprr,c/2);
					//IloRange rang(*env,IloNum(-IloInfinity),1);
					//mod->remove(*range);
					range->add(rang);
					mod->add(rang);
					//mod->add(*range);
					Exprr.end();
					flag=1;
					TempCuts.GetNext(pos);
				}
			}
		}
		
		for(i=0;i<=Col_numb;i++)
		{
			free(Clique[i]);
			Clique[i]=NULL;
		}
		free(Clique);
		Clique=NULL;

		for(i=0;i<=col_num;i++)
		{
			free(La[i]);
			La[i]=NULL;
		}
		free(La);
		La=NULL;

		for(i=0;i<=solu_numb;i++)
		{
			free(Cij[i]);
			Cij[i]=NULL;
		}
		free(Cij);
		La=NULL;

		for(i=0;i<=solu_numb;i++)
		{
			free(A[i]);
			A[i]=NULL;
		}
		free(A);
		La=NULL;


		cplex.clearModel();
		cplex.clear();
		cplex.end();
		Make.end();
		CC.end();
		mod2.removeAllProperties();
		mod2.end();
		env2.end();
	}

	delete []Col;
	return flag;
}
int GenerateCutss(int NumbInCol,int old_col_num,int NEW_COL_NUM,int possible,IloEnv *env, IloModel *mod, IloNumVarArray *Lamada, IloNumArray *C,IloRangeArray *range,IloObjective *cost,int b,int c)
{
	int flag=0;
	int       solu_numb,i;
	int       *Col;
	float   Flow[C_MAX+1][C_MAX+1];
	CList<AddCut,AddCut> TempCuts;

	/*记录模型中各个列在列池中的编号*/
	Col =new int[NEW_COL_NUM-old_col_num+NumbInCol+1];
	for(i=1;i<=NumbInCol;i++)
		Col[i]=ColInParent[i];

	for(i=NumbInCol+1;i<=NEW_COL_NUM - old_col_num +NumbInCol;i++)//??
		Col[i]=old_col_num+ i - NumbInCol;

	/*Here We Will Get The Solution of Linear Model*/     
	if(possible!=-1)
	{
/*		FILE *fp1;
		if((fp1=fopen("middle_result.txt","wt+"))==NULL)
		{
			printf("Cannot open file strike any key exit 7!");
			exit(1);
		}*/
		solu_numb=0;            //非零解个数
		int n=col_num;          //总列数
		CList<double,double> solution;   //非零解值
		CList <int,int>soluindex;     //非零解序号
		CList <CFurnace,CFurnace> temp_list;  //记录非零解所对应的组炉方案

		int **Clique;    //记录各个团约束，每个板卷所对应的行中非零变量，矩阵大小是Col_numb*col_num
		Clique=(int **)malloc((Col_numb+1)*sizeof(int*));
		for (i=0;i<=Col_numb;i++)
		{
			Clique[i]=(int *)malloc((n+1)*sizeof(int));
		}

		for (int cc=0;cc<=Col_numb;cc++)//初始化
		{
			for(int dd=0;dd<=Col_numb;dd++)
			{
				Flow[cc][dd]=0;
			}
			for(int ee=0;ee<n;ee++)
				Clique[cc][ee]=0;
		}

		int **La; //记录模型中A矩阵的系数
		La=(int **)malloc((col_num+1)*sizeof(int*));
		for (i=0;i<=col_num;i++)
		{
			La[i]=(int *)malloc((Col_numb+2)*sizeof(int));
		}

		POSITION pos=pitlist.GetHeadPosition();//
		POSITION sol=soln_of_each_node.GetHeadPosition();
		for(i=1;i<=n;i++)
		{
			CFurnace Now_pit=pitlist.GetAt(pos);
			double solu=soln_of_each_node.GetAt(sol);
			while(Now_pit.No!=Col[i])              //在列池中找到当前模型中所涉及的列
			{
				//				Now_pit.ingot.RemoveAll();
				pitlist.GetNext(pos);
				Now_pit=pitlist.GetAt(pos);
			}
			//			fprintf(fp1,"\nf:%d	%d	s:%f	c:%f	cen:%d	",i-1,Now_pit.index,solu,Now_pit.cost,Now_pit.center);
			if(Round(solu,6)>(double)EPS)          //记录模型中最优解值不为零的列
			{
//				fprintf(fp1,"\nf:%d	%d	%d	s:%f	c:%f	cen:%d	",solu_numb+1,i-1,Now_pit.index,solu,Now_pit.cost,Now_pit.center);
				solu_numb++;   
				solution.AddTail(solu);
				soluindex.AddTail(i);
				temp_list.AddTail(Now_pit);
			}

			for(int co=1;co<=Now_pit.lenth;co++)   //记录模型A矩阵中的系数
			{
				La[i][Now_pit.coils_in[co]]=1;
				//				fprintf(fp1,"%d	",Now_pit.coils_in[co]);
				if(Round(solu,6)>(double)EPS)
				{
					Clique[Now_pit.coils_in[co]][solu_numb]=1;
					
//					fprintf(fp1,"%d	",Now_pit.coils_in[co]);
				}
			}
			soln_of_each_node.GetNext(sol);
			Now_pit.ingot.RemoveAll();

		}


		int **A;  //求割模型中系数矩阵,记录的是两个方案之间的连线
		float **Cij;//求割模型中费用系数

		A=(int **)malloc((solu_numb+1)*sizeof(int*));
		for (i=0;i<=solu_numb;i++)
		{
			A[i]=(int *)malloc((solu_numb+2)*sizeof(int));
		}

		Cij=(float **)malloc((solu_numb+1)*sizeof(float*));
		for (i=0;i<=solu_numb;i++)
		{
			Cij[i]=(float *)malloc((solu_numb+2)*sizeof(float));
		}

		for(i=0;i<=solu_numb;i++)//初始化
			for (int j=0;j<=solu_numb;j++)
			{
				Cij[i][j]=0;
				A[i][j]=0;
			}

			for(int i1=1;i1<=Col_numb;i1++)//给参数赋值
			{
				for(int k=1;k<=solu_numb;k++)
				{
					for(int k1=k+1;k1<=solu_numb;k1++)
					{
						if(Clique[i1][k]==1&&Clique[i1][k1]==1)
						{
							A[k][k1]=1;
							A[k1][k]=1;
							Cij[k][k1]=solution.GetAt(solution.FindIndex(k-1))+solution.GetAt(solution.FindIndex(k1-1));
							Cij[k1][k]=Cij[k][k1];
						}
					}
				}
			}
/*			for(i=1;i<=Col_numb;i++)
			{
				fprintf(fp1,"\nCliaue	%d:	",i);
				for(int j=1;j<=solu_numb;j++)
				{
					if(Clique[i][j]>0)
						fprintf(fp1,"%d	",j);
				}
			}
			fclose(fp1);*/



			int sol1,sol2,sol3;
			for(int sol1=1;sol1<=solu_numb;sol1++)
			{
				CFurnace temppit=temp_list.GetAt(temp_list.FindIndex(sol1-1));
				float vv=solution.GetAt(solution.FindIndex(sol1-1));  //解的值
				if(vv>=1)
					continue;
				float fra_value;
				int total=temppit.lenth;
				for (int i=1;i<=total;i++)
				{
					
					for (int j=i+1;j<=total;j++)
					{
						for(int k=1;k<=Col_numb;k++)
						{
							fra_value=vv;
							int a=temppit.coils_in[i];
							int b=temppit.coils_in[j];
							if (k!=a&&k!=b)
							{
								bool is_find=false;
								int a_number=0;
								int b_number=0;
								for(sol2=1;sol2<=solu_numb;sol2++)
								{
									if(sol2!=sol1)
									{
										if(Clique[a][sol2]==1&&Clique[k][sol2]==1&&(a_number==0||b_number!=0))
										{
											fra_value+=solution.GetAt(solution.FindIndex(sol2-1));
											a_number++;
										}
										else if(Clique[b][sol2]==1&&Clique[k][sol2]==1)
										{
											fra_value+=solution.GetAt(solution.FindIndex(sol2-1));
											b_number++;
										}
										else if(Clique[b][sol2]==1&&Clique[a][sol2]==1)
										{
											fra_value+=solution.GetAt(solution.FindIndex(sol2-1));
										}
									}
								}
								if(Round(fra_value-1.2,6)>EPS&&b_number!=0&&a_number!=0)
								{
									if(a>b)
									{ 
										int d=a;a=b;b=d;
									}
									if(a>k)
									{ 
										int d=a;a=k;k=d;
									}
									if(b>k)
									{ 
										int d=b;b=k;k=d;
									}

									bool repeat_flag=false;
									POSITION pos=TempCuts.GetHeadPosition();
									while (pos)
									{
										struct AddCut ff=TempCuts.GetAt(pos);
										if(ff.coil[0][0]==a&&ff.coil[1][0]==b&&ff.coil[2][0]==k)
										{
											repeat_flag=true;
											break;
										}
										TempCuts.GetNext(pos);
									}
									if(repeat_flag==true)
										continue;
									struct AddCut st;
									st.coil[0][0]=a;
									st.coil[0][1]=b;

									st.coil[1][0]=b;
									st.coil[1][1]=k;

									st.coil[2][0]=k;
									st.coil[2][1]=a;

									st.number=3;
									st.right=c/2;

									printf("%d	%d	%d	%f",a,b,k,fra_value);

									TempCuts.AddTail(st);
									CurrentCutList.AddTail(st);

									IloExpr Exprr(*env);
									int to=soln_of_each_node.GetCount();
									for(int i=0;i<to;i++)
									{
										for (int k=0;k<st.number;k++)
										{
											if(La[i+1][st.coil[k][1]]==1&&La[i+1][st.coil[k][0]]==1)
											{
												Exprr+=(*Lamada)[i];
												//printf("%d	",i);
												break;
											}
										}
									}
									//printf("\n");
									IloRange rang(*env,IloNum(-IloInfinity),Exprr,c/2);
									range->add(rang);
									mod->add(rang);
									Exprr.end();
									flag=1;
								}
							}
						}
					}
				}
			}


			//建立模型，求解割
			for(i=0;i<=Col_numb;i++)
			{
				free(Clique[i]);
				Clique[i]=NULL;
			}
			free(Clique);
			Clique=NULL;

			for(i=0;i<=col_num;i++)
			{
				free(La[i]);
				La[i]=NULL;
			}
			free(La);
			La=NULL;

			for(i=0;i<=solu_numb;i++)
			{
				free(Cij[i]);
				Cij[i]=NULL;
			}
			free(Cij);
			La=NULL;

			for(i=0;i<=solu_numb;i++)
			{
				free(A[i]);
				A[i]=NULL;
			}
			free(A);
			La=NULL;



	}

	delete []Col;
	return flag;
}
int Get_Number_Upper()
{
	IloEnv env2;
	IloModel mod2(env2);
	IloNumVarArray22 Xij(env2);              //决策变量板卷放到炉子中
	IloExpr CC(env2);                         //目标函数，最大化，总的装包费用

	for (int i = 0; i <=Col_numb; i++)
	{
		Xij.add(IloNumVarArray(env2,Fur_numb+1,0,1)); //初始化决策变量
	}

	//加入约束
	/*能力约束*/
	
	for (int j=1;j<=Fur_numb;j++)
	{
		IloExpr availExpr(env2);
		for (int i=1;i<=Col_numb;i++)
		{
			availExpr+=(coil_inf[i].m_nWidth+70)*Xij[i][j];
			for(int k=1;k<=Col_numb;k++)
			{
				if(C2[i][k]>MAX)
					mod2.add(Xij[i][j]+Xij[k][j]<=1);
			}
			if(C1[i][j]>MAX)
				Xij[i][j]==0;
		}
		mod2.add(availExpr<=4700-70);
		
		availExpr.end();
	}
	for(int i=1;i<=Col_numb;i++)
	{
		IloExpr availExpr(env2);
		for (int j=1;j<=Fur_numb;j++)
		{
			availExpr+=Xij[i][j];
		}
		mod2.add(availExpr<=1);
		availExpr.end();
	}



	/*加入目标函数*/
	for(int j=1;j<=Fur_numb;j++)
	{
		for(int i=1;i<=Col_numb;i++)
			CC+=Xij[i][j];
	}

	mod2.add(IloMaximize(env2, CC));//最大化问题

	IloCplex cplex(env2);//采用CPLEX求解
	cplex.setParam(IloCplex::AdvInd, 0);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
//	cplex.setParam(IloCplex::MIPDisplay,0);
	//			cplex.setParam(IloCplex::SolnPoolAGap,0.5);
	cplex.extract(mod2);
	cplex.exportModel("number_upper.lp");
	int dd=cplex.solve();
	float number=cplex.getObjValue();
	int a=number;
	cplex.clearModel();
	cplex.clear();
	cplex.end();
	Xij.end();
	CC.end();
	mod2.removeAllProperties();
	mod2.end();
	env2.end();
	number_upperbound=a;
	return a;
}

void Variable_Fixing()
{


	double *price;
	price=new double[Col_numb+1];

	FILE *fpp;
	if((fpp=fopen("Variable_Fixing.txt","wt+"))==NULL)
	{
		printf("Cannot open file strike any key exit 7!");
		exit(1);
	}
	float aa[C_MAX];
	for(int i=1;i<=Col_numb;i++)
		aa[i]=0;

	for(int k=1;k<=Col_numb;k++)
	{
		for(int j=1;j<=Fur_numb;j++)
		{
			float dream=coil_inf[k].m_nWeight+coil_inf[k].PRI-C1[k][j];
			if(dream>aa[k])
				aa[k]=dream;
			if(C1[k][j]>C_MAX*1000||fur_coil[j][k]==4||median[k]==4)//
				continue;
			float temp_value=coil_inf[k].m_nWeight+coil_inf[k].PRI-C1[k][j]-dual_var[k+Fur_numb-1];
			for(int l=1;l<=Fur_numb;l++)
			{
				if(temp_value>eij[k][l])
					eij[k][l]=temp_value;
			}
			for (int l=1;l<=Col_numb;l++)
			{
				if(temp_value>eij[l][j])
					eij[l][j]=temp_value;
			}
			

		}
	}

	float dual_value=0;
	for(int i=1;i<=Col_numb;i++)
	{
		dual_value+=aa[i];
	}

	dual_value=dual_value+((float)NHB)*dual_var[0];
	dual_value=dual_value+((float)HHB)*dual_var[1];
	dual_value=dual_value+((float)NHS)*dual_var[2];
	dual_value = dual_value + ((float)HHS)*dual_var[3];

	fprintf(fpp,"\n");
	fprintf(fpp,"median\n");

	for(int i=1;i<=Col_numb;i++)//Xij
	{
		for(int j=1;j<=Fur_numb;j++)
		{
			if(C1[i][j]>C_MAX*1000||fur_coil[j][i]==4||median[i]==4||dual_var[j-1]<0||Ykk[i][j]==4)//
				continue;
			float a=Round(dkk[i][j]-dual_var[j-1],3);
			float b=result->lower_bound-Round(dual_value,3);
			if(Round(dkk[i][j]-dual_var[j-1],3)-result->lower_bound+Round(opt_soln_of_eachnode,3)<-0.01)//opt_soln_of_eachnode
			{
				Ykk[i][j]=4;
				fprintf(fpp,"k: %d	j:%d\n",i,j);
			}
		}
	}
	fclose(fpp);
	delete []price;
}
float Imp_Upbound(CList<CFurnace,CFurnace> &Flist,float cur_value,int *UCoil)
{
    if(sett<=10)
	    goto L2;
	float imp_value=cur_value;
	POSITION f_pos;
	f_pos=Flist.GetHeadPosition();
	while (f_pos)
	{
		CFurnace fur=Flist.GetAt(f_pos);
		BOOL flag1=false;
		for(int kk=1;kk<=fur.lenth;kk++)
		{
			for (int ll=1;ll<=fur.lenth;ll++)
			{
				if(arc_matrix[fur.coils_in[ll]][fur.coils_in[kk]]==3)
				{
					flag1=true;
					break;
				}
			}
			if(flag1==true)
			break;
		}
		if (flag1==true)
		{
			Flist.GetNext(f_pos);
			continue;
		}
		float total_height=0;
		for (int l=1;l<=fur.lenth;l++)
		{
			total_height+=coil_inf[fur.coils_in[l]].m_nWidth+70;
		}
		total_height=total_height-70;
		BOOL is_imp=false;
		POSITION co_pos=fur.ingot.GetHeadPosition();
		while(co_pos)
		{
			COb ob=fur.ingot.GetAt(co_pos);
			bool is_find=false;
			int i;
			for (i=1;i<=fur.lenth;i++)
			{
				if (fur.coils_in[i]==ob.index)
				{
					is_find=true;
					break;
				}
			}
			if(is_find==false)
			{

				fur.ingot.GetNext(co_pos);
				continue;
			}
			float most=0;
			int most_index=-1;
			POSITION max_pos;
			
			for (int j=1;j<=Col_numb;j++)
			{
				CFurnace temp=fur;
				if(fur.coils_in[i]==22&&j==7)
					int a=0;
				BOOL flag2=false;
				for(int kk=1;kk<=fur.lenth;kk++)
				{
					if(arc_matrix[fur.coils_in[kk]][j]==4)
						flag2=true;
				}
 				if(UCoil[j]==0&&(C2[fur.center][j]<MAX)&&(C1[j][fur.index]<MAX)&&(total_height-coil_inf[fur.coils_in[i]].m_nWidth+coil_inf[j].m_nWidth<=4700)&&(fur.coils_in[i]!=j)&&(fur.coils_in[i]!=fur.center)&&(flag2==false))
				{
					temp.coils_in[i]=j;
					POSITION i_pos=temp.ingot.GetHeadPosition();
					while (i_pos)
					{
						COb ob=temp.ingot.GetAt(i_pos);
						if(ob.index==fur.coils_in[i])
						{
							is_find=true;
							break;
						}
						temp.ingot.GetNext(i_pos);
					}
					if(is_find==true)
					{
						temp.ingot.SetAt(i_pos,coil_inf[j]);
					}

					float before=imp_value;
					float after=imp_value-(coil_inf[fur.coils_in[i]].m_nWeight+coil_inf[fur.coils_in[i]].PRI-C2[fur.center][fur.coils_in[i]]-C1[fur.coils_in[i]][fur.index])+(coil_inf[j].m_nWeight+coil_inf[j].PRI-C2[fur.center][j]-C1[j][fur.index]);
					if(after>before&&after>most)
					{
						most=after;
						most_index=j;
						
					}
				}
			}
			if(most_index!=-1)
			{
			    int pre_index=fur.coils_in[i];
				fur.coils_in[i]=most_index;
				fur.ingot.SetAt(co_pos,coil_inf[most_index]);
				fur.cost=GetValue(fur); //获得列所对应的目标函数
				Flist.SetAt(f_pos,fur);
				imp_value=most;
				UCoil[most_index]=1;
				UCoil[pre_index]=0;
				is_imp=true;
//				break;
			}
			fur.ingot.GetNext(co_pos);
		}
		if(is_imp)
			pit_pool.AddTail(fur);
		Flist.GetNext(f_pos);
	}
	if(imp_value>cur_value)
	{
		printf("%f\n",imp_value);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n");
	}
	return imp_value;
L2:	return 0;
}