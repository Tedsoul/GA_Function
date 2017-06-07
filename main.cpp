#include <cstdio>
#include <cstdlib> 
#include <cmath>     
#include <ctime> 
#include <limits>  
#include <iostream> 
#include<fstream>
#include<string>
#include<queue>
#include<vector>
#include"myfun.h"
void test_For(string s1, double(*Fun)(double present[], int DIM), double low, double high);
//遗传算法 GA

using namespace std;
using   std::numeric_limits;

#define Div 30							//维度            *****需要更改******   
#define size_population 1000		//种群大小
#define EVO_maxGen 100				//迭代次数
#define numOfRun_Fun 30			//运行一个函数的次数

#define rnd(low, uper) ((rand() / (double)RAND_MAX)*((uper)-(low)) + (low))

#define possibility_Crossover 0.8		//交配概率 
#define possibility_Mutation 0.05		//变异概率 最佳0.03【此时交配是0.8】
#define sizeOfTournament 0.2			//竞标赛规模
double   lowBound[Div], upperBound[Div];  //种群中个体的范围，   *****需要更改******  



class Individual{
public:
	double X[Div];
	double fitness_now;

	double(*myFun)(double present[], int DIM) = f1;
public:
	//初始化个体，随机生成定义域内坐标，并设置历史最好坐标为这个坐标，同时计算适应值
	void init(){
		for (int i = 0; i < Div; i++){
			X[i] = rnd(lowBound[i], upperBound[i]);
		}
		fitness_now = fitness();
	}

	void setFun(double(*tmp)(double[], int)){
		this->myFun = tmp;
	}

	double fitness(){
		fitness_now = myFun(X, Div);
		return fitness_now;
	}

	void copy_From(Individual &tmp){
		fitness_now = tmp.fitness_now;

		for (int j = 0; j < Div; j++){
			X[j] = tmp.X[j];
		}
	}
	//单点交叉
	void crossover_With(Individual &tmp){
		int cross_pos = rnd(1, Div - 2);
		for (int i = 0; i < cross_pos; i++){
			double tmp_double = X[i];
			X[i] = tmp.X[i];
			tmp.X[i] = tmp_double;
		}
	}
};

class GeneticAlgrothim{
public:
	Individual population[size_population];
	Individual newPopulation[size_population];
	double MyBest;
	double MyBest_pos[Div];


	void  selection_Bet(){
		/*轮盘赌选择，GA本来是求适应值最大的，
		如果求适应值最小的话，注意【如果只改大于小于号】，轮盘赌的概率就会反向，变成好的解概率小，
		如果直接将适应值变成倒数，又要注意，适应值为负数的情况，不是单纯的无穷大跟0的关系。
		*/
		double fitness_Sum = 0;
		for (int i = 0; i < size_population; i++)
			fitness_Sum += population[i].fitness_now;
		
		//possibility
		double * possibility = new double[size_population];
		for (int i = 0; i < size_population; i++){
			possibility[i] = population[i].fitness_now / fitness_Sum;
		}
		//生成新种群
		for (int i = 0; i < size_population; i++)
		{
			//select
			double num_rnd = rnd(0, 1);
			double tmp_possible = 0;
			int index;
			for (int j = 0; j < size_population; j++){
				index = j;
				tmp_possible += possibility[j];
				if (num_rnd <= tmp_possible)break;
			}

			newPopulation[i].copy_From(population[index]);
		}
		delete[] possibility;
	}

	void  selection_tournament(){//锦标赛选择方法 tournament selection model
		
		//生成新种群
		for (int i = 0; i < size_population; i++){
			//select
			int selectIndex = (int)rnd(0, size_population - 1);
			int local_Best = selectIndex;
			for (int j = 1; j < (sizeOfTournament * size_population); j++){
				selectIndex = (int)rnd(0, size_population - 1);
				if ((population[selectIndex].fitness_now) < (population[local_Best].fitness_now))//*************
					local_Best = selectIndex;
			}
			newPopulation[i].copy_From(population[local_Best]);
		}
	}


	/*
	void  selection_theBest(){//最佳个体保存法1 找出最佳与最差。2 若最佳比记录最佳好，则更新吉尼斯记录。 3 把最差的换成最佳
	update_MyBest();
	//最差
	double tmp = (population[0].fitness());
	int worstIndex = 0;

	for (int i = 0; i < size_population; i++){
	if ((population[i].fitness()) > tmp){//***************************
	tmp = population[i].fitness();
	worstIndex = i;
	}
	}
	//保存
	for (int i = 0; i < Div; i++)
	population[worstIndex].X[i] = MyBest_pos[i];

	population[worstIndex].fitness_now = MyBest;

	}
	*/
	void  crossover(){
		queue<int> cross;
		for (int i = 0; i < size_population; i++){
			if (rnd(0, 1) < possibility_Crossover)
				cross.push(i);
		}

		for (int i = 0; i < (int)(cross.size() / 2); i++)
		{
			int q1 = cross.front(); cross.pop();
			int q2 = cross.front(); cross.pop();
			newPopulation[q1].crossover_With(newPopulation[q2]);
		}

	}

	void  crossover_Random(){
		vector<int> *cross = new vector<int>();
		for (int i = 0; i < size_population; i++){
			if (rnd(0, 1) < possibility_Crossover)
				cross->push_back(i);
		}
		//洗牌
		int index, tmp;
		for (int i = (cross->size() - 1); i>0; i--)
		{
			index = rand() % i;
			tmp = (*cross)[i];
			(*cross)[i] = (*cross)[index];
			(*cross)[index] = tmp;
		}
		//交叉
		for (int i = 0; i < (int)(cross->size() / 2); i++)
		{
			int q1 = (*cross)[i];
			int q2 = (*cross)[(cross->size() - 1 - i)];
			newPopulation[q1].crossover_With(newPopulation[q2]);
		}
		cross->clear();
		delete cross;
	}

	void  mutation(){
		for (int i = 0; i < size_population; i++){
			for (int j = 0; j < Div; j++){
				if (rnd(0, 1) < possibility_Mutation)population[i].X[j] = rnd(lowBound[j], upperBound[j]);
			}
		}
	}
public:

	void copyNewPopulation(){
		for (int i = 0; i < size_population; i++)
			population[i].copy_From(newPopulation[i]);
	}
	void setBound(const double & low, const double & up){
		for (int i = 0; i < Div; i++){
			lowBound[i] = low;
			upperBound[i] = up;
		}
	}

	void update_MyBest(){//保存种群中适应值最小的染色体
		//评价适应值
		for (int i = 0; i < size_population; i++)
			population[i].fitness();

		double minFitness = (population[0].fitness());
		int minIndex = 0;

		for (int i = 0; i < size_population; i++){
			if ((population[i].fitness()) < minFitness){//***************************
				minFitness = population[i].fitness();
				minIndex = i;
			}
		}
		if (minFitness < MyBest){//***************************
			//保存
			for (int i = 0; i < Div; i++)
				MyBest_pos[i] = population[minIndex].X[i];
			MyBest = minFitness;
		}
	}
	void init(double lo, double hi){
		srand((unsigned)time(NULL));
		setBound(lo, hi);

		//init
		for (int i = 0; i < size_population; i++){
			population[i].init();
		}
		//保存种群中适应值最大的染色体
		double minFitness = (population[0].fitness());
		int minIndex = 0;

		for (int i = 0; i < size_population; i++){
			if ((population[i].fitness()) < minFitness){//***************************
				minFitness = population[i].fitness();
				minIndex = i;
			}
		}
		//保存
		for (int i = 0; i < Div; i++)
			MyBest_pos[i] = population[minIndex].X[i];
		MyBest = minFitness;
	}


	void run(){
		int gen = 0;
		//run之前运行init;init包含了 保存最佳染色体
		while (gen < EVO_maxGen){
			selection_tournament();

			crossover_Random();

			mutation();

			copyNewPopulation();

			update_MyBest();

			gen++;
		}
		cout << "gen:" << gen << endl;
		//for (int i = 0; i < size_population;i++)
		//cout << "The lowest value is" << population[i].pBest << endl;
		cout << "The lowest value is" << MyBest << endl;

		cout << endl;
	}

};
void test(){
	GeneticAlgrothim p;
	fstream file("D:\\data.txt", ios::out);
	file.close();
	double(*Fun[])(double present[], int DIM) = { f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 };

	test_For("f1", f1, -100, 100);
	test_For("f2", f2, -10, 10);
	test_For("f3", f3, -100, 100);
	test_For("f4", f4, -100, 100);
	test_For("f5", f5, -30, 30);
	test_For("f6", f6, -100, 100);
	test_For("f7", f7, -1.28, 1.28);
	test_For("f8", f8, -500, 500);
	test_For("f9", f9, -5.12, 5.12);
	test_For("f10", f10, -32, 32);
	test_For("f11", f11, -600, 600);
	test_For("f12", f12, -50, 50);
	test_For("f13", f13, -50, 50);
}
void test_For(string s1, double(*Fun)(double present[], int DIM), double low, double high){
	GeneticAlgrothim p;
	for (int i = 0; i < Div; i++){
		p.population[i].myFun = Fun;
	}

	double store[numOfRun_Fun] = { 0 };

	fstream file1("D:\\data.txt", ios::app);
	file1 << "Now,Function  " << s1 << endl;
	file1.close();

	for (int i = 0; i < numOfRun_Fun; i++){
		p.init(low, high);
		p.run();
		fstream file2("D:\\data.txt", ios::app);
		file2 << "No. " << (i + 1) << "  The lowest value is  " << p.MyBest << endl;
		file2.close();

		store[i] = p.MyBest;
	}
	fstream file3("D:\\data.txt", ios::app);
	file3 << "ave : " << ave(store, numOfRun_Fun) << endl << "Standard Deviation:  " << standardDeviation(store, numOfRun_Fun) << endl;
	file3.close();

}
void testOnlyOne(){
	GeneticAlgrothim p;
	fstream file("D:\\data.txt", ios::out);
	file.close();
	for (int i = 0; i < Div; i++){
		p.population[i].myFun = f13;
	}

	double store[30] = { 0 };

	for (int i = 0; i < 30; i++){
		p.init(-50, 50);
		p.run();
		fstream file2("D:\\data.txt", ios::app);
		//file2 << "No. " << (i + 1) << "  The lowest value is  " << p.gBest << endl;
		file2.close();

		//store[i] = p.gBest;
	}
	fstream file3("D:\\data.txt", ios::app);
	file3 << "ave : " << ave(store, 30) << endl << "Standard Deviation:  " << standardDeviation(store, 30) << endl;
	file3.close();

}
void printOnly(){
	GeneticAlgrothim p;
	for (int i = 0; i < Div; i++){
		p.population[i].myFun = f1;
	}
	for (int i = 0; i < 30; i++){
		p.init(-100, 100);
		p.run();
	}
}
int main(){
	//testOnlyOne();
	test();
}
