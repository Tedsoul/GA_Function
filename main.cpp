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
//�Ŵ��㷨 GA for jiga

using namespace std;
using   std::numeric_limits;

#define Div 30							//ά��            *****��Ҫ����******   
#define size_population 1000		//��Ⱥ��С
#define EVO_maxGen 100				//��������
#define numOfRun_Fun 30			//����һ�������Ĵ���

#define rnd(low, uper) ((rand() / (double)RAND_MAX)*((uper)-(low)) + (low))

#define possibility_Crossover 0.8		//������� 
#define possibility_Mutation 0.05		//������� ���0.03����ʱ������0.8��
#define sizeOfTournament 0.2			//��������ģ
double   lowBound[Div], upperBound[Div];  //��Ⱥ�и���ķ�Χ��   *****��Ҫ����******  



class Individual{
public:
	double X[Div];
	double fitness_now;

	double(*myFun)(double present[], int DIM) = f1;
public:
	//��ʼ�����壬������ɶ����������꣬��������ʷ�������Ϊ������꣬ͬʱ������Ӧֵ
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
	//���㽻��
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
		/*���̶�ѡ��GA����������Ӧֵ���ģ�
		�������Ӧֵ��С�Ļ���ע�⡾���ֻ�Ĵ���С�ںš������̶ĵĸ��ʾͻᷴ�򣬱�ɺõĽ����С��
		���ֱ�ӽ���Ӧֵ��ɵ�������Ҫע�⣬��ӦֵΪ��������������ǵ�����������0�Ĺ�ϵ��
		*/
		double fitness_Sum = 0;
		for (int i = 0; i < size_population; i++)
			fitness_Sum += population[i].fitness_now;
		
		//possibility
		double * possibility = new double[size_population];
		for (int i = 0; i < size_population; i++){
			possibility[i] = population[i].fitness_now / fitness_Sum;
		}
		//��������Ⱥ
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

	void  selection_tournament(){//������ѡ�񷽷� tournament selection model
		
		//��������Ⱥ
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
	void  selection_theBest(){//��Ѹ��屣�淨1 �ҳ��������2 ����ѱȼ�¼��Ѻã�����¼���˹��¼�� 3 �����Ļ������
	update_MyBest();
	//���
	double tmp = (population[0].fitness());
	int worstIndex = 0;

	for (int i = 0; i < size_population; i++){
	if ((population[i].fitness()) > tmp){//***************************
	tmp = population[i].fitness();
	worstIndex = i;
	}
	}
	//����
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
		//ϴ��
		int index, tmp;
		for (int i = (cross->size() - 1); i>0; i--)
		{
			index = rand() % i;
			tmp = (*cross)[i];
			(*cross)[i] = (*cross)[index];
			(*cross)[index] = tmp;
		}
		//����
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

	void update_MyBest(){//������Ⱥ����Ӧֵ��С��Ⱦɫ��
		//������Ӧֵ
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
			//����
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
		//������Ⱥ����Ӧֵ����Ⱦɫ��
		double minFitness = (population[0].fitness());
		int minIndex = 0;

		for (int i = 0; i < size_population; i++){
			if ((population[i].fitness()) < minFitness){//***************************
				minFitness = population[i].fitness();
				minIndex = i;
			}
		}
		//����
		for (int i = 0; i < Div; i++)
			MyBest_pos[i] = population[minIndex].X[i];
		MyBest = minFitness;
	}


	void run(){
		int gen = 0;
		//run֮ǰ����init;init������ �������Ⱦɫ��
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