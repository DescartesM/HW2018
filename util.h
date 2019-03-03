/*
 * util.h
 *
 *  Created on: 2018年4月8日
 *      Author: customer
 */

#ifndef UTIL_H_
#define UTIL_H_
#include "lib_io.h"
#include "pinv.h"

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <string.h>
#include <string>
#include <vector>
#include <ctime>
#include <unordered_map>
using namespace std;
struct VM{
	string vm_name;
	int vm_cpu;
	int vm_memory;
	int vm_buy;
};//virtual machine struct

class Server {
public:
	vector<VM> vms;
	Server(int cpu, int mem);
	bool put_flavor(VM vm);
	double get_cpu_usage_rate();
	double get_mem_usage_rate();
	string printans(string str, Server server);
private:
	int total_mem;
	int total_cpu;
	int free_mem;
	int free_cpu;
	unordered_map<string, int> vm_count;
};

struct PS{
	int id;
	int ps_cpu;
	int ps_memory;
};//physical server struct

//box node
typedef struct bnode{
    int id;
    int rm_cpu;
    int rm_mem;
    int have_vm[16];
    struct bnode *next;
}BNode;

void get_predict_factor(vector<vector<int>> vm_count_per_day, string flavour, vector<double>&factor, int day_count);

void get_vm_count_per_day(vector<vector<int>>&vm_count_per_day, char * data[MAX_DATA_NUM], int data_num, int& day_count);

vector<vector<float>> get_multi_line_predict_factor(vector<vector<int>> vm_count_per_day, string flavour,vector<float>& train_mean, float sigma, int day_count);
vector<float> Gaussian_kernal(vector<std::vector<float>> &X, float sigma);
void Gaussian_kernal_2(vector<std::vector<float>> &X, float sigma, vector<float>train_mean);
vector<vector<float>> add_feature(vector<vector<float>> X, vector<float> add);

int zhishu_smooth(vector<vector<int>> vm_count_per_day, string flavour,int day_count);

bool cmp_mem(VM a, VM b);
bool cmp_cpu(VM a, VM b);
bool cmp_buy(VM a, VM b);
vector<string> split(const string& src, string separate_character);
tm StringTotm(string str);
string fillingBox(vector<VM> vm,int n, int physical_cpu, int physical_mem);
string SA_fillbox(vector<VM> vm, int n, int physical_cpu, int physical_mem, bool flag);
string int_to_string(int i);
#endif /* UTIL_H_ */
