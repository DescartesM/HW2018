/*
 * util.cpp
 *
 *  Created on: 2018年4月8日
 *      Author: customer
 */

#include "util.h"
//各种功能性函数+预测函数+放置函数
//int to string
string int_to_string(int i)
{
    char ch[10];
    sprintf(ch, "%d", i);
    string s(ch);

    return s;
}
//server class func
Server::Server(int cpu, int mem) :total_cpu(cpu), total_mem(mem) {
	free_cpu = cpu;
	free_mem = mem;
}
bool Server::put_flavor(VM vm) {
	if (free_cpu >= vm.vm_cpu && free_mem >= vm.vm_memory) {
		free_cpu -= vm.vm_cpu;
		free_mem -= vm.vm_memory;
		unordered_map<string,int>::iterator key;
		key = vm_count.find(vm.vm_name);
		if(key!=vm_count.end())
				vm_count[vm.vm_name] ++;
		else
			vm_count[vm.vm_name] = 1;
		vms.push_back(vm);
		return true;
	}
	return false;
}
string Server::printans(string str, Server server){
	unordered_map<string,int>::iterator iter = server.vm_count.begin();
	for(;iter!=server.vm_count.end();iter++)
		str = str + " "+ iter->first + " " + int_to_string(iter->second);
	return str;
}
double Server::get_cpu_usage_rate() {
	return (1-free_cpu/(total_cpu*1.0));
}
double Server::get_mem_usage_rate() {
	return (1-free_mem/(total_mem*1.0));
}
//sa algorithm
string SA_fillbox(vector<VM> vm, int n, int physical_cpu, int physical_mem, bool flag){
	double used_server_num = vm.size() + 1;//相当于temp_max
	vector<Server> res_servers;
	double T_init = 100.0;  //SA初始温度
	double T_min = 1;   //SA终止温度
	double d_speed = 0.9999; //SA下降速度
	vector<int> gods_dice;//上帝骰子，制造变化
	for (unsigned int i = 0; i < vm.size(); i++) {
		gods_dice.push_back(i);
	}
	//core part
	int T = T_init;
	while (T > T_min) {
			random_shuffle(gods_dice.begin(), gods_dice.end());//STL 自带洗牌算法，乱序
			vector<VM> new_vm_copy = vm;
			swap(new_vm_copy[gods_dice[0]], new_vm_copy[gods_dice[1]]);
			//FF Algorithm
			vector<Server> servers;
			Server first_server(physical_cpu, physical_mem);
			servers.push_back(first_server);
			unsigned int i;
			for(i=0;i<new_vm_copy.size();i++){
				auto iter = servers.begin();
				for (; iter != servers.end(); ++iter) {
					if (iter->put_flavor(new_vm_copy[i])) {
						break;
					}
				}
				if (iter == servers.end()) {//装不下了，打开新箱子
					Server new_server(physical_cpu, physical_mem);
					new_server.put_flavor(new_vm_copy[i]);
					servers.push_back(new_server);
				}
			}
			double server_num;
			//cpu=0.mem=1
			if (flag == 0)
				server_num = servers.size() - 1 + servers.rbegin()->get_cpu_usage_rate();
			else
				server_num = servers.size() - 1 + servers.rbegin()->get_mem_usage_rate();
			//如果分数更低，则保存结果
			if (server_num < used_server_num) {
				used_server_num = server_num;
				res_servers = servers;
				vm = new_vm_copy;
			}
			//如果分数更高，则以一定概率保存结果，防止优化陷入局部最优解
			else {
				if (exp((used_server_num - server_num) / T) > rand() / RAND_MAX) {
					used_server_num = server_num;
					res_servers = servers;
					vm = new_vm_copy;
				}
			}
			T = d_speed * T;  //一次循环结束，温度降低
	}
	unsigned int Bnum = res_servers.size();
	string placement_str;
	unsigned int i;
	for(i=0;i<Bnum;i++){
		placement_str = placement_str + int_to_string(i+1);
		placement_str = res_servers[i].printans(placement_str, res_servers[i]);
		if(i != Bnum-1)
			placement_str = placement_str + "\n";
	}
	string res;
	res = int_to_string(Bnum) + "\n" + placement_str;
  return res;
}
//ffd algorithm
string fillingBox(vector<VM> vm,int n, int physical_cpu, int physical_mem)
{
		string placement_str;
    BNode *h=NULL;
    BNode *p=NULL;
    BNode *tail=NULL;
    int i;
    for(i=0;i<n;i++)
        {
    	int idx =1;
    		//遍历所有开了的箱子，找到第一个能放下去的箱子
    	for(p=h;p&&(p->rm_cpu<vm[i].vm_cpu||p->rm_mem<vm[i].vm_memory);p=p->next)
    		{
    			idx++; // used as physical server's ID
    		}
    	if(p==NULL)
            {
				//需要开新箱子
		    p=(BNode *)malloc(sizeof(BNode));
		    p->id = idx;
		    p->rm_cpu = physical_cpu;
		    p->rm_mem = physical_mem;
			  p->next=NULL;
			  for(int i=0;i<16;i++)
				      p->have_vm[i] = 0;
        if(!h)
                  //因为无头结点需要开新箱子
          h=tail=p;
        else
                 {
                    //因为现有的所有箱子剩余容积均不够需开新箱子,update tail pointer
         tail=tail->next=p;
                 }
            }
            //向确定要放入物品的箱子（不论新旧）放入物品
    	int vm_id = 0;
    	if(vm[i].vm_name.length() > 7)
    			vm_id = 10 + (int)(vm[i].vm_name[7]- '0');
    	else
    			vm_id = (int)(vm[i].vm_name[6]- '0');
     p->rm_cpu -= vm[i].vm_cpu;
     p->rm_mem -= vm[i].vm_memory;
     p->have_vm[vm_id]++;
     //printf("put vm flavor%d in Physical machine %d\n",vm_id, idx);
    }
	int Bnum = 0;
	while((h)){
		//get box num
		Bnum++;
		placement_str = placement_str + int_to_string(h->id);
		for(int i=1;i<16;i++){
				if(h->have_vm[i] != 0)
						placement_str = placement_str + " flavor"+ int_to_string(i) + " " + int_to_string(h->have_vm[i]);
		}
		if(h->next)
				placement_str = placement_str + "\n";
		h = h->next;
	}

	string res;
	res = int_to_string(Bnum) + "\n" + placement_str;
  return res;
}
//string to tm
tm StringTotm(string str)//只精确到天数，方便之后做统计
{
    char *cha = (char*)str.data();             // 将string转换成char*。
    tm tm_;                                    // 定义tm结构体。
    int year, month, day, hour, minute, second;// 定义时间的各个int临时变量。
    sscanf(cha, "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second);// 将string存储的日期时间，转换为int临时变量。
    tm_.tm_year = year - 1900;                 // 年，由于tm结构体存储的是从1900年开始的时间，所以tm_year为int临时变量减去1900。
    tm_.tm_mon = month - 1;                    // 月，由于tm结构体的月份存储范围为0-11，所以tm_mon为int临时变量减去1。
    tm_.tm_mday = day;                         // 日。
    tm_.tm_hour = 0;                        // 时。
    tm_.tm_min = 0;                       // 分。
    tm_.tm_sec = 0;                       // 秒。
    tm_.tm_isdst = 0;                          // 非夏令时。
    return tm_;                                 // 返回值。
}
//split string
vector<string> split(const string& src, string separate_character)
{
    vector<string> strs;

  int separate_characterLen = separate_character.size();//分割字符串的长度,这样就可以支持如“,,”多字符串的分隔符
    int lastPosition = 0,index = -1;
    while (-1 != (index = src.find(separate_character,lastPosition)))
    {
        strs.push_back(src.substr(lastPosition,index - lastPosition));
        lastPosition = index + separate_characterLen;
    }
    string lastString = src.substr(lastPosition);//截取最后一个分隔符后的内容
    if (!lastString.empty())
        strs.push_back(lastString);//如果最后一个分隔符后还有内容就入队
    return strs;
}
bool cmp_cpu(VM a, VM b){
	//primary key is cpu
	if(a.vm_cpu < b.vm_cpu)
		return false;
	else if (a.vm_cpu == b.vm_cpu){
		if(a.vm_memory < b.vm_memory)
			return false;
		else return true;
	}
	return true;
}
bool cmp_mem(VM a, VM b){
	//primary key is memory
	if(a.vm_memory < b.vm_memory)
			return false;
		else if (a.vm_memory == b.vm_memory){
			if(a.vm_cpu < b.vm_cpu)
				return false;
			else return true;
		}
		return true;
}
bool cmp_buy(VM a, VM b){
	if(a.vm_buy < b.vm_buy)
			return false;
	else return true;
}
//获得每天每种flavor的销售额
void get_vm_count_per_day(vector<vector<int>>&vm_count_per_day, char * data[MAX_DATA_NUM], int data_num, int& day_count){
	int day_idx = 0;
	int vm_idx = 0;
	vector<string> temp =split(data[0], "\t");
	tm start_time_tm = StringTotm(temp[2]);
	int start_time = mktime(&start_time_tm);
	for(int i=0;i<data_num;i++){
			vector<string> data_info =split(data[i], "\t");
			tm tm_time = StringTotm(data_info[2]);
			int temp_time = mktime(&tm_time);
			day_idx=(temp_time- start_time)/86400 +1;
			if(data_info[1].length() > 7)
					vm_idx = (int)(data_info[1][6]- '0')*10 + (int)(data_info[1][7]- '0');
			else
					vm_idx = (int)(data_info[1][6]- '0');
			vm_count_per_day[vm_idx][day_idx]++;

	}
	day_count = day_idx;
}
//获得局部线性回归的参数
void get_predict_factor(vector<vector<int>> vm_count_per_day, string flavour, vector<double>&factor, int day_count){
	int flavour_idx =0;
	if(flavour.length() > 7)
		flavour_idx = (int)(flavour[6]- '0')*10 + (int)(flavour[7]- '0');
	else
		flavour_idx = (int)(flavour[6]- '0');
	unsigned int i=0;
	int week_num = day_count / 7;
	int dropout = day_count % 7;
	std::vector<std::vector<float>> y(week_num,std::vector<float>(1,0));
	int idx =0;
	for(i=dropout;i<day_count;i+=7){
		int week_sum =0;
		for(int j= i;j<i+7;j++){
			week_sum += vm_count_per_day[flavour_idx][j+1];
		}
		y[idx][0] = week_sum;
		idx++;
	}
	//vector<int> x(day_count,1);   // x轴 - 时间轴
	vector<int> x(week_num,1);
	for(unsigned int i=0;i<x.size();++i)
		x[i] = i+1;
	std::vector<std::vector<float>> newx = addcols(x);
	/*
	//std::vector<std::vector<float>> y(day_count,std::vector<float>(1,0));
	for(int i=0;i<day_count;i++){
		y[i][0] = vm_count_per_day[flavour_idx][i+1];
		if(y[i][0] >= 15)
			y[i][0] /=1.5;
	}
	*/
	std::vector<std::vector<float>> theta_vec;//(2, std::vector<float>(1,0));
	float  pinvtoler = 1.e-6;
	float tau = 1;
	std::vector<std::vector<float>> dst;
	std::vector<std::vector<float>> dst2;
	std::vector<std::vector<float>> dst3; // 用于后面更新权重

	transpose(newx, dst);  // dst是转置的结果2x150 x'

	dst2 = matrix_mul(dst, newx);

	std::vector<std::vector<float>> pinv1;
	pinv(dst2, pinv1, pinvtoler); // pinv1: 2x2

	theta_vec = matrix_mul(matrix_mul(pinv1, dst), y);

	tau = 2.2;//best for now

	std::vector<std::vector<float>> y_est(1,std::vector<float>(newx.size(),0));
	std::vector<float> w_ii(newx.size(), 0);   // 150x1
	std::vector<std::vector<float>> W;
	for (unsigned int ii = 0; ii < newx.size(); ++ii)
	{
		for (unsigned int j = 0; j < newx.size(); ++j)
		{
			w_ii[j] = std::exp(-std::pow((newx[ii][1]) - newx[j][1],2) / (2 * tau*tau));
		}
		W = diag(w_ii);
		pinv(matrix_mul(matrix_mul(dst,W),newx), dst3, pinvtoler);
		theta_vec = matrix_mul(matrix_mul(matrix_mul(dst3, dst), W), y);
		for (unsigned int m = 0; m < newx[0].size(); ++m)
		{
			y_est[0][ii] += newx[ii][m] * theta_vec[m][0];
		}
	}
	//output factors
	factor[0] = theta_vec[0][0];
	factor[1] = theta_vec[1][0];
	//unit test
	/*
	for (int i = 0; i < y_est.size(); ++i)
		for (int j = 0; j < y_est[0].size(); ++j)
			std::cout <<"flavour"<<flavour_idx<<"\t"<< j<<"\t"<<y_est[i][j]<<"\t"<< y[j][0]<<std::endl;
	*/
}
vector<float> Gaussian_kernal(vector<std::vector<float>> &X, float sigma){
		int rows = X.size();
		int cols =X[0].size();
		//get mean
		int i, j;
		vector<float>train_mean(cols, 0);
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				train_mean[j] += X[i][j] / (float)rows;
			}
		}
		vector<std::vector<float>> X_decen;
		X_decen = X;
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				X_decen[i][j] -= train_mean[j];
			}
		}
		//X_decen第二列除1000 ??
		for(i=0;i<rows;i++){
			X_decen[i][1] /=1000;
		}
		float d =0.5/(sigma)*(sigma);
		vector<std::vector<float>> a;
		a = X;
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				a[i][j] = X_decen[i][j]*X_decen[i][j]*d;
				X[i][j] = exp(a[i][j]); //return X back
			}
		}
		return train_mean;
}
void Gaussian_kernal_2(vector<std::vector<float>> &X, float sigma, vector<float>train_mean){
		int rows = X.size();
		int cols =X[0].size();
		//get mean
		int i, j;
		vector<std::vector<float>> X_decen;
		X_decen = X;
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				X_decen[i][j] -= train_mean[j];
			}
		}
		//X_decen第二列除1000 ??
		for(i=0;i<rows;i++){
			X_decen[i][1] /=1000;
		}
		float d =0.5/(sigma)*(sigma);
		vector<std::vector<float>> a;
		a = X;
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				a[i][j] = X_decen[i][j]*X_decen[i][j]*d;
				X[i][j] = exp(a[i][j]); //return X back
			}
		}
}
vector<vector<float>> add_feature(vector<vector<float>> X, vector<float> add){
	if(X.size() != add.size())
	{
		printf("feature size not match!");
		return X;
	}
	vector<vector<float>> new_X;
	new_X = X;
	for(unsigned int i=0;i<X.size();i++){
		new_X[i].push_back(add[i]);
	}
	return new_X;
}
vector<vector<float>> get_multi_line_predict_factor(vector<vector<int>> vm_count_per_day, string flavour,vector<float>& train_mean, float sigma, int day_count){
		int flavour_idx =0;
		if(flavour.length() > 7)
			flavour_idx = (int)(flavour[6]- '0')*10 + (int)(flavour[7]- '0');
		else
			flavour_idx = (int)(flavour[6]- '0');
		unsigned int i=0;
		int week_num = day_count / 7;
		int dropout = day_count % 7;
		std::vector<std::vector<float>> y(week_num,std::vector<float>(1,0));
		int idx =0;
		for(int i=dropout;i<day_count;i+=7){
			int week_sum =0;
			for(int j= i;j<i+7;j++){
				week_sum += vm_count_per_day[flavour_idx][j+1];
			}
			y[idx][0] = week_sum;
			idx++;
		}
		vector<int> x(week_num,1);   // x轴 - 时间轴
		for(i=0;i<x.size();++i)
			x[i] = i+1;
		std::vector<std::vector<float>> newx = addcols(x);
		//add feature function
		vector<float> temp(newx.size(),0);
		temp[0] = 0;
		for(i=1;i<newx.size();i++){
			temp[i] = y[i-1][0] / 7.0;
		}
		newx =add_feature(newx, temp);
		for(i=1;i<newx.size();i++){
			temp[i] = y[i-1][0];
		}
		newx =add_feature(newx, temp);
		/*
		for(int i=0;i<day_count;i++){
			y[i][0] = vm_count_per_day[flavour_idx][i+1];
			if(y[i][0] >= 15)
				y[i][0] /=1.5;
		}
		*/
		float  pinvtoler = 1.e-12;
		std::vector<std::vector<float>> dst;
		std::vector<std::vector<float>> dst2;
		std::vector<std::vector<float>> dst3;
		train_mean = Gaussian_kernal(newx, sigma);

		transpose(newx, dst);  // dst是 newx'
		dst2 = matrix_mul(newx, dst); // dst2是 newx'*newx

		std::vector<std::vector<float>> pinv1;
		pinv(dst2, pinv1, pinvtoler); // pinv1
		std::vector<std::vector<float>> b;
		std::vector<std::vector<float>> w;
		b = matrix_mul(pinv1, y);
		transpose(b, dst3);// dst3是 b'
		w =matrix_mul(dst3, newx);
		return w;
}

int zhishu_smooth(vector<vector<int>> vm_count_per_day, string flavour,int day_count){
	//get week
	int flavour_idx =0;
	if(flavour.length() > 7)
		flavour_idx = (int)(flavour[6]- '0')*10 + (int)(flavour[7]- '0');
	else
		flavour_idx = (int)(flavour[6]- '0');
	unsigned int i=0;
	int week_num = day_count / 7;
	int dropout = day_count % 7;
	vector<float> y;
	for(int i=dropout;i<day_count;i+=7){
		int week_sum =0;
		for(int j= i;j<i+7;j++){
			week_sum += vm_count_per_day[flavour_idx][j+1];
		}
		y.push_back(week_sum);
	}
	float mean_val =(y[0] + y[1]+ y[2])/3.0;
	vector<float> S;
	S.push_back(mean_val);
	float alpha = 0.3;
	for(i =0;i<y.size();i++){
		float temp_S;
		temp_S =alpha*y[i] + (1-alpha)*S[i];
		S.push_back(temp_S);
	}
	return round(S[week_num]);
}


