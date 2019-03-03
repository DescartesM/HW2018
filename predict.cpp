/*
 * predict.cpp
 *
 *  Created on: 2018年3月11日
 *      Author: customer
 */
#include "predict.h"
//你要完成的功能总入口
void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
	// 需要输出的内容
	// char * result_file = (char *)"17\n\n0 8 0 20";
	/**********************1.initial (data & info)*************************/
	std::stringstream ss;
	string my_temp_str = info[0];
	vector<string> physical_info =split(my_temp_str, " ");
	int physical_cpu = atoi(physical_info[0].c_str());
	int physical_memory = atoi(physical_info[1].c_str())*1024;
	int vm_num = atoi(info[2]);
	vector<VM>vm;
	vm.resize(vm_num);
	for(int i=0;i<vm_num;i++){
		string temp_string =info[i+3];
		vector<string>temp_vm_info = split(temp_string, " ");
		vm[i].vm_name = temp_vm_info[0];
		vm[i].vm_cpu = atoi(temp_vm_info[1].c_str());
		vm[i].vm_memory = atoi(temp_vm_info[2].c_str());
	}
	string kind = info[vm_num+4];
	string start_time_predict_str = info[vm_num+6];
	string end_time_predict_str = info[vm_num+7];;
	tm start_time_predict = StringTotm(start_time_predict_str);
	tm end_time_predict = StringTotm(end_time_predict_str);
	int delta_predict_time = (mktime(&end_time_predict) - mktime(&start_time_predict))/86400;
	vector<string> start_time_train_str =split(data[0], "\t");
	tm start_time_train_tm = StringTotm(start_time_train_str[2]);
	vector<string> end_time_train_str =split(data[data_num-1], "\t");
	tm end_time_train_tm = StringTotm(end_time_train_str[2]);
	int delta_train_time = (mktime(&end_time_train_tm) - mktime(&start_time_train_tm))/86400 +1;//和预测部分不同，要加1；
	/************************2.predict part***************************/

	int day_count = 0;
	//记录每天每种flavour的销量(类型和天数都从1开始)
	vector<vector<int>>vm_count_per_day(100);
	for(unsigned int i=0;i<vm_count_per_day.size();i++)
			vm_count_per_day[i].resize(100);
	get_vm_count_per_day(vm_count_per_day, data, data_num, day_count);

	/*
	//对每一种虚拟机进行预测（用的是指数平滑）用一个星期
	for(unsigned int i=0;i<vm.size();i++){
		vm[i].vm_buy =zhishu_smooth(vm_count_per_day, vm[i].vm_name, day_count);
		vm[i].vm_buy = vm[i].vm_buy*delta_predict_time/7;
	}
	*/
	//对每一种虚拟机进行预测（用的是局部线性加权）
	for(unsigned int i=0;i<vm.size();i++){
			vector<double>pre_factor(2);
			get_predict_factor(vm_count_per_day, vm[i].vm_name, pre_factor, day_count);
			int day = day_count / 7 +1;
			vm[i].vm_buy = round(pre_factor[0] + pre_factor[1]*day);
			if(vm[i].vm_buy < 0) vm[i].vm_buy =100;
			/*
			vm[i].vm_buy = 0;
			double day_buy=0.0;
			for(int day = delta_train_time+1;day<=delta_train_time + delta_predict_time;day++){
				day_buy += (pre_factor[0] + pre_factor[1]*day);
			}
			vm[i].vm_buy = round(day_buy);
			if(vm[i].vm_buy < 0) vm[i].vm_buy =0; // 防止负数报错
			*/
	}

	//对每一种虚拟机进行预测（多元线性回归）
	/*
	for(unsigned int i=0;i<vm.size();i++){
		float sigma = 0.8;
		unsigned int feature_num = 4;//feature num (now=3, const_part && time && test_feature)
		vector<float>train_mean(feature_num);
		vector<vector<float>>w;
		w = get_multi_line_predict_factor(vm_count_per_day, vm[i].vm_name, train_mean, sigma, day_count);
		int flavour_idx =0;
		if(vm[i].vm_name.length() > 7)
			flavour_idx = (int)(vm[i].vm_name[6]- '0')*10 + (int)(vm[i].vm_name[7]- '0');
		else
			flavour_idx = (int)(vm[i].vm_name[6]- '0');
		vector<vector<float>>feature_input(1);

		feature_input[0].push_back(1);//常数b
		feature_input[0].push_back(day_count+1);//时间序号
		float my_mean = 0;
		for(int j = day_count-1;j>day_count-8;j--)
				my_mean +=vm_count_per_day[flavour_idx][j+1]/7.0;
		feature_input[0].push_back(my_mean);
		feature_input[0].push_back(my_mean*7);

		Gaussian_kernal_2(feature_input, sigma, train_mean);
		int ans_bug =0;
		ans_bug =w[0][0]*feature_input[0][0] + w[0][1]*feature_input[0][1] + w[0][2]*feature_input[0][2] + w[0][3]*feature_input[0][3];
		//vm[i].vm_buy =yjn[0][0];
		vm[i].vm_buy = ans_bug;
	}
*/
	//average predict
	/*
	int start_time_new = mktime(&start_time_predict) - 86400*7;//距离预测开始前7天
	vector<int>vm_count(100,0);
	for(int i=0;i<data_num;i++){
			vector<string> data_info =split(data[i], "\t");
			tm tm_time = StringTotm(data_info[2]);
			int temp = mktime(&tm_time);
			if(temp >= start_time_new){
					int vm_idx = 0;
					if(data_info[1].length() > 7)
						vm_idx = (int)(data_info[1][6]- '0')*10 + (int)(data_info[1][7]- '0');
					else
						vm_idx = (int)(data_info[1][6]- '0');
					vm_count[vm_idx]++;
			}
	}
	for(int i=0;i<vm_num;i++){
        int idx = 0;//the idx of flavorX
        if(vm[i].vm_name.length()>7)
        		idx = (int)(vm[i].vm_name[6]-'0')*10+(int)(vm[i].vm_name[7]-'0');
        else
        		idx = (int)(vm[i].vm_name[6]-'0');
        vm[i].vm_buy = vm_count[idx]/7.0*(delta_predict_time);
	}
	*/
	/************************3.placement part*************************/
	string ans;
	int vm_sum = 0;
	//sao cao zuo in placement
	/*
	if (kind == "CPU\r\n"){
		int vm_predict_sum = 0;
		for(int i=0;i<vm_num;i++){
			vm_predict_sum += vm[i].vm_buy*vm[i].vm_cpu;
			}
		if((vm_predict_sum%physical_cpu)/physical_cpu*1.0 < 0.15){
			int duoyu = vm_predict_sum%physical_cpu;
			sort(vm.begin(), vm.end(), cmp_buy);
			int idx=0;
			while(duoyu){
				vm[idx].vm_buy--;
				duoyu -=vm[idx].vm_cpu;
				idx++;
			}
		}
	}
*/
	for(int i=0;i<vm_num;i++){
		vm_sum +=vm[i].vm_buy;
	}
	ans = int_to_string(vm_sum) + "\n";
	for(int i=0;i<vm_num;i++){
		ans = ans + vm[i].vm_name + ' '+ int_to_string(vm[i].vm_buy) + "\n";
	}
	ans = ans +"\n";
	bool CPUorMEM;
	if (kind == "CPU\r\n")
		//sort(vm.begin(), vm.end(), cmp_cpu);
		CPUorMEM = false;
	else
		CPUorMEM = true;
		//sort(vm.begin(), vm.end(), cmp_mem);

//SA + First-Fit Algorithm
	vector<VM>new_vm;
	for(int i=0;i<vm_num;i++)
		for(int j=0;j<vm[i].vm_buy;j++){
			VM temp_vm;
			temp_vm.vm_cpu = vm[i].vm_cpu;
			temp_vm.vm_name = vm[i].vm_name;
			temp_vm.vm_memory = vm[i].vm_memory;
			new_vm.push_back(temp_vm);
	}
	string ans2;
	ans2=SA_fillbox(new_vm,vm_sum, physical_cpu, physical_memory, CPUorMEM);
	ans = ans + ans2;
	const char* result_file = ans.c_str();

	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	write_result(result_file, filename);
}




