/*
 * predict.h
 *
 *  Created on: 2018年3月11日
 *      Author: customer
 */

#ifndef __ROUTE_H__
#define __ROUTE_H__
#include "util.h"


using namespace std;

void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename);
vector<string> split(const string& src, string separate_character);
tm StringTotm(string str);



#endif /* PREDICT_H_ */
