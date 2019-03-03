/*
 * pinv.hpp
 *
 *  Created on: 2018年4月7日
 *      Author: customer
 */

#ifndef PINV_H_
#define PINV_H_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <float.h>

std::vector<std::vector<float>> matrix_mul(const std::vector<std::vector<float>>& mat1, const std::vector<std::vector<float>>& mat2);

int transpose(const std::vector<std::vector<float>>& src, std::vector<std::vector<float>>& dst);

void JacobiSVD(std::vector<std::vector<float>>& At,
	std::vector<std::vector<float>>& _W, std::vector<std::vector<float>>& Vt);

int svd(const std::vector<std::vector<float>>& matSrc,
	std::vector<std::vector<float>>& matD, std::vector<std::vector<float>>& matU, std::vector<std::vector<float>>& matVt);

std::vector<std::vector<float>> addcols(std::vector<int> x);

std::vector<std::vector<float>> diag(std::vector<float> x);

void print_matrix(const std::vector<std::vector<float>>& mat);

int pinv(const std::vector<std::vector<float>>& src, std::vector<std::vector<float>>& dst, float tolerance);

#endif /* PINV_H_ */
