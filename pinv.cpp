/*
 * pinv.cpp
 *
 *  Created on: 2018年4月8日
 *      Author: customer
 */
//矩阵求逆反相关函数，取自github，目前还没有修改
#ifndef PINV_CPP_
#define PINV_CPP_
#include "pinv.h"
std::vector<std::vector<float>> matrix_mul(const std::vector<std::vector<float>>& mat1, const std::vector<std::vector<float>>& mat2)
{
	std::vector<std::vector<float>> result;
	int m1 = mat1.size(), n1 = mat1[0].size();
	int m2 = mat2.size(), n2 = mat2[0].size();
	if (n1 != m2) {
		fprintf(stderr, "mat dimension dismatch\n");
		return result;
	}

	result.resize(m1);
	for (int i = 0; i < m1; ++i) {
		result[i].resize(n2, (float)0);
	}

	for (int y = 0; y < m1; ++y) {
		for (int x = 0; x < n2; ++x) {
			for (int t = 0; t < n1; ++t) {
				result[y][x] += mat1[y][t] * mat2[t][x];
			}
		}
	}

	return result;
}


// ================================= 矩阵转置 =================================
int transpose(const std::vector<std::vector<float>>& src, std::vector<std::vector<float>>& dst)
{
	int m = src.size();
	int n = src[0].size();

	dst.resize(n);
	for (int i = 0; i < n; ++i) {
		dst[i].resize(m);
	}

	for (int y = 0; y < n; ++y) {
		for (int x = 0; x < m; ++x) {
			dst[y][x] = src[x][y];
		}
	}

	return 0;
}

// ================================= 矩阵奇异值分解 =================================
void JacobiSVD(std::vector<std::vector<float>>& At,
	std::vector<std::vector<float>>& _W, std::vector<std::vector<float>>& Vt)
{
	double minval = FLT_MIN;
	float eps = (float)(FLT_EPSILON * 2);
	const int m = At[0].size();
	const int n = _W.size();
	const int n1 = m; // urows
	std::vector<double> W(n, 0.);

	for (int i = 0; i < n; i++) {
		double sd{0.};
		for (int k = 0; k < m; k++) {
			float t = At[i][k];
			sd += (double)t*t;
		}
		W[i] = sd;

		for (int k = 0; k < n; k++)
			Vt[i][k] = 0;
		Vt[i][i] = 1;
	}

	int max_iter = std::max(m, 30);
	for (int iter = 0; iter < max_iter; iter++) {
		bool changed = false;
		float c, s;

		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
				float *Ai = At[i].data(), *Aj = At[j].data();
				double a = W[i], p = 0, b = W[j];

				for (int k = 0; k < m; k++)
					p += (double)Ai[k] * Aj[k];

				if (std::abs(p) <= eps * std::sqrt((double)a*b))
					continue;

				p *= 2;
				double beta = a - b, gamma = hypot((float)p, beta);
				if (beta < 0) {
					double delta = (gamma - beta)*0.5;
					s = (float)std::sqrt(delta / gamma);
					c = (float)(p / (gamma*s * 2));
				} else {
					c = (float)std::sqrt((gamma + beta) / (gamma * 2));
					s = (float)(p / (gamma*c * 2));
				}

				a = b = 0;
				for (int k = 0; k < m; k++) {
					float t0 = c*Ai[k] + s*Aj[k];
					float t1 = -s*Ai[k] + c*Aj[k];
					Ai[k] = t0; Aj[k] = t1;

					a += (double)t0*t0; b += (double)t1*t1;
				}
				W[i] = a; W[j] = b;

				changed = true;

				float *Vi = Vt[i].data(), *Vj = Vt[j].data();

				for (int k = 0; k < n; k++) {
					float t0 = c*Vi[k] + s*Vj[k];
					float t1 = -s*Vi[k] + c*Vj[k];
					Vi[k] = t0; Vj[k] = t1;
				}
			}
		}

		if (!changed||W[n-1]<1e-70)
			break;
	}

	for (int i = 0; i < n; i++) {
		double sd{ 0. };
		for (int k = 0; k < m; k++) {
			float t = At[i][k];
			sd += (double)t*t;
		}
		W[i] = std::sqrt(sd);
	}

	for (int i = 0; i < n - 1; i++) {
		int j = i;
		for (int k = i + 1; k < n; k++) {
			if (W[j] < W[k])
				j = k;
		}
		if (i != j) {
			std::swap(W[i], W[j]);

			for (int k = 0; k < m; k++)
				std::swap(At[i][k], At[j][k]);

			for (int k = 0; k < n; k++)
				std::swap(Vt[i][k], Vt[j][k]);
		}
	}

	for (int i = 0; i < n; i++)
		_W[i][0] = (float)W[i];

	srand(time(nullptr));

	for (int i = 0; i < n1; i++) {
		double sd = i < n ? W[i] : 0;

		for (int ii = 0; ii < 100 && sd <= minval; ii++) {
			// if we got a zero singular value, then in order to get the corresponding left singular vector
			// we generate a random vector, project it to the previously computed left singular vectors,
			// subtract the projection and normalize the difference.
			const float val0 = (float)(1. / m);
			for (int k = 0; k < m; k++) {
				unsigned int rng = rand() % 4294967295; // 2^32 - 1
				float val = (rng & 256) != 0 ? val0 : -val0;
				At[i][k] = val;
			}
			for (int iter = 0; iter < 2; iter++) {
				for (int j = 0; j < i; j++) {
					sd = 0;
					for (int k = 0; k < m; k++)
						sd += At[i][k] * At[j][k];
					float asum = 0;
					for (int k = 0; k < m; k++) {
						float t = (float)(At[i][k] - sd*At[j][k]);
						At[i][k] = t;
						asum += std::abs(t);
					}
					asum = asum > eps * 100 ? 1 / asum : 0;
					for (int k = 0; k < m; k++)
						At[i][k] *= asum;
				}
			}

			sd = 0;
			for (int k = 0; k < m; k++) {
				float t = At[i][k];
				sd += (double)t*t;
			}
			sd = std::sqrt(sd);
		}

		float s = (float)(sd > minval ? 1 / sd : 0.);
		for (int k = 0; k < m; k++)
			At[i][k] *= s;
	}
}

// matSrc为原始矩阵，支持非方阵，matD存放奇异值，matU存放左奇异向量，matVt存放转置的右奇异向量
int svd(const std::vector<std::vector<float>>& matSrc,
	std::vector<std::vector<float>>& matD, std::vector<std::vector<float>>& matU, std::vector<std::vector<float>>& matVt)
{
	int m = matSrc.size();
	int n = matSrc[0].size();
	for (const auto& sz : matSrc) {
		if (n != (int)sz.size()) {
			fprintf(stderr, "matrix dimension dismatch\n");
			return -1;
		}
	}

	bool at = false;
	if (m < n) {
		std::swap(m, n);
		at = true;
	}

	matD.resize(n);
	for (int i = 0; i < n; ++i) {
		matD[i].resize(1, (float)0);
	}
	matU.resize(m);
	for (int i = 0; i < m; ++i) {
		matU[i].resize(m, (float)0);
	}
	matVt.resize(n);
	for (int i = 0; i < n; ++i) {
		matVt[i].resize(n, (float)0);
	}
	std::vector<std::vector<float>> tmp_u = matU, tmp_v = matVt;

	std::vector<std::vector<float>> tmp_a, tmp_a_;
	if (!at)
		transpose(matSrc, tmp_a);
	else
		tmp_a = matSrc;

	if (m == n) {
		tmp_a_ = tmp_a;
	} else {
		tmp_a_.resize(m);
		for (int i = 0; i < m; ++i) {
			tmp_a_[i].resize(m, (float)0);
		}
		for (int i = 0; i < n; ++i) {
			tmp_a_[i].assign(tmp_a[i].begin(), tmp_a[i].end());
		}
	}
	JacobiSVD(tmp_a_, matD, tmp_v);

	if (!at) {
		transpose(tmp_a_, matU);
		matVt = tmp_v;
	} else {
		transpose(tmp_v, matU);
		matVt = tmp_a_;
	}

	return 0;
}


// 一维变二维
std::vector<std::vector<float>> addcols(std::vector<int> x)
{
	unsigned int num = x.size();
	std::vector<std::vector<float>> newx(num, std::vector<float>(2, 1));
	for (unsigned int i = 0; i<num; ++i)
	{
		newx[i][1] = x[i];
	}
	return newx;
}

// diag() - 输入向量，输出是对角矩阵  // ok!
std::vector<std::vector<float>> diag(std::vector<float> x)
{
	unsigned int num = x.size();
	unsigned int i, j;
	std::vector<std::vector<float>> w(num, std::vector<float>(num, 0));  // w初始化为零
	for (i = 0; i<num; ++i)
	{
		for (j = 0; j<num; ++j)
		{
			if (i == j)  // 对角线上的元素为向量x的元素
			{
				w[i][j] = float(x[i]);
			}
		}
	}
	return w;
}

// ================================= 输出矩阵值 =================================
void print_matrix(const std::vector<std::vector<float>>& mat)
{
	int rows = mat.size();
	for (int y = 0; y < rows; ++y) {
		for (unsigned int x = 0; x < mat[y].size(); ++x) {
			fprintf(stderr, "  %f  ", mat[y][x]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
}

// ================================= 求伪逆矩阵 =================================
int pinv(const std::vector<std::vector<float>>& src, std::vector<std::vector<float>>& dst, float tolerance)
{
	std::vector<std::vector<float>> D, U, Vt;
	if (svd(src, D, U, Vt) != 0) {
		fprintf(stderr, "singular value decomposition fail\n");
		return -1;
	}

	int m = src.size();
	int n = src[0].size();

	std::vector<std::vector<float>> Drecip, DrecipT, Ut, V;

	transpose(Vt, V);
	transpose(U, Ut);

	if (m < n)
		std::swap(m, n);

	Drecip.resize(n);
	for (int i = 0; i < n; ++i) {
		Drecip[i].resize(m, (float)0);

		if (D[i][0] > tolerance)
			Drecip[i][i] = 1.0f / D[i][0];
	}

	if (src.size() < src[0].size())
		transpose(Drecip, DrecipT);
	else
		DrecipT = Drecip;

	std::vector<std::vector<float>> tmp = matrix_mul(V, DrecipT);
	dst = matrix_mul(tmp, Ut);

	return 0;
}

/*
int main()
{
	//test_pseudoinverse();

	std::vector<int> x(150,1);   // x轴 - 时间轴
	for(int i=0;i<x.size();++i)
	{
	x[i] = i+1;
	}
	std::vector<std::vector<float>> newx = addcols(x);
	std::vector<std::vector<float>> y(150,std::vector<float>(1,0));
	y[40][0] = 2;
	y[45][0] = 1;
	y[57][0] = 2;
	y[60][0] = 1;
	y[76][0] = 1;
	y[82][0] = 1;
	y[88][0] = 1;
	y[89][0] = 2;
	y[91][0] = 1;
	y[97][0] = 1;
	y[113][0] = 1;
	y[117][0] = 1;
	y[119][0] = 1;
	y[125][0] = 1;
	y[132][0] = 1;
	y[139][0] = 1;
	y[140][0] = 2;
	y[141][0] = 2;y[145][0] = 1;y[147][0] = 2;y[149][0] = 1;

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

	tau = 0.5;

	std::vector<std::vector<float>> y_est(1,std::vector<float>(newx.size(),0));
	std::vector<float> w_ii(150, 0);   // 150x1
	std::vector<std::vector<float>> W;
	for (int ii = 0; ii < newx.size(); ++ii)
	{
		for (int j = 0; j < newx.size(); ++j)
		{
			w_ii[j] = std::exp(-std::pow((newx[ii][1]) - newx[j][1],2) / (2 * tau*tau));
		}
		W = diag(w_ii);
		pinv(matrix_mul(matrix_mul(dst,W),newx), dst3, pinvtoler);
		theta_vec = matrix_mul(matrix_mul(matrix_mul(dst3, dst), W), y);
		for (int m = 0; m < newx[0].size(); ++m)
		{
			y_est[0][ii] += newx[ii][m] * theta_vec[m][0];
		}
	}
	for (int i = 0; i < y_est.size(); ++i)
		for (int j = 0; j < y_est[0].size(); ++j)
			std::cout << y_est[i][j] <<std::endl;
	return 0;
}
*/



#endif /* PINV_CPP_ */
