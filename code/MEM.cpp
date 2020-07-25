// Maximum Entropy Model (MEM)
// Jun-Seop Song

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdlib>
#include <cmath>

#define N 9
#define alpha 0.1
#define MAX_ITER 100000
#define tol 0.00001
#define PAIRWISE_MEM

const int M = 1 << N;
int b[N + 1]; // b[i] = 2^i;
double h[N], J[N][N];
double E[M], P[M]; // E: negative energy, P: e^E/sum(e^E)
double s[N], ss[N][N], s_m[N], ss_m[N][N];
double err;


void input()
{
	FILE *fid;
	int i, j;

	fid = fopen("s.txt", "r");
	for (i = 0; i < N; ++i) fscanf(fid, "%lf", &s[i]);
	fclose(fid);

	fid = fopen("ss.txt", "r");
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			fscanf(fid, "%lf", &ss[i][j]);
	fclose(fid);

	for (i = 0; i <= N; ++i) b[i] = 1 << i;
}


void calcProbability()
{
	int i, j, k;
	double sum = 0.0;

	// Calculate E
	for (k = 0; k < M; ++k)
	{
		E[k] = 0.0;

		for (i = 0; i < N; ++i)
		{
			if (k & b[i]) E[k] += h[i];
		}

#ifdef PAIRWISE_MEM
		for (i = 0; i < N; ++i)
		{
			for (j = i + 1; j < N; ++j)
			{
				if ((k & b[i]) && (k & b[j])) E[k] += J[i][j];
			}
		}
#endif
	}

	// Calculate P
	for (k = 0; k < M; ++k)
	{
		P[k] = exp(E[k]);
		sum += P[k];
	}
	for (k = 0; k < M; ++k) P[k] /= sum;
}


void updateParameter()
{
	int i, j, k;
	double delta;

	// Calculate s_m
	for (i = 0; i < N; ++i)
	{
		s_m[i] = 0.0;
		for (k = 0; k < M; ++k)
		{
			if (k & b[i]) s_m[i] += P[k];
		}
	}

#ifdef PAIRWISE_MEM
	// Calculate ss_m
	for (i = 0; i < N; ++i)
	{
		for (j = i + 1; j < N; ++j)
		{
			ss_m[i][j] = 0.0;
			for (k = 0; k < M; ++k)
			{
				if ((k & b[i]) && (k & b[j])) ss_m[i][j] += P[k];
			}
		}
	}
#endif

	// Update h, J
	err = 0.0;
	for (i = 0; i < N; ++i)
	{
		delta = alpha * log(s[i] / s_m[i]);
		h[i] += delta;
		if (fabs(delta) > err) err = fabs(delta);
	}
#ifdef PAIRWISE_MEM
	for (i = 0; i < N; ++i)
	{
		for (j = i + 1; j < N; ++j)
		{
			delta = alpha * log(ss[i][j] / ss_m[i][j]);
			J[i][j] += delta;
			if (fabs(delta) > err) err = fabs(delta);
		}
	}
#endif
}


void iteration()
{
	int i;

	for (i = 1; i <= MAX_ITER; ++i)
	{
		calcProbability();
		updateParameter();
		printf("%d    %lf\n", i, err);
		if (err < tol) break;
	}
}


void output()
{
	FILE *fid;
	int i, j;

	fid = fopen("h.txt", "w");
	for (i = 0; i < N; ++i) fprintf(fid, "%.10lf\n", h[i]);
	fclose(fid);

#ifdef PAIRWISE_MEM
	fid = fopen("J.txt", "w");
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			if (i < j) fprintf(fid, "%.10lf ", J[i][j]);
			else if (i > j) fprintf(fid, "%.10lf ", J[j][i]);
			else fprintf(fid, "0 ");
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
#endif

	fid = fopen("s_m.txt", "w");
	for (i = 0; i < N; ++i) fprintf(fid, "%.10lf\n", s_m[i]);
	fclose(fid);

#ifdef PAIRWISE_MEM
	fid = fopen("ss_m.txt", "w");
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			if (i < j) fprintf(fid, "%.10lf ", ss_m[i][j]);
			else if (i > j) fprintf(fid, "%.10lf ", ss_m[j][i]);
			else fprintf(fid, "0 ");
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
#endif

	fid = fopen("E.txt", "w");
	for (i = 0; i < M; ++i) fprintf(fid, "%.10lf\n", -E[i]);
	fclose(fid);

#ifdef PAIRWISE_MEM
	fid = fopen("P2.txt", "w");
#else
	fid = fopen("P1.txt", "w");
#endif
	for (i = 0; i < M; ++i) fprintf(fid, "%.10lf\n", P[i]);
	fclose(fid);
}


int main()
{
	input();
	iteration();
	output();
	system("pause");

	return 0;
}