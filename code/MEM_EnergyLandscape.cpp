// Energy Landscape Analysis
// Jun-Seop Song

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>

using namespace std;

#define N 9
#define min(a,b) (a)<(b)?(a):(b)

const int M = 1 << N;
int b[N + 1]; // b[i] = 2^i;
double E[M];
vector<int> minimaState;
vector<double> minimaEnergy;
int numMinima;
int basin[M], *basinCnt;
double **energyBarrier; // energy barrier between two local minima
bool isRemoved[M], isVisited[M];
vector<pair<double, int>> transitionState;


void input()
{
	FILE *fid;
	int i;

	fid = fopen("E.txt", "r");
	for (i = 0; i < M; ++i) fscanf(fid, "%lf", &E[i]);
	fclose(fid);

	for (i = 0; i <= N; ++i) b[i] = 1 << i;
	for (i = 0; i < M; ++i) basin[i] = -1;
}


void findLocalMinima()
{
	int i, j;
	bool isMin;

	printf("[Local Minima]\n");

	for (i = 0; i < M; ++i)
	{
		isMin = true;
		for (j = 0; j < N; ++j)
		{
			if (i & b[j])
			{
				if (E[i] > E[i - b[j]])
				{
					isMin = false; break;
				}
			}
			else
			{
				if (E[i] > E[i + b[j]])
				{
					isMin = false; break;
				}
			}
		}

		if (isMin)
		{
			printf("%d %d %lf\n", numMinima, i, E[i]);
			minimaState.push_back(i);
			minimaEnergy.push_back(E[i]);
			basin[i] = numMinima++;
		}
	}
}


int moveLowerEnergy(int i)
{
	int j, minAdjState;
	double minAdjEnergy = E[i];

	if (basin[i] != -1) return basin[i];

	for (j = 0; j < N; ++j)
	{
		if (i & b[j])
		{
			if (minAdjEnergy > E[i - b[j]])
			{
				minAdjState = i - b[j];
				minAdjEnergy = E[i - b[j]];
			}
		}
		else
		{
			if (minAdjEnergy > E[i + b[j]])
			{
				minAdjState = i + b[j];
				minAdjEnergy = E[i + b[j]];
			}
		}
	}

	basin[i] = moveLowerEnergy(minAdjState);
	return basin[i];
}


void calcBasinSize()
{
	int i;

	for (i = 0; i < M; ++i)
		if (basin[i] == -1) moveLowerEnergy(i);

	basinCnt = new int[numMinima];
	fill(basinCnt, basinCnt + numMinima, 0);
	for (i = 0; i < M; ++i) basinCnt[basin[i]]++;
}


bool dfs(int i, int target)
{
	int j;

	if (isVisited[i] || isRemoved[i]) return 0;
	isVisited[i] = 1;

	if (i == target) return 1;

	for (j = 0; j < N; ++j)
	{
		if (i & b[j])
		{
			if (dfs(i - b[j], target)) return 1;
		}
		else
		{
			if (dfs(i + b[j], target)) return 1;
		}
	}

	return 0;
}


void calcEnergyBarrier() // disconnectivity graph
{
	priority_queue<pair<double, int>> pq; // (energy, state)
	pair<double, int> Eth; // energy threshold
	int i, j;
	int cntConnected = numMinima * (numMinima - 1) / 2;

	printf("\n[Energy Barrier]\n");

	bool **isConnected = new bool*[numMinima];
	for (i = 0; i < numMinima; ++i) isConnected[i] = new bool[numMinima];
	for (i = 0; i < numMinima; ++i) for (j = 0; j < numMinima; ++j) isConnected[i][j] = true;

	energyBarrier = new double*[numMinima];
	for (i = 0; i < numMinima; ++i) energyBarrier[i] = new double[numMinima];
	for (i = 0; i < numMinima; ++i) for (j = 0; j < numMinima; ++j) energyBarrier[i][j] = 0.0;

	for (i = 0; i < M; ++i) pq.push(make_pair(E[i], i));

	while (cntConnected) // until all pairs are disconnected
	{
		Eth = pq.top(); pq.pop();
		isRemoved[Eth.second] = 1;

		// Determine if i, jth local minima are path-connected.
		for (i = 0; i < numMinima; ++i)
		{
			for (j = i + 1; j < numMinima; ++j)
			{
				if (isConnected[i][j])
				{
					fill(isVisited, isVisited + M, false);
					if (!dfs(minimaState[i], minimaState[j])) // disconnected
					{
						printf("E=%lf  %d %d\n", Eth.first, i, j);
						if (transitionState.empty() || transitionState.back() != Eth) transitionState.push_back(Eth);
						isConnected[i][j] = 0;
						cntConnected--;
						energyBarrier[i][j] = energyBarrier[j][i] = min(Eth.first - minimaEnergy[i], Eth.first - minimaEnergy[j]);
					}
				}
			}
		}
	}
}


void output()
{
	FILE *fid;
	int i, j;

	fid = fopen("LocalMinima.txt", "w");
	for (i = 0; i < numMinima; ++i)
	{
		fprintf(fid, "%d ", i); // #
		fprintf(fid, "%d ", minimaState[i]); // minima state
		for (j = 0; j < N; ++j) fprintf(fid, "%d", (bool)(minimaState[i] & b[j])); // minima state (binary)
		fprintf(fid, " %.10lf", minimaEnergy[i]); // minima energy
		fprintf(fid, " %.10lf\n", (double)basinCnt[i] / (double)M); // basin size
	}
	fclose(fid);

	fid = fopen("EnergyBarrier.txt", "w");
	for (i = 0; i < numMinima; ++i)
	{
		for (j = 0; j < numMinima; ++j) fprintf(fid, "%.10lf ", energyBarrier[i][j]);
		fprintf(fid, "\n");
	}
	fclose(fid);

	fid = fopen("TransitionState.txt", "w");
	for (i = 0; i < transitionState.size(); ++i)
	{
		fprintf(fid, "%d ", transitionState[i].second); // minima state
		for (j = 0; j < N; ++j) fprintf(fid, "%d", (bool)(transitionState[i].second & b[j])); // minima state (binary)
		fprintf(fid, " %.10lf\n", transitionState[i].first);
	}
	fclose(fid);
}


int main()
{
	input();
	findLocalMinima();
	calcBasinSize();
	calcEnergyBarrier();
	output();
	system("pause");

	return 0;
}