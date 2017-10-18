#pragma once
#include <chrono>
#include <iostream>
#include <vector>
using namespace std;

class ChronoTimer
{
private:
	std::string name;
	int ntimers;
	std::vector<double> time;
	std::vector<int> count;
	std::vector< std::chrono::time_point<std::chrono::steady_clock> > start;
	std::chrono::time_point<std::chrono::steady_clock> last;

public:
	ChronoTimer(const std::string &name, int ntimers = 1)
	{
		this->name = name;
		this->ntimers = ntimers;
		reset();
	}

	void reset()
	{
		time.resize(ntimers);
		count.resize(ntimers);
		start.resize(ntimers);
	}

	void tic(int i = 0)
	{
		start[i] = chrono::steady_clock::now();
		last = start[i];
	}

	void toc(int i = 0)
	{

		auto end = chrono::steady_clock::now();
		
		//auto diff = end - start[i];
		auto diff = end - last;
		last = end;
		time[i] = chrono::duration<double, std::nano>(diff).count();
		++count[i];

	}

	void print() const
	{
		if (ntimers == 1) {
			printf("%s    % 8d %0.9f\n", name.c_str(), count[0], time[0] * 1e-9);
		}
		else {
			printf("%s[0] % 8d %0.9f\n", name.c_str(), count[0], time[0] * 1e-9);
			for (int i = 1; i < ntimers; ++i) {
				printf("%s[%d] % 8d %0.9f\n", name.c_str(), i, count[i], time[i] * 1e-9);
			}
		}
	}
};