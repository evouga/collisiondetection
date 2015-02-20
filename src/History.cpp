#include "History.h"
#include <set>
#include <iostream>

using namespace Eigen;
using namespace std;

History::History(const VectorXd &qstart)
{
	int nverts = qstart.size()/3;
	for(int i=0; i<nverts; i++)
	{
		HistoryEntry firstentry;
		firstentry.time = 0;
		firstentry.pos = qstart.segment<3>(3*i);
		history_.push_back(vector<HistoryEntry>());
		history_[i].push_back(firstentry);
	}
}

void History::addHistory(int vert, double time, const Eigen::Vector3d &pos)
{
	HistoryEntry newentry;
	newentry.time = time;
	newentry.pos = pos;
	history_[vert].push_back(newentry);
}

void History::finishHistory(const Eigen::VectorXd &qend)
{
	int nverts = qend.size()/3;
	for(int i=0; i<nverts; i++)
	{
		HistoryEntry lastentry;
		lastentry.time = 1.0;
		lastentry.pos = qend.segment<3>(3*i);
		history_[i].push_back(lastentry);
	}
}

int History::countHistoryEntries()
{
	int total = 0;
	for(vector<vector<HistoryEntry> >::iterator it = history_.begin(); it != history_.end(); ++it)
		total += it->size();
	return total;
}

void History::getPosAtTime(int vert, double time, Eigen::Vector3d &pos, int &idx) const
{
	return getPosAtTime(vert, time, pos, idx, 0, history_[vert].size()-1);
}

void History::getPosAtTime(int vert, double time, Eigen::Vector3d &pos, int &idx, int start, int end) const
{
	const vector<HistoryEntry> &hist = history_[vert];
	if(end-start < 10)
	{
		int next = start;
		while(next < (int)hist.size() && hist[next].time <= time)
			next++;
		int prev = next-1;
		if(next == (int)hist.size())
		{
			pos = hist[prev].pos;
			idx = prev;
			return;
		}
		double dt = hist[next].time - hist[prev].time;
		double alpha = (time - hist[prev].time)/dt;
		pos = (1.0 - alpha)*hist[prev].pos + alpha*hist[next].pos;
		idx = prev;
		return;
	}
	int mid = (start+end)/2;
	if(hist[mid].time < time)
		return getPosAtTime(vert, time, pos, idx, mid, end);
	else
		return getPosAtTime(vert, time, pos, idx, start, mid);
}

double History::computeMinimumGap() const
{
	double result = 1.0;
	int numverts = history_.size();
	for(int i=0; i<numverts; i++)
	{
		int nentries = history_[i].size();
		for(int j=1; j<nentries; j++)
		{
			double gap = history_[i][j].time - history_[i][j-1].time;
			result = min(gap, result);
		}
	}
	return result;
}

void History::stitchCommonHistory(const std::vector<int> &verts, std::vector<StitchedEntry> &stitchedHistory) const
{
	// Assumes there is always a history entry at time 0 and 1
	int nverts = verts.size();
	stitchedHistory.clear();

	vector<vector<HistoryEntry>::const_iterator> its;
	for(int i=0; i<nverts; i++)
		its.push_back(history_[verts[i]].begin());

	double curtime = 0;
	while(curtime <= 1.0)
	{
		double newtime = numeric_limits<double>::infinity();
		StitchedEntry newentry;
		newentry.time = curtime;
		newentry.pos.resize(nverts);
		for(int i=0; i<nverts; i++)
		{
			vector<HistoryEntry>::const_iterator next = its[i];
			++next;
			while(next != history_[verts[i]].end() && next->time <= curtime)
			{
				++next;
				++its[i];
			}
			Vector3d oldpos = its[i]->pos;
			if(next == history_[verts[i]].end())
				newentry.pos[i] = oldpos;
			else
			{
				Vector3d newpos = next->pos;
				double dt = next->time - its[i]->time;
				double a = curtime - its[i]->time;
				double b = next->time - curtime;
				newentry.pos[i] = (b*oldpos + a*newpos)/dt;
				newtime = min(newtime, next->time);
			}
		}
		stitchedHistory.push_back(newentry);
		curtime = newtime;
	}
}
