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

const Eigen::Vector3d History::getPosAtTime(int vert, double time) const
{
	vector<HistoryEntry>::const_iterator prev = history_[vert].begin();
	for(vector<HistoryEntry>::const_iterator it = history_[vert].begin(); it != history_[vert].end() && it->time <= time; ++it)
		prev = it;
	vector<HistoryEntry>::const_iterator next = prev + 1;
	if(next == history_[vert].end())
		return prev->pos;
	assert(prev->time <= time && next->time > time);
	double dt = next->time - prev->time;
	if(dt == 0)
		return next->pos;
	double alpha = (time-prev->time)/dt;
	return (1.0-alpha)*prev->pos + alpha*next->pos;
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
