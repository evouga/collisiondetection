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
	set<double> times;
	int nverts = verts.size();
	for(int i=0; i<nverts; i++)
	{
		for(vector<HistoryEntry>::const_iterator it = history_[verts[i]].begin(); it != history_[verts[i]].end(); ++it)
			times.insert(it->time);
	}
	stitchedHistory.clear();
	for(set<double>::const_iterator it = times.begin(); it != times.end(); ++it)
	{
		StitchedEntry newentry;
		newentry.time = *it;
		for(int i=0; i<nverts; i++)
			newentry.pos.push_back(getPosAtTime(verts[i], *it));
		stitchedHistory.push_back(newentry);
	}
}

void History::stitchCommonHistoryInterval(const std::vector<int> &verts, std::vector<StitchedEntry> &stitchedHistory, double mint, double maxt) const
{
	set<double> times;
	int nverts = verts.size();
	for(int i=0; i<nverts; i++)
	{
		for(vector<HistoryEntry>::const_iterator it = history_[verts[i]].begin(); it != history_[verts[i]].end(); ++it)
			times.insert(it->time);
	}
	stitchedHistory.clear();
	for(set<double>::const_iterator it = times.begin(); it != times.end(); ++it)
	{
		set<double>::const_iterator next = it;
		++next;
		if(next != times.end() && *next < mint)
			continue;
		StitchedEntry newentry;
		newentry.time = *it;
		for(int i=0; i<nverts; i++)
			newentry.pos.push_back(getPosAtTime(verts[i], *it));
		stitchedHistory.push_back(newentry);
		if(*it > maxt)
			break;
	}
}
