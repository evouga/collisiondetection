#ifndef HISTORY_H
#define HISTORY_H

#include <Eigen/Core>
#include <vector>

struct HistoryEntry
{
	double time;
	Eigen::Vector3d pos;

	bool operator==(const HistoryEntry &other)
	{
		return time == other.time && pos == other.pos;
	}
};

struct StitchedEntry
{
	double time;
	std::vector<Eigen::Vector3d> pos;
};

class History
{
public:
	History(const Eigen::VectorXd &qstart);
	
	void addHistory(int vert, double time, const Eigen::Vector3d &pos, const History *oldhistory, bool checkOldHistory);
	void finishHistory(const Eigen::VectorXd &qend);

	int countHistoryEntries();
	const Eigen::Vector3d getPosAtTime(int vert, double time) const;
	void stitchCommonHistory(const std::vector<int> &verts, std::vector<StitchedEntry> &stitchedHistory) const;
	const std::vector<HistoryEntry> &getVertexHistory(int vert) const {return history_[vert];}

private:
	std::vector<std::vector<HistoryEntry> > history_;
};

#endif
