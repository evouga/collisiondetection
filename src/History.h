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
	
	void addHistory(int vert, double time, const Eigen::Vector3d &pos);
	void finishHistory(const Eigen::VectorXd &qend);

	int countHistoryEntries();
	void getPosAtTime(int vert, double time, Eigen::Vector3d &pos, int &idx) const;
	void stitchCommonHistory(const std::vector<int> &verts, std::vector<StitchedEntry> &stitchedHistory) const;
	const std::vector<HistoryEntry> &getVertexHistory(int vert) const {return history_[vert];}
	double computeMinimumGap() const;
private:
	bool atEnd(const std::vector<std::vector<HistoryEntry>::const_iterator> &its, const std::vector<int> &verts) const;
	void getPosAtTime(int vert, double time, Eigen::Vector3d &pos, int &idx, int begin, int end) const;
	
	std::vector<std::vector<HistoryEntry> > history_;
};

#endif
