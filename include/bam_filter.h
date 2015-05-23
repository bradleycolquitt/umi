#ifndef BAM_FILTER_H
#define BAM_FILTER_H

class BamRecord;

#include <bam_hash.h>

class StrandHash
{
    public:
        bool update(BamRecord * record);

        unordered_map<int, unique_ptr<PositionHash> >::iterator begin()
        {
            return position_map.begin();
        }

        unordered_map<int, unique_ptr<PositionHash> >::iterator end()
        {
            return position_map.end();
        }

    private:
        unordered_map<int, unique_ptr<PositionHash> > position_map;
};

#endif
