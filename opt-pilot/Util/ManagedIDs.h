#include <cstddef>
#include <queue>

/**
 *  \brief Simple class to manage an ID pool.
 *
 *  Previously freed ID's are redistributed on following requests.
 */
class ManagedIDs {

public:

    ManagedIDs() : next_free_(0)
    {}

    /// return next free ID
    size_t nextID() {

        size_t id = 0;

        if(freeids_.empty()) {
            id = next_free_;
            next_free_++;
        } else {
            id = freeids_.front();
            freeids_.pop();
        }

        return id;
    }


    /// free previously allocated ID
    void freeID(size_t id) {

        if(id == next_free_ - 1)
            next_free_--;
        else
            freeids_.push(id);
    }


private:

    /// queue to handle freed ID's
    std::queue<size_t> freeids_;

    /// next free ID
    size_t next_free_;

};
