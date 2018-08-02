#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <map>
#include <string>
#include <iostream>

template<typename T>
class Statistics {

public:

    Statistics(std::string name) : stat_name_(name) {}
    ~Statistics() {}

    void registerStatistic(std::string name, T initial_value = 0) {
        std::pair<statistics_iterator_t, bool> statistic_position;
        statistic_position = statistics_.insert(std::pair<std::string, T>(name, initial_value));

        if(statistic_position.second == false)
            std::cout << "Statistic " << statistic_position.first->second << " already exists!" << std::endl;
    }

    void changeStatisticBy(std::string name, T change_by_value) {
        statistics_iterator_t name_at;
        name_at = statistics_.find(name);

        if(name_at != statistics_.end())
            statistics_[name] += change_by_value;
        else
            std::cout << "Statistic " << name << " not registered!" << std::endl;
    }

    T getStatisticValue(std::string name) {
        return statistics_[name];
    }

    void dumpStatistics() {
        std::cout << "Statistics: " << stat_name_ << std::endl;

        T sum = 0;
        for(std::pair<std::string, T> stat : statistics_) {
            sum += stat.second;
            std::cout << "\t" << stat.first << " = " << stat.second << std::endl;
        }

        std::cout << "_________________________" << std::endl;
        std::cout << "Total: " << sum << std::endl;
    }

    void dumpStatistics(std::ostringstream &stream) {
        stream << "Statistics: " << stat_name_ << std::endl;

        T sum = 0;
        for (std::pair<std::string, T> stat : statistics_) {
            sum += stat.second;
            stream << "\t" << stat.first << " = " << stat.second << std::endl;
        }

        stream << "_________________________" << std::endl;
        stream << "Total: " << sum << std::endl;
    }


private:

    typedef typename std::map<std::string, T> statistics_t;
    typedef typename std::map<std::string, T>::iterator statistics_iterator_t;

    statistics_t statistics_;
    std::string stat_name_;

};

#endif
