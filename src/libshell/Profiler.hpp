//
//  Profiler.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/25/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Profiler_hpp
#define Profiler_hpp

#include <Timer.hpp>

class Profiler
{
    std::map<std::string, Timer<> *> timers;
    std::string lastTimerName;
public:
    
    Profiler(): timers()
    {}
    
    void push_start(const std::string & sTimerName)
    {
        const auto it = timers.find(sTimerName);
        
        // create a new timer if does not exists -- else we start again
        if(it == timers.end())
        {
            Timer<> * timer = new Timer<>();
            timers[sTimerName] = timer;
        }
        
        timers.at(sTimerName)->start();
        lastTimerName = sTimerName;
    }
    
    void pop_stop()
    {
        pop_stop(lastTimerName);
    }
    
    void pop_stop(const std::string & sTimerName)
    {
        const auto it = timers.find(sTimerName);
        assert(it != timers.end());
        
        if(it != timers.end())
            it->second->stop();
    }
    
    void printSummary() const
    {
        Real totalAllTimers = 0.0;
        for(const auto & timer : timers)
        {
            totalAllTimers += timer.second->getTotalTime_s();
        }
        for(const auto & timer : timers)
        {
            const std::string timer_name = timer.first;
            const auto nSamples = timer.second->getSamples();
            const Real totalTime = timer.second->getTotalTime_s();
            const double avgTime = (nSamples == 0 ? -1 : totalTime / nSamples);
            printf("[%20s]: \t %02.0f%%\t%03.3e s\t%03.3f s\t(%lu samples)\n",
                   timer_name.c_str(), 100*totalTime/totalAllTimers, avgTime, totalTime,  nSamples);
        }
        printf("[Total time]: \t%f\n", totalAllTimers);
    }
    
    void printSummaryToFile(FILE * f) const
    {
        assert(f!=NULL);
        
        Real totalAllTimers = 0.0;
        for(const auto & timer : timers)
        {
            totalAllTimers += timer.second->getTotalTime_s();
        }
        for(const auto & timer : timers)
        {
            const std::string timer_name = timer.first;
            const auto nSamples = timer.second->getSamples();
            const Real totalTime = timer.second->getTotalTime_s();
            const Real avgTime = (nSamples == 0 ? -1 : totalTime / nSamples);
            fprintf(f, "[%20s]: \t %02.0f%%\t%03.3e s\t%03.3f s\t(%lu samples)\n",
                    timer_name.c_str(), 100*totalTime/totalAllTimers, avgTime, totalTime,  nSamples);
        }
        fprintf(f, "[Total time]: \t%f\n", totalAllTimers);
    }
    
    void reset()
    {
        for(auto & timer : timers)
            timer.second->reset();
    }
    
    void clear()
    {
        for(auto & timer : timers)
        {
            delete (timer.second);
            timer.second = nullptr;
        }
        timers.clear();
    }
    
    ~Profiler()
    {
        clear();
    }
};

#endif /* Profiler_hpp */
