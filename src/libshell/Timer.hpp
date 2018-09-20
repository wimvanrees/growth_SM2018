//
//  Timer.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/25/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//
#ifndef Timer_hpp
#define Timer_hpp

#include "common.hpp"
#include <chrono>


/*! \class Timer
 * \brief Simple timer using c++11 chrono
 * \tparam one of std::chrono::duration, such as std::chrono::seconds or std::chrono::milliseconds
 *
 * Usage: Timer<std::chrono::seconds> or Timer<std::chrono::milliseconds>
 */
template<typename T = std::chrono::seconds>
class Timer
{
    typedef std::chrono::high_resolution_clock tClock;
    typedef std::chrono::duration<long double, typename T::period> tDuration;
    
    long unsigned int nSamples;
    tDuration totalDuration;
    std::chrono::time_point<tClock> t_start, t_end;
    
public:
    Timer():
    nSamples(0),
    totalDuration(tDuration(0))
    {}
    
    void start()
    {
        t_start = tClock::now();
    }
    
    void stop()
    {
        t_end = tClock::now();
        totalDuration = totalDuration + std::chrono::duration_cast<tDuration>(t_end - t_start);
        nSamples++;
    }
    
    void reset()
    {
        nSamples=0;
        totalDuration = tDuration(0);
    }
    
    long double getTotalTime_s() const
    {
        return std::chrono::duration_cast<std::chrono::duration<long double, std::chrono::seconds::period>>(totalDuration).count();
    }
    
    long double totalTime_ms() const
    {
        return std::chrono::duration_cast<std::chrono::duration<long double, std::chrono::milliseconds::period>>(totalDuration).count();
    }
    
    long unsigned int getSamples() const
    {
        return nSamples;
    }
};

#endif /* Timer_hpp */
