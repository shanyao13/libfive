/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once
#include <future>
#include <list>
#include <vector>

namespace libfive {

    //ProgressHandler 类的主要功能是处理任务进度的跟踪和更新，并在需要时通知用户或系统进度的变化。
    // 它通过多线程机制定期检查进度状态，并可以根据需要在不同阶段（Phase）中更新进度信息
class ProgressHandler
{
public:
    /*
     *  Constructor
     */
    ProgressHandler();
    virtual ~ProgressHandler();

    /*  Implement this function in a subclass to send progress updates.
     *
     *  The default implementation does nothing, but needs to exist
     *  so that we don't fail if the thread tries to call it during
     *  the destructor. */
    virtual void progress(double d);

    /*  Called by workers to update current progress */
    void tick(uint64_t i=1);

    /*  start() kicks off the worker thread for a progress operation
     *  with the given number of phases and per-phase weights. */
    void start(const std::vector<unsigned>& weights);
    void nextPhase(uint64_t total);

    /*  finish() is called by the same place that calls start(),
     *  and is responsible for joining the worker thread. */
    void finish();

protected:
    /*  This function is run in an async thread.  It checks the progress
     *  regularly and calls progress() when the value changes. */
    void run();

    struct Phase {
        Phase(unsigned weight) :
            weight(weight), total(0), counter(0)
        {
            /* nothing to do here */
        }
        unsigned weight;
        uint64_t total;
        std::atomic<uint64_t> counter;
    };
    std::list<Phase> phases;
    std::list<Phase>::iterator current_phase;
    unsigned total_weight;

    /*  Flag used to stop the progress thread in the destructor */
    std::atomic_bool done;

    /*  Used to stop the worker thread early, since it otherwise only
     *  updates at a fixed speed (e.g. 200 Hz). */
    std::timed_mutex timed_mut;
    std::mutex phase_mut;

    /*  This is our async update thread */
    std::future<void> future;
};

}   // namespace libfive
