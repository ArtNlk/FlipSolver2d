#include "threadpool.h"
#include <cmath>
#include <iostream>
#include <mutex>
#include <thread>

ThreadPool::ThreadPool():
    m_threads(),
    m_running(true),
    m_workingThreadCount(0)
{
    unsigned int size = std::thread::hardware_concurrency();
    //size = 1;
    for (size_t i = 0; i < size; i++) {
        m_threads.emplace_back(&ThreadPool::threadFunc,this);
    }
    std::cout << "Threadpool created with size " << m_threads.size() << std::endl;
}

ThreadPool::~ThreadPool()
{
    m_running = false;
    m_cv.notify_all();
    for (std::thread& t : m_threads)
    {
        t.join();
    }
}

ThreadPool *ThreadPool::i()
{
    static ThreadPool instance;
    return &instance;
}

std::vector<Range> ThreadPool::splitRange(size_t length, size_t minSize, size_t jobsPerThread)
{
    return splitRange(Range(0,length),minSize, jobsPerThread);
}

std::vector<Range> ThreadPool::splitRange(Range r, size_t minSize, size_t jobsPerThread)
{
    if(r.size() == 0)
    {
        return std::vector<Range>();
    }

    size_t numParts = m_threads.size() * jobsPerThread;
    if(r.size() < numParts)
    {
        numParts = r.size();
    }
    std::vector<Range> parts;
    parts.reserve(numParts);
    size_t partSize = std::floor(r.size() / numParts);
    partSize = std::max(minSize,partSize);
    size_t remainder = 0;
    if(partSize * numParts < r.size())
    {
        remainder = r.size() - (partSize * numParts);
    }
    size_t currentStart = r.start;
    size_t currentEnd = r.start + partSize;
    currentEnd = std::min(r.end,currentEnd);
    for (size_t i = 0; i < numParts; i++) {
        if (i < remainder) {
            currentEnd++;
        }
        parts.emplace_back(currentStart, currentEnd);
        currentStart = currentEnd;
        currentEnd = std::min(currentStart + partSize, r.end);
        if(currentEnd == currentStart) break;
    }
    parts.shrink_to_fit();
    return parts;
}

void ThreadPool::wait()
{
    std::unique_lock<std::mutex> lock(m_mutex);
    //std::cout << "entered wait"<< std::endl;
//    while(m_workingThreadCount != 0 || !m_tasks.empty())
//    {
//        std::this_thread::yield();
//    }
    m_poolDone.wait(lock, [this] {
        //std::cout << "wakeup" << std::endl;
        //std::cout << m_running << std::endl;
        //std::cout << m_tasks.empty() << std::endl;
        return !m_running || (m_workingThreadCount == 0 && m_tasks.empty());
    });
}

int ThreadPool::activeThreadCount()
{
    return m_workingThreadCount;
}

int ThreadPool::queueLength()
{
    std::unique_lock lock(m_mutex);
    int length = m_tasks.size();
    return length;
}

size_t ThreadPool::threadCount()
{
    return m_threads.size();
}
