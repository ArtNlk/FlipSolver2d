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
}

ThreadPool *ThreadPool::i()
{
    static ThreadPool instance;
    return &instance;
}

std::vector<Range> ThreadPool::splitRange(unsigned int length, unsigned int minSize)
{
    return splitRange(0,length,minSize);
}

std::vector<Range> ThreadPool::splitRange(unsigned int start, unsigned int end, unsigned int minSize)
{
    if(start == end)
    {
        return std::vector<Range>();
    }
    unsigned int rangeSize = end - start;
    unsigned int numParts = m_threads.size();
    if(rangeSize < numParts)
    {
        numParts = rangeSize;
    }
    std::vector<Range> parts;
    parts.reserve(numParts);
    unsigned int partSize = std::floor(rangeSize / numParts);
    partSize = std::max(minSize,partSize);
    unsigned int remainder = 0;
    if(partSize * numParts < rangeSize)
    {
        remainder = rangeSize - (partSize * numParts);
    }
    unsigned int currentStart = start;
    unsigned int currentEnd = start + partSize;
    currentEnd = std::min(end,currentEnd);
    for (unsigned int i = 0; i < numParts; i++) {
        if (i < remainder) {
            currentEnd++;
        }
        parts.emplace_back(Range{currentStart, currentEnd});
        currentStart = currentEnd;
        currentEnd = std::min(currentStart + partSize, end);
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
