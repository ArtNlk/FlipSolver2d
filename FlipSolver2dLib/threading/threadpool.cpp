#include "threadpool.h"
#include <thread>

ThreadPool::ThreadPool(unsigned int size):
    m_threads(size),
    m_tasks(),
    m_cv(),
    m_poolDone(),
    m_mutex(),
    m_running(true),
    m_workingThreadCount(0)
{
    if(size == 0)
    {
        size = 1;
    }
    for (size_t i = 0; i < size; i++) {
        m_threads.emplace_back(&ThreadPool::threadFunc,this);
    }
}

ThreadPool::~ThreadPool()
{
    m_running = false;
    m_cv.notify_all();
}

std::vector<Range> ThreadPool::splitRange(unsigned int length, unsigned int minSize)
{
    return splitRange(0,length,minSize);
}

std::vector<Range> ThreadPool::splitRange(unsigned int start, unsigned int end, unsigned int minSize)
{
        unsigned int rangeSize = end - start;
        unsigned int numParts = m_threads.size();
        if(rangeSize < numParts)
        {
            numParts = rangeSize;
        }
        std::vector<Range> parts;
        parts.reserve(numParts);
        unsigned int partSize = std::max(rangeSize / numParts, minSize);
        unsigned int remainder = rangeSize % partSize;
        unsigned int currentStart = start;
        unsigned int currentEnd = start + partSize;
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
//            std::this_thread::yield();
//    }
    m_poolDone.wait(lock, [this] {
        //std::cout << "wakeup" << std::endl;
        //std::cout << m_running << std::endl;
        //std::cout << m_tasks.empty() << std::endl;
        return !m_running || (m_workingThreadCount == 0 && m_tasks.empty());
    });
}
